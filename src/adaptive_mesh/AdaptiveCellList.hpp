/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * Shadowfax is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Shadowfax is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with Shadowfax. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file AdaptiveCellList.hpp
 *
 * @brief Specialized list for mesh evolution algorithm
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ADAPTIVECELLLIST_HPP
#define ADAPTIVECELLLIST_HPP

#include "AdaptiveMeshUtils.hpp"         // for get_periodic_position
#include "DelCont.hpp"                   // for DelCont
#include "Error.hpp"                     // for my_exit
#include "MPIGlobal.hpp"                 // for size, rank, sendbuffer, etc
#include "MPIMethods.hpp"                // for MyMPI_Pack, MyMPI_Gather, etc
#include "RestartFile.hpp"               // for RestartFile
#include "StateVector.hpp"               // for StateVector
#include "Vec.hpp"                       // for Vec
#include "utilities/Cuboid.hpp"          // for Cuboid
#include "utilities/GasParticle.hpp"     // for GasParticle
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include <algorithm>                     // for min
#include <cstddef>                       // for NULL
#include <ext/alloc_traits.h>
#include <iostream>       // for operator<<, basic_ostream, etc
#include <unordered_map>  // for unordered_map, etc
#include <utility>        // for pair
#include <vector>         // for vector

/**
 * @brief Specialized list for AdaptiveCell which can distinguish between
 * normal cells, ghost cells and orphan cells.
 *
 * The idea of this class is to provide a general way to store cells that
 * can be either normal cells, ghost cells or orphan cells and that can
 * be indexed using a single integer index. Using a single index requires us
 * to somehow encode extra information in this index and this class handles
 * this.
 * Outside of this class, we can just use the indices as if this class
 * is just a std::vector. However, we can also ask this class if a certain
 * index represents a normal cell, a ghost or an orphan. How we distinguish
 * between these three types is an internal business.
 */
template <class Cell> class AdaptiveCellList {
  private:
    /*! @brief Internal list of normal cells */
    std::vector<Cell*> _cells;

    /*! @brief Internal list of ghost cells */
    std::vector<Cell*> _ghosts;

    /*! @brief Internal list of orphan cells */
    std::vector<Cell*> _orphans;

    /*! @brief Internal list of flags indicating whether a certain orphan exists
     *  on this process or not */
    std::vector<std::vector<int> > _orphanflags;

  public:
    /**
     * @brief Constructor.
     *
     * Initializes a list with the given size for normal cells. Adds a dummy
     * orphan to the orphan list. This way, orphans can be addressed as
     * negative indices (where we just take the opposite of the index to locate
     * the orphan in the internal list), normal cells are indexed as positive
     * integers in the range 0 <= i < given size, while ghosts have positive
     * indices >= given size.
     *
     * @param normal_cell_size The number of normal cells that will be stored
     * in this list.
     */
    inline AdaptiveCellList(unsigned int normal_cell_size) {
        _cells.resize(normal_cell_size, NULL);
        _orphans.push_back(NULL);
        _orphanflags.resize(MPIGlobal::size);
    }

    inline ~AdaptiveCellList() {
        //        for(unsigned int i = 0; i < _cells.size(); i++){
        //            delete _cells[i];
        //        }
        //        for(unsigned int i = 0; i < _ghosts.size(); i++){
        //            delete _ghosts[i];
        //        }
        //        for(unsigned int i = 0; i < _orphans.size(); i++){
        //            delete _orphans[i];
        //        }
    }

    /**
     * @brief Add a normal cell to the AdaptiveCellList, with the given index
     *
     * @param index Index that the given cell should get in the AdaptiveCellList
     * @param normal_cell Pointer to the cell that will be stored. The
     * AdaptiveCellList takes ownership of the cell and will delete it upon
     * destruction
     */
    inline void add_normal_cell(int index, Cell* normal_cell) {
        _cells[index] = normal_cell;
    }

    /**
     * @brief Add the given ghost cell to the end of the AdaptiveCellList
     *
     * @param ghost_cell Pointer to the cell that will be stored. The
     * AdaptiveCellList takes ownership of the cell and will delete it upon
     * destruction.
     * @return The integer index with which the ghost cell can be retrieved
     */
    inline int add_ghost_cell(Cell* ghost_cell) {
        _ghosts.push_back(ghost_cell);
        return _ghosts.size() - 1 + _cells.size();
    }

    /**
     * @brief Add an orphan cell to the beginning of the AdaptiveCellList
     *
     * @param orphan_cell Pointer to the cell that will be stored. The
     * AdaptiveCellList takes ownership of the cell and will delete it upon
     * destruction.
     * @return The (negative) integer index with which the orphan cell can be
     * retrieved
     */
    inline int add_orphan_cell(Cell* orphan_cell) {
        _orphans.push_back(orphan_cell);
        unsigned int oid = orphan_cell->get_original();
        int rank = orphan_cell->get_rank();
        if(_orphanflags[rank].size() <= oid) {
            _orphanflags[rank].resize(oid + 1, 0);
            _orphanflags[rank][oid] = -(_orphans.size() - 1);
        }
        return -(_orphans.size() - 1);
    }

    /**
     * @brief Get the index of the orphan corresponding to the cell with the
     * given index on the process with the given rank
     *
     * If the corresponding orphan does not (yet) exist, index 0 is returned.
     *
     * @param rank Rank of the process holding the original copy of the cell
     * @param oid Index of the original copy of the cell on its home process
     * @return Index of the orphan if it exists, 0 otherwise
     */
    inline int get_orphan(int rank, unsigned int oid) {
        if(_orphanflags[rank].size() <= oid) {
            return 0;
        }
        return _orphanflags[rank][oid];
    }

    /**
     * @brief Check if the given index corresponds to a ghost cell
     *
     * Ghost cells have positive indices >= _cells.size().
     *
     * @param index The index which is queried
     * @return True if the given index corresponds to a ghost, false otherwise
     */
    inline bool is_ghost(int index) {
        return (index >= 0 && index >= (int)_cells.size());
    }

    /**
     * @brief Check if the given index corresponds to an orphan cell
     *
     * Orphan cells have strictly negative indices
     *
     * @param index The index which is queried
     * @return True if the given index corresponds to an orphan, false otherwise
     */
    inline bool is_orphan(int index) {
        return index < 0;
    }

    /**
     * @brief Check if the given index corresponds to a normal cell
     *
     * Normall cells have positive indices < _cells.size().
     *
     * @param index The index which is queried
     * @return True if the given index corresponds to a normal cell, false
     * otherwise
     */
    inline bool is_normal(int index) {
        return (index >= 0 && index < (int)_cells.size());
    }

    /**
     * @brief Access the orphan with the given index
     *
     * This function does not check whether the given index corresponds to a
     * valid orphan and might give rise to segmentation faults.
     *
     * @param index A valid negative integer index of an orphan cell
     * @return A reference to the pointer that is stored in the list
     */
    inline Cell*& get_orphan(int index) {
        return _orphans[-index];
    }

    /**
     * @brief Access the ghost with the given index
     *
     * This function does not check whether the given index corresponds to a
     * valid ghost and might give rise to segmentation faults.
     *
     * @param index A valid positive integer index of a ghost cell
     * @return A reference to the pointer that is stored in the list
     */
    inline Cell*& get_ghost(int index) {
        return _ghosts[index - _cells.size()];
    }

    /**
     * @brief Access the normal cell with the given index
     *
     * This function does not check whether the given index corresponds to a
     * valid normal cell and might give rise to segmentation fault.
     *
     * @param index A valid positive integer index of a normal cell
     * @return A reference to the pointer that is stored in the list
     */
    inline Cell*& get_normal(int index) {
        return _cells[index];
    }

    /**
     * @brief Access the cell with the given index
     *
     * This function uses the value of the index (i) to discriminate between
     *  - normal cells (0 <= i < _cells.size())
     *  - ghost cells (_cells.size() <= i)
     *  - orphan cells (i < 0)
     * If the given index does not correspond to any of these types, this
     * function will still give rise to segmentation faults.
     *
     * @param index A valid index of a cell of any type
     * @return A reference to the pointer that is stored in the list
     */
    inline Cell*& get_cell(int index) {
        if(is_orphan(index)) {
            return get_orphan(index);
        } else {
            if(is_ghost(index)) {
                return get_ghost(index);
            } else {
                return get_normal(index);
            }
        }
    }

    /**
     * @brief Operator wrapper around AdaptiveCellList::get_cell()
     *
     * This function uses the value of the index (i) to discriminate between
     *  - normal cells (0 <= i < _cells.size())
     *  - ghost cells (_cells.size() <= i)
     *  - orphan cells (i < 0)
     * If the given index does not correspond to any of these types, this
     * function will still give rise to segmentation faults.
     *
     * @param index A valid index of a cell of any type
     * @return A reference to the pointer that is stored in the list
     */
    inline Cell*& operator[](int index) {
        return get_cell(index);
    }

    /**
     * @brief Get the number of normal cells in the list
     *
     * @return The unsigned integer number of normal cells in the list
     */
    inline unsigned int normalsize() {
        return _cells.size();
    }

    /**
     * @brief Get the number of ghost cells in the list
     *
     * @return The unsigned integer number of ghost cells in the list
     */
    inline unsigned int ghostsize() {
        return _ghosts.size();
    }

    /**
     * @brief Access the last ghost cell in the list
     *
     * @return A reference to the last ghost cell that was added to the list
     */
    inline Cell*& ghostback() {
        return _ghosts.back();
    }

    /**
     * @brief Get the number of orphan cells in the list
     *
     * @return The unsigned integer number of orphans cells in the list
     */
    inline unsigned int orphansize() {
        return _orphans.size();
    }

    /**
     * @brief Reorder the normal cells in the list using the given list of new
     * indices
     *
     * The method first makes an internal copy of the old normal cell list and
     * then uses the given list of indices to assign new pointers to the
     * internal normal cell list, based on the old values stored.
     *
     * @param new_ids A std::vector of unsigned integer new indices for the
     * normal cells in the list
     */
    inline void reorder(std::vector<unsigned int>& new_ids) {
        std::vector<Cell*> old_values(_cells);
        for(unsigned int i = 0; i < _cells.size(); i++) {
            _cells[new_ids[i]] = old_values[i];
        }
    }

    /**
     * @brief Remove all orphans from the list, except those with the given
     * indices
     *
     * The method first makes an internal copy of the old orphan cell list and
     * then resizes the internal orphan list to the size of the given index
     * list + 1. The first position is blocked to allow easy indexing of the
     * new list. We then use the given list of indices to assign new values to
     * the internal orphan list, based on the copied old orphan cell list.
     *
     * @param to_keep A std::vector of negative integer indices corresponding to
     * the orphans that should not be removed
     */
    inline void update_orphans(std::vector<int>& to_keep) {
        std::vector<Cell*> orphanclone(_orphans);
        _orphans.resize(to_keep.size() + 1);
        _orphans[0] = NULL;
        for(unsigned int i = 0; i < to_keep.size(); i++) {
            _orphans[i + 1] = orphanclone[-to_keep[i]];
        }
    }

    /**
     * @brief Iterator to loop over ghost cell members of the AdaptiveCellList
     */
    class ghostiterator {
      private:
        /*! @brief Index of the ghost the iterator is currently pointing to */
        unsigned int _current;

        /*! @brief Reference to the AdaptiveCellList over which we iterate */
        AdaptiveCellList& _list;

      public:
        /**
         * @brief Constructor.
         *
         * @param current Index of the ghost the iterator is currently pointing
         * to
         * @param list Reference to the AdaptiveCellList over which we iterate
         */
        inline ghostiterator(unsigned int current, AdaptiveCellList& list)
                : _current(current), _list(list) {}

        /**
         * @brief Dereference operator
         *
         * @return Reference to the ghost cell the iterator is currently
         * pointing to
         */
        inline Cell*& operator*() {
            return _list._ghosts[_current];
        }

        /**
         * @brief Pointer operator
         *
         * @return Reference to the ghost cell the iterator is currently
         * pointing to
         */
        inline Cell*& operator->() {
            return _list._ghosts[_current];
        }

        /**
         * @brief Increment operator
         *
         * @return Reference to the incremented ghostiterator
         */
        inline ghostiterator& operator++() {
            ++_current;
            return *this;
        }

        /**
         * @brief Increment operator
         *
         * @return New instance of the ghostiterator which is incremented
         */
        inline ghostiterator operator++(int) {
            ghostiterator clone(*this);
            ++_current;
            return clone;
        }

        /**
         * @brief Compare ghostiterators
         *
         * @param it ghostiterator to compare with
         * @return True if both ghostiterators iterate over the same
         * AdaptiveCellList and have the same current position, false otherwise
         */
        inline bool operator==(ghostiterator it) {
            return (_current == it._current && &_list == &it._list);
        }

        /**
         * @brief Compare ghostiterators
         *
         * @param it ghostiterator to compare with
         * @return True if the ghostiterators are not equal as determined by
         * ghostiterator::operator==, false otherwise
         */
        inline bool operator!=(ghostiterator it) {
            return !(*this == it);
        }

        /**
         * @brief Get the index in the AdaptiveCellList of the ghost cell this
         * iterator is currently pointing to
         *
         * This is not the same index as the one that is stored internally, but
         * the one that can be used in e.g. AdaptiveCellList::get_ghost() to
         * access the ghost cell
         *
         * @return Integer index of the ghost cell in the AdaptiveCellList
         */
        inline int index() {
            return _list._cells.size() + _current;
        }
    };

    /**
     * @brief Get an AdaptiveCellList::ghostiterator to the first ghost in the
     * list
     *
     * @return An AdaptiveCellList::ghostiterator pointing to the first ghost
     * cell in the list
     */
    inline ghostiterator ghostbegin() {
        return ghostiterator(0, *this);
    }

    /**
     * @brief Get an AdaptiveCellList::ghostiterator to the end of the ghost
     * list
     *
     * @return An AdaptiveCellList::ghostiterator pointing to the end of the
     * ghost cell list
     */
    inline ghostiterator ghostend() {
        return ghostiterator(_ghosts.size(), *this);
    }

    /**
     * @brief Iterator to loop over orphan cell members of the AdaptiveCellList
     */
    class orphaniterator {
      private:
        /*! @brief Index of the orphan the iterator is currently pointing to */
        unsigned int _current;

        /*! @brief Reference to the AdaptiveCellList over which we iterate */
        AdaptiveCellList& _list;

      public:
        /**
         * @brief Constructor.
         *
         * @param current Index of the orphan the iterator is currently pointing
         * to
         * @param list Reference to the AdaptiveCellList over which we iterate
         */
        inline orphaniterator(unsigned int current, AdaptiveCellList& list)
                : _current(current), _list(list) {}

        /**
         * @brief Dereference operator
         *
         * @return Reference to the orphan the iterator is currently pointing to
         */
        inline Cell*& operator*() {
            return _list._orphans[_current];
        }

        /**
         * @brief Pointer operator
         *
         * @return Reference to the orphan the iterator is currently pointing to
         */
        inline Cell*& operator->() {
            return _list._orphans[_current];
        }

        /**
         * @brief Increment operator
         *
         * @return Reference to the incremented orphaniterator
         */
        inline orphaniterator& operator++() {
            ++_current;
            return *this;
        }

        /**
         * @brief Increment operator
         *
         * @return New instance of the orphaniterator which is incremented
         */
        inline orphaniterator operator++(int) {
            orphaniterator clone(*this);
            ++_current;
            return clone;
        }

        /**
         * @brief Compare orphaniterators
         *
         * @param it orphaniterator to compare with
         * @return True if both orphaniterators iterate over the same
         * AdaptiveCellList and their internal state is equal, false otherwise
         */
        inline bool operator==(orphaniterator it) {
            return (_current == it._current && &_list == &it._list);
        }

        /**
         * @brief Compare orphaniterators
         *
         * @param it orphaniterator to compare with
         * @return True if the orphaniterators are not equal, as determined by
         * orphaniterator::operator==, false otherwise
         */
        inline bool operator!=(orphaniterator it) {
            return !(*this == it);
        }

        /**
         * @brief Get the index in the AdaptiveCellList of the orphan cell this
         * iterator is currently pointing to
         *
         * This is not the same index as the one that is stored internally, but
         * the one that can be used in e.g. AdaptiveCellList::get_ghost() to
         * access the ghost cell
         *
         * @return Integer index of the orphan in the AdaptiveCellList
         */
        inline int index() {
            return -_current;
        }
    };

    /**
     * @brief Get an AdaptiveCellList::orphaniterator to the first orphan in the
     * list
     *
     * @return An AdaptiveCellList::orphaniterator to the first orphan in the
     * list
     */
    inline orphaniterator orphanbegin() {
        // we start counting from 1, since we reserved the first spot
        return orphaniterator(1, *this);
    }

    /**
     * @brief Get an AdaptiveCellList::orphaniterator to the end of the orphan
     * list
     *
     * @return An AdaptiveCellList::orphaniterator to the end of the orphan list
     */
    inline orphaniterator orphanend() {
        return orphaniterator(_orphans.size(), *this);
    }

    /**
     * @brief Iterator to loop over normal cell members of the AdaptiveCellList
     */
    class normaliterator {
      private:
        /*! @brief Index of the normal cell the iterator is currently pointing
         *  to */
        unsigned int _current;

        /*! @brief Reference to the AdaptiveCellList over which we iterate */
        AdaptiveCellList& _list;

      public:
        /**
         * @brief Constructor.
         *
         * @param current Index of the normal cell the iterator is currently
         * pointing to
         * @param list Reference to the AdaptiveCellList over which we iterate
         */
        inline normaliterator(unsigned int current, AdaptiveCellList& list)
                : _current(current), _list(list) {}

        /**
         * @brief Dereference operator
         *
         * @return Reference to the normal cell the iterator is currently
         * pointing to
         */
        inline Cell*& operator*() {
            return _list._cells[_current];
        }

        /**
         * @brief Pointer operator
         *
         * @return Reference to the normal cell the iterator is currently
         * pointing to
         */
        inline Cell*& operator->() {
            return _list._cells[_current];
        }

        /**
         * @brief Increment operator
         *
         * @return Reference to the incremented normaliterator
         */
        inline normaliterator& operator++() {
            ++_current;
            return *this;
        }

        // this version corresponds to iterator++
        // as you can see, much more work is done here
        // conclusion: use ++iterator instead if you don't need the return
        // value!!
        /**
         * @brief Increment operator
         *
         * @return New instance of then normaliterator which is incremented
         */
        inline normaliterator operator++(int) {
            normaliterator clone(*this);
            ++_current;
            return clone;
        }

        /**
         * @brief Compare normaliterators
         *
         * @param it normaliterator to compare with
         * @return True if both normaliterators iterate over the same
         * AdaptiveCellList and their internal states are equal, false otherwise
         */
        inline bool operator==(normaliterator it) {
            return (_current == it._current && &_list == &it._list);
        }

        /**
         * @brief Compare normaliterators
         *
         * @param it normaliterator to compare with
         * @return True if the normaliterators are not equal, as determined by
         * normaliterator::operator==, false otherwise
         */
        inline bool operator!=(normaliterator it) {
            return !(*this == it);
        }

        /**
         * @brief Get the index of the normal cell in the AdaptiveCellList
         *
         * For normal cells, this happens to be the same index as the internal
         * state of the iterator.
         *
         * @return The index of the normal cell in the AdaptiveCellList
         */
        inline int index() {
            return _current;
        }
    };

    /**
     * @brief Get an AdaptiveCellList::normaliterator to the first normal cell
     * in the list
     *
     * @return An AdaptiveCellList::normaliterator to the first normal cell in
     * the list
     */
    inline normaliterator normalbegin() {
        return normaliterator(0, *this);
    }

    /**
     * @brief Get an AdaptiveCellList::normaliterator to the end of the normal
     * cell list
     *
     * @return An AdaptiveCellList::normaliterator to the end of the normal cell
     * list
     */
    inline normaliterator normalend() {
        return normaliterator(_cells.size(), *this);
    }

    /**
     * @brief Dump the AdaptiveCellList to the given RestartFile, using the
     * given list of particles for the ordering
     *
     * @param rfile RestartFile to which the internal data are dumped
     * @param particles std::vector of particles which will be used to determine
     * the order in which normal cells are dumped to the RestartFile
     */
    template <class P>
    inline void dump(RestartFile& rfile, std::vector<P*>& particles) {
        unsigned int vsize = _cells.size();
        rfile.write(vsize);
        // we need to dump the cells in particle order, since we need to restore
        // the particle pointer in the cell when opening the restart file
        for(unsigned int i = 0; i < vsize; i++) {
            _cells[particles[i]->get_local_id()]->dump(rfile);
        }
        vsize = _ghosts.size();
        rfile.write(vsize);
        for(unsigned int i = 0; i < vsize; i++) {
            _ghosts[i]->dump(rfile);
        }
        vsize = _orphans.size() - 1;
        rfile.write(vsize);
        for(unsigned int i = 1; i < vsize + 1; i++) {
            _orphans[i]->dump(rfile);
        }
    }

    /**
     * @brief Restart constructor. Load an AdaptiveCellList from the given
     * RestartFile, using the given list of particles for ordering
     *
     * @param rfile RestartFile from which to read the internal data
     * @param particles std::vector of particles which will be used to determine
     * the order in which normal cells are added to the internal list
     */
    template <class P>
    inline AdaptiveCellList(RestartFile& rfile, std::vector<P*>& particles) {
        unsigned int vsize;
        rfile.read(vsize);
        _cells.resize(vsize);
        for(unsigned int i = 0; i < vsize; i++) {
            _cells[particles[i]->get_local_id()] =
                    new Cell(rfile, particles[i]);
        }
        rfile.read(vsize);
        _ghosts.resize(vsize);
        for(unsigned int i = 0; i < vsize; i++) {
            _ghosts[i] = new Cell(rfile);
        }
        rfile.read(vsize);
        _orphans.resize(vsize + 1);
        _orphans[0] = NULL;
        for(unsigned int i = 1; i < vsize + 1; i++) {
            _orphans[i] = new Cell(rfile);
        }
    }
};

/**
 * @brief Auxiliary class to obtain the position of a general neighbour cell
 *
 * In principle, the position can be obtained directly from the id of the
 * particle associated with the neighbour. It is however possible that the
 * neighbour is a ghost, either periodic or reflective, or that the neighbour
 * particle is on another MPI process. In these cases, the position needs to be
 * calculated or communicated first.
 */
class ParticleInfo {
  private:
    /*! @brief Rank of the process that stores the original particle */
    int _rank;

    /*! @brief Wall flag associated with the neighbour */
    int _wall;

    /*! @brief Id of the particle */
    unsigned long _id;

    /*! @brief Flag indicating if the position stored is reliable */
    bool _position_set;

    /*! @brief Position of the neighbour */
    Vec _position;

    /*! @brief Position of the neighbour in the interval [1,2] */
    Vec _position12;

    /*! @brief Local index of the cell/particle on its original process */
    unsigned int _locid;

    /*! @brief Primitive variables of the cell/particle */
    StateVector _W;

    /*! @brief Start time of the particle timestep on the integer timeline */
    unsigned long _starttime;

  public:
    /**
     * @brief Empty constructor
     *
     * Initializes an empty ParticleInfo instance
     */
    ParticleInfo()
            : _rank(-1), _wall(0), _id(0), _position_set(false), _locid(0) {}

    /**
     * @brief General constructor
     *
     * @param id Id of the associated particle
     * @param rank Rank of the MPI process where the original particle is stored
     * (defaults to the current process)
     * @param wall Wall flag associated with the neighour (defaults to 0: no
     * wall)
     */
    ParticleInfo(unsigned long id, int rank = MPIGlobal::rank, int wall = 0)
            : _rank(rank), _wall(wall), _id(id), _position_set(false),
              _locid(0) {}

    /**
     * @brief Get the id of the associated particle
     *
     * @return Id of the particle
     */
    unsigned long get_id() {
        return _id;
    }

    /**
     * @brief Set the id of the particle
     *
     * @param id Id of the particle
     */
    void set_id(unsigned long id) {
        _id = id;
    }

    /**
     * @brief Get the rank of the MPI process on which the particle is stored
     *
     * @return Rank of the owner process
     */
    int get_rank() {
        return _rank;
    }

    /**
     * @brief Set the rank of the owner process
     *
     * @param rank Rank of the owner process
     */
    void set_rank(int rank) {
        _rank = rank;
    }

    /**
     * @brief Test whether the stored position is reliable
     *
     * @return True if the position is reliable, false if it is not set or is
     * outdated
     */
    bool has_position() {
        return _position_set;
    }

    /**
     * @brief Set the position of the particle
     *
     * @param position Position of the particle
     * @param position12 Position of the particle in the interval [1,2]
     */
    void set_position(Vec position, Vec position12) {
        _position = position;
        _position12 = position12;
        _position_set = true;
    }

    /**
     * @brief Flag the position to be unreliable
     */
    void unset_position() {
        _position_set = false;
    }

    /**
     * @brief Get the particle position
     *
     * @warning This routine does not check if the position is reliable, this
     * should be done using ParticleInfo::has_position.
     *
     * @return The position associated with the neighbour
     */
    Vec get_position() {
        return _position;
    }

    /**
     * @brief Get the particle position in the interval [1,2]
     *
     * @return Particle position within the interval [1,2]
     */
    Vec get_position12() {
        return _position12;
    }

    /**
     * @brief Set the index of the cell or particle on the local process that
     * holds the original copy
     *
     * @param locid Index of the cell or particle on its original process
     */
    void set_locid(unsigned int locid) {
        _locid = locid;
    }

    /**
     * @brief Get the index of the cell or particle on the local process that
     * holds the original copy
     *
     * @return Index of the cell or particle on its original process
     */
    unsigned int get_locid() {
        return _locid;
    }

    /**
     * @brief Set the primitive variables of the particle
     *
     * @param W Primitive variables of the particle
     */
    void set_W(StateVector W) {
        _W = W;
    }

    /**
     * @brief Get the primitive variables of the particle
     *
     * @return Primitive variables of the particle
     */
    StateVector get_W() {
        return _W;
    }

    /**
     * @brief Set the start time of the particle timestep on the integer
     * timeline
     *
     * @param starttime Start time of the particle timestep
     */
    void set_starttime(unsigned long starttime) {
        _starttime = starttime;
    }

    /**
     * @brief Get the start time of the particle timestep on the integer
     * timeline
     *
     * @return Start time of the particle timestep
     */
    unsigned long get_starttime() {
        return _starttime;
    }
};

/**
 * @brief General hashmap that stores neighbour information
 *
 * Every neighbour is identified by a unique id, that is either a unique
 * particle id or a local ghost identifier. This class is responsible for making
 * sure that we can obtain a reliable position for every id, irrespective of
 * what kind of neighbour is actually considered or where the information is
 * stored.
 */
class NewCellList {
  private:
    /*! @brief Hashmap linking ids to ParticleInfo instances */
    std::unordered_map<unsigned long, ParticleInfo> _hashmap;

    /*! @brief Internal counter holding the last ghost id that was used */
    unsigned long _maxid;

    /*! @brief Cuboid specifying the dimensions of the simulation box */
    Cuboid _box;

  public:
    /**
     * @brief Constructor
     *
     * @param maxid Maximum id of all particles across all MPI processes. It is
     * assumed that ids larger than this can be safely used to reference ghost
     * copies of particles
     * @param box Cuboid specifying the dimensions of the simulation box
     */
    NewCellList(unsigned long maxid, Cuboid& box) : _maxid(maxid), _box(box) {}

    /**
     * @brief Add a ghost particle to the list
     *
     * @param original Id of the particle that is copied by the ghost
     * @param wall Flag indicating how the original coordinates of the particle
     * should be converted into actual coordinates
     * @param rank Rank of the process that holds the original particle
     * @return Unique id of the newly created ghost copy in the list
     */
    unsigned long add_ghost(unsigned long original, int wall, int rank = -1) {
        unsigned long newid = _maxid + original * 8 - 1 - wall;
        if(!_hashmap.count(newid)) {
            _hashmap[newid] = ParticleInfo(original, rank, -1 - wall);
        }
        return newid;
    }

    /**
     * @brief Get the ID of the given ghost copy of the given ID
     *
     * @param original ID of a particle
     * @param wall Flag indicating how the original coordinates of the particle
     * should be converted into actual coordinates
     * @return Unique id of the ghost copy
     */
    unsigned long get_ghost(unsigned long original, int wall) {
        unsigned long newid = _maxid + original * 8 - 1 - wall;
        return newid;
    }

    /**
     * @brief Signal that we want a reliable position for the neighbour with the
     * given id
     *
     * The routine checks if the hashmap already contains an entry for the given
     * particle. If not, it is assumed that the particle is on the current
     * MPI process.
     *
     * @param id Unique identifier of a neighbour
     * @param rank Rank of the MPI process that holds the original particle
     */
    void obtain_position(unsigned long id, int rank = -1) {
        if(!_hashmap.count(id)) {
            _hashmap[id] = ParticleInfo(id, rank);
        }
        _hashmap[id].unset_position();
    }

    /**
     * @brief Unset all recorded positions
     */
    void reset_positions() {
        for(auto it = _hashmap.begin(); it != _hashmap.end(); ++it) {
            it->second.unset_position();
        }
    }

    /**
     * @brief Make sure all particles in the list have reliable positions
     *
     * @param particles ParticleVector holding the local particles
     * @param box Cuboid used to convert coordinates to the [1,2] interval
     */
    void calculate_positions(ParticleVector& particles, Cuboid& box) {
        // do local particles and check for ghosts
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            GasParticle* p = particles.gas(i);
            unsigned long pid = p->id();
            if(!_hashmap[pid].has_position()) {
                Vec p12(p->get_position());
                for(unsigned int i = 0; i < ndim_; i++) {
                    p12[i] =
                            1. +
                            (p12[i] - box.get_anchor()[i]) / box.get_sides()[i];
                }
                _hashmap[pid].set_position(p->get_position(), p12);
                _hashmap[pid].set_locid(i);
                _hashmap[pid].set_rank(MPIGlobal::rank);
                _hashmap[pid].set_W(p->get_Wvec());
                _hashmap[pid].set_starttime(p->get_starttime());
            }
            for(int wall = 0; wall < 8; wall++) {
                unsigned long gid = _maxid + pid * 8 + wall;
                if(_hashmap.count(gid)) {
                    if(!_hashmap[gid].has_position()) {
                        // have to apply wall
                        double pos[2];
                        pos[0] = p->get_position().x();
                        pos[1] = p->get_position().y();
                        Cuboid box = particles.get_container().get_cuboid();
                        AdaptiveMeshUtils::get_periodic_position(pos, -wall - 1,
                                                                 box);
#if ndim_ == 3
                        Vec cpos(pos[0], pos[1], 0.);
#else
                        Vec cpos(pos[0], pos[1]);
#endif
                        Vec p12(cpos);
                        for(unsigned int i = 0; i < ndim_; i++) {
                            p12[i] = 1. +
                                     (p12[i] - box.get_anchor()[i]) /
                                             box.get_sides()[i];
                        }
                        _hashmap[gid].set_position(cpos, p12);
                        _hashmap[gid].set_locid(i);
                        _hashmap[gid].set_rank(MPIGlobal::rank);
                        _hashmap[gid].set_W(p->get_Wvec());
                    }
                }
            }
        }

        if(MPIGlobal::size == 1) {
            // ready!
            return;
        }

        // flag incomplete particles on other processes
        std::vector<std::vector<unsigned long> > incomplete(MPIGlobal::size);
        std::vector<unsigned long> orphans;
        for(auto it = _hashmap.begin(); it != _hashmap.end(); ++it) {
            if(!it->second.has_position()) {
                unsigned long id = it->first;
                // we only request the real particle, not its copies
                // in fact, this means we can request the same particle multiple
                // times...
                id = get_id(id);
                int rank = -1;
                if(_hashmap.count(id)) {
                    rank = _hashmap[id].get_rank();
                }
                if(rank == MPIGlobal::rank) {
                    // particle is flagged to be local, but isn't: it moved to
                    // another process
                    // flag it as unfound
                    it->second.set_rank(-1);
                    rank = -1;
                }
                if(rank >= 0) {
                    incomplete[rank].push_back(id);
                } else {
                    // if we do not know the rank, we will have to find it by a
                    // collective search
                    orphans.push_back(id);
                }
            }
        }

        // send requests to other processes
        std::vector<MPI_Request> reqs((MPIGlobal::size - 1) * 2,
                                      MPI_REQUEST_NULL);
        // Offsets in the MPIGlobal buffers. For every external process, we
        // reserve 2 blocks in each buffer to allow for non-blocking
        // communication of buffers (if we would use a single buffer, then
        // subsequent sends would overwrite it)
        std::vector<unsigned int> buffers((MPIGlobal::size - 1) * 2);
        unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
        for(unsigned int i = 0; i < buffers.size(); i++) {
            buffers[i] = i * bufsize;
        }

        // two vector components, an unsigned index and a byte flag
        unsigned int typesize = 2 * sizeof(double) + sizeof(unsigned int) + 1;
        unsigned int maxsize = bufsize / typesize;
        if(!(bufsize % typesize)) {
            maxsize--;
        }

        vector<int> allflag(MPIGlobal::size - 1);
        vector<unsigned int> numsend(MPIGlobal::size - 1);
        vector<unsigned int> sendsizes(MPIGlobal::size - 1);
        // we keep track of the total number of messages that should be sent
        // and received. We send and receive at least one message for every
        // process. For every message that is sent, an answer has to be
        // received. For incomplete sends, we need to send additional
        // messages. For incomplete receives, we have to receive additional
        // messages.
        int numtoreceive = MPIGlobal::size - 1;
        int numtosend = MPIGlobal::size - 1;
        int numrecv = 0;
        int numsent = 0;
        for(int i = 0; i < MPIGlobal::size - 1; i++) {
            sendsizes[i] =
                    incomplete[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                            .size();

            int send_pos = 0;
            for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]);
                si++) {
                MyMPI_Pack(
                        &incomplete[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                   [si],
                        1, MPI_UNSIGNED_LONG,
                        &MPIGlobal::sendbuffer[buffers[2 * i]], bufsize,
                        &send_pos);
            }
            allflag[i] = (sendsizes[i] <= maxsize);
            if(!allflag[i]) {
                numtosend++;
            }
            numsend[i] = 1;
            // add continuation signal
            MyMPI_Pack(&allflag[i], 1, MPI_INT,
                       &MPIGlobal::sendbuffer[buffers[2 * i]], bufsize,
                       &send_pos);

            MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[2 * i]], send_pos,
                        MPI_PACKED, (MPIGlobal::rank + 1 + i) % MPIGlobal::size,
                        0, &reqs[2 * i]);
            numsent++;
            numtoreceive++;
        }

        // receive loop: receive 2 messages from every external process
        vector<unsigned int> numreceived(MPIGlobal::size - 1, 0);
        while(numrecv < numtoreceive || numsent < numtosend) {
            MPI_Status status;
            if(numsent < numtosend) {
                for(int j = 0; j < MPIGlobal::size - 1; j++) {
                    if(!allflag[j]) {
                        int flag;
                        MyMPI_Test(&reqs[j], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[j] * maxsize;
                                si < std::min((numsend[j] + 1) * maxsize,
                                              sendsizes[j]);
                                si++) {
                                MyMPI_Pack(
                                        &incomplete[(MPIGlobal::rank + 1 + j) %
                                                    MPIGlobal::size][si],
                                        1, MPI_UNSIGNED_LONG,
                                        &MPIGlobal::sendbuffer[buffers[2 * j]],
                                        bufsize, &send_pos);
                            }
                            allflag[j] = (sendsizes[j] <=
                                          (numsend[j] + 1) * maxsize);
                            if(!allflag[j]) {
                                numtosend++;
                            }
                            numsend[j]++;
                            // add continuation signal
                            MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[2 * j]],
                                       bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[2 * j]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    0, &reqs[2 * j]);
                            numsent++;
                            numtoreceive++;
                        }
                    }
                }
            }

            if(numrecv < numtoreceive) {
                int index;
                int tag;
                int recv_pos;
                int send_pos;

                // wait for a message from any source and with any tag to
                // arrive
                MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
                numrecv++;
                // query the MPI_Status object to retrieve tag, source and
                // size
                tag = status.MPI_TAG;
                index = status.MPI_SOURCE;
                int nelements;
                MyMPI_Get_count(&status, MPI_PACKED, &nelements);
                // select an index in the buffers to use for receiving and
                // sending
                int freebuffer;
                if(index == MPIGlobal::size - 1) {
                    freebuffer = 2 * MPIGlobal::rank;
                } else {
                    freebuffer = 2 * index;
                }

                if(tag == 0) {
                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                               nelements, MPI_PACKED, index, 0, &status);
                    recv_pos = 0;
                    send_pos = 0;
                    while(recv_pos < nelements - 4) {
                        unsigned long id;
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &id, 1,
                                MPI_UNSIGNED_LONG);
                        Vec position;
                        unsigned int locid = 0;
                        char has_position = false;
                        if(_hashmap.count(id)) {
                            if(_hashmap[id].has_position()) {
                                position = _hashmap[id].get_position();
                                locid = _hashmap[id].get_locid();
                                has_position = true;
                            }
                        }
                        MyMPI_Pack(
                                &has_position, 1, MPI_BYTE,
                                &MPIGlobal::sendbuffer[buffers[freebuffer + 1]],
                                bufsize, &send_pos);
                        MyMPI_Pack(
                                &locid, 1, MPI_UNSIGNED,
                                &MPIGlobal::sendbuffer[buffers[freebuffer + 1]],
                                bufsize, &send_pos);
                        MyMPI_Pack(
                                &position[0], 2, MPI_DOUBLE,
                                &MPIGlobal::sendbuffer[buffers[freebuffer + 1]],
                                bufsize, &send_pos);
                    }
                    int flag;
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &flag, 1, MPI_INT);
                    if(!flag) {
                        numtoreceive++;
                    }

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[freebuffer + 1]],
                                send_pos, MPI_PACKED, index, 1,
                                &reqs[freebuffer + 1]);
                }

                if(tag == 1) {
                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer + 1]],
                               nelements, MPI_PACKED, index, 1, &status);
                    unsigned int j = numreceived[freebuffer / 2];
                    recv_pos = 0;
                    while(recv_pos < nelements) {
                        unsigned long id = incomplete[index][j];
                        char has_position;
                        unsigned int locid;
                        Vec position;
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer + 1]],
                                nelements, &recv_pos, &has_position, 1,
                                MPI_BYTE);
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer + 1]],
                                nelements, &recv_pos, &locid, 1, MPI_UNSIGNED);
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer + 1]],
                                nelements, &recv_pos, &position[0], 2,
                                MPI_DOUBLE);
                        if(has_position) {
                            Vec p12(position);
                            for(unsigned int i = 0; i < ndim_; i++) {
                                p12[i] = 1. +
                                         (p12[i] - box.get_anchor()[i]) /
                                                 box.get_sides()[i];
                            }
                            _hashmap[id].set_position(position, p12);
                            _hashmap[id].set_locid(locid);
                            _hashmap[id].set_rank(index);
                            for(int wall = 0; wall < 8; wall++) {
                                unsigned long gid = _maxid + id * 8 + wall;
                                if(_hashmap.count(gid)) {
                                    if(!_hashmap[gid].has_position()) {
                                        // have to apply wall
                                        double pos[2];
                                        pos[0] = position.x();
                                        pos[1] = position.y();
                                        Cuboid box = particles.get_container()
                                                             .get_cuboid();
                                        AdaptiveMeshUtils::
                                                get_periodic_position(
                                                        pos, -wall - 1, box);
#if ndim_ == 3
                                        Vec cpos(pos[0], pos[1], 0.);
#else
                                        Vec cpos(pos[0], pos[1]);
#endif
                                        Vec p12(cpos);
                                        for(unsigned int i = 0; i < ndim_;
                                            i++) {
                                            p12[i] = 1. +
                                                     (p12[i] -
                                                      box.get_anchor()[i]) /
                                                             box.get_sides()[i];
                                        }
                                        _hashmap[gid].set_position(cpos, p12);
                                        _hashmap[gid].set_locid(locid);
                                        _hashmap[gid].set_rank(index);
                                    }
                                }
                            }
                        } else {
                            orphans.push_back(id);
                        }
                        j++;
                    }
                    numreceived[freebuffer / 2] = j;
                }
            }
        }
        vector<MPI_Status> status((MPIGlobal::size - 1) * 2);
        MyMPI_Waitall((MPIGlobal::size - 1) * 2, &reqs[0], &status[0]);
        MyMPI_Barrier();

        // most of the particles now should have accurate positions
        // the only problem are particles for which we do not know on which
        // process they reside
        // we have to do a collective search for those
        unsigned int numorphans = orphans.size();
        unsigned int globnumorphans[MPIGlobal::size];
        MyMPI_Allgather(&numorphans, 1, MPI_UNSIGNED, globnumorphans, 1,
                        MPI_UNSIGNED);
        for(int irank = 0; irank < MPIGlobal::size; irank++) {
            if(globnumorphans[irank]) {
                std::vector<unsigned long> to_search(globnumorphans[irank]);
                if(irank == MPIGlobal::rank) {
                    to_search = orphans;
                }
                MyMPI_Bcast(&to_search[0], globnumorphans[irank],
                            MPI_UNSIGNED_LONG, irank);
                std::vector<char> has_position(globnumorphans[irank], false);
                std::vector<unsigned int> locid(globnumorphans[irank], false);
                std::vector<Vec> position(globnumorphans[irank]);
                for(unsigned int i = 0; i < to_search.size(); i++) {
                    if(_hashmap.count(to_search[i])) {
                        if(_hashmap[to_search[i]].has_position()) {
                            has_position[i] = true;
                            locid[i] = _hashmap[to_search[i]].get_locid();
                            position[i] = _hashmap[to_search[i]].get_position();
                        }
                    }
                }
                std::vector<char> position_flags;
                std::vector<unsigned int> locids;
                std::vector<Vec> positions;
                if(irank == MPIGlobal::rank) {
                    position_flags.resize(globnumorphans[irank] *
                                          MPIGlobal::size);
                    locids.resize(globnumorphans[irank] * MPIGlobal::size);
                    positions.resize(globnumorphans[irank] * MPIGlobal::size);
                }
                MyMPI_Gather(&has_position[0], globnumorphans[irank], MPI_BYTE,
                             &position_flags[0], globnumorphans[irank],
                             MPI_BYTE, irank);
                MyMPI_Gather(&locid[0], globnumorphans[irank], MPI_UNSIGNED,
                             &locids[0], globnumorphans[irank], MPI_UNSIGNED,
                             irank);
                MyMPI_Gather(&position[0], globnumorphans[irank] * 2,
                             MPI_DOUBLE, &positions[0],
                             globnumorphans[irank] * 2, MPI_DOUBLE, irank);
                if(irank == MPIGlobal::rank) {
                    for(int jrank = 0; jrank < MPIGlobal::size; jrank++) {
                        for(unsigned int i = 0; i < globnumorphans[irank];
                            i++) {
                            if(position_flags[jrank * globnumorphans[irank] +
                                              i]) {
                                Vec p12(positions
                                                [jrank * globnumorphans[irank] +
                                                 i]);
                                for(unsigned int i = 0; i < ndim_; i++) {
                                    p12[i] = 1. +
                                             (p12[i] - box.get_anchor()[i]) /
                                                     box.get_sides()[i];
                                }
                                _hashmap[orphans[i]].set_position(
                                        positions
                                                [jrank * globnumorphans[irank] +
                                                 i],
                                        p12);
                                _hashmap[orphans[i]].set_locid(
                                        locids[jrank * globnumorphans[irank] +
                                               i]);
                                _hashmap[orphans[i]].set_rank(jrank);
                                for(int wall = 0; wall < 8; wall++) {
                                    unsigned long gid =
                                            _maxid + orphans[i] * 8 + wall;
                                    if(_hashmap.count(gid)) {
                                        if(!_hashmap[gid].has_position()) {
                                            // have to apply wall
                                            double pos[2];
                                            pos[0] =
                                                    positions
                                                            [jrank * globnumorphans
                                                                             [irank] +
                                                             i].x();
                                            pos[1] =
                                                    positions
                                                            [jrank * globnumorphans
                                                                             [irank] +
                                                             i].y();
                                            Cuboid box =
                                                    particles.get_container()
                                                            .get_cuboid();
                                            AdaptiveMeshUtils::
                                                    get_periodic_position(
                                                            pos, -wall - 1,
                                                            box);
#if ndim_ == 3
                                            Vec cpos(pos[0], pos[1], 0.);
#else
                                            Vec cpos(pos[0], pos[1]);
#endif
                                            Vec p12(cpos);
                                            for(unsigned int i = 0; i < ndim_;
                                                i++) {
                                                p12[i] = 1. +
                                                         (p12[i] -
                                                          box.get_anchor()[i]) /
                                                                 box.get_sides()
                                                                         [i];
                                            }
                                            _hashmap[gid].set_position(cpos,
                                                                       p12);
                                            _hashmap[gid].set_locid(
                                                    locids[jrank * globnumorphans
                                                                           [irank] +
                                                           i]);
                                            _hashmap[gid].set_rank(jrank);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Get the position associated with the neighbour with the given id
     *
     * @param id Unique identifier of a neighbour
     * @return Position of the neighbour
     */
    Vec get_position(unsigned long id) {
        if(id >= _maxid) {
            int wall = get_wall(id);
            id = get_id(id);
            if(!_hashmap.count(id)) {
                cerr << "No position for id = " << id << "!" << endl;
                my_exit();
            }
            double pos[2];
            pos[0] = _hashmap[id].get_position().x();
            pos[1] = _hashmap[id].get_position().y();
            AdaptiveMeshUtils::get_periodic_position(pos, wall, _box);
#if ndim_ == 3
            Vec cpos(pos[0], pos[1], 0.);
#else
            Vec cpos(pos[0], pos[1]);
#endif
            return cpos;
        } else {
            if(!_hashmap.count(id)) {
                cerr << "No position for id = " << id << "!" << endl;
                my_exit();
            }
            return _hashmap[id].get_position();
        }
    }

    /**
     * @brief Get the original id corresponding to the given id
     *
     * @param id Unique ID of a particle or ghost copy
     * @return Original ID of the particle
     */
    unsigned long get_id(unsigned long id) {
        if(id < _maxid) {
            return id;
        }
        id = (id - _maxid) >> 3;
        return id;
    }

    /**
     * @brief Get the wall flag for the given particle ID
     *
     * @param id ID of a particle (can be a ghost copy)
     * @return Wall flag for the given ID
     */
    int get_wall(unsigned long id) {
        //        if(!_hashmap.count(id)){
        //            cerr << "No wall for id = " << id << "!" << endl;
        //            my_exit();
        //        }
        if(id < _maxid) {
            return 0;
        } else {
            unsigned long pid = (id - _maxid) >> 3;
            return _maxid + 8 * pid - 1 - id;
        }
    }

    /**
     * @brief Get the primitive variables for the particle with given ID
     *
     * @param id ID of a particle (can be a ghost copy)
     * @return Primitive variables of the (original) particle
     */
    StateVector get_W(unsigned long id) {
        if(id < _maxid) {
            return _hashmap[id].get_W();
        } else {
            unsigned long pid = (id - _maxid) >> 3;
            return _hashmap[pid].get_W();
        }
    }

    /**
     * @brief Get the start time of the timestep of the given particle
     *
     * @param id ID of a particle (can be a ghost copy)
     * @return Start time of the timestep of the (original) particle
     */
    unsigned long get_starttime(unsigned long id) {
        if(id < _maxid) {
            return _hashmap[id].get_starttime();
        } else {
            unsigned long pid = (id - _maxid) >> 3;
            return _hashmap[pid].get_starttime();
        }
    }

    /**
     * @brief Check if the given first id is a ghost copy of the given second id
     *
     * @param cid Possible ghost copy
     * @param id Original ID
     * @return True if the first id is equal to or a ghost copy of the second
     */
    bool is_copy(unsigned long cid, unsigned long id) {
        if(cid == id) {
            return true;
        }
        cid = (cid - _maxid) >> 3;
        return cid == id;
    }

    /**
     * @brief Check if the particle with the given id is local
     *
     * @param id Unique ID of a particle
     * @return True if the original particle is on the local MPI process
     */
    bool is_local(unsigned long id) {
        //        if(!_hashmap.count(id)){
        //            cerr << "No particle for id = " << id << "!" << endl;
        //            my_exit();
        //        }
        id = get_id(id);
        return _hashmap[id].get_rank() == MPIGlobal::rank;
    }

    /**
     * @brief Get the rank of the process that holds the original particle
     *
     * @param id Unique ID of a particle
     * @return Rank of the process that holds the original particle
     */
    int get_rank(unsigned long id) {
        //        if(!_hashmap.count(id)){
        //            cerr << "No particle for id = " << id << "!" << endl;
        //            my_exit();
        //        }
        id = get_id(id);
        return _hashmap[id].get_rank();
    }

    /**
     * @brief Get the index of the cell or particle on the original process
     *
     * @param id Unique ID of a particle or ghost copy
     * @return Local index of the cell or particle on its original process
     */
    unsigned int get_locid(unsigned long id) {
        id = get_id(id);
        if(!_hashmap.count(id)) {
            cerr << "No locid for id = " << id << "!" << endl;
            my_exit();
        }
        return _hashmap[id].get_locid();
    }

    /**
     * @brief Dump the list to a RestartFile
     *
     * @param rfile RestartFile to write to
     */
    void dump(RestartFile& rfile) {
        _box.dump(rfile);
        rfile.write(_maxid);
    }

    /**
     * @brief Restart constructor
     *
     * @param rfile RestartFile to read from
     */
    NewCellList(RestartFile& rfile) : _box(rfile) {
        rfile.read(_maxid);
    }
};

#endif
