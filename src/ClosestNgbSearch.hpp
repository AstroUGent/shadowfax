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
 * @file ClosestNgbSearch.hpp
 *
 * @brief TreeWalker implementation for closest neighbour searching
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef CLOSESTNGBSEARCH_HPP
#define CLOSESTNGBSEARCH_HPP

#include "MPIMethods.hpp"
#include "Vec.hpp"  // for Vec
#include "utilities/StarParticle.hpp"
#include "utilities/TreeWalker.hpp"
#include <vector>  // for vector

class GasParticle;
class Leaf;
class PseudoNode;
class TreeNode;

/**
 * @brief TreeWalker implementation used to find the closest neighbour of a
 * given GasParticle
 */
class ClosestNgbSearch : public PeriodicTreeWalker {
  private:
    /*! @brief StarParticle for which the closest neighbour search is
     *  performed */
    StarParticle* _star;

    /*! @brief Center of the closest neighbour search */
    Vec _center;

    /*! @brief Radius squared of the closest neighbour search */
    double _radius;

    /*! @brief Current closest particle */
    GasParticle* _closest;

  public:
    /**
     * @brief Auxiliary class used to communicate information between MPI
     * processes during the closest neighbour search tree walk
     */
    class Export {
      private:
        /*! @brief StarParticle for which the closest neighbour search is
         *  performed */
        StarParticle* _star;

        /*! @brief Center of the closest neighbour search */
        Vec _center;

        /*! @brief Radius squared of the closest neighbour search */
        double _radius2;

      public:
        /**
         * @brief Constructor
         *
         * @param star StarParticle for which the closest neighbour search is
         * performed
         * @param center Center of the closest neighbour search
         * @param radius2 Radius squared of the closest neighbour search
         */
        Export(StarParticle* star, Vec center, double radius2) {
            _star = star;
            _center = center;
            _radius2 = radius2;
        }

        /**
         * @brief Add the relevant data to the given communication buffer
         *
         * @param buffer Buffer to write to
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_center[0], ndim_, MPI_DOUBLE, buffer, bufsize,
                       position);
            MyMPI_Pack(&_radius2, 1, MPI_DOUBLE, buffer, bufsize, position);
        }

        /**
         * @brief Read the response of the external communication from the given
         * communication buffer
         *
         * @param buffer Buffer to read from
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void unpack_data(void* buffer, int bufsize, int* position) {
            double radius2;
            int rank;
            GasParticle* closest;
            MyMPI_Unpack(buffer, bufsize, position, &radius2, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &rank, 1, MPI_INT);
            // a pointer should be an unsigned 64-bit integer on 64-bit systems
            MyMPI_Unpack(buffer, bufsize, position, &closest, 1,
                         MPI_UNSIGNED_LONG);
            if(radius2 < _star->get_closest_radius2()) {
                _star->set_closest_gasparticle(closest, radius2, rank);
            }
        }
    };

    /**
     * @brief Auxiliary class used to communicate particle information between
     * MPI processes during the gravity treewalk.
     */
    class Import {
      private:
        /*! @brief Center of the closest neighbour search */
        Vec _center;

        /*! @brief Radius squared of the closest neighbour search */
        double _radius2;

        /*! @brief Pointer to the closest GasParticle */
        GasParticle* _closest;

      public:
        /**
         * @brief Constructor. Initialize the import based on the given
         * communication stream
         *
         * @param buffer Buffer to read from
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        Import(void* buffer, int bufsize, int* position) {
            MyMPI_Unpack(buffer, bufsize, position, &_center[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_radius2, 1, MPI_DOUBLE);
        }

        /**
         * @brief Get the center of the closest neighbour search
         *
         * @return Center of the closest neighbour search
         */
        Vec get_center() {
            return _center;
        }

        /**
         * @brief Get the radius squared of the closest neighbour search
         *
         * @return Radius squared of the closest neighbour search
         */
        double get_radius2() {
            return _radius2;
        }

        /**
         * @brief Set the pointer to the closest GasParticle
         *
         * @param closest Pointer to the closest GasParticle
         */
        void set_closest(GasParticle* closest) {
            _closest = closest;
        }

        /**
         * @brief Set the most recent value of the search radius squared
         *
         * @param radius2 Most recent value of search radius squared
         */
        void set_radius2(double radius2) {
            _radius2 = radius2;
        }

        /**
         * @brief Add the relevant data to the given communication buffer to
         * send them back to the original process
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            int rank = MPIGlobal::rank;
            MyMPI_Pack(&_radius2, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&rank, 1, MPI_INT, buffer, bufsize, position);
            MyMPI_Pack(&_closest, 1, MPI_UNSIGNED_LONG, buffer, bufsize,
                       position);
        }
    };

    ClosestNgbSearch(StarParticle* star, Vec center, double radius);
    ClosestNgbSearch(Import& import);
    void set_position(Vec position);
    Vec get_position();

    GasParticle* get_closest();

    void increase_radius();

    bool splitnode(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);

    bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table);
    void periodicpseudonodeaction(PseudoNode* node, EwaldTable& ewald_table);
    void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table);

    void after_walk();

    void after_walk(Import& import);
};

#endif  // CLOSESTNGBSEARCH_HPP
