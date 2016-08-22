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
 * @file TenetGrid.hpp
 *
 * @brief Tenet grid: stores the cells and provides convenient access functions
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TENETGRID_HPP
#define TENETGRID_HPP

#include "TenetCell.hpp"         // for TenetCell
#include "Vec.hpp"               // for Vec, operator*, operator+
#include "utilities/Cuboid.hpp"  // for Cuboid
#include <cstddef>               // for NULL
#include <vector>                // for vector

/**
 * @brief Abstraction of the grid used by Tenet
 *
 * The TenetGrid contains the geometrical information about the grid, such as
 * cell sizes, cell origins, cell neighbour relations...
 * It also defines an iterator to facilitate grid access.
 */
class TenetGrid {
  private:
    /*! @brief Cells contained in the grid */
    std::vector<TenetCell*> _cells;

    /*! @brief Dimensions of the grid in number of cells in every direction */
    unsigned int _dimensions[ndim_];

    /*! @brief Physical width of a single (cubic) cell of the grid */
    double _cellwidth;

    /*! @brief Physical volume of a single cell of the grid */
    double _cellvolume;

  public:
    /**
     * @brief Constructor
     *
     * Allocates the cells and initializes cell centers.
     *
     * @param dimensions Dimensions of the grid in number of cells in every
     * direction
     * @param box Physical dimensions of the box containing the grid
     */
    inline TenetGrid(unsigned int* dimensions, Cuboid& box) {
        _dimensions[0] = dimensions[0];
        _dimensions[1] = dimensions[1];
        unsigned int size = _dimensions[0] * _dimensions[1];
#if ndim_ == 3
        _dimensions[2] = dimensions[2];
        size *= _dimensions[2];
#endif
        _cells.resize(size, NULL);
        Vec dx = box.get_sides();
        dx[0] /= dimensions[0];
        dx[1] /= dimensions[1];
#if ndim_ == 3
        dx[2] /= dimensions[2];
#endif
        for(unsigned int ix = 0; ix < _dimensions[0]; ++ix) {
            for(unsigned int iy = 0; iy < _dimensions[1]; ++iy) {
#if ndim_ == 2
                Vec i(ix + 0.5, iy + 0.5);
                Vec x = box.get_anchor() + i * dx;
                _cells[ix * _dimensions[0] + iy] = new TenetCell(x);
#else
                for(unsigned int iz = 0; iz < _dimensions[2]; ++iz) {
                    Vec i(ix + 0.5, iy + 0.5, iz + 0.5);
                    Vec x = box.get_anchor() + i * dx;
                    _cells[ix * _dimensions[0] * _dimensions[1] +
                           iy * _dimensions[1] + iz] = new TenetCell(x);
                }
#endif
            }
        }
        _cellwidth = dx[0];
        _cellvolume = dx[0] * dx[1];
#if ndim_ == 3
        _cellvolume *= dx[2];
#endif
    }

    /**
     * @brief Destructor
     *
     * Deallocates the cells.
     */
    inline ~TenetGrid() {
        for(unsigned int i = 0; i < _cells.size(); i++) {
            delete _cells[i];
        }
    }

    /**
     * @brief Get the physical width of a single (cubic) cell of the grid
     *
     * @return Width of a single cell
     */
    inline double get_cellwidth() {
        return _cellwidth;
    }

    /**
     * @brief Get the volume of a single cell of the grid
     *
     * @return Volume of a single cell
     */
    inline double get_cellvolume() {
        return _cellvolume;
    }

    /**
     * @brief Access operator
     *
     * @param i Index of a cell in the grid
     * @return Reference to the cell pointer
     */
    inline TenetCell*& operator[](unsigned int i) {
        return _cells[i];
    }

    /**
     * @brief Size of the grid in number of cells
     *
     * @return Number of cells in the grid
     */
    inline unsigned int size() {
#if ndim_ == 3
        return _dimensions[0] * _dimensions[1] * _dimensions[2];
#else
        return _dimensions[0] * _dimensions[1];
#endif
    }

    /**
     * @brief Iterator to loop over the cells of the grid
     *
     * The iterator also supports access to the neighbouring cells of the
     * current cell.
     */
    class iterator {
      private:
        /*! @brief Separate counters for the different directions, allowing
         * faster neigbour access */
        unsigned int _i[ndim_];

        /*! @brief Current index of the iterator in the TenetGrid */
        unsigned int _index;

        /*! @brief Reference to the grid over which we iterate */
        TenetGrid& _grid;

      public:
        /**
         * @brief Constructor.
         *
         * @param grid Reference to the TenetGrid over which we iterate
         * @param index Starting index of the iterator
         */
        inline iterator(TenetGrid& grid, unsigned int index = 0) : _grid(grid) {
            _i[0] = 0;
            _i[1] = 0;
#if ndim_ == 3
            _i[2] = 0;
#endif
            _index = index;
        }

        /**
         * @brief Dereference operator
         *
         * @return Reference to the cell the iterator is currently pointing to
         */
        inline TenetCell*& operator*() {
            return _grid._cells[_index];
        }

        /**
         * @brief Pointer operator
         *
         * @return Reference to the cell the iterator is currently pointing to
         */
        inline TenetCell*& operator->() {
            return _grid._cells[_index];
        }

        /**
         * @brief Increment operator
         *
         * We do not implement the iterator++, as this copies the entire
         * iterator.
         *
         * @return Reference to the incremented iterator
         */
        inline iterator& operator++() {
            ++_index;
#if ndim_ == 3
            ++_i[2];
            if(_i[2] == _grid._dimensions[2]) {
                _i[2] = 0;
                ++_i[1];
                if(_i[1] == _grid._dimensions[1]) {
                    _i[1] = 0;
                    ++_i[0];
                }
            }
#else
            ++_i[1];
            if(_i[1] == _grid._dimensions[1]) {
                _i[1] = 0;
                ++_i[0];
            }
#endif
            return *this;
        }

        /**
         * @brief Compare iterators
         *
         * @param it iterator to compare with
         * @return True if both iterators iterate over the same TenetGrid and
         * their internal states are equal, false otherwise
         */
        inline bool operator==(iterator it) {
            return (_index == it._index && &_grid == &it._grid);
        }

        /**
         * @brief Compare iterators
         *
         * @param it iterator to compare with
         * @return True if the iterators are not equal, as determined by
         * iterator::operator==, false otherwise
         */
        inline bool operator!=(iterator it) {
            return !(*this == it);
        }

        /**
         * @brief Get the index of the current cell the iterator is pointing to
         *
         * @return Index of the current cell
         */
        inline unsigned int index() {
            return _index;
        }

#if ndim_ == 3
        /**
         * @brief Get the index of the cell to the left of the cell the iterator
         * is currently pointing to
         *
         * @return Index of cell left of current cell
         */
        inline unsigned int index_left() {
            unsigned int ix = _i[0];
            if(ix) {
                --ix;
            } else {
                ix = _grid._dimensions[0] - 1;
            }
            return ix * _grid._dimensions[0] * _grid._dimensions[1] +
                   _i[1] * _grid._dimensions[1] + _i[2];
        }

        /**
         * @brief Get a reference to the cell to the left of the cell the
         * iterator is currently pointing to
         *
         * @return Cell to the left of current cell
         */
        inline TenetCell*& cell_left() {
            return _grid._cells[index_left()];
        }

        /**
         * @brief Get the index of the cell to the right of the cell the
         * iterator
         * is currently pointing to
         *
         * @return Index of cell right of current cell
         */
        inline unsigned int index_right() {
            unsigned int ix = _i[0];
            ++ix;
            if(ix == _grid._dimensions[0]) {
                ix = 0;
            }
            return ix * _grid._dimensions[0] * _grid._dimensions[1] +
                   _i[1] * _grid._dimensions[1] + _i[2];
        }

        /**
         * @brief Get a reference to the cell to the right of the cell the
         * iterator is currently pointing to
         *
         * @return Cell to the right of current cell
         */
        inline TenetCell*& cell_right() {
            return _grid._cells[index_right()];
        }

        /**
         * @brief Get the index of the cell in front of the cell the iterator is
         * currently pointing to
         *
         * @return Index of cell in front of current cell
         */
        inline unsigned int index_front() {
            unsigned int iy = _i[1];
            if(iy) {
                --iy;
            } else {
                iy = _grid._dimensions[1] - 1;
            }
            return _i[0] * _grid._dimensions[0] * _grid._dimensions[1] +
                   iy * _grid._dimensions[1] + _i[2];
        }

        /**
         * @brief Get a reference to the cell in front of the cell the iterator
         * is currently pointing to
         *
         * @return Cell in front of current cell
         */
        inline TenetCell*& cell_front() {
            return _grid._cells[index_front()];
        }

        /**
         * @brief Get the index of the cell at the back of the cell the iterator
         * is currently pointing to
         *
         * @return Index of cell at back of current cell
         */
        inline unsigned int index_back() {
            unsigned int iy = _i[1];
            ++iy;
            if(iy == _grid._dimensions[1]) {
                iy = 0;
            }
            return _i[0] * _grid._dimensions[0] * _grid._dimensions[1] +
                   iy * _grid._dimensions[1] + _i[2];
        }

        /**
         * @brief Get a reference to the cell at the back of the cell the
         * iterator is currently pointing to
         *
         * @return Cell tat back of current cell
         */
        inline TenetCell*& cell_back() {
            return _grid._cells[index_back()];
        }

        /**
         * @brief Get the index of the cell below the cell the iterator is
         * currently pointing to
         *
         * @return Index of cell below current cell
         */
        inline unsigned int index_below() {
            unsigned int iz = _i[2];
            if(iz) {
                --iz;
            } else {
                iz = _grid._dimensions[2] - 1;
            }
            return _i[0] * _grid._dimensions[0] * _grid._dimensions[1] +
                   _i[1] * _grid._dimensions[1] + iz;
        }

        /**
         * @brief Get a reference to the cell below the cell the iterator is
         * currently pointing to
         *
         * @return Cell below current cell
         */
        inline TenetCell*& cell_below() {
            return _grid._cells[index_below()];
        }

        /**
         * @brief Get the index of the cell above the cell the iterator is
         * currently pointing to
         *
         * @return Index of cell above current cell
         */
        inline unsigned int index_above() {
            unsigned int iz = _i[2];
            ++iz;
            if(iz == _grid._dimensions[2]) {
                iz = 0;
            }
            return _i[0] * _grid._dimensions[0] * _grid._dimensions[1] +
                   _i[1] * _grid._dimensions[1] + iz;
        }

        /**
         * @brief Get a reference to the cell above the cell the iterator is
         * currently pointing to
         *
         * @return Cell above current cell
         */
        inline TenetCell*& cell_above() {
            return _grid._cells[index_above()];
        }
#else
        /**
         * @brief Get the index of the cell to the left of the cell the iterator
         * is currently pointing to
         *
         * @return Index of cell left of current cell
         */
        inline unsigned int index_left() {
            unsigned int ix = _i[0];
            if(ix) {
                --ix;
            } else {
                ix = _grid._dimensions[0] - 1;
            }
            return ix * _grid._dimensions[0] + _i[1];
        }

        /**
         * @brief Get a reference to the cell to the left of the cell the
         * iterator is currently pointing to
         *
         * @return Cell to the left of current cell
         */
        inline TenetCell*& cell_left() {
            return _grid._cells[index_left()];
        }

        /**
         * @brief Get the index of the cell to the right of the cell the
         * iterator is currently pointing to
         *
         * @return Index of cell right of current cell
         */
        inline unsigned int index_right() {
            unsigned int ix = _i[0];
            ++ix;
            if(ix == _grid._dimensions[0]) {
                ix = 0;
            }
            return ix * _grid._dimensions[0] + _i[1];
        }

        /**
         * @brief Get a reference to the cell to the right of the cell the
         * iterator is currently pointing to
         *
         * @return Cell to the right of current cell
         */
        inline TenetCell*& cell_right() {
            return _grid._cells[index_right()];
        }

        /**
         * @brief Get the index of the cell above the cell the iterator is
         * currently pointing to
         *
         * @return Index of cell above current cell
         */
        inline unsigned int index_above() {
            unsigned int iy = _i[1];
            ++iy;
            if(iy == _grid._dimensions[1]) {
                iy = 0;
            }
            return _i[0] * _grid._dimensions[0] + iy;
        }

        /**
         * @brief Get a reference to the cell above the cell the iterator is
         * currently pointing to
         *
         * @return Cell above current cell
         */
        inline TenetCell*& cell_above() {
            return _grid._cells[index_above()];
        }

        /**
         * @brief Get the index of the cell below the cell the iterator is
         * currently pointing to
         *
         * @return Index of cell below current cell
         */
        inline unsigned int index_below() {
            unsigned int iy = _i[1];
            if(iy) {
                --iy;
            } else {
                iy = _grid._dimensions[1] - 1;
            }
            return _i[0] * _grid._dimensions[0] + iy;
        }

        /**
         * @brief Get a reference to the cell below the cell the iterator is
         * currently pointing to
         *
         * @return Cell below current cell
         */
        inline TenetCell*& cell_below() {
            return _grid._cells[index_below()];
        }
#endif

        /**
         * @brief Iterator used to loop over the faces of a single cell
         */
        class face_iterator {
          private:
            /*! @brief Index of the current face */
            unsigned char _index;

            /*! @brief Reference to an iterator pointing to the current cell */
            TenetGrid::iterator& _iterator;

          public:
            /**
             * @brief Constructor
             *
             * @param iterator Reference to an iterator pointing to the current
             * cell
             * @param index Index of the initial face
             */
            inline face_iterator(TenetGrid::iterator& iterator,
                                 unsigned char index = 0)
                    : _iterator(iterator) {
                _index = index;
            }

            /**
             * @brief Get the midpoint of the face in scaled coordinates
             *
             * @return Midpoint of the face
             */
            inline Vec get_midpoint() {
#if ndim_ == 3
                switch(_index) {
                    case 0:
                        // LEFT
                        return Vec(-1., 0., 0.);
                    case 1:
                        // RIGHT
                        return Vec(1., 0., 0.);
                    case 2:
                        // FRONT
                        return Vec(0., -1., 0.);
                    case 3:
                        // BACK
                        return Vec(0., 1., 0.);
                    case 4:
                        // BOTTOM
                        return Vec(0., 0., -1.);
                    case 5:
                        // TOP
                        return Vec(0., 0., 1.);
                }
#else
                switch(_index) {
                    case 0:
                        // BOTTOM
                        return Vec(0., -1.);
                    case 1:
                        // TOP
                        return Vec(0., 1.);
                    case 2:
                        // LEFT
                        return Vec(-1., 0.);
                    case 3:
                        // RIGHT
                        return Vec(1., 0.);
                }
#endif
                return Vec();
            }

            /**
             * @brief Index of the neighbouring cell that shares the face with
             * the current cell
             *
             * @return Index of neighbouring cell
             */
            inline unsigned int get_index_R() {
#if ndim_ == 3
                switch(_index) {
                    case 0:
                        // LEFT
                        return _iterator.index_left();
                    case 1:
                        // RIGHT
                        return _iterator.index_right();
                    case 2:
                        // FRONT
                        return _iterator.index_front();
                    case 3:
                        // BACK
                        return _iterator.index_back();
                    case 4:
                        // BOTTOM
                        return _iterator.index_below();
                    case 5:
                        // TOP
                        return _iterator.index_above();
                }
#else
                switch(_index) {
                    case 0:
                        // BOTTOM
                        return _iterator.index_below();
                    case 1:
                        // TOP
                        return _iterator.index_above();
                    case 2:
                        // LEFT
                        return _iterator.index_left();
                    case 3:
                        // RIGHT
                        return _iterator.index_right();
                }
#endif
                return 0;
            }

            /**
             * @brief Get the normal vector to the current face
             *
             * @return Normal vector
             */
            inline Vec get_normal() {
#if ndim_ == 3
                switch(_index) {
                    case 0:
                        // LEFT
                        return Vec(-1., 0., 0.);
                    case 1:
                        // RIGHT
                        return Vec(1., 0., 0.);
                    case 2:
                        // FRONT
                        return Vec(0., -1., 0.);
                    case 3:
                        // BACK
                        return Vec(0., 1., 0.);
                    case 4:
                        // BOTTOM
                        return Vec(0., 0., -1.);
                    case 5:
                        // TOP
                        return Vec(0., 0., 1.);
                }
#else
                switch(_index) {
                    case 0:
                        // BOTTOM
                        return Vec(0., -1.);
                    case 1:
                        // TOP
                        return Vec(0., 1.);
                    case 2:
                        // LEFT
                        return Vec(-1., 0.);
                    case 3:
                        // RIGHT
                        return Vec(1., 0.);
                }
#endif
                return Vec();
            }

            /**
             * @brief Increment operator
             *
             * We do not implement the iterator++, as this copies the entire
             * iterator.
             *
             * @return Reference to the incremented iterator
             */
            inline face_iterator& operator++() {
                ++_index;
                return *this;
            }

            /**
             * @brief Compare iterators
             *
             * @param it iterator to compare with
             * @return True if both iterators iterate over the same
             * TenetGrid::iterator and their internal states are equal, false
             * otherwise
             */
            inline bool operator==(face_iterator it) {
                return (_index == it._index && &_iterator == &it._iterator);
            }

            /**
             * @brief Compare iterators
             *
             * @param it iterator to compare with
             * @return True if the iterators are not equal, as determined by
             * iterator::operator==, false otherwise
             */
            inline bool operator!=(face_iterator it) {
                return !(*this == it);
            }
        };

        /**
         * @brief Get a face_iterator to the beginning of the current cell list
         * of faces
         *
         * @return face_iterator to first face
         */
        inline face_iterator face_begin() {
            return face_iterator(*this);
        }

        /**
         * @brief Get a face_iterator to the end of the current list of faces
         *
         * @return face_iterator to end of faces
         */
        inline face_iterator face_end() {
            return face_iterator(*this, 2 * ndim_);
        }
    };

    /**
     * @brief Get an iterator to the beginning of the TenetGrid
     *
     * @return Iterator pointing to first cell of the grid
     */
    inline iterator begin() {
        return iterator(*this);
    }

    /**
     * @brief Get an iterator to the end of the TenetGrid
     *
     * @return Iterator pointing to the end of the grid
     */
    inline iterator end() {
        return iterator(*this, size());
    }

    /**
     * @brief Iterator that loops over all faces of the grid once
     */
    class advanced_face_iterator {
      private:
        /*! @brief Reference to the TenetGrid over which we iterate */
        TenetGrid& _grid;

        /*! @brief Current position of the iterator */
        unsigned int _index;

        /*! @brief Index of the current face of the current cell */
        unsigned char _face_index;

        /*! @brief Separate counters for the different directions, allowing
         * faster neigbour access */
        unsigned int _i[ndim_];

      public:
        /**
         * @brief Constructor
         *
         * @param grid Reference to the TenetGrid over which we iterate
         * @param index Initial cell of the iterator
         * @param face_index Initial face of the initial cell
         */
        inline advanced_face_iterator(TenetGrid& grid, unsigned int index = 0,
                                      unsigned char face_index = 0)
                : _grid(grid) {
            _index = index;
            _face_index = face_index;
            _i[0] = 0;
            _i[1] = 0;
#if ndim_ == 3
            _i[2] = 0;
#endif
        }

        /**
         * @brief Increment operator
         *
         * @return Reference to the incremented advanced_face_iterator
         */
        inline advanced_face_iterator& operator++() {
            ++_face_index;
            if(_face_index == ndim_) {
                _face_index = 0;
                ++_index;
#if ndim_ == 3
                ++_i[2];
                if(_i[2] == _grid._dimensions[2]) {
                    _i[2] = 0;
                    ++_i[1];
                    if(_i[1] == _grid._dimensions[1]) {
                        _i[1] = 0;
                        ++_i[0];
                    }
                }
#else
                ++_i[1];
                if(_i[1] == _grid._dimensions[1]) {
                    _i[1] = 0;
                    ++_i[0];
                }
#endif
            }
            return *this;
        }

#if ndim_ == 3
        /**
         * @brief Get the index of the cell to the right of the cell the
         * iterator
         * is currently pointing to
         *
         * @return Index of cell right of current cell
         */
        inline unsigned int index_right() {
            unsigned int ix = _i[0];
            ++ix;
            if(ix == _grid._dimensions[0]) {
                ix = 0;
            }
            return ix * _grid._dimensions[0] * _grid._dimensions[1] +
                   _i[1] * _grid._dimensions[1] + _i[2];
        }

        /**
         * @brief Get the index of the cell at the back of the cell the iterator
         * is currently pointing to
         *
         * @return Index of cell at back of current cell
         */
        inline unsigned int index_back() {
            unsigned int iy = _i[1];
            ++iy;
            if(iy == _grid._dimensions[1]) {
                iy = 0;
            }
            return _i[0] * _grid._dimensions[0] * _grid._dimensions[1] +
                   iy * _grid._dimensions[1] + _i[2];
        }

        /**
         * @brief Get the index of the cell above the cell the iterator is
         * currently pointing to
         *
         * @return Index of cell above current cell
         */
        inline unsigned int index_above() {
            unsigned int iz = _i[2];
            ++iz;
            if(iz == _grid._dimensions[2]) {
                iz = 0;
            }
            return _i[0] * _grid._dimensions[0] * _grid._dimensions[1] +
                   _i[1] * _grid._dimensions[1] + iz;
        }
#else
        /**
         * @brief Get the index of the cell to the right of the cell the
         * iterator is currently pointing to
         *
         * @return Index of cell right of current cell
         */
        inline unsigned int index_right() {
            unsigned int ix = _i[0];
            ++ix;
            if(ix == _grid._dimensions[0]) {
                ix = 0;
            }
            return ix * _grid._dimensions[0] + _i[1];
        }

        /**
         * @brief Get the index of the cell above the cell the iterator is
         * currently pointing to
         *
         * @return Index of cell above current cell
         */
        inline unsigned int index_above() {
            unsigned int iy = _i[1];
            ++iy;
            if(iy == _grid._dimensions[1]) {
                iy = 0;
            }
            return _i[0] * _grid._dimensions[0] + iy;
        }
#endif

        /**
         * @brief Get the index of the cell on the left of the face
         *
         * @return Index of left cell
         */
        inline unsigned int get_index_L() {
            return _index;
        }

        /**
         * @brief Get the index of the cell on the right of the face
         *
         * @return Index of right cell
         */
        inline unsigned int get_index_R() {
#if ndim_ == 3
            switch(_face_index) {
                case 0:
                    // RIGHT
                    return index_right();
                case 1:
                    // BACK
                    return index_back();
                case 2:
                    // TOP
                    return index_above();
            }
#else
            switch(_face_index) {
                case 0:
                    // RIGHT
                    return index_right();
                case 1:
                    // TOP
                    return index_above();
            }
#endif
            return 0;
        }

        /**
         * @brief Get the normal vector to the current face
         *
         * @return Normal vector
         */
        inline Vec get_normal() {
#if ndim_ == 3
            switch(_face_index) {
                case 0:
                    // RIGHT
                    return Vec(1., 0., 0.);
                case 1:
                    // BACK
                    return Vec(0., 1., 0.);
                case 2:
                    // TOP
                    return Vec(0., 0., 1.);
            }
#else
            switch(_face_index) {
                case 0:
                    // RIGHT
                    return Vec(1., 0.);
                case 1:
                    // TOP
                    return Vec(0., 1.);
            }
#endif
            return Vec();
        }

        /**
         * @brief Get the midpoint of the face in scaled coordinates
         *
         * @return Midpoint of the face
         */
        inline Vec get_midpoint() {
#if ndim_ == 3
            switch(_face_index) {
                case 0:
                    // RIGHT
                    return Vec(1., 0., 0.);
                case 1:
                    // BACK
                    return Vec(0., 1., 0.);
                case 2:
                    // TOP
                    return Vec(0., 0., 1.);
            }
#else
            switch(_face_index) {
                case 0:
                    // RIGHT
                    return Vec(1., 0.);
                case 1:
                    // TOP
                    return Vec(0., 1.);
            }
#endif
            return Vec();
        }

        /**
         * @brief Comparison operator
         *
         * @param it advanced_face_iterator to compare with
         * @return True if both advanced_face_iterators are at the same position
         */
        inline bool operator==(advanced_face_iterator it) {
            return (_index == it._index && &_grid == &it._grid);
        }

        /**
         * @brief Comparison operator
         *
         * @param it advanced_face_iterator to compare with
         * @return True if both advanced_face_iterators are at different
         * positions
         */
        inline bool operator!=(advanced_face_iterator it) {
            return !(*this == it);
        }
    };

    /**
     * @brief Get an advanced_face_iterator to the first face of the grid
     *
     * @return advanced_face_iterator to beginning of grid
     */
    inline advanced_face_iterator faces_begin() {
        return advanced_face_iterator(*this);
    }

    /**
     * @brief Get an advanced_face_iterator to the end of the grid
     *
     * @return advanced_face_iterator to end of grid
     */
    inline advanced_face_iterator faces_end() {
        return advanced_face_iterator(*this, size());
    }
};

#endif  // TENETGRID_HPP
