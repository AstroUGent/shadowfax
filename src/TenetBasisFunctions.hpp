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
 * @file TenetBasisFunctions.hpp
 *
 * @brief Basis functions used for the Tenet representation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TENETBASISFUNCTIONS_HPP
#define TENETBASISFUNCTIONS_HPP

#include "Legendre.hpp"
#include "TenetCellWeights.hpp"
#include "Vec.hpp"
#include <vector>

/**
 * @brief Basis functions used for the Tenet representation
 */
class TenetBasisFunctions {
  private:
    /*! @brief Legendre polynomials used as basis functions */
    Legendre _legendre;

    /*! @brief Roots of the Legendre polynomials used for Gauss-Legendre
     * quadrature */
    std::vector<double> _zeros;

    /*! @brief Weights used for Gauss-Legendre quadrature */
    std::vector<double> _weights;

    /*! @brief Maximum order of the basis functions */
    unsigned char _order;

    /*! @brief Indices of the order 1 basis functions */
    unsigned char _o1[ndim_];

    /*! @brief Total number of basis functions */
    unsigned char _size;

  public:
    /**
     * @brief Constructor
     *
     * @param order Maximum order of the basis functions
     */
    TenetBasisFunctions(unsigned char order) : _legendre(order) {
        _zeros = _legendre.get_zeros();
        _weights = _legendre.get_weights();
        _order = order;
#if ndim_ == 3
        _size = (order + 1) * (order + 2) * (order + 3) / 6;
        _o1[0] = order * (order + 1) / 2 + order + 1;
        _o1[1] = order + 1;
        _o1[2] = 1;
#else
        _size = (order + 1) * (order + 2) / 2;
        _o1[0] = order + 1;
        _o1[1] = 1;
#endif
    }

    /**
     * @brief Get an empty TenetCellWeights instance with the correct size
     *
     * @return Empty TenetCellWeigths
     */
    TenetCellWeights get_empty_weights() {
        return TenetCellWeights(_size, _o1);
    }

    /**
     * @brief iterator to loop over basis functions
     */
    class iterator {
      private:
        /*! @brief Reference to the TenetBasisFunctions over which we iterate */
        TenetBasisFunctions& _basis;

        /*! @brief Index of the basis function the iterator is currently
         * pointing to */
        unsigned char _index;

        /*! @brief Order of x-direction basis function component */
        unsigned char _u;

        /*! @brief Order of y-direction basis function component */
        unsigned char _v;

#if ndim_ == 3
        /*! @brief Order of z-direction basis function component */
        unsigned char _w;
#endif

      public:
        /**
         * @brief Constructor
         *
         * @param basis Reference to the TenetBasisFunctions over which we
         * iterate
         * @param index Index of teh basis function the iterator is initially
         * pointing to
         */
        inline iterator(TenetBasisFunctions& basis, unsigned char index = 0)
                : _basis(basis) {
            _index = index;
            _u = 0;
            _v = 0;
#if ndim_ == 3
            _w = 0;
#endif
        }

        /**
         * @brief Increment operator
         *
         * @return Reference to the incremented iterator
         */
        inline iterator& operator++() {
            ++_index;
#if ndim_ == 3
            ++_w;
            if(_w > _basis._order - _u - _v) {
                _w = 0;
                ++_v;
                if(_v > _basis._order - _u) {
                    _v = 0;
                    ++_u;
                }
            }
#else
            ++_v;
            if(_v > _basis._order - _u) {
                _v = 0;
                ++_u;
            }
#endif
            return *this;
        }

        /**
         * @brief Evaluation of the current basis function
         *
         * @param x Coordinates for which the basis function is evaluated
         * @return Value of the basis function at the given coordinate
         */
        inline double operator()(Vec x) {
#if ndim_ == 3
            return _basis._legendre.polynomial1d_scaled(x[0], _u) *
                   _basis._legendre.polynomial1d_scaled(x[1], _v) *
                   _basis._legendre.polynomial1d_scaled(x[2], _w);
#else
            return _basis._legendre.polynomial1d_scaled(x[0], _u) *
                   _basis._legendre.polynomial1d_scaled(x[1], _v);
#endif
        }

#if ndim_ == 3
        /**
         * @brief Derivative of the current basis function in the x-direction
         *
         * @param x Coordinates for which the basis function is evaluated
         * @return Value of the derivative of the basis function at the given
         * coordinate
         */
        inline double prime_x(Vec x) {
            return _basis._legendre.polynomial1d_scaled_prime(x[0], _u) *
                   _basis._legendre.polynomial1d_scaled(x[1], _v) *
                   _basis._legendre.polynomial1d_scaled(x[2], _w);
        }

        /**
         * @brief Derivative of the current basis function in the y-direction
         *
         * @param x Coordinates for which the basis function is evaluated
         * @return Value of the derivative of the basis function at the given
         * coordinate
         */
        inline double prime_y(Vec x) {
            return _basis._legendre.polynomial1d_scaled(x[0], _u) *
                   _basis._legendre.polynomial1d_scaled_prime(x[1], _v) *
                   _basis._legendre.polynomial1d_scaled(x[2], _w);
        }

        /**
         * @brief Derivative of the current basis function in the z-direction
         *
         * @param x Coordinates for which the basis function is evaluated
         * @return Value of the derivative of the basis function at the given
         * coordinate
         */
        inline double prime_z(Vec x) {
            return _basis._legendre.polynomial1d_scaled(x[0], _u) *
                   _basis._legendre.polynomial1d_scaled(x[1], _v) *
                   _basis._legendre.polynomial1d_scaled_prime(x[2], _w);
        }
#else
        /**
         * @brief Derivative of the current basis function in the x-direction
         *
         * @param x Coordinates for which the basis function is evaluated
         * @return Value of the derivative of the basis function at the given
         * coordinate
         */
        inline double prime_x(Vec x) {
            return _basis._legendre.polynomial1d_scaled_prime(x[0], _u) *
                   _basis._legendre.polynomial1d_scaled(x[1], _v);
        }

        /**
         * @brief Derivative of the current basis function in the y-direction
         *
         * @param x Coordinates for which the basis function is evaluated
         * @return Value of the derivative of the basis function at the given
         * coordinate
         */
        inline double prime_y(Vec x) {
            return _basis._legendre.polynomial1d_scaled(x[0], _u) *
                   _basis._legendre.polynomial1d_scaled_prime(x[1], _v);
        }
#endif

        /**
         * @brief Index of the basis function to which the iterator is currently
         * pointing
         *
         * @return Current index
         */
        inline unsigned char index() {
            return _index;
        }

        /**
         * @brief Comparison operator
         *
         * @param it iterator to compare with
         * @return True if both iterators point to the same TenetBasisFunctions
         * and are at the same index
         */
        inline bool operator==(iterator it) {
            return (_index == it._index && &_basis == &it._basis);
        }

        /**
         * @brief Comparison operator
         *
         * @param it iterator to compare with
         * @return True if both iterators are at different indices or point to
         * different TenetBasisFunctions
         */
        inline bool operator!=(iterator it) {
            return !(*this == it);
        }
    };

    /**
     * @brief Get an iterator to the first basis function
     *
     * @return iterator to first basis function
     */
    inline iterator begin() {
        return iterator(*this);
    }

    /**
     * @brief Get an iterator to the end of the basis functions
     *
     * @return iterator to end of basis functions
     */
    inline iterator end() {
        return iterator(*this, _size);
    }

    /**
     * @brief Iterator used to carry out a volume quadrature over a cell
     */
    class cell_quadrature {
      private:
        /*! @brief Reference to the TenetBasisFunctions */
        TenetBasisFunctions& _basis;

        /*! @brief Index in the x-direction */
        unsigned char _ix;

        /*! @brief Index in the y-direction */
        unsigned char _iy;

#if ndim_ == 3
        /*! @brief Index in the z-direction */
        unsigned char _iz;
#endif

        /*! @brief Index of the current state of the iterator */
        unsigned char _index;

      public:
        /**
         * @brief Constructor
         *
         * @param basis Reference to the TenetBasisFunctions
         * @param index Initial index of the iterator
         */
        inline cell_quadrature(TenetBasisFunctions& basis,
                               unsigned char index = 0)
                : _basis(basis) {
            _index = index;
            _ix = 0;
            _iy = 0;
#if ndim_ == 3
            _iz = 0;
#endif
        }

        /**
         * @brief Get the current evaluation point of the iterator in scaled
         * cell coordinates
         *
         * @return Current evaluation point
         */
        inline Vec get_ksi() {
#if ndim_ == 3
            return Vec(_basis._zeros[_ix], _basis._zeros[_iy],
                       _basis._zeros[_iz]);
#else
            return Vec(_basis._zeros[_ix], _basis._zeros[_iy]);
#endif
        }

        /**
         * @brief Get the weight of the current evaluation point of the iterator
         *
         * @return Current weight
         */
        inline double get_weight() {
#if ndim_ == 3
            return _basis._weights[_ix] * _basis._weights[_iy] *
                   _basis._weights[_iz];
#else
            return _basis._weights[_ix] * _basis._weights[_iy];
#endif
        }

        /**
         * @brief Increment operator
         *
         * @return Reference to the incremented iterator
         */
        inline cell_quadrature& operator++() {
            ++_index;
#if ndim_ == 3
            ++_iz;
            if(_iz == _basis._zeros.size()) {
                _iz = 0;
                ++_iy;
                if(_iy == _basis._zeros.size()) {
                    _iy = 0;
                    ++_ix;
                }
            }
#else
            ++_iy;
            if(_iy == _basis._zeros.size()) {
                _iy = 0;
                ++_ix;
            }
#endif
            return *this;
        }

        /**
         * @brief Comparison operator
         *
         * @param it cell_quadrature to compare with
         * @return True if both cell_quadratures have the same index
         */
        inline bool operator==(cell_quadrature it) {
            return (_index == it._index && &_basis == &it._basis);
        }

        /**
         * @brief Comparison operator
         *
         * @param it cell_quadrature to compare with
         * @return True if both cell_quadratures have different indices
         */
        inline bool operator!=(cell_quadrature it) {
            return !(*this == it);
        }
    };

    /**
     * @brief Get a cell_quadrature pointing to the first evaluation point of
     * a volume cell quadrature
     *
     * @return cell_quadrature to begin of quadrature
     */
    inline cell_quadrature cell_quadrature_begin() {
        return cell_quadrature(*this);
    }

    /**
     * @brief Get a cell_quadrature pointing to the end of the volume cell
     * quadrature
     *
     * @return cell_quadrature to end of quadrature
     */
    inline cell_quadrature cell_quadrature_end() {
#if ndim_ == 3
        return cell_quadrature(*this,
                               _zeros.size() * _zeros.size() * _zeros.size());
#else
        return cell_quadrature(*this, _zeros.size() * _zeros.size());
#endif
    }

    /**
     * @brief Iterator used to perform face quadrature
     */
    class face_quadrature {
      private:
        /*! @brief Reference to the TenetBasisFunctions */
        TenetBasisFunctions& _basis;

        /*! @brief Reference to the midpoint of the face */
        Vec& _midpoint;

        /*! @brief Unit vectors in the face */
        Vec _face[ndim_ - 1];

        /*! @brief Index in the first direction inside the face */
        unsigned char _ix;
#if ndim_ == 3
        /*! @brief Index in the second direection inside the face */
        unsigned char _iy;
#endif
        /*! @brief Current position of the iterator */
        unsigned char _index;

      public:
        /**
         * @brief Constructor
         *
         * @param basis Reference to the TenetBasisFunctions
         * @param midpoint Reference to the midpoint of the face
         * @param normal Reference to the normal to the face
         * @param index Initial position of the iterator
         */
        face_quadrature(TenetBasisFunctions& basis, Vec& midpoint, Vec& normal,
                        unsigned char index = 0)
                : _basis(basis), _midpoint(midpoint) {
            _index = index;
            _ix = 0;
#if ndim_ == 3
            _iy = 0;
            if(normal.x()) {
                _face[0] = Vec(0., 1., 0.);
                _face[1] = Vec(0., 0., 1.);
            } else {
                if(normal.y()) {
                    _face[0] = Vec(1., 0., 0.);
                    _face[1] = Vec(0., 0., 1.);
                } else {
                    _face[0] = Vec(1., 0., 0.);
                    _face[1] = Vec(0., 1., 0.);
                }
            }
#else
            if(normal.x()) {
                _face[0] = Vec(0., 1.);
            } else {
                _face[0] = Vec(1., 0.);
            }
#endif
        }

        /**
         * @brief Get the left quadrature point
         *
         * @return Left quadrature point
         */
        inline Vec get_ksiL() {
#if ndim_ == 3
            return _face[0] * _basis._zeros[_ix] +
                   _face[1] * _basis._zeros[_iy] + _midpoint;
#else
            return _face[0] * _basis._zeros[_ix] + _midpoint;
#endif
        }

        /**
         * @brief Get the right quadrature point
         *
         * @return Right quadrature point
         */
        inline Vec get_ksiR() {
#if ndim_ == 3
            return _face[0] * _basis._zeros[_ix] +
                   _face[1] * _basis._zeros[_iy] - _midpoint;
#else
            return _face[0] * _basis._zeros[_ix] - _midpoint;
#endif
        }

        /**
         * @brief Get the weight of this point in the quadrature
         *
         * @return Weight
         */
        inline double get_weight() {
#if ndim_ == 3
            return _basis._weights[_ix] * _basis._weights[_iy];
#else
            return _basis._weights[_ix];
#endif
        }

        /**
         * @brief Addition operator
         *
         * @return Reference to the incremented face_quadrature
         */
        inline face_quadrature& operator++() {
            ++_index;
#if ndim_ == 3
            ++_iy;
            if(_iy == _basis._zeros.size()) {
                _iy = 0;
                ++_ix;
            }
#else
            ++_ix;
#endif
            return *this;
        }

        /**
         * @brief Comparison operator
         *
         * @param it face_quadrature to compare with
         * @return True if both face_quadratures are at the same position
         */
        inline bool operator==(face_quadrature it) {
            return (_index == it._index && &_midpoint == &it._midpoint &&
                    &_basis == &it._basis);
        }

        /**
         * @brief Comparison operator
         *
         * @param it face_quadrature to compare with
         * @return True if both face_quadratures are at different positions
         */
        inline bool operator!=(face_quadrature it) {
            return !(*this == it);
        }
    };

    /**
     * @brief Get a face_quadrature iterator to the beginning of a face
     *
     * @param midpoint Midpoint of the face
     * @param normal Face normal
     * @return face_quadrature pointing to beginning of face
     */
    inline face_quadrature face_quadrature_begin(Vec& midpoint, Vec& normal) {
        return face_quadrature(*this, midpoint, normal);
    }

    /**
     * @brief Get a face_quadrature iterator to the end of a face
     *
     * @param midpoint Midpoint of the face
     * @param normal Face normal
     * @return face_quadrature pointing to end of face
     */
    inline face_quadrature face_quadrature_end(Vec& midpoint, Vec& normal) {
#if ndim_ == 3
        return face_quadrature(*this, midpoint, normal,
                               _zeros.size() * _zeros.size());
#else
        return face_quadrature(*this, midpoint, normal, _zeros.size());
#endif
    }
};

#endif  // TENETBASISFUNCTIONS_HPP
