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
 * @file TenetCellWeights.hpp
 *
 * @brief Weights associated with the hydrodynamical variables inside a
 * TenetCell.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TENETCELLWEIGHTS_HPP
#define TENETCELLWEIGHTS_HPP

#include "StateVector.hpp"
#include <vector>

/**
 * @brief Weights for a single TenetCell
 *
 * Stores the weights and provides utility functions.
 */
class TenetCellWeights {
  private:
    /*! @brief Weights associated with the cell */
    std::vector<StateVector> _weights;

    /*! @brief Indices of the order 1 basis function weights in the internal
     * vector */
    unsigned char _o1[ndim_];

  public:
    /**
     * @brief Empty constructor
     *
     * Should not be used, but is provided to make it possible to use
     * TenetCellWeights inside std::vector
     */
    inline TenetCellWeights() {
        _o1[0] = 0;
        _o1[1] = 0;
#if ndim_ == 3
        _o1[2] = 0;
#endif
    }

    /**
     * @brief Constructor
     *
     * @param size Number of weights to be stored
     * @param o1 Indices of the first order weights, used for slope limiting
     */
    inline TenetCellWeights(unsigned char size, unsigned char* o1) {
        _o1[0] = o1[0];
        _o1[1] = o1[1];
#if ndim_ == 3
        _o1[2] = o1[2];
#endif
        _weights.resize(size, 0.);
    }

    /**
     * @brief Get the first order weights in the given direction
     *
     * @param index Direction (x = 0, y = 1, ...)
     * @return Weights in the given direction
     */
    inline StateVector& get_first_order(unsigned char index) {
        return _weights[_o1[index]];
    }

    /**
     * @brief Get the zeroth order weights
     *
     * @return Zeroth order weights
     */
    inline StateVector& get_zeroth_order() {
        return _weights[0];
    }

    /**
     * @brief Slope limit the cell weights by using the given first order
     * limited weights
     *
     * This resets all but the zeroth order weights to zero, and then sets the
     * first order weights to the given weights. The calling function should
     * check that slope limiting is indeed necessary, as this method does not
     * check if the given first order weigths are equal to the old first order
     * weights.
     *
     * @param o1 First order slope limited weights
     */
    inline void slope_limit(StateVector* o1) {
        // do not reset first element, which is the zeroth order weight
        std::fill(_weights.begin() + 1, _weights.end(), 0.);
        _weights[_o1[0]] = o1[0];
        _weights[_o1[1]] = o1[1];
#if ndim_ == 3
        _weights[_o1[2]] = o1[2];
#endif
    }

    /**
     * @brief Access operator
     *
     * @param index Index of a basis function in the cell
     * @return Weights associated with that basis function
     */
    inline StateVector& operator[](unsigned char index) {
        return _weights[index];
    }

    /**
     * @brief Function that adds the given TenetCellWeights to the current
     * weights
     *
     * @param a TenetCellWeights to add
     * @return Reference to the current weights
     */
    inline TenetCellWeights& operator+=(TenetCellWeights a) {
        for(unsigned char i = 0; i < _weights.size(); ++i) {
            _weights[i] += a._weights[i];
        }
        return *this;
    }

    /**
     * @brief Multiplication operator
     *
     * @param s Scalar to multiply with
     * @return Reference to the TenetCellWeights with multiplied components
     */
    inline TenetCellWeights& operator*=(double s) {
        for(unsigned char i = 0; i < _weights.size(); ++i) {
            _weights[i] *= s;
        }
        return *this;
    }
};

#endif  // TENETCELLWEIGHTS_HPP
