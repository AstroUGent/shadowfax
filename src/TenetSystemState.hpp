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
 * @file TenetSystemState.hpp
 *
 * @brief State of the Tenet grid, used by the Runge-Kutta integration scheme
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TENETSYSTEMSTATE_HPP
#define TENETSYSTEMSTATE_HPP

#include "TenetCellWeights.hpp"
#include <vector>

/**
 * @brief System state for the Tenet grid
 */
class TenetSystemState {
  private:
    /*! @brief Internal vector of weights for every cell in the grid */
    std::vector<TenetCellWeights> _weights;

  public:
    /**
     * @brief Constructor
     *
     * @param size Size of the grid
     */
    inline TenetSystemState(unsigned int size) {
        _weights.resize(size);
    }

    /**
     * @brief Size of the internal weight vector
     *
     * @return Number of cells in the grid
     */
    inline unsigned int size() {
        return _weights.size();
    }

    /**
     * @brief Access operator
     *
     * @param index Index of a cell in the grid
     * @return Weights of the requested cell
     */
    inline TenetCellWeights& operator[](unsigned int index) {
        return _weights[index];
    }

    /**
     * @brief Addition operator
     *
     * @param a TenetSystemState to add to this TenetSystemState
     * @return Reference to the incremented TenetSystemState
     */
    inline TenetSystemState& operator+=(TenetSystemState a) {
        for(unsigned int i = 0; i < _weights.size(); i++) {
            _weights[i] += a._weights[i];
        }
        return *this;
    }

    /**
     * @brief Multiplication operator
     *
     * @param s Scalar to multiply with the components of this TenetSystemState
     * @return Reference to the resulting TenetSystemState
     */
    inline TenetSystemState& operator*=(double s) {
        for(unsigned int i = 0; i < _weights.size(); ++i) {
            _weights[i] *= s;
        }
        return *this;
    }
};

/**
 * @brief Free multiplication operator for TenetSystemState and scalar
 *
 * @param a TenetSystemState
 * @param s Scalar
 * @return TenetSystemState multiplied by scalar
 */
inline TenetSystemState operator*(TenetSystemState a, double s) {
    return a *= s;
}

/**
 * @brief Free multiplication operator for scalar and TenetSystemState
 *
 * @param s Scalar
 * @param a TenetSystemState
 * @return TenetSystemState multiplied by scalar
 */
inline TenetSystemState operator*(double s, TenetSystemState a) {
    return a *= s;
}

/**
 * @brief Free addition operator for two TenetSystemState instances
 *
 * @param a TenetSystemState A
 * @param b TenetSystemState B
 * @return Sum of A and B
 */
inline TenetSystemState operator+(TenetSystemState a, TenetSystemState b) {
    return a += b;
}

#endif  // TENETSYSTEMSTATE_HPP
