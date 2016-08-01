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
 * @file TenetMinModSlopeLimiter.hpp
 *
 * @brief TenetSlopeLimiter implementation containing a simple conserved
 * variables min-mod slope limiter
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TENETMINMODSLOPELIMITER_HPP
#define TENETMINMODSLOPELIMITER_HPP

#include "StateVector.hpp"
#include "TenetSlopeLimiter.hpp"
#include "TenetSystemState.hpp"

/**
 * @brief TenetSlopeLimiter implementation containing a simple conserved
 * variables min-mod slope limiter
 */
class TenetMinModSlopeLimiter : public TenetSlopeLimiter {
  private:
    /*! @brief Parameter setting the strength of the limiter */
    double _beta;

    /**
     * @brief Minmod slope limiter
     *
     * @param a StateVector a
     * @param b StateVector b
     * @param c StateVector c
     * @return minmod of a, b and c
     */
    static StateVector minmod(const StateVector& a, const StateVector& b,
                              const StateVector& c) {
        StateVector minmod;
        for(unsigned int i = 0; i < ndim_ + 2; i++) {
            if(a[i] < 0. && b[i] < 0. && c[i] < 0.) {
                minmod[i] = std::max(a[i], std::max(b[i], c[i]));
            }
            if(a[i] > 0. && b[i] > 0. && c[i] > 0.) {
                minmod[i] = std::min(a[i], std::min(b[i], c[i]));
            }
        }
        return minmod;
    }

    /**
     * @brief Check if the given statevectors are equal
     *
     * @param a StateVector a
     * @param b StateVector b
     * @return True if all components of a equal those of b
     */
    static bool statevector_equal(const StateVector& a, const StateVector& b) {
        bool equal = true;
        for(unsigned int i = 0; i < ndim_ + 2; i++) {
            equal &= (a[i] == b[i]);
        }
        return equal;
    }

  public:
    /**
     * @brief Constructor
     *
     * @param grid Reference to the TenetGrid on which the slope limiter acts
     * @param beta Parameter setting the strength of the limiter
     */
    TenetMinModSlopeLimiter(TenetGrid& grid, double beta = 0.5)
            : TenetSlopeLimiter(grid) {
        _beta = beta;
    }

    /**
     * @brief Limit the given TenetSystemState weights
     *
     * @param state TenetSystemState to limit
     */
    virtual void limit(TenetSystemState& state) {
#if ndim_ == 3
        for(auto it = _grid.begin(); it != _grid.end(); ++it) {
            unsigned int i = it.index();
            StateVector weightsL = state[it.index_left()].get_zeroth_order();
            StateVector weightsR = state[it.index_right()].get_zeroth_order();
            StateVector weightsT = state[it.index_above()].get_zeroth_order();
            StateVector weightsB = state[it.index_below()].get_zeroth_order();
            StateVector weightsF = state[it.index_front()].get_zeroth_order();
            StateVector weightsA = state[it.index_back()].get_zeroth_order();
            StateVector weights_O0 = state[i].get_zeroth_order();
            StateVector weights_O1[3];
            weights_O1[0] = state[i].get_first_order(0);
            weights_O1[1] = state[i].get_first_order(1);
            weights_O1[2] = state[i].get_first_order(2);
            StateVector new_weights[3];
            new_weights[0] = minmod(sqrt(3.) * weights_O1[0],
                                    _beta * (weights_O0 - weightsL),
                                    _beta * (weightsR - weights_O0)) /
                             sqrt(3.);
            new_weights[1] = minmod(sqrt(3.) * weights_O1[1],
                                    _beta * (weights_O0 - weightsA),
                                    _beta * (weightsF - weights_O0)) /
                             sqrt(3.);
            new_weights[2] = minmod(sqrt(3.) * weights_O1[2],
                                    _beta * (weights_O0 - weightsB),
                                    _beta * (weightsT - weights_O0)) /
                             sqrt(3.);
            if(!statevector_equal(new_weights[0], weights_O1[0]) ||
               !statevector_equal(new_weights[1], weights_O1[1]) ||
               !statevector_equal(new_weights[2], weights_O1[2])) {
                state[i].slope_limit(new_weights);
            }
        }
#else
        for(auto it = _grid.begin(); it != _grid.end(); ++it) {
            unsigned int i = it.index();
            StateVector weightsL = state[it.index_left()].get_zeroth_order();
            StateVector weightsR = state[it.index_right()].get_zeroth_order();
            StateVector weightsT = state[it.index_above()].get_zeroth_order();
            StateVector weightsB = state[it.index_below()].get_zeroth_order();
            StateVector weights_O0 = state[i].get_zeroth_order();
            StateVector weights_O1[2];
            weights_O1[0] = state[i].get_first_order(0);
            weights_O1[1] = state[i].get_first_order(1);
            StateVector new_weights[2];
            new_weights[0] = minmod(sqrt(3.) * weights_O1[0],
                                    _beta * (weights_O0 - weightsL),
                                    _beta * (weightsR - weights_O0)) /
                             sqrt(3.);
            new_weights[1] = minmod(sqrt(3.) * weights_O1[1],
                                    _beta * (weights_O0 - weightsB),
                                    _beta * (weightsT - weights_O0)) /
                             sqrt(3.);
            if(!statevector_equal(new_weights[0], weights_O1[0]) ||
               !statevector_equal(new_weights[1], weights_O1[1])) {
                state[i].slope_limit(new_weights);
            }
        }
#endif
    }
};

#endif  // TENETMINMODSLOPELIMITER_HPP
