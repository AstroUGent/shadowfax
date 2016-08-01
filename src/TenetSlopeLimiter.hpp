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
 * @file TenetSlopeLimiter.hpp
 *
 * @brief General interface for Tenet slope limiters
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TENETSLOPELIMITER_HPP
#define TENETSLOPELIMITER_HPP

class TenetGrid;
class TenetSystemState;

/**
 * @brief General interface for Tenet slope limiters
 */
class TenetSlopeLimiter {
  protected:
    /*! @brief Reference to the TenetGrid on which the slope limiter acts */
    TenetGrid& _grid;

  public:
    /**
     * @brief Constructor
     *
     * @param grid Reference to the TenetGrid on which the slope limiter acts
     */
    TenetSlopeLimiter(TenetGrid& grid) : _grid(grid) {}

    /**
     * @brief Limit the weights of the given TenetSystemState
     *
     * @param state TenetSystemState on which to act
     */
    virtual void limit(TenetSystemState& state) = 0;
};

#endif  // TENETSLOPELIMITER_HPP
