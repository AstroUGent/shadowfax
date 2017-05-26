/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ExternalGravity.hpp
 *
 * @brief External gravity support.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#ifndef HEAD_EXTERNALGRAVITY
#define HEAD_EXTERNALGRAVITY

#include "Vec.hpp"

/**
 * @brief External gravity support.
 */
class ExternalGravity {
  private:
    /*! @brief Surface density of the disc .*/
    double _surface_density;

    /*! @brief Scale height of the disc. */
    double _scale_height;

    /*! @brief Position of the disc along the x-axis. */
    double _x_disc;

    /*! @brief Dynamical time of the system. */
    double _dynamical_time;

    /*! @brief Time over which to grow the disk in units of the dynamical time.
     */
    double _growth_time;

  public:
    ExternalGravity();

    Vec get_acceleration(const Vec& position, double time) const;
};

#endif  // HEAD_EXTERNALGRAVITY
