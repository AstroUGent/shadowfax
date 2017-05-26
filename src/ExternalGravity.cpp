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
 * @file ExternalGravity.cpp
 *
 * @brief External gravity support.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "ExternalGravity.hpp"

#include <cmath>

/**
 * @brief Constructor.
 */
ExternalGravity::ExternalGravity() {
    _surface_density = 10.;
    _scale_height = 100.;
    _x_disc = 200.;
    _growth_time = 5.;
    _dynamical_time = 1.;
}

/**
 * @brief Get the acceleration at the given position, at the given time.
 *
 * @param position Position.
 * @return Acceleration.
 */
Vec ExternalGravity::get_acceleration(const Vec& position, double time) const {
    const double dx = position.x() - _x_disc;

    double reduction_factor = 1.;
    const double growth_time = _growth_time * _dynamical_time;
    if(time < growth_time) {
        reduction_factor = time / growth_time;
    }

    double x_accel = reduction_factor * 2. * M_PI * _surface_density *
                     std::tanh(std::abs(dx) / _scale_height);
    if(dx > 0.) {
        x_accel = -x_accel;
    }
    return Vec(x_accel, 0., 0.);
}
