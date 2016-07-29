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
 * @file RungeKuttaFlux.hpp
 *
 * @brief Interface for Runge-Kutta fluxes
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef RUNGEKUTTAFLUX_HPP
#define RUNGEKUTTAFLUX_HPP

/**
 * @brief Interface for Runge-Kutta fluxes
 *
 * For a differential equation of the form
 * \f[
 *   \frac{\text{d}}{\text{d}t} x(t) + R(x, t) = 0,
 * \f]
 * the Runge-Kutta flux is \f$ R(x,t) \f$.
 *
 * Implementations of this interface should define a method get_flux() that
 * returns the evaluated value of the flux, based on the given time and variable
 * value. The variable value is a template argument to this class.
 */
template <typename SystemState> class RungeKuttaFlux {
  public:
    /**
     * @brief Get the Runge-Kutta flux for the given time and system state
     *
     * @param t Time
     * @param x System state
     * @return Runge-Kutta flux
     */
    virtual SystemState get_flux(double t, SystemState& x) = 0;
};

#endif  // RUNGEKUTTAFLUX_HPP
