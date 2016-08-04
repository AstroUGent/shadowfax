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
 * @file RungeKutta.hpp
 *
 * @brief Runge-Kutta integrator
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef RUNGEKUTTA_HPP
#define RUNGEKUTTA_HPP

#include <iostream>  // for operator<<, cerr, endl, etc

template <typename SystemState> class RungeKuttaFlux;

/**
 * @brief Runge-Kutta integrator
 */
template <class SystemState> class RungeKutta {
  private:
    /*! @brief Order of the integrator */
    unsigned char _order;

    /*! @brief RungeKuttaFlux used to calculate fluxes */
    RungeKuttaFlux<SystemState>& _flux;

    /*! @brief Current state of the system */
    SystemState& _state;

  public:
    /**
     * @brief Constructor
     *
     * @param order Order of the integrator
     * @param flux RungeKuttaFlux used to calculate fluxes
     * @param x0 Initial state of the system
     */
    RungeKutta(unsigned char order, RungeKuttaFlux<SystemState>& flux,
               SystemState& x0)
            : _flux(flux), _state(x0) {
        _order = order;
    }

    /**
     * @brief Integrate the system forward in time over the given time interval
     *
     * @param dt Time interval
     * @return System state at current time + dt
     */
    SystemState integrate(double dt) {
        SystemState xnp1 = _state;
        if(_order == 2) {
            SystemState temp = _state;
            SystemState k1 = _flux.get_flux(0., temp);
            temp = _state + dt * k1;
            SystemState k2 = _flux.get_flux(0., temp);
            xnp1 += 0.5 * dt * (k1 + k2);
            return xnp1;
        }
        if(_order == 3) {
            SystemState temp = _state;
            SystemState k1 = _flux.get_flux(0., temp);
            temp = _state + dt * k1;
            SystemState k2 = _flux.get_flux(0., temp);
            temp = _state + 0.25 * dt * (k1 + k2);
            SystemState k3 = _flux.get_flux(0., temp);
            xnp1 += dt * ((1. / 6.) * k1 + (1. / 6.) * k2 + (2. / 3.) * k3);
            return xnp1;
        }
        if(_order == 4) {
            SystemState temp = _state;
            SystemState k1 = _flux.get_flux(0., temp);
            temp = _state + 0.39175222700392 * dt * k1;
            SystemState k2 = _flux.get_flux(0., temp);
            temp = _state +
                   dt * (0.21766909633821 * k1 + 0.36841059262959 * k2);
            SystemState k3 = _flux.get_flux(0., temp);
            temp = _state +
                   dt * (0.08269208670950 * k1 + 0.13995850206999 * k2 +
                         0.25189177424738 * k3);
            SystemState k4 = _flux.get_flux(0., temp);
            temp = _state +
                   dt * (0.06796628370320 * k1 + 0.11503469844438 * k2 +
                         0.20703489864929 * k3 + 0.54497475021237 * k4);
            SystemState k5 = _flux.get_flux(0., temp);
            xnp1 += dt * (0.14681187618661 * k1 + 0.24848290924556 * k2 +
                          0.10425883036650 * k3 + 0.27443890091960 * k4 +
                          0.22600748319395 * k5);
            return xnp1;
        }
        std::cerr << "Order " << _order << " methods are not supported!"
                  << std::endl;
        return xnp1;
    }
};

#endif  // RUNGEKUTTA_HPP
