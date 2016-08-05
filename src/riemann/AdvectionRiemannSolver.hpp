/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file AdvectionRiemannSolver.hpp
 *
 * @brief Riemann solver that solves the advection equation instead of the
 * Euler equation: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ADVECTIONRIEMANNSOLVER_HPP
#define ADVECTIONRIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"    // for RiemannSolver
#include "StateVector.hpp"      // for StateVector
#include "utilities/Timer.hpp"  // for Timer
#include <string>               // for string

class RestartFile;

/**
 * @brief Riemann solver that solves the advection equation instead of the Euler
 * equations
 */
class AdvectionRiemannSolver : public RiemannSolver {
  private:
    /*! @brief Counts the number of times the solver was called */
    unsigned long _counter;

    /*! @brief Timer quantifying the time spent in the solver */
    Timer _timer;

    /*! @brief Adiabatic index of the fluid */
    double _gamma;

  public:
    AdvectionRiemannSolver(double gamma = 5. / 3.);
    ~AdvectionRiemannSolver();

    /**
     * @brief Get the identifying tag for this Solver
     *
     * @return "ADVE", because this is an AdvectionRiemannSolver
     */
    virtual std::string tag() {
        return "NOPR";
    }

    virtual StateVector solve(StateVector& WL, StateVector& WR, Vec& n,
                              double& mach0, bool reflective = false);
    virtual StateVector solve_for_flux(StateVector& WL, StateVector& WR, Vec& n,
                                       Vec& v, bool reflective = false);
    virtual double get_soundspeed(const StateVector& W);
    virtual void test();
    virtual StateVector get_Q(double volume, const StateVector& W);
    virtual StateVector get_W(double volume, StateVector& Q, bool use_energy);
    virtual StateVector get_flux(const Vec& v, unsigned int index,
                                 const StateVector& W);
    virtual StateVector get_time_derivative(const StateVector& W,
                                            const StateVector* gradients);
    virtual double get_gamma();

    unsigned long get_neval();

    virtual void dump(RestartFile& rfile);
    AdvectionRiemannSolver(RestartFile& rfile);
};

#endif  // ADVECTIONRIEMANNSOLVER_HPP
