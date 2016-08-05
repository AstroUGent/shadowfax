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
 * @file CollelaGlazRiemannSolver.hpp
 *
 * @brief Collela-Glaz RiemannSolver: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "RiemannSolver.hpp"  // for RiemannSolver
#include "StateVector.hpp"    // for StateVector

/**
 * @brief Generalized Solver for advanced equations of state
 *
 * Not working (yet).
 */
class CollelaGlazRiemannSolver : public RiemannSolver {
  private:
    /*! @brief Tolerance value */
    double _tolerance;

    double dp_du(const StateVector& W);
    double dp_drho(const StateVector& W);

    double calculate_p_star(const StateVector& WL, const StateVector& WR);
    double get_wavespeed(double pstar, const StateVector& W, double gammaS,
                         double gammaO, double GammaS, double GammaO);
    double get_energy(const StateVector& W);

  public:
    CollelaGlazRiemannSolver(double tolerance = 1.e-8);
    ~CollelaGlazRiemannSolver() {}

    virtual StateVector solve(StateVector& WL, StateVector& WR, Vec& n,
                              double& mach, bool reflective = false);

    double get_soundspeed(const StateVector& W);

    void test();

    StateVector get_Q(double volume, const StateVector& W);
    StateVector get_W(double volume, StateVector& Q, bool use_energy);
    StateVector get_flux(const Vec& v, unsigned int index,
                         const StateVector& W);

    /**
     * @brief Dummy adiabatic index
     *
     * @return 0, because a constant adiabatic index is not defined for this
     * type of Solver
     */
    double get_gamma() {
        return 0.;
    }
};
