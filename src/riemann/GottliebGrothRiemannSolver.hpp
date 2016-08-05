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
 * @file GottliebGrothRiemannSolver.hpp
 *
 * @brief Gottlieb-Groth RiemannSolver: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "RiemannSolver.hpp"  // for RiemannSolver
#include "StateVector.hpp"    // for StateVector

/**
 * @brief Solver that implements a Gottlieb-Groth generalized Riemann solver
 *
 * Does not work (yet).
 */
class GottliebGrothRiemannSolver : public RiemannSolver {
  private:
    /*! @brief Tolerance value */
    double _tolerance;

    /*! @brief Adiabatic index of the fluid */
    double _gamma;

    double fpstarL(double ustari, double uL, double gammaL, double aL,
                   double pL, double CL);
    double fpstarR(double ustari, double uR, double gammaR, double aR,
                   double pR, double CR);
    double fpstarLprime(double ustari, double uL, double gammaL, double aL,
                        double pL, double CL);
    double fpstarRprime(double ustari, double uR, double gammaR, double aR,
                        double pR, double CR);

    double get_energy(const StateVector& W);

  public:
    GottliebGrothRiemannSolver(double tolerance = 1.e-8,
                               double gamma = 5. / 3.);
    ~GottliebGrothRiemannSolver() {}

    virtual StateVector solve(StateVector& WL, StateVector& WR, Vec& n,
                              double& mach, bool reflective = false);

    double get_soundspeed(const StateVector& W);

    void test();

    StateVector get_Q(double volume, const StateVector& W);
    StateVector get_W(double volume, StateVector& Q, bool use_energy);
    StateVector get_flux(const Vec& v, unsigned int index,
                         const StateVector& W);

    /**
     * @brief Get the adiabatic index of the fluid
     *
     * @return Adiabatic index of the fluid
     */
    double get_gamma() {
        return _gamma;
    }
};
