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
 * @file HLLCRiemannSolver.hpp
 *
 * @brief HLLC RiemannSolver: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HLLCRIEMANNSOLVER_HPP
#define HLLCRIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"  // for RiemannSolver
#include "StateVector.hpp"    // for StateVector
#include <string>             // for string

class RestartFile;
class Vec;

/**
 * @brief HLLC RiemannSolver
 */
class HLLCRiemannSolver : public RiemannSolver {
  private:
    /*! @brief Adiabatic index */
    double _gamma;

    /*! @brief Number of evaluations */
    unsigned long _neval;

    StateVector vacuum_flux(StateVector& WL, StateVector& WR, Vec& n, Vec& v,
                            double uL, double uR, double aL, double aR);

  public:
    HLLCRiemannSolver(double gamma);
    virtual ~HLLCRiemannSolver();

    /**
     * @brief Tag for RiemannSolver
     *
     * @return "HLLC"
     */
    virtual std::string tag() {
        return "HLLC";
    }

    virtual StateVector solve(StateVector& WL, StateVector& WR, Vec& n,
                              double& mach, bool reflective = false);

    StateVector solve_for_flux(StateVector& WL, StateVector& WR, Vec& n, Vec& v,
                               bool reflective = false);

    virtual StateVector get_time_derivative(const StateVector& W,
                                            const StateVector* gradients);

    virtual double get_soundspeed(const StateVector& W);

    void test();
  
    virtual StateVector get_Q(double volume, const StateVector& W);

    virtual StateVector get_W(double volume, StateVector& Q,
                              bool use_energy = true);

    virtual StateVector get_flux(const Vec& v, unsigned int index,
                                 const StateVector& W);

    virtual double get_gamma();

    virtual void dump(RestartFile& rfile);
    HLLCRiemannSolver(RestartFile& rfile);

    virtual unsigned long get_neval();
};

#endif  // HLLCRIEMANNSOLVER_HPP
