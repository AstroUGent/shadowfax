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
 * @file ApproximateSolver.hpp
 *
 * @brief Two Rarefaction Riemann Solver (TRRS): header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef APPROXIMATESOLVER_HPP
#define APPROXIMATESOLVER_HPP

#include "RiemannSolver.hpp"

class RestartFile;

/**
 * @brief Two Rarefaction Riemann Solver
 *
 * Approximate Riemann solver that assumes that both the left and the right
 * intermediate state in the solution of the Riemann problem are rarefaction
 * waves. This allows us to skip the iterative procedure for finding p and
 * speeds up the solution of the Riemann problem.
 */
class TRRSSolver : public RiemannSolver {
  private:
    /*! \brief Adiabatic index \f$\gamma{}\f$ of the gas */
    double _gamma;
    /*! \brief \f$\frac{\gamma{}+1}{2\gamma{}}\f$ */
    double _gp1d2g;
    /*! \brief \f$\frac{\gamma{}-1}{2\gamma{}}\f$ */
    double _gm1d2g;
    /*! \brief \f$\frac{\gamma{}-1}{\gamma{}+1}\f$ */
    double _gm1dgp1;
    /*! \brief \f$\frac{2}{\gamma{}+1}\f$ */
    double _tdgp1;
    /*! \brief \f$\frac{2}{\gamma{}-1}\f$ */
    double _tdgm1;
    /*! \brief \f$\frac{\gamma{}-1}{2}\f$ */
    double _gm1d2;
    /*! \brief \f$\frac{2\gamma{}}{\gamma{}-1}\f$ */
    double _tgdgm1;
    /*! \brief \f$\frac{1}{\gamma}\f$ */
    double _ginv;

    /*! \brief Counts the number of times the solver was called */
    unsigned long _counter;

    /*! \brief Timer quantifying the time spent in the solver */
    Timer _timer;

    double get_energy(const StateVector& W);

  public:
    TRRSSolver(double gamma);
    ~TRRSSolver();

    /**
     * @brief Get the identifying tag for this Solver
     *
     * @return "TRRS", because this is a TRRSSolver
     */
    virtual std::string tag() { return "TRRS"; }

    virtual StateVector solve(StateVector& WL, StateVector& WR, Vec& n,
                              double& mach0, bool reflective = false);
    virtual StateVector solve_for_flux(StateVector& WL, StateVector& WR, Vec& n,
                                       Vec& v, bool reflective = false);
    virtual StateVector get_time_derivative(const StateVector& W,
                                            const StateVector* gradients);
    virtual double get_soundspeed(const StateVector& W);
    virtual void test();
    virtual StateVector get_Q(double volume, const StateVector& W);
    virtual StateVector get_W(double volume, StateVector& Q, bool use_energy);
    virtual StateVector get_flux(const Vec& v, unsigned int index,
                                 const StateVector& W);
    virtual double get_gamma();

    unsigned long get_neval();

    virtual void dump(RestartFile& rfile);
    TRRSSolver(RestartFile& rfile);
};

#endif
