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
 * @file ExactRiemannSolver.hpp
 *
 * @brief Exact Riemann solver: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef EXACTRIEMANNSOLVER_HPP
#define EXACTRIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"
#include "StateVector.hpp"
#include "../src/utilities/Timer.hpp"

class RestartFile;

/**
 * @brief Solver implementation that represents an exact Riemann solver
 *
 * The solver finds the pressure in the intermediate state by performing a
 * mixed Brent/Newthon-Raphson iteration. This pressure is then used to
 * find the exact solution of the Riemann problem.
 */
class ExactRiemannSolver : public RiemannSolver{
private:
    /*! \brief Adiabatic index \f$\gamma\f$ of the fluid */
    double _gamma;

    /*! \brief \f$\frac{\gamma+1}{2\gamma}\f$ */
    double _gp1d2g;

    /*! \brief \f$\frac{\gamma-1}{2\gamma}\f$ */
    double _gm1d2g;

    /*! \brief \f$\frac{\gamma-1}{\gamma+1}\f$ */
    double _gm1dgp1;

    /*! \brief \f$\frac{2}{\gamma+1}\f$ */
    double _tdgp1;

    /*! \brief \f$\frac{2}{\gamma-1}\f$ */
    double _tdgm1;

    /*! \brief \f$\frac{\gamma-1}{2}\f$ */
    double _gm1d2;

    /*! \brief \f$\frac{2\gamma}{\gamma-1}\f$ */
    double _tgdgm1;

    /*! \brief \f$\frac{1}{\gamma}\f$ */
    double _ginv;

    /*! \brief Tolerance used to decide when the pressure iteration is
     *  converged */
    double _tolerance;

    /*! \brief Cutoff used to distinguish between a Newton-Raphson iteration
     *  and Brent's method */
    double _cutoff;

    /*! \brief Total number of Brent Riemann evaluations */
    unsigned long _counterBrent;

    /*! \brief Total number of Riemann solver evaluations */
    unsigned long _counterTotal;

    /*! \brief Total number of vacuum Riemann solver evaluations */
    unsigned long _vacuum;

    /*! \brief Timer for time spent in Riemann solver evaluations */
    Timer _timer;


    StateVector solve_vacuum(StateVector& WL, StateVector& WR, double vL,
                             double vR, double aL, double aR, Vec &n);
    double fb(double p, StateVector& Wb, double a);
    double gb(double p, StateVector& Wb);
    double fprimeb(double p, StateVector& Wb, double a);
    double fprime(double p, StateVector& WL, StateVector& WR, double aL,
                  double aR);
    double f(double p, StateVector& WL, StateVector& WR, double vL, double vR,
             double aL, double aR);
    double guess_p(StateVector& WL, StateVector& WR, double vL, double vR,
                   double aL, double aR);
    double BrentsMethodSolve(double lowerLimit, double upperLimit,
                             double errorTol, StateVector& WL, StateVector& WR,
                             double vL, double vR, double aL, double aR);

public:
    ExactRiemannSolver(double gamma = 1.66667, double tolerance = 1.e-8,
                  double cutoff = -5.);
    ~ExactRiemannSolver();

    /**
     * @brief Tag for ExactRiemannSolver
     *
     * @return "EXAC"
     */
    virtual std::string tag(){
        return "EXAC";
    }

    StateVector solve(StateVector& WL, StateVector& WR, Vec &n, double& mach,
                      bool reflective = false);

    /**
     * @brief Return the adiabatic index of the fluid
     *
     * @return Adiabatic index of the fluid
     */
    double get_gamma(){
        return _gamma;
    }

    double get_soundspeed(const StateVector& W);

    double get_energy(const StateVector& W);

    void test();

    StateVector get_Q(double volume, const StateVector &W);
    StateVector get_W(double volume, StateVector &Q, bool use_energy = true);
    StateVector get_flux(const Vec &v, unsigned int index,
                         const StateVector &W);

    unsigned long get_neval();

    void dump(RestartFile &rfile);
    ExactRiemannSolver(RestartFile &rfile);
};

#endif // EXACTRIEMANNSOLVER_HPP
