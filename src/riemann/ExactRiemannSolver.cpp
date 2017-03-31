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
 * @file ExactRiemannSolver.cpp
 *
 * @brief Exact Riemann solver implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 * @author Robbert Verbeke (robbert.verbeke@ugent.be)
 */
#include "ExactRiemannSolver.hpp"
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include "ProgramLog.hpp"
#include "RestartFile.hpp"
#include "utilities/HelperFunctions.hpp"
#include <cmath>
#include <iostream>
using namespace std;

/**
 * @brief Constructor
 *
 * Initializes all quantities that are derived from \f$\gamma\f$.
 *
 * @param gamma Adiabatic index
 * @param tolerance Tolerance value used to decide when the pressure iteration
 * is converged
 * @param cutoff Cutoff to distinguish between a Newton-Raphson iteration and
 * Brent's method
 */
ExactRiemannSolver::ExactRiemannSolver(double gamma, double tolerance,
                                       double cutoff) {
    _gamma = gamma;
    // related quantities:
    _gp1d2g = 0.5 * (_gamma + 1.) / _gamma;  // gamma plus 1 divided by 2 gamma
    _gm1d2g = 0.5 * (_gamma - 1.) / _gamma;  // gamma minus 1 divided by 2 gamma
    _gm1dgp1 = (_gamma - 1.) /
               (_gamma + 1.);      // gamma minus 1 divided by gamma plus 1
    _tdgp1 = 2. / (_gamma + 1.);   // two divided by gamma plus 1
    _tdgm1 = 2. / (_gamma - 1.);   // two divided by gamma minus 1
    _gm1d2 = 0.5 * (_gamma - 1.);  // gamma minus 1 divided by 2
    _tgdgm1 = 2. * _gamma /
              (_gamma - 1.);  // two times gamma divided by gamma minus 1
    _ginv = 1. / _gamma;      // gamma inverse
    _tolerance = tolerance;
    _cutoff = cutoff;
    _counterBrent = 0;
    _counterTotal = 0;
    _vacuum = 0;

    LOGS("Initialized ExactRiemannSolver");
}

/**
 * @brief Destructor. Print status information to the stdout
 */
ExactRiemannSolver::~ExactRiemannSolver() {
    if(MPIGlobal::size > 1) {
        unsigned long counterBrent, counterTotal, vacuum;
        MyMPI_Reduce(&_counterBrent, &counterBrent, 1, MPI_LONG, MPI_SUM, 0);
        MyMPI_Reduce(&_counterTotal, &counterTotal, 1, MPI_LONG, MPI_SUM, 0);
        MyMPI_Reduce(&_vacuum, &vacuum, 1, MPI_LONG, MPI_SUM, 0);
        _counterBrent = counterBrent;
        _counterTotal = counterTotal;
        _vacuum = vacuum;
    }
    cout << "### Riemann solver statistics ###" << endl;
    cout << HelperFunctions::human_readable_counter(_counterBrent)
         << " Riemann solver evaluations using Brent's method" << endl;
    cout << HelperFunctions::human_readable_counter(_counterTotal)
         << " Riemann solver evaluations" << endl;
    cout << HelperFunctions::human_readable_counter(_vacuum)
         << " vacuum Riemann solver evaluations" << endl;
    cout << "Total time spent in Riemann solver: " << _timer.value()
         << " seconds" << endl;

    LOGS("ExactRiemannSolver destructed");
}

/**
 * @brief Solve the Riemann problem with the given left and right state
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param n Normal vector to the interface
 * @param mach Maximal mach-number of all shocks in the left state
 * @param reflective Flag indicating if the right state should be a
 * reflective copy of the left state
 * @return Solution of the Riemann problem
 */
StateVector ExactRiemannSolver::solve(StateVector& WL, StateVector& WR, Vec& n,
                                      double& mach, bool reflective) {
    _timer.start();
    _counterTotal++;

#if ndim_ == 3
    double vL = WL[1] * n[0] + WL[2] * n[1] + WL[3] * n[2];
    double vR = WR[1] * n[0] + WR[2] * n[1] + WR[3] * n[2];
#else
    double vL = WL[1] * n[0] + WL[2] * n[1];
    double vR = WR[1] * n[0] + WR[2] * n[1];
#endif

    if(reflective) {
        WR = WL;
        vR = -vL;
    }

    double aL = get_soundspeed(WL);
    double aR = get_soundspeed(WR);

    StateVector solution;
    if(!WL.rho() || !WR.rho()) {
        solution = solve_vacuum(WL, WR, vL, vR, aL, aR, n);
        _timer.stop();
        // a shockwave can never be adjacent to a region of vacuum (dixit Toro)
        // hence, the mach number for the shockwaves in this system is 0
        return solution;
    }

    if(2. * aL / (_gamma - 1.) + 2. * aR / (_gamma - 1.) <= vR - vL) {
        solution = solve_vacuum(WL, WR, vL, vR, aL, aR, n);
    } else {
        double p = 0.;
        double pguess = guess_p(WL, WR, vL, vR, aL, aR);
        double fp = f(p, WL, WR, vL, vR, aL, aR);
        double fpguess = f(pguess, WL, WR, vL, vR, aL, aR);
        if(fp < _cutoff) {
            while(fabs(p - pguess) > _tolerance * 0.5 * (p + pguess)) {
                p = pguess;
                pguess = pguess - fpguess / fprime(pguess, WL, WR, aL, aR);
                if(pguess < 0.) {
                    pguess = 0.;
                }
                fpguess = f(pguess, WL, WR, vL, vR, aL, aR);
            }
            p = pguess;
        } else {
            if(fp * fpguess >= 0.) {
                // Newton-Raphson until convergence or until usable interval is
                // found to use Brent's method
                while(fabs(p - pguess) > _tolerance * 0.5 * (p + pguess) &&
                      fpguess < 0.) {
                    p = pguess;
                    pguess = pguess - fpguess / fprime(pguess, WL, WR, aL, aR);
                    fpguess = f(pguess, WL, WR, vL, vR, aL, aR);
                }
            }
            // As soon as there is a usable interval: use Brent's method
            if(fabs(p - pguess) > _tolerance * 0.5 * (p + pguess) &&
               fpguess > 0.) {
                p = 0.;
                _counterBrent++;
                p = BrentsMethodSolve(p, pguess, _tolerance, WL, WR, vL, vR, aL,
                                      aR);
            } else {
                p = pguess;
            }
        }

        double u = 0.5 * (vL + vR) + 0.5 * (fb(p, WR, aR) - fb(p, WL, aL));

        // set mach number
        if(p > WL.p()) {
            // left shock
            mach = std::max(mach, sqrt(_gp1d2g * p / WL.p() + _gm1d2g));
        }
        if(p > WR.p()) {
            // right shock
            mach = std::max(mach, sqrt(_gp1d2g * p / WR.p() + _gm1d2g));
        }

        double vhalf;
        if(u < 0) {
            solution = WR;
            double pdpR = p / WR.p();
            if(p > WR.p()) {
                // shockwave
                double SR = vR + aR * sqrt(_gp1d2g * pdpR + _gm1d2g);
                if(SR > 0) {
                    solution.set_rho(WR.rho() * (pdpR + _gm1dgp1) /
                                     (_gm1dgp1 * pdpR + 1.));
                    solution.set_p(p);
                    vhalf = u - vR;
                } else {
                    // solution = WR
                    vhalf = 0.;
                }
            } else {
                // rarefaction wave
                double SHR = vR + aR;
                if(SHR > 0) {
                    double STR = u + aR * pow(pdpR, _gm1d2g);
                    if(STR <= 0) {
                        solution.set_rho(
                                WR.rho() *
                                pow(_tdgp1 - _gm1dgp1 / aR * vR, _tdgm1));
                        vhalf = _tdgp1 * (-aR + _gm1d2 * vR) - vR;
                        solution.set_p(WR.p() * pow(_tdgp1 - _gm1dgp1 / aR * vR,
                                                    _tgdgm1));
                    } else {
                        solution.set_rho(WR.rho() * pow(pdpR, _ginv));
                        solution.set_p(p);
                        vhalf = u - vR;
                    }
                } else {
                    // solution = WR
                    vhalf = 0.;
                }
            }
        } else {
            solution = WL;
            double pdpL = p / WL.p();
            if(p > WL.p()) {
                // shockwave
                double SL = vL - aL * sqrt(_gp1d2g * pdpL + _gm1d2g);
                if(SL < 0) {
                    solution.set_rho(WL.rho() * (pdpL + _gm1dgp1) /
                                     (_gm1dgp1 * pdpL + 1.));
                    solution.set_p(p);
                    vhalf = u - vL;
                } else {
                    // solution = WL
                    vhalf = 0.;
                }
            } else {
                // rarefaction wave
                double SHL = vL - aL;
                if(SHL < 0) {
                    double STL = u - aL * pow(pdpL, _gm1d2g);
                    if(STL > 0) {
                        solution.set_rho(
                                WL.rho() *
                                pow(_tdgp1 + _gm1dgp1 / aL * vL, _tdgm1));
                        vhalf = _tdgp1 * (aL + _gm1d2 * vL) - vL;
                        solution.set_p(WL.p() * pow(_tdgp1 + _gm1dgp1 / aL * vL,
                                                    _tgdgm1));
                    } else {
                        solution.set_rho(WL.rho() * pow(pdpL, _ginv));
                        vhalf = u - vL;
                        solution.set_p(p);
                    }
                } else {
                    // solution = WL
                    vhalf = 0.;
                }
            }
        }

        solution[1] += vhalf * n[0];
        solution[2] += vhalf * n[1];
#if ndim_ == 3
        solution[3] += vhalf * n[2];
#endif
    }

    _timer.stop();

    // sanity check
    bool check = (solution.rho() < 0.);
    check |= (solution.rho() != solution.rho());
    check |= (solution.vx() != solution.vx());
    check |= (solution.vy() != solution.vy());
#if ndim_ == 3
    check |= (solution.vz() != solution.vz());
#endif
    check |= (solution.p() < 0.);
    check |= (solution.p() != solution.p());
    if(check) {
        throw RiemannSolverException(WL, WR, solution);
    }
    return solution;
}

/**
 * @brief Solve the Riemann problem in the presence of vacuum
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param vL Left velocity in the frame of the interface
 * @param vR Right velocity in the frame of the interface
 * @param aL Soundspeed of the left state
 * @param aR Soundspeed of the right state
 * @param n Normal vector to the interface
 * @return Solution of the Riemann problem
 */
StateVector ExactRiemannSolver::solve_vacuum(StateVector& WL, StateVector& WR,
                                             double vL, double vR, double aL,
                                             double aR, Vec& n) {
    _vacuum++;
    StateVector solution;
    if(!WR.rho() && !WL.rho()) {
        return solution;
    }
    double vhalf;
    if(!WR.rho()) {
        solution = WL;
        // vacuum right state
        if(vL < aL) {
            double SL = vL + 2. * aL / (_gamma - 1.);
            if(SL > 0.) {
                solution.set_rho(
                        WL.rho() *
                        pow(2. / (_gamma + 1.) +
                                    (_gamma - 1.) / (_gamma + 1.) / aL * vL,
                            2. / (_gamma - 1.)));
                vhalf = 2. / (_gamma + 1.) * (aL + 0.5 * (_gamma - 1.) * vL) -
                        vL;
                solution.set_p(
                        WL.p() *
                        pow(2. / (_gamma + 1.) +
                                    (_gamma - 1.) / (_gamma + 1.) / aL * vL,
                            2. * _gamma / (_gamma - 1.)));
            } else {
                // vacuum
                solution.reset();
                vhalf = 0.;
            }
        } else {
            // solution = WL
            vhalf = 0.;
        }
    } else {
        if(!WL.rho()) {
            solution = WR;
            // vacuum left state
            if(-aR < vR) {
                double SR = vR - 2. * aR / (_gamma - 1.);
                if(SR >= 0.) {
                    // vacuum
                    solution.reset();
                    vhalf = 0.;
                } else {
                    solution.set_rho(
                            WR.rho() *
                            pow(2. / (_gamma + 1.) -
                                        (_gamma - 1.) / (_gamma + 1.) / aR * vR,
                                2. / (_gamma - 1.)));
                    vhalf = 2. / (_gamma + 1.) *
                                    (-aR + 0.5 * (_gamma - 1.) * vR) -
                            vR;
                    solution.set_p(
                            WR.p() *
                            pow(2. / (_gamma + 1.) -
                                        (_gamma - 1.) / (_gamma + 1.) / aR * vR,
                                2. * _gamma / (_gamma - 1.)));
                }
            } else {
                // solution = WR
                vhalf = 0.;
            }
        } else {
            // vacuum generation
            double SR = vR - 2. * aR / (_gamma - 1.);
            double SL = vL + 2. * aL / (_gamma - 1.);
            if(SR > 0. && SL < 0.) {
                solution.reset();
                vhalf = 0.;
            } else {
                if(SL >= 0.) {
                    solution = WL;
                    if(aL > vL) {
                        solution.set_rho(WL.rho() *
                                         pow(2. / (_gamma + 1.) +
                                                     (_gamma - 1.) /
                                                             (_gamma + 1.) /
                                                             aL * vL,
                                             2. / (_gamma - 1.)));
                        vhalf = 2. / (_gamma + 1.) *
                                        (aL + 0.5 * (_gamma - 1.) * vL) -
                                vL;
                        solution.set_p(WL.p() *
                                       pow(2. / (_gamma + 1.) +
                                                   (_gamma - 1.) /
                                                           (_gamma + 1.) / aL *
                                                           vL,
                                           2. * _gamma / (_gamma - 1.)));
                    } else {
                        // solution = WL
                        vhalf = 0.;
                    }
                } else {
                    solution = WR;
                    if(-aR < vR) {
                        solution.set_rho(WR.rho() *
                                         pow(2. / (_gamma + 1.) -
                                                     (_gamma - 1.) /
                                                             (_gamma + 1.) /
                                                             aR * vR,
                                             2. / (_gamma - 1.)));
                        vhalf = 2. / (_gamma + 1.) *
                                        (-aR + 0.5 * (_gamma - 1.) * vR) -
                                vR;
                        solution.set_p(WR.p() *
                                       pow(2. / (_gamma + 1.) -
                                                   (_gamma - 1.) /
                                                           (_gamma + 1.) / aR *
                                                           vR,
                                           2. * _gamma / (_gamma - 1.)));
                    } else {
                        // solution = WR
                        vhalf = 0.;
                    }
                }
            }
        }
    }

    solution[1] += vhalf * n[0];
    solution[2] += vhalf * n[1];
#if ndim_ == 3
    solution[3] += vhalf * n[2];
#endif

    // sanity check
    bool check = (solution.rho() < 0.);
    check |= (solution.rho() != solution.rho());
    check |= (solution.vx() != solution.vx());
    check |= (solution.vy() != solution.vy());
#if ndim_ == 3
    check |= (solution.vz() != solution.vz());
#endif
    check |= (solution.p() < 0.);
    check |= (solution.p() != solution.p());
    if(check) {
        throw RiemannSolverException(WL, WR, solution);
    }
    return solution;
}

/**
 * @brief Find the pressure in the intermediate state by using Brent's method
 *
 * @param lowerLimit Pressure value lower than intermediate pressure
 * @param upperLimit Pressure value higher than intermediate pressure
 * @param errorTol Tolerance to decide when the iteration is converged
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param vL Left velocity in the frame of the interface
 * @param vR Right velocity in the frame of the interface
 * @param aL Soundspeed of the left state
 * @param aR Soundspeed of the right state
 * @return Pressure in the intermediate state
 */
double ExactRiemannSolver::BrentsMethodSolve(double lowerLimit,
                                             double upperLimit, double errorTol,
                                             StateVector& WL, StateVector& WR,
                                             double vL, double vR, double aL,
                                             double aR) {
    double a = lowerLimit;
    double b = upperLimit;
    double c = 0.;
    double d = 1e230;

    double fa = f(a, WL, WR, vL, vR, aL, aR);
    double fb = f(b, WL, WR, vL, vR, aL, aR);

    double fc = 0.;
    double s = 0.;
    double fs = 0.;

    // if f(a) f(b) >= 0 then error-exit
    if(fa * fb >= 0) {
        // cout << "Error exit" << endl;
        if(fa < fb)
            return a;
        else
            return b;
    }

    // if |f(a)| < |f(b)| then swap (a,b) end if
    if(fabs(fa) < fabs(fb)) {
        double tmp = a;
        a = b;
        b = tmp;
        tmp = fa;
        fa = fb;
        fb = tmp;
    }

    c = a;
    fc = fa;
    bool mflag = true;

    while(!(fb == 0) && (fabs(a - b) > errorTol * (a + b) * 0.5)) {
        if((fa != fc) && (fb != fc))
            // Inverse quadratic interpolation
            s = a * fb * fc / (fa - fb) / (fa - fc) +
                b * fa * fc / (fb - fa) / (fb - fc) +
                c * fa * fb / (fc - fa) / (fc - fb);
        else
            // Secant Rule
            s = b - fb * (b - a) / (fb - fa);

        double tmp2 = (3. * a + b) / 4.;
        if(!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b))) ||
           (mflag && (fabs(s - b) >= (fabs(b - c) / 2.))) ||
           (!mflag && (fabs(s - b) >= (fabs(c - d) / 2.))) ||
           (mflag && (fabs(b - c) < errorTol * (b + c) * 0.5)) ||
           (!mflag && (fabs(c - d) < errorTol * (c + d) * 0.5))) {
            s = (a + b) / 2.;
            mflag = true;
        } else {
            mflag = false;
        }
        fs = f(s, WL, WR, vL, vR, aL, aR);
        d = c;
        c = b;
        fc = fb;
        if(fa * fs < 0.) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        // if |f(a)| < |f(b)| then swap (a,b) end if
        if(fabs(fa) < fabs(fb)) {
            double tmp = a;
            a = b;
            b = tmp;
            tmp = fa;
            fa = fb;
            fb = tmp;
        }
    }
    return b;
}

/**
 * @brief Riemann fL or fR function
 *
 * @param p (Temporary) pressure in the intermediate state
 * @param Wb Left or right StateVector
 * @param a Soundspeed of the left or right state
 * @return Riemann fL or fR function value
 */
double ExactRiemannSolver::fb(double p, StateVector& Wb, double a) {
    double fval = 0.;
    if(p > Wb.p()) {
        double A = _tdgp1 / Wb.rho();
        double B = _gm1dgp1 * Wb.p();
        fval = (p - Wb.p()) * sqrt(A / (p + B));
    } else {
        fval = _tdgm1 * a * (pow(p / Wb.p(), _gm1d2g) - 1.);
    }
    return fval;
}

/**
 * @brief Riemann gL or gR function
 *
 * @param p (Temporary) pressure in the intermediate state
 * @param Wb Left or right StateVector
 * @return Value of the Riemann gL or gR function
 */
double ExactRiemannSolver::gb(double p, StateVector& Wb) {
    double A = _tdgp1 / Wb.rho();
    double B = _gm1dgp1 * Wb.p();
    return sqrt(A / (p + B));
}

/**
 * @brief Derivative of Riemann fL or fR function
 *
 * @param p (Temporary) pressure in the intermediate state
 * @param Wb Left or right StateVector
 * @param a Soundpeed of the left or right state
 * @return Value of the derivative of the Riemann fL or fR function
 */
double ExactRiemannSolver::fprimeb(double p, StateVector& Wb, double a) {
    double fval = 0.;
    if(p > Wb.p()) {
        double A = _tdgp1 / Wb.rho();
        double B = _gm1dgp1 * Wb.p();
        fval = (1. - 0.5 * (p - Wb.p()) / (B + p)) * sqrt(A / (p + B));
    } else {
        fval = 1. / (Wb.rho() * a) * pow(p / Wb.p(), -_gp1d2g);
    }
    return fval;
}

/**
 * @brief Derivative of the Riemann f function
 *
 * @param p (Temporary) pressure in the intermediate state
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param aL Soundspeed of the left state
 * @param aR Soundspeed of the right state
 * @return Value of the derivative of the Riemann f function
 */
double ExactRiemannSolver::fprime(double p, StateVector& WL, StateVector& WR,
                                  double aL, double aR) {
    return fprimeb(p, WL, aL) + fprimeb(p, WR, aR);
}

/**
 * @brief Riemann f function
 *
 * @param p (Temporary) pressure in the intermediate state
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param vL Left velocity in the frame of the interface
 * @param vR Right velocity in the frame of the interface
 * @param aL Soundspeed of the left state
 * @param aR Soundspeed of the right state
 * @return Value of the Riemann f function
 */
double ExactRiemannSolver::f(double p, StateVector& WL, StateVector& WR,
                             double vL, double vR, double aL, double aR) {
    return fb(p, WL, aL) + fb(p, WR, aR) + (vR - vL);
}

/**
 * @brief Calculate a good starting estimate for the intermediate pressure
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param vL Left velocity in the frame of the interface
 * @param vR Right velocity in the frame of the interface
 * @param aL Soundspeed of the left state
 * @param aR Soundspeed of the right state
 * @return Starting estimate for the pressure iteration
 */
double ExactRiemannSolver::guess_p(StateVector& WL, StateVector& WR, double vL,
                                   double vR, double aL, double aR) {
    double pguess;
    double pmin = min(WL.p(), WR.p());
    double pmax = max(WL.p(), WR.p());
    double qmax = pmax / pmin;
    double ppv = 0.5 * (WL.p() + WR.p()) -
                 0.125 * (vR - vL) * (WL.rho() + WR.rho()) * (aL + aR);
    ppv = max(_tolerance * 0.5 * (WL.p() + WR.p()), ppv);
    if(qmax <= 2. && pmin <= ppv && ppv <= pmax) {
        pguess = ppv;
    } else {
        if(ppv < pmin) {
            // two rarefactions
            pguess = pow((aL + aR - _gm1d2 * (vR - vL)) /
                                 (aL / pow(WL.p(), _gm1d2g) +
                                  aR / pow(WR.p(), _gm1d2g)),
                         _tgdgm1);
        } else {
            // two shocks
            pguess = (gb(ppv, WL) * WL.p() + gb(ppv, WR) * WR.p() - vR + vL) /
                     (gb(ppv, WL) + gb(ppv, WR));
        }
    }
    // Toro: "Not that approximate solutions may predict, incorrectly, a
    // negative value for pressure (...). Thus in order to avoid negative guess
    // values we introduce the small positive constant _tolerance"
    pguess = max(_tolerance * 0.5 * (WL.p() + WR.p()), pguess);
    return pguess;
}

/**
 * @brief Get the soundspeed corresponding to the given StateVector
 *
 * @param W StateVector
 * @return Soundspeed corresponding to the StateVector
 */
double ExactRiemannSolver::get_soundspeed(const StateVector& W) {
    if(W.rho()) {
        return sqrt(_gamma * W.p() / W.rho());
    } else {
        return 0.;
    }
}

/**
 * @brief Calculate the energy for the given primitive variables
 *
 * @param W StateVector of primitive variables
 * @return Total energy
 */
double ExactRiemannSolver::get_energy(const StateVector& W) {
#if ndim_ == 3
    return 0.5 * W.rho() *
                   (W.vx() * W.vx() + W.vy() * W.vy() + W.vz() * W.vz()) +
           W.p() / (_gamma - 1.);
#else
    return 0.5 * W.rho() * (W.vx() * W.vx() + W.vy() * W.vy()) +
           W.p() / (_gamma - 1.);
#endif
}

/**
 * @brief Get the primitive variables corresponding to the given conserved
 * variables and the given volume
 *
 * @param volume Volume of a cell
 * @param Q Conserved variables in a cell
 * @param use_energy Flag indicating if the energy or entropy should be used to
 * convert total energy to pressure
 * @return Primitive variable StateVector
 */
StateVector ExactRiemannSolver::get_W(double volume, StateVector& Q,
                                      bool use_energy) {
    StateVector W;
    if(!Q.m()) {
        return W;
    }
    W.set_rho(Q.m() / volume);
    W.set_vx(Q.px() / Q.m());
    W.set_vy(Q.py() / Q.m());
#if ndim_ == 3
    W.set_vz(Q.pz() / Q.m());
    W.set_p((_gamma - 1.) *
            (Q.e() -
             0.5 * (Q.px() * Q.px() + Q.py() * Q.py() + Q.pz() * Q.pz()) /
                     Q.m()) /
            volume);
#else
    W.set_p((_gamma - 1.) *
            (Q.e() - 0.5 * (Q.px() * Q.px() + Q.py() * Q.py()) / Q.m()) /
            volume);
#endif
    for(unsigned int i = 0; i < NUM_PAQ; i++) {
        W.paq(i) = Q.paq(i) / Q.m();
    }
    if(use_energy) {
        // reset entropy
        Q.set_paq(Q.m() * W.p() / pow(W.rho(), _gamma));
        W.set_paq(Q.paq() / Q.m());
    } else {
        W.set_paq(Q.paq() / Q.m());
        W.set_p(W.paq() * pow(W.rho(), _gamma));
        Q.set_e(get_energy(W) * volume);
    }
    if(W.p() < 1.e-30) {
        W.set_p(1.e-30);
#if ndim_ == 3
        Q.set_e(0.5 * (Q.px() * Q.px() + Q.py() * Q.py() + Q.pz() * Q.pz()) /
                        Q.m() +
                W.p() * volume / (_gamma - 1.));
#else
        Q.set_e(0.5 * (Q.px() * Q.px() + Q.py() * Q.py()) / Q.m() +
                W.p() * volume / (_gamma - 1.));
#endif
    }

    if(W.rho() < 1.e-30) {
        W.reset();
    }

    return W;
}

/**
 * @brief Get the conserved variables corresponding to the given primitive
 * variables and volume
 *
 * @param volume Volume of a cell
 * @param W Primitive variables in a cell
 * @return Conserved variables StateVector
 */
StateVector ExactRiemannSolver::get_Q(double volume, const StateVector& W) {
    StateVector Q;
    Q.set_m(W.rho() * volume);
    Q.set_px(Q.m() * W.vx());
    Q.set_py(Q.m() * W.vy());
#if ndim_ == 3
    Q.set_pz(Q.m() * W.vz());
    Q.set_e((0.5 * W.rho() *
                     (W.vx() * W.vx() + W.vy() * W.vy() + W.vz() * W.vz()) +
             W.p() / (_gamma - 1.)) *
            volume);
#else
    Q.set_e((0.5 * W.rho() * (W.vx() * W.vx() + W.vy() * W.vy()) +
             W.p() / (_gamma - 1.)) *
            volume);
#endif
    // Ai = Pi/rhoi^gamma
    Q.set_paq(Q.m() * W.p() / pow(W.rho(), _gamma));
    for(unsigned int i = 0; i < NUM_PAQ; i++) {
        Q.paq(i) = Q.m() * W.paq(i);
    }

    // sanity check results
    if(Q.m() < 0. || Q.e() < 0.) {
        throw ConservedVariablesException(W, Q, volume);
    }

    return Q;
}

/**
 * @brief Get the given component of the flux matrix corresponding to the given
 * primitive variables at an interface with given velocity
 *
 * @param v Velocity of the interface
 * @param index unsigned integer index of the requested component
 * @param W Primitive variables at the interface
 * @return Flux through the interface
 */
StateVector ExactRiemannSolver::get_flux(const Vec& v, unsigned int index,
                                         const StateVector& W) {
    StateVector F;
    if(!W.rho()) {
        return F;
    }
    F[0] = W.rho() * (W[1 + index] - v[index]);
    F[1] = W.rho() * (W[1 + index] - v[index]) * W.vx();
    F[2] = W.rho() * (W[1 + index] - v[index]) * W.vy();
#if ndim_ == 3
    F[3] = W.rho() * (W[1 + index] - v[index]) * W.vz();
#endif
    F[1 + index] += W.p();
    F[ndim_ + 1] =
            (W[1 + index] - v[index]) * get_energy(W) + W.p() * W[1 + index];
    F.set_paq(W.rho() * (W[1 + index] - v[index]) * W.paq());
    F.set_dye(W.rho() * (W[1 + index] - v[index]) * W.dye());
    for(unsigned int i = 0; i < NUM_PAQ; i++) {
        F.paq(i) = W.rho() * (W[1 + index] - v[index]) * W.paq(i);
    }
    return F;
}

/**
 * @brief Get the number of Riemann solver evaluations
 *
 * @return Number of Riemann solver evaluations
 */
unsigned long ExactRiemannSolver::get_neval() {
    return _counterTotal;
}

/**
 * @brief Solve the Riemann problem and return fluxes instead of a half state
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param n Normal to the interface
 * @param v Interface velocity
 * @param reflective Flag indicating if the right state should be a
 * reflective copy of the left state
 * @return Flux along the interface normal
 */
StateVector ExactRiemannSolver::solve_for_flux(StateVector& WL, StateVector& WR,
                                               Vec& n, Vec& v,
                                               bool reflective) {
    double maxmach;
    StateVector Whalf = solve(WL, WR, n, maxmach, reflective);
    Whalf += v;

    StateVector fluxes[ndim_];
    fluxes[0] = get_flux(v, 0, Whalf);
    fluxes[1] = get_flux(v, 1, Whalf);
#if ndim_ == 3
    fluxes[2] = get_flux(v, 2, Whalf);
#endif

    StateVector flux;
#if ndim_ == 3
    flux = fluxes[0] * n[0] + fluxes[1] * n[1] + fluxes[2] * n[2];
#else
    flux = fluxes[0] * n[0] + fluxes[1] * n[1];
#endif

    return flux;
}

/**
 * @brief Convert gradients to time derivatives using the Euler equations in
 * primitive form
 *
 * These can be used in the MUSCL-Hancock prediction step.
 *
 * @param W Primitive variables
 * @param gradients Spatial derivatives of the primitive variables
 * @return Time derivatives of the primitive variables
 */
StateVector ExactRiemannSolver::get_time_derivative(
        const StateVector& W, const StateVector* gradients) {

    StateVector dWdt;

#if ndim_ == 3
    // there is a mistake in Toro (2009): the element rho in the first row of
    // matrices 3.91 and 3.92 should move 1 and 2 cells to the right
    // respectively
    dWdt[0] = -W.rho() * (gradients[0].vx() + gradients[1].vy() +
                          gradients[2].vz()) -
              W.vx() * gradients[0].rho() - W.vy() * gradients[1].rho() -
              W.vz() * gradients[2].rho();
    dWdt[1] = -W.vx() * gradients[0].vx() - W.vy() * gradients[1].vx() -
              W.vz() * gradients[2].vx() - gradients[0].p() / W.rho();
    dWdt[2] = -W.vx() * gradients[0].vy() - W.vy() * gradients[1].vy() -
              W.vz() * gradients[2].vy() - gradients[1].p() / W.rho();
    dWdt[3] = -W.vx() * gradients[0].vz() - W.vy() * gradients[1].vz() -
              W.vz() * gradients[2].vz() - gradients[2].p() / W.rho();
    dWdt[4] = -_gamma * W.p() * (gradients[0].vx() + gradients[1].vy() +
                                 gradients[2].vz()) -
              W.vx() * gradients[0].p() - W.vy() * gradients[1].p() -
              W.vz() * gradients[2].p();
#else
    dWdt[0] = -W.vx() * gradients[0].rho() - W.vy() * gradients[1].rho() -
              W.rho() * (gradients[0].vx() + gradients[1].vy());
    dWdt[1] = -W.vx() * gradients[0].vx() - W.vy() * gradients[1].vx() -
              gradients[0].p() / W.rho();
    dWdt[2] = -W.vx() * gradients[1].vy() - W.vy() * gradients[1].vy() -
              gradients[1].p() / W.rho();
    dWdt[3] = -_gamma * W.p() * (gradients[0].vx() + gradients[1].vy()) -
              W.vx() * gradients[0].p() - W.vy() * gradients[1].p();
#endif
    for(unsigned int i = 0; i < NUM_PAQ; i++) {
        dWdt.paq(i) =
                -W.vx() * gradients[0].paq(i) - W.vy() * gradients[1].paq(i);
#if ndim_ == 3
        dWdt.paq(i) -= W.vz() * gradients[2].paq(i);
#endif
    }

    return dWdt;
}

/**
 * @brief Dump the solver to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void ExactRiemannSolver::dump(RestartFile& rfile) {
    _timer.dump(rfile);

    rfile.write(_gamma);
    rfile.write(_gp1d2g);
    rfile.write(_gm1d2g);
    rfile.write(_gm1dgp1);
    rfile.write(_tdgp1);
    rfile.write(_tdgm1);
    rfile.write(_gm1d2);
    rfile.write(_tgdgm1);
    rfile.write(_ginv);
    rfile.write(_tolerance);
    rfile.write(_cutoff);

    rfile.write(_counterBrent);
    rfile.write(_counterTotal);
    rfile.write(_vacuum);

    LOGS("RiemannSolver dumped");
}

/**
 * @brief Restart constructor. Initialize the solver based on the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 */
ExactRiemannSolver::ExactRiemannSolver(RestartFile& rfile) : _timer(rfile) {
    rfile.read(_gamma);
    rfile.read(_gp1d2g);
    rfile.read(_gm1d2g);
    rfile.read(_gm1dgp1);
    rfile.read(_tdgp1);
    rfile.read(_tdgm1);
    rfile.read(_gm1d2);
    rfile.read(_tgdgm1);
    rfile.read(_ginv);
    rfile.read(_tolerance);
    rfile.read(_cutoff);

    rfile.read(_counterBrent);
    rfile.read(_counterTotal);
    rfile.read(_vacuum);

    LOGS("RiemannSolver restarted");
}
