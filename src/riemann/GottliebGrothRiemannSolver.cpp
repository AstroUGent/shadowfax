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
 * @file GottliebGrothRiemannSolver.cpp
 *
 * @brief Gottlieb-Groth solver for the Riemann problem
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "GottliebGrothRiemannSolver.hpp"
#include "Error.hpp"  // for my_exit
#include "Vec.hpp"    // for Vec
#include <cmath>      // for pow, sqrt, fabs
#include <iostream>   // for operator<<, basic_ostream, etc
using namespace std;

/**
 * @brief Constructor
 *
 * @param tolerance Tolerance value used internally
 * @param gamma Adiabatic index of the fluid
 */
GottliebGrothRiemannSolver::GottliebGrothRiemannSolver(double tolerance,
                                                       double gamma) {
    _tolerance = tolerance;
    _gamma = gamma;
}

/**
 * @brief Solve the Riemann problem with given left and right state
 *
 * @warning Not implemented (yet)!
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param n Interface normal
 * @param mach Maximal mach-number of all shocks in the left state
 * @param reflective Flag indicating if the interface is reflective
 * @return  Solution of the Riemann problem
 */
StateVector GottliebGrothRiemannSolver::solve(StateVector& WL, StateVector& WR,
                                              Vec& n, double& mach,
                                              bool reflective) {
    StateVector Wstar;

    double aL = get_soundspeed(WL);
    double aR = get_soundspeed(WR);
    double gammaL = _gamma;
    double gammaR = _gamma;
    double CL = gammaL * WL.p() / aL;
    double CR = gammaR * WR.p() / aR;

    double uLtilda = WL.vx() + 2. / (gammaL - 1.) * aL;
    double uRtilda = WR.vx() - 2. / (gammaR - 1.) * aR;
    double sigma = gammaL;
    if(WL.p() < WR.p()) {
        sigma = gammaR;
    }
    double z = (gammaL - 1.) / (gammaR - 1.) * aR / aL *
               pow(WL.p() / WR.p(), 0.25 * (sigma - 1.) / sigma);
    double ustar0 = (uLtilda * z + uRtilda) / (1. + z);

    double ustari = ustar0;

    double pstarL = fpstarL(ustari, WL.vx(), gammaL, aL, WL.p(), CL);
    double pstarR = fpstarR(ustari, WR.vx(), gammaR, aR, WR.p(), CR);
    while(fabs(1. - pstarL / pstarR) > _tolerance) {
        double pstarLprime =
                fpstarLprime(ustari, WL.vx(), gammaL, aL, WL.p(), CL);
        double pstarRprime =
                fpstarRprime(ustari, WR.vx(), gammaR, aR, WR.p(), CR);
        ustari -= (pstarL - pstarR) / (pstarLprime - pstarRprime);
        pstarL = fpstarL(ustari, WL.vx(), gammaL, aL, WL.p(), CL);
        pstarR = fpstarR(ustari, WR.vx(), gammaR, aR, WR.p(), CR);
    }

    double u = ustari;
    double p = pstarL;

    if(u < 0) {
        Wstar = WR;
        if(p > WR.p()) {
            double SR = WR.vx() +
                        sqrt(_gamma * WR.p() / WR.rho()) *
                                sqrt(0.5 * (_gamma + 1.) / _gamma * p / WR.p() +
                                     0.5 * (_gamma - 1.) / _gamma);
            if(SR > 0) {
                Wstar.set_rho(
                        WR.rho() *
                        (p / WR.p() + (_gamma - 1.) / (_gamma + 1.)) /
                        ((_gamma - 1.) / (_gamma + 1.) * p / WR.p() + 1.));
                Wstar.set_p(p);
                Wstar.set_vx(u);
            }
            // else : solution = WR
        } else {
            double SHR = WR.vx() + aR;
            if(SHR > 0) {
                double STR =
                        u + aR * pow(p / WR.p(), 0.5 * (_gamma - 1.) / _gamma);
                if(STR <= 0) {
                    Wstar.set_rho(WR.rho() *
                                  pow(2. / (_gamma + 1.) -
                                              (_gamma - 1.) / (_gamma + 1.) /
                                                      aR * WR.vx(),
                                      2. / (_gamma - 1.)));
                    Wstar.set_vx(2. / (_gamma + 1.) *
                                 (-aR + 0.5 * (_gamma - 1.) * WR.vx()));
                    Wstar.set_p(WR.p() *
                                pow(2. / (_gamma + 1.) -
                                            (_gamma - 1.) / (_gamma + 1.) / aR *
                                                    WR.vx(),
                                    2. * _gamma / (_gamma - 1.)));
                } else {
                    Wstar.set_rho(WR.rho() * pow(p / WR.p(), 1. / _gamma));
                    Wstar.set_p(p);
                    Wstar.set_vx(u);
                }
            }
            // else : solution = WR
        }
    } else {
        Wstar = WL;
        if(p > WL.p()) {
            double SL = WL.vx() -
                        aL * sqrt(0.5 * (_gamma + 1.) / _gamma * p / WL.p() +
                                  0.5 * (_gamma - 1.) / _gamma);
            if(SL < 0) {
                Wstar.set_rho(
                        WL.rho() *
                        (p / WL.p() + (_gamma - 1.) / (_gamma + 1.)) /
                        ((_gamma - 1.) / (_gamma + 1.) * p / WL.p() + 1.));
                Wstar.set_p(p);
                Wstar.set_vx(u);
            }
            // else : solution = WL
        } else {
            double SHL = WL.vx() - aL;
            if(SHL < 0) {
                double STL =
                        u - aL * pow(p / WL.p(), 0.5 * (_gamma - 1.) / _gamma);
                if(STL > 0) {
                    Wstar.set_rho(WL.rho() *
                                  pow(2. / (_gamma + 1.) +
                                              (_gamma - 1.) / (_gamma + 1.) /
                                                      aL * WL.vx(),
                                      2. / (_gamma - 1.)));
                    Wstar.set_vx(2. / (_gamma + 1.) *
                                 (aL + 0.5 * (_gamma - 1.) * WL.vx()));
                    Wstar.set_p(WL.p() *
                                pow(2. / (_gamma + 1.) +
                                            (_gamma - 1.) / (_gamma + 1.) / aL *
                                                    WL.vx(),
                                    2. * _gamma / (_gamma - 1.)));
                } else {
                    Wstar.set_rho(WL.rho() * pow(p / WL.p(), 1. / _gamma));
                    Wstar.set_vx(u);
                    Wstar.set_p(p);
                }
            }
            // else : solution = WL
        }
    }

    return Wstar;
}

/**
 * @brief Auxiliary function, no idea what it does
 *
 * @warning Probably does not work!
 *
 * @param ustari Velocity in the intermediate state
 * @param uL Velocity in the left state
 * @param gammaL Adiabatic index of the left state
 * @param aL Soundspeed of the left state
 * @param pL Pressure of the left state
 * @param CL Shock speed in the left state
 * @return Value of the function
 */
double GottliebGrothRiemannSolver::fpstarL(double ustari, double uL,
                                           double gammaL, double aL, double pL,
                                           double CL) {
    if(ustari > uL) {
        double astarL = aL - 0.5 * (gammaL - 1.) * (ustari - uL);
        return pL * pow(astarL / aL, 2. * gammaL / (gammaL - 1.));
    } else {
        double WL = 0.25 * (gammaL + 1.) * (ustari - uL) / aL;
        WL -= sqrt(1. + WL * WL);
        return pL + CL * (ustari - uL) * WL;
    }
}

/**
 * @brief Derivative of GottliebGrothRiemannSolver::fpstarL()
 *
 * @warning Probably does not work!
 *
 * @param ustari Velocity in the intermediate state
 * @param uL Velocity in the left state
 * @param gammaL Adiabatic index of the left state
 * @param aL Soundspeed of the left state
 * @param pL Pressure of the left state
 * @param CL Shock speed in the left state
 * @return Value of the derivative
 */
double GottliebGrothRiemannSolver::fpstarLprime(double ustari, double uL,
                                                double gammaL, double aL,
                                                double pL, double CL) {
    if(ustari > uL) {
        double astarL = aL - 0.5 * (gammaL - 1.) * (ustari - uL);
        return -gammaL * pL * pow(astarL / aL, 2. * gammaL / (gammaL - 1.)) /
               astarL;
    } else {
        double WL = 0.25 * (gammaL + 1.) * (ustari - uL) / aL;
        WL -= sqrt(1. + WL * WL);
        double WL2 = WL * WL;
        return 2. * CL * WL * WL2 / (1. + WL2);
    }
}

/**
 * @brief Auxiliary function, no idea what it does
 *
 * @warning Probably does not work!
 *
 * @param ustari Velocity in the intermediate state
 * @param uR Velocity in the right state
 * @param gammaR Adiabatic index of the right state
 * @param aR Soundspeed of the right state
 * @param pR Pressure of the right state
 * @param CR Shock speed in the right state
 * @return Value of the function
 */
double GottliebGrothRiemannSolver::fpstarR(double ustari, double uR,
                                           double gammaR, double aR, double pR,
                                           double CR) {
    if(ustari < uR) {
        double astarR = aR + 0.5 * (gammaR - 1.) * (ustari - uR);
        return pR * pow(astarR / aR, 2. * gammaR / (gammaR - 1.));
    } else {
        double WR = 0.25 * (gammaR + 1.) * (ustari - uR) / aR;
        WR += sqrt(1. + WR * WR);
        return pR + CR * (ustari - uR) * WR;
    }
}

/**
 * @brief Derivative of GottliebGrothRiemannSolver::fpstarR()
 *
 * @warning Probably does not work!
 *
 * @param ustari Velocity in the intermediate state
 * @param uR Velocity in the right state
 * @param gammaR Adiabatic index of the right state
 * @param aR Soundspeed of the right state
 * @param pR Pressure of the right state
 * @param CR Shock speed in the right state
 * @return Value of the derivative
 */
double GottliebGrothRiemannSolver::fpstarRprime(double ustari, double uR,
                                                double gammaR, double aR,
                                                double pR, double CR) {
    if(ustari < uR) {
        double astarR = aR + 0.5 * (gammaR - 1.) * (ustari - uR);
        return gammaR * pR * pow(astarR / aR, 2. * gammaR / (gammaR - 1.)) /
               astarR;
    } else {
        double WR = 0.25 * (gammaR + 1.) * (ustari - uR) / aR;
        WR += sqrt(1. + WR * WR);
        double WR2 = WR * WR;
        return 2. * CR * WR * WR2 / (1. + WR2);
    }
}

/**
 * @brief Get the soundspeed corresponding to the given primitive StateVector
 *
 * @param W StateVector containing primitive variables
 * @return Soundspeed
 */
double GottliebGrothRiemannSolver::get_soundspeed(const StateVector& W) {
    return sqrt(_gamma * W.p() / W.rho());
}

/**
 * @brief Method to test the solver on a set of problems with known solution
 */
void GottliebGrothRiemannSolver::test() {
    cout << "Testing the Gottlieb-Groth Riemann solver" << endl;
    double rhoL[6] = {1., 1., 1., 1., 5.99924, 1.};
    double rhoR[6] = {0.125, 1., 1., 1., 5.99242, 1.};
    double uL[6] = {0., -2., 0., 0., 19.5975, -1.};
    double uR[6] = {0., 2., 0., 0., -6.19633, 1.};
    double pL[6] = {1., 0.4, 1000., 0.01, 460.894, 1.e-6};
    double pR[6] = {0.1, 0.4, 0.01, 100., 46.0950, 1.0005e-6};
    double rhosol[6] = {0.47969, 0.00617903, 0.615719, 0.61577, 12.743, 0.};
    double usol[6] = {0.841194, 0., 18.2812, -5.78011, 8.56045, 0.};
    double psol[6] = {0.293945, 8.32249e-05, 445.626, 44.5687, 1841.82, 0.};
    for(unsigned int i = 0; i < 6; i++) {
        cout << "Test problem " << i + 1 << endl;
        StateVector WL, WR;
        WL.set_rho(rhoL[i]);
        WL.set_vx(uL[i]);
        WL.set_vy(2.);
        WL.set_p(pL[i]);
        WR.set_rho(rhoR[i]);
        WR.set_vx(uR[i]);
        WR.set_vy(1.);
        WR.set_p(pR[i]);
        cout << "WL: " << WL.rho() << " " << WL.vx() << " " << WL.p() << "\t("
             << WL.vy() << ")" << endl;
        cout << "WR: " << WR.rho() << " " << WR.vx() << " " << WR.p() << "\t("
             << WR.vy() << ")" << endl;
        double mach;
        Vec n;
        n[0] = 1.;
        StateVector solution = solve(WL, WR, n, mach);
        cout << "solution : " << solution.rho() << " " << solution.vx() << " "
             << solution.p() << "\t(" << solution.vy() << ")" << endl;
        cout << "should be: " << rhosol[i] << " " << usol[i] << " " << psol[i]
             << endl;
    }
}

/**
 * @brief Get the conserved variables corresponding to the given primitive
 * variables and cell volume
 *
 * @param volume Volume of the cell
 * @param W Primitive variables in the cell
 * @return Conserved variables StateVector
 */
StateVector GottliebGrothRiemannSolver::get_Q(double volume,
                                              const StateVector& W) {
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
    Q.set_paq(W.paq() * Q.m());
    return Q;
}

/**
 * @brief Get the primitive variables corresponding to the given conserved
 * variables and cell volume
 *
 * @param volume Volume in the cell
 * @param Q Conserved variables in the cell
 * @param use_energy Flag indicating if we should use the energy or entropy
 * formalism to convert energy to pressure
 * @return Primitive variables StateVector
 */
StateVector GottliebGrothRiemannSolver::get_W(double volume, StateVector& Q,
                                              bool use_energy) {
    StateVector W;
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
    W.set_paq(Q.paq() / Q.m());
    if(W.p() < 0.) {
        cerr << "A negative value for the pressure was obtained. We better "
                "stop."
             << endl;
        my_exit();
    }
    return W;
}

/**
 * @brief Get the energy corresponding to the given primitive variables
 *
 * @param W Primitive variables StateVector
 * @return Energy
 */
double GottliebGrothRiemannSolver::get_energy(const StateVector& W) {
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
 * @brief Get the given component of the flux through a face with given velocity
 * and given primitive variables at the interface
 *
 * @param v Velocity of the face
 * @param index Index of requested flux component
 * @param W Primitive variables at the interface
 * @return Component of the flux through the interface
 */
StateVector GottliebGrothRiemannSolver::get_flux(const Vec& v,
                                                 unsigned int index,
                                                 const StateVector& W) {
    StateVector F;
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
    return F;
}
