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
 * @file HLLCRiemannSolver.cpp
 *
 * @brief HLLC RiemannSolver: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "HLLCRiemannSolver.hpp"
#include "Error.hpp"        // for my_exit
#include "RestartFile.hpp"  // for RestartFile
#include "Vec.hpp"          // for Vec
#include <algorithm>        // for max
#include <cmath>            // for pow, sqrt, fabs
#include <iostream>         // for operator<<, basic_ostream, etc

/**
 * @brief Solve the vacuum Riemann problem for the flux
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param n Interface normal
 * @param v Interface velocity
 * @param uL Left state velocity along interface normal
 * @param uR Right state velocity along interface normal
 * @param aL Left sound speed
 * @param aR Right sound speed
 * @return Deboosted flux in the lab frame
 */
StateVector HLLCRiemannSolver::vacuum_flux(StateVector& WL, StateVector& WR,
                                           Vec& n, Vec& v, double uL, double uR,
                                           double aL, double aR) {
    // STEP 1: solve Riemann problem
    StateVector solution;

    // Left and right state vacuum has already been treated

    double vhalf;
    if(!WR.rho()) {
        solution = WL;
        // vacuum right state
        if(aL > uL) {
            double SL = uL + 2. * aL / (_gamma - 1.);
            if(SL > 0.) {
                solution.set_rho(
                        WL.rho() *
                        pow(2. / (_gamma + 1.) +
                                    (_gamma - 1.) / (_gamma + 1.) / aL * uL,
                            2. / (_gamma - 1.)));
                vhalf = 2. / (_gamma + 1.) * (aL + 0.5 * (_gamma - 1.) * uL) -
                        uL;
                solution.set_p(
                        WL.p() *
                        pow(2. / (_gamma + 1.) +
                                    (_gamma - 1.) / (_gamma + 1.) / aL * uL,
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
            if(-uR < aR) {
                double SR = uR - 2. * aR / (_gamma - 1.);
                if(SR >= 0.) {
                    // vacuum
                    solution.reset();
                    vhalf = 0.;
                } else {
                    solution.set_rho(
                            WR.rho() *
                            pow(2. / (_gamma + 1.) -
                                        (_gamma - 1.) / (_gamma + 1.) / aR * uR,
                                2. / (_gamma - 1.)));
                    vhalf = 2. / (_gamma + 1.) *
                                    (-aR + 0.5 * (_gamma - 1.) * uR) -
                            uR;
                    solution.set_p(
                            WR.p() *
                            pow(2. / (_gamma + 1.) -
                                        (_gamma - 1.) / (_gamma + 1.) / aR * uR,
                                2. * _gamma / (_gamma - 1.)));
                }
            } else {
                // solution = WR
                vhalf = 0.;
            }
        } else {
            // vacuum generation
            double SR = uR - 2. * aR / (_gamma - 1.);
            double SL = uL + 2. * aL / (_gamma - 1.);
            if(SR > 0. && SL < 0.) {
                solution.reset();
                vhalf = 0.;
            } else {
                if(SL >= 0.) {
                    solution = WL;
                    if(uL < aL) {
                        solution.set_rho(WL.rho() *
                                         pow(2. / (_gamma + 1.) +
                                                     (_gamma - 1.) /
                                                             (_gamma + 1.) /
                                                             aL * uL,
                                             2. / (_gamma - 1.)));
                        vhalf = 2. / (_gamma + 1.) *
                                        (aL + 0.5 * (_gamma - 1.) * uL) -
                                uL;
                        solution.set_p(WL.p() *
                                       pow(2. / (_gamma + 1.) +
                                                   (_gamma - 1.) /
                                                           (_gamma + 1.) / aL *
                                                           uL,
                                           2. * _gamma / (_gamma - 1.)));
                    } else {
                        // solution = WL
                        vhalf = 0.;
                    }
                } else {
                    solution = WR;
                    if(-uR < aR) {
                        solution.set_rho(WR.rho() *
                                         pow(2. / (_gamma + 1.) -
                                                     (_gamma - 1.) /
                                                             (_gamma + 1.) /
                                                             aR * uR,
                                             2. / (_gamma - 1.)));
                        vhalf = 2. / (_gamma + 1.) *
                                        (-aR + 0.5 * (_gamma - 1.) * uR) -
                                uR;
                        solution.set_p(WR.p() *
                                       pow(2. / (_gamma + 1.) -
                                                   (_gamma - 1.) /
                                                           (_gamma + 1.) / aR *
                                                           uR,
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

    StateVector flux;

    if(!solution.rho()) {
        return flux;
    }

    double vsol = solution[1] * n[0] + solution[2] * n[1];
    double v2 = solution[1] * solution[1] + solution[2] * solution[2];
#if ndim_ == 3
    vsol += solution[3] * n[2];
    v2 += solution[3] * solution[3];
#endif
    flux[0] = solution.rho() * vsol;
    flux[1] = solution.rho() * solution.vx() * vsol + solution.p() * n[0];
    flux[2] = solution.rho() * solution.vy() * vsol + solution.p() * n[1];
#if ndim_ == 3
    flux[3] = solution.rho() * solution.vz() * vsol + solution.p() * n[2];
#endif
    double e = solution.p() / (_gamma - 1.) / solution.rho() + 0.5 * v2;
    flux[ndim_ + 1] = solution.rho() * e * vsol + solution.p() * vsol;

    // deboost to lab frame
    // we add the flux contribution due to the movement of the interface
    // the density flux is unchanged
    // we add the extra velocity flux due to the absolute motion of the fluid
    // similarly, we need to add the energy fluxes due to the absolute motion
    flux[1] += v[0] * flux[0];
    flux[2] += v[1] * flux[0];
#if ndim_ == 3
    flux[3] += v[2] * flux[0];
    flux[4] += v[0] * flux[1] + v[1] * flux[2] + v[2] * flux[3] +
               0.5 * v.norm2() * flux[0];
#else
    flux[3] += v[0] * flux[1] + v[1] * flux[2] + 0.5 * v.norm2() * flux[0];
#endif

    return flux;
}

/**
 * @brief Constructor
 *
 * @param gamma Adiabatic index
 */
HLLCRiemannSolver::HLLCRiemannSolver(double gamma) {
    _gamma = gamma;
    _neval = 0;
}

/**
 * @brief Destructor
 *
 * Print some general statistics
 */
HLLCRiemannSolver::~HLLCRiemannSolver() {
    cout << _neval << " Riemann solver evaluations" << endl;
}

/**
 * @brief Dummy classical solver
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param n Interface normal
 * @param mach Mach number
 * @param reflective Reflective flag
 * @return Dummy solution: average of left and right state
 */
StateVector HLLCRiemannSolver::solve(StateVector& WL, StateVector& WR, Vec& n,
                                     double& mach, bool reflective) {
    cerr << "This solver should not be used!" << endl;
    my_exit();
    return 0.5 * (WL + WR);
}

/**
 * @brief Solve for the flux directly
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param n Interface normal
 * @param v Interface velocity
 * @param reflective Flag indicating if the right state should be a
 * reflective copy of the left state
 * @return Flux along the interface normal
 */
StateVector HLLCRiemannSolver::solve_for_flux(StateVector& WL, StateVector& WR,
                                              Vec& n, Vec& v, bool reflective) {
    _neval++;

    // Handle vacuum
    if(!WL.rho() && (reflective || !WR.rho())) {
        return StateVector();
    }

// STEP 0: obtain velocity in interface frame and apply boundary conditions
#if ndim_ == 3
    double uL = WL[1] * n[0] + WL[2] * n[1] + WL[3] * n[2];
    double uR = WR[1] * n[0] + WR[2] * n[1] + WR[3] * n[2];
#else
    double uL = WL[1] * n[0] + WL[2] * n[1];
    double uR = WR[1] * n[0] + WR[2] * n[1];
#endif
    if(reflective) {
        WR = WL;
        uR = -uL;
    }
    double aL = get_soundspeed(WL);
    double aR = get_soundspeed(WR);

    // Handle vacuum: vacuum does not require iteration and is always exact
    if(!WL.rho() || !WR.rho()) {
        return vacuum_flux(WL, WR, n, v, uL, uR, aL, aR);
    }
    if(2. * aL / (_gamma - 1.) + 2. * aR / (_gamma - 1.) < fabs(uL - uR)) {
        return vacuum_flux(WL, WR, n, v, uL, uR, aL, aR);
    }

    // STEP 1: pressure estimate
    double rhobar = 0.5 * (WL.rho() + WR.rho());
    double abar = 0.5 * (aL + aR);
    double pPVRS = 0.5 * (WL.p() + WR.p()) - 0.5 * (uR - uL) * rhobar * abar;
    double pstar = std::max(0., pPVRS);

    // STEP 2: wave speed estimates
    // all these speeds are along the interface normal, since uL and uR are
    double qL = 1.;
    if(pstar > WL.p()) {
        qL = sqrt(1. + 0.5 * (_gamma + 1.) / _gamma * (pstar / WL.p() - 1.));
    }
    double qR = 1.;
    if(pstar > WR.p()) {
        qR = sqrt(1. + 0.5 * (_gamma + 1.) / _gamma * (pstar / WR.p() - 1.));
    }
    double SL = uL - aL * qL;
    double SR = uR + aR * qR;
    double Sstar = (WR.p() - WL.p() + WL.rho() * uL * (SL - uL) -
                    WR.rho() * uR * (SR - uR)) /
                   (WL.rho() * (SL - uL) - WR.rho() * (SR - uR));

    // STEP 3: HLLC flux in a frame moving with the interface velocity
    StateVector flux;
    if(Sstar >= 0.) {
        // flux FL
        flux[0] = WL.rho() * uL;
        // these are the actual correct fluxes in the boosted lab frame
        // (not rotated to interface frame)
        flux[1] = WL.rho() * WL.vx() * uL + WL.p() * n[0];
        flux[2] = WL.rho() * WL.vy() * uL + WL.p() * n[1];
        double v2 = WL.vx() * WL.vx() + WL.vy() * WL.vy();
#if ndim_ == 3
        flux[3] = WL.rho() * WL.vz() * uL + WL.p() * n[2];
        v2 += WL.vz() * WL.vz();
#endif
        double eL = WL.p() / (_gamma - 1.) / WL.rho() + 0.5 * v2;
        flux[ndim_ + 1] = WL.rho() * eL * uL + WL.p() * uL;
        if(SL < 0.) {
            // add flux FstarL
            StateVector UstarL;
            UstarL[0] = 1.;
            // we need UstarL in the lab frame:
            // subtract the velocity in the interface frame from the lab frame
            // velocity and then add Sstar in interface frame
            UstarL[1] = WL.vx() + (Sstar - uL) * n[0];
            UstarL[2] = WL.vy() + (Sstar - uL) * n[1];
#if ndim_ == 3
            UstarL[3] = WL.vz() + (Sstar - uL) * n[2];
#endif
            UstarL[ndim_ + 1] =
                    eL +
                    (Sstar - uL) * (Sstar + WL.p() / (WL.rho() * (SL - uL)));
            UstarL *= WL.rho() * (SL - uL) / (SL - Sstar);
            flux[0] += SL * (UstarL[0] - WL.rho());
            flux[1] += SL * (UstarL[1] - WL.rho() * WL.vx());
            flux[2] += SL * (UstarL[2] - WL.rho() * WL.vy());
#if ndim_ == 3
            flux[3] += SL * (UstarL[3] - WL.rho() * WL.vz());
#endif
            flux[ndim_ + 1] += SL * (UstarL[ndim_ + 1] - WL.rho() * eL);
        }
    } else {
        // flux FR
        flux[0] = WR.rho() * uR;
        flux[1] = WR.rho() * WR.vx() * uR + WR.p() * n[0];
        flux[2] = WR.rho() * WR.vy() * uR + WR.p() * n[1];
        double v2 = WR.vx() * WR.vx() + WR.vy() * WR.vy();
#if ndim_ == 3
        flux[3] = WR.rho() * WR.vz() * uR + WR.p() * n[2];
        v2 += WR.vz() * WR.vz();
#endif
        double eR = WR.p() / (_gamma - 1.) / WR.rho() + 0.5 * v2;
        flux[ndim_ + 1] = WR.rho() * eR * uR + WR.p() * uR;
        if(SR > 0.) {
            // add flux FstarR
            StateVector UstarR;
            UstarR[0] = 1.;
            UstarR[1] = WR.vx() + (Sstar - uR) * n[0];
            UstarR[2] = WR.vy() + (Sstar - uR) * n[1];
#if ndim_ == 3
            UstarR[3] = WR.vz() + (Sstar - uR) * n[2];
#endif
            UstarR[ndim_ + 1] =
                    eR +
                    (Sstar - uR) * (Sstar + WR.p() / (WR.rho() * (SR - uR)));
            UstarR *= WR.rho() * (SR - uR) / (SR - Sstar);
            flux[0] += SR * (UstarR[0] - WR.rho());
            flux[1] += SR * (UstarR[1] - WR.rho() * WR.vx());
            flux[2] += SR * (UstarR[2] - WR.rho() * WR.vy());
#if ndim_ == 3
            flux[3] += SR * (UstarR[3] - WR.rho() * WR.vz());
#endif
            flux[ndim_ + 1] += SR * (UstarR[ndim_ + 1] - WR.rho() * eR);
        }
    }

    // deboost to lab frame
    // we add the flux contribution due to the movement of the interface
    // the density flux is unchanged
    // we add the extra velocity flux due to the absolute motion of the fluid
    // similarly, we need to add the energy fluxes due to the absolute motion
    flux[1] += v[0] * flux[0];
    flux[2] += v[1] * flux[0];
#if ndim_ == 3
    flux[3] += v[2] * flux[0];
    flux[4] += v[0] * flux[1] + v[1] * flux[2] + v[2] * flux[3] +
               0.5 * v.norm2() * flux[0];
#else
    flux[3] += v[0] * flux[1] + v[1] * flux[2] + 0.5 * v.norm2() * flux[0];
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
StateVector HLLCRiemannSolver::get_time_derivative(
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

    return dWdt;
}

/**
 * @brief Get the sound speed
 *
 * @param W StateVector
 * @return Sound speed
 */
double HLLCRiemannSolver::get_soundspeed(const StateVector& W) {
    return sqrt(_gamma * W.p() / W.rho());
}

/**
 * @brief Test the Riemann solver on 6 problems with known solutions
 */
void HLLCRiemannSolver::test() {

  cout <<  "This test is empty ! Implement me !" << endl;

}

/**
 * @brief Convert primitive to conserved variables
 *
 * @warning Dummy!
 *
 * @param volume Volume
 * @param W Primitive variables
 * @return Conserved variables
 */
StateVector HLLCRiemannSolver::get_Q(double volume, const StateVector& W) {
    StateVector Q;
    Q[0] = W.rho() * volume;
    Q[1] = Q.m() * W.vx();
    Q[2] = Q.m() * W.vy();
    double v2 = W.vx() * W.vx() + W.vy() * W.vy();
#if ndim_ == 3
    Q[3] = Q.m() * W.vz();
    v2 += W.vz() * W.vz();
#endif
    Q[ndim_ + 1] = volume * W.p() / (_gamma - 1.) + 0.5 * Q.m() * v2;
    return Q;
}

/**
 * @brief Convert conserved to primitive variables
 *
 * @warning Dummy!
 *
 * @param volume Volume
 * @param Q Conserved variables
 * @param use_energy Energy flag
 * @return Primitive variables
 */
StateVector HLLCRiemannSolver::get_W(double volume, StateVector& Q,
                                     bool use_energy) {
    StateVector W;
    // handle vacuum
    if(!Q[0]) {
        return W;
    }
    W[0] = Q.m() / volume;
    W[1] = Q.px() / Q.m();
    W[2] = Q.py() / Q.m();
    double p2 = Q.px() * Q.px() + Q.py() * Q.py();
#if ndim_ == 3
    W[3] = Q.pz() / Q.m();
    p2 += Q.pz() * Q.pz();
#endif
    W[ndim_ + 1] = (_gamma - 1.) * (Q.e() - 0.5 * p2 / Q.m()) / volume;

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
        W.set_rho(0.);
        W.set_vx(0.);
        W.set_vy(0.);
#if ndim_ == 3
        W.set_vz(0.);
#endif
        W.set_p(0.);
        W.set_paq(0.);
    }
    return W;
}

/**
 * @brief Get the flux
 *
 * @warning Dummy!
 *
 * @param v Interface velocity
 * @param index Required index
 * @param W Primitive variables
 * @return Flux
 */
StateVector HLLCRiemannSolver::get_flux(const Vec& v, unsigned int index,
                                        const StateVector& W) {
    cerr << "This function should not be used!" << endl;
    my_exit();
    return W;
}

/**
 * @brief Get the adiabatic index
 *
 * @return Adiabatic index
 */
double HLLCRiemannSolver::get_gamma() {
    return _gamma;
}

/**
 * @brief Dump the HLLCRiemannSolver to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void HLLCRiemannSolver::dump(RestartFile& rfile) {
    rfile.write(_gamma);
    rfile.write(_neval);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
HLLCRiemannSolver::HLLCRiemannSolver(RestartFile& rfile) {
    rfile.read(_gamma);
    rfile.read(_neval);
}

/**
 * @brief Get the number of HLLCRiemannSolver evaluations
 *
 * @return Number of evaluations
 */
unsigned long HLLCRiemannSolver::get_neval() {
    return _neval;
}
