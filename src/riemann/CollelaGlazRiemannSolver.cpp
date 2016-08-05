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
 * @file CollelaGlazRiemannSolver.cpp
 *
 * @brief Collela-Glaz Riemann solver implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "CollelaGlazRiemannSolver.hpp"
#include <algorithm>  // for max, min
#include <cmath>      // for fabs, sqrt
#include <iostream>   // for operator<<, basic_ostream, etc
using namespace std;

/**
 * @brief Constructor
 *
 * @param tolerance Tolerance value
 */
CollelaGlazRiemannSolver::CollelaGlazRiemannSolver(double tolerance) {
    _tolerance = tolerance;
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
StateVector CollelaGlazRiemannSolver::solve(StateVector& WL, StateVector& WR,
                                            Vec& n, double& mach,
                                            bool reflective) {
    StateVector W;
    W[ndim_ + 1] = calculate_p_star(WL, WR);
    return W;
}

/**
 * @brief Derivative of the pressure w.r.t. the thermal energy
 *
 * @warning Not implemented (yet)!
 *
 * @param W Primitive variables
 * @return Derivative of the pressure
 */
double CollelaGlazRiemannSolver::dp_du(const StateVector& W) {
    return 0.;
}

/**
 * @brief Derivative of the pressure w.r.t. the density
 *
 * @warning Not implemented (yet)!
 *
 * @param W Primitive variables
 * @return Derivative of the pressure
 */
double CollelaGlazRiemannSolver::dp_drho(const StateVector& W) {
    return 0.;
}

/**
 * @brief Get the soundspeed corresponding to the given StateVector
 *
 * @warning Not correctly implemented (yet)!
 *
 * @param W Primitive variables
 * @return Soundspeed
 */
double CollelaGlazRiemannSolver::get_soundspeed(const StateVector& W) {
    // c^2 = 1/rho^2*p*dp/du + dp/drho
    //    return sqrt(W.p()*dp_du(W)/W.rho()/W.rho() + dp_drho(W));
    // polytropic EOS c^2 = 1.4 p/rho
    return sqrt(1.4 * W.p() / W.rho());
    //    return sqrt(5.*W.p()/3./W.rho());
}

/**
 * @brief Test the Riemann solver
 */
void CollelaGlazRiemannSolver::test() {
    cout << "Testing the Riemann solver" << endl;
    double rhoL[5] = {1., 1., 1., 1., 5.99924};
    double rhoR[5] = {0.125, 1., 1., 1., 5.99242};
    double uL[5] = {0., -2., 0., 0., 19.5975};
    double uR[5] = {0., 2., 0., 0., -6.19633};
    double pL[5] = {1., 0.4, 1000., 0.01, 460.894};
    double pR[5] = {0.1, 0.4, 0.01, 100., 46.0950};
    //    double rhosol[5] = {0.47969, 0.00617903, 0.615719, 0.61577, 12.743};
    //    double usol[5] = {0.841194, 0., 18.2812, -5.78011, 8.56045};
    //    double psol[5] = {0.293945, 8.32249e-05, 445.626, 44.5687, 1841.82};
    double psol[5] = {0.30313, 0.00189, 460.894, 46.095, 1691.64};
    for(unsigned int i = 0; i < 5; i++) {
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
        double pstar = calculate_p_star(WL, WR);
        cout << "pstar: " << pstar << ", should be " << psol[i] << "\t"
             << (fabs(pstar - psol[i]) / psol[i]) << endl;
        //        StateVector solution = solve(WL, WR);
        //        cout << "solution : " << solution.rho() << " " <<
        // solution.vx() << " "
        //        << solution.p() << "\t(" << solution.vy() << ")" << endl;
        //        cout << "should be: " << rhosol[i] << " " << usol[i] << " " <<
        // psol[i]
        //        << endl;
    }
}

/**
 * @brief Get the conserved variables corresponding to the given primitive
 * variables and volume
 *
 * @warning Not implemented!
 *
 * @param volume Volume of a cell
 * @param W Primitive variables of a cell
 * @return Conserved variables StateVector
 */
StateVector CollelaGlazRiemannSolver::get_Q(double volume,
                                            const StateVector& W) {
    StateVector Q;
    return Q;
}

/**
 * @brief Get the primitive variables corresponding to the given conserved
 * variables and volume
 *
 * @warning Not implemented!
 *
 * @param volume Volume in a cell
 * @param Q Conserved variables in a cell
 * @param use_energy Flag indicating if the energy or entropy formulation should
 * be used to convert total energy to pressure
 * @return
 */
StateVector CollelaGlazRiemannSolver::get_W(double volume, StateVector& Q,
                                            bool use_energy) {
    StateVector W;
    return W;
}

/**
 * @brief Get the given component of the flux matrix through the interface with
 * the given primitive variables and velocity
 *
 * @warning Not implemented!
 *
 * @param v Velocity of the interface
 * @param index unsigned integer index of the requested component
 * @param W Primitive variables at the interface
 * @return Flux through the interface
 */
StateVector CollelaGlazRiemannSolver::get_flux(const Vec& v, unsigned int index,
                                               const StateVector& W) {
    StateVector flux;
    return flux;
}

/**
 * @brief Calculate the pressure in the intermediate state
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @return Pressure in the intermediate state
 */
double CollelaGlazRiemannSolver::calculate_p_star(const StateVector& WL,
                                                  const StateVector& WR) {
    // secant method
    // the initial guesses are shamelessly stolen from a piece of code on
    // http://www.astro.sunysb.edu/mzingale/phy688/riemann.f90
    double pstar[2] = {0.};
    double ustarR[2] = {0.};
    double ustarL[2] = {0.};

    // small gamma
    double gammaL = WL.p() / WL.rho() / get_energy(WL) + 1.;
    double gammaR = WR.p() / WR.rho() / get_energy(WR) + 1.;

    double CR = WR.rho() * get_soundspeed(WR);
    double CL = WL.rho() * get_soundspeed(WL);

    //    pstar[0] = WR.p() - WL.p() - CR*(WR.vx() - WL.vx());
    //    pstar[0] = WL.p() + pstar[0]*(CL/(CL+CR));
    //    pstar[0] = std::max(_tolerance, pstar[0]);
    pstar[0] = 0.5 * (WR.p() + WL.p());

    // capital Gamma
    double GammaL = CL * CL / WL.rho() / WL.p();
    double GammaR = CR * CR / WR.rho() / WR.p();

    double SL = get_wavespeed(pstar[0], WL, gammaL, gammaR, GammaL, GammaR);
    double SR = get_wavespeed(pstar[0], WR, gammaR, gammaL, GammaR, GammaL);

    ustarR[0] = WR.vx() + (pstar[0] - WR.p()) / SR;
    ustarL[0] = WL.vx() - (pstar[0] - WL.p()) / SL;

    pstar[1] = WR.p() - WL.p() - SR * (WR.vx() - WL.vx());
    pstar[1] = WL.p() + pstar[1] * SL / (SL + SR);
    pstar[1] = std::max(_tolerance, pstar[1]);

    SL = get_wavespeed(pstar[1], WL, gammaL, gammaR, GammaL, GammaR);
    SR = get_wavespeed(pstar[1], WR, gammaR, gammaL, GammaR, GammaL);

    ustarR[1] = WR.vx() + (pstar[1] - WR.p()) / SR;
    ustarL[1] = WL.vx() - (pstar[1] - WL.p()) / SL;

    if(pstar[0] == pstar[1]) {
        pstar[0] = 0.5 * (WL.p() + WR.p());
        ustarR[0] = 0.5 * (WL.vx() + WR.vx());
        ustarL[0] = ustarR[0];
    }

    unsigned int loopcount = 0;
    while(fabs(pstar[0] - pstar[1]) >=
          0.5 * _tolerance * (pstar[0] + pstar[1])) {
        double ptemp = pstar[1];
        double uRtemp = ustarR[1];
        double uLtemp = ustarL[1];

        pstar[1] = pstar[1] -
                   (ustarR[1] - ustarL[1]) * (fabs(pstar[1] - pstar[0]) /
                                              (fabs(ustarL[1] - ustarL[0]) +
                                               fabs(ustarR[1] - ustarR[0])));

        SR = get_wavespeed(pstar[1], WR, gammaR, gammaL, GammaR, GammaL);
        SL = get_wavespeed(pstar[1], WL, gammaL, gammaR, GammaL, GammaR);

        ustarR[1] = WR.vx() + (pstar[1] - WR.p()) / SR;
        ustarL[1] = WL.vx() - (pstar[1] - WL.p()) / SL;

        pstar[0] = ptemp;
        ustarR[0] = uRtemp;
        ustarL[0] = uLtemp;
        loopcount++;
    }

    cout << "used " << loopcount << " iterations" << endl;

    return pstar[1];
}

/**
 * @brief Get the wave speed in the intermediate state
 *
 * @warning No idea what this function does...
 *
 * @param pstar Pressure in the intermediate state
 * @param W Primitive variables
 * @param gammaS Unknown parameter
 * @param gammaO Unknown parameter
 * @param GammaS Unknown parameter
 * @param GammaO Unknown parameter
 * @return Wave speed in the intermediate state
 */
double CollelaGlazRiemannSolver::get_wavespeed(double pstar,
                                               const StateVector& W,
                                               double gammaS, double gammaO,
                                               double GammaS, double GammaO) {
    double gamma = 0.5 * (gammaS + gammaO);
    gamma = gammaS +
            2. * (1. - 2. * gamma / (GammaS + GammaO)) * (gamma - 1.) *
                    (pstar - W.p()) / (pstar + W.p());
    // check the min and max values for the specific equations you want to use
    // (maybe set these as parameters)
    gamma = std::max(1., std::min(gamma, 1.66667));

    //    if(pstar > W.p()){
    //        double S = (pstar-W.p())*(pstar + 0.5*(gamma-1.)*(pstar+W.p()))/
    //    (pstar - (gamma-1.)/(gammaS-1.)*W.p())*W.rho();
    //        return sqrt(S);
    //    } else {
    //        double S = 2.*sqrt(gamma*W.p()/W.rho())/(gamma-1.)*
    //    (pow(pstar/W.p(), 0.5*(gamma-1.)/gamma)-1.);
    //        return (pstar-W.p())/S;
    //    }
    double S = (pstar - W.p()) *
               (pstar + 0.5 * (gamma - 1.) * (pstar + W.p())) /
               (pstar - (gamma - 1.) / (gammaS - 1.) * W.p()) * W.rho();
    return sqrt(S);
}

/**
 * @brief Get the total energy corresponding to the given primitive variables
 *
 * @param W Primitive variables
 * @return Total energy
 */
double CollelaGlazRiemannSolver::get_energy(const StateVector& W) {
    // we temporarily use a polytropic EOS p = 0.4 rho u <=> u = 5/2 p/rho
    return 2.5 * W.p() / W.rho();
    //    return 1.5*W.p()/W.rho();
}
