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
 * @file AdvectionRiemannSolver.cpp
 *
 * @brief Riemann solver that solves the advection equation instead of the
 * Euler equation: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdvectionRiemannSolver.hpp"
#include "RestartFile.hpp"  // for RestartFile
#include "Vec.hpp"          // for Vec
#include <iostream>         // for operator<<, basic_ostream, etc
using namespace std;

/**
 * @brief Constructor
 */
AdvectionRiemannSolver::AdvectionRiemannSolver(double gamma) {
    _counter = 0;
    _gamma = gamma;
}

/**
 * @brief Destructor
 */
AdvectionRiemannSolver::~AdvectionRiemannSolver() {
    cout << _counter << " fake Riemann solver evaluations" << endl;
}

/**
 * @brief Solve the Riemann problem with given left and right states
 *
 * @param WL Left state
 * @param WR Right state
 * @param n Interface normal vector
 * @param mach0 Maximum mach number of all shocks involved (ignored for this
 * solver)
 * @param reflective Flag indicating whether the interface corresponds to a
 * reflective boundary
 * @return StateVector corresponding to the sampled Riemann solution at the
 * interface, in this case just the average of both states
 */
StateVector AdvectionRiemannSolver::solve(StateVector& WL, StateVector& WR,
                                          Vec& n, double& mach0,
                                          bool reflective) {
    _counter++;
    _timer.start();
    StateVector Wstar;
    if(reflective) {
        Wstar = WL;
    } else {
        double ustar = 0.5 * (WL.vx() + WR.vx());
        if(ustar > 0.) {
            Wstar = WL;
        } else {
            Wstar = WR;
        }
    }
    _timer.stop();
    return Wstar;
}

/**
 * @brief Solve the Riemann problem and return the flux directly
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param n Interface normal
 * @param v Interface velocity
 * @param reflective Flag indicating if the right state should be a
 * reflective copy of the left state
 * @return Flux along the interface normal
 */
StateVector AdvectionRiemannSolver::solve_for_flux(StateVector& WL,
                                                   StateVector& WR, Vec& n,
                                                   Vec& v, bool reflective) {
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
StateVector AdvectionRiemannSolver::get_time_derivative(
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
 * @brief Get the sound speed corresponding to the given state
 *
 * 0 in this case, since the fluid has no pressure and hence no travelling waves
 *
 * @param W Input StateVector
 * @return Soundspeed
 */
double AdvectionRiemannSolver::get_soundspeed(const StateVector& W) {
    return sqrt(_gamma * W.p() / W.rho());
}

/**
 * @brief Test the Riemann solver
 *
 * Does nothing.
 */
void AdvectionRiemannSolver::test() {
    // do nothing
}

/**
 * @brief Get the conserved variables corresponding to the given primitive
 * variables and volume
 *
 * @param volume Volume of the cell
 * @param W Primitive variables inside the cell
 * @return Conserved variables inside the cell
 */
StateVector AdvectionRiemannSolver::get_Q(double volume, const StateVector& W) {
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
    //    Q.set_paq(W.paq()*Q.m());
    // Ai = Pi/rhoi^gamma
    Q.set_paq(Q.m() * W.p() / pow(W.rho(), _gamma));

    // sanity check results
    if(Q.m() < 0. || Q.e() < 0.) {
        throw ConservedVariablesException(W, Q, volume);
    }

    return Q;
}

/**
 * @brief Get the primitive variables corresponding to the given conserved
 * variables and volume
 *
 * @param volume Volume of the cell
 * @param Q Conserved variables inside the cell
 * @param use_energy Flag indicating whether we should use the energy or entropy
 * formulation to calculate the pressure (ignored for this solver)
 * @return Primitive variables inside the cell
 */
StateVector AdvectionRiemannSolver::get_W(double volume, StateVector& Q,
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
    if(use_energy) {
        // reset entropy
        Q.set_paq(Q.m() * W.p() / pow(W.rho(), _gamma));
        W.set_paq(Q.paq() / Q.m());
    } else {
        W.set_paq(Q.paq() / Q.m());
        W.set_p(W.paq() * pow(W.rho(), _gamma));

#if ndim_ == 3
        Q.set_e(0.5 * W.rho() *
                        (W.vx() * W.vx() + W.vy() * W.vy() + W.vz() * W.vz()) +
                W.p() / (_gamma - 1.) * volume);
#else
        Q.set_e(0.5 * W.rho() * (W.vx() * W.vx() + W.vy() * W.vy()) +
                W.p() / (_gamma - 1.) * volume);
#endif
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

    // sanity check results
    //    if(W.rho() < 0. || W.p() < 0.){
    //        throw PrimitiveVariablesException(W, Q, volume);
    //    }
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
 * @brief Get the requested component of the flux of the conserved variables
 * through an interface with given velocity and primitive variables
 *
 * @param v Velocity of the interface
 * @param index Component of the flux that is requested
 * @param W Primitive variables at the interface
 * @return Requested component of the flux of the conserved variables through
 * the interface
 */
StateVector AdvectionRiemannSolver::get_flux(const Vec& v, unsigned int index,
                                             const StateVector& W) {
    StateVector F;
    F[0] = W.rho() * (W[1 + index] - v[index]);
    F[1] = W.rho() * (W[1 + index] - v[index]) * W.vx();
    F[2] = W.rho() * (W[1 + index] - v[index]) * W.vy();
#if ndim_ == 3
    F[3] = W.rho() * (W[1 + index] - v[index]) * W.vz();
    double e = 0.5 * W.rho() *
               (W.vx() * W.vx() + W.vy() * W.vy() + W.vz() * W.vz());
#else
    double e = 0.5 * W.rho() * (W.vx() * W.vx() + W.vy() * W.vy());
#endif
    F[ndim_ + 1] = (W[1 + index] - v[index]) * e;
    F.set_paq(W.rho() * (W[1 + index] - v[index]) * W.paq());
    return F;
}

/**
 * @brief Get the adiabatic index of the fluid
 *
 * @return The adiabiatic index of the fluid
 */
double AdvectionRiemannSolver::get_gamma() {
    return _gamma;
}

/**
 * @brief Get the number of Riemann solver evaluations since the start of the
 * simulation
 *
 * @return Number of solver evaluations
 */
unsigned long AdvectionRiemannSolver::get_neval() {
    return _counter;
}

/**
 * @brief Dump the solver to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void AdvectionRiemannSolver::dump(RestartFile& rfile) {
    _timer.dump(rfile);
    rfile.write(_counter);
    rfile.write(_gamma);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
AdvectionRiemannSolver::AdvectionRiemannSolver(RestartFile& rfile)
        : _timer(rfile) {
    rfile.read(_counter);
    rfile.read(_gamma);
}
