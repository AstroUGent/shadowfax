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
 * @file NoPressureRiemannSolver.cpp
 *
 * @brief Riemann solver for a fluid without pressure: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "NoPressureRiemannSolver.hpp"
#include "ProgramLog.hpp"   // for LOGS
#include "RestartFile.hpp"  // for RestartFile
#include "Vec.hpp"          // for Vec, operator*, operator/
#include <iostream>         // for operator<<, basic_ostream, etc
using namespace std;

/**
 * @brief Constructor
 */
NoPressureRiemannSolver::NoPressureRiemannSolver() {
    _counter = 0;

    LOGS("NoPressureRiemannSolver initialized");
}

/**
 * @brief Destructor
 */
NoPressureRiemannSolver::~NoPressureRiemannSolver() {
    cout << _counter << " fake Riemann solver evaluations" << endl;

    LOGS("NoPressureRiemannSolver destructed");
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
StateVector NoPressureRiemannSolver::solve(StateVector& WL, StateVector& WR,
                                           Vec& n, double& mach0,
                                           bool reflective) {
    _counter++;
    _timer.start();
    StateVector Wstar;
    if(reflective) {
        Wstar = WL;
    } else {
        Wstar = 0.5 * (WL + WR);
    }
    _timer.stop();
    return Wstar;
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
StateVector NoPressureRiemannSolver::solve_for_flux(StateVector& WL,
                                                    StateVector& WR, Vec& n,
                                                    Vec& v, bool reflective) {
    double maxmach;
    StateVector Whalf = solve(WL, WR, n, maxmach, reflective);

    StateVector fluxes[3];
    fluxes[0] = get_flux(v, 0, Whalf);
    fluxes[1] = get_flux(v, 1, Whalf);
    fluxes[2] = get_flux(v, 2, Whalf);

    StateVector flux;
#if ndim_ == 3
    flux[0] = fluxes[0][0] * n[0] + fluxes[1][0] * n[1] + fluxes[2][0] * n[2];
    flux[1] = fluxes[0][1] * n[0] + fluxes[1][1] * n[1] + fluxes[2][1] * n[2];
    flux[2] = fluxes[0][2] * n[0] + fluxes[1][2] * n[1] + fluxes[2][2] * n[2];
    flux[3] = fluxes[0][3] * n[0] + fluxes[1][3] * n[1] + fluxes[2][3] * n[2];
    flux[4] = fluxes[0][4] * n[0] + fluxes[1][4] * n[1] + fluxes[2][4] * n[2];
#else
    flux[0] = fluxes[0][0] * n[0] + fluxes[1][0] * n[1];
    flux[1] = fluxes[0][1] * n[0] + fluxes[1][1] * n[1];
    flux[2] = fluxes[0][2] * n[0] + fluxes[1][2] * n[1];
    flux[3] = fluxes[0][3] * n[0] + fluxes[1][3] * n[1];
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
StateVector NoPressureRiemannSolver::get_time_derivative(
        const StateVector& W, const StateVector* gradients) {
    return StateVector();
}

/**
 * @brief Get the sound speed corresponding to the given state
 *
 * 0 in this case, since the fluid has no pressure and hence no travelling waves
 *
 * @param W Input StateVector
 * @return Soundspeed
 */
double NoPressureRiemannSolver::get_soundspeed(const StateVector& W) {
    return 0.;
}

/**
 * @brief Test the Riemann solver
 *
 * Does nothing.
 */
void NoPressureRiemannSolver::test() {
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
StateVector NoPressureRiemannSolver::get_Q(double volume,
                                           const StateVector& W) {
    double m = volume * W.rho();
#if ndim_ == 3
    Vec v(W.vx(), W.vy(), W.vz());
#else
    Vec v(W.vx(), W.vy());
#endif
    Vec p = m * v;
    double e = 0.5 * m * v.norm2();
#if ndim_ == 3
    StateVector Q(m, p.x(), p.y(), p.z(), e);
#else
    StateVector Q(m, p.x(), p.y(), e);
#endif
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
StateVector NoPressureRiemannSolver::get_W(double volume, StateVector& Q,
                                           bool use_energy) {
    if(!Q.m()) {
        // vacuum
        return StateVector();
    }
    double rho = Q.m() / volume;
#if ndim_ == 3
    Vec p(Q.px(), Q.py(), Q.pz());
#else
    Vec p(Q.px(), Q.py());
#endif
    Vec v = p / Q.m();
#if ndim_ == 3
    StateVector W(rho, v.x(), v.y(), v.z(), 0.);
#else
    StateVector W(rho, v.x(), v.y(), 0.);
#endif
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
StateVector NoPressureRiemannSolver::get_flux(const Vec& v, unsigned int index,
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
 * @return 0., since the fluid has no pressure and an adiabatic index is hence
 * meaningless
 */
double NoPressureRiemannSolver::get_gamma() {
    return 0.;
}

/**
 * @brief Get the number of Riemann solver evaluations since the start of the
 * simulation
 *
 * @return Number of solver evaluations
 */
unsigned long NoPressureRiemannSolver::get_neval() {
    return _counter;
}

/**
 * @brief Dump the solver to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void NoPressureRiemannSolver::dump(RestartFile& rfile) {
    _timer.dump(rfile);
    rfile.write(_counter);

    LOGS("NoPressureRiemannSolver dumped");
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
NoPressureRiemannSolver::NoPressureRiemannSolver(RestartFile& rfile)
        : _timer(rfile) {
    rfile.read(_counter);

    LOGS("NoPressureRiemannSolver restarted");
}
