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
 * @file ApproximateSolver.cpp
 *
 * @brief Two Rarefaction Riemann Solver (TRRS): implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ApproximateSolver.hpp"
#include "RestartFile.hpp"
#include "ProgramLog.hpp"
#include "utilities/HelperFunctions.hpp"
#include <iostream>
using namespace std;

/**
 * @brief Get the energy corresponding to the given primitive variables
 *
 * @param W Primitive variables StateVector
 * @return Energy
 */
double TRRSSolver::get_energy(const StateVector &W){
#if ndim_==3
    return 0.5*W.rho()*(W.vx()*W.vx() + W.vy()*W.vy() + W.vz()*W.vz()) +
            W.p()/(_gamma-1.);
#else
    return 0.5*W.rho()*(W.vx()*W.vx() + W.vy()*W.vy()) + W.p()/(_gamma-1.);
#endif
}

/**
 * @brief Constructor
 *
 * @param gamma Adiabatic index \f$\gamma\f$ of the gas
 */
TRRSSolver::TRRSSolver(double gamma) : _gamma(gamma) {
    _gp1d2g = 0.5*(_gamma+1.)/_gamma; // gamma plus 1 divided by 2 gamma
    _gm1d2g = 0.5*(_gamma-1.)/_gamma; // gamma minus 1 divided by 2 gamma
    _gm1dgp1 = (_gamma-1.)/(_gamma+1.); // gamma minus 1 divided by gamma plus 1
    _tdgp1 = 2./(_gamma+1.); // two divided by gamma plus 1
    _tdgm1 = 2./(_gamma-1.); // two divided by gamma minus 1
    _gm1d2 = 0.5*(_gamma-1.); // gamma minus 1 divided by 2
    _tgdgm1 = 2.*_gamma/(_gamma-1.); // two times gamma divided by gamma minus 1
    _ginv = 1./_gamma; // gamma inverse

    _counter = 0;

    LOGS("TRRSSolver initialized");
}

/**
 * @brief Destructor
 *
 * Print the number of calls and the total time spent in the solver to the
 * stdout.
 */
TRRSSolver::~TRRSSolver(){
    cout << "TRRS Riemann solver: "
         << HelperFunctions::human_readable_counter(_counter)
         << " evaluations in " << _timer.value() << " seconds." << endl;

    LOGS("TRRSSolver destructed");
}

/**
 * @brief Solve the Riemann problem in between the given left and right state
 *
 * @param WL Left primitive variables StateVector
 * @param WR Right primitive variables StateVector
 * @param n Normal vector to the interface
 * @param mach0 Reference to the variable in which the maximal mach number
 * encountered shall be stored
 * @param reflective Flag indicating if the right state should be a
 * reflective copy of the left state
 * @return Primitive variables StateVector that is the solution of the stated
 * Riemann problem
 */
StateVector TRRSSolver::solve(StateVector &WL, StateVector &WR, Vec &n,
                              double &mach0, bool reflective){
    _counter++;
    _timer.start();

    StateVector Whalf;

#if ndim_==3
    double vL = WL[1]*n[0] + WL[2]*n[1] + WL[3]*n[2];
    double vR = WR[1]*n[0] + WR[2]*n[1] + WR[3]*n[2];
#else
    double vL = WL[1]*n[0] + WL[2]*n[1];
    double vR = WR[1]*n[0] + WR[2]*n[1];
#endif

    if(reflective){
        WR = WL;
        vR = -vL;
    }

    double aL = get_soundspeed(WL);
    double aR = get_soundspeed(WR);

    double PLR = pow(WL.p()/WR.p(), _gm1d2g);
    double ustar = ( PLR*vL/aL + vR/aR + _tdgm1*(PLR-1.) ) /
            ( PLR/aL + 1./aR );
    double pstar = 0.5*( WL.p()*pow(1.+_gm1d2/aL*(vL-ustar), _tgdgm1) +
                         WR.p()*pow(1.+_gm1d2/aR*(ustar-vR), _tgdgm1) );

    double vhalf;
    if(ustar < 0.){
        Whalf = WR;
        double pdpR = pstar/WR.p();
        // always a rarefaction wave, that's the approximation
        double SHR = vR + aR;
        if(SHR > 0){
            double STR = ustar + aR*pow(pdpR, _gm1d2g);
            if(STR <= 0){
                Whalf.set_rho(WR.rho()*pow(_tdgp1 - _gm1dgp1/aR*vR,
                                           _tdgm1));
                vhalf = _tdgp1*(-aR + _gm1d2*vR) - vR;
                Whalf.set_p(WR.p()*pow(_tdgp1 - _gm1dgp1/aR*vR, _tgdgm1));
            } else {
                Whalf.set_rho(WR.rho()*pow(pdpR, _ginv));
                Whalf.set_p(pstar);
                vhalf = ustar - vR;
            }
        } else {
            // Whalf = WR
            vhalf = 0.;
        }
    } else {
        Whalf = WL;
        double pdpL = pstar/WL.p();
        // rarefaction wave
        double SHL = vL - aL;
        if(SHL < 0){
            double STL = ustar - aL*pow(pdpL, _gm1d2g);
            if(STL > 0){
                Whalf.set_rho(WL.rho()*pow(_tdgp1 + _gm1dgp1/aL*vL,
                                           _tdgm1));
                vhalf = _tdgp1*(aL + _gm1d2*vL) - vL;
                Whalf.set_p(WL.p()*pow(_tdgp1 + _gm1dgp1/aL*vL, _tgdgm1));
            } else {
                Whalf.set_rho(WL.rho()*pow(pdpL, _ginv));
                vhalf = ustar - vL;
                Whalf.set_p(pstar);
            }
        } else {
            // Whalf = WL
            vhalf = 0.;
        }
    }

    Whalf[1] += vhalf*n[0];
    Whalf[2] += vhalf*n[1];
#if ndim_==3
    Whalf[3] += vhalf*n[2];
#endif

    _timer.stop();
    return Whalf;
}

/**
 * @brief Get the soundspeed corresponding to the given primitive variables
 *
 * @param W Primitive variables StateVector
 * @return Soundspeed
 */
double TRRSSolver::get_soundspeed(const StateVector &W){
    return sqrt(_gamma*W.p()/W.rho());
}

/**
 * @brief Test the Riemann solver on a set of problems with known solutions
 *
 * Since this is an approximate solver, it will fail on most of these problems.
 */
void TRRSSolver::test(){
    cout << "Testing the Riemann solver" << endl;
    double rhoL[6] = {1., 1., 1., 1., 5.99924, 1.};
    double rhoR[6] = {0.125, 1., 1., 1., 5.99242, 1.};
    double uL[6] = {0., -2., 0., 0., 19.5975, -1.};
    double uR[6] = {0., 2., 0., 0., -6.19633, 1.};
    double pL[6] = {1., 0.4, 1000., 0.01, 460.894, 1.e-6};
    double pR[6] = {0.1, 0.4, 0.01, 100., 46.0950, 1.0005e-6};
    double rhosol[6] = {0.47969, 0.00617903, 0.615719, 0.61577, 12.743, 0.};
    double usol[6] = {0.841194, 0., 18.2812, -5.78011, 8.56045, 0.};
    double psol[6] = {0.293945, 8.32249e-05, 445.626, 44.5687, 1841.82, 0.};
    for(unsigned int i = 0; i < 6; i++){
        cout << "Test problem " << i+1 << endl;
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
 * variables and volume
 *
 * @param volume Volume of the cell
 * @param W Primitive variables StateVector of the cell
 * @return Conserved variables StateVector
 */
StateVector TRRSSolver::get_Q(double volume, const StateVector &W){
    StateVector Q;
    Q.set_m(W.rho()*volume);
    Q.set_px(Q.m()*W.vx());
    Q.set_py(Q.m()*W.vy());
#if ndim_==3
    Q.set_pz(Q.m()*W.vz());
    Q.set_e((0.5*W.rho()*(W.vx()*W.vx() + W.vy()*W.vy() + W.vz()*W.vz()) +
             W.p()/(_gamma-1.))*volume);
#else
    Q.set_e((0.5*W.rho()*(W.vx()*W.vx() + W.vy()*W.vy()) +
             W.p()/(_gamma-1.))*volume);
#endif
    Q.set_paq(Q.m()*W.p()/pow(W.rho(), _gamma));
    return Q;
}

/**
 * @brief Get the primitive variables corresponding to the given conserved
 * variables and volume
 *
 * @param volume Volume of the cell
 * @param Q Conserved variables StateVector of the cell
 * @param use_energy Flag indicating if the pressure should be calculated using
 * the entropy (false) or the energy (true) formalism
 * @return Primitive variables StateVector
 */
StateVector TRRSSolver::get_W(double volume, StateVector &Q, bool use_energy){
    StateVector W;
    W.set_rho(Q.m()/volume);
    W.set_vx(Q.px()/Q.m());
    W.set_vy(Q.py()/Q.m());
#if ndim_==3
    W.set_vz(Q.pz()/Q.m());
    W.set_p((_gamma-1.)*(Q.e() - 0.5*(Q.px()*Q.px() + Q.py()*Q.py() +
                                      Q.pz()*Q.pz())/Q.m())/volume);
#else
    W.set_p((_gamma-1.)*(Q.e() - 0.5*(Q.px()*Q.px() +
                                      Q.py()*Q.py())/Q.m())/volume);
#endif
    if(use_energy){
        // reset entropy
        Q.set_paq(Q.m()*W.p()/pow(W.rho(), _gamma));
        W.set_paq(Q.paq()/Q.m());
    } else {
        W.set_paq(Q.paq()/Q.m());
        W.set_p(W.paq()*pow(W.rho(), _gamma));
        Q.set_e(get_energy(W)*volume);
    }
    return W;
}

/**
 * @brief Get the given component of the flux corresponding to the given
 * primitive variables
 *
 * @param v Velocity of the interface
 * @param index Requested component of the flux
 * @param W Primitive variables at the interface
 * @return Flux StateVector
 */
StateVector TRRSSolver::get_flux(const Vec &v, unsigned int index,
                                 const StateVector &W){
    StateVector F;
    F[0] = W.rho()*(W[1+index]-v[index]);
    F[1] = W.rho()*(W[1+index]-v[index])*W.vx();
    F[2] = W.rho()*(W[1+index]-v[index])*W.vy();
#if ndim_==3
    F[3] = W.rho()*(W[1+index]-v[index])*W.vz();
#endif
    F[1+index] += W.p();
    F[ndim_+1] = (W[1+index]-v[index])*get_energy(W) + W.p()*W[1+index];
    F.set_paq(W.rho()*(W[1+index]-v[index])*W.paq());
    return F;
}

/**
 * @brief Get the adiabatic index
 *
 * @return Adiabatic index
 */
double TRRSSolver::get_gamma(){
    return _gamma;
}

/**
 * @brief Get number of Riemann solver evaluations
 *
 * @return Number of Riemann solver evaluations
 */
unsigned long TRRSSolver::get_neval(){
    return _counter;
}

/**
 * @brief Dump the solver to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void TRRSSolver::dump(RestartFile &rfile){
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

    rfile.write(_counter);

    LOGS("TRRSSolver dumped");
}

/**
 * @brief Restart constructor. Initialize the solver from the given RestartFile
 *
 * @param rfile RestartFile to read from
 */
TRRSSolver::TRRSSolver(RestartFile &rfile) : _timer(rfile){
    rfile.read(_gamma);
    rfile.read(_gp1d2g);
    rfile.read(_gm1d2g);
    rfile.read(_gm1dgp1);
    rfile.read(_tdgp1);
    rfile.read(_tdgm1);
    rfile.read(_gm1d2);
    rfile.read(_tgdgm1);
    rfile.read(_ginv);

    rfile.read(_counter);

    LOGS("TRRSSolver restarted");
}
