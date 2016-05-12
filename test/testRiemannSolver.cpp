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
 * @file testRiemannSolver.cpp
 *
 * @brief Unit test for the Riemann solver
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "StateVector.hpp"                 // for StateVector, operator*
#include "Vec.hpp"                         // for Vec
#include "myAssert.hpp"                    // for assert_values_equal
#include "riemann/ExactRiemannSolver.hpp"  // for ExactRiemannSolver
#include "riemann/RiemannSolver.hpp"       // for RiemannSolver
#include <iostream>

/**
 * @brief Test the given RiemannSolver on a given problem and verify the result
 *
 * @param rhoL Left density
 * @param vxL Left velocity
 * @param pL Left pressure
 * @param rhoR Right density
 * @param vxR Right velocity
 * @param pR Right pressure
 * @param rho Theoretical resulting density
 * @param vx Theoretical resulting velocity
 * @param p Theoretical resulting pressure
 * @param solver RiemannSolver to test
 */
void test_RiemannSolver(double rhoL, double vxL, double pL, double rhoR,
                        double vxR, double pR, double rho, double vx, double p,
                        RiemannSolver& solver) {
    StateVector WL, WR;

    WL.set_rho(rhoL);
    WL.set_vx(vxL);
    WL.set_p(pL);

    WR.set_rho(rhoR);
    WR.set_vx(vxR);
    WR.set_p(pR);

    Vec n;
    n[0] = 1.;
    double mach;

    StateVector Whalf = solver.solve(WL, WR, n, mach);

    assert_values_equal(Whalf.rho(), rho, "Density incorrect!");
    assert_values_equal(Whalf.vx(), vx, "Velocity incorrect!");
    assert_values_equal(Whalf.p(), p, "Pressure incorrect!");
}

/**
 * @brief Test the given RiemannSolver on the given problem with given flux
 * solution
 *
 * @param rhoL Left density
 * @param vxL Left velocity
 * @param pL Left pressure
 * @param rhoR Right density
 * @param vxR Right velocity
 * @param pR Right pressure
 * @param Frho Theoretical resulting density flux
 * @param Fvx Theoretical resulting velocity flux
 * @param Fp Theoretical resulting pressure flux
 * @param solver RiemannSolver to test
 */
void test_RiemannFlux(double rhoL, double vxL, double pL, double rhoR,
                      double vxR, double pR, double Frho, double Fvx, double Fp,
                      RiemannSolver& solver) {
    StateVector WL, WR;

    WL.set_rho(rhoL);
    WL.set_vx(vxL);
    WL.set_p(pL);

    WR.set_rho(rhoR);
    WR.set_vx(vxR);
    WR.set_p(pR);

    Vec n;
    n[0] = 1.;
    Vec v;

    StateVector flux = solver.solve_for_flux(WL, WR, n, v);

    cerr << flux.rho() << " =?= " << Frho << endl;
    cerr << flux.vx() << " =?= " << Fvx << endl;
    cerr << flux.p() << " =?= " << Fp << endl;
    //    assert_values_equal(flux.rho(), Frho, "Density flux incorrect!");
    //    assert_values_equal(flux.vx(), Fvx, "Velocity flux incorrect!");
    //    assert_values_equal(flux.p(), Fp, "Pressure flux incorrect!");
}

/**
 * @brief Run a small 1D fixed grid Sod shock test and output result to stdout
 */
void test_Sod() {
    ExactRiemannSolver solver(5. / 3.);

    StateVector Qs[100];
    StateVector Ws[100];

    for(unsigned int i = 0; i < 100; i++) {
        if(i < 50) {
            Ws[i].set_rho(0.);
            Ws[i].set_vx(0.);
            Ws[i].set_p(0.);
        } else {
            Ws[i].set_rho(1.);
            Ws[i].set_vx(0.);
            Ws[i].set_p(1.);
        }
        Qs[i] = solver.get_Q(0.01, Ws[i]);
    }

#if ndim_ == 3
    Vec n(1., 0.25, 0.5);
#else
    Vec n(1., 0.5);
#endif
    // normalize n
    n /= n.norm();
    Vec v;
    //    v[0] = 0.1;
    Vec vzero;
    double dt = 0.1 / 64.;
#define PERIODIC
    for(unsigned int t = 0; t < 64; t++) {
        for(unsigned int i = 0; i < 100; i++) {
#ifdef PERIODIC
            StateVector WL = Ws[i];
            WL -= v;
            unsigned int ip = i + 1;
            if(ip == 100) {
                ip = 0;
            }
            StateVector WR = Ws[ip];
            WR -= v;
            StateVector flux = solver.solve_for_flux(WL, WR, n, v, false);
            Qs[i] -= dt * flux;
            Qs[ip] += dt * flux;
#else
            StateVector WL = Ws[i];
            if(!i) {
                WL -= vzero;
                StateVector flux =
                        solver.solve_for_flux(WL, WL, n, vzero, true);
                Qs[i] += dt * flux;
            }
            if(i == 99) {
                WL -= vzero;
                StateVector flux =
                        solver.solve_for_flux(WL, WL, n, vzero, true);
                Qs[i] -= dt * flux;
            } else {
                WL -= v;
                StateVector WR = Ws[i + 1];
                WR -= v;
                StateVector flux = solver.solve_for_flux(WL, WR, n, v, false);
                Qs[i] -= dt * flux;
                Qs[i + 1] += dt * flux;
            }
#endif
        }
        for(unsigned int i = 0; i < 100; i++) {
            Ws[i] = solver.get_W(0.01, Qs[i], true);
        }
    }
    for(unsigned int i = 0; i < 100; i++) {
        cout << (i + 0.5) * 0.01 << "\t" << Ws[i].rho() << endl;
    }
}

/**
 * @brief Test the Riemann solver(s)
 *
 * Currently only tests the HLLCRiemannSolver, without really testing the
 * results
 *
 * Also runs the small fixed grid 1D sod test
 *
 * @param argc Number of command line arguments (ignored)
 * @param argv Command line arguments (ignored)
 * @return 0 on succes. Aborts otherwise
 */
int main(int argc, char** argv) {
    //    ExactRiemannSolver solver;
    //    HLLCRiemannSolver solver(5. / 3.);
    ExactRiemannSolver solver;
    StateVector WL(0.125, 0., 0., 0., 0.1);
    StateVector WR(0.5625, 0., 0., 0., 0.55);
    Vec n(1., 0., 0.);
    Vec v;
    StateVector flux = solver.solve_for_flux(WL, WR, n, v);
    cout << flux[0] << " " << flux[1] << " " << flux[2] << " " << flux[3] << " "
         << flux[4] << endl;

    //    test_RiemannSolver(1., 0., 1., 0.125, 0., 0.1, 0.47969, 0.841194,
    // 0.293945,
    //                       solver);
    //    test_RiemannSolver(1., -2., 0.4, 1., 2., 0.4, 0.00617903, 0.,
    // 8.32249e-05,
    //                       solver);
    //    test_RiemannSolver(1., 0., 1000., 1., 0., 0.01, 0.615719, 18.2812,
    // 445.626,
    //                       solver);
    //    test_RiemannSolver(1., 0., 0.01, 1., 0., 100., 0.61577, -5.78011,
    // 44.5687,
    //                       solver);
    //    test_RiemannSolver(5.99924, 19.5975, 460.894, 5.99242, -6.19633,
    // 46.0950,
    //                       12.743, 8.56045, 1841.82, solver);
    //    test_RiemannSolver(1., -1., 1.e-6, 1., 1., 1.0005e-6, 0., 0., 0.,
    // solver);

    //    test_RiemannFlux(1., 0., 1., 0.125, 0., 0.1, 0.403512, 0.633377,
    //    0.760924,
    //                     solver);
    //    test_RiemannFlux(1., -2., 0.4, 1., 2., 0.4, 0., 8.32249e-05, 0.,
    //    solver);
    //    test_RiemannFlux(1., 0., 1000., 1., 0., 0.01, 11.2561, 651.4, 22247.3,
    //                     solver);
    //    test_RiemannFlux(1., 0., 0.01, 1., 0., 100., -3.55922, 65.1414,
    //    -703.485,
    //                     solver);
    //    test_RiemannFlux(5.99924, 19.5975, 460.894, 5.99242, -6.19633,
    //    46.0950,
    //                     109.086, 2775.65, 43413.9, solver);
    //    test_RiemannFlux(1., -1., 1.e-6, 1., 1., 1.0005e-6, 0., 0., 0.,
    //    solver);

    //    test_Sod();

    return 0;
}
