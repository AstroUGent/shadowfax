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
 * @file Tenet.cpp
 *
 * @brief Implementation of Kevin Schaal's Tenet: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Tenet.hpp"
#include "Legendre.hpp"                 // for Legendre
#include "MPIMethods.hpp"               // for MyMPI_Finalize, MyMPI_Init
#include "RungeKutta.hpp"               // for RungeKutta
#include "StateVector.hpp"              // for StateVector, operator*
#include "TenetBasisFunctions.hpp"      // for TenetBasisFunctions, etc
#include "TenetCell.hpp"                // for TenetCell
#include "TenetCellWeights.hpp"         // for TenetCellWeights
#include "TenetGrid.hpp"                // for TenetGrid::iterator, etc
#include "TenetMinModSlopeLimiter.hpp"  // for TenetMinModSlopeLimiter
#include "TenetRungeKuttaFlux.hpp"      // for TenetRungeKuttaFlux
#include "TenetSystemState.hpp"         // for TenetSystemState, operator*, etc
#include "Vec.hpp"                      // for Vec, operator-, operator*, etc
#include "riemann/RiemannSolver.hpp"    // for RiemannSolver
#include "riemann/RiemannSolverFactory.hpp"  // for RiemannSolverFactory
#include "utilities/Cuboid.hpp"              // for Cuboid
#include <algorithm>                         // for min
#include <iostream>  // for operator<<, basic_ostream, etc
#include <math.h>    // for fabs, sqrt
#include <stdlib.h>  // for atoi
#include <vector>    // for vector
using namespace std;

/**
 * @brief Main program
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
Tenet::Tenet(int argc, char** argv) {
    unsigned char order = 2;
    if(argc > 1) {
        order = atoi(argv[1]);
    }

    cout << "Using order " << (unsigned int)order << " polynomials" << endl;

    Legendre legendre(order);
    TenetBasisFunctions basis(order);

    std::vector<double> zeros = legendre.get_zeros();
    std::vector<double> weights = legendre.get_weights();

#if ndim_ == 3
    unsigned int num1D = 45;
    unsigned int numcell = num1D * num1D * num1D;
    double cellwidth = 1. / num1D;
#else
    unsigned int num1D = 40;
    unsigned int numcell = num1D * num1D;
    double cellwidth = 1. / num1D;
#endif
    TenetSystemState cellweights(numcell);
    Vec anchor;
#if ndim_ == 3
    Vec sides(1., 1., 1.);
    unsigned int dims[3] = {num1D, num1D, num1D};
#else
    Vec sides(1., 1.);
    unsigned int dims[2] = {num1D, num1D};
#endif
    Cuboid box(anchor, sides);
    TenetGrid grid(dims, box);
    numcell = grid.size();

#if ndim_ == 3
    Vec center(0.5, 0.5, 0.5);
#else
    Vec center(0.5, 0.5);
#endif
    for(auto it = grid.begin(); it != grid.end(); ++it) {
        unsigned int i = it.index();
        cellweights[i] = basis.get_empty_weights();
        for(auto bit = basis.begin(); bit != basis.end(); ++bit) {
            for(auto cit = basis.cell_quadrature_begin();
                cit != basis.cell_quadrature_end(); ++cit) {
                Vec ksi = cit.get_ksi();
                Vec p = 0.5 * cellwidth * ksi + it->get_center();
                StateVector uksi;
                if((p - center).norm2() < 0.0625) {
                    uksi.set_rho(1.);
                    uksi.set_p(1.5);
                } else {
                    uksi.set_rho(0.125);
                    uksi.set_p(0.15);
                }
                double phij = bit(ksi);
                cellweights[i][bit.index()] += uksi * phij * cit.get_weight();
            }
#if ndim_ == 3
            cellweights[i][bit.index()] *= 0.125;
#else
            cellweights[i][bit.index()] *= 0.25;
#endif
        }
    }

    double gamma = 5. / 3.;
    RiemannSolver* solver =
            RiemannSolverFactory::generate("HLLC", gamma, 1.e-8, 5.);
    Vec vface;
    // time integration
    double t = 0.;
    double dt_fix = 0.001;
    Vec cellcenter;
    TenetRungeKuttaFlux<TenetSystemState> flux(grid, legendre, basis, solver);
    RungeKutta<TenetSystemState> integrator(order + 1, flux, cellweights);
    TenetMinModSlopeLimiter limiter(grid);
    double CFL = 0.2;
    while(t < dt_fix) {
        cout << "t = " << t << endl;
        double dt = dt_fix;
        double dt_min = dt_fix;
        for(auto it = grid.begin(); it != grid.end(); ++it) {
            StateVector ucell = it->get_variables(
                    cellcenter, cellweights[it.index()], basis);
            StateVector Wcell;
#if ndim_ == 3
            Wcell[0] = ucell[0];
            Wcell[1] = ucell[1] / ucell[0];
            Wcell[2] = ucell[2] / ucell[0];
            Wcell[3] = ucell[3] / ucell[0];
            Wcell[4] = (gamma - 1.) *
                       (ucell[4] -
                        0.5 * (ucell[1] * ucell[1] + ucell[2] * ucell[2] +
                               ucell[3] * ucell[3]) /
                                ucell[0]);
#else
            Wcell[0] = ucell[0];
            Wcell[1] = ucell[1] / ucell[0];
            Wcell[2] = ucell[2] / ucell[0];
            Wcell[3] = (gamma - 1.) *
                       (ucell[3] -
                        0.5 * (ucell[1] * ucell[1] + ucell[2] * ucell[2]) /
                                ucell[0]);
#endif
            double c_s = sqrt(gamma * Wcell.p() / Wcell.rho());
            double vsum = fabs(Wcell[1]) + fabs(Wcell[2]);
#if ndim_ == 3
            vsum += fabs(Wcell[3]);
#endif
            vsum += ndim_ * c_s;
            vsum /= cellwidth;
            double dt_cell = CFL / (2. * order + 1.) / vsum;
            dt_min = std::min(dt_cell, dt_min);
        }
        while(dt > dt_min) {
            dt *= 0.5;
        }
        cout << "dt = " << dt << endl;

        limiter.limit(cellweights);

        cellweights = integrator.integrate(dt);
        t += dt;
    }

    ofstream plot("plot.txt");
    for(auto it = grid.begin(); it != grid.end(); ++it) {
        unsigned int i = it.index();
        StateVector ucell =
                it->get_variables(cellcenter, cellweights[i], basis);
        plot << it->get_center().x() << "\t" << it->get_center().y();
#if ndim_ == 3
        plot << "\t" << it->get_center().z();
#endif
        plot << "\t" << (it->get_center() - center).norm() << "\t"
             << ucell.rho() << endl;
    }

    delete solver;
}

/**
 * @brief Entrance point for the Tenet program
 *
 * Calls the Tenet constructor.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    MyMPI_Init(&argc, &argv);

    Tenet(argc, argv);

    return MyMPI_Finalize();
}
