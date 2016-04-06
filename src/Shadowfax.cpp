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
 * @file Shadowfax.cpp
 *
 * @brief Main simulation program
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "DelTess.hpp"
#include "ExArith.h"
#include "MPIMethods.hpp"
#include "Plotter.hpp"
#include "SidePrograms.hpp"
#include "Simulation.hpp"
#include <getopt.h>
#include <iostream>
using namespace std;

/**
 * @brief Sideprograms that can be called using the corresponding command line
 * argument
 */
enum SideProgramTypes {
    /*! \brief Plot sideprogram */
    SIDEPROGRAM_PLOT = 0,
    /*! \brief ICMaker sideprogram (no longer supported and replaced by an
     *  separate program) */
    SIDEPROGRAM_IC,
    /*! \brief Make a grid of the given points and return the cell areas (2D
     *  only) */
    SIDEPROGRAM_AREA,
    /*! \brief Make a grid of the given points and return the inverse volumes
     *  as a local density measure (3D only) */
    SIDEPROGRAM_DENSITY,
    /*! \brief Calculate masses for the given snapshot */
    SIDEPROGRAM_MASS,
    /*! \brief Calculate gravitational potential energies for the given
     *  snapshot */
    SIDEPROGRAM_EPOT,
    /*! \brief Sort the given snapshot in Hilbert space filling order */
    SIDEPROGRAM_SORT,
    /*! \brief Test the Delaunay construction algorithms */
    SIDEPROGRAM_DELTEST,
    /*! \brief Test the exact arithmetics routines */
    SIDEPROGRAM_EXARITH
};

/**
 * @brief Check if a command line argument calls a sideprogram and start the
 * corresponding program if applicable
 *
 * Supported sideprograms:
 *  - plot: make an image of the given snapshot file
 *  - ic: deprecated and replaced by the icmaker program (will still flag
 *    however)
 *  - area: calculate the areas of the 2D Voronoi cells for a given set of
 *    points
 *  - density: calculate the inverse volumes of the 3D Voronoi cells for a
 *    given set of points
 *  - mass: calculate masses by constructing the Voronoi grid for the given
 *    snapshot
 *  - potential: calculate gravitational potential energies for the given
 *    snapshot
 *  - sort: sort the given points in space filling Hilbert order
 *  - delaunaycheck: check the Delaunay construction algorithms and exit
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return 1 if a sideprogram was flagged and run, 0 if the main simulation
 * program should be started
 */
unsigned int bootstrap(int argc, char** argv) {
    if(argc < 2) {
        return 0;
    }

    static struct option long_options[] = {
            {"plot", no_argument, NULL, SIDEPROGRAM_PLOT},
            {"ic", no_argument, NULL, SIDEPROGRAM_IC},
            {"area", no_argument, NULL, SIDEPROGRAM_AREA},
            {"density", no_argument, NULL, SIDEPROGRAM_DENSITY},
            {"mass", no_argument, NULL, SIDEPROGRAM_MASS},
            {"epot", no_argument, NULL, SIDEPROGRAM_EPOT},
            {"sort", no_argument, NULL, SIDEPROGRAM_SORT},
            {"delaunaycheck", no_argument, NULL, SIDEPROGRAM_DELTEST},
            {"exact_arithmetics", no_argument, NULL, SIDEPROGRAM_EXARITH},
            {0, 0, 0, 0}};

    int c;
    opterr = 0;
    while((c = getopt_long(argc, argv, "", long_options, NULL)) != -1) {
        // the extra {} need to be present for this code to compile, see
        // http://stackoverflow.com/questions/2351936/create-an-object-in-
        // switch-case
        switch(c) {
            case SIDEPROGRAM_PLOT: {
                // plot snapshot
                Plotter pl(argc, argv);
                return 1;
            }
            case SIDEPROGRAM_IC: {
                // make ic
                cerr << "This functionality is no longer supported by the main "
                        "Shadowfax program. Use the program icmakerNd instead!"
                     << endl;
                return 1;
            }
            case SIDEPROGRAM_AREA: {
                AreaCalculator ac(argc, argv);
                return 1;
            }
            case SIDEPROGRAM_DENSITY: {
                // density code for Cloet-Osselaer et al. 2014
                DensityCalculator dc(argc, argv);
                return 1;
            }
            case SIDEPROGRAM_MASS: {
                MassCalculator mc(argc, argv);
                return 1;
            }
            case SIDEPROGRAM_EPOT: {
                PotentialCalculator pc(argc, argv);
                return 1;
            }
            case SIDEPROGRAM_SORT: {
                HilbertSorter hs(argc, argv);
                return 1;
            }
            case SIDEPROGRAM_DELTEST: {
                // check the delaunay algorithms
                DelTess delaunay(NULL, 0);
                delaunay.check_methods();
                return 1;
            }
            case SIDEPROGRAM_EXARITH: {
                ExactArithmetic::test_predicates();
                return 1;
            }
            case '?':
                return 0;
            default:
                return 0;
        }
    }
    return 0;
}

/**
 * @brief Entrance point of the main shadowfax simulaton program
 *
 * Checks if a sideprogram should be run by calling bootstrap(). If not, starts
 * the main simulation program by calling the Simulation constructor.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    MyMPI_Init(&argc, &argv);

    if(!bootstrap(argc, argv)) {
        // no specific action required. Run simulation
        Simulation sim(argc, argv);
    }

    return MyMPI_Finalize();
}

/**
 * \mainpage
 * \author Bert Vandenbroucke \n
 *         Vakgroep Fysica en Sterrenkunde \n
 *         Krijgslaan 281 \n
 *         9000 Gent \n
 *         Belgium \n
 *         bert.vandenbroucke@ugent.be
 *
 *
 * \section purpose Purpose of the program
 *
 * Shadowfax is a parallel moving mesh hydrodynamical integration code that can
 * be used to solve a wide range of hydrodynamical problems.
 *
 * The Euler equations of hydrodynamics are discretized on a co-moving voronoi
 * grid and are solved using a second order MUSCL-Hancock integration scheme.
 */
