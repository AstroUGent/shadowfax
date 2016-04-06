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
 * @file RiemannSolverFactory.hpp
 *
 * @brief Factory to load and dump Solver instances without having to worry
 * about their actual type
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef RIEMANNSOLVERFACTORY_HPP
#define RIEMANNSOLVERFACTORY_HPP

#include "ApproximateSolver.hpp"
#include "ExactRiemannSolver.hpp"
#include <istream>

#include "RestartFile.hpp"

/**
 * @brief Factory to create and dump Solver instances to a RestartFile
 *
 * Currently supports:
 *  - RiemannSolver
 *  - TRRSSolver
 */
class RiemannSolverFactory {
  public:
    /**
     * @brief Load a Solver from the given RestartFile
     *
     * We first read the tag from the RestartFile and then call the restart
     * constructor for the corresponding Solver implementation.
     *
     * @param rfile RestartFile to read from
     * @return Solver instance
     */
    static RiemannSolver* load(RestartFile& rfile) {
        std::string tag;
        rfile.read(tag);
        if(tag == "EXAC") {
            return new ExactRiemannSolver(rfile);
        }
        if(tag == "TRRS") {
            return new TRRSSolver(rfile);
        }
        return NULL;
    }

    /**
     * @brief Dump the given Solver to the given RestartFile
     *
     * We first write an identifying tag to the RestartFile and then call the
     * dump method of the solver itself.
     *
     * @param rfile RestartFile to write to
     * @param solver RiemannSolver to dump
     */
    static void dump(RestartFile& rfile, RiemannSolver* solver) {
        // it is important we make a local copy of the tag, since otherwise it
        // is not written to the restart file correctly!
        std::string tag = solver->tag();
        rfile.write(tag);
        solver->dump(rfile);
    }

    /**
     * @brief Generate a Riemann solver of the given type
     *
     * @param name Type of Riemann solver to generate
     * @param gamma Adiabatic index of the gas
     * @param tolerance Tolerance value used to decide when the iterative
     * method to find the pressure is converged (exact solver only)
     * @param cutoff Cutoff value to distinguish between Newton-Raphson and
     * Brent's method in the iterative pressure finding procedure (exact solver
     * only)
     * @return A pointer to a Solver instance
     */
    static RiemannSolver* generate(std::string name, double gamma,
                                   double tolerance, double cutoff) {
        if(name == "Exact") {
            return new ExactRiemannSolver(gamma, tolerance, cutoff);
        }
        if(name == "TRRS") {
            return new TRRSSolver(gamma);
        }
        std::cerr << "Error! Unknown Riemann solver type: " << name << "!"
                  << std::endl;
        my_exit();
        return NULL;
    }
};

#endif
