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

#include "AdvectionRiemannSolver.hpp"   // for AdvectionRiemannSolver
#include "ApproximateSolver.hpp"        // for TRRSSolver
#include "Error.hpp"                    // for my_exit
#include "ExactRiemannSolver.hpp"       // for ExactRiemannSolver
#include "HLLCRiemannSolver.hpp"        // for HLLCRiemannSolver
#include "NoPressureRiemannSolver.hpp"  // for NoPressureRiemannSolver
#include "ParameterFile.hpp"
#include "RestartFile.hpp"
#include "RiemannSolver.hpp"  // for RiemannSolver
#include <cstddef>            // for NULL
#include <string>             // for operator==, basic_string, etc

#define RIEMANNSOLVERFACTORY_DEFAULT_TYPE "Exact"
#define RIEMANNSOLVERFACTORY_DEFAULT_GAMMA (5. / 3.)
#define RIEMANNSOLVERFACTORY_DEFAULT_TOLERANCE 1.e-8
#define RIEMANNSOLVERFACTORY_DEFAULT_CUTOFF -5.

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
     * @brief Load a RiemannSolver from the given RestartFile
     *
     * We first read the tag from the RestartFile and then call the restart
     * constructor for the corresponding RiemannSolver implementation.
     *
     * @param rfile RestartFile to read from
     * @return RiemannSolver instance
     */
    static RiemannSolver* load(RestartFile& rfile) {
        std::string tag;
        rfile.read(tag);
        if(tag == "EXAC") {
            return new ExactRiemannSolver(rfile);
        }
        if(tag == "HLLC") {
            return new HLLCRiemannSolver(rfile);
        }
        if(tag == "TRRS") {
            return new TRRSSolver(rfile);
        }
        if(tag == "ADVE") {
            return new AdvectionRiemannSolver(rfile);
        }
        if(tag == "NOPR") {
            return new NoPressureRiemannSolver(rfile);
        }
        return NULL;
    }

    /**
     * @brief Dump the given RiemannSolver to the given RestartFile
     *
     * We first write an identifying tag to the RestartFile and then call the
     * dump method of the RiemannSolver itself.
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
     * @brief Generate a RiemannSolver of the given type
     *
     * @param name Type of RiemannSolver to generate
     * @param gamma Adiabatic index of the gas
     * @param tolerance Tolerance value used to decide when the iterative
     * method to find the pressure is converged (exact solver only)
     * @param cutoff Cutoff value to distinguish between Newton-Raphson and
     * Brent's method in the iterative pressure finding procedure (exact solver
     * only)
     * @return A pointer to a RiemannSolver instance
     */
    static RiemannSolver* generate(std::string name, double gamma,
                                   double tolerance, double cutoff) {
        if(name == "Exact") {
            return new ExactRiemannSolver(gamma, tolerance, cutoff);
        }
        if(name == "HLLC") {
            return new HLLCRiemannSolver(gamma);
        }
        if(name == "TRRS") {
            return new TRRSSolver(gamma);
        }
        if(name == "NoPressure") {
            return new NoPressureRiemannSolver();
        }
        if(name == "Advection") {
            return new AdvectionRiemannSolver(gamma);
        }
        std::cerr << "Error! Unknown Riemann solver type: " << name << "!"
                  << std::endl;
        my_exit();
        return NULL;
    }

    /**
     * @brief Generate a RiemannSolver based on the given parameters
     *
     * @param parameters ParameterFile containing the desired parameter values
     * @return A pointer to a RiemannSolver instance
     */
    static RiemannSolver* generate(ParameterFile* parameters) {
        std::string name = parameters->get_parameter<std::string>(
                "RiemannSolver.Type", RIEMANNSOLVERFACTORY_DEFAULT_TYPE);
        double gamma = parameters->get_parameter<double>(
                "Hydro.Gamma", RIEMANNSOLVERFACTORY_DEFAULT_GAMMA);
        double tolerance = parameters->get_parameter<double>(
                "RiemannSolver.Tolerance",
                RIEMANNSOLVERFACTORY_DEFAULT_TOLERANCE);
        double cutoff = parameters->get_parameter<double>(
                "RiemannSolver.CutOff", RIEMANNSOLVERFACTORY_DEFAULT_CUTOFF);
        return generate(name, gamma, tolerance, cutoff);
    }
};

#endif
