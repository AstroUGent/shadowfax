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
 * @file RiemannSolver.hpp
 *
 * @brief General definitions for Riemann solver classes: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef RIEMANNSOLVER_HPP
#define RIEMANNSOLVER_HPP

#include "StateVector.hpp"
#include "../src/utilities/Timer.hpp"

class RestartFile;

/**
 * @brief General interface for Riemann solvers
 */
class RiemannSolver{
public:
    // if we would not declare this destructor, it would not be virtual
    // and the destructors of deriving classes would not be called
    virtual ~RiemannSolver(){}

    /**
     * @brief Get a std::string describing the specific implementation of
     * RiemannSolver
     *
     * @return An implementation specific identifying tag for the Solver
     */
    virtual std::string tag()=0;

    /**
     * @brief Solve the Riemann problem with given left and right states
     *
     * @param WL Left StateVector
     * @param WR Right StateVector
     * @param n Normal vector to the interface
     * @param mach Mach-number of the strongest shock for the left state
     * @param reflective Flag indicating if the right state should be a
     * reflective copy of the left state
     * @return Solution of the Riemann problem
     */
    virtual StateVector solve(StateVector& WL, StateVector& WR, Vec& n,
                              double& mach, bool reflective = false)=0;

    /**
     * @brief Get the hydrodynamical soundspeed for the given StateVector
     *
     * @param W Given StateVector
     * @return Soundspeed of the fluid with the given state
     */
    virtual double get_soundspeed(const StateVector& W)=0;

    /**
     * @brief Test the Riemann solver on a set of problems with known solutions
     */
    virtual void test()=0;

    /**
     * @brief Calculate conserved quantities from the given primitive variables
     * and given volume
     *
     * @param volume Volume of a cell
     * @param W StateVector in a cell
     * @return Conserved quantities StateVector
     */
    virtual StateVector get_Q(double volume, const StateVector& W)=0;

    /**
     * @brief Calculate primitive variables from the given conserved variables
     * and volume
     *
     * @param volume Volume of a cell
     * @param Q Conserved quantities StateVector of a cell
     * @param use_energy Flag indicating if we should use the energy or entropy
     * formulation to convert total energy to pressure
     * @return Primitive variables StateVector
     */
    virtual StateVector get_W(double volume, StateVector& Q,
                              bool use_energy = true)=0;

    /**
     * @brief Get the given component of the flux matrix for the given state and
     * given interface velocity
     *
     * @param v Velocity of the interface
     * @param index unsigned integer index of the requested flux component
     * @param W StateVector at the interface
     * @return StateVector of fluxes through the interface
     */
    virtual StateVector get_flux(const Vec& v, unsigned int index,
                                 const StateVector& W)=0;

    /**
     * @brief Get the adiabatic index used in this Solver
     *
     * @return The adiabatic index
     */
    virtual double get_gamma()=0;

    /**
     * @brief Dump the Solver to the given RestartFile
     *
     * @param rfile RestartFile to write to
     */
    virtual void dump(RestartFile &rfile)=0;

    /**
     * @brief Get the number of Riemann solver evaluations
     *
     * @return Number of Riemann solver evaluations
     */
    virtual unsigned long get_neval()=0;
};

/**
 * @brief Exception thrown when something goes wrong during the update of the
 * primitive variables
 */
class PrimitiveVariablesException : public std::exception{
private:
    /*! \brief Derived primitive variables */
    StateVector _W;

    /*! \brief Conserved variables */
    StateVector _Q;

    /*! \brief Volume */
    double _V;

public:
    PrimitiveVariablesException(StateVector W, StateVector Q, double V);

    virtual const char* what() const throw();
};

/**
 * @brief Exception thrown when something goes wrong during the conversion of
 * primitive variables to convserved variables
 */
class ConservedVariablesException : public std::exception{
private:
    /*! \brief Primitive variables */
    StateVector _W;

    /*! \brief Derived conserved variables */
    StateVector _Q;

    /*! \brief Volume */
    double _V;

public:
    ConservedVariablesException(StateVector W, StateVector Q, double V);

    virtual const char* what() const throw();
};

/**
 * @brief Exception thrown when the RiemannSolver outputs unphysical results
 */
class RiemannSolverException : public std::exception{
private:
    /*! \brief Left StateVector */
    StateVector _WL;

    /*! \brief Right StateVector */
    StateVector _WR;

    /*! \brief Riemann solver solution */
    StateVector _Wstar;

public:
    RiemannSolverException(StateVector WL, StateVector WR, StateVector Wstar);

    virtual const char* what() const throw();
};

#endif // RIEMANNSOLVER_HPP
