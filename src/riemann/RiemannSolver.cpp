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
 * @file RiemannSolver.cpp
 *
 * @brief RiemannSolver exceptions
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "RiemannSolver.hpp"
#include "ProgramLog.hpp"

/**
 * @brief Constructor
 *
 * @param W Derived primitive variables
 * @param Q Conserved variables
 * @param V Volume
 */
PrimitiveVariablesException::PrimitiveVariablesException(StateVector W,
                                                         StateVector Q,
                                                         double V) {
    _W = W;
    _Q = Q;
    _V = V;
}

/**
 * @brief Human readable error message for this Exception
 * @return Human readable C-string
 */
const char* PrimitiveVariablesException::what() const throw() {
    cerr << "Error! Negative density/pressure after calculation of primitive "
            "quantities!"
         << endl;
    cerr << "Q: " << _Q.m() << "\t" << _Q.px() << "\t" << _Q.py() << "\t";
#if ndim_ == 3
    cerr << _Q.pz() << "\t";
#endif
    cerr << _Q.e() << endl;
    cerr << "W: " << _W.rho() << "\t" << _W.vx() << "\t" << _W.vy() << "\t";
#if ndim_ == 3
    cerr << _W.vz() << "\t";
#endif
    cerr << _W.p() << endl;
    cerr << "V: " << _V << endl;
    return "PrimitiveVariablesException";
}

/**
 * @brief Constructor
 *
 * @param W Primitive variables
 * @param Q Derived conserved variables
 * @param V Volume
 */
ConservedVariablesException::ConservedVariablesException(StateVector W,
                                                         StateVector Q,
                                                         double V) {
    _W = W;
    _Q = Q;
    _V = V;

    LOGE("ConservedVariablesException raised!");
}

/**
 * @brief Human readable error message
 * @return Human readable C-string
 */
const char* ConservedVariablesException::what() const throw() {
    cerr << "Error! Negative mass/energy after calculation of conserved "
            "quantities!"
         << endl;
    cerr << "W: " << _W.rho() << "\t" << _W.vx() << "\t" << _W.vy() << "\t";
#if ndim_ == 3
    cerr << _W.vz() << "\t";
#endif
    cerr << _W.p() << endl;
    cerr << "Q: " << _Q.m() << "\t" << _Q.px() << "\t" << _Q.py() << "\t";
#if ndim_ == 3
    cerr << _Q.pz() << "\t";
#endif
    cerr << _Q.e() << endl;
    cerr << "V: " << _V << endl;
    return "ConservedVariablesException";
}

/**
 * @brief Constructor
 *
 * @param WL Left StateVector
 * @param WR Right StateVector
 * @param Wstar Solution of the Riemann problem
 */
RiemannSolverException::RiemannSolverException(StateVector WL, StateVector WR,
                                               StateVector Wstar) {
    _WL = WL;
    _WR = WR;
    _Wstar = Wstar;

    LOGE("RiemannSolverException raised!");
}

/**
 * @brief Human readable error message
 * @return Human readable C-string
 */
const char* RiemannSolverException::what() const throw() {
    cerr << "Error! Unacceptable solution to Riemann problem!" << endl;
    cerr << "WL: " << _WL.rho() << "\t" << _WL.vx() << "\t" << _WL.vy() << "\t";
#if ndim_ == 3
    cerr << _WL.vz() << "\t";
#endif
    cerr << _WL.p() << endl;
    cerr << "WR: " << _WR.rho() << "\t" << _WR.vx() << "\t" << _WR.vy() << "\t";
#if ndim_ == 3
    cerr << _WR.vz() << "\t";
#endif
    cerr << _WR.p() << endl;
    cerr << "Wstar: " << _Wstar.rho() << "\t" << _Wstar.vx() << "\t"
         << _Wstar.vy() << "\t";
#if ndim_ == 3
    cerr << _Wstar.vz() << "\t";
#endif
    cerr << _Wstar.p() << endl;
    return "RiemannSolverException";
}
