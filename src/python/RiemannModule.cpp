/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2016 Bert Vandenbroucke
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
 * @file RiemannModule.cpp
 *
 * @brief Expose the exact Riemann solver to Python
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "riemann/ApproximateSolver.hpp"
#include "riemann/ExactRiemannSolver.hpp"
#include "riemann/RiemannSolver.hpp"

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/python/module.hpp>
using namespace boost::python;

/**
 * @brief Wrapper around ExactRiemannSolver::solve
 *
 * @param solver Active RiemannSolver instance used to solve the Riemann problem
 * @param WL Python list containing the left state hydrodynamical variables
 * @param WR Python list containing the right state hydrodynamical variables
 * @return Python list containing the solution hydrodynamical variables
 */
template <typename Solver> list python_solve(Solver& solver, list WL, list WR) {

    double WLarr[ndim_ + 2];
    for(unsigned int i = 0; i < ndim_ + 2; i++) {
        WLarr[i] = extract<double>(WL[i]);
    }
    double WRarr[ndim_ + 2];
    for(unsigned int i = 0; i < ndim_ + 2; i++) {
        WRarr[i] = extract<double>(WR[i]);
    }

#if ndim_ == 3
    StateVector WLvec(WLarr[0], WLarr[1], WLarr[2], WLarr[3], WLarr[4]);
    StateVector WRvec(WRarr[0], WRarr[1], WRarr[2], WRarr[3], WRarr[4]);
#else
    StateVector WLvec(WLarr[0], WLarr[1], WLarr[2], WLarr[3]);
    StateVector WRvec(WRarr[0], WRarr[1], WRarr[2], WRarr[3]);
#endif
    Vec normal;
    normal[0] = 1.;
    double maxmach;
    StateVector Whalfvec = solver.solve(WLvec, WRvec, normal, maxmach);

    list Whalf;
    Whalf.append(Whalfvec[0]);
    Whalf.append(Whalfvec[1]);
    Whalf.append(Whalfvec[2]);
    Whalf.append(Whalfvec[3]);
#if ndim_ == 3
    Whalf.append(Whalfvec[4]);
#endif
    return Whalf;
}

#if ndim_ == 3
/*! @brief Register the libriemann Python module */
BOOST_PYTHON_MODULE(libpython_riemann3d)
#else
BOOST_PYTHON_MODULE(libpython_riemann2d)
#endif
{
    class_<ExactRiemannSolver>("ExactRiemannSolver", init<double>())
            .def("solve", &python_solve<ExactRiemannSolver>);

    class_<TRRSSolver>("TRRSSolver", init<double>())
            .def("solve", &python_solve<TRRSSolver>);
}
