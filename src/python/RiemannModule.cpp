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
 * @file RiemannModule.cpp
 *
 * @brief Expose the exact Riemann solver to Python
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "riemann/ExactRiemannSolver.hpp"

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
using namespace boost::python;

#if ndim_ == 3
/*! @brief Register the libriemann Python module */
BOOST_PYTHON_MODULE(libpython_riemann3d)
#else
BOOST_PYTHON_MODULE(libpython_riemann2d)
#endif
{
#if ndim_ == 3
    class_<StateVector>("StateVector",
                        init<double, double, double, double, double>())
            .def("rho", &StateVector::rho)
            .def("vx", &StateVector::vx)
            .def("vy", &StateVector::vy)
            .def("vz", &StateVector::vz)
            .def("p", &StateVector::p)
            .def("set", static_cast<void (StateVector::*)(
                                double, double, double, double, double)>(
                                &StateVector::set));
#else
    class_<StateVector>("StateVector", init<double, double, double, double>())
            .def("rho", &StateVector::rho)
            .def("vx", &StateVector::vx)
            .def("vy", &StateVector::vy)
            .def("p", &StateVector::p)
            .def("set",
                 static_cast<void (StateVector::*)(double, double, double,
                                                   double)>(&StateVector::set));
#endif

    class_<ExactRiemannSolver>("ExactRiemannSolver", init<double>())
            .def("solve", &ExactRiemannSolver::solve);
}
