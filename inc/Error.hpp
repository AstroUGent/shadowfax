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
 * @file Error.hpp
 *
 * @brief Custom abort macro that prints out useful information to the stderr
 *
 * Based on a similar macro in SWIFT.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ERROR_HPP
#define ERROR_HPP

#include <iostream>
#include <MPIGlobal.hpp>

#define my_exit() { \
    std::cerr << MPIGlobal::rank << ": " << __FILE__ << "::" << __FUNCTION__ \
    << ":" << __LINE__ << std::endl; MPI_Abort(MPI_COMM_WORLD, -1); }

#endif // ERROR_HPP
