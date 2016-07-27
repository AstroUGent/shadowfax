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
 * @file HelperFunctions.hpp
 *
 * @brief Various auxiliary functions: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HELPERFUNCIONS_HPP
#define HELPERFUNCIONS_HPP

#include <cstdlib>
#include <string>

/**
 * @brief General auxiliary functions
 */
namespace HelperFunctions {
std::string human_readable_bytes(unsigned int bytes);
unsigned int machine_readable_bytes(std::string input);

std::string human_readable_counter(unsigned long counter);

std::string human_readable_long(unsigned long num);

std::string make_hdf5_file(std::string name);

std::string get_absolute_path(std::string relative_path);

/**
 * @brief Generate a random double precision floating point value in the range
 * [0, 1[
 *
 * @return A uniform random double
 */
inline double rand_double() {
    return ((double)rand()) / ((double)RAND_MAX);
}

/**
 * @brief Generate a random value from the given range
 *
 * @param min Lower limit of the range
 * @param max Upper limit of the range
 * @return Random value in the range
 */
inline double rand_value(double min, double max) {
    return min + rand_double() * (max - min);
}

/**
 * @brief Get the sign of the argument or 0 when the argument is zero
 *
 * Code taken from
 * http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-
 * signum-sgn-in-c-c
 *
 * @param val Value
 * @return -1, 0 or 1, depending on the sign of the value
 */
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}
}

#endif  // HELPERFUNCIONS_HPP
