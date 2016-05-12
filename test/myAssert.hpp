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
 * @file myAssert.hpp
 *
 * @brief Custom macro to assert if floating point values are equal
 *
 * Useful for checking if the results of analytic test problems are reproduced.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#define my_assert(expression, message)                                   \
    if(!(expression)) {                                                  \
        std::cerr << message << std::endl;                               \
        std::cerr << __FILE__ << "::" << __FUNCTION__ << ":" << __LINE__ \
                  << std::endl;                                          \
        abort();                                                         \
    }

#define assert_values_equal(a, b, message)               \
    if(b != 0.) {                                        \
        if(fabs(a - b) / fabs(b) > 1.e-4) {              \
            std::cerr << message << std::endl;           \
            std::cerr << a << " =/= " << b << std::endl; \
            abort();                                     \
        }                                                \
    } else {                                             \
        if(fabs(a - b) > 1.e-4) {                        \
            std::cerr << message << std::endl;           \
            std::cerr << a << " =/= " << b << std::endl; \
            abort();                                     \
        }                                                \
    }
