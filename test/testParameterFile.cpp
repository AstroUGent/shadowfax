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
 * @file testParameterFile.cpp
 *
 * @brief Unit test for the ParameterFile class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "MPIMethods.hpp"
#include "ParameterFile.hpp"
#include "io/Unit.hpp"
#include "myAssert.hpp"

/**
 * @brief Unit test for the ParameterFile class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char** argv) {
    std::string filenames[2] = {"test.ini", "paramtest.yml"};
    for(unsigned int i = 0; i < 2; ++i) {
        ParameterFile params(filenames[i]);

        my_assert(params.get_parameter<double>("Test.test_parameter_float",
                                               0.) == 42.,
                  "Floating point parameter reading fails!");
        my_assert(params.get_parameter<int>("Test.test_parameter_int", 0) == 42,
                  "Integer parameter reading fails!");
        my_assert(
                params.get_parameter<std::string>("Test.test_parameter_string",
                                                  "None") == "test_string",
                "String parameter reading fails!");
        my_assert(params.get_parameter<bool>("Test.test_parameter_bool",
                                             true) == false,
                  "Boolean parameter reading fails!");

        Unit unit_length("length", "m", 1.);
        double value = params.get_quantity("Test.test_physical_value",
                                           unit_length, "0. m");
        my_assert(value == 1., "Physical value parameter reading fails!");
    }

    return 0;
}
