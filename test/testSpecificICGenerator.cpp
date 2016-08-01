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
 * @file testSpecificICGenerator.cpp
 *
 * @brief Unit test for the SpecificICGenerator class
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "SpecificICGenerator.hpp"
#include <iostream>

/**
 * @brief SpecificICGenerator test
 *
 * Tests the functionality of the SpecificICGenerator class
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {

    for(double x = 0.; x < 1.; x += 0.01) {
        std::cout << x << "\t" << SpecificICGenerator::gplummer(x) << std::endl;
    }

    return 0;
}
