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
 * @file testUnit.cpp
 *
 * @brief Unit test for the Unit class
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "io/Unit.hpp"
#include "io/UnitDefinitions.hpp"
#include "myAssert.hpp"
#include <iostream>
#include <string>

/**
 * @brief Unit test
 *
 * Tests the functionality of the Unit class
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    Unit unit_length("length", "m", 1.);
    Unit unit_velocity("length/time", "m/s", 1.);
    Unit unit_time("time", "s", 1.);

    Unit unit_time2 = unit_length / unit_velocity;

    if(unit_time.get_quantity().compare(unit_time2.get_quantity())) {
        std::cerr << unit_time.get_quantity() << "\t"
                  << unit_time2.get_quantity() << std::endl;
        std::cerr << unit_time.get_quantity().compare(unit_time2.get_quantity())
                  << std::endl;
        std::cerr << "Units are not compatible!" << std::endl;
        exit(1);
    }

    Unit unit_length2 = unit_velocity * unit_time;

    if(unit_length.get_quantity().compare(unit_length2.get_quantity())) {
        std::cerr << unit_length.get_quantity() << "\t"
                  << unit_length2.get_quantity() << std::endl;
        std::cerr << unit_length.get_quantity().compare(
                             unit_length2.get_quantity())
                  << std::endl;
        std::cerr << "Units are not compatible!" << std::endl;
        exit(1);
    }

    Unit unit_length3 = UnitDefinitions::get_unit("m");
    if(unit_length3.get_quantity() != "length") {
        std::cerr << "Wrong quantity!" << std::endl;
        exit(1);
    }

    return 0;
}
