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
 * @file testCosmology.cpp
 *
 * @brief Unit test for the cosmological integration
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Cosmology.hpp"
#include <iostream>

/**
 * @brief Unit test for the Cosmology class
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {

    double H0 = 2.e-18;  // per second
    Cosmology cosmo(H0, 0.28, 0.72);

    for(double x = 0.01; x < 1.; x += 0.01) {
        std::cout << x << "\t" << cosmo.get_velocity_factor(0.01, x) << "\t"
                  << cosmo.get_acceleration_factor(0.01, x) << std::endl;
    }

    return 0;
}
