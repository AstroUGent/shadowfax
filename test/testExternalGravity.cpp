/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testExternalGravity.cpp
 *
 * @brief Unit test for the ExternalGravity class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "ExternalGravity.hpp"

#include <iostream>

/**
 * @brief Unit test for the ExternalGravity class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char** argv) {
    ExternalGravity external_gravity;

    for(unsigned int i = 0; i < 1000; ++i) {
        Vec pos((i + 0.5) * 0.001 * 400., 0., 0.);
        Vec a = external_gravity.get_acceleration(pos, 6.);
        std::cout << pos.x() << "\t" << a.x() << std::endl;
    }

    return 0;
}
