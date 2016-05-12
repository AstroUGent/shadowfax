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
 * @file Physics.cpp
 *
 * @brief Physical constants used in the simulation: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Physics.hpp"
#include "io/PhysicalConstant.hpp"
#include "io/UnitSet.hpp"
using namespace std;

/**
 * @brief Constructor
 *
 * Initialize the physical constants in SI-units and store their values in
 * simulation units
 *
 * @param units UnitSet used internally during the simulation
 * @param real_units Flag indicating whether physical values or idealized values
 * should be used for the physical constants
 */
Physics::Physics(UnitSet& units, bool real_units) {
    Unit Gunit("length*length*length/mass/time/time", "Gunit", 1.);
    double Gval;

    if(real_units) {
        // hard coded values in SI-units
        // gravitational constant
        Gval = 6.674e-11;
    } else {
        // hard coded idealized values (1 in most cases)
        // gravitational constant
        Gval = 1.;
    }
    PhysicalConstant G(Gval, Gunit);
    _G = G.get_value(units);
}

/**
 * @brief Get the value of the gravitational constant in simulation units
 *
 * @return The value of the gravitational constant in simulation units
 */
double Physics::get_gravitational_constant() {
    return _G;
}
