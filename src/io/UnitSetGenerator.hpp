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
 * @file UnitSetGenerator.hpp
 *
 * @brief Factory for UnitSet instances
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef UNITSETGENERATOR_HPP
#define UNITSETGENERATOR_HPP

#include <iostream>
#include <string>

#include "Error.hpp"
#include "UnitSet.hpp"

/**
  * \brief Class that generates a UnitSet based on an input string
  */
class UnitSetGenerator {
  public:
    /**
     * @brief Generate a UnitSet with the given name
     *
     * Currently supported:
     *  - SI units ("SI")
     *  - CGS units ("CGS")
     *  - galactic units ("galactic"): length in kpc, mass in solar mass, time
     *    in Gyr
     *
     * If another name is given, the program will abort.
     *
     * @param name Name of the UnitSet to generate
     * @return Pointer to a UnitSet, should be deleted by the caller
     */
    static UnitSet* generate(std::string name) {
        if(name == "SI") {
            Unit unit_length("length", "m", 1.);
            Unit unit_mass("mass", "kg", 1.);
            Unit unit_time("time", "s", 1.);
            return new UnitSet(unit_length, unit_mass, unit_time);
        }
        if(name == "CGS") {
            Unit unit_length("length", "cm", 0.01);
            Unit unit_mass("mass", "g", 0.001);
            Unit unit_time("time", "s", 1.);
            return new UnitSet(unit_length, unit_mass, unit_time);
        }
        if(name == "galactic") {
            Unit unit_length("length", "kpc", 3.08567758e19);
            Unit unit_mass("mass", "Msol", 1.9891e30);
            Unit unit_time("time", "Gyr", 3.154e16);
            return new UnitSet(unit_length, unit_mass, unit_time);
        }
        std::cerr << "Error! Unknown UnitSet: " << name << "!" << std::endl;
        my_exit();
        return NULL;
    }
};

#endif  // UNITSETGENERATOR_HPP
