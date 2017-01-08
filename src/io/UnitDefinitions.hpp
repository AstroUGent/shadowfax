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
 * @file UnitDefinitions.hpp
 *
 * @brief Definitions of units that can be used in the parameter and output
 * files.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef UNITDEFINITIONS_HPP
#define UNITDEFINITIONS_HPP

#include "Unit.hpp"

/**
 * @brief Definitions of units that can be used in the parameter and output
 * files.
 */
class UnitDefinitions {
  public:
    /**
     * @brief Get the unit corresponding to the given name.
     *
     * @param name Name of a unit.
     * @return Unit instance with the correct quantity and SI value for that
     * unit.
     */
    inline static Unit get_unit(std::string name) {
        // units are ordered alphabetically, irrespective of their quantity
        if(name == "cm") {
            return Unit("length", "cm", 0.01);
        } else if(name == "g") {
            return Unit("mass", "g", 0.001);
        } else if(name == "Gyr") {
            return Unit("time", "Gyr", 3.154e16);
        } else if(name == "kg") {
            return Unit("mass", "kg", 1.);
        } else if(name == "kpc") {
            return Unit("length", "kpc", 3.08567758e19);
        } else if(name == "m") {
            return Unit("length", "m", 1.);
        } else if(name == "Msol") {
            return Unit("mass", "Msol", 1.9891e30);
        } else if(name == "s") {
            return Unit("time", "s", 1.);
        } else if(name == "") {
            // this one is always last
            return Unit("dimensionless", "", 1.);
        } else {
            std::cerr << "Unknown unit: " << name << std::endl;
            abort();
            return Unit();
        }
    }
};

#endif  // UNITDEFINITIONS_HPP
