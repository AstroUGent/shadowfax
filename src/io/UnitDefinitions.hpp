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
        if(name == "amu") {
            // check value!
            return Unit("mass", "amu", 1.660539040e-27);
        } else if(name == "cm") {
            return Unit("length", "cm", 0.01);
        } else if(name == "erg") {
            return Unit("length*length*mass/time/time", "erg", 1.e7);
        } else if(name == "g") {
            return Unit("mass", "g", 0.001);
        } else if(name == "Gyr") {
            return Unit("time", "Gyr", 3.154e16);
        } else if(name == "h") {
            return Unit("time", "h", 3600.);
        } else if(name == "J") {
            return Unit("length*length*mass/time/time", "J", 1.);
        } else if(name == "K") {
            return Unit("temperature", "K", 1.);
        } else if(name == "kg") {
            return Unit("mass", "kg", 1.);
        } else if(name == "kg^-1") {
            return Unit("1/mass", "kg^-1", 1.);
        } else if(name == "kg m^-3") {
            return Unit("mass/length/length/length", "kg m^-3", 1.);
        } else if(name == "kpc") {
            return Unit("length", "kpc", 3.08567758e19);
        } else if(name == "m") {
            return Unit("length", "m", 1.);
        } else if(name == "Msol") {
            return Unit("mass", "Msol", 1.9891e30);
        } else if(name == "s") {
            return Unit("time", "s", 1.);
        } else if(name == "yr") {
            return Unit("time", "yr", 3.15576e7);
        } else if(name == "") {
            // this one is always last
            return Unit("dimensionless", "", 1.);
        } else {
            std::cerr << "Unknown unit: " << name << std::endl;
            abort();
            return Unit();
        }
    }

    /**
     * @brief Get the quantity represented by the given quantity, converting
     * quantity aliases to their proper basic quantity representations.
     *
     * @param quantity Quantity, can be either a basic quantity expression or an
     * aliased quantity (expressions containing aliased quantities are not
     * supported).
     * @return Basic quantity representation of the given quantity. This
     * representation only contains the five basic SI quantities: length, mass,
     * time, temperature and current.
     */
    inline static std::string get_quantity(std::string quantity) {
        if(quantity == "energy") {
            return "length*length*mass/time/time";
        } else if(quantity == "mass_density") {
            return "mass/length/length/length";
        } else {
            // no alias found, return the input string
            return quantity;
        }
    }
};

#endif  // UNITDEFINITIONS_HPP
