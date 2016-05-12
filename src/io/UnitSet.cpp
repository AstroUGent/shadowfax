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
 * @file UnitSet.cpp
 *
 * @brief Set of units: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "UnitSet.hpp"
#include <iostream>
using namespace std;

/**
  * \brief Default constructor.
  *
  * The idea would be to make a few predefined sets to make life easier,
  * but currently the only predefined set are the default SI units
  *
  * @param type UnitType specifying the predefined set, currently only
  * UNITS_DEFAULT
  */
UnitSet::UnitSet(UnitType type) {
    _unit_position = Unit("length", "m", 1.);
    _unit_velocity = Unit("length/time", "m/s", 1.);
    _unit_density = Unit("mass/length/length/length", "kg/m^3", 1.);
    _unit_pressure = Unit("mass/length/time/time", "Pascal", 1.);
    _unit_time = Unit("time", "s", 1.);
    _unit_mass = Unit("mass", "kg", 1.);
}

/**
  * \brief Construct a UnitSet from a given length, mass and time Unit
  *
  * The other units are derived from these using multiplication and division.
  * As a result, all units are mutually consistent.
  *
  * @param unit_length Length Unit
  * @param unit_mass Mass Unit
  * @param unit_time Time Unit
  */
UnitSet::UnitSet(Unit unit_length, Unit unit_mass, Unit unit_time) {
    _unit_position = unit_length;
    _unit_velocity = unit_length / unit_time;
    _unit_density = unit_mass / unit_length / unit_length / unit_length;
    _unit_pressure = unit_mass / unit_length / unit_time / unit_time;
    _unit_time = unit_time;
    _unit_mass = unit_mass;
}

/**
  * @brief et the length unit
  *
  * @returns The length Unit
  */
Unit UnitSet::get_length_unit() {
    return _unit_position;
}

/**
  * @brief Get the density unit
  *
  * @returns The density Unit
  */
Unit UnitSet::get_density_unit() {
    return _unit_density;
}

/**
  * @brief Get the velocity unit
  *
  * @returns The velocity Unit
  */
Unit UnitSet::get_velocity_unit() {
    return _unit_velocity;
}

/**
  * @brief Get the pressure unit
  *
  * @returns The pressure Unit
  */
Unit UnitSet::get_pressure_unit() {
    return _unit_pressure;
}

/**
  * @brief Get the time unit
  *
  * @returns The time Unit
  */
Unit UnitSet::get_time_unit() {
    return _unit_time;
}

/**
  * @brief Get the mass unit
  *
  * @returns The mass Unit
  */
Unit UnitSet::get_mass_unit() {
    return _unit_mass;
}

/**
  * @brief Get the unit for the given quantity
  *
  * This function only works if the quantity is expressed with the basic
  * SI-quantities length, mass and time, using multiplication and division. The
  * expression is broken down and the unit is then constructed using
  * multiplications and divisions involving the basic SI-quantities.
  *
  * @param quantity Quantity for which to obtain the Unit
  * @returns The Unit of the given quantity in the system of units specified by
  * this UnitSet
  */
Unit UnitSet::get_unit(string quantity) {
    unsigned int pos = 0;
    string names[3] = {"length", "mass", "time"};
    Unit units[3] = {_unit_position, _unit_mass, _unit_time};
    Unit unit;
    bool multiply = true;
    while(pos < quantity.length()) {
        unsigned int i = 0;
        while(i < 3 && quantity.find(names[i], pos) > pos) {
            i++;
        }
        if(i == 3) {
            if(quantity.find("*", pos) == pos) {
                multiply = true;
            } else {
                multiply = false;
            }
            pos += 1;
        } else {
            if(unit.get_name() == "dimensionless") {
                unit = units[i];
            } else {
                if(multiply) {
                    unit *= units[i];
                } else {
                    unit /= units[i];
                }
            }
            pos += names[i].length();
        }
    }
    return unit;
}

/**
 * @brief Dump the unit set to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void UnitSet::dump(RestartFile& rfile) {
    _unit_position.dump(rfile);
    _unit_density.dump(rfile);
    _unit_velocity.dump(rfile);
    _unit_pressure.dump(rfile);
    _unit_time.dump(rfile);
    _unit_mass.dump(rfile);
}

/**
 * @brief Restart constructor. Initialize the unit set from the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 */
UnitSet::UnitSet(RestartFile& rfile)
        : _unit_position(rfile), _unit_density(rfile), _unit_velocity(rfile),
          _unit_pressure(rfile), _unit_time(rfile), _unit_mass(rfile) {}
