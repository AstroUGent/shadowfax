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
 * @file UnitSet.hpp
 *
 * @brief Set of units: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef UNITSET_HPP
#define UNITSET_HPP

#include "Unit.hpp"  // for Unit
#include <string>    // for string

class RestartFile;

/**
  * @brief Collection of units which are mutually consistent
  *
  * In other words: if you divide e.g. the length unit associated with the set
  * by the time unit associated with the set, you arrive at the velocity unit
  * associated with the set.
  */
class UnitSet {
  private:
    /*! @brief Length Unit associated with this UnitSet */
    Unit _unit_position;

    /*! @brief Density Unit associated with this UnitSet */
    Unit _unit_density;

    /*! @brief Velocity Unit associated with this UnitSet */
    Unit _unit_velocity;

    /*! @brief Pressure Unit associated with this UnitSet */
    Unit _unit_pressure;

    /*! @brief Time Unit associated with this UnitSet */
    Unit _unit_time;

    /*! @brief Mass Unit associated with this UnitSet */
    Unit _unit_mass;

    /*! @brief Acceleration Unit associated with this UnitSet */
    Unit _unit_acceleration;

    /*! @brief Energy Unit associated with this UnitSet */
    Unit _unit_energy;

  public:
    UnitSet();
    UnitSet(Unit unit_length, Unit unit_mass, Unit unit_time);

    ~UnitSet() {}

    Unit get_length_unit();
    Unit get_density_unit();
    Unit get_velocity_unit();
    Unit get_pressure_unit();
    Unit get_time_unit();
    Unit get_mass_unit();
    Unit get_acceleration_unit();
    Unit get_energy_unit();

    Unit get_unit(std::string quantity);

    void dump(RestartFile& rfile);
    UnitSet(RestartFile& rfile);
};

#endif  // UNITSET_HPP
