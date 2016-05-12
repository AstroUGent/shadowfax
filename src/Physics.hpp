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
 * @file Physics.hpp
 *
 * @brief Physical constants used in the simulation: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PHYSICS_HPP
#define PHYSICS_HPP

class UnitSet;

/**
 * @brief Physical constants used in the simulation
 *
 * Physical constants have a constant, hard-coded value in SI units. We store
 * them in simulation units, which might require unit conversion.
 */
class Physics {
  private:
    /*! \brief The gravitational constant G */
    double _G;

  public:
    Physics(UnitSet& units, bool real_units = true);

    double get_gravitational_constant();
};

#endif  // PHYSICS_HPP
