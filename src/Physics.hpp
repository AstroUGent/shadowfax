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

class ParameterFile;
class RestartFile;
class UnitSet;

#define PHYSICS_DEFAULT_MEANMOLWEIGHT "1.219512195 amu"

/**
 * @brief Physical constants used in the simulation
 *
 * Physical constants have a constant, hard-coded value in SI units. We store
 * them in simulation units, which might require unit conversion.
 */
class Physics {
  private:
    /*! @brief The gravitational constant G */
    double _G;
    /*! @brief The Hubble constant, assuming a value of 100 km/s/Mpc */
    double _H0;
    /*! @brief Mean molecular weight */
    double _mean_mol_weight;

  public:
    Physics(UnitSet& units, double mean_mol_weight, bool real_units = true);
    Physics(UnitSet& units, ParameterFile* parameters);

    double get_gravitational_constant();
    double get_Hubble_constant();
    double get_mean_mol_weight();

    void dump(RestartFile& rfile);
    Physics(RestartFile& rfile);
};

#endif  // PHYSICS_HPP
