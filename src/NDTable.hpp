/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 *                    Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file NDTable.hpp
 *
 * @brief ND Rectangular table: header
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#ifndef NDTABLE_HPP
#define NDTABLE_HPP
#include "io/Unit.hpp"
#include "io/UnitConverter.hpp"
#include "io/UnitSet.hpp"
#include <array>
#include <dirent.h>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

class RestartFile;

/**
  * @brief 5D table that is rectangular
  *
  * Used for interpolating from Fe, Mg, z, T and nH to cooling rate.
  * The values are saved in a 5D vector of doubles.
  */
class FiveDTable {
  private:
    /*! @brief Values of Fe for which the table has values */
    vector<double> _Fe_values;
    /*! @brief Values of Mg for which the table has values */
    vector<double> _Mg_values;
    /*! @brief Values of redshift z for which the table has values */
    vector<double> _redshift_values;
    /*! @brief Values of nH for which the table has values */
    vector<double> _nH_values;
    /*! @brief Values of T for which the table has values */
    vector<double> _T_values;
    /*! @brief Ints that determine if the table is flattened along an axis */
    vector<int> _collapsed;
    /*! @brief 5D vector of doules that contain the cooling rate values */
    vector<vector<vector<vector<vector<double>>>>> _table;
    /*! @brief array used to store values in interpolation */
    array<double, 10> _axisranges;

    double linear_interp(vector<double> value);

  public:
    FiveDTable(vector<string> filenames, UnitSet* simulation_units);
    double get_value(vector<double> value);

    double get_T_max();
    double get_z_max();

    void dump(RestartFile& rfile);
    FiveDTable(RestartFile& rfile);
};

#endif  // NDTABLE_HPP
