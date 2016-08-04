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
 * @file NDIrregTable.hpp
 *
 * @brief 3D non-rectangular Table: header
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#ifndef NDIRREGTABLE_HPP
#define NDIRREGTABLE_HPP

#include <array>
#include <string>   // for string
#include <utility>  // for pair
#include <vector>   // for vector

class RestartFile;
class UnitSet;

/**
  * @brief 3D table that is non-rectangular along its third axis
  *
  * Used for interpolating from Fe, Mg and n to nH. These tables are non-
  * rectangular along the n-axis and thus need special treatment.
  * The values are saved in a 3D vector of pair of doubles.
  * Each pair contains the value of n and nH respectively.
  */
class ThreeDIrregTable {
  private:
    /*! @brief Values of Fe for which the table has values */
    std::vector<double> _Fe_values;
    /*! @brief Values of Mg for which the table has values */
    std::vector<double> _Mg_values;
    /*! @brief Ints that determine if the table is flattened along an axis */
    std::vector<int> _collapsed;
    /*! @brief 3D vector of pair that contain the n & nH values */
    std::vector<std::vector<std::vector<std::pair<double, double>>>> _table;
    /*! @brief array used to store values in interpolation */
    std::array<double, 12> _axisranges;

    double linear_interp(std::vector<double> value);
    int find_n_in_table(double value, int a, int b);

  public:
    ThreeDIrregTable(std::vector<std::string> filenames,
                     UnitSet* simulation_units);

    double get_value(std::vector<double> value);

    void dump(RestartFile& rfile);
    ThreeDIrregTable(RestartFile& rfile);
};

#endif  // NDIRREGTABLE_HPP
