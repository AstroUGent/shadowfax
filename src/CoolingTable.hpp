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
 * @file CoolingTable.hpp
 *
 * @brief Cooling Table: header
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#ifndef COOLINGTABLE_HPP
#define COOLINGTABLE_HPP

#include <string>  // for string
#include <vector>  // for vector

class FiveDTable;
class RestartFile;
class ThreeDIrregTable;
class UnitSet;

/**
  * @brief Class that calculates cooling rate for a given Fe, Mg, z, T and n
  *
  * implementation can depend on available cooling tables
  */
class CoolingTable {
  private:
    /*! @brief Value of cooling for Fe, Mg, z, T and nH */
    FiveDTable* _coolingtable;
    /*! @brief Value of nH for Fe, Mg and n */
    ThreeDIrregTable* _ntable;
    /*! @brief Previous density value */
    double _n_prev;
    /*! @brief Previous Fe value */
    double _Fe_prev;
    /*! @brief Previous Mg value */
    double _Mg_prev;
    /*! @brief Previous hydrogen density value */
    double _nH_prev;

    /*! @brief Maximum temperature value that is tabulated */
    double _T_max;
    /*! @brief Maximum redshift value that is tabulated */
    double _z_max;

    std::vector<std::string> get_file_list(std::string directory);

  public:
    CoolingTable(std::string directory, UnitSet* simulation_units);
    ~CoolingTable();

    double get_value(double Fe, double Mg, double rs, double n, double T);

    void dump(RestartFile& rfile);
    CoolingTable(RestartFile& rfile);
};

#endif  // COOLINGTABLE_HPP
