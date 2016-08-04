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
 * @brief 5D rectangular Table: implementation
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#include "NDTable.hpp"
#include "RestartFile.hpp"
#include "io/Unit.hpp"           // for operator/, Unit, operator*
#include "io/UnitConverter.hpp"  // for UnitConverter
#include "io/UnitSet.hpp"        // for UnitSet
#include <algorithm>             // for sort, unique, upper_bound
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <float.h>   // for DBL_MAX
#include <iostream>  // for ifstream
#include <iterator>  // for distance
#include <math.h>    // for pow
#include <sstream>
#include <string>  // for basic_string, string
#include <vector>  // for vector, vector<>::iterator
using namespace std;

inline int find_in_set1(vector<double>& values, double value) {
    // check if set is not 1 value
    int s = values.size();
    if(s > 1) {
        int res = distance(values.begin(),
                           upper_bound(values.begin(), values.end(), value)) -
                  1;
        if(res < 0) {
            return 0;
        } else if(res < s - 1) {
            return res;
        } else {
            return s - 2;
        }
    } else {
        return 0;
    }
}

inline double axis_range1(vector<double>& values, double value, int offset,
                          int collapse, int index) {
    // table collapsed on 1 axis?
    switch(collapse) {
        case 0:
            switch(offset) {
                case 0:
                    return values.at(index + 1) - value;
                case 1:
                    return value - values.at(index);
                default:
                    return -1.;
            }
        case 1:
            return 1.;

        default:
            return -1.;
    }
}

inline double width(vector<double>& values, int collapsed, int index) {
    switch(collapsed) {
        case 0:
            return values.at(index + 1) - values.at(index);
        case 1:
            return 1.;
        default:
            return -1.;
    }
}

/**
 * @brief Constructor
 *
 * Loads in the tables from the files and puts them in the vector
 *
 * @param filenames The list of filenames of the cooling tables
 * @param simulation_units The UnitSet used in the simulation
 */
FiveDTable::FiveDTable(vector<string> filenames, UnitSet* simulation_units) {
    _Fe_values = vector<double>();
    _Mg_values = vector<double>();
    _redshift_values = vector<double>();
    _nH_values = vector<double>();
    _T_values = vector<double>();

    // units
    Unit unit_mass("mass", "g", 0.001);
    Unit unit_length("length", "cm", 0.01);
    Unit unit_time("time", "s", 1.);

    // erg
    Unit unit_energy =
            unit_mass * unit_length * unit_length / unit_time / unit_time;
    Unit unit_cooling =
            unit_energy / unit_time / unit_length / unit_length / unit_length;

    Unit sim_energy = simulation_units->get_mass_unit() *
                      simulation_units->get_length_unit() *
                      simulation_units->get_length_unit() /
                      simulation_units->get_time_unit() /
                      simulation_units->get_time_unit();
    Unit sim_cooling = sim_energy / simulation_units->get_time_unit() /
                       simulation_units->get_length_unit() /
                       simulation_units->get_length_unit() /
                       simulation_units->get_length_unit();

    // UnitConverter
    UnitConverter cool_conv(unit_cooling, sim_cooling);

    vector<string> strs;
    for(unsigned int i = 0; i < filenames.size(); i++) {
        boost::split(strs, filenames[i], boost::is_any_of("_"));
        _Fe_values.push_back(stod(strs[strs.size() - 4]));
        _Mg_values.push_back(stod(strs[strs.size() - 3]));
        _redshift_values.push_back(stod(strs[strs.size() - 2]));
        _nH_values.push_back(stod(strs[strs.size() - 1]));
    }
    // read the Temperature values from the first file in the list
    ifstream infile(filenames[0]);
    string line;
    double T;
    // get rid of the first 2 lines, they don't contain the correct info
    std::getline(infile, line);
    std::getline(infile, line);
    while(std::getline(infile, line)) {
        std::istringstream iss(line);
        if(!(iss >> T)) {
            break;  // error
        } else {
            _T_values.push_back(T);
        }
    }
    int zero_T;
    if(!_T_values.at(0)) {
        zero_T = 0;
    } else {
        _T_values.push_back(0.);
        zero_T = 1;
    }
    // sort and erase duplicates
    sort(_Fe_values.begin(), _Fe_values.end());
    _Fe_values.erase(unique(_Fe_values.begin(), _Fe_values.end()),
                     _Fe_values.end());

    sort(_Mg_values.begin(), _Mg_values.end());
    _Mg_values.erase(unique(_Mg_values.begin(), _Mg_values.end()),
                     _Mg_values.end());

    sort(_redshift_values.begin(), _redshift_values.end());
    _redshift_values.erase(
            unique(_redshift_values.begin(), _redshift_values.end()),
            _redshift_values.end());

    sort(_nH_values.begin(), _nH_values.end());
    _nH_values.erase(unique(_nH_values.begin(), _nH_values.end()),
                     _nH_values.end());

    sort(_T_values.begin(), _T_values.end());
    _T_values.erase(unique(_T_values.begin(), _T_values.end()),
                    _T_values.end());

    _collapsed = vector<int>({0, 0, 0, 0, 0});
    // determines if table is collapsed on certain axes
    if(_Fe_values.size() == 1) {
        _collapsed[0] = 1;
    } else {
        _collapsed[0] = 0;
    }
    if(_Mg_values.size() == 1) {
        _collapsed[1] = 1;
    } else {
        _collapsed[1] = 0;
    }
    if(_redshift_values.size() == 1) {
        _collapsed[2] = 1;
    } else {
        _collapsed[2] = 0;
    }
    if(_nH_values.size() == 1) {
        _collapsed[3] = 1;
    } else {
        _collapsed[3] = 0;
    }
    if(_T_values.size() == 1) {
        _collapsed[4] = 1;
    } else {
        _collapsed[4] = 0;
    }

    // make a 5-dimensional vector of the correct size and fill it with 0's
    _table = vector<vector<vector<vector<vector<double>>>>>(_Fe_values.size());
    for(unsigned int i = 0; i < _Fe_values.size(); i++) {
        vector<vector<vector<vector<double>>>> row;
        _table.push_back(row);
        for(unsigned int j = 0; j < _Mg_values.size(); j++) {
            vector<vector<vector<double>>> row;
            _table.at(i).push_back(row);
            for(unsigned int k = 0; k < _redshift_values.size(); k++) {
                vector<vector<double>> row;
                _table.at(i)[j].push_back(row);
                for(unsigned int l = 0; l < _nH_values.size(); l++) {
                    vector<double> row;
                    _table.at(i)[j][k].push_back(row);
                    for(unsigned int m = 0; m < _T_values.size(); m++) {
                        _table.at(i)[j][k][l].push_back(0.);
                    }
                }
            }
        }
    }
    // fill in the vector with values from the files
    vector<string>::iterator it2;
    int i, j, k, l;
    double res;
    int n = 0;

    for(it2 = filenames.begin(); it2 != filenames.end(); ++it2) {
        boost::split(strs, *it2, boost::is_any_of("_"));
        // find indices for table; set is sorted and unique-valued (which we
        // need) but doesn't have indices
        i = distance(_Fe_values.begin(),
                     find(_Fe_values.begin(), _Fe_values.end(),
                          stod(strs[strs.size() - 4])));
        j = distance(_Mg_values.begin(),
                     find(_Mg_values.begin(), _Mg_values.end(),
                          stod(strs[strs.size() - 3])));
        k = distance(_redshift_values.begin(),
                     find(_redshift_values.begin(), _redshift_values.end(),
                          stod(strs[strs.size() - 2])));
        l = distance(_nH_values.begin(),
                     find(_nH_values.begin(), _nH_values.end(),
                          stod(strs[strs.size() - 1])));
        ifstream infile(*it2);

        // get rid of the first 2 lines, they don't contain the correct info
        std::getline(infile, line);
        std::getline(infile, line);
        n = zero_T;
        while(std::getline(infile, line)) {
            std::istringstream iss(line);
            if(!(iss >> T >> res)) {
                break;  // error
            } else {
                _table.at(i)[j][k][l][n] = log10(cool_conv.convert(res));
                //_table.at(i)[j][k][l][n] = _cool_conv.convert(res);
                n++;
            }
        }
        if(zero_T == 1) {
            // fill in log(0) for 0 K, which is -infinity; approximate by
            // smallest negative double
            _table.at(i)[j][k][l][0] = -DBL_MAX;
        }
    }
    _axisranges = array<double, 10>();
}

/**
 * @brief Return the value for the requested coordinates
 *
 * Return the interpolated value of cooling rate for the given Fe, Mg, z, T and
 * nH
 *
 * @param value a vector of the requested values of Fe, Mg, z, T and nH
 *
 * @return Interpolated value of the cooling rate
 */
double FiveDTable::get_value(vector<double> value) {
    // can easily be switched to other method when necesary without changing
    // other files
    return this->linear_interp(value);
}

/**
 * @brief Get maximal temperature value that is tabulated
 *
 * @return Highest temperature value that is tabulated
 */
double FiveDTable::get_T_max() {
    return _T_values.back();
}

/**
 * @brief Get maximal redshift value that is tabulated
 *
 * @return Highest redshift value that is tabulated
 */
double FiveDTable::get_z_max() {
    return _redshift_values.back();
}

/**
 * @brief Linear interpolator
 *
 * Calculates the requested value of the cooling rate using multilinear,
 *  volume based interpolation.
 *
 * @param value a vector of the requested values of Fe, Mg, z, T and nH
 *
 * @return Lineary interpolated value of the cooling rate
 */
double FiveDTable::linear_interp(vector<double> value) {
    /*GENERAL IDEA:
     * point lies withing 5-cuboid with 32 corners
     * planes determined by the point's coords split hypercuboid into 32
     * subcuboids
     * each value of the table of a corner has a contribution of
     * value*(volume of subcuboid it's part of)/(total volume of cuboid)
     * to calculate this:
     * 1. We find the 2 "points" on each axis between which our points lie
     * 2. We iterate over all corners, calculating value*subvolume
     * 3. We divide by the total volume
     */

    // find points; this index and this index+1
    int i = find_in_set1(_Fe_values, value[0]);
    int j = find_in_set1(_Mg_values, value[1]);
    int k = find_in_set1(_redshift_values, value[2]);
    int l = find_in_set1(_nH_values, value[3]);
    int m = find_in_set1(_T_values, value[4]);

    // calculate and save axis ranges
    for(int a = 0; a <= 1 - _collapsed[0]; a++) {
        _axisranges[0 + 5 * a] =
                axis_range1(_Fe_values, value[0], a, _collapsed[0], i);
    }
    for(int b = 0; b <= 1 - _collapsed[1]; b++) {
        _axisranges[1 + 5 * b] =
                axis_range1(_Mg_values, value[1], b, _collapsed[1], j);
    }
    for(int c = 0; c <= 1 - _collapsed[2]; c++) {
        _axisranges[2 + 5 * c] =
                axis_range1(_redshift_values, value[2], c, _collapsed[2], k);
    }
    for(int d = 0; d <= 1 - _collapsed[3]; d++) {
        _axisranges[3 + 5 * d] =
                axis_range1(_nH_values, value[3], d, _collapsed[3], l);
    }
    for(int e = 0; e <= 1 - _collapsed[4]; e++) {
        _axisranges[4 + 5 * e] =
                axis_range1(_T_values, value[4], e, _collapsed[4], m);
    }

    // iterate over corners
    double res = 0.;
    for(int a = 0; a <= 1 - _collapsed[0]; ++a) {
        for(int b = 0; b <= 1 - _collapsed[1]; ++b) {
            for(int c = 0; c <= 1 - _collapsed[2]; ++c) {
                for(int d = 0; d <= 1 - _collapsed[3]; ++d) {
                    for(int e = 0; e <= 1 - _collapsed[4]; ++e) {
                        res += (_axisranges[5 * a] * _axisranges[1 + 5 * b] *
                                _axisranges[2 + 5 * c] *
                                _axisranges[3 + 5 * d] *
                                _axisranges[4 + 5 * e] *
                                _table.at(i + a)[j + b][k + c][l + d][m + e]);
                    }
                }
            }
        }
    }
    // calculate full volume of hypercuboid
    double vol_tot = width(_Fe_values, _collapsed[0], i) *
                     width(_Mg_values, _collapsed[1], j) *
                     width(_redshift_values, _collapsed[2], k) *
                     width(_nH_values, _collapsed[3], l) *
                     width(_T_values, _collapsed[4], m);

    return pow(10., res / vol_tot);
}

/**
 * @brief Dump the object to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void FiveDTable::dump(RestartFile& rfile) {
    rfile.write(_Fe_values);
    rfile.write(_Mg_values);
    rfile.write(_redshift_values);
    rfile.write(_nH_values);
    rfile.write(_T_values);
    rfile.write(_collapsed);
    rfile.write(_table);
    rfile.write(_axisranges);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
FiveDTable::FiveDTable(RestartFile& rfile) {
    rfile.read(_Fe_values);
    rfile.read(_Mg_values);
    rfile.read(_redshift_values);
    rfile.read(_nH_values);
    rfile.read(_T_values);
    rfile.read(_collapsed);
    rfile.read(_table);
    rfile.read(_axisranges);
}
