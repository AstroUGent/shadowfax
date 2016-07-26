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
 * @file NDIrregTable.cpp
 *
 * @brief 3D non-rectangular Table: implementation
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#include "NDIrregTable.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <fstream>
#include <math.h>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

/**
 * @brief Get the index of the given value in the given array
 *
 * Notice the obsolete usage of "inline" in a source file. Only works in header
 * files!
 *
 * @param values Array with values
 * @param value Value to search
 * @return Index of the value in the array
 */
inline int find_in_set(vector<double>& values, double value) {
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

/**
 * @brief No idea what this function does
 *
 * @param values Array of values
 * @param collapsed Unknown flag
 * @param index Unknown index
 * @return No idea
 */
inline double width(vector<double>& values, int collapsed, int index) {
    switch(collapsed) {
        case 0:
            return values.at(index + 1) - values.at(index);
        case 1:
            return 1.;
        default:
            // should probably crash?
            return -1.;
    }
}

/**
 * @brief No idea what this function does
 *
 * @param values Array of values
 * @param value Value
 * @param offset Unknown offset
 * @param collapsed Unknown flag
 * @param index Unknown index
 * @return No idea
 */
inline double axis_range(vector<double>& values, double value, int offset,
                         int collapsed, int index) {
    // table collapsed on 1 axis?
    switch(collapsed) {
        case 0:
            switch(offset) {
                case 0:
                    return values.at(index + 1) - value;
                case 1:
                    return value - values.at(index);
                default:
                    // should probably crash
                    return -1.;
            }
        case 1:
            return 1.;
        default:
            // should probably crash
            return -1.;
    }
}

/**
 * @brief No idea what this function does
 *
 * @param values Array of values
 * @param value Value
 * @param offset Unknown offset
 * @param collapsed Unknown flag
 * @param index Unknown index
 * @return No idea
 */
inline double axis_range_sq(vector<double>& values, double value, int offset,
                            int collapsed, int index) {
    // table collapsed on 1 axis?
    switch(collapsed) {
        case 0: {
            double dx = value - values.at(index);
            switch(offset) {
                case 0: {
                    double dd = values.at(index + 1) - value;
                    return dd * dd - dx * dx;
                }
                case 1:
                    return dx * dx;
                default:
                    // should probably crash
                    return -1.;
            }
        }
        case 1:
            return 1.;

        default:
            // should probably crash
            return -1.;
    }
}

/**
 * @brief Function to find the index for the density when Fe and Mg indices are
 * known
 * Since the table is non-rectangular along this axis the indexes of Fe and Mg
 * need to be determined before to then find in which cell the n were are
 * interpolating resides
 *
 * @param value The density to be found
 * @param a Index of the Fe value
 * @param b Index of the Mg value
 *
 * @return Index of the value of n directly below the entered one
 */
int ThreeDIrregTable::find_n_in_table(double value, int a, int b) {
    double x, y;
    switch(_collapsed[0]) {
        case 0:
            x = _axisranges[3];
            break;
        case 1:
            x = 0.;
            break;
        default:
            // should probably crash
            x = -1.;
            break;
    }
    switch(_collapsed[1]) {
        case 0:
            y = _axisranges[4];
            break;
        case 1:
            y = 0.;
            break;
        default:
            // should probably crash
            y = -1.;
            break;
    }
    // binary tree
    int max = _table.at(a)[b].size() - 1;
    int sizem = max;
    int min = 0;
    int med;
    while(min < max - 1) {
        // the usage of floor seems unnecessary...
        // (int-int)/2 will still be an int, and you get floor for free
        med = floor(min + ((max - min) / 2));
        if(value <
           (1. - x) * (1. - y) * _table.at(a)[b][med].first +
                   x * (1. - y) *
                           _table.at(a + 1 - _collapsed[0])[b][med].first +
                   y * (1. - x) *
                           _table.at(a)[b + 1 - _collapsed[1]][med].first +
                   x * y *
                           _table.at(a + 1 -
                                     _collapsed[0])[b + 1 - _collapsed[1]]
                                                   [med].first) {
            max = med;
        } else {
            min = med;
        }
    }
    if(min < sizem) {
        return min;
    } else {
        return sizem - 1;
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
ThreeDIrregTable::ThreeDIrregTable(vector<string> filenames,
                                   UnitSet* simulation_units) {
    _Fe_values = vector<double>();
    _Mg_values = vector<double>();

    // units
    Unit unit_mass("mass", "g", 0.001);
    Unit unit_length("length", "cm", 0.01);
    Unit unit_density = unit_mass / unit_length / unit_length / unit_length;
    UnitConverter dens_conv(unit_density, simulation_units->get_density_unit());

    vector<string> strs;
    for(unsigned int i = 0; i < filenames.size(); i++) {
        boost::split(strs, filenames[i], boost::is_any_of("_"));
        _Fe_values.push_back(stod(strs[strs.size() - 4]));
        _Mg_values.push_back(stod(strs[strs.size() - 3]));
    }
    sort(_Fe_values.begin(), _Fe_values.end());
    _Fe_values.erase(unique(_Fe_values.begin(), _Fe_values.end()),
                     _Fe_values.end());

    sort(_Mg_values.begin(), _Mg_values.end());
    _Mg_values.erase(unique(_Mg_values.begin(), _Mg_values.end()),
                     _Mg_values.end());
    double N;
    // determines if table is collapsed on certain axes
    _collapsed = vector<int>({0, 0});
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

    // make a 3-dimensional vector of the correct size and fill it with 0's
    _table = vector<vector<vector<pair<double, double>>>>(_Fe_values.size());
    for(unsigned int i = 0; i < _Fe_values.size(); i++) {
        vector<vector<pair<double, double>>> row;
        _table.push_back(row);
        for(unsigned int j = 0; j < _Mg_values.size(); j++) {
            vector<pair<double, double>> row;
            _table.at(i).push_back(row);
        }
    }

    // fill in the vector with values from the files
    vector<string>::iterator it2;
    int i, j;
    string line;
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

        ifstream infile(*it2);

        std::getline(infile, line);
        std::istringstream iss(line);
        iss >> N;
        _table.at(i)[j].push_back(make_pair(log10(dens_conv.convert(N)),
                                            stod(strs[strs.size() - 1])));
    }
    // sort density values
    for(unsigned int i = 0; i < _Fe_values.size(); i++) {
        for(unsigned int j = 0; j < _Mg_values.size(); j++) {
            sort(_table.at(i)[j].begin(), _table.at(i)[j].end(),
                 boost::bind(&pair<double, double>::first, _1) <
                         boost::bind(&pair<double, double>::first, _2));
            // remove duplicates (from redshift tables)
            _table.at(i)[j].erase(
                    unique(_table.at(i)[j].begin(), _table.at(i)[j].end()),
                    _table.at(i)[j].end());
        }
    }

    _axisranges = array<double, 12>();
}

/**
 * @brief Return the value for the requested coordinates
 *
 * Return the interpolated value of nH for the given Fe, Mg and density
 *
 * @param value a vector of the requested values of Fe, Mg and density
 *
 * @return Interpolated value of the hydrogen density
 */
double ThreeDIrregTable::get_value(vector<double> value) {
    // can easily be switched to other method when necesary without changing
    // other files
    return this->linear_interp(value);
}

/**
 * @brief Linear interpolator
 *
 * Calculates the requested value of nH using multilinear, volume based
 *interpolation.
 * The volume are calculated via coordinate remapping since the table is
 * non-rectangular along the n-axis.
 *
 * @param value a vector of the requested values of Fe, Mg and density
 *
 * @return Linearly interpolated value of the hydrogen density
 */
double ThreeDIrregTable::linear_interp(vector<double> value) {
    /*GENERAL IDEA:
     * point lies withing truncated rectangular prism with 8 corners
     * planes determined by the point's coords split prism into 8
     * subvolumes
     * each value of the table of a corner has a contribution of
     * value*(volume of subcuboid it's part of)/(total volume of cuboid)
     * to calculate this:
     * 1. We calculate the several constants needed for calculating suvolume
     * 2. We iterate over all corners, calculating value*subvolume
     * 3. We divide by the total volume
     * The subvolumes can be found by using a coordinate map from our prism
     * to the unit cube and then calculating the Jacobian and integrating
     * in unit-cube-space
     */

    // find points; this index and this index+1
    // we need the ones of Fe and Mg first
    int i = find_in_set(_Fe_values, value[0]);
    int j = find_in_set(_Mg_values, value[1]);

    // save axis ranges of Fe and Mg in an array
    double temp = width(_Fe_values, _collapsed[0], i);
    for(int a = 0; a <= 1 - _collapsed[0]; ++a) {
        _axisranges[3 * a] =
                axis_range(_Fe_values, value[0], a, _collapsed[0], i) / temp;
        _axisranges[3 * (a + 2)] =
                axis_range_sq(_Fe_values, value[0], a, _collapsed[0], i) /
                (temp * temp);
    }
    temp = width(_Mg_values, _collapsed[1], j);
    for(int b = 0; b <= 1 - _collapsed[1]; ++b) {
        _axisranges[1 + 3 * b] =
                axis_range(_Mg_values, value[1], b, _collapsed[1], j) / temp;
        _axisranges[1 + 3 * (b + 2)] =
                axis_range_sq(_Mg_values, value[1], b, _collapsed[1], j) /
                (temp * temp);
    }

    // then we can find the index for n
    int k = find_n_in_table(log10(value[2]), i, j);

    // calculate constants needed in the volume calculation
    int di = 1 - _collapsed[0];
    int dj = 1 - _collapsed[1];
    double a0 = _table.at(i)[j][k].first;
    double a1 = _table.at(i + di)[j][k].first - _table.at(i)[j][k].first;
    double a2 = _table.at(i)[j + dj][k].first - _table.at(i)[j][k].first;
    double a3 = _table.at(i)[j][k + 1].first - _table.at(i)[j][k].first;
    double a4 = _table.at(i)[j][k].first + _table.at(i + di)[j + dj][k].first -
                _table.at(i + di)[j][k].first - _table.at(i)[j + dj][k].first;
    double a5 = _table.at(i)[j][k].first + _table.at(i)[j + dj][k + 1].first -
                _table.at(i)[j][k + 1].first - _table.at(i)[j + dj][k].first;
    double a6 = _table.at(i)[j][k].first + _table.at(i + di)[j][k + 1].first -
                _table.at(i + di)[j][k].first - _table.at(i)[j][k + 1].first;
    double a7 = _table.at(i)[j][k + 1].first +
                _table.at(i + di)[j + dj][k + 1].first +
                _table.at(i + di)[j][k].first + _table.at(i)[j + dj][k].first -
                _table.at(i + di)[j][k + 1].first -
                _table.at(i)[j + dj][k + 1].first -
                _table.at(i + di)[j + dj][k].first - _table.at(i)[j][k].first;

    // calculate the axis ranges for n
    temp = (log10(value[2]) - a0 - a1 * value[0] - a2 * value[1] -
            a4 * value[0] * value[1]) /
           (a3 + a5 * value[1] + a6 * value[0] + a7 * value[0] * value[1]);
    _axisranges[2] = 1. - temp;
    _axisranges[5] = temp;

    // iterate over corners
    double res = 0.;
    for(int a = 0; a <= 1 - _collapsed[0]; ++a) {
        for(int b = 0; b <= 1 - _collapsed[1]; ++b) {
            for(int c = 0; c <= 1; ++c) {
                res += (_axisranges[2 + 3 * c] *
                        (a3 * _axisranges[3 * a] * _axisranges[1 + 3 * b] +
                         a5 * _axisranges[3 * a] *
                                 _axisranges[1 + 3 * (b + 2)] / 2 +
                         a6 * _axisranges[3 * (a + 2)] *
                                 _axisranges[1 + 3 * b] / 2 +
                         a7 * _axisranges[3 * (a + 2)] *
                                 _axisranges[1 + 3 * (b + 2)] / 4) *
                        _table.at(i + a)[j + b][k + c].second);
            }
        }
    }
    return res / abs(a3 + 0.5 * a5 + 0.5 * a6 + 0.25 * a7);
}

/**
 * @brief Dump the object to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void ThreeDIrregTable::dump(RestartFile& rfile) {
    rfile.write(_Fe_values);
    rfile.write(_Mg_values);
    rfile.write(_collapsed);
    rfile.write(_table);
    rfile.write(_axisranges);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
ThreeDIrregTable::ThreeDIrregTable(RestartFile& rfile) {
    rfile.read(_Fe_values);
    rfile.read(_Mg_values);
    rfile.read(_collapsed);
    rfile.read(_table);
    rfile.read(_axisranges);
}
