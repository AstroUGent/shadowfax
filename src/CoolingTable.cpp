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
 * @file CoolingTable.cpp
 *
 * @brief Cooling Table: implementation
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#include "CoolingTable.hpp"
#include "CoolingLocation.hpp"
#include "Error.hpp"         // for my_exit
#include "NDIrregTable.hpp"  // for ThreeDIrregTable
#include "NDTable.hpp"       // for FiveDTable
#include "RestartFile.hpp"
#include <algorithm>  // for min
#include <cstddef>    // for NULL
#include <dirent.h>   // for closedir, opendir, readdir, etc
#include <iostream>   // for operator<<, basic_ostream, etc
#include <string>     // for string, allocator, etc
#include <vector>     // for vector
using namespace std;

/**
 * @brief Cooling Table file scanner
 *
 * Makes a list of all the files in the folder /coolingtables and returns it
 *
 * @param directory Folder containing the cooling tables on the local system
 * @return List of cooling table files
 */
vector<string> CoolingTable::get_file_list(string directory) {
    if(directory == "") {
        directory = COOLING_LOCATION;
    }
    DIR* dir;
    struct dirent* ent;
    vector<string> result;
    cout << endl
         << "Looking for cooling tables in " << directory << "..." << endl;
    int n = 0;
    if((dir = opendir(directory.c_str())) != NULL) {
        while((ent = readdir(dir)) != NULL) {
            string str(ent->d_name);
            if(str.find(".rates") != string::npos) {
                result.push_back(directory + str);
                n++;
            }
        }
        closedir(dir);
        cout << n << " tables found." << endl;
    } else {
        cout << "None found." << endl;
    }
    cout << endl;
    return result;
}
/**
 * @brief Constructor
 *
 * Creates the 3D and 5D table from the given list of files
 *
 * @param directory Folder containing the cooling tables on the system
 * @param simulation_units The UnitSet used in the simulation
 */
CoolingTable::CoolingTable(string directory, UnitSet* simulation_units) {
    vector<string> filenames = get_file_list(directory);
    _coolingtable = new FiveDTable(filenames, simulation_units);
    _ntable = new ThreeDIrregTable(filenames, simulation_units);
    // must be out of real range
    _n_prev = -1.;
    _Fe_prev = -200.;
    _Mg_prev = -200.;
    _nH_prev = -1.;

    _T_max = _coolingtable->get_T_max();
    _z_max = _coolingtable->get_z_max();
}

/**
 * @brief Destructor.
 */
CoolingTable::~CoolingTable() {
    delete _coolingtable;
    delete _ntable;
}

/**
 * @brief Return the value for the requested coordinates
 *
 * Return the interpolated value of cooling for the given Fe, Mg, z, T and
 *density
 *
 * @param Fe Iron abundance
 * @param Mg Magnesium abundance
 * @param z Redshift
 * @param n the density
 * @param T the temperature
 *
 * @return Interpolated value of the cooling rate
 */
double CoolingTable::get_value(double Fe, double Mg, double z, double n,
                               double T) {
    // can easily be switched to other method when necessary without changing
    // other files
    // currently hardcoded, change this later
    if(Fe < -99.) {
        cerr << "Coolingtable received Fe/H lower than minimum" << endl;
        my_exit();
    }

    // Sanity checks on temperature value
    // There is no cooling for T = 0
    if(!T) {
        return 0.;
    }
    if(T < 0.) {
        cerr << "Coolingtable received negative temperature" << endl;
        my_exit();
    }
    // we should linearly extrapolate for larger values...
    T = std::min(T, _T_max);

    // Sanity checks on redshift value
    if(z < 0.) {
        cerr << "Coolingtable received negative redshift" << endl;
        my_exit();
    }
    // this is correct: we want to use the z = 11 cooling tables for higher
    // redshifts; all of these have no UVB and are hence redshift independent
    z = std::min(z, _z_max);

    // actual interpolation
    double result;
    // first we try to find the hydrogen density...
    // no need to interpolate 3D again if same Fe, Mg and n as last time
    if(Fe != _Fe_prev || Mg != _Mg_prev || n != _n_prev) {
        _nH_prev = _ntable->get_value(vector<double>({Fe, Mg, n}));
    }
    // ...then we interpolate on the full cooling table
    result = _coolingtable->get_value(vector<double>({Fe, Mg, z, _nH_prev, T}));
    _Mg_prev = Mg;
    _Fe_prev = Fe;
    _n_prev = n;
    return result;
}

/**
 * @brief Dump the object to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void CoolingTable::dump(RestartFile& rfile) {
    _coolingtable->dump(rfile);
    _ntable->dump(rfile);
    rfile.write(_n_prev);
    rfile.write(_Fe_prev);
    rfile.write(_Mg_prev);
    rfile.write(_nH_prev);
    rfile.write(_T_max);
    rfile.write(_z_max);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
CoolingTable::CoolingTable(RestartFile& rfile) {
    _coolingtable = new FiveDTable(rfile);
    _ntable = new ThreeDIrregTable(rfile);
    rfile.read(_n_prev);
    rfile.read(_Fe_prev);
    rfile.read(_Mg_prev);
    rfile.read(_nH_prev);
    rfile.read(_T_max);
    rfile.read(_z_max);
}
