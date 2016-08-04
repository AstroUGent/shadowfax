/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *                    Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 *               2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ParameterFile.cpp
 *
 * @brief Parameter file: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ParameterFile.hpp"
#include "Error.hpp"                           // for my_exit
#include "RestartFile.hpp"                     // for RestartFile
#include <boost/property_tree/ini_parser.hpp>  // for read_ini
#include <fstream>   // for operator<<, basic_ostream, etc
#include <iostream>  // for cout, cerr
#include <utility>   // for pair
using namespace std;

/**
 * @brief Print the contents of the parameter file to the stdout
 */
void ParameterFile::print_contents() {
    for(auto it = _parameters.begin(); it != _parameters.end(); ++it) {
        cout << "### " << it->first.data() << " ###" << endl;
        for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            cout << it2->first.data() << ": " << it2->second.data() << endl;
        }
        cout << endl;
    }
}

/**
 * @brief Constructor
 *
 * We try to open the file with the given name and initialize a
 * boost::property_tree based on its contents, assuming it is written in
 * .ini-format.
 *
 * We then parse all parameters and initialize them internally. Except for the
 * total simulation time, which needs to be specified, all parameters have
 * reasonable default values and can be omitted from the parameter file.
 *
 * @param name Filename of the parameter file
 */
ParameterFile::ParameterFile(std::string name) {
    ifstream file(name.c_str());
    if(!file) {
        cerr << "Cannot read parameterfile \"" << name << "\"!" << endl;
        my_exit();
    }

    boost::property_tree::ini_parser::read_ini(name, _parameters);

    // print out some info
    cout << "Read parameters from " << name << endl;

    print_contents();
}

#define PARAMETERFILE_RESTART_HEADERFLAG 0
#define PARAMETERFILE_RESTART_PROPERTYFLAG 1
#define PARAMETERFILE_RESTART_ENDFLAG 2

/**
 * @brief Dump the parameterfile to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void ParameterFile::dump(RestartFile& rfile) {
    int headerflag = PARAMETERFILE_RESTART_HEADERFLAG;
    int propertyflag = PARAMETERFILE_RESTART_PROPERTYFLAG;
    int endflag = PARAMETERFILE_RESTART_ENDFLAG;
    string header;
    string property;
    string value;
    for(auto it = _parameters.begin(); it != _parameters.end(); ++it) {
        rfile.write(headerflag);
        header = it->first.data();
        rfile.write(header);
        for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            rfile.write(propertyflag);
            property = it2->first.data();
            rfile.write(property);
            value = it2->second.data();
            rfile.write(value);
        }
    }
    rfile.write(endflag);
}

/**
 * @brief Restart constructor. Initialize the restartfile from the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 */
ParameterFile::ParameterFile(RestartFile& rfile) {
    int val;
    rfile.read(val);
    string header;
    while(val != PARAMETERFILE_RESTART_ENDFLAG) {
        if(val == PARAMETERFILE_RESTART_HEADERFLAG) {
            // new header
            rfile.read(header);
        } else {
            string name;
            string value;
            rfile.read(name);
            rfile.read(value);
            _parameters.put(header + string(".") + name, value);
        }
        rfile.read(val);
    }

    print_contents();
}
