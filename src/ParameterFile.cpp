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
#include "Error.hpp"        // for my_exit
#include "RestartFile.hpp"  // for RestartFile
#include "YMLFile.hpp"
#include <boost/property_tree/ini_parser.hpp>  // for read_ini
#include <fstream>   // for operator<<, basic_ostream, etc
#include <iostream>  // for cout, cerr
#include <string>
#include <utility>  // for pair
using namespace std;

/**
 * @brief Print the contents of the parameter file to the stdout
 */
void ParameterFile::print_contents() {
    if(_yml_file) {
        _yml_file->print_contents(cout);
    } else {
        for(auto it = _parameters.begin(); it != _parameters.end(); ++it) {
            cout << "### " << it->first.data() << " ###" << endl;
            for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                cout << it2->first.data() << ": " << it2->second.data() << endl;
            }
            cout << endl;
        }
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
ParameterFile::ParameterFile(std::string name) : _yml_file(nullptr) {
    ifstream file(name.c_str());
    if(!file) {
        cerr << "Cannot read parameterfile \"" << name << "\"!" << endl;
        my_exit();
    }

    // check type
    std::string extension = name.substr(name.size() - 3, 3);
    if(extension == "ini") {
        cout << ".ini parameter file" << endl;
        boost::property_tree::ini_parser::read_ini(name, _parameters);
    } else if(extension == "yml") {
        cout << ".yml parameter file" << endl;
        _yml_file = new YMLFile(name);
    } else {
        cerr << "Unknown parameter file type: ." << extension << endl;
        my_exit();
    }

    // print out some info
    cout << "Read parameters from " << name << endl;

    print_contents();
}

/**
 * @brief Destructor.
 *
 * Free YMLFile (if necessary).
 */
ParameterFile::~ParameterFile() {
    if(_yml_file) {
        delete _yml_file;
    }
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
    int propertyflag = PARAMETERFILE_RESTART_PROPERTYFLAG;
    if(_yml_file) {
        std::string tag = "YML";
        rfile.write(tag);
        for(auto it = _yml_file->begin(); it != _yml_file->end(); ++it) {
            rfile.write(propertyflag);
            rfile.write(it.get_key());
            rfile.write(it.get_value());
        }
    } else {
        std::string tag = "INI";
        rfile.write(tag);
        int headerflag = PARAMETERFILE_RESTART_HEADERFLAG;
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
    }
    int endflag = PARAMETERFILE_RESTART_ENDFLAG;
    rfile.write(endflag);
}

/**
 * @brief Restart constructor. Initialize the restartfile from the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 */
ParameterFile::ParameterFile(RestartFile& rfile) {
    std::string tag;
    rfile.read(tag);
    if(tag == "YML") {
        _yml_file = new YMLFile();
        int val;
        rfile.read(val);
        while(val != PARAMETERFILE_RESTART_ENDFLAG) {
            std::string key;
            std::string value;
            rfile.read(key);
            rfile.read(value);
            _yml_file->add_value(key, value);
            rfile.read(val);
        }
    } else {
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
    }

    print_contents();
}
