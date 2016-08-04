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
 * @file AsciiInput.cpp
 *
 * @brief ASCII file reader
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AsciiInput.hpp"
#include "Block.hpp"   // for Block
#include "Error.hpp"   // for my_exit
#include "Header.hpp"  // for Header
#include "sstream"
#include <iostream>  // for operator<<, basic_istream, etc
#include <string>    // for string, operator<<, getline, etc
#include <vector>    // for vector
using namespace std;

/**
  * @brief Construct an AsciiInput that reads the file with the given name
  *
  * The file is opened and the data is written to an internal vector. If
  * present, the names of the columns are read in from a commented out line.
  *
  * @param filename Name of the ASCII file
  */
AsciiInput::AsciiInput(std::string filename) {
    ifstream is(filename.c_str());
    if(!is) {
        cerr << "Error! Cannot open file \"" << filename << "\"!" << endl;
        my_exit();
    }
    string line;
    double val;
    while(getline(is, line)) {
        if(line[0] == '#') {
            line.erase(0, 1);
            read_column_names(line);
            continue;
        }
        vector<double> row;
        istringstream linestream(line);
        while(linestream >> val) {
            row.push_back(val);
        }
        _data.push_back(row);
    }
}

/**
  * @brief Read a block from the file
  *
  * Every name in the header of the Block is looked up in the internal list of
  * column names and the requested columns are then written line by line to the
  * Block. If a column is requested that is not part of the file, this results
  * in an error. If no column names were present in the file, we only check the
  * number of columns and it is assumed that the columns have the same order as
  * the elements in the header of the Block
  *
  * @param block Block to fill with data
  * @param npart Number of particles to read (actually not used in this
  * function)
  */
void AsciiInput::read(Block& block, unsigned int npart) {
    if(_column_names.size()) {
        vector<string> headers = block.get_headers();
        vector<unsigned int> indices(headers.size());
        for(unsigned int i = headers.size(); i--;) {
            unsigned int j = 0;
            while(j < _column_names.size() &&
                  headers[i].compare(_column_names[j])) {
                j++;
            }
            if(j == _column_names.size()) {
                cerr << "Requested column '" << headers[i]
                     << "' which is not "
                        "part of the Ascii file"
                     << endl;
                my_exit();
            }
            indices[i] = j;
        }
        for(unsigned int i = 0; i < _data.size(); i++) {
            vector<double> data(indices.size());
            for(unsigned int j = indices.size(); j--;) {
                data[j] = _data[i][indices[j]];
            }
            block.add_data(data);
        }
    } else {
        if(_data[0].size() != block.get_headers().size()) {
            cerr << "Number of colums in ASCII file is incompatible with "
                    "requested number of columns. Error!"
                 << endl;
            my_exit();
        }
        for(unsigned int i = 0; i < _data.size(); i++) {
            block.add_data(_data[i]);
        }
    }
}

/**
  * @brief Read the header of the file
  *
  * For now, we only set the number of gas particles, since no other data can be
  * present in the file.
  *
  * @param header The Header to fill
  */
void AsciiInput::read_header(Header& header) {
    header.set_ngaspart(_data.size());
}

/**
  * @brief Read column names from a string
  *
  * The names are expected to be single strings (without internal whitespace)
  * and separated by whitespace (spaces or tabs)
  *
  * @param line String containing column names separated by whitespace
  */
void AsciiInput::read_column_names(string line) {
    istringstream row(line);
    string name;
    while(row >> name) {
        _column_names.push_back(name);
    }
}
