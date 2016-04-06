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
 * @file AsciiOutput.cpp
 *
 * @brief ASCII file writer
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Block.hpp"
#include "Header.hpp"
#include "Output.hpp"
using namespace std;

/**
  * \brief Construct an AsciiOutput to write a file with the given name
  *
  * @param filename Name of the file to write. The file is overwritten if it
  * already exists
  */
AsciiOutput::AsciiOutput(string filename) { _filename = filename; }

/**
  * \brief Destructor. Does the actual writing of the file
  *
  * We first write the column names to a comment line and then we write the
  * data.
  */
AsciiOutput::~AsciiOutput() {
    ofstream file(_filename.c_str());
    write_column_names(file);
    for(unsigned int i = 0; i < _data.size(); i++) {
        for(unsigned int j = 0; j < _data[i].size(); j++) {
            file << _data[i][j];
            if(j < _data[i].size() - 1) {
                file << "\t";
            } else {
                file << "\n";
            }
        }
    }
}

/**
  * \brief Write the column names to the file
  *
  * The line starts with a # and column names are separated by a tab.
  *
  * @param file ofstream to write to
  */
void AsciiOutput::write_column_names(ofstream& file) {
    file << "# ";
    for(unsigned int i = 0; i < _column_names.size(); i++) {
        file << _column_names[i];
        if(i < _column_names.size() - 1) {
            file << "\t";
        } else {
            file << "\n";
        }
    }
}

/**
  * \brief Write the given Header to the file
  *
  * This function is only provided to implement the Output interface and does
  * not do anything (yet?)
  *
  * @param header Header to write to the file
  */
void AsciiOutput::write_header(Header& header) {
    // do nothing, since an ascii-file does not contain a header (yet?)
}

/**
  * \brief Write the given Block to the file
  *
  * The data is actually written to internal vectors and only written out in the
  * end (in the destructor). Column names in the header of the Block are
  * appended to the internal vector of column names and the data are appended in
  * the same order to every line of the internal data vector.
  *
  * @param block Block to write to the file
  */
void AsciiOutput::write(Block& block) {
    vector<string> headers = block.get_headers();
    for(unsigned int i = 0; i < headers.size(); i++) {
        _column_names.push_back(headers[i]);
    }
    if(_data.size() < block.number_of_lines()) {
        _data.resize(block.number_of_lines());
    }
    for(unsigned int i = 0; i < block.number_of_lines(); i++) {
        vector<double> line = block.get_line(i);
        for(unsigned int j = 0; j < line.size(); j++) {
            _data[i].push_back(line[j]);
        }
    }
}
