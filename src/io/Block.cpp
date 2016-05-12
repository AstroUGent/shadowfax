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
 * @file Block.cpp
 *
 * @brief Data block for old input/output: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Block.hpp"
#include "Error.hpp"
#include <cstdlib>
#include <iostream>
using namespace std;

/**
  * @brief Constructor
  *
  * Initialize column headers and data dimensions. Resize the internal
  * std::vector to the given size.
  *
  * @param name Name of the Block
  * @param headers Column headers for the data
  * @param dimensions Dimensions of the data columns (the number of values
  * stored in every cell of the column)
  * @param units Units of the data columns
  * @param size Size of the data: the number of cells in every column that will
  * be filled later on
  */
Block::Block(string name, vector<string>& headers,
             vector<unsigned int>& dimensions, vector<Unit>& units,
             unsigned int size) {
    _name = name;
    _data.resize(headers.size());
    for(unsigned int i = 0; i < headers.size(); i++) {
        _headers.push_back(headers[i]);
        _dimensions.push_back(dimensions[i]);
        _data[i].resize(size * dimensions[i]);
        _units.push_back(units[i]);
    }
    _emptyline = 0;
}

/**
  * @brief Fill the entire data buffer at once using the given
  * std::vector<std::vector>
  *
  * @param data A data buffer with dimensions equal to or smaller than the
  * dimensions of the internal buffer.
  */
void Block::add_data(std::vector<std::vector<double> >& data) {
    for(unsigned int i = data.size(); i--;) {
        for(unsigned int j = data[i].size(); j--;) {
            _data[i][j] = data[i][j];
        }
    }
}

/**
  * @brief Add a row of data to the end of the internal buffer
  *
  * The end of the internal buffer is
  *  - the first row with uninitialized data,
  * or, if no unitialized rows exist,
  *  - a newly created row at the end of the internal buffer
  *
  * @param data A row of data with the same dimensions as the number of columns
  * in the Block
  */
void Block::add_data(std::vector<double>& data) {
    if(_emptyline == _data[0].size()) {
        for(unsigned int i = _data.size(); i--;) {
            _data[i].push_back(0.);
        }
    }
    for(unsigned int i = _data.size(); i--;) {
        _data[i][_emptyline] = data[i];
    }
    _emptyline++;
}

/**
  * @brief Add a column of data to the Block
  *
  * @param data A column of data with unspecified dimensions
  * @param y The index of the column in the internal data buffer (or the index
  * of the column in the headers)
  */
void Block::add_column(vector<double>& data, unsigned int y) {
    _data[y] = data;
}

/**
  * @brief Add a row of data to the end of the internal buffer
  *
  * The end of the internal buffer is
  *  - the first row with uninitialized data,
  * or, if no unitialized rows exist,
  *  - a newly created row at the end of the internal buffer
  *
  * @param data A row of data with the same dimensions as the number of columns
  * in the Block (and taking into account multidimensional columns)
  */
void Block::add_row(vector<double>& data) {
    if(_emptyline == _data[0].size()) {
        for(unsigned int i = _data.size(); i--;) {
            for(unsigned int j = _dimensions[i]; j--;) {
                _data[i].push_back(0.);
            }
        }
    }
    unsigned int index = 0;
    for(unsigned int i = 0; i < _data.size(); i++) {
        for(unsigned int j = 0; j < _dimensions[i]; j++) {
            _data[i][_emptyline * _dimensions[i] + j] = data[index++];
        }
    }
    _emptyline++;
}

/**
  * @brief Add data to the specified cell in the internal data buffer
  *
  * @param x Index of the cell in the column (so a row number)
  * @param y Index of the column (so a column number)
  * @param data The data to add to the cell
  */
void Block::add_data(unsigned int x, unsigned int y, double data) {
    if(y >= _headers.size()) {
        // error!
        cout << "Error in Block.cpp: trying to add element in column that does "
                "not exist! (x="
             << x << ")" << endl;
        my_exit();
    }
    if(x >= _data[0].size()) {
        for(unsigned int i = _headers.size(); i--;) {
            _data[i].push_back(0.);
        }
    }
    _data[y][x] = data;
}

/**
  * @brief Access the internal data buffer
  *
  * @return A reference to the internal data buffer
  */
const vector<vector<double> >& Block::get_data() {
    return _data;
}

/**
  * @brief Access a column of the internal data buffer
  *
  * @param index The index of the column in the internal data buffer (or the
  * index of the column in the header list)
  * @return A c-type void pointer to the buffer underlying the specified column,
  * which can be used to read from or write to the column
  */
const void* Block::get_buffer(unsigned int index) {
    return &_data[index][0];
}

/**
  * @brief Get the Unit associated with a column in the Block
  *
  * @param index The index of the column in the header list
  * @return A reference to the unit assiociated to the given column
  */
Unit& Block::get_unit(unsigned int index) {
    return _units[index];
}

/**
  * @brief Access the name of the Block
  *
  * @return A reference to the name of the Block
  */
const string& Block::get_name() {
    return _name;
}

/**
  * @brief Access the header list
  *
  * @return A reference to the header list of the Block
  */
const vector<string>& Block::get_headers() {
    return _headers;
}

/**
  * @brief Access the dimensions list. The dimension of a given column is the
  * number of data values inside one cell of the column
  *
  * @return A reference to the dimensions list of the Block
  */
const vector<unsigned int>& Block::get_dimensions() {
    return _dimensions;
}

/**
  * @brief Get information on the number of cells in the Block
  *
  * @return The number of cells in one column of the Block
  */
unsigned int Block::number_of_lines() {
    return _data[0].size();
}

/**
  * @brief Access a row of the Block
  *
  * @param index The valid index of a row in the internal buffer
  * @return A std::vector containing the data in the specified row of the Block
  */
const vector<double> Block::get_line(unsigned int index) {
    vector<double> line(_headers.size());
    for(unsigned int i = _data.size(); i--;) {
        line[i] = _data[i][index];
    }
    return line;
}

/**
  * @brief Get information on the number of cells in the Block
  *
  * @return The number of cells in a single column of the Block, taking into
  * account the number of dimensions of the columns
  */
unsigned int Block::get_size() {
    return _data[0].size() / _dimensions[0];
}
