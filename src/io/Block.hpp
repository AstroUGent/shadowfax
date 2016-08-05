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
 * @file Block.hpp
 *
 * @brief Data block for old input/output: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "Unit.hpp"
#include <string>
#include <vector>

/**
  * @brief A block of data to be written to or read from a snapshot file
  *
  * A Block represents a chunk of data that can be written to a snapshot file
  * or read from it. The data is stored in columns, every column has a header
  * which holds the name of the data in that column.
  * The Block stores its data in an internal std::vector<std::vector>. The
  * internal data can be accessed through a void pointer, as is customary in
  * standard c writing and reading routines.
  */
class Block {
  private:
    /*! @brief The name of the Block */
    std::string _name;

    /*! @brief The header list of the Block */
    std::vector<std::string> _headers;

    /*! @brief The dimensions list of the Block (the dimensions of a given
     *  column are the number of values inside one cell of that column) */
    std::vector<unsigned int> _dimensions;

    /*! @brief The internal data buffer of the Block */
    std::vector<std::vector<double> > _data;

    /*! @brief The index of the first uninitialized row in the internal data
     *  buffer. If no uninitialized rows exist, this is equal to the column
     *  length */
    unsigned int _emptyline;

    /*! @brief A vector with a Unit for every column in the Block */
    std::vector<Unit> _units;

  public:
    Block(std::string name, std::vector<std::string>& headers,
          std::vector<unsigned int>& dimensions, std::vector<Unit>& units,
          unsigned int size = 0);
    void add_data(std::vector<std::vector<double> >& data);
    void add_data(std::vector<double>& data);
    void add_data(unsigned int x, unsigned int y, double data);
    void add_column(std::vector<double>& data, unsigned int y);
    void add_row(std::vector<double>& data);

    const std::vector<std::vector<double> >& get_data();
    const void* get_buffer(unsigned int index);
    Unit& get_unit(unsigned int index);
    const std::string& get_name();
    const std::vector<std::string>& get_headers();
    const std::vector<unsigned int>& get_dimensions();

    unsigned int number_of_lines();
    const std::vector<double> get_line(unsigned int index);

    unsigned int get_size();
};

#endif  // BLOCK_HPP
