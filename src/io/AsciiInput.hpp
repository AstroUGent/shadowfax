/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file AsciiInput.hpp
 *
 * @brief Old interface for snapshot readers: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ASCIIINPUT_HPP
#define ASCIIINPUT_HPP

#include "Input.hpp"
#include <string>  // for string
#include <vector>  // for vector

class Block;
class Header;

/**
  * @brief Read in data from an ASCII file
  *
  * The data is expected to be in columns. The first line of the file can
  * optionally contain column names that are linked to the data. No Header
  * information can be present in the file. All data is read in during
  * construction and is stored in internal vectors.
  */
class AsciiInput : public Input {
  private:
    /*! @brief Internal data vector */
    std::vector<std::vector<double> > _data;

    /*! @brief Internal vector of column names */
    std::vector<std::string> _column_names;

    void read_column_names(std::string line);

  public:
    AsciiInput(std::string filename);
    ~AsciiInput() {}

    void read(Block& block, unsigned int npart);
    void read_header(Header& header);
};

#endif  // ASCIIINPUT_HPP
