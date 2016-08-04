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
 * @file AsciiOutput.hpp
 *
 * @brief Old interface for snapshot writers: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ASCIIOUTPUT_HPP
#define ASCIIOUTPUT_HPP

#include "Output.hpp"
#include <fstream>
#include <string>
#include <vector>

class Block;
class Header;

/**
  * @brief Output data to an ASCII file
  *
  * All blocks are written to one large table in the order they are written to
  * the AsciiOutput. Unit information is discarded and the Header is ignored.
  * Column names are written to the first line, which is preceded by a #. All
  * data is stored internally and the actual file is written by the
  * deconstructor.
  */
class AsciiOutput : public Output {
  private:
    /*! @brief Internal data vector */
    std::vector<std::vector<double> > _data;

    /*! @brief Internal vector with column names */
    std::vector<std::string> _column_names;

    /*! @brief Name of the output file */
    std::string _filename;

    void write_column_names(std::ofstream& file);

  public:
    AsciiOutput(std::string filename);
    ~AsciiOutput();

    void write(Block& block);
    void write_header(Header& header);
};

#endif  // ASCIIOUTPUT_HPP
