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
 * @file FileInput.hpp
 *
 * @brief Old interface for snapshot readers: header
 *
 * This interface is still used by ShadowfaxSnapshotReader.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef FILEINPUT_HPP
#define FILEINPUT_HPP

#include "Input.hpp"
#include <string>  // for string

class Block;
class Header;

/**
  * @brief Read in data from the default HDF5 snapshot format
  */
class FileInput : public Input {
  private:
    /*! @brief Name of the HDF5 file that is being read */
    std::string _filename;

  public:
    FileInput(std::string filename);
    ~FileInput() {}

    void read(Block& block, unsigned int npart);
    void read_header(Header& header);
};

#endif  // FILEINPUT_HPP
