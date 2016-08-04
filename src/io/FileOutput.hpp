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
 * @file FileOutput.hpp
 *
 * @brief Old interface for snapshot writers: header
 *
 * This interface is still used by ShadowfaxSnapshotWriter.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef FILEOUTPUT_HPP
#define FILEOUTPUT_HPP

#include "Output.hpp"
#include <string>

class Block;
class Header;

/**
  * @brief Output a default HDF5 snapshot
  *
  * Every Block written to the FileOutput corresponds to a separate group in the
  *  HDF5 file. Header information is written to a separate group called Header.
  * Columns correspond to separate datasets inside a group and Unit information
  * is provided for every dataset through HDF5 attributes.
  */
class FileOutput : public Output {
  private:
    /*! @brief Name of the HDF5 file (has to end with .hdf5) */
    std::string _filename;

  public:
    FileOutput(std::string filename);
    ~FileOutput() {}
    void write(Block& block);
    void write_header(Header& header);
};

#endif  // OUTPUT_HPP
