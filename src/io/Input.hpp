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
 * @file Input.hpp
 *
 * @brief Old interface for snapshot readers: header
 *
 * This interface is still used by ShadowfaxSnapshotReader.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef INPUT_HPP
#define INPUT_HPP

#include <hdf5.h>
#include <string>
#include <fstream>
#include <vector>

class Block;
class Header;

/**
  * \brief Interface for data input
  *
  * Data is organized in blocks which basically are tables with column names and
  * rows of data. Extra information is written to the file through a Header.
  */
class Input{
public:
    /**
      * \brief Fill the given Block with data
      *
      * The column names have to be pre-specified in the header of the Block and
      * only those data are accessible. Since many input methods gain
      * performance by pre-allocating memory, the size of the data that has to
      * be filled can be specified.
      *
      * @param block Block to be filled
      * @param npart Number of particles to read information for
      */
    virtual void read(Block& block, unsigned int npart)=0;

    /**
      * \brief Fill the given Header with data
      *
      * Contrary to the function that reads a Block, the implementation itself
      * determines what data are actually read in.
      *
      * @param header Header to fill
      */
    virtual void read_header(Header& header)=0;
};

/**
  * \brief Read in data from the default HDF5 snapshot format
  */
class FileInput : public Input{
private:
    /*! \brief Name of the HDF5 file that is being read */
    std::string _filename;

public:
    FileInput(std::string filename);
    ~FileInput(){}

    void read(Block& block, unsigned int npart);
    void read_header(Header &header);
};

/**
  * \brief Read in data from an ASCII file
  *
  * The data is expected to be in columns. The first line of the file can
  * optionally contain column names that are linked to the data. No Header
  * information can be present in the file. All data is read in during
  * construction and is stored in internal vectors.
  */
class AsciiInput : public Input{
private:
    /*! \brief Internal data vector */
    std::vector< std::vector<double> > _data;

    /*! \brief Internal vector of column names */
    std::vector<std::string> _column_names;

    void read_column_names(std::string line);

public:
    AsciiInput(std::string filename);
    ~AsciiInput(){}

    void read(Block &block, unsigned int npart);
    void read_header(Header &header);
};

#endif // INPUT_HPP
