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
 * @file Output.hpp
 *
 * @brief Old interface for snapshot writers: header
 *
 * This interface is still used by ShadowfaxSnapshotWriter.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <string>
#include <fstream>
#include <vector>

class Block;
class Header;

/**
  * \brief Interface for data output
  *
  * Data is organized in blocks, which are a kind of tables with column names
  * and rows of data. Extra information is provided through a Header.
  */
class Output{
public:
    /**
      * \brief Output a Block
      *
      * All information in the header of the Block and all actual data is
      * guaranteed to be outputted. Extra information (e.g. units) may or may
      * not be outputted depending on the implementation.
      *
      * @param block Block to write out
      */
    virtual void write(Block& block)=0;

    /**
      * \brief Output a Header
      *
      * What information is outputted is implementation dependant.
      *
      * @param header Header to write out
      */
    virtual void write_header(Header& header)=0;
};

/**
  * \brief Output a default HDF5 snapshot
  *
  * Every Block written to the FileOutput corresponds to a separate group in the
  *  HDF5 file. Header information is written to a separate group called Header.
  * Columns correspond to separate datasets inside a group and Unit information
  * is provided for every dataset through HDF5 attributes.
  */
class FileOutput : public Output{
private:
    /*! \brief Name of the HDF5 file (has to end with .hdf5) */
    std::string _filename;

public:
    FileOutput(std::string filename);
    ~FileOutput(){}
    void write(Block& block);
    void write_header(Header &header);
};

/**
  * \brief Output data to an ASCII file
  *
  * All blocks are written to one large table in the order they are written to
  * the AsciiOutput. Unit information is discarded and the Header is ignored.
  * Column names are written to the first line, which is preceded by a #. All
  * data is stored internally and the actual file is written by the
  * deconstructor.
  */
class AsciiOutput : public Output{
private:
    /*! \brief Internal data vector */
    std::vector< std::vector<double> > _data;

    /*! \brief Internal vector with column names */
    std::vector<std::string> _column_names;

    /*! \brief Name of the output file */
    std::string _filename;

    void write_column_names(std::ofstream& file);

public:
    AsciiOutput(std::string filename);
    ~AsciiOutput();

    void write(Block &block);
    void write_header(Header &header);
};

#endif // OUTPUT_HPP
