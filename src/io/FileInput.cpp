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
 * @file FileInput.cpp
 *
 * @brief Old interface for snapshot readers: implementation
 *
 * This interface is still used by ShadowfaxSnapshotReader.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "FileInput.hpp"
#include "Block.hpp"          // for Block
#include "HDF5tools.hpp"      // for read_attribute, BOOL, etc
#include "Header.hpp"         // for Header
#include "MPIGlobal.hpp"      // for rank, size
#include "MPIMethods.hpp"     // for MyMPI_Bcast, MyMPI_Barrier
#include "Unit.hpp"           // for Unit
#include "UnitConverter.hpp"  // for UnitConverter
#include <cstdlib>            // for NULL
#include <iostream>           // for operator<<, basic_ostream, etc
#include <sstream>
#include <string>  // for string, operator<<
#include <vector>  // for vector
using namespace std;

/**
  * @brief Construct a FileInput for the file with the given name
  *
  * No attempt is made to actually open the file and no check on the existence
  * of the file is present (yet?).
  *
  * @param filename Name of the file (should end with .hdf5)
  */
FileInput::FileInput(string filename) {
    _filename = filename;
}

/**
  * @brief Read information from a HDF5-file and store the contents in the
  * correct fields in the given Block.
  *
  * The Block initially contains a list with a Unit for every variable. The
  * variables in the snapshot can also contain a Unit. When the variable is read
  * in, a check is performed on compatibility of the units. If they are
  * compatible, the variable is converted from the snapshot Unit to the
  * requested Unit in the Block. If the snapshot does not contain units, we
  * assume that the snapshot units and the requested units are equal.
  *
  * The C++ HDF5-API does not offer functionality to check if attributes exist.
  * We therefore use the C API, which does offer this functionality.
  *
  * @param block Block to fill with data
  * @param npart Number of lines that should be read
  */
void FileInput::read(Block& block, unsigned int npart) {
#ifdef PYTHON_MODULE
    int rank = 0;
    int size = 1;
#else
    int rank = MPIGlobal::rank;
    int size = MPIGlobal::size;
#endif
    // if npart is a multiple of world.size(), we want nice blocks. If it is
    // not, we should make sure that all particles are read
    unsigned int npart_local = npart / size + ((npart % size) > 0);
    unsigned int npart_other = rank * npart_local;
    while(npart_other + npart_local > npart) {
        npart_local--;
    }
    for(int i = 0; i < size; i++) {
        if(i == rank) {
            hid_t file =
                    H5Fopen(_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            herr_t status;
            hid_t group = H5Gopen(file, block.get_name().c_str());
            for(unsigned int j = 0; j < block.get_headers().size(); j++) {
                stringstream name;
                name << "/" << block.get_name() << "/"
                     << block.get_headers()[j];
                hid_t dataset = H5Dopen(file, name.str().c_str());
                Unit unit;
                if(H5Aexists(dataset, "unit_quantity")) {
                    string unit_quantity = HDF5tools::read_attribute_string(
                            dataset, "unit_quantity");
                    string unit_name = HDF5tools::read_attribute_string(
                            dataset, "unit_name");
                    HDF5tools::read_attribute(dataset, "unit_SI_value",
                                              HDF5types::DOUBLE,
                                              unit.SI_value());
                    unit = Unit(unit_quantity, unit_name, unit.get_SI_value());
                } else {
                    unit = block.get_unit(j);
                }
                UnitConverter uc(block.get_unit(j), unit);
                hid_t space = H5Dget_space(dataset);
                hsize_t dims2[1] = {npart_local};
                hsize_t offset[1] = {npart_other};
                // for some reason, we need to set up a memory space, otherwise
                // H5Sclose(space) results in a SEGFAULT on parallel systems...
                hid_t memspace = H5Screate_simple(1, dims2, NULL);
                vector<double> data(npart_local, 0.);
                status = H5Sselect_hyperslab(space, H5S_SELECT_SET, offset,
                                             NULL, dims2, NULL);
                status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, space,
                                 H5P_DEFAULT, &data[0]);
                status = H5Sclose(space);
                status = H5Sclose(memspace);
                status = H5Dclose(dataset);
                for(unsigned int k = data.size(); k--;) {
                    data[k] = uc.convert(data[k]);
                }
                block.add_column(data, j);
            }
            status = H5Gclose(group);
            status = H5Fclose(file);

            if(status < 0) {
                std::cerr << "ERROR!" << std::endl;
            }
        }
#ifndef PYTHON_MODULE
        MyMPI_Barrier();
#endif
    }
}

/**
  * @brief Fill the given Header with the header information of the file
  *
  * All information is read in and then a consistency check is performed to
  * check if the compiled version of the code is compatible with the information
  * in the Header.
  *
  * @param header Header to fill
  */
void FileInput::read_header(Header& header) {
#ifdef PYTHON_MODULE
    int rank = 0;
#else
    int rank = MPIGlobal::rank;
#endif
    if(!rank) {
        hid_t flag = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fclose_degree(flag, H5F_CLOSE_SEMI);
        hid_t file = H5Fopen(_filename.c_str(), H5F_ACC_RDONLY, flag);
        herr_t status;
        hid_t group = H5Gopen(file, "/Header");

        HDF5tools::read_attribute(group, "npart", HDF5types::UINT,
                                  header.get_npartspec());
        HDF5tools::read_attribute(group, "ndim", HDF5types::UINT,
                                  header.get_ndim());
        HDF5tools::read_attribute(group, "box", HDF5types::DOUBLE,
                                  header.get_box());
        HDF5tools::read_attribute(group, "time", HDF5types::DOUBLE,
                                  header.get_time());
        HDF5tools::read_attribute(group, "periodic", HDF5types::BOOL,
                                  header.get_periodic());
        HDF5tools::read_attribute(group, "second order", HDF5types::BOOL,
                                  header.get_second_order());
        HDF5tools::read_attribute(group, "static", HDF5types::BOOL,
                                  header.get_static());
        HDF5tools::read_attribute(group, "global timestep", HDF5types::BOOL,
                                  header.get_global_timestep());
        HDF5tools::read_attribute(group, "adiabatic index", HDF5types::DOUBLE,
                                  header.get_gamma());
        HDF5tools::read_attribute(group, "gravity", HDF5types::BOOL,
                                  header.get_gravity());
        if(header.ndmpart() && header.gravity()) {
            HDF5tools::read_attribute(group, "hsoft", HDF5types::DOUBLE,
                                      header.get_hsoft());
        }

        status = H5Gclose(group);
        status = H5Fclose(file);

        if(status < 0) {
            std::cerr << "ERROR!" << std::endl;
        }

        header.check_makeflags();
    }
#ifndef PYTHON_MODULE
    // we do not have the MPIGlobal buffer yet at this point,
    // so we have to provide our own buffer
    int bufsize = sizeof(Header);
    char* buffer = new char[bufsize];
    int send_pos = 0;
    if(!rank) {
        header.pack_data(buffer, bufsize, &send_pos);
    }
    MyMPI_Bcast(&send_pos, 1, MPI_INT, 0);
    MyMPI_Bcast(buffer, send_pos, MPI_PACKED, 0);
    int recv_pos = 0;
    header = Header(buffer, send_pos, &recv_pos);
    delete[] buffer;
#endif
}
