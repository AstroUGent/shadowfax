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
 * @file FileOutput.cpp
 *
 * @brief Old interface for snapshot writers: implementation
 *
 * This interface is still used by ShadowfaxSnapshotWriter.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Output.hpp"
#include "Block.hpp"
#include "Header.hpp"
#include "HDF5tools.hpp"
#include <vector>
#include <iostream>
#include <sstream>
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include <hdf5.h>
using namespace std;

/**
  * \brief Construct a FileOutput with the given name
  *
  * We check if the file can be created and in the process override a possibly
  * existing file with the same name.
  *
  * @param filename Name of the file (should end with .hdf5)
  */
FileOutput::FileOutput(string filename){
    int rank = MPIGlobal::local_rank;
    // check if we can open the file (and immediately overwrite it if it already
    // exists)
    if(!rank){
        hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                               H5P_DEFAULT);
        herr_t status;
        status = H5Fclose(file);

        if(status < 0){
            std::cerr << "ERROR!" << std::endl;
        }
    }
    _filename = filename;
}

/**
  * \brief Write a Block to the HDF5 file.
  *
  * The corresponding dataset in the file is first created by the process with
  * rank 0. Since we do not use parallel HDF5, we have to close and reopen the
  * file for every process. This induces some overhead, but there seems to be no
  *  other option...
  *
  * @param block The Block the write to the file
  */
void FileOutput::write(Block& block){
    unsigned int blocksize = block.get_size();
    int rank = MPIGlobal::rank;
    int size = MPIGlobal::size;
    int lrank = MPIGlobal::local_rank;
    vector<unsigned int> blocksizes(size);
    MyMPI_Allgather(&blocksize, 1, MPI_UNSIGNED, &blocksizes[0], 1,
            MPI_UNSIGNED);
    for(int i = 0; i < size; i++){
        if(i == lrank){
            if(!lrank){
                blocksize = 0;
                for(unsigned int j = size; j--;){
                    blocksize += blocksizes[j];
                }
                hid_t file = H5Fopen(_filename.c_str(), H5F_ACC_RDWR,
                                     H5P_DEFAULT);
                herr_t status;
                hid_t group = H5Gcreate(file, block.get_name().c_str(), -1);
                for(unsigned int k = 0; k < block.get_headers().size(); k++){
                    unsigned int rank = 1;
                    if(block.get_dimensions()[k] > 1){
                        rank = 2;
                    }
                    hsize_t dims[2] = {blocksize, block.get_dimensions()[k]};
                    stringstream name;
                    name << "/" << block.get_name() << "/"
                         << block.get_headers()[k];
                    hid_t filespace = H5Screate_simple(rank, dims, NULL);
                    hid_t dataset = H5Dcreate(group, name.str().c_str(),
                                              H5T_NATIVE_DOUBLE, filespace,
                                              H5P_DEFAULT);
                    Unit& unit = block.get_unit(k);

                    HDF5tools::write_attribute_string(dataset, "unit_quantity",
                                                      unit.quantity());
                    HDF5tools::write_attribute_string(dataset, "unit_name",
                                                      unit.name());
                    HDF5tools::write_attribute_scalar(dataset, "unit_SI_value",
                                                      HDF5types::DOUBLE,
                                                      unit.SI_value());

                    hsize_t dims2[2] = {blocksizes[0],
                                        block.get_dimensions()[k]};
                    hsize_t offset[2] = {0, 0};
                    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                                 offset, NULL, dims2, NULL);
                    hid_t memspace = H5Screate_simple(rank, dims2, NULL);
                    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace,
                                      filespace, H5P_DEFAULT,
                                      block.get_buffer(k));

                    status = H5Sclose(memspace);
                    status = H5Sclose(filespace);
                    status = H5Dclose(dataset);
                }
                status = H5Gclose(group);
                status = H5Fclose(file);

                if(status < 0){
                    std::cerr << "ERROR!" << std::endl;
                }
            } else {
                hid_t file = H5Fopen(_filename.c_str(), H5F_ACC_RDWR,
                                     H5P_DEFAULT);
                herr_t status;
                hid_t group = H5Gopen(file, block.get_name().c_str());
                blocksize = 0;
                for(unsigned int j = rank; j--;){
                    blocksize += blocksizes[j];
                }
                for(unsigned int k = 0; k < block.get_headers().size(); k++){
                    unsigned int rank = 1;
                    if(block.get_dimensions()[k] > 1){
                        rank = 2;
                    }
                    stringstream name;
                    name << "/" << block.get_name() << "/"
                         << block.get_headers()[k];
                    hid_t dataset = H5Dopen(file, name.str().c_str());
                    hid_t filespace = H5Dget_space(dataset);
                    hsize_t dims2[2] = {blocksizes[i],
                                        block.get_dimensions()[k]};
                    hsize_t offset[2] = {blocksize, 0};
                    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                                 offset, NULL, dims2, NULL);
                    hid_t memspace = H5Screate_simple(rank, dims2, NULL);
                    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace,
                                      filespace, H5P_DEFAULT,
                                      block.get_buffer(k));

                    status = H5Sclose(memspace);
                    status = H5Sclose(filespace);
                    status = H5Dclose(dataset);
                }
                status = H5Gclose(group);
                status = H5Fclose(file);

                if(status < 0){
                    std::cerr << "ERROR!" << std::endl;
                }
            }
        }
        // we have to wait for the previous process to close the file (because
        // opening a file that is already opened does result in errors)
        MyMPI_Barrier();
    }
}

/**
  * \brief Write a Header to the file
  *
  * Only the process with rank 0 writes the header.
  *
  * @param header The Header to write to the file
  */
void FileOutput::write_header(Header &header){
    int rank = MPIGlobal::local_rank;
    // only process 0 writes the header
    if(!rank){
        hid_t file = H5Fopen(_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        herr_t status;
        hid_t group = H5Gcreate(file, "Header", -1);

        HDF5tools::write_attribute_array(group, 2, "npart", HDF5types::UINT,
                                         header.get_npartspec());
        HDF5tools::write_attribute_scalar(group, "ndim", HDF5types::UINT,
                                          header.get_ndim());
        HDF5tools::write_attribute_array(group, ndim_+ndim_, "box",
                                         HDF5types::DOUBLE, header.get_box());
        HDF5tools::write_attribute_scalar(group, "time", HDF5types::DOUBLE,
                                          header.get_time());
        HDF5tools::write_attribute_scalar(group, "periodic", HDF5types::BOOL,
                                          header.get_periodic());
        HDF5tools::write_attribute_scalar(group, "second order",
                                          HDF5types::BOOL,
                                          header.get_second_order());
        HDF5tools::write_attribute_scalar(group, "static", HDF5types::BOOL,
                                          header.get_static());
        HDF5tools::write_attribute_scalar(group, "global timestep",
                                          HDF5types::BOOL,
                                          header.get_global_timestep());
        HDF5tools::write_attribute_scalar(group, "adiabatic index",
                                          HDF5types::DOUBLE,
                                          header.get_gamma());
        HDF5tools::write_attribute_scalar(group, "gravity", HDF5types::BOOL,
                                          header.get_gravity());
        if(header.ndmpart() && header.gravity()){
            HDF5tools::write_attribute_scalar(group, "hsoft", HDF5types::DOUBLE,
                                              header.get_hsoft());
        }

        status = H5Gclose(group);
        status = H5Fclose(file);

        if(status < 0){
            std::cerr << "ERROR!" << std::endl;
        }
    }
}
