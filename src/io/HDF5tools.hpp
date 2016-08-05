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
 * @file HDF5tools.hpp
 *
 * @brief Auxiliary functions to work with the HDF5 library
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HDF5TOOLS_HPP
#define HDF5TOOLS_HPP

#include "Error.hpp"
#include <hdf5.h>
#include <iostream>  // for operator<<, basic_ostream, etc
#include <stddef.h>  // for NULL
#include <string>    // for operator<<, char_traits, etc
#include <vector>    // for vector

/**
  * @brief Convenient aliases for the HDF5 datatypes
  */
namespace HDF5types {
/*! @brief Boolean representation: a native 32-bit unsigned integer */
static const hid_t BOOL __attribute__((used)) = H5T_NATIVE_UINT32;
/*! @brief Double precision floating point representation: a native 64-bit
 *  double */
static const hid_t DOUBLE __attribute__((used)) = H5T_NATIVE_DOUBLE;
/*! @brief Unsigned integer representation: a native 32-bit unsigned integer */
static const hid_t UINT __attribute__((used)) = H5T_NATIVE_UINT32;
/*! @brief Single precision floating point representation: a native 32-bit
 *  float */
static const hid_t FLOAT __attribute__((used)) = H5T_NATIVE_FLOAT;
/*! @brief Unsigned long representation: a native 64-bit unsigned integer */
static const hid_t ULONG __attribute__((used)) = H5T_NATIVE_UINT64;
/*! @brief Integer representation: a native 32-bit integer */
static const hid_t INT __attribute__((used)) = H5T_NATIVE_INT32;
}

/**
  * @brief Convenient aliases for HDF5 constants
  */
namespace HDF5constants {
/*! @brief The default property list used if properties are required */
static const hid_t DEFAULT_PROPERTY_LIST __attribute__((used)) = H5P_DEFAULT;
/*! @brief Overwrite existing files when creating a new HDF5 file */
static const hid_t OVERWRITE_EXISTING_FILES __attribute__((used)) =
        H5F_ACC_TRUNC;
/*! @brief Open a HDF5 file in read/write modus */
static const hid_t READWRITE __attribute__((used)) = H5F_ACC_RDWR;
/*! @brief Open a HDF5 file in read only modus */
static const hid_t READONLY __attribute__((used)) = H5F_ACC_RDONLY;
/*! @brief Use the entire dataspace of the dataset or attribute as file space or
 *  memory space during read/write operations */
static const hid_t COMPLETE_DATASPACE __attribute__((used)) = H5S_ALL;

/*! @brief Intialize a new selection when selecting hyperslabs of dataspaces */
static const H5S_seloper_t SET_NEW_SELECTION __attribute__((used)) =
        H5S_SELECT_SET;
/*! @brief Select a continuous block in the dataspace during hyperslab
 *  selection */
static const hsize_t* CONTINUOUS_SELECTION __attribute__((used)) = NULL;
/*! @brief Return single element blocks when selecting hyperslabs of
 *  dataspaces */
static const hsize_t* SINGLE_ELEMENT_BLOCKS __attribute__((used)) = NULL;
}

/**
  * @brief Simplified wrappers around low level HDF5 API functions
  */
namespace HDF5tools {
// if you remove any of the "inline" qualifiers, the code below will
// fail to compile!!

// write functions

/**
  * @brief Write a scalar attribute to the specified group or dataset
  *
  * The method creates the appropriate memory space and then creates the
  * attribute. Both the memory space and the attribute are also closed again.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the attribute
  * @param type HDF5 id for the type of the attribute (can be an alias from
  * HDF5types)
  * @param value C-style void pointer to the data that should be written to the
  * attribute
  */
inline void write_attribute_scalar(hid_t group, std::string name, hid_t type,
                                   void* value) {
    // create scalar dataspace
    hid_t attspace = H5Screate(H5S_SCALAR);
    if(attspace < 0) {
        std::cerr << "Error! Failed to create dataspace for scalar attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // create attribute
    hid_t att = H5Acreate(group, name.c_str(), type, attspace, H5P_DEFAULT);
    if(att < 0) {
        std::cerr << "Error! Failed to create scalar attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // write attribute
    herr_t status = H5Awrite(att, type, value);
    if(status < 0) {
        std::cerr << "Error! Failed to write scalar attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(attspace);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace for scalar attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(att);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }
}

/**
  * @brief Write an array attribute to the specified group or dataset
  *
  * The method creates the appropriate memory space and then creates the
  * attribute. Both the memory space and the attribute are also closed again.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param length Number of elements in the array attribute
  * @param name Name of the attribute
  * @param type HDF5 id for the type of the attribute (can be an alias from
  * HDF5types)
  * @param value C-style void pointer to the data that should be written to the
  * attribute
  */
inline void write_attribute_array(hid_t group, unsigned int length,
                                  std::string name, hid_t type, void* value) {
    hsize_t dims[1] = {length};

    // create dataspace
    hid_t attspace = H5Screate_simple(1, dims, NULL);
    if(attspace < 0) {
        std::cerr << "Error! Failed to create dataspace for array attribute \""
                  << name << "\" of length " << length << std::endl;
        my_exit();
    }

    // create attribute
    hid_t att = H5Acreate(group, name.c_str(), type, attspace, H5P_DEFAULT);
    if(att < 0) {
        std::cerr << "Error! Failed to create array attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // write attribute
    herr_t status = H5Awrite(att, type, value);
    if(status < 0) {
        std::cerr << "Error! Failed to write array attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(attspace);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace for array attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(att);
    if(status < 0) {
        std::cerr << "Error! Failed to close array attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }
}

/**
  * @brief Write a string attribute to the specified group or dataset
  *
  * Writing variable size strings requires the creation of the appropriate
  * HDF5 data type. The function then creates the appropriate memory space
  * and writes the attribute.
  * Both the attribute and the memory space are also closed again.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the attribute
  * @param value C-type void pointer to the data that should be written to the
  * attribute
  */
inline void write_attribute_string(hid_t group, std::string name,
                                   const void* value) {
    // create C-string datatype
    hid_t strtype = H5Tcopy(H5T_C_S1);
    if(strtype < 0) {
        std::cerr << "Error! Failed to copy C-string datatype for attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // set datatype length to variable
    herr_t status = H5Tset_size(strtype, H5T_VARIABLE);
    if(status < 0) {
        std::cerr << "Error! Failed to set size of C-string datatype for "
                     "attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // create dataspace
    hid_t attspace = H5Screate(H5S_SCALAR);
    if(attspace < 0) {
        std::cerr << "Error! Failed to create dataspace for string attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // create attribute
    hid_t att = H5Acreate(group, name.c_str(), strtype, attspace, H5P_DEFAULT);
    if(att < 0) {
        std::cerr << "Error! Failed to create string attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // write attribute
    status = H5Awrite(att, strtype, value);
    if(status < 0) {
        std::cerr << "Error! Failed to write string attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // close string type
    status = H5Tclose(strtype);
    if(status < 0) {
        std::cerr << "Error! Failed to close C-string datatype of attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(attspace);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace for scalar attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(att);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }
}

/**
  * @brief Write an array attribute to the specified group or dataset
  *
  * The method creates the appropriate memory space and then creates the
  * attribute. Both the memory space and the attribute are also closed again.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the attribute
  * @param type HDF5 id for the type of the attribute (can be an alias from
  * HDF5types)
  * @param data Reference to the std::vector containing the data to write
  */
template <typename T>
inline void write_attribute_array(hid_t group, std::string name, hid_t type,
                                  std::vector<T>& data) {
    hsize_t dims[1] = {data.size()};

    // create dataspace
    hid_t attspace = H5Screate_simple(1, dims, NULL);
    if(attspace < 0) {
        std::cerr << "Error! Failed to create dataspace for array attribute \""
                  << name << "\" of length " << data.size() << std::endl;
        my_exit();
    }

    // create attribute
    hid_t att = H5Acreate(group, name.c_str(), type, attspace,
                          HDF5constants::DEFAULT_PROPERTY_LIST);
    if(att < 0) {
        std::cerr << "Error! Failed to create array attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // write attribute
    herr_t status = H5Awrite(att, type, &data[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to write array attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(attspace);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace for array attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(att);
    if(status < 0) {
        std::cerr << "Error! Failed to close array attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }
}

/**
  * @brief Write an array attribute to the specified group or dataset
  *
  * The method creates the appropriate memory space and then creates the
  * attribute. Both the memory space and the attribute are also closed again.
  *
  * This version writes a single element and hence creates an array with length
  * 1.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the attribute
  * @param type HDF5 id for the type of the attribute (can be an alias from
  * HDF5types)
  * @param data Reference to the single data value to write
  */
template <typename T>
inline void write_attribute_array(hid_t group, std::string name, hid_t type,
                                  T& data) {
    hsize_t dims[1] = {1};

    // create dataspace
    hid_t attspace = H5Screate_simple(1, dims, NULL);
    if(attspace < 0) {
        std::cerr << "Error! Failed to create dataspace for array attribute \""
                  << name << "\" of length " << 1 << std::endl;
        my_exit();
    }

    // create attribute
    hid_t att = H5Acreate(group, name.c_str(), type, attspace,
                          HDF5constants::DEFAULT_PROPERTY_LIST);
    if(att < 0) {
        std::cerr << "Error! Failed to create array attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // write attribute
    herr_t status = H5Awrite(att, type, &data);
    if(status < 0) {
        std::cerr << "Error! Failed to write array attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(attspace);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace for array attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(att);
    if(status < 0) {
        std::cerr << "Error! Failed to close array attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }
}

/**
  * @brief Write a scalar dataset to the specified group
  *
  * The method creates the dataset, writes the data to it and then closes it
  * again.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the dataset
  * @param type HDF5 id for the type of the dataset (can be an alias from
  * HDF5types)
  * @param data Reference to the std::vector containing the data to write to the
  * new dataset
  */
template <typename T>
inline void write_dataset_scalar(hid_t group, std::string name, hid_t type,
                                 std::vector<T>& data) {
    hsize_t dims[1] = {data.size()};

    // create dataspace
    hid_t filespace = H5Screate_simple(1, dims, NULL);
    if(filespace < 0) {
        std::cerr << "Error! Failed to create dataspace for scalar dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // create dataset
    hid_t dataset = H5Dcreate(group, name.c_str(), type, filespace,
                              HDF5constants::DEFAULT_PROPERTY_LIST);
    if(dataset < 0) {
        std::cerr << "Error! Failed to create scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // write dataset
    herr_t status =
            H5Dwrite(dataset, type, HDF5constants::COMPLETE_DATASPACE,
                     filespace, HDF5constants::DEFAULT_PROPERTY_LIST, &data[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to write scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataset
    status = H5Dclose(dataset);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar dataset \"" << name << "\""
                  << std::endl;
    }
}

/**
  * @brief Write a vector dataset to the specified group
  *
  * The method creates the dataset, writes the data to it and then closes it
  * again. The data is assumed to have a length that is a multiple of 3, wherein
  * 3 consecutive elements are considered to be the x-, y- and z-components of a
  * 3-element vector. They will be written to separate columns in the dataset.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the dataset
  * @param type HDF5 id for the type of the attribute (can be an alias from
  * HDF5types)
  * @param data Reference to the std::vector containing the data to write to the
  * new dataset
  */
template <typename T>
inline void write_dataset_vector(hid_t group, std::string name, hid_t type,
                                 std::vector<T>& data) {
    hsize_t dims[2] = {data.size() / 3, 3};

    // create dataspace
    hid_t filespace = H5Screate_simple(2, dims, NULL);
    if(filespace < 0) {
        std::cerr << "Error! Failed to create dataspace for vector dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // create dataset
    hid_t dataset = H5Dcreate(group, name.c_str(), type, filespace,
                              HDF5constants::DEFAULT_PROPERTY_LIST);
    if(dataset < 0) {
        std::cerr << "Error! Failed to create vector dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // write dataset
    herr_t status =
            H5Dwrite(dataset, type, HDF5constants::DEFAULT_PROPERTY_LIST,
                     filespace, HDF5constants::DEFAULT_PROPERTY_LIST, &data[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to write vector dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataset
    status = H5Dclose(dataset);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar dataset \"" << name << "\""
                  << std::endl;
    }
}

/**
 * @brief Create an extendable scalar dataset with the given name, type and size
 *
 * The dataset can be filled later on using
 * HDF5tools::write_dataset_scalar_chunk()
 *
 * @param group HDF5 id of a group
 * @param name Name of the new dataset to create
 * @param type HDF5 id for the type of the dataset
 * @param size Size of the data that will be written to the dataset in total
 * @return HDF5 id for the newly created dataset
 */
inline hid_t create_dataset_scalar(hid_t group, std::string name, hid_t type,
                                   unsigned int size) {
    hsize_t dims[1] = {size};

    // create dataspace
    hid_t filespace = H5Screate_simple(1, dims, NULL);
    if(filespace < 0) {
        std::cerr << "Error! Failed to create dataspace for scalar dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // create dataset
    hid_t dataset = H5Dcreate(group, name.c_str(), type, filespace,
                              HDF5constants::DEFAULT_PROPERTY_LIST);
    if(dataset < 0) {
        std::cerr << "Error! Failed to create scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    herr_t status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    return dataset;
}

/**
 * @brief Create an extendable vector dataset with the given name, type and size
 *
 * The dataset can be filled later on using
 * HDF5tools::write_dataset_vector_chunk()
 *
 * @param group HDF5 id of a group
 * @param name Name of the new dataset to create
 * @param type HDF5 id for the type of the dataset
 * @param size Size of the data that will be written to the dataset in total
 * @return HDF5 id for the newly created dataset
 */
inline hid_t create_dataset_vector(hid_t group, std::string name, hid_t type,
                                   unsigned int size) {
    hsize_t dims[2] = {size, 3};

    // create dataspace
    hid_t filespace = H5Screate_simple(2, dims, NULL);
    if(filespace < 0) {
        std::cerr << "Error! Failed to create dataspace for scalar dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // create dataset
    hid_t dataset = H5Dcreate(group, name.c_str(), type, filespace,
                              HDF5constants::DEFAULT_PROPERTY_LIST);
    if(dataset < 0) {
        std::cerr << "Error! Failed to create scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    herr_t status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    return dataset;
}

/**
 * @brief Write a chunk of data to an existing scalar dataset
 *
 * @param dataset HDF5 id of a dataset
 * @param type HDF5 id of the dataset type
 * @param offset Offset of the data in the dataset
 * @param data Data to write
 */
template <typename T>
inline void write_dataset_scalar_chunk(hid_t dataset, hid_t type,
                                       unsigned int offset,
                                       std::vector<T>& data) {
    // get filespace
    hid_t filespace = H5Dget_space(dataset);
    if(filespace < 0) {
        std::cerr << "Error! Failed to obtain dataspace of scalar dataset!"
                  << std::endl;
        my_exit();
    }

    // select hyperslab in filespace
    hsize_t dims[1] = {data.size()};
    hsize_t offs[1] = {offset};
    herr_t status =
            H5Sselect_hyperslab(filespace, HDF5constants::SET_NEW_SELECTION,
                                offs, HDF5constants::CONTINUOUS_SELECTION, dims,
                                HDF5constants::SINGLE_ELEMENT_BLOCKS);
    if(status < 0) {
        std::cerr << "Error! Unable to select hyperslab of dataset filespace!"
                  << std::endl;
        my_exit();
    }

    // create memory space
    hid_t memspace = H5Screate_simple(1, dims, NULL);
    if(memspace < 0) {
        std::cerr << "Error! Failed to create memory space to write chunk to "
                     "scalar dataset!"
                  << std::endl;
        my_exit();
    }

    // write to dataset
    status = H5Dwrite(dataset, type, memspace, filespace,
                      HDF5constants::DEFAULT_PROPERTY_LIST, &data[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to write chunk to scalar dataset!"
                  << std::endl;
        my_exit();
    }

    // close memory space
    status = H5Sclose(memspace);
    if(status < 0) {
        std::cerr << "Error! Failed to close memory space for scalar dataset!"
                  << std::endl;
        my_exit();
    }

    // close filespace
    status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close filespace of scalar dataset!"
                  << std::endl;
        my_exit();
    }
}

/**
 * @brief Write a chunk of data to an existing vector dataset
 *
 * @param dataset HDF5 id of a dataset
 * @param type HDF5 id of the dataset type
 * @param offset Offset of the data in the dataset
 * @param data Data to write
 */
template <typename T>
inline void write_dataset_vector_chunk(hid_t dataset, hid_t type,
                                       unsigned int offset,
                                       std::vector<T>& data) {
    // get filespace
    hid_t filespace = H5Dget_space(dataset);
    if(filespace < 0) {
        std::cerr << "Error! Failed to obtain dataspace of vector dataset!"
                  << std::endl;
        my_exit();
    }

    // select hyperslab in filespace
    hsize_t dims[2] = {data.size() / 3, 3};
    hsize_t offs[2] = {offset, 0};
    herr_t status =
            H5Sselect_hyperslab(filespace, HDF5constants::SET_NEW_SELECTION,
                                offs, HDF5constants::CONTINUOUS_SELECTION, dims,
                                HDF5constants::SINGLE_ELEMENT_BLOCKS);
    if(status < 0) {
        std::cerr << "Error! Unable to select hyperslab of dataset filespace!"
                  << std::endl;
        my_exit();
    }

    // create memory space
    hid_t memspace = H5Screate_simple(2, dims, NULL);
    if(memspace < 0) {
        std::cerr << "Error! Failed to create memory space to write chunk to "
                     "vector dataset!"
                  << std::endl;
        my_exit();
    }

    // write to dataset
    status = H5Dwrite(dataset, type, memspace, filespace,
                      HDF5constants::DEFAULT_PROPERTY_LIST, &data[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to write chunk to vector dataset!"
                  << std::endl;
        my_exit();
    }

    // close memory space
    status = H5Sclose(memspace);
    if(status < 0) {
        std::cerr << "Error! Failed to close memory space for vector dataset!"
                  << std::endl;
        my_exit();
    }

    // close filespace
    status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close filespace of vector dataset!"
                  << std::endl;
        my_exit();
    }
}

/// read functions

/**
  * @brief Read an array or scalar attribute from the specified group or dataset
  *
  * The attribute is opened, read and closed again.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the attribute. Errors are thrown if there is no
  * attribute with this name
  * @param type HDF5 id for the type of the attribute (can be an alias from
  * HDF5types)
  * @param value C-type void pointer to the memory location where the data
  * should be written to
  */
inline void read_attribute(hid_t group, std::string name, hid_t type,
                           void* value) {
    // open attribute
    hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
    if(attr < 0) {
        std::cerr << "Error! Failed to open attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // read attribute
    herr_t status = H5Aread(attr, type, value);
    if(status < 0) {
        std::cerr << "Error! Failed to read attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(attr);
    if(status < 0) {
        std::cerr << "Error! Failed to close attribute \"" << name << "\""
                  << std::endl;
    }
}

/**
  * @brief Read a string attribute from the specified group or dataset
  *
  * Handling variable size strings with HDF5 requires the creation
  * of a specific string data type. The attribute is opened and the
  * data is written to a char buffer. This buffer is then converted
  * to a string. The attribute is closed together with the memory
  * space that was needed to read it.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the attribute. Errors are thrown if there is no
  * attribute with this name
  * @return The string contained in the attribute
  */
inline std::string read_attribute_string(hid_t group, std::string name) {
    // open attribute
    hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
    if(attr < 0) {
        std::cerr << "Error! Failed to open string attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    char* data[1];

    // create C-string datatype
    hid_t strtype = H5Tcopy(H5T_C_S1);
    if(strtype < 0) {
        std::cerr << "Error! Failed to copy C-string datatype for attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // set datatype length to variable
    herr_t status = H5Tset_size(strtype, H5T_VARIABLE);
    if(status < 0) {
        std::cerr << "Error! Failed to set size of C-string datatype for "
                     "attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // read attribute
    status = H5Aread(attr, strtype, data);
    if(status < 0) {
        std::cerr << "Error! Failed to read string attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close string type
    status = H5Tclose(strtype);
    if(status < 0) {
        std::cerr << "Error! Failed to close C-string datatype for attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(attr);
    if(status < 0) {
        std::cerr << "Error! Failed to close string attribute \"" << name
                  << "\"" << std::endl;
    }

    // make string
    std::string attstring(data[0]);

    // free buffer
    free(data[0]);

    return attstring;
}

/**
  * @brief Read a scalar attribute from the specified group or dataset
  *
  * The attribute is opened, read and closed again and its value is returned.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the attribute. Errors are thrown if there is no
  * attribute with this name
  * @param type HDF5 id for the type of the attribute (can be an alias from
  * HDF5types)
  * @return The value of the requested attribute
  */
template <typename T>
inline T read_attribute_scalar(hid_t group, std::string name, hid_t type) {
    // open attribute
    hid_t attr =
            H5Aopen(group, name.c_str(), HDF5constants::DEFAULT_PROPERTY_LIST);
    if(attr < 0) {
        std::cerr << "Error! Failed to open scalar attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // read attribute
    T value;
    herr_t status = H5Aread(attr, type, &value);
    if(status < 0) {
        std::cerr << "Error! Failed to read scalar attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(attr);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    return value;
}

/**
  * @brief Read a vector attribute from the specified group or dataset
  *
  * The attribute is opened, read and closed again and its value is returned.
  * If a single attribute is found, it is returned as a std::vector with a
  * single element.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the attribute. Errors are thrown if there is no
  * attribute with this name
  * @param type HDF5 id for the type of the attribute (can be an alias from
  * HDF5types)
  * @return std::vector containing the data in the attribute
  */
template <typename T>
inline std::vector<T> read_attribute_vector(hid_t group, std::string name,
                                            hid_t type) {
    // open attribute
    hid_t attr =
            H5Aopen(group, name.c_str(), HDF5constants::DEFAULT_PROPERTY_LIST);
    if(attr < 0) {
        std::cerr << "Error! Failed to open vector attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // open attribute dataspace
    hid_t space = H5Aget_space(attr);
    if(space < 0) {
        std::cerr << "Error! Failed to open dataspace of vector attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // query dataspace extents
    hsize_t size[1];
    hsize_t maxsize[1];
    int ndim = H5Sget_simple_extent_dims(space, size, maxsize);
    if(!ndim) {
        size[0] = 1;
    }
    if(ndim < 0) {
        std::cerr
                << "Error! Unable to query extent of dataspace of attribute \""
                << name << "\"" << std::endl;
        my_exit();
    }

    // read attribute
    std::vector<T> value(size[0]);
    herr_t status = H5Aread(attr, type, &value[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to read vector attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(space);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace of vector attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(attr);
    if(status < 0) {
        std::cerr << "Error! Failed to close vector attribute \"" << name
                  << "\"" << std::endl;
    }

    return value;
}

/**
  * @brief Read a scalar dataset from the specified group
  *
  * The dataset is opened, its dimensions are read and used to create a
  * std::vector in which the dataset contents are stored.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the dataset.
  * @param type HDF5 id for the type of the dataset (can be an alias from
  * HDF5types)
  * @return std::vector containing the data in the dataset
  */
template <typename T>
inline std::vector<T> read_dataset_scalar(hid_t group, std::string name,
                                          hid_t type) {
    // open dataset
    hid_t dataset = H5Dopen(group, name.c_str());
    if(dataset < 0) {
        std::cerr << "Error! Failed to open scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // open dataset dataspace
    hid_t filespace = H5Dget_space(dataset);
    if(filespace < 0) {
        std::cerr << "Error! Failed to open dataspace of scalar dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // query dataspace extents
    hsize_t size[1];
    hsize_t maxsize[1];
    int ndim = H5Sget_simple_extent_dims(filespace, size, maxsize);
    if(ndim < 0) {
        std::cerr << "Error! Unable to query extent of dataspace of dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // read dataset
    std::vector<T> data(size[0]);
    herr_t status = H5Dread(dataset, type, HDF5constants::COMPLETE_DATASPACE,
                            HDF5constants::COMPLETE_DATASPACE,
                            HDF5constants::DEFAULT_PROPERTY_LIST, &data[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to read scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace of scalar dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close dataset
    status = H5Dclose(dataset);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    return data;
}

/**
  * @brief Read a scalar dataset chunk from the specified group
  *
  * The dataset is opened, its dimensions are read and used to create a
  * std::vector in which the dataset contents are stored.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param offset Offset of the data to be read
  * @param nlines Number of lines to be read
  * @param name Name of the dataset.
  * @param type HDF5 id for the type of the dataset (can be an alias from
  * HDF5types)
  * @return std::vector containing the data in the dataset
  */
template <typename T>
inline std::vector<T> read_dataset_scalar_chunk(hid_t group,
                                                unsigned int offset,
                                                unsigned int nlines,
                                                std::string name, hid_t type) {
    // open dataset
    hid_t dataset = H5Dopen(group, name.c_str());
    if(dataset < 0) {
        std::cerr << "Error! Failed to open scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // open dataset dataspace
    hid_t filespace = H5Dget_space(dataset);
    if(filespace < 0) {
        std::cerr << "Error! Failed to open dataspace of scalar dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // select hyperslab in filespace
    hsize_t dims[2] = {nlines, 1};
    hsize_t offs[2] = {offset, 0};
    herr_t status =
            H5Sselect_hyperslab(filespace, HDF5constants::SET_NEW_SELECTION,
                                offs, HDF5constants::CONTINUOUS_SELECTION, dims,
                                HDF5constants::SINGLE_ELEMENT_BLOCKS);
    if(status < 0) {
        std::cerr << "Error! Unable to select hyperslab of dataset filespace!"
                  << std::endl;
        my_exit();
    }

    // create memory space
    hid_t memspace = H5Screate_simple(2, dims, NULL);
    if(memspace < 0) {
        std::cerr << "Error! Failed to create memory space for scalar dataset "
                     "chunk"
                  << std::endl;
        my_exit();
    }

    // read dataset
    std::vector<T> data(nlines);
    status = H5Dread(dataset, type, memspace, filespace,
                     HDF5constants::DEFAULT_PROPERTY_LIST, &data[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to read scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close memory space
    status = H5Sclose(memspace);
    if(status < 0) {
        std::cerr << "Error! Failed to close memory space for scalar dataset "
                     "chunk"
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace of scalar dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close dataset
    status = H5Dclose(dataset);
    if(status < 0) {
        std::cerr << "Error! Failed to close scalar dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    return data;
}

/**
  * @brief Read a vector dataset from the specified group
  *
  * The dataset is opened, its dimensions are read and used to create a
  * std::vector in which the dataset contents are stored. No assumptions
  * have to be made about the size of the vectors contained in the dataset.
  * Vector elements are read row-by-row and then pasted together, such that
  * the first n elements will be the n components of the first vector and
  * so on.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param name Name of the dataset.
  * @param type HDF5 id for the type of the dataset (can be an alias from
  * HDF5types)
  * @return std::vector containing the data in the dataset
  */
template <typename T>
inline std::vector<T> read_dataset_vector(hid_t group, std::string name,
                                          hid_t type) {
    // open dataset
    hid_t dataset = H5Dopen(group, name.c_str());
    if(dataset < 0) {
        std::cerr << "Error! Failed to open vector dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // open dataspace
    hid_t filespace = H5Dget_space(dataset);
    if(filespace < 0) {
        std::cerr << "Error! Failed to open dataspace of vector dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // query dataspace extents
    hsize_t size[2];
    hsize_t maxsize[2];
    int ndim = H5Sget_simple_extent_dims(filespace, size, maxsize);
    if(ndim < 0) {
        std::cerr << "Error! Unable to query extent of dataspace of dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // read dataset
    std::vector<T> data(size[0] * size[1]);
    herr_t status = H5Dread(dataset, type, HDF5constants::COMPLETE_DATASPACE,
                            HDF5constants::COMPLETE_DATASPACE,
                            HDF5constants::DEFAULT_PROPERTY_LIST, &data[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to read vector dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace of vector dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close dataset
    status = H5Dclose(dataset);
    if(status < 0) {
        std::cerr << "Error! Failed to close vector dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    return data;
}

/**
  * @brief Read a vector dataset chunk from the specified group
  *
  * The dataset is opened, its dimensions are read and used to create a
  * std::vector in which the dataset contents are stored. No assumptions
  * have to be made about the size of the vectors contained in the dataset.
  * Vector elements are read row-by-row and then pasted together, such that
  * the first n elements will be the n components of the first vector and
  * so on.
  *
  * @param group HDF5 id of either a dataset or a group
  * @param offset Offset of the data to be read
  * @param nlines Number of lines to be read
  * @param name Name of the dataset.
  * @param type HDF5 id for the type of the dataset (can be an alias from
  * HDF5types)
  * @return std::vector containing the data in the dataset
  */
template <typename T>
inline std::vector<T> read_dataset_vector_chunk(hid_t group,
                                                unsigned int offset,
                                                unsigned int nlines,
                                                std::string name, hid_t type) {
    // open dataset
    hid_t dataset = H5Dopen(group, name.c_str());
    if(dataset < 0) {
        std::cerr << "Error! Failed to open vector dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // open dataspace
    hid_t filespace = H5Dget_space(dataset);
    if(filespace < 0) {
        std::cerr << "Error! Failed to open dataspace of vector dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // select hyperslab in filespace
    hsize_t dims[2] = {nlines, 3};
    hsize_t offs[2] = {offset, 0};
    herr_t status =
            H5Sselect_hyperslab(filespace, HDF5constants::SET_NEW_SELECTION,
                                offs, HDF5constants::CONTINUOUS_SELECTION, dims,
                                HDF5constants::SINGLE_ELEMENT_BLOCKS);
    if(status < 0) {
        std::cerr << "Error! Unable to select hyperslab of dataset filespace!"
                  << std::endl;
        my_exit();
    }

    // create memory space
    hid_t memspace = H5Screate_simple(2, dims, NULL);
    if(memspace < 0) {
        std::cerr << "Error! Failed to create memory space for scalar dataset "
                     "chunk"
                  << std::endl;
        my_exit();
    }

    // read dataset
    std::vector<T> data(nlines * 3);
    status = H5Dread(dataset, type, memspace, filespace,
                     HDF5constants::DEFAULT_PROPERTY_LIST, &data[0]);
    if(status < 0) {
        std::cerr << "Error! Failed to read vector dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // close memory space
    status = H5Sclose(memspace);
    if(status < 0) {
        std::cerr << "Error! Failed to close memory space for scalar dataset "
                     "chunk"
                  << std::endl;
        my_exit();
    }

    // close dataspace
    status = H5Sclose(filespace);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace of vector dataset \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close dataset
    status = H5Dclose(dataset);
    if(status < 0) {
        std::cerr << "Error! Failed to close vector dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    return data;
}

/// Other utility functions

/**
 * @brief Determine if the attribute from the given group with the given name
 * has the given type
 *
 * @param group HDF5 id of either a dataset or a group
 * @param name Name of the attribute. Errors are thrown if there is no attribute
 * with this name
 * @param type HDF5 id for a type we want to compare to the actual attribute
 * type
 * @return True if the given type is compatible with the type of the attribute,
 * false otherwise
 */
inline bool attribute_has_type(hid_t group, std::string name, hid_t type) {
    // open attribute
    hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
    if(attr < 0) {
        std::cerr << "Error! Failed to open attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // get attribute type
    hid_t atype = H5Aget_type(attr);
    if(atype < 0) {
        std::cerr << "Error! Unable to determine type of attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // close attribute
    herr_t status = H5Aclose(attr);
    if(status < 0) {
        std::cerr << "Error! Failed to close attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // compare attribute type with given type
    htri_t equal = H5Tequal(atype, type);
    if(equal < 0) {
        std::cerr << "Error! Unable to compare type of attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // close attribute type
    status = H5Tclose(atype);
    if(status < 0) {
        std::cerr << "Error! Failed to close type of attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    return equal > 0;
}

/**
 * @brief Determine if the dataset from the given group with the given name has
 * the given type
 *
 * @param group HDF5 id of a group
 * @param name Name of the dataset
 * @param type HDF5 id for a type we want to compare to the actual attribute
 * type
 * @return True if the given type is compatible with the type of the dataset,
 * false otherwise
 */
inline bool dataset_has_type(hid_t group, std::string name, hid_t type) {
    // open dataset
    hid_t dataset = H5Dopen(group, name.c_str());
    if(dataset < 0) {
        std::cerr << "Error! Failed to open dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // get dataset type
    hid_t dtype = H5Dget_type(dataset);
    if(dtype < 0) {
        std::cerr << "Error! Unable to determine type of dataset \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // compare dataset type with given type
    htri_t equal = H5Tequal(dtype, type);
    if(equal < 0) {
        std::cerr << "Error! Unable to compare type of dataset \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // close datatype
    herr_t status = H5Tclose(dtype);
    if(status < 0) {
        std::cerr << "Error! Failed to close datatype of dataset \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // close dataset
    status = H5Dclose(dataset);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataset \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    return equal > 0;
}

/**
  * @brief Determine the size of the attribute of the given group with the given
  * name
  *
  * @param group HDF5 id of a group or dataset
  * @param name Name of the attribute
  * @return The length of the attribute. If the attribute is a scalar, length 1
  * is returned.
  */
inline unsigned int get_attribute_size(hid_t group, std::string name) {
    // open attribute
    hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
    if(attr < 0) {
        std::cerr << "Error! Failed to open attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    // get attribute dataspace
    hid_t space = H5Aget_space(attr);
    if(space < 0) {
        std::cerr << "Error! Failed to open dataspace of attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // get dataspace extents
    hsize_t size[1];
    hsize_t maxsize[1];
    int ndim = H5Sget_simple_extent_dims(space, size, maxsize);
    if(ndim < 0) {
        std::cerr << "Error! Unable to query extent of dataspace of "
                     "attribute \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    // close dataspace
    herr_t status = H5Sclose(space);
    if(status < 0) {
        std::cerr << "Error! Failed to close dataspace of attribute \"" << name
                  << "\"" << std::endl;
        my_exit();
    }

    // close attribute
    status = H5Aclose(attr);
    if(status < 0) {
        std::cerr << "Error! Failed to close attribute \"" << name << "\""
                  << std::endl;
        my_exit();
    }

    if(!ndim) {
        size[0] = 1;
    }

    return size[0];
}

/**
  * @brief Check if the given group has a child with the given name
  *
  * @param group HDF5 id of a group
  * @param name Name of a dataset or attribute
  * @return True if the given name is a child of group, false otherwise
  */
inline bool exists(hid_t group, std::string name) {
    // query group
    htri_t link = H5Lexists(group, name.c_str(), H5P_DEFAULT);
    if(link < 0) {
        std::cerr << "Error! Unable to query existence of object with name \""
                  << name << "\"" << std::endl;
        my_exit();
    }

    return link > 0;
}

/**
 * @brief Turn off default HDF5 error handling.
 *
 * This allows us to use our own more useful error messages without the output
 * being cluttered with meaningless HDF5 errors.
 */
inline void turn_off_error_handling() {
    herr_t status = H5Eset_auto(NULL, NULL);
    if(status < 0) {
        std::cerr << "Error! Unable to turn off default HDF5 error handling!"
                  << std::endl;
        my_exit();
    }
}

}  // namespace HDF5tools

#endif  // HDF5TOOLS_HPP
