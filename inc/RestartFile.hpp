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
 * @file RestartFile.hpp
 *
 * @brief Restart file support: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef RESTARTFILE_HPP
#define RESTARTFILE_HPP

#include <fstream>
#include <string>
#include <vector>

#include "Vec.hpp"
#include "StateVector.hpp"

/**
 * @brief Wrapper around the std::ofstream or std::ifstream that is used to dump
 * a Simulation to a restart-file or the read it in again
 *
 * Upon construction the RestartFile writes useful system and simulation info to
 * the ofstream or reads it from the ifstream.
 * The class also provides template functions to write and read the most common
 * data types to and from the restart-file. Other classes either write/read
 * their members using these functions or invoke the dump method or restart
 * constructor of their members. No explicit reading or writing should be done
 * outside this class.
 */
class RestartFile{
private:
    /*! \brief Output file stream if this RestartFile is used to dump data */
    std::ofstream _ofile;

    /*! \brief Input file stream if this RestartFile is used to restart */
    std::ifstream _ifile;

    std::string get_padded_digit(int digit, int limit);

    std::string get_filename(std::string prefix, int rank, int size,
                             std::string suffix = "");
    std::string get_foldername(std::string prefix, int rank, int size,
                               std::string suffix = "");

public:
    RestartFile(std::string outputdir, double simtime);
    RestartFile(std::string filename);

    /**
     * @brief Write a standard datatype value to the stream
     *
     * @param value Value to write to the stream. Can be any type that can be
     * casted to a binary array with a size known at compile time
     */
    template<typename T> void write(T value){
        _ofile.write(reinterpret_cast<char*>(&value), sizeof(T));
    }

    /**
     * @brief RestartFile::write() specialization for booleans
     *
     * The boolean is first converted to a byte and then written to the stream.
     *
     * @param value Boolean value to write to the stream
     */
    void write(bool value){
        char boolean = value;
        _ofile.write(&boolean, 1);
    }

    /**
     * @brief RestartFile::write() specialization for Vec
     *
     * We know the internal structure of the Vec and use this knowledge to
     * directly write its internal double array to the stream.
     *
     * @param value Vec to write to the stream
     */
    void write(Vec &value){
        _ofile.write(reinterpret_cast<char*>(&value[0]), ndim_*sizeof(double));
    }

    /**
     * @brief RestartFile::write() specialization for StateVector
     *
     * We know the internal structure of the StateVector and use this knowledge
     * to directly write its internal double array to the stream.
     *
     * @param value StateVector to write to the stream
     */
    void write(StateVector &value){
        _ofile.write(reinterpret_cast<char*>(&value[0]),
                     (ndim_+2)*sizeof(double));
    }

    /**
     * @brief RestartFile::write() specialization for std::vector
     *
     * By using a recursive implementation, things like
     * write(vector< vector<int> >) actually work!
     *
     * @param vec std::vector containing any other datatype that can be written
     * to the stream
     */
    template<typename T> void write(std::vector<T> &vec){
        unsigned int vsize = vec.size();
        write(vsize);
        for(unsigned int i = 0; i < vsize; i++){
            write(vec[i]);
        }
    }

    /**
     * @brief RestartFile::write() specialization for std::string
     *
     * A string is under the hood a binary array with a length. We first write
     * this length to the stream and then the raw array.
     *
     * @param value std::string to write to the stream
     */
    void write(std::string &value){
        unsigned int ssize = value.size();
        write(ssize);
        _ofile.write(value.c_str(), ssize);
    }

    /**
     * @brief Write the given array with the given size to the stream
     *
     * The datatype of the array can be any standard type that supports casting
     * to a byte array and with a fixed size at compile time.
     *
     * @param array Array to write to the stream
     * @param size unsigned integer length of the array
     */
    template<typename T> void write(T *array, unsigned int size){
        _ofile.write(reinterpret_cast<char*>(array), sizeof(T)*size);
    }

    /**
     * @brief Write the given boolean array with the given size to the stream
     *
     * The template method does not work for boolean arrays, since a boolean
     * cannot be simply casted to a char.
     *
     * @param array Array to write to the stream
     * @param size unsigned integer length of the array
     */
    void write(bool *array, unsigned int size){
        char *byte_array = new char[size];
        for(unsigned int i = 0; i < size; i++){
            byte_array[i] = array[i];
        }
        _ofile.write(byte_array, size);
        delete [] byte_array;
    }

    /**
     * @brief Read a value with a standard datatype from the stream
     *
     * The datatype can be any standard type that supports casting to a byte
     * array and with a fixed size at compile time.
     *
     * @param value Value to read from the stream
     */
    template<typename T> void read(T &value){
        _ifile.read(reinterpret_cast<char*>(&value), sizeof(value));
    }

    /**
     * @brief RestartFile::read() specialization for Vec
     *
     * We use our knowledge of the internal structure of Vec to directly read
     * in its internal double array.
     *
     * @param value Vec to read from the stream
     */
    void read(Vec &value){
        _ifile.read(reinterpret_cast<char*>(&value[0]), ndim_*sizeof(double));
    }

    /**
     * @brief RestartFile::read() specialization for StateVector
     *
     * We use our knowledge of the internal structure of StateVector to directly
     * read in its internal double array.
     *
     * @param value StateVector to read from the stream
     */
    void read(StateVector &value){
        _ifile.read(reinterpret_cast<char*>(&value[0]),
                    (ndim_+2)*sizeof(double));
    }

    /**
     * @brief RestartFile::read() specialization for std::vector
     *
     * By using a recursive definition, things like
     * read(std::vector< std::vector<int> >) actually work!
     *
     * @param vec std::vector to read from the stream
     */
    template<typename T> void read(std::vector<T> &vec){
        unsigned int vsize;
        read(vsize);
        vec.resize(vsize);
        for(unsigned int i = 0; i < vsize; i++){
            T element;
            read(element);
            vec[i] = element;
        }
    }

    /**
     * @brief RestartFile::read() specialization for std::string
     *
     * This one is quite complicated. We first read the length of the string
     * from the stream and then initialize a temporary byte array with this
     * length +1. We store the byte array from the stream in this temporary
     * buffer and set the last element of the buffer to the string terminating
     * character. We then initialize a std::string based on the temporary
     * buffer, free the temporary buffer and return the std::string.
     *
     * @param value std::string to read from the stream
     */
    void read(std::string &value){
        unsigned int ssize;
        read(ssize);
        char* stringval = new char[ssize+1];
        _ifile.read(stringval, ssize);
        stringval[ssize] = '\0';
        value = string(stringval);
        delete [] stringval;
    }

    /**
     * @brief Read an array with the given size from the stream
     *
     * The datatype of the array can be any standard datatype that supports
     * casting to a byte array and has a fixed size at compile time.
     *
     * @param array Array to read from the stream
     * @param size Size of the array
     */
    template<typename T> void read(T *array, unsigned int size){
        _ifile.read(reinterpret_cast<char*>(array), sizeof(T)*size);
    }

    /**
     * @brief Read a boolean array with the given size from the stream
     *
     * We have to provide a custom implementation for boolean arrays, since a
     * boolean cannot be casted to a char.
     *
     * @param array Array to read from the stream
     * @param size Size of the array
     */
    void read(bool *array, unsigned int size){
        char *byte_array = new char[size];
        _ifile.read(byte_array, size);
        for(unsigned int i = 0; i < size; i++){
            array[i] = byte_array[i] > 0;
        }
        delete [] byte_array;
    }
};

#endif // RESTARTFILE_HPP
