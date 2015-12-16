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
 * @file HelperFunctions.cpp
 *
 * @brief Various auxiliary functions: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "HelperFunctions.hpp"
#include <algorithm>
#include <sstream>
#include <vector>
using namespace std;

/**
 * @brief Convert the given byte size to a human readable string using KB, MB...
 *
 * Since the byte size is an integer, more than GB is impossible.
 *
 * @param bytes Integer size in bytes
 * @return String representation of the byte size
 */
string HelperFunctions::human_readable_bytes(unsigned int bytes){
    stringstream output;
    unsigned int i = 0;
    double number = (double)bytes;
    while((bytes >> 10)){
        bytes >>= 10;
        i++;
    }
    // funny detail: when you only allow unsigned int sizes, you cannot have
    // more than GB
    string names[4] = {"bytes", "KB", "MB", "GB"};
    unsigned int sizes[4] = {1, 1<<10, 1<<20, 1<<30};
    number /= ((double)sizes[i]);
    output << number << " " << names[i];
    return output.str();
}

/**
 * @brief Convert the given string to an integer byte size
 *
 * @param input String containing a number and a byte unit (KB, MB...)
 * @return Size in bytes
 */
unsigned int HelperFunctions::machine_readable_bytes(std::string input){
    unsigned int bytes = 1;
    string names[4] = {"bytes", "KB", "MB", "GB"};
    unsigned int sizes[4] = {1, 1<<10, 1<<20, 1<<30};
    for(unsigned int i = 0; i < 4; i++){
        if(input.find(names[i]) != string::npos){
            bytes = sizes[i];
        }
    }

    istringstream inputstream(input);
    double number;
    inputstream >> number;
    return (number*bytes);
}

/**
 * @brief Convert the given long integer to a human readable string with comma
 * seperated mutiples of 1000
 *
 * @param counter Unsigned long integer counter value
 * @return Human readable counter string
 */
string HelperFunctions::human_readable_counter(unsigned long counter){
    stringstream output;
    if(counter){
        vector <unsigned int> counterparts;
        while(counter){
            unsigned int part = counter%1000;
            counterparts.push_back(part);
            counter /= 1000;
        }
        reverse(counterparts.begin(), counterparts.end());
        for(unsigned int i = 0; i < counterparts.size()-1; i++){
            if(i){
                output.fill('0');
                output.width(3);
            }
            output << counterparts[i] << ",";
        }
        if(counterparts.size() > 1){
            output.fill('0');
            output.width(3);
        }
        output << counterparts.back();
        counterparts.clear();
    } else {
        output << "0";
    }
    return output.str();
}

/**
 * @brief Convert the given long integer to a power of two to make it more human
 * readable
 *
 * @param num Unsigned long integer value
 * @return Human readable power of two string
 */
string HelperFunctions::human_readable_long(unsigned long num){
    stringstream output;
    if(num){
        unsigned int exponent = 60;
        while(!(num>>exponent)){
            exponent--;
        }
        output << "2^" << exponent;
        unsigned long curnum = 1;
        curnum <<= exponent;
        num -= curnum;
        while(num){
            exponent = 60;
            while(!(num>>exponent)){
                exponent--;
            }
            output << "+2^" << exponent;
            curnum = 1;
            curnum <<= exponent;
            num -= curnum;
        }
    } else {
        output << num;
    }
    return output.str();
}

/**
 * @brief Make sure the given filename ends with .hdf5
 *
 * If it already has the correct extension, we do nothing. If not, we add the
 * extension to the name.
 *
 * @param name Name of the file
 * @return Name of the file with a guaranteed .hdf5 extension
 */
std::string HelperFunctions::make_hdf5_file(std::string name){
    if(name.find(".hdf5") != name.size()-5){
        return name + string(".hdf5");
    } else {
        return name;
    }
}
