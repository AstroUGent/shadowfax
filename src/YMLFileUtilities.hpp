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
 * @file YMLFileUtilities.hpp
 *
 * @brief General functions that are used by the YMLFile class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef YMLFILEUTILITIES_HPP
#define YMLFILEUTILITIES_HPP

#include "Error.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <stdexcept>
#include <string>

/**
 * @brief Utility functions that are not really related to a single class.
 */
namespace YMLFileUtilities {

/**
 * @brief Get a random double precision floating point value in between 0 and 1.
 *
 * @return Random uniform double precision floating point value.
 */
inline double random_double() {
    return ((double)rand()) / ((double)RAND_MAX);
}

/**
 * @brief Split a string of the form [str1, str2, str3] into its parts.
 *
 * @param value std::string having the form mentioned above.
 * @param str1 Variable to store the first part in.
 * @param str2 Variable to store the second part in.
 * @param str3 Variable to store the third part in.
 */
inline void split_string(const std::string& value, std::string& str1,
                         std::string& str2, std::string& str3) {
    size_t pos1 = value.find('[') + 1;
    size_t pos2 = value.find(',', pos1);
    str1 = value.substr(pos1, pos2 - pos1);
    pos1 = pos2 + 1;
    pos2 = value.find(',', pos1);
    str2 = value.substr(pos1, pos2 - pos1);
    pos1 = pos2 + 1;
    pos2 = value.find(']', pos1);
    str3 = value.substr(pos1, pos2 - pos1);
}

/**
 * @brief Convert the given string to a variable of the given template type.
 *
 * @param value std::string value.
 * @return Variable of the given template type containing the parsed contents of
 * the given std::string.
 */
template <typename _datatype_> _datatype_ convert(const std::string& value);

/**
 * @brief Convert the given string to a double precision floating point value.
 *
 * @param value std::string value.
 * @return Double precision floating point stored in the string.
 */
template <> inline double convert<double>(const std::string& value) {
    char* str_end;
    double dvalue = strtod(value.c_str(), &str_end);
    if(str_end == value.c_str()) {
        std::cerr << "Error converting \"" << value
                  << "\" to a floating point value!" << std::endl;
        my_exit();
    }
    return dvalue;
}

/**
 * @brief Convert the given string to an integer value.
 *
 * @param value std::string value.
 * @return Integer stored in the string.
 */
template <> inline int convert<int>(const std::string& value) {
    char* str_end;
    int ivalue = strtol(value.c_str(), &str_end, 0);
    if(str_end == value.c_str()) {
        std::cerr << "Error converting \"" << value << "\" to an integer value!"
                  << std::endl;
        my_exit();
    }
    return ivalue;
}

/**
 * @brief Convert the given string to an unsigned integer value.
 *
 * @param value std::string value.
 * @return Unsigned integer stored in the string.
 */
template <>
inline unsigned int convert<unsigned int>(const std::string& value) {
    char* str_end;
    unsigned int ivalue = strtol(value.c_str(), &str_end, 0);
    if(str_end == value.c_str()) {
        std::cerr << "Error converting \"" << value
                  << "\" to an unsigned integer value!" << std::endl;
        my_exit();
    }
    return ivalue;
}

/**
 * @brief Convert the given string to an unsigned char value.
 *
 * @param value std::string value.
 * @return Unsigned char stored in the string.
 */
template <>
inline unsigned char convert<unsigned char>(const std::string& value) {
    char* str_end;
    unsigned char ivalue = strtol(value.c_str(), &str_end, 0);
    if(str_end == value.c_str()) {
        //    cmac_error("Error converting \"%s\" to an unsigned char value!",
        //               value.c_str());
        my_exit();
    }
    return ivalue;
}

/**
 * @brief Convert the given string to a boolean value.
 *
 * The following string literals map to true: "true", "yes", "on", "y".
 * The following string literals map to false: "false", "no", "off", "n".
 * The string is converted to lowercase before it is parsed, so upper case or
 * mixed case versions, e.g. "True", "FALSE", "oFf" will also be correctly
 * parsed. All other string literals will result in an error.
 *
 * @param value std::string value.
 * @return True or false.
 */
template <> inline bool convert<bool>(const std::string& value) {
    std::string value_copy(value);
    // convert to lowercase
    std::transform(value_copy.begin(), value_copy.end(), value_copy.begin(),
                   ::tolower);
    // strip trailing whitespace
    unsigned int i = 0;
    while(value_copy[i] == ' ') {
        ++i;
    }
    value_copy = value_copy.substr(i);
    i = value_copy.size() - 1;
    while(value_copy[i] == ' ') {
        --i;
    }
    value_copy = value_copy.substr(0, i + 1);
    if(value_copy == "true" || value_copy == "yes" || value_copy == "on" ||
       value_copy == "y") {
        return true;
    } else if(value_copy == "false" || value_copy == "no" ||
              value_copy == "off" || value_copy == "n") {
        return false;
    } else {
        //    cmac_error("Error converting \"%s\" to a boolean value!",
        //               value_copy.c_str());
        my_exit();
        return false;
    }
}

/**
 * @brief Convert the given value to a std::string.
 *
 * @param value Value to convert.
 * @return std::string.
 */
template <typename _datatype_> std::string to_string(_datatype_ value) {
    std::stringstream sstream;
    sstream << value;
    return sstream.str();
}

/**
 * @brief to_string specialization for unsigned char values.
 *
 * We have to convert the unsigned char to an integer before outputting,
 * otherwise it is outputted as the character it is supposed to represent, which
 * yields garbage.
 *
 * @param value Value to convert.
 * @return std::string.
 */
template <> inline std::string to_string<unsigned char>(unsigned char value) {
    std::stringstream sstream;
    unsigned int ivalue = value;
    sstream << ivalue;
    return sstream.str();
}

/**
 * @brief to_string specialization for boolean values.
 *
 * @param value Bool value.
 * @return "true" or "false".
 */
template <> inline std::string to_string<bool>(bool value) {
    if(value) {
        return "true";
    } else {
        return "false";
    }
}

/**
 * @brief Split the given string containing a value and an associated unit into
 * a std::pair.
 *
 * @param svalue std::string containing a value - unit pair.
 * @return std::pair containing the value and unit.
 */
inline std::pair<double, std::string> split_value(const std::string& svalue) {
    size_t idx;
    double value;
    try {
        value = std::stod(svalue, &idx);
    } catch(std::invalid_argument e) {
        //    cmac_error("Error extracting value from \"%s\" unit-value pair!",
        //               svalue.c_str());
        my_exit();
    }

    while(svalue[idx] == ' ') {
        ++idx;
    }
    std::string unit = svalue.substr(idx);
    return make_pair(value, unit);
}

/**
 * @brief Get the index of the last element in the given ordered array that is
 * smaller than the given value.
 *
 * This routine uses bisection, and always returns a value in the range
 * [0, length-2], even if the given value is outside the given array.
 *
 * @param x Value to locate.
 * @param xarr Array in which to search.
 * @param length Length of the array.
 * @return Index of the last element in the ordered array that is smaller than
 * the given value, i.e. value is in between xarr[index] and xarr[index+1].
 */
inline unsigned int locate(double x, const double* xarr, unsigned int length) {
    unsigned int jl = 0;
    unsigned int ju = length + 1;
    while(ju - jl > 1) {
        unsigned int jm = (ju + jl) / 2;
        if(x > xarr[jm]) {
            jl = jm;
        } else {
            ju = jm;
        }
    }
    if(jl == length - 1) {
        --jl;
    }
    return jl;
}

/**
 * @brief Compose a filename made up by the given prefix and counter value,
 * appropriately zero padded.
 *
 * @param folder Folder to add to the filename.
 * @param prefix Prefix for the filename.
 * @param extension Extension for the filename.
 * @param counter Value of the counter.
 * @param padding Number of digits the counter should have.
 * @return std::string with format: "<prefix>XX<counter>XX.<extension>", where
 * the number of Xs is equal to padding.
 */
inline std::string compose_filename(const std::string& folder,
                                    const std::string& prefix,
                                    const std::string& extension,
                                    unsigned int counter,
                                    unsigned int padding) {
    std::stringstream namestring;
    if(!folder.empty()) {
        namestring << folder << "/";
    }
    namestring << prefix;
    namestring.fill('0');
    namestring.width(padding);
    namestring << counter;
    namestring << "." << extension;
    return namestring.str();
}

/**
 * @brief Get the absolute path corresponding to the given path.
 *
 * @param path Path, can be relative or absolute.
 * @return Absolute path.
 */
inline std::string get_absolute_path(std::string path) {
    // strip trailing / from path
    if(path[path.size() - 1] == '/') {
        path = path.substr(0, path.size() - 1);
    }

    char* absolute_path_ptr = realpath(path.c_str(), nullptr);
    if(absolute_path_ptr == nullptr) {
        //    cmac_error("Unable to resolve path \"%s\"!", path.c_str());
        my_exit();
    }
    std::string absolute_path(absolute_path_ptr);
    free(absolute_path_ptr);
    return absolute_path;
}

/**
 * @brief Get a time stamp of the form day/month/year, hour:minutes:seconds.
 *
 * @return Time stamp.
 */
inline std::string get_timestamp() {
    std::time_t timestamp = std::time(nullptr);
    std::tm* time = std::localtime(&timestamp);
    std::stringstream timestream;
    if(time->tm_mday < 10) {
        timestream << "0";
    }
    timestream << time->tm_mday << "/";
    if(time->tm_mon < 9) {
        timestream << "0";
    }
    // tm_mon counts from 0 to 11
    // tm_year is the number of years since 1900
    timestream << (time->tm_mon + 1) << "/" << (time->tm_year + 1900) << ", ";
    if(time->tm_hour < 10) {
        timestream << "0";
    }
    timestream << time->tm_hour << ":";
    if(time->tm_min < 10) {
        timestream << "0";
    }
    timestream << time->tm_min << ":";
    if(time->tm_sec < 10) {
        timestream << "0";
    }
    timestream << time->tm_sec;
    return timestream.str();
}

/**
 * @brief Convert a floating point time value in seconds to a human readable
 * time string.
 *
 * If the time is larger than a minute, it is displayed as Xm Ys, and so on for
 * larger subdivisions of time, up to years.
 *
 * @param time Double precision floating point time value.
 * @return std::string containing a human readable version of the time value.
 */
inline std::string human_readable_time(double time) {
    std::stringstream timestream;
    // 2^32 minutes is 8166 years. We can safely assume no internal timer will
    // ever reach that value
    unsigned int minutes = time / 60.;
    double seconds = time - 60. * minutes;
    if(minutes > 0) {
        unsigned int hours = minutes / 60;
        if(hours > 0) {
            minutes -= 60 * hours;
            unsigned int days = hours / 24;
            if(days > 0) {
                hours -= 24 * days;
                unsigned int years = days / 365;
                if(years > 0) {
                    days -= 365 * years;
                    timestream << years << "y ";
                }
                timestream << days << "d ";
            }
            timestream << hours << "h ";
        }
        timestream << minutes << "m ";
    }
    timestream << seconds << "s";
    return timestream.str();
}

/**
 * @brief Get the unit of \f$2^{10e}\f$ bytes.
 *
 * @param exponent Exponent \f$e\f$.
 * @return Name for \f$2^{10e}\f$ bytes: (\f$2^{10}\f$ bytes = KB, ...).
 */
inline std::string byte_unit(unsigned char exponent) {
    switch(exponent) {
        case 0:
            return "bytes";
        case 1:
            return "KB";
        case 2:
            return "MB";
        case 3:
            return "GB";
        case 4:
            return "TB";
        default:
            //    cmac_error("Name for 2^(10*%u) bytes was not implemented!",
            //    exponent);
            my_exit();
            return "";
    }
}

/**
 * @brief Convert the given number of bytes to a human readable string.
 *
 * @param bytes Number of bytes.
 * @return std::string containing the given number of bytes in "bytes", "KB",
 * "MB", "GB"...
 */
inline std::string human_readable_bytes(unsigned long bytes) {
    unsigned char sizecount = 0;
    double bytefloat = bytes;
    while((bytes >> 10) > 0) {
        bytes >>= 10;
        ++sizecount;
        bytefloat /= 1024.;
    }
    std::stringstream bytestream;
    bytefloat = std::round(100. * bytefloat) * 0.01;
    bytestream << bytefloat << " " << byte_unit(sizecount);
    return bytestream.str();
}
}

#endif  // YMLFILEUTILITIES_HPP
