/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ParameterFile.hpp
 *
 * @brief Parameter file: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PARAMETERFILE_HPP
#define PARAMETERFILE_HPP

#include "YMLFile.hpp"
#include "YMLFileUtilities.hpp"

#include <boost/property_tree/ptree.hpp>  // for basic_ptree
#include <sstream>  // for basic_stringbuf<>::int_type, etc
#include <string>   // for string, operator<

class RestartFile;
class Unit;
class UnitSet;

/**
 * @brief Abstraction of the parameter file that contains vital run information
 *
 * The parameter file is a .ini-file which contains parameters that set the code
 * behaviour at runtime. Every parameter has a default value which is used if
 * the parameter is not provided.
 *
 * A complete parameter file would look like this:
 * \verbinclude default.ini
 */
class ParameterFile {
  private:
    /*! @brief Property tree containing the contents of the .ini file */
    boost::property_tree::ptree _parameters;

    /*! @brief YMLFile containing the contents of the .yml file. */
    YMLFile* _yml_file;

    /*! @brief UnitSet containing the internal units. Used to convert physical
     *  quantities to the correct units. */
    UnitSet* _internal_units;

    void print_contents();

  public:
    ParameterFile(std::string name);

    ~ParameterFile();

    /**
     * @brief Get the value of the parameter with the given name
     *
     * @param name Name of the parameter
     * @param default_value Default value for the parameter
     * @return Value of the parameter
     */
    template <typename T> T get_parameter(std::string name, T default_value) {
        if(_yml_file) {
            // convert all '.' to ':' to comply with the YML syntax
            YMLFileUtilities::replace_char('.', ':', name);
            return _yml_file->get_value<T>(name, default_value);
        } else {
            // convert all ':' to '.' to comply with the boost::ini syntax
            YMLFileUtilities::replace_char(':', '.', name);
            return _parameters.get<T>(name, default_value);
        }
    }

    /**
     * @brief Get the value of the parameter with the given name
     *
     * This version uses a C style string.
     *
     * @param name Name of the parameter
     * @param default_value Default value for the parameter
     * @return Value of the parameter
     */
    template <typename T> T get_parameter(const char* name, T default_value) {
        return get_parameter<T>(std::string(name), default_value);
    }

    double get_quantity(std::string name, std::string quantity,
                        std::string default_value);

    double get_quantity(std::string name, Unit unit, std::string default_value);

    /**
     * @brief Check if the given parameter is true or false/not specified
     *
     * @param name Parameter name
     * @return True only if the parameter is specified and has a "true" value
     */
    inline bool check_parameter(const char* name) {
        return get_parameter<bool>(name, false);
    }

    void dump(RestartFile& rfile);
    ParameterFile(RestartFile& rfile);
};

#endif  // PARAMETERFILE_HPP
