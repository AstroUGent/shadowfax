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
 * @file UnitConverter.hpp
 *
 * @brief Convert between units
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef UNITCONVERTER_HPP
#define UNITCONVERTER_HPP

#include "Unit.hpp"  // for Unit
#include <iostream>  // for basic_ostream, operator<<, etc
#include <string>    // for operator<<, string

/**
  * @brief Convert from one Unit to another.
  *
  * This is of course only possible if both units are compatible, i.e. if the
  * quantities of both units are the same.
  */
class UnitConverter {
  private:
    /*! @brief Conversion factor that converts the SI-value of one Unit to that
     *  of another */
    double _conversion_factor;

  public:
    /**
      * @brief Constructor
      *
      * Checks if the units are compatible and calculates the conversion factor
      * between the units.
      *
      * @param in Original Unit, i.e. the Unit you feed to the converter
      * @param out New Unit, what comes out of the converter
      */
    UnitConverter(Unit in, Unit out) {
        if(in.get_quantity().compare(out.get_quantity())) {
            std::cerr << in.get_quantity() << "\t" << out.get_quantity()
                      << std::endl;
            std::cerr << in.get_quantity().compare(out.get_quantity())
                      << std::endl;
            std::cerr << "Units are not compatible!" << std::endl;
            my_exit();
        }
        _conversion_factor = in.get_SI_value() / out.get_SI_value();
    }

    ~UnitConverter() {}

    /**
      * @brief Convert the given value from the input Unit to the output Unit
      *
      * @param value Input value, in original Unit
      * @returns Converted value in the new Unit
      */
    double convert(const double value) {
        return value * _conversion_factor;
    }
};

#endif  // UNITCONVERTER_HPP
