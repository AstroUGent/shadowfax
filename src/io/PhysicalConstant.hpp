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
 * @file PhysicalConstant.hpp
 *
 * @brief Physical constant representation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PHYSICALCONSTANT_HPP
#define PHYSICALCONSTANT_HPP

#include "Unit.hpp"
#include "UnitConverter.hpp"
#include "UnitSet.hpp"

/**
  * \brief Representation of a physical constant.
  *
  * The constant has a value and a Unit, such that we can easily
  * convert from one set of units to another.
  */
class PhysicalConstant {
  private:
    /*! \brief Value of the physical constant */
    double _value;

    /*! \brief Unit of the physical constant */
    Unit _unit;

  public:
    /**
      * \brief Construct a physical constant with a Unit and a value
      *
      * @param value Value of the physical constant
      * @param unit Unit of the physical constant
      */
    inline PhysicalConstant(double value, Unit unit)
            : _value(value), _unit(unit) {}

    /**
      * \brief Get the raw value of the physical constant
      *
      * This function totally disregards the Unit of the physical constant, so
      * it is unsafe to use.
      *
      * @return Raw value of the physical constant in its original Unit
      */
    inline double get_value() { return _value; }

    /**
      * \brief Get the raw value of the Unit associated to this physical
      * constant
      *
      * This function can be used together with the get_value() function to make
      * it more safe.
      *
      * @return The original Unit associated with the physical constant
      */
    inline Unit get_unit() { return _unit; }

    /**
      * \brief Get the value of the physical constant in the specified system of
      * units
      *
      * The function determines what the Unit for the physical constant is in
      * the specified UnitSet by making use of the quantity of its original
      * Unit. The value of the physical constant is then converted using a
      * UnitConverter and this value is returned.
      *
      * @param unitset UnitSet specifying the system of units in which you want
      * to use the physical constant
      * @return The value of the physical constant in the desired system of
      * units
      */
    inline double get_value(UnitSet unitset) {
        Unit unit = unitset.get_unit(_unit.get_quantity());
        UnitConverter converter(_unit, unit);
        return converter.convert(_value);
    }
};

#endif  // PHYSICALCONSTANT_HPP
