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
 * @file Unit.hpp
 *
 * @brief Unit support
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef UNIT_HPP
#define UNIT_HPP

#include <istream>
#include <ostream>
#include <string>

#include "RestartFile.hpp"

/**
  * \brief Basic unit abstraction.
  *
  * A unit is defined by a physical quantity, an arbitrary unit name and the
  * SI-value of that unit.
  *
  * For example: kiloparsecs:
  *     Unit("length", "kpc", 3.08567758e19)
  *
  * So the idea is that once you know which quantity you have, you take the
  * basic SI unit for that quantity and then the given value is the value of the
  * new unit in basic SI units. For the example above, the basic unit is meters.
  *
  * As a reminder, the seven basic SI quantities and their respective units are:
  *  - length: meters (m)
  *  - mass: kilogram (kg)
  *  - time: seconds (s)
  *  - temperature: Kelvin (K)
  *  - electric current: Ampere (A)
  *  - luminous intensity: candela (cd)
  *  - amount of substance: mole (mol)
  */
class Unit {
  protected:
    /*! \brief Name of the Unit */
    std::string _name;

    /*! \brief Quantity of the Unit, expressed in basic SI quantities */
    std::string _quantity;

    /*! \brief SI value of the Unit in basic SI units for the quantity of this
     *  Unit */
    double _SI_value;

  public:
    /**
      * @brief Default constructor
      *
      * An empty Unit has name "dimensionless" and value 1. (which makes sure
      * that conversion will not change the value of the quantity).
      */
    Unit() : _name("dimensionless"), _quantity(""), _SI_value(1.) {}

    /**
     * @brief Constructor
     *
     * @param quantity SI-quantity or combination that makes clear what quantity
     * this unit represents
     * @param name Name of the unit, that might be meaningful to a subset of
     * people
     * @param SI_value Value of the unit if you express 1 \a unit \a name in
     * SI-units
     */
    Unit(std::string quantity, std::string name, double SI_value)
            : _name(name), _quantity(quantity), _SI_value(SI_value) {}

    ~Unit() {}

    /**
      * @brief Get the SI value of this unit in the basic SI unit for that
      * quantity
      *
      * @returns The SI value of the Unit
      */
    inline double get_SI_value() { return _SI_value; }

    /**
      * @brief Get the quantity associated with this Unit
      *
      * The quantity should be expressed in terms of the seven basic SI units
      * (see above) by means of multiplication and division. For convenience,
      * quantities in this expression should be ordered: first all quantities
      * that are multiplied (or 1 if there is no such quantity), then the ones
      * that are divided. These two groups then should be ordered
      * alphabetically. Powers are not allowed and every quantity that is
      * divided should be preceded by a backslash (/). All quantities that are
      * multiplied (except for the first one) should be preceded by a *.
      *
      * @returns The name of the quantity of the Unit
      */
    inline const std::string get_quantity() { return _quantity; }

    /**
      * @brief Get the name of the Unit
      *
      * e.g. "meters" or "kilogram"
      *
      * Symbols are also allowed. This field is ignored in all operations and is
      * provided to make the Unit more human readable. If your custom length
      * Unit e.g. is kiloparsec, this is where you put "kiloparsec" or "kpc".
      *
      * @returns The name of the Unit
      */
    inline const std::string get_name() { return _name; }

    /**
      * @brief void pointer to the SI value needed for C-style writing routines
      *
      * @returns A void pointer to the internal _SI_value
      */
    inline void* SI_value() { return &_SI_value; }

    /**
      * @brief void pointer to the name needed for C-style writing routines
      *
      * @returns A void pointer to the internal _name
      */
    inline const void* name() { return &_name; }

    /**
      * @brief void pointer to the quantity needed for C-style writing routines
      *
      * @returns A void pointer to the internal _quantity
      */
    inline const void* quantity() { return &_quantity; }

    /**
      * @brief Divide this Unit by another Unit
      *
      * The name and quantity are just put together with a backslash in between.
      * The SI value is divided.
      *
      * @param u Unit to divide this Unit by
      * @returns Reference to this Unit
      */
    inline Unit& operator/=(Unit u) {
        _quantity += "/" + u._quantity;
        _name += "/" + u._name;
        _SI_value /= u._SI_value;
        return *this;
    }

    /**
      * @brief Multiply this Unit with another Unit
      *
      * The name and quantity are just put together with a * in between.
      * The SI value is multiplied.
      *
      * @param u Unit to multiply this Unit with
      * @returns Reference to this Unit
      */
    inline Unit& operator*=(Unit u) {
        _quantity += "*" + u._quantity;
        _name += "*" + u._name;
        _SI_value *= u._SI_value;
        return *this;
    }

    /**
     * @brief Dump the unit to the given RestartFile
     *
     * @param rfile RestartFile to write to
     */
    inline void dump(RestartFile& rfile) {
        rfile.write(_name);
        rfile.write(_quantity);
        rfile.write(_SI_value);
    }

    /**
     * @brief Restart constructor. Initialize the unit from the given
     * RestartFile
     *
     * @param rfile RestartFile to write to
     */
    inline Unit(RestartFile& rfile) {
        rfile.read(_name);
        rfile.read(_quantity);
        rfile.read(_SI_value);
    }
};

/**
  * @brief Divide two units
  *
  * @param a Unit in numerator of the division
  * @param b Unit in denominator of the division
  * @returns A new Unit which is the result of the division
  */
inline Unit operator/(Unit a, Unit b) { return a /= b; }

/**
  * @brief Multiply two units
  *
  * @param a Unit one
  * @param b Unit two
  * @returns A new Unit which is the result of the multiplication
  */
inline Unit operator*(Unit a, Unit b) { return a *= b; }

#endif  // UNIT_HPP
