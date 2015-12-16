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
 * @file Cuboid.hpp
 *
 * @brief Representation of a Cuboid, a box with an origin and different side
 * lengths in every dimension
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef CUBOID_HPP
#define CUBOID_HPP

#include "Vec.hpp"
#include "RestartFile.hpp"

/**
 * @brief Representation of a 2D or 3D cuboid
 *
 * The cuboid consists of two Vecs: an anchor that gives the corner of the
 * cuboid with the lowest coordinates and a side Vec with the lengths of
 * the 2 or 3 sides.
 */
class Cuboid{
private:
    /*! \brief Bottom left (front) corner of the box */
    Vec _anchor;

    /*! \brief Sides of the box */
    Vec _sides;

public:
    inline Cuboid() {}

    /**
     * @brief Constructor
     *
     * @param anchor Bottom left (front) corner of the box
     * @param sides Sides of the box
     */
    inline Cuboid(Vec anchor, Vec sides) : _anchor(anchor), _sides(sides) {}

    /**
     * @brief Set the extents of the box
     *
     * @param anchor Bottom left (front) corner of the box
     * @param sides Sides of the box
     */
    inline void set(Vec anchor, Vec sides){
        _anchor = anchor;
        _sides = sides;
    }

    /**
     * @brief Get the bottom left (front) corner of the box
     * @return Bottom left (front) corner of the box
     */
    inline Vec get_anchor(){
        return _anchor;
    }

    /**
     * @brief Get the sides of the box
     * @return Sides of the box
     */
    inline Vec get_sides(){
        return _sides;
    }

    /**
     * @brief Dump the box to the given RestartFile
     *
     * @param rfile RestartFile to write to
     */
    inline void dump(RestartFile &rfile){
        rfile.write(_anchor);
        rfile.write(_sides);
    }

    /**
     * @brief Restart constructor. Initialize the box based on the given
     * RestartFile
     *
     * @param rfile RestartFile to read from
     */
    inline Cuboid(RestartFile &rfile){
        rfile.read(_anchor);
        rfile.read(_sides);
    }
};

#endif // CUBOID_HPP
