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
 * @file Box.hpp
 *
 * @brief Representation of a cubic box
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef BOX_HPP
#define BOX_HPP

#include "Vec.hpp"

/**
 * @brief A 2D square or 3D cube.
 *
 * Stores an anchor and a side. A Box is similar to a Cuboid, but with a side
 * that is equal for all dimensions.
 */
class Box {
  private:
    /*! \brief Bottom left (front) corner of the box */
    Vec _anchor;

    /*! \brief Side length of the box */
    double _side;

  public:
    /**
     * @brief Constructor. Initialize a box width zero extents
     */
    inline Box() : _side(0.) {}

    /**
     * @brief Constructor
     *
     * @param anchor Bottom left (front) corner of the box
     * @param side Side length of the box
     */
    inline Box(Vec anchor, double side) : _anchor(anchor), _side(side) {}

    /**
     * @brief Get the bottom left (front) corner of the box
     *
     * @return Bottom left (front) corner of the box
     */
    inline Vec get_anchor() {
        return _anchor;
    }

    /**
     * @brief Get the side length of the box
     *
     * @return Side length of the box
     */
    inline double get_side() {
        return _side;
    }
};

#endif  // BOX_HPP
