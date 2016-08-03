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
 * @file CubicBox.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief Containers for the Delaunay tesselation: general header
 *
 * The CubicBox represents a cubic (or square) box with given origin and side.
 */
#ifndef CUBICBOX_HPP
#define CUBICBOX_HPP

#include "DelCont.hpp"
#include "Vec.hpp"
#include "utilities/Cuboid.hpp"
#include <vector>

class Particle;
class VorGen;

/**
 * @brief A cubic box
 *
 * The edges of the box are parallel to the coordinate axes.
 */
class CubicBox : public DelCont {
  private:
    /*! @brief The origin of the CubicBox (the coordinates of the center of the
     *  box) */
    Vec _origin;

    /*! @brief The side of the CubicBox (the side of the box) */
    double _side;

    /*! @brief The bitwidth of the box (1 with n 0's, where n is 60 divided by
     *  the number of dimensions; n = 20 for 3D and n = 30 for 2D) */
    unsigned long _bitwidth;

  public:
    CubicBox(Vec origin, double side);
    CubicBox();
    virtual ~CubicBox() {}
    bool inside(VorGen* point, double radius);
    bool inside(Vec pntpos);
    std::vector<double> get_bounding_tetrahedron();
    void get_bounding_box(double* box);
    double get_box_width();
    unsigned long get_key(Vec& coords);
    void get_ngb_keys(Vec& coords, double radius, unsigned long* keys);
    void keep_inside(Particle* p);
    void closest_copy(Vec& v);

    Cuboid get_cuboid();
};

#endif  // CUBICBOX_HPP
