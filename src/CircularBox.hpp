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
 * @file CircularBox.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief Containers for the Delaunay tesselation: general header
 *
 * The CircularBox represents a round or spheric volume with given central point
 * and radius.
 */
#ifndef CIRCULARBOX_HPP
#define CIRCULARBOX_HPP

#include "DelCont.hpp"  // for DelCont
#include "Vec.hpp"
#include "utilities/Cuboid.hpp"
#include <vector>

class Particle;
class VorGen;

/**
 * @brief A spherical volume
 *
 * A sphere in 3D or a circle in 2D.
 */
class CircularBox : public DelCont {
  private:
    /*! @brief The origin of the sphere or circle */
    Vec _origin;

    /*! @brief The radius of the sphere or circle */
    double _radius;

    /*! @brief The radius squared of the sphere or circle (for convenience) */
    double _radius2;

  public:
    CircularBox(Vec origin, double radius);
    virtual ~CircularBox() {}
    bool inside(VorGen* point, double radius);
    bool inside(Vec pntpos);
    std::vector<double> get_bounding_tetrahedron();
    void get_bounding_box(double* box);
    double get_box_width();
    unsigned long get_key(Vec& coords);

    /**
     * @brief Empty
     * @param coords
     * @param radius
     * @param keys
     */
    void get_ngb_keys(Vec& coords, double radius, unsigned long* keys) {}
    void keep_inside(Particle* p);

    /**
     * @brief Empty
     * @param v
     */
    void closest_copy(Vec& v) {}

    Cuboid get_cuboid();
};

#endif  // CIRCULARBOX_HPP
