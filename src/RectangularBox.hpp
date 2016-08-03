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
 * @file RectangularBox.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief Containers for the Delaunay tesselation: general header
 */
#ifndef RECTANGULARBOX_HPP
#define RECTANGULARBOX_HPP

#include "DelCont.hpp"
#include "Vec.hpp"
#include "utilities/Cuboid.hpp"
#include <vector>

class Particle;
class VorGen;
class RestartFile;

/**
 * @brief DelCont implementation consisting of a rectangular shaped box
 *
 * The box can have different side lengths in different dimensions.
 */
class RectangularBox : public DelCont {
  private:
    /*! @brief Center of the cuboid */
    Vec _center;

    /*! @brief Sides of the cuboid */
    Vec _sides;

    /*! @brief Maximum size of an integer coordinate in one of the dimensions */
    unsigned long _bitwidth;

    /*! @brief Largest side length of the cuboid */
    double _maxside;

  public:
    RectangularBox(Vec center, Vec sides);
    RectangularBox();
    virtual ~RectangularBox() {}

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

    void dump(RestartFile& rfile);
    RectangularBox(RestartFile& rfile);
};

#endif  // RECTANGULARBOX_HPP
