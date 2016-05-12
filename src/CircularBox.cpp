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
 * \file CircularBox.cpp
 *
 * @brief Circular box implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "DelCont.hpp"
#include "VorGen.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

/**
 * Initialize a CircularBox with given origin and radius.
 * @param origin Coordinates of the origin
 * @param radius The radius of the sphere
 */
CircularBox::CircularBox(Vec origin, double radius) {
    _origin = origin;
    _radius = radius;
    // storing the radius squared makes internal checks faster
    _radius2 = radius * radius;
}

/**
 * @brief Determine whether the given coordinates define a point inside the
 * sphere.
 *
 * Points on the boundary of the sphere are defined to be outside the sphere.
 *
 * @param pntpos A Vec holding the coordinates of a point
 * @return True if the radius of the given coordinates is smaller than the
 * radius of the sphere
 */
bool CircularBox::inside(Vec pntpos) {
    double radius = 0.0;
    for(unsigned i = ndim_; i--;) {
        pntpos[i] -= _origin[i];
        radius += pntpos[i] * pntpos[i];
    }
    return radius < _radius2;
}

/**
 * @brief Determine whether the volume defined by the given Point and radius
 * around the Point is completely inside the CircularBox.
 *
 * Points on the boundary are defined to be outside.
 *
 * @param point Point that is the center of the spherical region that has to be
 * tested
 * @param radius Radius of the spherical region to be tested
 * @return True if the complete region specified by the Point and radius are
 * strictly inside the sphere
 */
bool CircularBox::inside(VorGen* point, double radius) {
    double r2 = 0;
    for(unsigned int i = ndim_; i--;) {
        double diff = point->pos(i) - _origin[i];
        r2 += diff * diff;
    }
    return radius + sqrt(r2) < _radius;
}

/**
 * @brief Get the coordinates of the Simplex that encompasses the complete
 * volume inside the CircularBox
 *
 * To ensure that mirror points are also inside this Simplex, the boundary
 * points are chosen with a serious margin.
 *
 * @return A vector with coordinate-vectors of 4 points.
 */
vector<double> CircularBox::get_bounding_tetrahedron() {
    vector<double> coords((ndim_ + 1) * ndim_);
    double r = 4 * _radius;
#if ndim_ == 3
    // equilateral tetrahedron symmetrically enclosing a sphere with
    // radius==4*_radius at its center
    coords[0] = _origin[0];
    coords[1] = _origin[1] + 3 * r;
    coords[2] = _origin[2];

    coords[3] = _origin[0] + sqrt(8.0) * r;
    coords[4] = _origin[1] - r;
    coords[5] = _origin[2];

    coords[6] = _origin[0] - sqrt(2.0) * r;
    coords[7] = _origin[1] - r;
    coords[8] = _origin[2] + sqrt(6.0) * r;

    coords[9] = _origin[0] - sqrt(2.0) * r;
    coords[10] = _origin[1] - r;
    coords[11] = _origin[2] - sqrt(6.0) * r;
#else
    // equilateral triangle symmetrically enclosing a circle with
    // radius==4*_radius at its center
    double side = sqrt(3.0) * r;
    coords[0] = _origin[0] - side;
    coords[1] = _origin[1] - r;

    coords[2] = _origin[0] + side;
    coords[3] = _origin[1] - r;

    coords[4] = _origin[0];
    coords[5] = _origin[1] + 2 * r;
#endif
    return coords;
}

/**
  * @brief Retrieve the origin and side of a box containing the entire sphere or
  * circle.
  *
  * @param box A 3- or 4-element array to store the radius and side in
  */
void CircularBox::get_bounding_box(double* box) {
    for(unsigned int i = ndim_; i--;) {
        box[i] = _origin[i] - _radius;
    }
    box[ndim_] = 2. * _radius;
}

/**
 * @brief Return the width of the box encompassing the entire circle/sphere.
 *
 * This is needed for tree construction.
 *
 * @return 2 times the radius of the sphere
 */
double CircularBox::get_box_width() {
    return 2 * _radius;
}

/**
  * @brief Calculate the key for the given coordinates.
  *
  * \attention Not implemented (yet)!
  *
  * @param coords A Vec containing the coordinates of a particle inside the
  * CircularBox
  * @return The hilbert key for the given coordinates, taking into account the
  * size of the internal volume
  */
unsigned long CircularBox::get_key(Vec& coords) {
    return 0;
}

/**
  * @brief Keep the given particle inside the CircularBox, using periodicity.
  *
  * \attention Not implemented (yet)!
  *
  * @param p A particle, not necessarily inside the CircularBox
  */
void CircularBox::keep_inside(Particle* p) {
    // do nothing
}

/**
 * @brief Get the cuboid that encompasses the entire CircularBox
 *
 * @warning Dummy function with dummy return value
 *
 * @return Cuboid specifying the dimensions of the surrouning box
 */
Cuboid CircularBox::get_cuboid() {
    Cuboid cuboid;
    return cuboid;
}
