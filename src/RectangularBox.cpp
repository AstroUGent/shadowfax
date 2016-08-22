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
 * @file RectangularBox.cpp
 *
 * @brief Rectangular box implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "RectangularBox.hpp"
#include "ProgramLog.hpp"          // for LOGS
#include "RestartFile.hpp"         // for RestartFile
#include "VorGen.hpp"              // for VorGen
#include "utilities/Hilbert.hpp"   // for get_key
#include "utilities/Particle.hpp"  // for Particle
#include <algorithm>               // for max
#include <cmath>                   // for fabs, sqrt
using namespace std;

/**
 * @brief Constructor
 *
 * @param center Center of the box
 * @param sides Side lengths of the box
 */
RectangularBox::RectangularBox(Vec center, Vec sides)
        : _center(center), _sides(sides) {
    _bitwidth = 1;
    _bitwidth <<= (60 / ndim_);
    _maxside = std::max(sides[0], sides[1]);
#if ndim_ == 3
    _maxside = std::max(_maxside, sides[2]);
#endif

    LOGS("RectangularBox created from origin and sides");
}

/**
 * @brief Empty construcctor
 */
RectangularBox::RectangularBox() {
    _bitwidth = 1;
    _bitwidth <<= (60 / ndim_);
    _maxside = 0.;

    LOGS("Empty RectangularBox created");
}

/**
 * @brief Check if the sphere centered around the given point and with the given
 * radius is inside the box
 *
 * @param point VorGen around which the sphere is centered
 * @param radius Radius of the sphere
 * @return True if the sphere is inside the box or touches its border, false
 * otherwise
 */
bool RectangularBox::inside(VorGen* point, double radius) {
    Vec pntpos = point->get_position();
    for(unsigned int i = 0; i < ndim_; i++) {
        double test = _sides[i] / 2;
        if(fabs(_center[i] - pntpos[i]) > test - radius) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Check if the given Vec lies inside the box
 *
 * @param pntpos Vec specifying a set of coordinates
 * @return True if the coordinates lie inside the box, false otherwise
 */
bool RectangularBox::inside(Vec pntpos) {
    for(unsigned int i = 0; i < ndim_; i++) {
        double test = _sides[i] / 2;
        if(fabs(_center[i] - pntpos[i]) > test) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Get the vertices of a simplex encompassing the entire box
 *
 * @return Coordinates of the vertices of a simplex that encompasses the box
 */
vector<double> RectangularBox::get_bounding_tetrahedron() {
    vector<double> coords((ndim_ + 1) * ndim_);
    double side = 0.;
    for(unsigned int i = 0; i < ndim_; i++) {
        side = std::max(side, _sides[i]);
    }
#if ndim_ == 3
    // r = 4 * radius of sphere that contains cube
    double r = sqrt(12.) * side;
    coords[0] = _center[0];
    coords[1] = _center[1] + 3 * r;
    coords[2] = _center[2];

    coords[3] = _center[0] + sqrt(8.0) * r;
    coords[4] = _center[1] - r;
    coords[5] = _center[2];

    coords[6] = _center[0] - sqrt(2.0) * r;
    coords[7] = _center[1] - r;
    coords[8] = _center[2] + sqrt(6.0) * r;

    coords[9] = _center[0] - sqrt(2.0) * r;
    coords[10] = _center[1] - r;
    coords[11] = _center[2] - sqrt(6.0) * r;
#else
    // r = 4 * radius of circle that contains square
    double r = sqrt(8.) * side;
    double l = sqrt(3.) * r;
    coords[0] = _center[0] - l;
    coords[1] = _center[1] - r;

    coords[2] = _center[0] + l;
    coords[3] = _center[1] - r;

    coords[4] = _center[0];
    coords[5] = _center[1] + 2 * r;
#endif
    return coords;
}

/**
 * @brief Get the origin and extents of the box
 *
 * @param box Array to store the result in
 */
void RectangularBox::get_bounding_box(double* box) {
    for(int i = 0; i < ndim_; i++) {
        box[i] = _center[i] - 0.5 * _sides[i];
        box[ndim_ + i] = _sides[i];
    }
}

/**
 * @brief Get the width of the box
 *
 * @warning This quantity is not defined for this type of box and this function
 * should not be used!
 *
 * @return 0, since this function should not be used
 */
double RectangularBox::get_box_width() {
    // problem!!!!
    return 0.;
}

/**
 * @brief Get the Hilbert-key for the given coordinates
 *
 * @param coords Vec specifying a set of coordinates
 * @return The Hilbert-key for the given coordinates
 */
unsigned long RectangularBox::get_key(Vec& coords) {
    unsigned long bits[ndim_] = {0};
    for(unsigned int i = ndim_; i--;) {
        bits[i] = ((coords[i] - _center[i] + 0.5 * _maxside) / _maxside) *
                  _bitwidth;
    }
    return HB::get_key(bits, 60);
}

/**
 * @brief Dummy method
 *
 * @warning Not (to be) used.
 *
 * @param coords XXX
 * @param radius XXX
 * @param keys XXX
 */
void RectangularBox::get_ngb_keys(Vec& coords, double radius,
                                  unsigned long* keys) {
    // not used, so not implemented...
}

/**
 * @brief Keep the given Particle inside the box by periodically trimming its
 * coordinates
 *
 * @param p Particle that resides inside or outside the box
 */
void RectangularBox::keep_inside(Particle* p) {
    for(unsigned int i = 0; i < ndim_; i++) {
        if(p->pos(i) < _center[i] - 0.5 * _sides[i]) {
            p->get_position()[i] += _sides[i];
        }
        // if both boundaries are inclusive, this gives rise to problems
        // (because the boundaries are their own periodic copies)
        // hence the >=
        if(p->pos(i) >= _center[i] + 0.5 * _sides[i]) {
            p->get_position()[i] -= _sides[i];
        }
    }
}

/**
 * @brief Get the closest copy of the given coordinates to the origin (0,0,0)
 *
 * @param v Given coordinates
 */
void RectangularBox::closest_copy(Vec& v) {
    for(unsigned int i = 0; i < ndim_; i++) {
        if(v[i] > 0.5 * _sides[i]) {
            v[i] -= _sides[i];
        }
        if(v[i] < -0.5 * _sides[i]) {
            v[i] += _sides[i];
        }
    }
}

/**
 * @brief Get the cuboid with the extents of this box
 *
 * @return Cuboid specifying the dimensions of the box
 */
Cuboid RectangularBox::get_cuboid() {
#if ndim_ == 3
    Vec origin(_center[0], _center[1], _center[2]);
    Vec sides(_sides[0], _sides[1], _sides[2]);
#else
    Vec origin(_center[0], _center[1]);
    Vec sides(_sides[0], _sides[1]);
#endif
    Cuboid cuboid(origin - 0.5 * sides, sides);
    return cuboid;
}

/**
 * @brief Dump the box to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void RectangularBox::dump(RestartFile& rfile) {
    rfile.write(_center);
    rfile.write(_sides);
    rfile.write(_bitwidth);
    rfile.write(_maxside);

    LOGS("RectangularBox dumped");
}

/**
 * @brief Restart constructor. Initialize the box from the given RestartFile
 *
 * @param rfile RestartFile to read from
 */
RectangularBox::RectangularBox(RestartFile& rfile) {
    rfile.read(_center);
    rfile.read(_sides);
    rfile.read(_bitwidth);
    rfile.read(_maxside);

    LOGS("RectangularBox restarted");
}
