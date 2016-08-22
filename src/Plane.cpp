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
 * @file Plane.cpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief 3D plane: implementation
 */
#if ndim_ == 3

#include "Plane.hpp"
#include "Line.hpp"    // for Line
#include "VorGen.hpp"  // for VorGen
#include <cmath>       // for fabs, sqrt
using namespace std;

/**
 * @brief Construct a plane through the given point coordinates and with the
 * given normal
 *
 * The normal does not have to be normalized, this is done in this method.
 *
 * @param nx x-component of the normal vector of the plane
 * @param ny y-component of the normal vector of the plane
 * @param nz z-component of the normal vector of the plane
 * @param x x-coodinate of a reference point on the plane
 * @param y y-coodinate of a reference point on the plane
 * @param z z-coodinate of a reference point on the plane
 */
Plane::Plane(double nx, double ny, double nz, double x, double y, double z) {
    double norm = sqrt(nx * nx + ny * ny + nz * nz);
    _nx = nx / norm;
    _ny = ny / norm;
    _nz = nz / norm;
    _x = x;
    _y = y;
    _z = z;
}

/**
 * @brief Construct a plane throught the given points
 *
 * @param point1 First VorGen
 * @param point2 Second VorGen
 * @param point3 Third VorGen
 */
Plane::Plane(VorGen* point1, VorGen* point2, VorGen* point3) {
    _x = point1->x();
    _y = point1->y();
    _z = point1->z();
    double a[3], b[3];
    a[0] = point2->x() - point1->x();
    a[1] = point2->y() - point1->y();
    a[2] = point2->z() - point1->z();
    b[0] = point3->x() - point1->x();
    b[1] = point3->y() - point1->y();
    b[2] = point3->z() - point1->z();
    _nx = a[1] * b[2] - a[2] * b[1];
    _ny = a[2] * b[0] - a[0] * b[2];
    _nz = a[0] * b[1] - a[1] * b[0];
    double norm = sqrt(_nx * _nx + _ny * _ny + _nz * _nz);
    _nx = _nx / norm;
    _ny = _ny / norm;
    _nz = _nz / norm;
}

/**
 * @brief Calculate the intersection of the Plane with the given Line
 *
 * @param line Line to intersect with
 * @return VorGen that is the intersection of the plane and the line
 */
VorGen Plane::intersect(Line* line) {
    double num, denom, d;
    num = (_x - line->x()) * _nx + (_y - line->y()) * _ny +
          (_z - line->z()) * _nz;
    denom = line->xdir() * _nx + line->ydir() * _ny + line->zdir() * _nz;
    d = num / denom;
    return VorGen(line->x() + d * line->xdir(), line->y() + d * line->ydir(),
                  line->z() + d * line->zdir());
}

/**
 * @brief Test whether the given Plane intersects this Plane
 *
 * @deprecated Inexact function that is not used for the moment
 *
 * @param plane Second Plane
 * @return True if the planes intersect, false otherwise
 */
bool Plane::intersect(Plane* plane) {
    double a[3];
    a[0] = _ny * plane->_nz - _nz * plane->_ny;
    a[1] = _nz * plane->_nx - _nx * plane->_nz;
    a[2] = _nx * plane->_ny - _ny * plane->_nx;
    // if a is the null vector, the planes coincide and there is no intersection
    bool result = true;
    if(!(a[0] == 0 && a[1] == 0 && a[2] == 0)) {
        result = false;
    }
    return result;
}

/**
 * @brief Determine the intersection point of this Plane and the 2 given planes
 *
 * @param plane1 Second Plane
 * @param plane2 Third plane
 * @return VorGen that is the intersection of the three planes
 */
VorGen Plane::intersect_planes(Plane* plane1, Plane* plane2) {
    double dot[3], num[3], denom;
    dot[0] = _x * _nx + _y * _ny + _z * _nz;
    dot[1] = plane1->_x * plane1->_nx + plane1->_y * plane1->_ny +
             plane1->_z * plane1->_nz;
    dot[2] = plane2->_x * plane2->_nx + plane2->_y * plane2->_ny +
             plane2->_z * plane2->_nz;
    num[0] = dot[0] * (plane1->_ny * plane2->_nz - plane1->_nz * plane2->_ny) +
             dot[1] * (plane2->_ny * _nz - plane2->_nz * _ny) +
             dot[2] * (_ny * plane1->_nz - _nz * plane1->_ny);
    num[1] = dot[0] * (plane1->_nz * plane2->_nx - plane1->_nx * plane2->_nz) +
             dot[1] * (plane2->_nz * _nx - plane2->_nx * _nz) +
             dot[2] * (_nz * plane1->_nx - _nx * plane1->_nz);
    num[2] = dot[0] * (plane1->_nx * plane2->_ny - plane1->_ny * plane2->_nx) +
             dot[1] * (plane2->_nx * _ny - plane2->_ny * _nx) +
             dot[2] * (_nx * plane1->_ny - _ny * plane1->_nx);
    denom = _nx * plane1->_ny * plane2->_nz + _ny * plane1->_nz * plane2->_nx +
            _nz * plane1->_nx * plane2->_ny - _nz * plane1->_ny * plane2->_nx -
            _ny * plane1->_nx * plane2->_nz - _nx * plane1->_nz * plane2->_ny;
    return VorGen(num[0] / denom, num[1] / denom, num[2] / denom);
}

/**
 * @brief Calculate the distance between this Plane and the given VorGen
 *
 * @param point VorGen to compute the distance to
 * @return Distance between the point and the plane
 */
double Plane::distance(VorGen* point) {
    return fabs(_nx * (point->x() - _x) + _ny * (point->y() - _y) +
                _nz * (point->z() - _z));
}

/**
 * @brief Calculate the distance between two planes
 *
 * Only makes sense when planes are parallel (which can be checked by the
 * function Plane::intersect(Plane*).
 *
 * @deprecated Function is not used for the moment
 *
 * @param plane Second plane
 * @return Distance between this plane and the second plane
 */
double Plane::distance(Plane* plane) {
    return fabs(_nx * (plane->_x - _x) + _ny * (plane->_y - _y) +
                _nz * (plane->_z - _z));
}

#endif
