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
 * @file Line.cpp
 *
 * @brief Geometrical line (segment): implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Line.hpp"
#include "VorGen.hpp"
#include "Plane.hpp"
#include "ExArith.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ostream>
#include <valarray>
#include <cmath>
using namespace std;

/**
 * @brief Constructor
 *
 * Construct the line segment through two points.
 *
 * @param point1 VorGen representing the first point
 * @param point2 VorGen representing the second point
 */
Line::Line(VorGen* point1, VorGen* point2){
    Vec p1 = point1->get_p12();
    Vec p2 = point2->get_p12();
    _x1 = p1.x();
    _x2 = p2.x();
    _y1 = p1.y();
    _y2 = p2.y();
    _xdir = _x2 - _x1;
    _ydir = _y2 - _y1;
#if ndim_==3
    _z1 = p1.z();
    _z2 = p2.z();
    _zdir = _z2 - _z1;
    double dirnorm = sqrt(_xdir*_xdir + _ydir*_ydir + _zdir*_zdir);
    _zdir /= dirnorm;
#else
    double dirnorm = sqrt(_xdir*_xdir + _ydir*_ydir);
#endif
    _xdir /= dirnorm;
    _ydir /= dirnorm;
}

#if ndim_==3
/**
 * @brief Constructor
 *
 * Construct the line through the given point and with the given direction. The
 * given direction does not have to be normalized, this is done in this
 * function.
 *
 * @param x x-coordinate of the point
 * @param y y-coordinate of the point
 * @param z z-coordinate of the point
 * @param xdir x-coordinate of the direction vector of the line
 * @param ydir y-coordinate of the direction vector of the line
 * @param zdir z-coordinate of the direction vector of the line
 */
Line::Line(double x, double y, double z, double xdir, double ydir, double zdir){
    _x1 = x;
    _y1 = y;
    _z1 = z;
    _xdir = xdir;
    _ydir = ydir;
    _zdir = zdir;
    double dirnorm = sqrt(_xdir*_xdir + _ydir*_ydir + _zdir*_zdir);
    _xdir /= dirnorm;
    _ydir /= dirnorm;
    _zdir /= dirnorm;
    _x2 = x + xdir;
    _y2 = y + ydir;
    _z2 = z + zdir;
}
#else
/**
 * @brief Constructor
 *
 * Construct the line through the given point and with the given direction. The
 * given direction does not have to be normalized, this is done in this
 * function.
 *
 * @param x x-coordinate of the point
 * @param y y-coordinate of the point
 * @param xdir x-coordinate of the direction vector of the line
 * @param ydir y-coordinate of the direction vector of the line
 */
Line::Line(double x, double y, double xdir, double ydir){
    _x1 = x;
    _y1 = y;
    _xdir = xdir;
    _ydir = ydir;
    double dirnorm = sqrt(_xdir*_xdir + _ydir*_ydir);
    _xdir /= dirnorm;
    _ydir /= dirnorm;
    _x2 = x + xdir;
    _y2 = y + ydir;
}
#endif

/**
 * @brief Get the midpoint of the line segment
 *
 * @return The coordinates of the midpoint of the line segment
 */
vector<double> Line::get_midpoint(){
    vector<double> x(ndim_);
    x[0] = (_x1+_x2)*0.5;
    x[1] = (_y1+_y2)*0.5;
#if ndim_==3
    x[2] = (_z1+_z2)*0.5;
#endif
    return x;
}

#if ndim_==3
/**
 * @brief Get the orthogonal bisector of the line segment
 *
 * @return Plane that is the orthogonal bisector of the line segment
 */
Plane Line::get_bisector(){
    vector<double> pos(3);
    pos = this->get_midpoint();
    return Plane(_xdir, _ydir, _zdir, pos[0], pos[1], pos[2]);
}
#endif

/**
 * @brief Print the line to the given stream
 *
 * @param stream std::ostream to write to
 */
void Line::print(ostream &stream){
    stream << "(" << _x1 << "," << _y1;
#if ndim_==3
    stream << "," << _z1;
#endif
    stream << ")\t(" << _x2 << "," << _y2;
#if ndim_==3
    stream << "," << _z2;
#endif
    stream << ")";
}

/**
 * @brief Get the x-coordinate of the first endpoint
 *
 * @return x-coordinate of the first endpoint
 */
double Line::x(){
    return _x1;
}

/**
 * @brief Get the x-coordinate of the direction vector of the line
 *
 * @return x-coordinate of the direction vector
 */
double Line::xdir(){
    return _xdir;
}

/**
 * @brief Get the y-coordinate of the first endpoint
 *
 * @return y-coordinate of the first endpoint
 */
double Line::y(){
    return _y1;
}

/**
 * @brief Get the y-coordinate of the direction vector of the line
 *
 * @return y-coordinate of the direction vector
 */
double Line::ydir(){
    return _ydir;
}

#if ndim_==3
/**
 * @brief Get the z-coordinate of the first endpoint
 *
 * @return z-coordinate of the first endpoint
 */
double Line::z(){
    return _z1;
}

/**
 * @brief Get the z-coordinate of the direction vector of the line
 *
 * @return z-coordinate of the direction vector
 */
double Line::zdir(){
    return _zdir;
}
#endif

/**
 * @brief Get the coordinates of the point on the line with given x-, y- or
 * z-coordinate
 *
 * @deprecated Not used
 *
 * @param coord char representing the coordinate ('x', 'y' or 'z')
 * @param value Value of the given coordinate
 * @return All coordinates of a point on the line with the given coordinate
 */
vector<double> Line::get_coordinates(char coord, double value){
    vector<double> result;
    double factor;
    if(coord == 'x'){
        factor = (value - _x1)/_xdir;
        result.push_back(value);
        result.push_back(_y1 + factor*_ydir);
#if ndim_==3
        result.push_back(_z1 + factor*_zdir);
#endif
    } else if(coord == 'y'){
        factor = (value - _y1)/_ydir;
        result.push_back(_x1 + factor*_xdir);
        result.push_back(value);
#if ndim_==3
        result.push_back(_z1 + factor*_zdir);
    } else {
        factor = (value - _z1)/_zdir;
        result.push_back(_x1 + factor*_xdir);
        result.push_back(_y1 + factor*_ydir);
        result.push_back(value);
#endif
    }
    return result;
}
