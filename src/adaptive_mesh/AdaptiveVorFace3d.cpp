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
 * @file AdaptiveVorFace3d.cpp
 *
 * @brief 3D Voronoi face, used for geometrical face operations: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdaptiveVorFace3d.hpp"
#include "AdaptiveMeshUtils.hpp"  // for get_periodic_position
#include "Error.hpp"              // for my_exit
#include "RestartFile.hpp"        // for RestartFile
#include "StateVector.hpp"        // for StateVector
#include <cmath>                  // for sqrt, fabs
#include <iostream>               // for cerr, cout
using namespace std;

/**
 * @brief Restart constructor. Initialize the face from the given RestartFile
 *
 * @param rfile RestartFile to read from
 */
AdaptiveVorFace3d::AdaptiveVorFace3d(RestartFile& rfile) {
    rfile.read(_pos);
}

/**
 * @brief Dump the face to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void AdaptiveVorFace3d::dump(RestartFile& rfile) {
    rfile.write(_pos);
}

/**
 * @brief Constructor.
 */
AdaptiveVorFace3d::AdaptiveVorFace3d() {}

/**
 * @brief Add the given vertex to the internal list
 *
 * @param vert Coordinates of the vertex to add
 */
void AdaptiveVorFace3d::add_vertexpoint(double* vert) {
    _pos.push_back(vert[0]);
    _pos.push_back(vert[1]);
    _pos.push_back(vert[2]);
}

/**
 * @brief Print the face to the given stream
 *
 * @param stream std::ostream to write to
 */
void AdaptiveVorFace3d::print(ostream& stream) {
    if(!_pos.size()) {
        return;
    }
    for(unsigned int i = 0; i < _pos.size(); i += 3) {
        stream << _pos[i] << "\t" << _pos[i + 1] << "\t" << _pos[i + 2] << "\n";
    }
    stream << _pos[0] << "\t" << _pos[1] << "\t" << _pos[2] << "\n\n\n";
}

/**
 * @brief Print the requested periodic copy of the face to the given stream
 *
 * @param stream std::ostream to write to
 * @param wall Boundary flag indicating which periodic copy to print
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorFace3d::print_copy(ostream& stream, int wall, Cuboid& box) {
    if(!_pos.size()) {
        return;
    }
    double pos[3];
    for(unsigned int i = 0; i < _pos.size(); i += 3) {
        pos[0] = _pos[i];
        pos[1] = _pos[i + 1];
        pos[2] = _pos[i + 2];
        if(wall < 0) {
            AdaptiveMeshUtils::get_periodic_position(pos, wall, box);
        }
        stream << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\n";
    }
    pos[0] = _pos[0];
    pos[1] = _pos[1];
    pos[2] = _pos[2];
    if(wall < 0) {
        AdaptiveMeshUtils::get_periodic_position(pos, wall, box);
    }
    stream << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\n\n\n";
}

/**
 * @brief Get the vertices of the face
 *
 * @return std::vector containing the coordinates of the vertices of the face
 */
vector<double> AdaptiveVorFace3d::get_vertices() {
    return _pos;
}

/**
 * @brief Calculate the midpoint and area of the face
 *
 * The midpoint is the area-weighted average of the midpoints of triangles that
 * constitute the face. The same triangles are used to calculate the total area.
 *
 * The triangles are formed by taking the first vertex and combining it with all
 * possible other consecutive vertices.
 *
 * @param midpoint Array to store the coordinates of the midpoint in
 * @return The total area of the face
 */
double AdaptiveVorFace3d::calculate_quantities(double* midpoint) {
    double area = 0.;
    midpoint[0] = 0.;
    midpoint[1] = 0.;
    midpoint[2] = 0.;
    for(unsigned int i = 3; i < _pos.size() - 3; i += 3) {
        double ab[3] = {_pos[i] - _pos[i + 3], _pos[i + 1] - _pos[i + 4],
                        _pos[i + 2] - _pos[i + 5]};
        double ac[3] = {_pos[i] - _pos[0], _pos[i + 1] - _pos[1],
                        _pos[i + 2] - _pos[2]};
        double ab2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        double ac2 = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
        double abac = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];
        // in principle abac^2 can not be larger than ab2*ac2 (because it is
        // formally equal to ab2*ac2*cos(enclosed angle)^2) the problem is that
        // roundoff error can sometimes make the term between brackets negative
        // this will only happen for small values of tri_area, so no problem if
        // we artificially keep it positive
        // notice a difference of a factor 0.5 with the expression below.
        // we do not include this factor because we divide by the total area
        // in the end anyway...
        double tri_area = sqrt(fabs(ab2 * ac2 - abac * abac));
        area += tri_area;
        midpoint[0] += (_pos[i] + _pos[i + 3] + _pos[0]) * tri_area;
        midpoint[1] += (_pos[i + 1] + _pos[i + 4] + _pos[1]) * tri_area;
        midpoint[2] += (_pos[i + 2] + _pos[i + 5] + _pos[2]) * tri_area;
    }
    midpoint[0] /= 3. * area;
    midpoint[1] /= 3. * area;
    midpoint[2] /= 3. * area;
    return 0.5 * area;
}

/**
 * @brief Calculate the midpoint of the face
 *
 * The midpoint is the area-weighted average of the midpoints of triangles that
 * constitute the face.
 *
 * The triangles are formed by taking the first vertex and combining it with all
 * possible other consecutive vertices.
 *
 * @param midpoint Array to store the coordinates of the midpoint in
 */
void AdaptiveVorFace3d::get_midpoint(double* midpoint) {
    double area = 0.;
    midpoint[0] = 0.;
    midpoint[1] = 0.;
    midpoint[2] = 0.;
    if(_pos.size() < 9) {
        cerr << "hmmm, _pos.size() = " << _pos.size() << endl;
    }
    for(unsigned int i = 3; i < _pos.size() - 3; i += 3) {
        double ab[3] = {_pos[i] - _pos[i + 3], _pos[i + 1] - _pos[i + 4],
                        _pos[i + 2] - _pos[i + 5]};
        double ac[3] = {_pos[i] - _pos[0], _pos[i + 1] - _pos[1],
                        _pos[i + 2] - _pos[2]};
        double ab2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        double ac2 = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
        double abac = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];
        // in principle abac^2 can not be larger than ab2*ac2 (because it is
        // formally equal to ab2*ac2*cos(enclosed angle)^2) the problem is that
        // roundoff error can sometimes make the term between brackets negative
        // this will only happen for small values of tri_area, so no problem if
        // we artificially keep it positive
        // notice a difference of a factor 0.5 with the expression below.
        // we do not include this factor because we divide by the total area
        // in the end anyway...
        double tri_area = sqrt(fabs(ab2 * ac2 - abac * abac));
        area += tri_area;
        midpoint[0] += (_pos[i] + _pos[i + 3] + _pos[0]) * tri_area;
        midpoint[1] += (_pos[i + 1] + _pos[i + 4] + _pos[1]) * tri_area;
        midpoint[2] += (_pos[i + 2] + _pos[i + 5] + _pos[2]) * tri_area;
    }
    midpoint[0] /= 3. * area;
    midpoint[1] /= 3. * area;
    midpoint[2] /= 3. * area;
}

/**
 * @brief Calculate the area of the face
 *
 * The total area is the sum of the areas of small triangles that constitute the
 * face.
 *
 * The triangles are formed by taking the first vertex and combining it with all
 * possible other consecutive vertices.
 *
 * @return The total area of the face
 */
double AdaptiveVorFace3d::get_area() {
    double area = 0.;
    for(unsigned int i = 3; i < _pos.size() - 3; i += 3) {
        double ab[3] = {_pos[i] - _pos[i + 3], _pos[i + 1] - _pos[i + 4],
                        _pos[i + 2] - _pos[i + 5]};
        double ac[3] = {_pos[i] - _pos[0], _pos[i + 1] - _pos[1],
                        _pos[i + 2] - _pos[2]};
        double ab2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        double ac2 = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
        double abac = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];
        // in principle abac^2 can not be larger than ab2*ac2 (because it is
        // formally equal to ab2*ac2*cos(enclosed angle)^2) the problem is that
        // roundoff error can sometimes make the term between brackets negative
        // this will only happen for small values of tri_area, so no problem if
        // we artificially keep it positive
        double tri_area = 0.5 * sqrt(fabs(ab2 * ac2 - abac * abac));
        area += tri_area;
    }
    if(area != area) {
        area = 0.;
    }
    return area;
}

/**
 * @brief Rotate the velocity components of the given StateVector to a reference
 * frame where the x-axis coincides with the face normal
 *
 * The y- and z-axis are chosen as arbitrary orthogonal directions in the plane
 * of the face.
 *
 * @param W StateVector to transform
 * @param left Coordinates of the left neighbour of the face
 * @param right Coordinates of the right neighbour of the face
 */
void AdaptiveVorFace3d::transform(StateVector& W, double* left, double* right) {
    double d[4];
    d[0] = left[0] - right[0];
    d[1] = left[1] - right[1];
    d[2] = left[2] - right[2];
    d[3] = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    // p is the angle between the projection of the face normal
    // in the xy-plane and the x-axis
    // t is the angle between the face normal and the z-axis
    // see https://en.wikipedia.org/wiki/Euler_angles:
    // p is the alpha Euler angle
    // t is the beta Euler angle
    // the gamma Euler angle does not have to be specified
    // (because we don't care about the y- and z-axis, as
    //  long as the x-axis coincides with the face normal)
    double cost = -d[2] / d[3];
    double sint = sqrt(1. - cost * cost);
    if(sint != sint) {
        cout << "it happens" << endl;
        my_exit();
    }
    double sinp, cosp;
    if(sint) {
        sinp = -d[1] / d[3] / sint;
        cosp = -d[0] / d[3] / sint;
    } else {
        sinp = 0.;
        cosp = 1.;
    }
    // coordinate transformation
    // this is conceptually done in two stages:
    //  * rotate around z-axis over p, x- and y-coordinate change
    //  * rotate around y-axis over t, x- and z-coordinate change
    double W1acc = W[1] * cosp + W[2] * sinp;
    double W1 = W1acc * sint + W[3] * cost;
    double W2 = -W[1] * sinp + W[2] * cosp;
    W[3] = -W1acc * cost + W[3] * sint;
    W[1] = W1;
    W[2] = W2;
}

/**
 * @brief Rotate the velocity components of the given StateVector back to the
 * original reference frame
 *
 * @param W StateVector to transform
 * @param left Coordinates of the left neighbour of the face
 * @param right Coordinates of the right neighbour of the face
 */
void AdaptiveVorFace3d::invtransform(StateVector& W, double* left,
                                     double* right) {
    double d[4];
    d[0] = left[0] - right[0];
    d[1] = left[1] - right[1];
    d[2] = left[2] - right[2];
    d[3] = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    double cost = -d[2] / d[3];
    double sint = sqrt(1. - cost * cost);
    double sinp, cosp;
    if(sint) {
        sinp = -d[1] / d[3] / sint;
        cosp = -d[0] / d[3] / sint;
    } else {
        sinp = 0.;
        cosp = 1.;
    }
    double W1acc = W[1] * sint - W[3] * cost;
    double W1 = W1acc * cosp - W[2] * sinp;
    double W2 = W1acc * sinp + W[2] * cosp;
    W[3] = W[1] * cost + W[3] * sint;
    W[1] = W1;
    W[2] = W2;
}
