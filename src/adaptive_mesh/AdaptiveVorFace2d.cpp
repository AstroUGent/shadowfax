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
 * @file AdaptiveVorFace2d.cpp
 *
 * @brief 2D Voronoi face, used for geometrical face operations: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdaptiveVorFace2d.hpp"
#include "RestartFile.hpp"
#include "StateVector.hpp"
#include <cmath>
#include <iostream>
using namespace std;

/**
 * @brief Constructor
 *
 * @param pos Coordinates of the midpoint of the face
 */
AdaptiveVorFace2d::AdaptiveVorFace2d(double* pos) {
    _mid[0] = pos[0];
    _mid[1] = pos[1];
    _pos[0] = 0.;
    _pos[1] = 0.;
    _pos[2] = 0.;
    _pos[3] = 0.;
    _full = true;
}

/**
 * @brief Deactivate the face
 *
 * @deprecated No longer used
 */
void AdaptiveVorFace2d::deactivate() {
    _full = false;
}

/**
 * @brief Check whether this face is active
 *
 * @deprecated No longer used
 *
 * @return True if the face is active, false otherwise
 */
bool AdaptiveVorFace2d::active() {
    return _full;
}

/**
 * @brief Set the velocity of the face
 *
 * @param v New velocity for the face
 */
void AdaptiveVorFace2d::set_v(double* v) {
    _v[0] = v[0];
    _v[1] = v[1];
}

/**
 * @brief Calculate and set the left vertex of the face
 *
 * The left vertex is the midpoint of the circumsphere through three mesh
 * generating points, A, B, and C.
 *
 * @param a Coordinates of mesh generator A
 * @param b Coordinates of mesh generator B
 * @param c Coordinates of mesh generator C
 * @param vert Array to store the resulting vertex in
 * @return True if the vertex is valid, false otherwise
 */
bool AdaptiveVorFace2d::add_Lvertex(double* a, double* b, double* c,
                                    double* vert) {
    double rel1[2], rel2[2];
    rel1[0] = b[0] - a[0];
    rel1[1] = b[1] - a[1];
    rel2[0] = c[0] - a[0];
    rel2[1] = c[1] - a[1];
    double D, rel1_2, rel2_2, R[2];
    D = 2. * (rel1[0] * rel2[1] - rel1[1] * rel2[0]);
    rel1_2 = rel1[0] * rel1[0] + rel1[1] * rel1[1];
    rel2_2 = rel2[0] * rel2[0] + rel2[1] * rel2[1];
    R[0] = (rel2[1] * rel1_2 - rel1[1] * rel2_2) / D;
    R[1] = (rel1[0] * rel2_2 - rel2[0] * rel1_2) / D;
    if(R[0] != R[0] || R[1] != R[1]) {
        cerr << "NaN vertex!" << endl;
        cerr << a[0] << "\t" << a[1] << endl;
        cerr << b[0] << "\t" << b[1] << endl;
        cerr << c[0] << "\t" << c[1] << endl;
        vert[0] = (a[0] + b[0] + c[0]) / 3.;
        vert[1] = (a[1] + b[1] + c[1]) / 3.;
    } else {
        vert[0] = a[0] + R[0];
        vert[1] = a[1] + R[1];
    }
    _pos[0] = vert[0];
    _pos[1] = vert[1];
    return !(R[0] != R[0] || R[1] != R[1]);
}

/**
 * @brief Set the right vertex of the face
 *
 * @param vert New coordinates for the right vertex
 */
void AdaptiveVorFace2d::set_Rvertex(double* vert) {
    _pos[2] = vert[0];
    _pos[3] = vert[1];
}

/**
 * @brief Print face information to the given stream
 *
 * @param stream std::ostream to write to
 */
void AdaptiveVorFace2d::print(ostream& stream) {
    stream << _pos[0] << "\t" << _pos[1] << "\n";
    stream << _pos[2] << "\t" << _pos[3] << "\n\n";
}

/**
 * @brief Get the vertices of the face
 *
 * @param L Array to store the coordinates of the left vertex in
 * @param R Array to store the coordinates of the right vertex in
 */
void AdaptiveVorFace2d::get_points(double* L, double* R) {
    L[0] = _pos[0];
    L[1] = _pos[1];
    R[0] = _pos[2];
    R[1] = _pos[3];
}

/**
 * @brief Set the vertices of the face
 *
 * @param L Array with new coordinates for the left vertex
 * @param R Array with new coordinates for the right vertex
 */
void AdaptiveVorFace2d::set_points(double* L, double* R) {
    _pos[0] = L[0];
    _pos[1] = L[1];
    _pos[2] = R[0];
    _pos[3] = R[1];
}

// void AdaptiveVorFace2d::save_restart(ostream& stream){
//    stream.write((char*) _v, 2*sizeof(double));
//    stream.write((char*) _pos, 4*sizeof(double));
//    stream.write((char*) _mid, 2*sizeof(double));
//    char full = _full;
//    stream.write(&full, 1);
//}

// AdaptiveVorFace2d::AdaptiveVorFace2d(istream& stream){
//    stream.read((char*) _v, 2*sizeof(double));
//    stream.read((char*) _pos, 4*sizeof(double));
//    stream.read((char*) _mid, 2*sizeof(double));
//    char full;
//    stream.read(&full, 1);
//    _full = full > 0;
//}

/**
 * @brief Calculate the area of the face
 *
 * If the area is smaller than an arbitrary small floating point value (1.e-13),
 * the area is set to zero.
 *
 * @return The area of the face
 */
double AdaptiveVorFace2d::get_area() {
    double d[2];
    d[0] = _pos[0] - _pos[2];
    d[1] = _pos[1] - _pos[3];
    double area = sqrt(d[0] * d[0] + d[1] * d[1]);
    if(area < 1.e-13) {
        area = 0.;
    }
    return area;
}

/**
 * @brief Get the midpoint of the face
 *
 * @param midpoint Array to store the result in
 */
void AdaptiveVorFace2d::get_midpoint(double* midpoint) {
    midpoint[0] = 0.5 * (_pos[0] + _pos[2]);
    midpoint[1] = 0.5 * (_pos[1] + _pos[3]);
}

/**
 * @brief Rotate the velocity components of the given StateVector to a reference
 * frame where the x-axis coincides with the face normal
 *
 * @param W StateVector to transform
 * @param left Coordinates of the left neighbour of the face
 * @param right Coordinates of the right neighbour of the face
 */
void AdaptiveVorFace2d::transform(StateVector& W, double* left, double* right) {
    double d[3];
    d[0] = left[0] - right[0];
    d[1] = left[1] - right[1];
    d[2] = sqrt(d[0] * d[0] + d[1] * d[1]);
    // t is the angle between the face normal and the x-axis
    double cost = -d[0] / d[2];
    double sint = -d[1] / d[2];
    double W1 = W[1] * cost + W[2] * sint;
    W[2] = -W[1] * sint + W[2] * cost;
    W[1] = W1;
}

/**
 * @brief Rotate the velocity components of the given StateVector back to the
 * original reference frame
 *
 * @param W StateVector to transform
 * @param left Coordinates of the left neighbour of the face
 * @param right Coordinates of the right neighbour of the face
 */
void AdaptiveVorFace2d::invtransform(StateVector& W, double* left,
                                     double* right) {
    double d[3];
    d[0] = left[0] - right[0];
    d[1] = left[1] - right[1];
    d[2] = sqrt(d[0] * d[0] + d[1] * d[1]);
    double cost = -d[0] / d[2];
    double sint = -d[1] / d[2];
    double W1 = W[1] * cost - W[2] * sint;
    W[2] = W[1] * sint + W[2] * cost;
    W[1] = W1;
}

/**
 * @brief Dump the face to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void AdaptiveVorFace2d::dump(RestartFile& rfile) {
    rfile.write(_v, 2);
    rfile.write(_pos, 4);
    rfile.write(_mid, 2);
    rfile.write(_full);
}

/**
 * @brief Restart constructor. Initialize the face from the given RestartFile
 *
 * @param rfile RestartFile to read from
 */
AdaptiveVorFace2d::AdaptiveVorFace2d(RestartFile& rfile) {
    rfile.read(_v, 2);
    rfile.read(_pos, 4);
    rfile.read(_mid, 2);
    rfile.read(_full);
}
