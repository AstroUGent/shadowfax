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
 * @file VorFace.cpp
 *
 * @brief Voronoi face: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "VorFace.hpp"
#include "Error.hpp"
#include "MPIGlobal.hpp"
#include "TimeLine.hpp"
#include "VorCell.hpp"
#include "VorGen.hpp"
#include "riemann/RiemannSolver.hpp"
#include "utilities/GasParticle.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
using namespace std;

/**
 * @brief Constructor
 *
 * Create a Voronoi face with the VorGen with the given index in the given list
 * as left generator.
 *
 * @param index Index of the generator in the DelTess VorGen list
 * @param points Reference to the DelTess VorGen list
 */
VorFace::VorFace(unsigned int index, vector<VorGen*>& points) {
    _left = points[index];
    _right = NULL;
    _vleft = points[index]->get_particle_id();
    _vright = 0;
    _area = 0.;
}

/**
 * @brief Write the face in plottable format to the given stream
 *
 * @param stream std::ostream to write to
 */
void VorFace::print(std::ostream& stream) {
    for(unsigned int i = 0; i < _vertices.size(); i++) {
        _vertices[i]->print(stream);
        stream << "\n";
    }
    _vertices[0]->print(stream);
    stream << "\n";
    stream << "\n";
    stream << "\n";
}

/**
 * @brief Write the face in POVRAY plottable format to the given stream
 *
 * @param stream std::ostream to write to
 */
void VorFace::print_pov_frame(ostream& stream) {
    Vec last = _vertices.back()->get_position();
    for(unsigned int i = 0; i < _vertices.size(); i++) {
        Vec current = _vertices[i]->get_position();
        if((current - last).norm() < 1.e-5) {
            continue;
        }
        stream << "sphere{<" << current[0] << "," << current[1] << ","
               << current[2] << ">,r\n texture { Frame } }" << endl;
        stream << "cylinder{<" << last[0] << "," << last[1] << "," << last[2]
               << ">,<" << current[0] << "," << current[1] << "," << current[2]
               << ">,r\n texture { Frame } }" << endl;
        last = current;
    }
}

/**
 * @brief Write the face in POVRAY plottable format to the given stream
 *
 * @param stream std::ostream to write to
 */
void VorFace::print_pov(ostream& stream) {
    Vec first = _vertices[0]->get_position();
    for(unsigned int i = 1; i < _vertices.size() - 1; i++) {
        Vec a = _vertices[i]->get_position();
        Vec b = _vertices[i + 1]->get_position();
        if((first - a).norm() < 1.e-9 || (first - b).norm() < 1.e-9 ||
           (a - b).norm() < 1.e-9) {
            continue;
        }
        stream << "triangle {\n <" << first[0] << ", " << first[1] << ", "
               << first[2] << ">, <" << a[0] << ", " << a[1] << ", " << a[2]
               << ">, <" << b[0] << ", " << b[1] << ", " << b[2] << ">" << endl;
        stream << " texture { Cell }\n }" << endl;
    }
}

/**
 * @brief Print the face in Leaflet plottable format to the given stream
 *
 * @param vstream std::ostream to write to
 * @param ox x-coordinate of the origin of the reference frame
 * @param oy y-coordinate of the origin of the reference frame
 * @param center Pointer to the generator of the cell currently being plotted
 */
void VorFace::print_leaflet(ostream& vstream, int ox, int oy, VorGen* center) {
    short coords[2];
    if(center == _left) {
        coords[0] =
                (short)(((int)(_vertices[0]->get_position().x() * 16384)) - ox);
        coords[1] =
                (short)(256 -
                        (((int)((_vertices[0]->get_position().y()) * 16384)) -
                         oy));
    } else {
        coords[0] =
                (short)(((int)(_vertices[1]->get_position().x() * 16384)) - ox);
        coords[1] =
                (short)(256 -
                        (((int)((_vertices[1]->get_position().y()) * 16384)) -
                         oy));
    }
    vstream.write((char*)&coords[0], 2 * sizeof(short));
}

/**
 * @brief Check if the face overlaps with the given geometrical box
 *
 * @warning Only works in 2D!
 *
 * @param box Origin and sides of a 2D box
 * @return True if the face overlaps with the box, false otherwise
 */
bool VorFace::overlap(double* box) {
    for(unsigned int i = 0; i < 2; i++) {
        Vec pos = _vertices[i]->get_position();
        if(pos.x() >= box[0] && pos.x() <= box[2] && pos.y() >= box[1] &&
           pos.y() <= box[3]) {
            return true;
        }
    }
    Vec pos = _vertices[0]->get_position();
    if(pos.x() >= box[0] && pos.x() <= box[2]) {
        pos = _vertices[1]->get_position();
        if(pos.y() >= box[1] && pos.y() <= box[3]) {
            return true;
        }
    } else {
        if(pos.y() >= box[1] && pos.y() <= box[3]) {
            pos = _vertices[1]->get_position();
            if(pos.x() >= box[0] && pos.x() <= box[2]) {
                return true;
            }
        }
    }
    return false;
}

/**
 * @brief Add the given vertex to the face
 *
 * @param point VorGen to be added to the end of the internal vertex list
 */
void VorFace::add_vertex(VorGen* point) {
    _vertices.push_back(point);
}

/**
 * @brief Get the vertices of the face
 *
 * @return Vertices of the face
 */
vector<VorGen*> VorFace::get_vertices() {
    return _vertices;
}

/**
 * @brief Add the given VorGen as face neighbour to the face
 *
 * This is only done to be able to construct the evolving mesh if the mesh
 * evolution algorithm is used.
 *
 * @param ngb VorGen to be added to the end of the face neighbour list
 */
void VorFace::add_facengb(VorGen* ngb) {
    _ngbs.push_back(ngb);
}

/**
 * @brief Get the face neighbours of the face
 * @return Face neighbours of the face
 */
vector<VorGen*> VorFace::get_facengbs() {
    return _ngbs;
}

/**
 * @brief Make the given VorGen the right generator of this face
 *
 * @param ngb Right generator of the face
 */
void VorFace::add_ngb(VorGen* ngb) {
    _right = ngb;
}

/**
 * @brief Get the right generator of the face
 *
 * @return Right generator of the face
 */
VorGen* VorFace::get_ngb() {
    return _right;
}

/**
 * @brief Get the left generator of the face
 *
 * @return Left generator of the face
 */
VorGen* VorFace::get_left() {
    return _left;
}

/**
 * @brief Set the index of the right generator of the face in the DelTess VorGen
 * list
 *
 * @param ngb Index of the right generator of the face
 */
void VorFace::add_ngb_id(unsigned int ngb) {
    _vright = ngb;
}

/**
 * @brief Get the index of the right generator of the face in the DelTess VorGen
 * list
 *
 * @return Index of the right generator
 */
unsigned int VorFace::get_ngb_id() {
    return _vright;
}

/**
 * @brief Get the geometrical midpoint of the face
 *
 * @return The midpoint of the face
 */
Vec& VorFace::get_midpoint() {
    return _midpoint;
}

/**
 * @brief Set the geometrical midpoint of the face to the given value
 *
 * @param midpoint Coordinates for the midpoint of the face
 */
void VorFace::set_midpoint(double* midpoint) {
    for(unsigned int i = ndim_; i--;) {
        _midpoint[i] = midpoint[i];
    }
}

/**
 * @brief Calculate the geometrical midpoint of the face
 *
 * For a 2D face, this is simply the midpoint of the segment between the two
 * vertices of the face.
 *
 * For a 3D face, the midpoint is the centroid of the polygon constituted by the
 * vertices. The midpoint is then the area weighted mean of the centroids of
 * small triangles consisting of the first vertex of the face and two
 * consecutive other vertices in the list.
 */
void VorFace::calculate_midpoint() {
#if ndim_ == 3
    double area = 0.;
    _midpoint[0] = 0.;
    _midpoint[1] = 0.;
    _midpoint[2] = 0.;
    for(unsigned int i = 1; i < _vertices.size() - 1; i++) {
        double ab[3] = {_vertices[i]->x() - _vertices[i + 1]->x(),
                        _vertices[i]->y() - _vertices[i + 1]->y(),
                        _vertices[i]->z() - _vertices[i + 1]->z()};
        double ac[3] = {_vertices[i]->x() - _vertices[0]->x(),
                        _vertices[i]->y() - _vertices[0]->y(),
                        _vertices[i]->z() - _vertices[0]->z()};
        double ab2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        double ac2 = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
        double abac = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];
        // in principle abac^2 can not be larger than ab2*ac2 (because it is
        // formally equal to ab2*ac2*cos(enclosed angle)^2)
        // the problem is that roundoff error can sometimes make the term
        // between brackets negative
        // this will only happen for small values of tri_area, so no problem if
        // we artificially keep it positive
        double tri_area = 0.5 * sqrt(fabs(ab2 * ac2 - abac * abac));
        area += tri_area;
        _midpoint[0] += (_vertices[i]->x() + _vertices[i + 1]->x() +
                         _vertices[0]->x()) *
                        tri_area;
        _midpoint[1] += (_vertices[i]->y() + _vertices[i + 1]->y() +
                         _vertices[0]->y()) *
                        tri_area;
        _midpoint[2] += (_vertices[i]->z() + _vertices[i + 1]->z() +
                         _vertices[0]->z()) *
                        tri_area;
    }
    _midpoint[0] /= 3. * area;
    _midpoint[1] /= 3. * area;
    _midpoint[2] /= 3. * area;
    if(_midpoint[0] != _midpoint[0]) {
        // do nothing, this problem is tackled elsewhere by setting the area
        // to 0 and not using faces with zero area
    }
#else
    _midpoint[0] = 0.5 * (_vertices[0]->x() + _vertices[1]->x());
    _midpoint[1] = 0.5 * (_vertices[0]->y() + _vertices[1]->y());
#endif
}

/**
 * @brief Get the geometrical area of the face
 *
 * @return The area of the face
 */
double VorFace::get_area() {
    return _area;
}

/**
 * @brief Set the geometrical area of the face to the given value
 *
 * @param area New value for the area of the face
 */
void VorFace::set_area(double area) {
    _area = area;
}

/**
 * @brief Calculate the geometrical area of the face
 *
 * For 2D, this is just the length of the segment between the two vertices of
 * the face.
 *
 * For 3D, the area is the area of the polygon constituted by the vertices of
 * the face. We calculate it as the sum of the areas of the small triangles
 * consisting of the first vertex and two consecutive other vertices in the
 * list.
 */
void VorFace::calculate_area() {
#if ndim_ == 3
    double area = 0.;
    // temporary fix: if the midpoint is NaN, then return 0
    if(_midpoint[0] != _midpoint[0]) {
        _area = 0.;
    } else {
        for(unsigned int i = 0; i < _vertices.size(); i++) {
            // calculate area triangle _vertices[i], _vertices[i+1], _midpoint
            double ab[3] = {_vertices[i]->x() -
                                    _vertices[(i + 1) % _vertices.size()]->x(),
                            _vertices[i]->y() -
                                    _vertices[(i + 1) % _vertices.size()]->y(),
                            _vertices[i]->z() -
                                    _vertices[(i + 1) % _vertices.size()]->z()};
            double ac[3] = {_vertices[i]->x() - _midpoint[0],
                            _vertices[i]->y() - _midpoint[1],
                            _vertices[i]->z() - _midpoint[2]};
            double ab2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
            double ac2 = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
            double abac = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];
            area += 0.5 * sqrt(fabs(ab2 * ac2 - abac * abac));
        }
        if(area != area) {
            area = 0.;
        }
        _area = area;
    }
#else
    double d[2];
    d[0] = _vertices[0]->x() - _vertices[1]->x();
    d[1] = _vertices[0]->y() - _vertices[1]->y();
    _area = sqrt(d[0] * d[0] + d[1] * d[1]);
#endif
}

/**
 * @brief Rotate the velocity components of the given StateVector to a reference
 * frame where the x-axis coincides with the segment between the given left and
 * right Vec
 *
 * The y- and z-axis are chosen as arbitrary orthogonal directions in the plane
 * of the face generated by the left and right Vec.
 *
 * @param W StateVector to transform
 * @param left Position of the left generator of the face
 * @param right Position of the right generator of the face
 */
void VorFace::transform(StateVector& W, const Vec& left, const Vec& right) {
#if ndim_ == 3
    // vector n is the new x-axis and is just the normalized vector pointing
    // from left to right
    double n[4];
    n[0] = right.x() - left.x();
    n[1] = right.y() - left.y();
    n[2] = right.z() - left.z();
    n[3] = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= n[3];
    n[1] /= n[3];
    n[2] /= n[3];
    // vector o is an orthonormal vector to n
    // we construct it as (-ny, nx, 0), except when n coincides with the z-axis
    double o[4];
    double p[3];
    if(n[0] || n[1]) {
        o[0] = -n[1];
        o[1] = n[0];
        o[2] = 0.;
        o[3] = sqrt(o[0] * o[0] + o[1] * o[1]);
        o[0] /= o[3];
        o[1] /= o[3];
        // p is orthonormal to both n and o. It is easy to check that this is
        // indeed the case.
        p[0] = -n[2] * o[1];
        p[1] = n[2] * o[0];
        p[2] = o[3];
    } else {
        o[0] = 0.;
        o[1] = 1.;
        o[2] = 0.;
        o[3] = 1.;
        // in this case, p is the negative x-axis
        // we could also use the positive x-axis, but then we switch the
        // orientation of the axes
        // which is of course fine, but still feels somewhat unethical :p
        p[0] = -1.;
        p[1] = 0.;
        p[2] = 0.;
    }
    Vec v(W[1], W[2], W[3]);
    W[1] = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
    W[2] = v[0] * o[0] + v[1] * o[1];
    W[3] = v[0] * p[0] + v[1] * p[1] + v[2] * p[2];
#else
    double d[3];
    d[0] = left.x() - right.x();
    d[1] = left.y() - right.y();
    d[2] = sqrt(d[0] * d[0] + d[1] * d[1]);
    // t is the angle between the face normal and the x-axis
    double cost = -d[0] / d[2];
    double sint = -d[1] / d[2];
    double W1 = W[1] * cost + W[2] * sint;
    W[2] = -W[1] * sint + W[2] * cost;
    W[1] = W1;
#endif
}

/**
 * @brief The inverse of VorFace::transform()
 *
 * Rotate the velocity components of the given StateVector back to the static
 * simulation reference frame.
 *
 * @param W StateVector to transform
 * @param left Position of the left generator of the face
 * @param right Position of the right generator of the face
 */
void VorFace::invtransform(StateVector& W, const Vec& left, const Vec& right) {
#if ndim_ == 3
    // vector n is the new x-axis and is just the normalized vector pointing
    // from left to right
    double n[4];
    n[0] = right.x() - left.x();
    n[1] = right.y() - left.y();
    n[2] = right.z() - left.z();
    n[3] = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= n[3];
    n[1] /= n[3];
    n[2] /= n[3];
    // vector o is an orthonormal vector to n
    // we construct it as (-ny, nx, 0), except when n coincides with the z-axis
    double o[4];
    double p[3];
    if(n[0] || n[1]) {
        o[0] = -n[1];
        o[1] = n[0];
        o[2] = 0.;
        o[3] = sqrt(o[0] * o[0] + o[1] * o[1]);
        o[0] /= o[3];
        o[1] /= o[3];
        // p is orthonormal to both n and o. It is easy to check that this is
        // indeed the case.
        p[0] = -n[2] * o[1];
        p[1] = n[2] * o[0];
        p[2] = o[3];
    } else {
        o[0] = 0.;
        o[1] = 1.;
        o[2] = 0.;
        o[3] = 1.;
        // in this case, p is the negative x-axis
        // we could also use the positive x-axis, but then we switch the
        // orientation of the axes
        // which is of course fine, but still feels somewhat unethical :p
        p[0] = -1.;
        p[1] = 0.;
        p[2] = 0.;
    }
    Vec v(W[1], W[2], W[3]);
    W[1] = v[0] * n[0] + v[1] * o[0] + v[2] * p[0];
    W[2] = v[0] * n[1] + v[1] * o[1] + v[2] * p[1];
    W[3] = v[0] * n[2] + v[2] * p[2];
#else
    double d[3];
    d[0] = left.x() - right.x();
    d[1] = left.y() - right.y();
    d[2] = sqrt(d[0] * d[0] + d[1] * d[1]);
    double cost = -d[0] / d[2];
    double sint = -d[1] / d[2];
    double W1 = W[1] * cost - W[2] * sint;
    W[2] = W[1] * sint + W[2] * cost;
    W[1] = W1;
#endif
}

/**
 * @brief Get the components of the normal to the face
 *
 * If t is the angle between the normal and the z-axis and p the angle between
 * the projection of the normal in the xy-plane and the x-axis, then the
 * components are \f$(\sin(t)\cos(p), \sin(t)\sin(p), \cos(t))\f$.
 *
 * @param angles Array to store the result in
 */
void VorFace::get_normal(double* angles) {
#if ndim_ == 3
    double d[4];
    d[0] = _left->x() - _right->x();
    d[1] = _left->y() - _right->y();
    d[2] = _left->z() - _right->z();
    d[3] = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    angles[0] = -d[0] / d[3];
    angles[1] = -d[1] / d[3];
    angles[2] = -d[2] / d[3];
#else
    double d[3];
    d[0] = _left->x() - _right->x();
    d[1] = _left->y() - _right->y();
    d[2] = sqrt(d[0] * d[0] + d[1] * d[1]);
    angles[0] = -d[0] / d[2];
    angles[1] = -d[1] / d[2];
#endif
}

/**
 * @brief Get the velocity of the face
 *
 * @return Velocity of the face
 */
Vec& VorFace::get_v() {
    return _v;
}

/**
 * @brief Set the velocity of the face based on the velocities of the left and
 * right generator
 *
 * @param particles Obsolete parameter, no idea why it is here
 */
void VorFace::set_v(ParticleVector& particles) {
#ifdef STATIC
// do nothing
#else
#if ndim_ == 3
    if(_right->get_particle() != NULL) {
        double rRL[4];
        rRL[0] = _right->x() - _left->x();
        rRL[1] = _right->y() - _left->y();
        rRL[2] = _right->z() - _left->z();
        rRL[3] = rRL[0] * rRL[0] + rRL[1] * rRL[1] + rRL[2] * rRL[2];
        double fac =
                (_left->get_particle()->vx() - _right->get_particle()->vx()) *
                (_midpoint[0] - 0.5 * (_left->x() + _right->x()));
        fac += (_left->get_particle()->vy() - _right->get_particle()->vy()) *
               (_midpoint[1] - 0.5 * (_left->y() + _right->y()));
        fac += (_left->get_particle()->vz() - _right->get_particle()->vz()) *
               (_midpoint[2] - 0.5 * (_left->z() + _right->z()));
        fac /= rRL[3];
        _v[0] = 0.5 * (_left->get_particle()->vx() +
                       _right->get_particle()->vx()) +
                fac * rRL[0];
        _v[1] = 0.5 * (_left->get_particle()->vy() +
                       _right->get_particle()->vy()) +
                fac * rRL[1];
        _v[2] = 0.5 * (_left->get_particle()->vz() +
                       _right->get_particle()->vz()) +
                fac * rRL[2];
    } else {
        _v[0] = 0.;
        _v[1] = 0.;
        _v[2] = 0.;
    }
#else
    if(_right->get_particle() != NULL) {
        double rRL[3];
        rRL[0] = _right->x() - _left->x();
        rRL[1] = _right->y() - _left->y();
        rRL[2] = rRL[0] * rRL[0] + rRL[1] * rRL[1];
        double fac =
                (_left->get_particle()->vx() - _right->get_particle()->vx()) *
                (_midpoint[0] - 0.5 * (_left->x() + _right->x()));
        fac += (_left->get_particle()->vy() - _right->get_particle()->vy()) *
               (_midpoint[1] - 0.5 * (_left->y() + _right->y()));
        fac /= rRL[2];
        _v[0] = 0.5 * (_left->get_particle()->vx() +
                       _right->get_particle()->vx()) +
                fac * rRL[0];
        _v[1] = 0.5 * (_left->get_particle()->vy() +
                       _right->get_particle()->vy()) +
                fac * rRL[1];
    } else {
        double rRL[3];
        rRL[0] = _right->x() - _left->x();
        rRL[1] = _right->y() - _left->y();
        rRL[2] = rRL[0] * rRL[0] + rRL[1] * rRL[1];
        double vRx, vRy;
        StateVector Wdummy(0., _left->get_particle()->vx(),
                           _left->get_particle()->vy(), 0.);
        transform(Wdummy, _left->get_position(), _right->get_position());
        Wdummy[1] = -Wdummy[1];
        invtransform(Wdummy, _left->get_position(), _right->get_position());
        vRx = Wdummy[1];
        vRy = Wdummy[2];
        double fac = (_left->get_particle()->vx() - vRx) *
                     (_midpoint[0] - 0.5 * (_left->x() + _right->x()));
        fac += (_left->get_particle()->vy() - vRy) *
               (_midpoint[1] - 0.5 * (_left->y() + _right->y()));
        fac /= rRL[2];
        _v[0] = 0.5 * (_left->get_particle()->vx() + vRx) + fac * rRL[0];
        _v[1] = 0.5 * (_left->get_particle()->vy() + vRy) + fac * rRL[1];
    }
#endif
#endif
}

#ifndef ICMAKER
/**
 * @brief Calculate the flux through the face using the given Riemann Solver and
 * update the conserved variables of the left and right cell
 *
 * @param timeline TimeLine of the simulation
 * @param solver Solver used to solve the Riemann problem at the interface
 */
void VorFace::calculate_flux(TimeLine& timeline, RiemannSolver& solver) {
    // if both particles are inactive, we don't consider the face
    if(!get_area()) {
        return;
    }
    if(_right->get_particle() == NULL &&
       _left->get_particle()->get_starttime() != timeline.get_integertime()) {
        return;
    }
    if(_right->get_particle() != NULL &&
       _right->get_particle()->get_starttime() != timeline.get_integertime() &&
       _left->get_particle()->get_starttime() != timeline.get_integertime()) {
        return;
    }

    // transform the boundaries:
    //   apply a velocity correction for the movement of the face
    //   the axes are rotated such that the midpoints of both cells lie on the
    //   positive x-axis
    //   the current cell acts as left state, the neighbour as right state
    // solve the riemann problem between left and right state
    //   the y'-component of the velocity is passively evolved
    // transform the resulting flux back to the original reference frame
    double dt = timeline.get_realtime(_left->get_particle()->get_timestep());
    if(_right->get_particle() != NULL) {
        dt = std::min(dt, timeline.get_realtime(
                                  _right->get_particle()->get_timestep()));
    }

    StateVector WL = _left->get_particle()->get_Wvec();
    StateVector WR;
    if(_right->get_particle() != NULL) {
        WR = _right->get_particle()->get_Wvec();
    }
    WL -= _v;
    WR -= _v;

    StateVector WLrec, WRrec;
    StateVector Wtemp;
    GasParticle* cell;
    VorGen* point;
    double d[ndim_];
    // right state
    if(_right->get_particle()) {
        Wtemp = WR;
        cell = _right->get_particle();
        point = _right;

        double dtp = dt;
        // add gravitational force
        if(timeline.has_gravity()) {
            for(unsigned int m = ndim_; m--;) {
                Wtemp[1 + m] +=
                        0.5 * dtp * cell->get_gravitational_acceleration()[m];
            }
        }
        Vec centroid = cell->get_centroid();
        // correction is +1*boxside, 0 or -1*boxside, depending on whether the
        // point is a periodic copy of the cell
        for(unsigned int n = ndim_; n--;) {
            double correction = cell->pos(n) - point->pos(n);
            d[n] = _midpoint[n] - centroid[n] + correction;
        }
        StateVector gradients[ndim_];
        cell->get_gradients(gradients);
        StateVector deltaW =
                0.5 * dtp * solver.get_time_derivative(Wtemp, gradients);
        Wtemp = reconstruct(Wtemp, gradients, d);
        Vec centdiff = _left->get_particle()->get_centroid() -
                       _right->get_particle()->get_centroid();
        Wtemp = limit(Wtemp, WR, WL, d, centdiff.norm());
        // sanity checks
        if(Wtemp[0] < 0. || Wtemp[ndim_ + 1] < 0.) {
            Wtemp = WR;
        }

        for(unsigned int m = ndim_ + 2; m--;) {
            Wtemp[m] += deltaW[m];
        }
        // sanity checks
        if(Wtemp[0] < 0. || Wtemp[ndim_ + 1] < 0.) {
            Wtemp -= deltaW;
        }

        // we cannot override the old WR, since we use it in the slope limiter
        // for WL
        WRrec = Wtemp;
    }

    // left state
    Wtemp = WL;
    cell = _left->get_particle();
    point = _left;

    double dtp = dt;
    // add gravitational force
    if(timeline.has_gravity()) {
        for(unsigned int m = ndim_; m--;) {
            Wtemp[1 + m] +=
                    0.5 * dtp * cell->get_gravitational_acceleration()[m];
        }
    }
    Vec centroid = cell->get_centroid();
    // correction is +1*boxside, 0 or -1*boxside, depending on whether the point
    // is a periodic copy of the cell
    for(unsigned int n = ndim_; n--;) {
        double correction = cell->pos(n) - point->pos(n);
        d[n] = _midpoint[n] - centroid[n] + correction;
    }
    StateVector gradientsL[ndim_];
    cell->get_gradients(gradientsL);
    StateVector deltaW =
            0.5 * dtp * solver.get_time_derivative(Wtemp, gradientsL);
    Wtemp = reconstruct(Wtemp, gradientsL, d);
    if(_right->get_particle()) {
        Vec centdiff = _left->get_particle()->get_centroid() -
                       _right->get_particle()->get_centroid();
        Wtemp = limit(Wtemp, WL, WR, d, centdiff.norm());
    }
    // sanity checks
    if(Wtemp[0] < 0. || Wtemp[ndim_ + 1] < 0.) {
        Wtemp = WL;
    }

    for(unsigned int m = ndim_ + 2; m--;) {
        Wtemp[m] += deltaW[m];
    }
    // sanity checks
    if(Wtemp[0] < 0. || Wtemp[ndim_ + 1] < 0.) {
        Wtemp -= deltaW;
    }

    WLrec = Wtemp;

    // now it is safe to override the old W's
    // if a state is vacuum, the time integration will fail (due to the 1/rho
    // term) and we keep the old value
    if(WL.rho()) {
        WL = WLrec;
    }
    if(WR.rho()) {
        WR = WRrec;
    }

    StateVector Whalf;
    double maxmach = 0.;
    Vec n;
    double angles[ndim_];
    get_normal(angles);
    for(unsigned int i = 0; i < ndim_; ++i) {
        n[i] = angles[i];
    }

    _left->get_particle()->set_max_mach(
            std::max(_left->get_particle()->get_max_mach(), maxmach));
    if(_right->get_particle()) {
        _right->get_particle()->set_max_mach(
                std::max(_right->get_particle()->get_max_mach(), maxmach));
    }

    StateVector flux = solver.solve_for_flux(WL, WR, n, _v,
                                             _right->get_particle() == NULL);
    StateVector dQL = dt * get_area() * flux;
    StateVector dQR = -dt * get_area() * flux;

    // flux limiter
    // if one of the particles is inactive, it is possible that we still
    // overflux. We therefore add a time correction factor as well
    double timefac = 1.;
    if(_right->get_particle()) {
        double dtmax =
                timeline.get_realtime(_left->get_particle()->get_timestep());
        dtmax = std::max(
                dtmax,
                timeline.get_realtime(_right->get_particle()->get_timestep()));
        timefac = dt / dtmax;
    }
    double areafac = get_area() / _left->get_particle()->get_total_area();
    areafac *= timefac;
    if(fabs(dQL[0]) > areafac * _left->get_particle()->get_Qvec().m()) {
        double corfac =
                fabs(areafac * _left->get_particle()->get_Qvec().m() / dQL[0]);
        dQL *= corfac;
        dQR = -1. * dQL;
    }
    if(_right->get_particle()) {
        double areafac = get_area() / _right->get_particle()->get_total_area();
        areafac *= timefac;
        if(fabs(dQR[0]) > areafac * _right->get_particle()->get_Qvec().m()) {
            double corfac = fabs(
                    areafac * _right->get_particle()->get_Qvec().m() / dQR[0]);
            dQR *= corfac;
            dQL = -1. * dQR;
        }
    }

    // to explicitly conserve the conserved quantities, we have to make sure
    // that all fluxes that are exchanged are exactly the same (up to round off)
    // fluxes will not be exactly the same for an exchange of left and right,
    // since switching left and right changes the rotation matrices (which are
    // very sensitive to round off)
    // in most cases, we solve this problem by only treating every face once
    // for ghost cells, this is tricky, since the face needs to be present on
    // both processes and will be treated twice
    // here we make sure that these faces are only treated on one of the
    // processes by only treating it on the process with the lowest rank
    if((int)_right->get_process() != MPIGlobal::rank) {
        // if left is inactive, we do not exchange fluxes; the flux exchange
        // will be done on the other process
        if(_left->get_particle()->get_starttime() !=
           timeline.get_integertime()) {
            return;
        }
        // if right is a mirror particle (very unlikely, but possible), we have
        // to do the flux exchange
        // if right is inactive, we have to do the flux exchange on this process
        // if both cells are active, we choose the one with the lowest rank and
        // exchange the fluxes on that process
        if((int)_right->get_process() > MPIGlobal::rank ||
           !_right->get_particle() ||
           _right->get_particle()->get_starttime() !=
                   timeline.get_integertime()) {
            _left->get_particle()->increase_dQ(dQL);
            if(timeline.has_gravity()) {
                _left->get_particle()->increase_delta_E(
                        dQL[0] *
                        (_left->get_position() - _right->get_position()));
            }
            if(_right->get_particle()) {
                _right->get_particle()->increase_dQ(dQR);
                if(timeline.has_gravity()) {
                    _right->get_particle()->increase_delta_E(
                            dQR[0] *
                            (_right->get_position() - _left->get_position()));
                }
            }
        }
        return;
    }

    if(_right->get_particle() != NULL) {
        // if _left is inactive and _right is a ghost => do nothing
        // the flux exchange between _left and _right will happen on another
        // processor
        // the last condition makes sure that periodic ghosts are only processed
        // once, since the id of a periodic ghost cell is not 1
        if(_left->get_particle()->get_starttime() ==
                   timeline.get_integertime() ||
           _right->get_id() != 1 ||
           (_right->get_particle()->get_position() - _right->get_position())
                   .norm2()) {
            _left->get_particle()->increase_dQ(dQL);
            if(timeline.has_gravity()) {
                _left->get_particle()->increase_delta_E(
                        dQL[0] *
                        (_left->get_position() - _right->get_position()));
            }
        }
        // if the id of _right is 1, then _right is a ghost
        // there are two possiblities:
        //   * the cell of the periodic copy also exists; the face is
        //     represented twice and will be processed twice
        //   * the cell of the periodic copy does not exist; this face lies
        //     between two inactive particles and won't be processed
        // In the second case, we never get here. In the first case, we don't
        // increase _right's dQ, as we only want to account for this face once
        if(_right->get_id() != 1) {
            _right->get_particle()->increase_dQ(dQR);
            if(timeline.has_gravity()) {
                _right->get_particle()->increase_delta_E(
                        dQR[0] *
                        (_right->get_position() - _left->get_position()));
            }
        }
    } else {
        _left->get_particle()->increase_dQ(dQL);
        if(timeline.has_gravity()) {
            _left->get_particle()->increase_delta_E(
                    dQL[0] * (_left->get_position() - _right->get_position()));
        }
    }
    if(flux[0] != flux[0]) {
#if ndim_ == 3
        cerr << "position left: " << _left->x() << "\t" << _left->y() << "\t"
             << _left->z() << endl;
        cerr << "position right: " << _right->x() << "\t" << _right->y() << "\t"
             << _right->z() << endl;
        cerr << "volume left: " << _left->get_cell()->get_volume() << endl;
        cerr << "WL:" << endl;
        cerr << _left->get_particle()->get_Wvec()[0] << "\t"
             << _left->get_particle()->get_Wvec()[1] << "\t"
             << _left->get_particle()->get_Wvec()[2] << "\t"
             << _left->get_particle()->get_Wvec()[3] << "\t"
             << _left->get_particle()->get_Wvec()[4] << endl;
        cerr << "volume right: " << _right->get_cell()->get_volume() << endl;
        cerr << "WR:" << endl;
        cerr << _right->get_particle()->get_Wvec()[0] << "\t"
             << _right->get_particle()->get_Wvec()[1] << "\t"
             << _right->get_particle()->get_Wvec()[2] << "\t"
             << _right->get_particle()->get_Wvec()[3] << "\t"
             << _right->get_particle()->get_Wvec()[4] << endl;
        cerr << "v_face: " << _v[0] << "\t" << _v[1] << "\t" << _v[2] << endl;
        cerr << "Left gradients:" << endl;
        for(unsigned int i = 0; i < ndim_ + 2; i++) {
            cerr << _left->get_particle()->get_gradient(i, 0) << "\t"
                 << _left->get_particle()->get_gradient(i, 1) << "\t"
                 << _left->get_particle()->get_gradient(i, 2) << endl;
        }
        cerr << "Right gradients:" << endl;
        for(unsigned int i = 0; i < ndim_ + 2; i++) {
            cerr << _right->get_particle()->get_gradient(i, 0) << "\t"
                 << _right->get_particle()->get_gradient(i, 1) << "\t"
                 << _right->get_particle()->get_gradient(i, 2) << endl;
        }
#else
        cerr << "position left: " << _left->x() << "\t" << _left->y() << endl;
        cerr << "position right: " << _right->x() << "\t" << _right->y()
             << endl;
        cerr << "volume left: " << _left->get_cell()->get_volume() << endl;
        cerr << "WL:" << endl;
        cerr << _left->get_particle()->get_Wvec()[0] << "\t"
             << _left->get_particle()->get_Wvec()[1] << "\t"
             << _left->get_particle()->get_Wvec()[2] << "\t"
             << _left->get_particle()->get_Wvec()[3] << endl;
        cerr << "volume right: " << _right->get_cell()->get_volume() << endl;
        cerr << "WR:" << endl;
        cerr << _right->get_particle()->get_Wvec()[0] << "\t"
             << _right->get_particle()->get_Wvec()[1] << "\t"
             << _right->get_particle()->get_Wvec()[2] << "\t"
             << _right->get_particle()->get_Wvec()[3] << endl;
        cerr << "v_face: " << _v[0] << "\t" << _v[1] << endl;
        cerr << "Left gradients:" << endl;
        for(unsigned int i = 0; i < ndim_ + 2; i++) {
            cerr << _left->get_particle()->get_gradient(i, 0) << "\t"
                 << _left->get_particle()->get_gradient(i, 1) << endl;
        }
        cerr << "Right gradients:" << endl;
        for(unsigned int i = 0; i < ndim_ + 2; i++) {
            cerr << _right->get_particle()->get_gradient(i, 0) << "\t"
                 << _right->get_particle()->get_gradient(i, 1) << endl;
        }
        Vec d = _midpoint - _right->get_cell()->get_centroid();
        double deltarho = _right->get_particle()->get_gradient(0, 0) * d[0] +
                          _right->get_particle()->get_gradient(0, 1) * d[1];
        cerr << deltarho << "\t"
             << (_right->get_particle()->get_Wvec().rho() + deltarho) << endl;
#endif
        throw FluxCalculationException(WL, WR, MPIGlobal::rank, solver);
    }
}
#endif

/**
 * @brief Calculate the advective flux through this face using the give timestep
 *
 * @warning This function does not work!
 *
 * @param dt Current timestep of the simulation
 */
void VorFace::calculate_advection_flux(double dt) {
    if(!get_area()) {
        return;
    }
    // find the upwind direction
    StateVector Whalf;
#if ndim_ == 3
    StateVector test(0., _v[0], _v[1], _v[2], 0.);
#else
    StateVector test(0., _v[0], _v[1], 0.);
#endif
    transform(test, _left->get_position(), _right->get_position());
    if(test.vx() < 0.) {
        Whalf = _left->get_particle()->get_Wvec();
    } else {
        Whalf = _right->get_particle()->get_Wvec();
    }
    StateVector flux[ndim_];
    flux[0].set_rho(_v[0] * Whalf.rho());
    flux[1].set_rho(_v[1] * Whalf.rho());
    flux[0].set_vx(_v[0] * Whalf.vx());
    flux[1].set_vx(_v[1] * Whalf.vx());
    flux[0].set_vy(_v[0] * Whalf.vy());
    flux[1].set_vy(_v[1] * Whalf.vy());
    flux[0].set_p(_v[0] * Whalf.p());
    flux[1].set_p(_v[1] * Whalf.p());
#if ndim_ == 3
    flux[2].set_rho(_v[2] * Whalf.rho());
    flux[2].set_vx(_v[2] * Whalf.vx());
    flux[2].set_vy(_v[2] * Whalf.vy());
    flux[0].set_vz(_v[0] * Whalf.vz());
    flux[1].set_vz(_v[1] * Whalf.vz());
    flux[2].set_vz(_v[2] * Whalf.vz());
    flux[2].set_p(_v[2] * Whalf.p());
#endif

    double angles[ndim_];
    get_normal(angles);
    StateVector dQL;
    StateVector dQR;
    for(unsigned int i = ndim_; i--;) {
        dQL -= dt * get_area() * angles[i] * flux[i];
        dQR += dt * get_area() * angles[i] * flux[i];
    }
    if(_right->get_particle() != NULL) {
        if(_right->get_id() != 1) {
            _right->get_particle()->increase_dQ(dQR);
        }
    }
    _left->get_particle()->increase_dQ(dQL);
}

/**
 * @brief Get the triangles that make up the polygon that is this face
 *
 * @param positions std:vector to store the triangle positions in
 * @param connectivity std::vector to store the triangle connections in
 * @return The number of vertices added to the positions list
 */
unsigned int VorFace::get_triangles(std::vector<float>& positions,
                                    std::vector<int>& connectivity) {
    unsigned int pointsize = positions.size() / 3;
    for(unsigned int i = 0; i < _vertices.size(); i++) {
        positions.push_back(_vertices[i]->x());
        positions.push_back(_vertices[i]->y());
#if ndim_ == 3
        positions.push_back(_vertices[i]->z());
#else
        positions.push_back(0.);
#endif
    }
    unsigned int n = _vertices.size();
    connectivity.push_back(n);
    for(unsigned int i = 0; i < _vertices.size(); i++) {
        connectivity.push_back(pointsize + i);
    }
    return n;
}

/**
 * @brief Constructor
 *
 * @param WL Left StateVector
 * @param WR Rigth StateVector
 * @param rank Rank of the MPI process
 * @param solver RiemannSolver used to solve the Riemann problem at the
 * interface
 */
FluxCalculationException::FluxCalculationException(StateVector WL,
                                                   StateVector WR,
                                                   unsigned int rank,
                                                   RiemannSolver& solver)
        : _solver(solver) {
    _WL = WL;
    _WR = WR;
    _rank = rank;
}

/**
 * @brief Get a human-readable message describing this Exception
 *
 * @return A comprehendable C-string
 */
const char* FluxCalculationException::what() const throw() {
    cerr << "An error occured during the flux calculation on process " << _rank
         << "!\n\n";
    cerr << "Left state:\n(" << _WL[0] << ", " << _WL[1] << ", " << _WL[2]
         << ", " << _WL[3];
#if ndim_ == 3
    cerr << ", " << _WL[4];
#endif
    cerr << ")\n\n";
    cerr << "Right state:\n(" << _WR[0] << ", " << _WR[1] << ", " << _WR[2]
         << ", " << _WR[3];
#if ndim_ == 3
    cerr << ", " << _WR[4];
#endif
    cerr << ")\n\n";
    StateVector WL = _WL;
    StateVector WR = _WR;
    double mach;
    Vec n;
    n[0] = 1.;
    StateVector sol = _solver.solve(WL, WR, n, mach);
    cerr << "Star state:\n(" << sol[0] << ", " << sol[1] << ", " << sol[2]
         << ", " << sol[3];
#if ndim_ == 3
    cerr << ", " << sol[4];
#endif
    cerr << ")\n\n";
    return "Flux calculation error";
}
