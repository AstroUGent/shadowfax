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
 * @file Simplex.cpp
 *
 * @brief A general simplex in X dimensions. If X=2, this is a triangle, if X=3,
 * we call it a tetrahedron: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Simplex.hpp"
#include "ExArith.h"
#include "Line.hpp"
#include "Plane.hpp"
#include "VorGen.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <ostream>
#include <utility>
#include <vector>
using namespace std;

#if ndim_ == 2
/**
 * @brief Constructor
 *
 * @param id1 Index of the first vertex in the DelTess VorGen list
 * @param id2 Index of the second vertex in the DelTess VorGen list
 * @param id3 Index of the third vertex in the DelTess VorGen list
 */
Simplex::Simplex(unsigned int id1, unsigned int id2, unsigned int id3) {
    _vorgens[0] = id1;
    _vorgens[1] = id2;
    _vorgens[2] = id3;
    _ngbs[0] = 0;
    _ngbs[1] = 0;
    _ngbs[2] = 0;
    _ngbface[0] = -1;
    _ngbface[1] = -1;
    _ngbface[2] = -1;
    _midpoint_circumsphere = NULL;
}
#endif

#if ndim_ == 3
/**
* @brief Check if the given VorGen lies inside this Simplex
*
* To prevent roundoff error from skewing the result, we use arbitrary
* precision arithmetics to test whether the point lies above the planes
* through the vertices.
*
* To ensure speed, we use barycentric coordinates and only use the more slow
* (and more exact) plane test if the barycentric coordinates could give the
* wrong answer. This code is based on
* http://www.blackpawn.com/texts/pointinpoly/default.html, but
* adapted to 3 dimensions.
*
* It should be possible to speed this up even more, but this is not of our
* present concern
*
* The return value also indicates special cases:
*  - 0: point is outside tetrahedron
*  - 1: point is inside tetrahedron
*  - 2: point is on face of tetrahedron
*  - 3: point is on edge of tetrahedron
*
* @param point VorGen to be checked
* @param points Reference to the DelTess VorGen list
* @return Value in the range 0-3. 0 indicates a point outside the tetrahedron
*/
/*
 * Test whether the given Point is inside this Tetrahedron.
 * The resulting int indicates the case at hand:
 */
int Simplex::inside(VorGen* point, vector<VorGen*>& points) {
    VorGen* points0 = points[_vorgens[0]];
    VorGen* points1 = points[_vorgens[1]];
    VorGen* points2 = points[_vorgens[2]];
    VorGen* points3 = points[_vorgens[3]];
    Vec p0 = points0->get_p12();
    Vec p1 = points1->get_p12();
    Vec p2 = points2->get_p12();
    Vec p3 = points3->get_p12();
    Vec p4 = point->get_p12();
    int inside = 0;
    double abce, acde, adbe, bdce;
    abce = ExactArithmetic::orient3d(p0, p1, p2, p4);
    acde = ExactArithmetic::orient3d(p0, p2, p3, p4);
    adbe = ExactArithmetic::orient3d(p0, p3, p1, p4);
    bdce = ExactArithmetic::orient3d(p1, p3, p2, p4);
    if(abce <= 0 && acde <= 0 && adbe <= 0 && bdce <= 0) {
        inside = 1;
        if(abce == 0) {
            inside++;
        }
        if(acde == 0) {
            inside++;
        }
        if(adbe == 0) {
            inside++;
        }
        if(bdce == 0) {
            inside++;
        }
    }
    return inside;
}
#else
/**
 * @brief Check if the given point lies inside the triangle
 *
 * We always use arbitrary precision arithmetics for this.
 *
 * The return value indicates special cases:
 *  - 0 point lies outside triangle
 *  - 1 point lies inside triangle
 *  - 2 point lies on a side of the triangle
 *
 * In principle, return values 3 and 4 are possible, but in these cases the
 * point lies on a corner of the triangle, in which case it coincides with one
 * of the vertices, or on top of all the vertices, in which case the triangle is
 * a single point. These cases should not occur for well-behaved mesh
 * generators.
 *
 * @param point VorGen to test
 * @param points Reference to the DelTess VorGen list
 * @return Value in the range 0-2. 0 means a point outside the triangle
 */
int Simplex::inside(VorGen* point, vector<VorGen*>& points) {
    int inside = 0;
    Vec a = points[_vorgens[0]]->get_p12();
    Vec b = points[_vorgens[1]]->get_p12();
    Vec c = points[_vorgens[2]]->get_p12();
    Vec d = point->get_p12();
    double abd, acd, bcd;
    abd = ExactArithmetic::orient2d(a, b, d);
    acd = ExactArithmetic::orient2d(a, c, d);
    bcd = ExactArithmetic::orient2d(b, c, d);
    if(abd >= 0 && acd <= 0 && bcd >= 0) {
        inside++;
        // the code below make sure that
        //  (a) inside is odd if the point lies on multiple faces, we throw an
        //      error in this case
        //  (b) inside/2-1 gives the index of the ngb that shares the face on
        //      which d resides
        if(!abd) {
            inside += 5;
        }
        if(!acd) {
            inside += 3;
        }
        if(!bcd) {
            inside++;
        }
        if(inside > 1 && inside % 2) {
            cerr << "Error! Points have exactly the same coordinates!" << endl;
            my_exit();
        }
    }
    return inside;
}
#endif

#if ndim_ == 3
/**
 * @brief Get the midpoint of the circumsphere of the tetrahedron
 *
 * This is detemined by the intersection of 3 of the 4 orthogonal bisectors of
 * the edges of the tetrahedron.
 *
 * The midpoint is only calculated once and then stored for later use.
 *
 * @param points Reference to the DelTess VorGen list
 * @return Pointer to the midpoint of the circumsphere of the tetrahedron
 */
VorGen* Simplex::get_special_point(vector<VorGen*>& points) {
    // based on wikipedia formula
    if(_midpoint_circumsphere == NULL) {
        double r1[3], r2[3], r3[3];
        VorGen* points0 = points[_vorgens[0]];
        VorGen* points1 = points[_vorgens[1]];
        VorGen* points2 = points[_vorgens[2]];
        VorGen* points3 = points[_vorgens[3]];
        r1[0] = points1->x() - points0->x();
        r1[1] = points1->y() - points0->y();
        r1[2] = points1->z() - points0->z();
        r2[0] = points2->x() - points0->x();
        r2[1] = points2->y() - points0->y();
        r2[2] = points2->z() - points0->z();
        r3[0] = points3->x() - points0->x();
        r3[1] = points3->y() - points0->y();
        r3[2] = points3->z() - points0->z();
        double fac1 = r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2];
        double fac2 = r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2];
        double fac3 = r3[0] * r3[0] + r3[1] * r3[1] + r3[2] * r3[2];
        double R[3];
        R[0] = fac3 * (r1[1] * r2[2] - r1[2] * r2[1]) +
               fac2 * (r3[1] * r1[2] - r3[2] * r1[1]) +
               fac1 * (r2[1] * r3[2] - r2[2] * r3[1]);
        R[1] = fac3 * (r1[2] * r2[0] - r1[0] * r2[2]) +
               fac2 * (r3[2] * r1[0] - r3[0] * r1[2]) +
               fac1 * (r2[2] * r3[0] - r2[0] * r3[2]);
        R[2] = fac3 * (r1[0] * r2[1] - r1[1] * r2[0]) +
               fac2 * (r3[0] * r1[1] - r3[1] * r1[0]) +
               fac1 * (r2[0] * r3[1] - r2[1] * r3[0]);
        double V = 2. * (r1[0] * r2[1] * r3[2] + r1[1] * r2[2] * r3[0] +
                         r1[2] * r2[0] * r3[1] - r1[2] * r2[1] * r3[0] -
                         r2[2] * r3[1] * r1[0] - r3[2] * r1[1] * r2[0]);
        if(V < 1.e-13) {
            Vec position;
            position[0] = 0.25 * (points0->x() + points1->x() + points2->x() +
                                  points3->x());
            position[1] = 0.25 * (points0->y() + points1->y() + points2->y() +
                                  points3->y());
            position[2] = 0.25 * (points0->z() + points1->z() + points2->z() +
                                  points3->z());
            delete _midpoint_circumsphere;
            _midpoint_circumsphere =
                    new VorGen(position[0], position[1], position[2]);
        } else {
            double Vinv = 1. / V;
            _midpoint_circumsphere = new VorGen(points0->x() + R[0] * Vinv,
                                                points0->y() + R[1] * Vinv,
                                                points0->z() + R[2] * Vinv);
        }
    }
    return _midpoint_circumsphere;
}
#else
/**
 * @brief Get the midpoint of the circumcircle for this triangle
 *
 * The midpoint is calculated using formulae found on wikipedia.
 *
 * The midpoint is only calculated once and then stored for later use.
 *
 * @param points Reference to the DelTess VorGen list
 * @return A pointer to the midpoint of the circumcircle of the triangle
 */
VorGen* Simplex::get_special_point(vector<VorGen*>& points) {
    if(_midpoint_circumsphere == NULL) {
        VorGen* points0 = points[_vorgens[0]];
        VorGen* points1 = points[_vorgens[1]];
        VorGen* points2 = points[_vorgens[2]];
        Vec a = points1->get_position() - points0->get_position();
        Vec b = points2->get_position() - points0->get_position();
        double D = 2. * (a[0] * b[1] - a[1] * b[0]);
        double a2 = a.norm2();
        double b2 = b.norm2();
        double R[2];
        R[0] = (b[1] * a2 - a[1] * b2) / D;
        R[1] = (a[0] * b2 - b[0] * a2) / D;
        _midpoint_circumsphere =
                new VorGen(points0->x() + R[0], points0->y() + R[1]);
    }
    return _midpoint_circumsphere;
}
#endif

/**
 * @brief Print the tetrahedron (in ascii) to the given stream
 *
 * Format:
\verbatim
t:\t(c1,c2,c3)\t(d1,d2,d3)\t(e1,e2,e3)\t(f1,f2,f3)\n
\endverbatim
 * (with the ci, di, ei, fi doubles)
 *
 * This method calls the print method of the 4 points
 *
 * @param stream std::ostream to write to
 * @param points Reference to the DelTess VorGen list
 */
void Simplex::print(ostream& stream, vector<VorGen*>& points) {
    VorGen* points0 = points[_vorgens[0]];
    VorGen* points1 = points[_vorgens[1]];
    VorGen* points2 = points[_vorgens[2]];
#if ndim_ == 3
    VorGen* points3 = points[_vorgens[3]];
#endif
    points0->print(stream);
    stream << "\n";
    points1->print(stream);
    stream << "\n";
    points2->print(stream);
    stream << "\n";
    points0->print(stream);
#if ndim_ == 3
    stream << "\n";
    points3->print(stream);
    stream << "\n";
    points1->print(stream);
    stream << "\n";
    points2->print(stream);
    stream << "\n";
    points3->print(stream);
#endif
    stream << "\n\n";
}

/**
 * @brief Calculate the centroid of the simplex
 *
 * @param centroid Array to store the result in
 * @param points Reference to the DelTess VorGen list
 */
void Simplex::get_centroid(double* centroid, vector<VorGen*>& points) {
    for(unsigned int i = ndim_; i--;) {
        centroid[i] = 0;
        for(unsigned int j = ndim_ + 1; j--;) {
            centroid[i] += points[_vorgens[j]]->pos(i);
        }
        centroid[i] /= (ndim_ + 1);
    }
}

/**
 * @brief Get the sum of the indices of the vertices of this simplex
 *
 * This method can be used to test if the simplex contains the correct vertices
 * in an idealized setup. This is useful to test the mesh construction
 * algorithms.
 *
 * @param points Reference to the DelTess VorGen list
 * @return The sum of the indices of the vertices of this simplex
 */
unsigned int Simplex::get_idsum(vector<VorGen*>& points) {
    unsigned int sum = 0;
    for(unsigned int i = ndim_ + 1; i--;) {
        sum += points[_vorgens[i]]->get_id();
    }
    return sum;
}
