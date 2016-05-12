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
 * @file testDelTess.cpp
 *
 * @brief Unit test for the Delaunay tesselation algorithm
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "DelCont.hpp"           // for CubicBox, DelCont
#include "DelTess.hpp"           // for DelTess
#include "Error.hpp"             // for my_exit
#include "Simplex.hpp"           // for Simplex
#include "Vec.hpp"               // for Vec
#include "VorGen.hpp"            // for VorGen
#include "myAssert.hpp"          // for my_assert
#include "utilities/Cuboid.hpp"  // for Cuboid
#include <cstddef>               // for NULL
#include <iostream>              // for basic_ostream, char_traits, etc
#include <vector>                // for vector
#if ndim_ == 3
#include <list>  // for _List_iterator, list, etc
#endif
using namespace std;

/**
 * @brief Function to test the robustness of all point insertion algorithms
 *
 * The function adds points with a known id and then checks if all neighbour
 * relations are correct. In order to facilitate checks, the ids are multiples
 * of 2; the sum of ids will always be unique for every tetrahedron/triangle.
 */
#if ndim_ == 3
void DelTess::check_methods() {
    // first step: setup known large tetrahedron
    cout << "Setting up test volume..." << endl;
    Vec origin(0.5, 0.5, 0.5);
    _container = new CubicBox(origin, 1.);
    Vec anchor(-2., -2., -2.);
    Vec sides(10., 10., 10.);
    _box12 = Cuboid(anchor, sides);
    // manually add ghosts (because the ghost positions changed)
    _points.push_back(new VorGen(-2., -2., -2.));
    _points.push_back(new VorGen(8., 0., 0.));
    _points.push_back(new VorGen(0., 8., 0.));
    _points.push_back(new VorGen(0., 0., 8.));
    _points[0]->set_p12(get_p12(_points[0]->get_position()));
    _points[1]->set_p12(get_p12(_points[1]->get_position()));
    _points[2]->set_p12(get_p12(_points[2]->get_position()));
    _points[3]->set_p12(get_p12(_points[3]->get_position()));
    Simplex* tetrahedron = new Simplex(0, 1, 2, 3, _points);
    Simplex* ghosttetrahedron = new Simplex(0, 1, 2, 3, _points);
    _simplices.push_back(ghosttetrahedron);
    _simplices.push_back(tetrahedron);
    tetrahedron->set_next_check(0);
    tetrahedron->set_previous_check(0);
    _lastindex = 1;
    // set id's of ghost particles
    // the order and position of the points is determined by the CubicBox
    // we expect them to be:
    // 1 --> (-2., -2., -2.)
    // 2 --> ( 8.,  0.,  0.)
    // 4 --> ( 0.,  8.,  0.)
    // 8 --> ( 0.,  0.,  8.)
    _points[0]->set_id(1);
    _points[1]->set_id(2);
    _points[2]->set_id(4);
    _points[3]->set_id(8);
    // check for normal insertion algorithm
    cout << "Checking basic algorithm..." << endl;
    VorGen* A = new VorGen(2., 2., 2.);
    A->set_id(16);
    _points.push_back(A);
    add_point(4);
    int result[128][65];
    int faces[128][65];
    // g0g1g2A
    result[23][1] = 30;
    result[23][2] = 29;
    result[23][4] = 27;
    result[23][16] = -1;
    faces[23][1] = 8;
    faces[23][2] = 8;
    faces[23][4] = 8;
    faces[23][16] = -1;
    // g0g1g3A
    result[27][1] = 30;
    result[27][2] = 29;
    result[27][8] = 23;
    result[27][16] = -1;
    faces[27][1] = 4;
    faces[27][2] = 4;
    faces[27][8] = 4;
    faces[27][16] = -1;
    // g0g2g3A
    result[29][1] = 30;
    result[29][4] = 27;
    result[29][8] = 23;
    result[29][16] = -1;
    faces[29][1] = 2;
    faces[29][4] = 2;
    faces[29][8] = 2;
    faces[29][16] = -1;
    // g1g2g3A
    result[30][2] = 29;
    result[30][4] = 27;
    result[30][8] = 23;
    result[30][16] = -1;
    faces[30][2] = 1;
    faces[30][4] = 1;
    faces[30][8] = 1;
    faces[30][16] = -1;
    // we need to set the relations for these checks to work
    set_relations();
    for(unsigned int i = 0; i < 4; i++) {
        VorGen* P = _points[i];
        list<Simplex*> simplices = P->get_tetrahedra();
        for(list<Simplex*>::iterator s = simplices.begin();
            s != simplices.end(); s++) {
            unsigned int* Pidx = (*s)->get_vorgens();
            VorGen* PP[4] = {_points[Pidx[0]], _points[Pidx[1]],
                             _points[Pidx[2]], _points[Pidx[3]]};
            unsigned int sum = (*s)->get_idsum(_points);
            for(unsigned int i = 0; i < 4; i++) {
                if(!(*s)->get_ngb(i)) {
                    my_assert(result[sum][PP[i]->get_id()] == -1,
                              "incorrect NULL neighbour!");
                    my_assert(faces[sum][PP[i]->get_id()] == -1,
                              "NULL neighbour has ngbface!");
                } else {
                    my_assert(result[sum][PP[i]->get_id()] ==
                                      static_cast<int>(
                                              _simplices[(*s)->get_ngb(i)]
                                                      ->get_idsum(_points)),
                              "Incorrect neighbour!");
                    unsigned long ngbid =
                            _points[_simplices[(*s)->get_ngb(i)]->vorgen(
                                            (*s)->get_ngbface(i))]
                                    ->get_id();
                    my_assert(faces[sum][PP[i]->get_id()] ==
                                      static_cast<int>(ngbid),
                              "Incorrect ngbface!");
                }
            }
        }
        // reset the simplices for the next test
        P->reset_simplices();
    }

    // check for 3 to 6 flip algorithm
    cout << "Checking n to 2n flip algorithm..." << endl;
    VorGen* B = new VorGen(0.25, 0.25, 0.25);
    B->set_id(32);
    _points.push_back(B);
    add_point(5);
    // g1g2g3A
    result[30][2] = 60;
    result[30][4] = 58;
    result[30][8] = 54;
    result[30][16] = -1;
    faces[30][2] = 32;
    faces[30][4] = 32;
    faces[30][8] = 32;
    faces[30][16] = -1;
    // g0g1g2B
    result[39][1] = 54;
    result[39][2] = 45;
    result[39][4] = 43;
    result[39][32] = -1;
    faces[39][1] = 16;
    faces[39][2] = 8;
    faces[39][4] = 8;
    faces[39][32] = -1;
    // g0g2g3B
    result[45][1] = 60;
    result[45][4] = 43;
    result[45][8] = 39;
    result[45][32] = -1;
    faces[45][1] = 16;
    faces[45][4] = 2;
    faces[45][8] = 2;
    faces[45][32] = -1;
    // g0g1g3B
    result[43][1] = 58;
    result[43][2] = 45;
    result[43][8] = 39;
    result[43][32] = -1;
    faces[43][1] = 16;
    faces[43][2] = 4;
    faces[43][8] = 4;
    faces[43][32] = -1;
    // Ag1g2B
    result[54][16] = 39;
    result[54][2] = 60;
    result[54][4] = 58;
    result[54][32] = 30;
    faces[54][16] = 1;
    faces[54][2] = 8;
    faces[54][4] = 8;
    faces[54][32] = 8;
    // Ag2g3B
    result[60][16] = 45;
    result[60][4] = 58;
    result[60][8] = 54;
    result[60][32] = 30;
    faces[60][16] = 1;
    faces[60][4] = 2;
    faces[60][8] = 2;
    faces[60][32] = 2;
    // Ag1g3B
    result[58][16] = 43;
    result[58][2] = 60;
    result[58][8] = 54;
    result[58][32] = 30;
    faces[58][16] = 1;
    faces[58][2] = 4;
    faces[58][8] = 4;
    faces[58][32] = 4;
    set_relations();
    for(unsigned int i = 0; i < 4; i++) {
        VorGen* P = _points[i];
        list<Simplex*> simplices = P->get_tetrahedra();
        for(list<Simplex*>::iterator s = simplices.begin();
            s != simplices.end(); s++) {
            unsigned int* Pidx = (*s)->get_vorgens();
            VorGen* PP[4] = {_points[Pidx[0]], _points[Pidx[1]],
                             _points[Pidx[2]], _points[Pidx[3]]};
            unsigned int sum = (*s)->get_idsum(_points);
            for(unsigned int i = 0; i < 4; i++) {
                if(!(*s)->get_ngb(i)) {
                    my_assert(result[sum][PP[i]->get_id()] == -1,
                              "incorrect NULL neighbour!");
                    my_assert(faces[sum][PP[i]->get_id()] == -1,
                              "NULL neighbour has ngbface!");
                } else {
                    my_assert(result[sum][PP[i]->get_id()] ==
                                      static_cast<int>(
                                              _simplices[(*s)->get_ngb(i)]
                                                      ->get_idsum(_points)),
                              "Incorrect neigbour!");
                    unsigned long ngbid =
                            _points[_simplices[(*s)->get_ngb(i)]->vorgen(
                                            (*s)->get_ngbface(i))]
                                    ->get_id();
                    my_assert(faces[sum][PP[i]->get_id()] ==
                                      static_cast<int>(ngbid),
                              "Incorrect ngbface!");
                }
            }
        }
        P->reset_simplices();
    }

    // check for 2 to 6 flip algorithm
    cout << "Checking 2 to 6 flip algorithm..." << endl;
    VorGen* C = new VorGen(1., 1., 1.5);
    C->set_id(64);
    _points.push_back(C);
    add_point(6);
    // g0g1g3B
    result[43][1] = 106;
    result[43][2] = 45;
    result[43][8] = 39;
    result[43][32] = -1;
    faces[43][1] = 64;
    faces[43][2] = 4;
    faces[43][8] = 4;
    faces[43][32] = -1;
    // g0g2g3B
    result[45][1] = 108;
    result[45][4] = 43;
    result[45][8] = 39;
    result[45][32] = -1;
    faces[45][1] = 64;
    faces[45][4] = 2;
    faces[45][8] = 2;
    faces[45][32] = -1;
    // g2g3BC
    result[108][4] = 106;
    result[108][8] = 116;
    result[108][32] = 92;
    result[108][64] = 45;
    faces[108][4] = 2;
    faces[108][8] = 16;
    faces[108][32] = 16;
    faces[108][64] = 1;
    // g1g3BC
    result[106][2] = 108;
    result[106][8] = 114;
    result[106][32] = 90;
    result[106][64] = 43;
    faces[106][2] = 4;
    faces[106][8] = 16;
    faces[106][32] = 16;
    faces[106][64] = 1;
    // g2g3AC
    result[92][4] = 90;
    result[92][8] = 116;
    result[92][16] = 108;
    result[92][64] = 30;
    faces[92][4] = 2;
    faces[92][8] = 32;
    faces[92][16] = 32;
    faces[92][64] = 2;
    // g1g3AC
    result[90][2] = 92;
    result[90][8] = 114;
    result[90][16] = 106;
    result[90][64] = 30;
    faces[90][2] = 4;
    faces[90][8] = 32;
    faces[90][16] = 32;
    faces[90][64] = 4;
    // g2ABC
    result[116][4] = 114;
    result[116][16] = 108;
    result[116][32] = 92;
    result[116][64] = 54;
    faces[116][4] = 2;
    faces[116][16] = 8;
    faces[116][32] = 8;
    faces[116][64] = 2;
    // g1ABC
    result[114][2] = 116;
    result[114][16] = 106;
    result[114][32] = 90;
    result[114][64] = 54;
    faces[114][2] = 4;
    faces[114][16] = 8;
    faces[114][32] = 8;
    faces[114][64] = 4;
    // g1g2g3A
    result[30][2] = 92;
    result[30][4] = 90;
    result[30][8] = 54;
    result[30][16] = -1;
    faces[30][2] = 64;
    faces[30][4] = 64;
    faces[30][8] = 32;
    faces[30][16] = -1;
    // g1g2AB
    result[54][2] = 116;
    result[54][4] = 114;
    result[54][16] = 39;
    result[54][32] = 30;
    faces[54][2] = 64;
    faces[54][4] = 64;
    faces[54][16] = 1;
    faces[54][32] = 8;
    // g0g1g2B
    result[39][1] = 54;
    result[39][2] = 45;
    result[39][4] = 43;
    result[39][32] = -1;
    faces[39][1] = 16;
    faces[39][2] = 8;
    faces[39][4] = 8;
    faces[39][32] = -1;
    set_relations();
    for(unsigned int i = 0; i < 4; i++) {
        VorGen* P = _points[i];
        list<Simplex*> simplices = P->get_tetrahedra();
        for(list<Simplex*>::iterator s = simplices.begin();
            s != simplices.end(); s++) {
            unsigned int* Pidx = (*s)->get_vorgens();
            VorGen* PP[4] = {_points[Pidx[0]], _points[Pidx[1]],
                             _points[Pidx[2]], _points[Pidx[3]]};
            unsigned int sum = (*s)->get_idsum(_points);
            for(unsigned int i = 0; i < 4; i++) {
                if(!(*s)->get_ngb(i)) {
                    my_assert(result[sum][PP[i]->get_id()] == -1,
                              "incorrect NULL neighbour!");
                    my_assert(faces[sum][PP[i]->get_id()] == -1,
                              "NULL neighbour has ngbface!");
                } else {
                    my_assert(result[sum][PP[i]->get_id()] ==
                                      static_cast<int>(
                                              _simplices[(*s)->get_ngb(i)]
                                                      ->get_idsum(_points)),
                              "Incorrect neigbour!");
                    unsigned long ngbid =
                            _points[_simplices[(*s)->get_ngb(i)]->vorgen(
                                            (*s)->get_ngbface(i))]
                                    ->get_id();
                    my_assert(faces[sum][PP[i]->get_id()] ==
                                      static_cast<int>(ngbid),
                              "Incorrect ngbface!");
                }
            }
        }
        P->reset_simplices();
    }

    // we start over again for the flipping algorithms
    for(unsigned int i = 0; i < 4; i++) {
        delete _points[i];
    }
    _points.clear();
    delete A;
    delete B;
    delete C;
    for(unsigned int i = 0; i < _simplices.size(); i++) {
        if(_simplices[i] != NULL) {
            delete _simplices[i];
        }
    }
    _simplices.clear();
    // manually add ghosts (because the ghost positions changed)
    _points.push_back(new VorGen(-2., -2., -2.));
    _points.push_back(new VorGen(8., 0., 0.));
    _points.push_back(new VorGen(0., 8., 0.));
    _points.push_back(new VorGen(0., 0., 8.));
    _points[0]->set_p12(get_p12(_points[0]->get_position()));
    _points[1]->set_p12(get_p12(_points[1]->get_position()));
    _points[2]->set_p12(get_p12(_points[2]->get_position()));
    _points[3]->set_p12(get_p12(_points[3]->get_position()));
    tetrahedron = new Simplex(0, 1, 2, 3, _points);
    ghosttetrahedron = new Simplex(0, 1, 2, 3, _points);
    _simplices.push_back(ghosttetrahedron);
    _simplices.push_back(tetrahedron);
    _lastindex = 1;
    _points[0]->set_id(1);
    _points[1]->set_id(2);
    _points[2]->set_id(4);
    _points[3]->set_id(8);

    cout << "Insertion algorithms work. Moving on to the flipping algorithms..."
         << endl;
    // check on ordinary 2 to 3 flip algorithm
    cout << "Checking 2 to 3 flip..." << endl;
    A = new VorGen(0.5, 0.5, 0.5);
    A->set_id(16);
    _points.push_back(A);
    add_point(4);
    B = new VorGen(1., 0.9, 0.4);
    B->set_id(32);
    _points.push_back(B);
    add_point(5);
    // g0g2g3A
    result[29][1] = 60;
    result[29][4] = 27;
    result[29][8] = 53;
    result[29][16] = -1;
    faces[29][1] = 32;
    faces[29][4] = 2;
    faces[29][8] = 32;
    faces[29][16] = -1;
    // g0g1g3A
    result[27][1] = 58;
    result[27][2] = 29;
    result[27][8] = 51;
    result[27][16] = -1;
    faces[27][1] = 32;
    faces[27][2] = 4;
    faces[27][8] = 32;
    faces[27][16] = -1;
    // g0g2AB
    result[53][1] = 60;
    result[53][4] = 51;
    result[53][16] = 39;
    result[53][32] = 29;
    faces[53][1] = 8;
    faces[53][4] = 2;
    faces[53][16] = 2;
    faces[53][32] = 8;
    // g0g1AB
    result[51][1] = 58;
    result[51][2] = 53;
    result[51][16] = 39;
    result[51][32] = 27;
    faces[51][1] = 8;
    faces[51][2] = 4;
    faces[51][16] = 4;
    faces[51][32] = 8;
    // g0g1g2B
    result[39][1] = 46;
    result[39][2] = 53;
    result[39][4] = 51;
    result[39][32] = -1;
    faces[39][1] = 8;
    faces[39][2] = 16;
    faces[39][4] = 16;
    faces[39][32] = -1;
    // g2g3AB
    result[60][4] = 58;
    result[60][8] = 53;
    result[60][16] = 46;
    result[60][32] = 29;
    faces[60][4] = 2;
    faces[60][8] = 1;
    faces[60][16] = 2;
    faces[60][32] = 1;
    // g1g3AB
    result[58][2] = 60;
    result[58][8] = 51;
    result[58][16] = 46;
    result[58][32] = 27;
    faces[58][2] = 4;
    faces[58][8] = 1;
    faces[58][16] = 4;
    faces[58][32] = 1;
    // g1g2g3B
    result[46][2] = 60;
    result[46][4] = 58;
    result[46][8] = 39;
    result[46][32] = -1;
    faces[46][2] = 16;
    faces[46][4] = 16;
    faces[46][8] = 1;
    faces[46][32] = -1;
    set_relations();
    for(unsigned int i = 0; i < 4; i++) {
        VorGen* P = _points[i];
        list<Simplex*> simplices = P->get_tetrahedra();
        for(list<Simplex*>::iterator s = simplices.begin();
            s != simplices.end(); s++) {
            unsigned int* Pidx = (*s)->get_vorgens();
            VorGen* PP[4] = {_points[Pidx[0]], _points[Pidx[1]],
                             _points[Pidx[2]], _points[Pidx[3]]};
            unsigned int sum = (*s)->get_idsum(_points);
            for(unsigned int i = 0; i < 4; i++) {
                if(!(*s)->get_ngb(i)) {
                    my_assert(result[sum][PP[i]->get_id()] == -1,
                              "incorrect NULL neighbour!");
                    my_assert(faces[sum][PP[i]->get_id()] == -1,
                              "NULL neighbour has ngbface!");
                } else {
                    my_assert(result[sum][PP[i]->get_id()] ==
                                      static_cast<int>(
                                              _simplices[(*s)->get_ngb(i)]
                                                      ->get_idsum(_points)),
                              "Incorrect neigbour!");
                    unsigned long ngbid =
                            _points[_simplices[(*s)->get_ngb(i)]->vorgen(
                                            (*s)->get_ngbface(i))]
                                    ->get_id();
                    my_assert(faces[sum][PP[i]->get_id()] ==
                                      static_cast<int>(ngbid),
                              "Incorrect ngbface!");
                }
            }
        }
        P->reset_simplices();
    }

    // we start over again for the flipping algorithms
    for(unsigned int i = 0; i < 4; i++) {
        delete _points[i];
    }
    _points.clear();
    delete A;
    delete B;
    for(unsigned int i = 0; i < _simplices.size(); i++) {
        if(_simplices[i] != NULL) {
            delete _simplices[i];
        }
    }
    _simplices.clear();
    // manually add ghosts (because the ghost positions changed)
    _points.push_back(new VorGen(-2., -2., -2.));
    _points.push_back(new VorGen(8., 0., 0.));
    _points.push_back(new VorGen(0., 8., 0.));
    _points.push_back(new VorGen(0., 0., 8.));
    _points[0]->set_p12(get_p12(_points[0]->get_position()));
    _points[1]->set_p12(get_p12(_points[1]->get_position()));
    _points[2]->set_p12(get_p12(_points[2]->get_position()));
    _points[3]->set_p12(get_p12(_points[3]->get_position()));
    tetrahedron = new Simplex(0, 1, 2, 3, _points);
    ghosttetrahedron = new Simplex(0, 1, 2, 3, _points);
    _simplices.push_back(ghosttetrahedron);
    _simplices.push_back(tetrahedron);
    _lastindex = 1;
    _points[0]->set_id(1);
    _points[1]->set_id(2);
    _points[2]->set_id(4);
    _points[3]->set_id(8);

    // check on 4 to 4 flip
    cout << "Checking 4 to 4 flip..." << endl;
    A = new VorGen(0., 0., 0.);
    A->set_id(16);
    _points.push_back(A);
    add_point(4);
    B = new VorGen(0.1, 2.5, 0.1);
    B->set_id(32);
    _points.push_back(B);
    add_point(5);
    // g0g2g3B
    result[45][1] = 46;
    result[45][4] = 57;
    result[45][8] = 39;
    result[45][32] = -1;
    faces[45][1] = 2;
    faces[45][4] = 16;
    faces[45][8] = 2;
    faces[45][32] = -1;
    // g0g3AB
    result[57][1] = 58;
    result[57][8] = 51;
    result[57][16] = 45;
    result[57][32] = 27;
    faces[57][1] = 2;
    faces[57][8] = 2;
    faces[57][16] = 4;
    faces[57][32] = 2;
    // g0g1g2B
    result[39][1] = 46;
    result[39][2] = 45;
    result[39][4] = 51;
    result[39][32] = -1;
    faces[39][1] = 8;
    faces[39][2] = 8;
    faces[39][4] = 16;
    faces[39][32] = -1;
    // g0g1AB
    result[51][1] = 58;
    result[51][2] = 57;
    result[51][16] = 39;
    result[51][32] = 27;
    faces[51][1] = 8;
    faces[51][2] = 8;
    faces[51][16] = 4;
    faces[51][32] = 8;
    // g1g2g3B
    result[46][2] = 45;
    result[46][4] = 58;
    result[46][8] = 39;
    result[46][32] = -1;
    faces[46][2] = 1;
    faces[46][4] = 16;
    faces[46][8] = 1;
    faces[46][32] = -1;
    // g1g3AB
    result[58][2] = 57;
    result[58][8] = 51;
    result[58][16] = 46;
    result[58][32] = 27;
    faces[58][2] = 1;
    faces[58][8] = 1;
    faces[58][16] = 4;
    faces[58][32] = 1;
    // g0g1g3A
    result[27][1] = 58;
    result[27][2] = 57;
    result[27][8] = 51;
    result[27][16] = -1;
    faces[27][1] = 32;
    faces[27][2] = 32;
    faces[27][8] = 32;
    faces[27][16] = -1;
    set_relations();
    for(unsigned int i = 0; i < 4; i++) {
        VorGen* P = _points[i];
        list<Simplex*> simplices = P->get_tetrahedra();
        for(list<Simplex*>::iterator s = simplices.begin();
            s != simplices.end(); s++) {
            unsigned int* Pidx = (*s)->get_vorgens();
            VorGen* PP[4] = {_points[Pidx[0]], _points[Pidx[1]],
                             _points[Pidx[2]], _points[Pidx[3]]};
            unsigned int sum = (*s)->get_idsum(_points);
            for(unsigned int i = 0; i < 4; i++) {
                if(!(*s)->get_ngb(i)) {
                    my_assert(result[sum][PP[i]->get_id()] == -1,
                              "incorrect NULL neighbour!");
                    my_assert(faces[sum][PP[i]->get_id()] == -1,
                              "NULL neighbour has ngbface!");
                } else {
                    my_assert(result[sum][PP[i]->get_id()] ==
                                      static_cast<int>(
                                              _simplices[(*s)->get_ngb(i)]
                                                      ->get_idsum(_points)),
                              "Incorrect neigbour!");
                    unsigned long ngbid =
                            _points[_simplices[(*s)->get_ngb(i)]->vorgen(
                                            (*s)->get_ngbface(i))]
                                    ->get_id();
                    my_assert(faces[sum][PP[i]->get_id()] ==
                                      static_cast<int>(ngbid),
                              "Incorrect ngbface!");
                }
            }
        }
        P->reset_simplices();
    }

    // we start over again for the flipping algorithms
    for(unsigned int i = 0; i < 4; i++) {
        delete _points[i];
    }
    _points.clear();
    delete A;
    delete B;
    for(unsigned int i = 0; i < _simplices.size(); i++) {
        if(_simplices[i] != NULL) {
            delete _simplices[i];
        }
    }
    _simplices.clear();
    // manually add ghosts (because the ghost positions changed)
    _points.push_back(new VorGen(-2., -2., -2.));
    _points.push_back(new VorGen(8., 0., 0.));
    _points.push_back(new VorGen(0., 8., 0.));
    _points.push_back(new VorGen(0., 0., 8.));
    _points[0]->set_p12(get_p12(_points[0]->get_position()));
    _points[1]->set_p12(get_p12(_points[1]->get_position()));
    _points[2]->set_p12(get_p12(_points[2]->get_position()));
    _points[3]->set_p12(get_p12(_points[3]->get_position()));
    tetrahedron = new Simplex(0, 1, 2, 3, _points);
    ghosttetrahedron = new Simplex(0, 1, 2, 3, _points);
    _simplices.push_back(ghosttetrahedron);
    _simplices.push_back(tetrahedron);
    _lastindex = 1;
    _points[0]->set_id(1);
    _points[1]->set_id(2);
    _points[2]->set_id(4);
    _points[3]->set_id(8);

    // check on 3 to 2 flip
    // we also have a 2 to 3 flip in the process, but since this one is already
    // ok, this is no problem
    cout << "Checking 3 to 2 flip..." << endl;
    A = new VorGen(0., 0., 0.);
    A->set_id(16);
    _points.push_back(A);
    add_point(4);
    B = new VorGen(0.2, 2.5, 0.1);
    B->set_id(32);
    _points.push_back(B);
    add_point(5);
    // g0g2g3B
    result[45][1] = 46;
    result[45][4] = 57;
    result[45][8] = 39;
    result[45][32] = -1;
    faces[45][1] = 2;
    faces[45][4] = 16;
    faces[45][8] = 2;
    faces[45][32] = -1;
    // g0g3AB
    result[57][1] = 58;
    result[57][8] = 51;
    result[57][16] = 45;
    result[57][32] = 27;
    faces[57][1] = 2;
    faces[57][8] = 2;
    faces[57][16] = 4;
    faces[57][32] = 2;
    // g0g1g2B
    result[39][1] = 46;
    result[39][2] = 45;
    result[39][4] = 51;
    result[39][32] = -1;
    faces[39][1] = 8;
    faces[39][2] = 8;
    faces[39][4] = 16;
    faces[39][32] = -1;
    // g0g1AB
    result[51][1] = 58;
    result[51][2] = 57;
    result[51][16] = 39;
    result[51][32] = 27;
    faces[51][1] = 8;
    faces[51][2] = 8;
    faces[51][16] = 4;
    faces[51][32] = 8;
    // g1g2g3B
    result[46][2] = 45;
    result[46][4] = 58;
    result[46][8] = 39;
    result[46][32] = -1;
    faces[46][2] = 1;
    faces[46][4] = 16;
    faces[46][8] = 1;
    faces[46][32] = -1;
    // g1g3AB
    result[58][2] = 57;
    result[58][8] = 51;
    result[58][16] = 46;
    result[58][32] = 27;
    faces[58][2] = 1;
    faces[58][8] = 1;
    faces[58][16] = 4;
    faces[58][32] = 1;
    // g0g1g3A
    result[27][1] = 58;
    result[27][2] = 57;
    result[27][8] = 51;
    result[27][16] = -1;
    faces[27][1] = 32;
    faces[27][2] = 32;
    faces[27][8] = 32;
    faces[27][16] = -1;
    set_relations();
    for(unsigned int i = 0; i < 4; i++) {
        VorGen* P = _points[i];
        list<Simplex*> simplices = P->get_tetrahedra();
        for(list<Simplex*>::iterator s = simplices.begin();
            s != simplices.end(); s++) {
            unsigned int* Pidx = (*s)->get_vorgens();
            VorGen* PP[4] = {_points[Pidx[0]], _points[Pidx[1]],
                             _points[Pidx[2]], _points[Pidx[3]]};
            unsigned int sum = (*s)->get_idsum(_points);
            for(unsigned int i = 0; i < 4; i++) {
                if(!(*s)->get_ngb(i)) {
                    my_assert(result[sum][PP[i]->get_id()] == -1,
                              "incorrect NULL neighbour!");
                    my_assert(faces[sum][PP[i]->get_id()] == -1,
                              "NULL neighbour has ngbface!");
                } else {
                    my_assert(result[sum][PP[i]->get_id()] ==
                                      static_cast<int>(
                                              _simplices[(*s)->get_ngb(i)]
                                                      ->get_idsum(_points)),
                              "Incorrect neigbour!");
                    unsigned long ngbid =
                            _points[_simplices[(*s)->get_ngb(i)]->vorgen(
                                            (*s)->get_ngbface(i))]
                                    ->get_id();
                    my_assert(faces[sum][PP[i]->get_id()] ==
                                      static_cast<int>(ngbid),
                              "Incorrect ngbface!");
                }
            }
        }
        P->reset_simplices();
    }

    cout << "Checks finished, no problems encountered." << endl;

    delete _container;
}
#else
void DelTess::check_methods() {
    // first step: setup known large tetrahedron
    _points.resize(3, NULL);
    cout << "Setting up test volume..." << endl;
    Vec pos(0., 0.);
    _container = new CubicBox(pos, 1.);
    Vec anchor(-2., -2.);
    Vec sides(10., 10.);
    _box12 = Cuboid(anchor, sides);
    add_ghosts();
    // set id's of ghost particles
    // the order and position of the points is determined by the CubicBox
    // we expect them to be:
    // 1 --> (-2., -2.)
    // 2 --> ( 8.,  0.)
    // 4 --> ( 0.,  8.)
    _points[0]->set_id(1);
    _points[1]->set_id(2);
    _points[2]->set_id(4);
    // check for normal insertion algorithm
    cout << "Checking basic algorithm..." << endl;
    VorGen* A = new VorGen(2., 2.);
    A->set_id(8);
    _points.push_back(A);
    add_point(3);
    int result[64][64];
    int faces[64][33];
    // g0g1A
    result[11][1] = 14;
    result[11][2] = 13;
    result[11][8] = -1;
    faces[11][1] = 4;
    faces[11][2] = 4;
    faces[11][8] = -1;
    // g0g2A
    result[13][1] = 14;
    result[13][4] = 11;
    result[13][8] = -1;
    faces[13][1] = 2;
    faces[13][4] = 2;
    faces[13][8] = -1;
    // g1g2A
    result[14][2] = 13;
    result[14][4] = 11;
    result[14][8] = -1;
    faces[14][2] = 1;
    faces[14][4] = 1;
    faces[14][8] = -1;
    for(unsigned int i = 1; i < _simplices.size(); i++) {
        unsigned int* Pidx = _simplices[i]->get_vorgens();
        VorGen* P[3] = {_points[Pidx[0]], _points[Pidx[1]], _points[Pidx[2]]};
        unsigned int sum = _simplices[i]->get_idsum(_points);
        for(unsigned int j = 0; j < 3; j++) {
            if(!_simplices[i]->get_ngb(j)) {
                my_assert(result[sum][P[j]->get_id()] == -1,
                          "incorrect NULL neighbour!");
                my_assert(faces[sum][P[j]->get_id()] == -1,
                          "NULL neighbour has ngbface!");
            } else {
                my_assert(result[sum][P[j]->get_id()] ==
                                  static_cast<int>(
                                          _simplices[_simplices[i]->get_ngb(j)]
                                                  ->get_idsum(_points)),
                          "Incorrect neigbour!");
                unsigned long ngbid =
                        _points[_simplices[_simplices[i]->get_ngb(j)]->vorgen(
                                        _simplices[i]->get_ngbface(j))]
                                ->get_id();
                my_assert(faces[sum][P[j]->get_id()] == static_cast<int>(ngbid),
                          "Incorrect ngbface!");
            }
        }
    }

    // Next up: degenerate case
    cout << "Checking degenerate case of insertion algorithm..." << endl;
    VorGen* B = new VorGen(1., 1.);
    B->set_id(16);
    _points.push_back(B);
    add_point(4);
    // g0g2B
    result[21][1] = 28;
    result[21][4] = 19;
    result[21][16] = -1;
    faces[21][1] = 8;
    faces[21][4] = 2;
    faces[21][16] = -1;
    // g2AB
    result[28][4] = 26;
    result[28][8] = 21;
    result[28][16] = 14;
    faces[28][4] = 2;
    faces[28][8] = 1;
    faces[28][16] = 2;
    // g0g1B
    result[19][1] = 26;
    result[19][2] = 21;
    result[19][16] = -1;
    faces[19][1] = 8;
    faces[19][2] = 4;
    faces[19][16] = -1;
    // g1AB
    result[26][2] = 28;
    result[26][8] = 19;
    result[26][16] = 14;
    faces[26][2] = 4;
    faces[26][8] = 1;
    faces[26][16] = 4;
    // g1g2A
    result[14][2] = 28;
    result[14][4] = 26;
    result[14][8] = -1;
    faces[14][2] = 16;
    faces[14][4] = 16;
    faces[14][8] = -1;
    for(unsigned int i = 1; i < _simplices.size(); i++) {
        unsigned int* Pidx = _simplices[i]->get_vorgens();
        VorGen* P[3] = {_points[Pidx[0]], _points[Pidx[1]], _points[Pidx[2]]};
        unsigned int sum = _simplices[i]->get_idsum(_points);
        for(unsigned int j = 0; j < 3; j++) {
            if(!_simplices[i]->get_ngb(j)) {
                my_assert(result[sum][P[j]->get_id()] == -1,
                          "incorrect NULL neighbour!");
                my_assert(faces[sum][P[j]->get_id()] == -1,
                          "NULL neighbour has ngbface!");
            } else {
                my_assert(result[sum][P[j]->get_id()] ==
                                  static_cast<int>(
                                          _simplices[_simplices[i]->get_ngb(j)]
                                                  ->get_idsum(_points)),
                          "Incorrect neigbour!");
                unsigned long ngbid =
                        _points[_simplices[_simplices[i]->get_ngb(j)]->vorgen(
                                        _simplices[i]->get_ngbface(j))]
                                ->get_id();
                my_assert(faces[sum][P[j]->get_id()] == static_cast<int>(ngbid),
                          "Incorrect ngbface!");
            }
        }
    }

    // Last check: flipping algorithm
    cout << "Checking flip algorithm..." << endl;
    VorGen* C = new VorGen(0.6, 0.7);
    C->set_id(32);
    _points.push_back(C);
    add_point(5);
    // g0g2C
    result[37][1] = 52;
    result[37][4] = 35;
    result[37][32] = -1;
    faces[37][1] = 16;
    faces[37][4] = 2;
    faces[37][32] = -1;
    // g2BC
    result[52][4] = 50;
    result[52][16] = 37;
    result[52][32] = 28;
    faces[52][4] = 2;
    faces[52][16] = 1;
    faces[52][32] = 8;
    // g2AB
    result[28][4] = 26;
    result[28][8] = 52;
    result[28][16] = 14;
    faces[28][4] = 2;
    faces[28][8] = 32;
    faces[28][16] = 2;
    // g0g1C
    result[35][1] = 50;
    result[35][2] = 37;
    result[35][32] = -1;
    faces[35][1] = 16;
    faces[35][2] = 4;
    faces[35][32] = -1;
    // g1BC
    result[50][2] = 52;
    result[50][16] = 35;
    result[50][32] = 26;
    faces[50][2] = 4;
    faces[50][16] = 1;
    faces[50][32] = 8;
    // g1AB
    result[26][2] = 28;
    result[26][8] = 50;
    result[26][16] = 14;
    faces[26][2] = 4;
    faces[26][8] = 32;
    faces[26][16] = 4;
    // g1g2A
    result[14][2] = 28;
    result[14][4] = 26;
    result[14][8] = -1;
    faces[14][2] = 16;
    faces[14][4] = 16;
    faces[14][8] = -1;
    for(unsigned int i = 1; i < _simplices.size(); i++) {
        unsigned int* Pidx = _simplices[i]->get_vorgens();
        VorGen* P[3] = {_points[Pidx[0]], _points[Pidx[1]], _points[Pidx[2]]};
        unsigned int sum = _simplices[i]->get_idsum(_points);
        for(unsigned int j = 0; j < 3; j++) {
            if(!_simplices[i]->get_ngb(j)) {
                my_assert(result[sum][P[j]->get_id()] == -1,
                          "incorrect NULL neighbour!");
                my_assert(faces[sum][P[j]->get_id()] == -1,
                          "NULL neighbour has ngbface!");
            } else {
                my_assert(result[sum][P[j]->get_id()] ==
                                  static_cast<int>(
                                          _simplices[_simplices[i]->get_ngb(j)]
                                                  ->get_idsum(_points)),
                          "Incorrect neigbour!");
                unsigned long ngbid =
                        _points[_simplices[_simplices[i]->get_ngb(j)]->vorgen(
                                        _simplices[i]->get_ngbface(j))]
                                ->get_id();
                my_assert(faces[sum][P[j]->get_id()] == static_cast<int>(ngbid),
                          "Incorrect ngbface!");
            }
        }
    }

    cout << "Finished testing, all tests passed!" << endl;

    delete _container;
}
#endif

/**
 * @brief Delaunay tesselation algorithm test
 *
 * @param argc Number of command line arguments (ignored)
 * @param argv Command line arguments (ignored)
 * @return 0 on succes. Aborts otherwise
 */
int main(int argc, char** argv) {

    DelTess delaunay(NULL, 0);
    delaunay.check_methods();
    return 0;
}
