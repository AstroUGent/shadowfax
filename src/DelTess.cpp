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
 * @file DelTess.cpp
 *
 * @brief Delaunay tesselation implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "DelTess.hpp"
#include "DelCont.hpp"
#include "Error.hpp"
#include "ExArith.hpp"
#include "Line.hpp"
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include "ProgramLog.hpp"
#include "Simplex.hpp"
#include "VorCell.hpp"
#include "VorGen.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/HelperFunctions.hpp"
#include "utilities/Tree.hpp"
#include "utilities/TreeWalker.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>
using namespace std;

// size of the buffer to store simplices in during the find_simplex routine
// normally, a point is located inside a single simplex and 1 would suffice
// however, in degenerate cases, the point can lie on the common face of two
// simplices, or even on the common edge of N simplices (3D). The buffer should
// at least have size N. You reasonably expect N to be rather small, but for
// degenerate meshes (e.g. cartesian grids), N can get quite large (I know of
// at least 1 example where N > 100).
// Increasing this number can solve the "encountered an axis with more than N
// common tetrahedra" error.
#define FIND_SIMPLEX_BUFFER_SIZE 1000

/**
 * @brief Convert the given position to a position in the range [1,2] (for exact
 * geometrical tests)
 *
 * @param pos General coordinates inside the box of the DelTess
 * @return Rescaled coordinates in the range [1,2]
 */
Vec DelTess::get_p12(Vec pos) {
    Vec p12(pos);
    for(unsigned int i = 0; i < ndim_; i++) {
        p12[i] = 1. + (p12[i] - _box12.get_anchor()[i]) / _box12.get_sides()[i];
    }
    return p12;
}

/**
 * @brief Constructor
 *
 * Initialize DelaunayTesselation and add ghost tetrahedron encompassing the
 * entire domain of the tesselation
 *
 * @param container A Container that defines the volume in which all particles
 * that will be added reside
 * @param numpart Number of particles that will be added to the DelTess
 * @param periodic Bool specifying if the DelTess is embedded in a periodic
 * (true) or reflective (false) box
 * @param tolerance A tolerance used to discriminate between exact and inexact
 * geometric tests
 */
DelTess::DelTess(DelCont* container, unsigned int numpart, bool periodic,
                 double tolerance) {
    _container = container;
    _voronoi = NULL;
    _lastvorgen = 0;

    vector<double> box = _container->get_bounding_tetrahedron();
#if ndim_ == 3
    Vec anchor(box[0], box[1], box[2]);
    Vec sides(box[0], box[1], box[2]);
#else
    Vec anchor(box[0], box[1]);
    Vec sides(box[0], box[1]);
#endif
    for(unsigned int j = 1; j < ndim_ + 1; j++) {
        for(unsigned int i = 0; i < ndim_; i++) {
            anchor[i] = std::min(anchor[i], box[ndim_ * j + i]);
            sides[i] = std::max(sides[i], box[ndim_ * j + i]);
        }
    }
    sides = sides - anchor;
    sides[0] = std::max(sides[0], sides[1]);
#if ndim_ == 3
    sides[0] = std::max(sides[0], sides[2]);
    sides.set(sides[0], sides[0], sides[0]);
#else
    sides.set(sides[0], sides[0]);
#endif
    anchor -= 1.e-5 * sides;
    sides *= 1.00002;
    _box12 = Cuboid(anchor, sides);

    if(container != NULL) {
        _points.resize(numpart + ndim_ + 1, NULL);
        add_ghosts();
    }
    _tolerance = tolerance;
    _lastchecked = 1;
#if ndim_ == 3
    for(unsigned int i = 0; i < 1000; i++) {
        _free_array[i] = 0;
    }
    _free_size = 0;
#endif
    _periodic = periodic;

    LOGS("DelTess constructed");
}

/**
 * \brief Destructor. Remove all points, mirrors, existing simplices, ghosts,
 * ghostscells and exportcopies.
 */
DelTess::~DelTess() {
#if ndim_ == 3
    cout << "DelTess contains " << HelperFunctions::human_readable_counter(
                                           _lastvorgen + _mirrors.size())
         << " points and "
         << HelperFunctions::human_readable_counter(_lastindex + 1 - _free_size)
         << " simplices.\n";
    cout << "This is " << ((double)(_lastindex + 1 - _free_size) /
                           (double)(_lastvorgen + _mirrors.size()))
         << " simplices per point" << endl;
#else
    cout << "DelTess contains " << HelperFunctions::human_readable_counter(
                                           _lastvorgen + _mirrors.size())
         << " points and "
         << HelperFunctions::human_readable_counter(_simplices.size())
         << " simplices.\n";
    cout << "This is " << ((double)(_simplices.size()) /
                           (double)(_lastvorgen + _mirrors.size()))
         << " simplices per point" << endl;
#endif
// pointers on indices > _lastvorgen are also contained in _mirrors and are
// deleted there
// pointers on lower indices that were never initialized are NULL and hence
// can be deleted safely
#if ndim_ == 2
    for(unsigned int i = 0; i < _points.size(); i++) {
        delete _points[i];
    }
#else
    for(unsigned int i = 0; i < _points.size(); i++) {
        delete _points[i];
    }
// not necessary: see above
//    for(unsigned int i = 0; i < _mirrors.size(); i++){
//        delete _mirrors[i];
//    }
#endif
    for(unsigned int i = _simplices.size(); i--;) {
        delete _simplices[i];
    }
    for(unsigned int i = _ghosts.size(); i--;) {
        for(unsigned int j = _ghosts[i].size(); j--;) {
            delete _ghosts[i][j];
        }
    }
    for(unsigned int i = _exportcopies.size(); i--;) {
        for(unsigned int j = _exportcopies[i].size(); j--;) {
            delete _exportcopies[i][j];
        }
    }
    for(unsigned int i = _ghostcells.size(); i--;) {
        delete _ghostcells[i];
    }

    LOGS("DelTess destructed");
}

/**
 * @brief Add pointer to the VoronoiTesselation being build with this
 * DelaunayTesselation
 *
 * Only needed for brute force point insertion with reflective boundaries if all
 * other options fail.
 * Only really needed in the special case of a cartesian grid.
 *
 * @param voronoi A VoronoiTesselation constructed based on this
 * DelaunayTesselation
 */
void DelTess::add_voronoi_tesselation(VorTess* voronoi) { _voronoi = voronoi; }

/**
 * @brief Get a pointer to the DelCont defining the tesselation domain.
 *
 * @return The DelCont used by this DelTess
 */
DelCont* DelTess::get_container() { return _container; }

#if ndim_ == 3
void DelTess::add_ghosts() {
    vector<double> coords = _container->get_bounding_tetrahedron();
    _points[_points.size() - 4] = new VorGen(coords[0], coords[1], coords[2]);
    _points[_points.size() - 3] = new VorGen(coords[3], coords[4], coords[5]);
    _points[_points.size() - 2] = new VorGen(coords[6], coords[7], coords[8]);
    _points[_points.size() - 1] = new VorGen(coords[9], coords[10], coords[11]);
    _points[_points.size() - 4]->set_p12(
            get_p12(_points[_points.size() - 4]->get_position()));
    _points[_points.size() - 3]->set_p12(
            get_p12(_points[_points.size() - 3]->get_position()));
    _points[_points.size() - 2]->set_p12(
            get_p12(_points[_points.size() - 2]->get_position()));
    _points[_points.size() - 1]->set_p12(
            get_p12(_points[_points.size() - 1]->get_position()));
    Simplex* tetrahedron =
            new Simplex(_points.size() - 4, _points.size() - 3,
                        _points.size() - 2, _points.size() - 1, _points);
    Simplex* ghosttetrahedron =
            new Simplex(_points.size() - 4, _points.size() - 3,
                        _points.size() - 2, _points.size() - 1, _points);
    _simplices.push_back(ghosttetrahedron);
    _simplices.push_back(tetrahedron);

    _simplices[0]->set_previous_check(0);
    _simplices[0]->set_next_check(0);
    _simplices[1]->set_previous_check(0);
    _simplices[1]->set_next_check(0);
    _lastindex = 1;
}
#else
/**
 * @brief Add 3 or 4 ghost points and the simplex they form
 *
 * This simplex has to encompass the entire tesselation domain to make sure that
 * every point added later at least always lies inside this simplex. The exact
 * bounds of the simplex are determined by the DelCont.
 */
void DelTess::add_ghosts() {
    vector<double> coords = _container->get_bounding_tetrahedron();
    _points[_points.size() - 3] = new VorGen(coords[0], coords[1]);
    _points[_points.size() - 2] = new VorGen(coords[2], coords[3]);
    _points[_points.size() - 1] = new VorGen(coords[4], coords[5]);
    _points[_points.size() - 3]->set_p12(
            get_p12(_points[_points.size() - 3]->get_position()));
    _points[_points.size() - 2]->set_p12(
            get_p12(_points[_points.size() - 2]->get_position()));
    _points[_points.size() - 1]->set_p12(
            get_p12(_points[_points.size() - 1]->get_position()));
    Simplex* triangle = new Simplex(_points.size() - 3, _points.size() - 2,
                                    _points.size() - 1);
    // we need to block index 0 in the _simplices vector, this is the easiest
    // way
    Simplex* ghosttriangle = new Simplex(_points.size() - 3, _points.size() - 2,
                                         _points.size() - 1);
    _simplices.push_back(ghosttriangle);
    _simplices.push_back(triangle);

    ghosttriangle->set_next_check(0);
    ghosttriangle->set_previous_check(0);
    triangle->set_next_check(0);
    triangle->set_previous_check(0);
    _lastindex = 1;
}
#endif

/**
 * @brief Add particle to the DelaunayTesselation
 *
 * This function is only a wrapper around the function DelTess::add_point. It
 * creates a VorGen corresponding to the Particle and registers it with the
 * Particle. For later use, the VorGen gets id 1.
 *
 * @param particle The GasParticle to add to the tesselation
 * @param index Index of the particle which can later be used to access the
 * corresponding VorCell
 */
void DelTess::add_particle(GasParticle* particle, unsigned int index) {
#if ndim_ == 3
    _points[index] = new VorGen(particle->x(), particle->y(), particle->z());
#else
    _points[index] = new VorGen(particle->x(), particle->y());
#endif
    _lastvorgen++;
    VorGen* point = _points[index];
    particle->set_vorgen(point);
    particle->set_vorgenid(index);
    point->set_particle(particle);
    point->set_particle_id(index);
    point->set_id(1);
    add_point(index);
}

#if ndim_ == 3
void DelTess::add_point(unsigned int index) {
    _points[index]->set_p12(get_p12(_points[index]->get_position()));
    // STEP 1: locate point
    // n will never be this large. At least, I hope so...
    unsigned int simplex_2[FIND_SIMPLEX_BUFFER_SIZE] = {0};
    unsigned int numsimplex = find_simplex(_points[index], simplex_2);
    vector<Simplex*> tetrahedron;
    for(unsigned int i = 0; i < numsimplex; i++) {
        tetrahedron.push_back(_simplices[simplex_2[i]]);
    }
    // STEP 2: process
    // if multiple tetrahedra are present, we are dealing with a degenerate case
    if(numsimplex > 1) {
        // find out which case we're dealing with...
        if(numsimplex == 2) {
            // point lies on face of tesselation: 2 to 6 flip
            two_to_six_flip(index, simplex_2);

        } else {
            // tetrahedron on edge: n to 2n flip (with n = tetrahedron.size())
            n_to_2n_flip(index, simplex_2, numsimplex);
        }
    } else {
        // the normal case: one tetrahedron is split into 4
        unsigned int* vorgens = _simplices[simplex_2[0]]->get_vorgens();
        unsigned int ngbs[4];
        unsigned int ngbfaces[4];
        for(unsigned int i = 0; i < 4; i++) {
            ngbs[i] = tetrahedron[0]->get_ngb(i);
            ngbfaces[i] = tetrahedron[0]->get_ngbface(i);
        }
        remove_tetrahedron(tetrahedron[0]);
        Simplex* new_tetrahedra[4];
        new_tetrahedra[0] =
                new Simplex(vorgens[0], vorgens[1], vorgens[2], index, _points);
        new_tetrahedra[1] =
                new Simplex(vorgens[0], vorgens[1], index, vorgens[3], _points);
        new_tetrahedra[2] =
                new Simplex(vorgens[0], index, vorgens[2], vorgens[3], _points);
        new_tetrahedra[3] =
                new Simplex(index, vorgens[1], vorgens[2], vorgens[3], _points);

        _simplices[simplex_2[0]] = new_tetrahedra[0];
        unsigned int sids[4] = {simplex_2[0], 0, 0, 0};
        // a 3 to 2 flip generates empty entries in _simplices
        // we check here if we can fill these entries
        for(unsigned int i = 1; i < 4; i++) {
            if(_free_size) {
                unsigned int free = _free_array[--_free_size];
                _simplices[free] = new_tetrahedra[i];
                sids[i] = free;
            } else {
                _simplices.push_back(new_tetrahedra[i]);
                _lastindex++;
                sids[i] = _lastindex;
            }
        }
        // a self-reference means we are at the beginning of the stack
        _simplices[sids[0]]->set_previous_check(sids[0]);
        _simplices[sids[0]]->set_next_check(sids[1]);
        _simplices[sids[1]]->set_previous_check(sids[0]);
        _simplices[sids[1]]->set_next_check(sids[2]);
        _simplices[sids[2]]->set_previous_check(sids[1]);
        _simplices[sids[2]]->set_next_check(sids[3]);
        _simplices[sids[3]]->set_previous_check(sids[2]);
        // a self-reference means we have reached the end of the stack
        _simplices[sids[3]]->set_next_check(sids[3]);

        // last tetrahedron: 123p
        int ngbface = -1;
        int index1 = 0;
        if(ngbs[0]) {
            _simplices[ngbs[0]]->change_ngb(ngbfaces[0], sids[3], index1);
            _simplices[sids[3]]->add_ngb(index1, ngbs[0], ngbfaces[0]);
        }

        // 023p
        ngbface = -1;
        index1 = 1;
        if(ngbs[1]) {
            _simplices[ngbs[1]]->change_ngb(ngbfaces[1], sids[2], index1);
            _simplices[sids[2]]->add_ngb(index1, ngbs[1], ngbfaces[1]);
        }

        // 013p
        ngbface = -1;
        index1 = 2;
        if(ngbs[2]) {
            _simplices[ngbs[2]]->change_ngb(ngbfaces[2], sids[1], index1);
            _simplices[sids[1]]->add_ngb(index1, ngbs[2], ngbfaces[2]);
        }

        // 012p
        ngbface = -1;
        index1 = 3;
        if(ngbs[3]) {
            _simplices[ngbs[3]]->change_ngb(ngbfaces[3], sids[0], index1);
            _simplices[sids[0]]->add_ngb(index1, ngbs[3], ngbfaces[3]);
        }

        // inter-tetrahedra relations
        index1 = 0;
        ngbface = 1;
        _simplices[sids[3]]->add_ngb(ngbface, sids[2], index1);
        _simplices[sids[2]]->add_ngb(index1, sids[3], ngbface);

        index1 = 0;
        ngbface = 2;
        _simplices[sids[3]]->add_ngb(ngbface, sids[1], index1);
        _simplices[sids[1]]->add_ngb(index1, sids[3], ngbface);

        index1 = 0;
        ngbface = 3;
        _simplices[sids[3]]->add_ngb(ngbface, sids[0], index1);
        _simplices[sids[0]]->add_ngb(index1, sids[3], ngbface);

        index1 = 1;
        ngbface = 2;
        _simplices[sids[2]]->add_ngb(ngbface, sids[1], index1);
        _simplices[sids[1]]->add_ngb(index1, sids[2], ngbface);

        index1 = 1;
        ngbface = 3;
        _simplices[sids[2]]->add_ngb(ngbface, sids[0], index1);
        _simplices[sids[0]]->add_ngb(index1, sids[2], ngbface);

        index1 = 2;
        ngbface = 3;
        _simplices[sids[1]]->add_ngb(ngbface, sids[0], index1);
        _simplices[sids[0]]->add_ngb(index1, sids[1], ngbface);
    }
    // STEP 3: test new tetrahedra for Delaunayhood
    unsigned int lastindex = simplex_2[0];
    while(_simplices[lastindex]->get_next_check()) {
        bool again = check_simplex(lastindex, index);
        while(again) {
            again = check_simplex(lastindex, index);
        }
        _simplices[lastindex]->set_previous_check(0);
        // we need a temporary variable, might as well update _lastchecked at
        // the same time :)
        _lastchecked = _simplices[lastindex]->get_next_check();
        // suppose lastindex contained a self-reference, then this line also
        // ensures that the loop will end :)
        _simplices[lastindex]->set_next_check(0);
        lastindex = _lastchecked;
        // we know for sure that nextindex is the new first element in the list,
        // so we set its previous to a self-reference
        _simplices[lastindex]->set_previous_check(0);
    }
    // garbage collection. If omitted, you're bound to run out of memory at some
    // point
    clear_remove_stack();
}
#else
/**
 * @brief Add a VorGen to the DelTess
 *
 * This happens in three distinct steps:
 *   - first the VorGen is located inside an existing Simnplex (or multiple)
 *   - according to the case at hand, the Simplex is split into 4 (general
 *     case), 6 or 2n new simplices (3 or 4 in 2D)
 *   - the new simplices are tested for Delaunayhood and if needed replaced
 *
 * To accelerate the first step, a first guess of the Simplex can be given
 * as the second argument of this method. If not given, this defaults to the
 * last Simplex that was tested for Delaunayhood.
 *
 * @param index The index of the VorGen to add in the internal _points vector
 */
void DelTess::add_point(unsigned int index) {
    _points[index]->set_p12(get_p12(_points[index]->get_position()));
    unsigned int simplex[100] = {0};
    unsigned int numsimplex = find_simplex(_points[index], simplex);
    if(numsimplex > 1) {
        Simplex* triangle[2] = {_simplices[simplex[0]], _simplices[simplex[1]]};
        // degenerate case
        unsigned int* vorgens0 = triangle[0]->get_vorgens();
        unsigned int* vorgens1 = triangle[1]->get_vorgens();
        unsigned int ngb[2][3];
        unsigned int ngbfaces[2][3];
        // the indices of the points that are NOT in common
        unsigned int p0 = 0, p1 = 0;
        for(unsigned int i = 0; i < 3; i++) {
            ngb[0][i] = triangle[0]->get_ngb(i);
            ngbfaces[0][i] = triangle[0]->get_ngbface(i);
            if(ngb[0][i] == simplex[1]) {
                p0 = i;
                p1 = ngbfaces[0][i];
            }
            ngb[1][i] = triangle[1]->get_ngb(i);
            ngbfaces[1][i] = triangle[1]->get_ngbface(i);
        }
        remove_tetrahedron(triangle[0]);
        remove_tetrahedron(triangle[1]);
        Simplex* new_triangles[4];
        // we might as well already order them counterclockwise
        new_triangles[0] =
                new Simplex(index, vorgens0[p0], vorgens0[(p0 + 1) % 3]);
        new_triangles[1] =
                new Simplex(index, vorgens0[(p0 + 2) % 3], vorgens0[p0]);
        new_triangles[2] =
                new Simplex(index, vorgens1[p1], vorgens1[(p1 + 1) % 3]);
        new_triangles[3] =
                new Simplex(index, vorgens1[(p1 + 2) % 3], vorgens1[p1]);

        _simplices[simplex[0]] = new_triangles[0];
        _simplices[simplex[1]] = new_triangles[1];
        _simplices.push_back(new_triangles[2]);
        _simplices.push_back(new_triangles[3]);

        new_triangles[0]->set_previous_check(simplex[0]);
        new_triangles[0]->set_next_check(simplex[1]);
        new_triangles[1]->set_previous_check(simplex[0]);
        new_triangles[1]->set_next_check(_simplices.size() - 2);
        new_triangles[2]->set_previous_check(simplex[1]);
        new_triangles[2]->set_next_check(_simplices.size() - 1);
        new_triangles[3]->set_previous_check(_simplices.size() - 2);
        new_triangles[3]->set_next_check(_simplices.size() - 1);

        // neighbours
        int ngbface = -1;
        // you know this one is zero...
        int index1 = 0;
        if(ngb[0][(p0 + 2) % 3]) {
            // this does nothing, but we need to set the face somehow
            _simplices[ngb[0][(p0 + 2) % 3]]->change_ngb(
                    ngbfaces[0][(p0 + 2) % 3], simplex[0], index1);
        }
        new_triangles[0]->add_ngb(index1, ngb[0][(p0 + 2) % 3],
                                  ngbfaces[0][(p0 + 2) % 3]);

        ngbface = -1;
        index1 = 0;
        if(ngb[0][(p0 + 1) % 3]) {
            _simplices[ngb[0][(p0 + 1) % 3]]->change_ngb(
                    ngbfaces[0][(p0 + 1) % 3], simplex[1], index1);
        }
        new_triangles[1]->add_ngb(index1, ngb[0][(p0 + 1) % 3],
                                  ngbfaces[0][(p0 + 1) % 3]);

        ngbface = -1;
        index1 = 0;
        if(ngb[1][(p1 + 2) % 3]) {
            _simplices[ngb[1][(p1 + 2) % 3]]->change_ngb(
                    ngbfaces[1][(p1 + 2) % 3], _simplices.size() - 2, index1);
        }
        new_triangles[2]->add_ngb(index1, ngb[1][(p1 + 2) % 3],
                                  ngbfaces[1][(p1 + 2) % 3]);

        ngbface = -1;
        index1 = 0;
        if(ngb[1][(p1 + 1) % 3]) {
            _simplices[ngb[1][(p1 + 1) % 3]]->change_ngb(
                    ngbfaces[1][(p1 + 1) % 3], _simplices.size() - 1, index1);
        }
        new_triangles[3]->add_ngb(index1, ngb[1][(p1 + 1) % 3],
                                  ngbfaces[1][(p1 + 1) % 3]);

        index1 = 2;
        ngbface = 1;
        new_triangles[0]->add_ngb(ngbface, _simplices.size() - 1, index1);
        new_triangles[3]->add_ngb(index1, simplex[0], ngbface);

        index1 = 1;
        ngbface = 2;
        new_triangles[1]->add_ngb(ngbface, _simplices.size() - 2, index1);
        new_triangles[2]->add_ngb(index1, simplex[1], ngbface);

        index1 = 1;
        ngbface = 2;
        new_triangles[0]->add_ngb(ngbface, simplex[1], index1);
        new_triangles[1]->add_ngb(index1, simplex[0], ngbface);

        index1 = 1;
        ngbface = 2;
        new_triangles[2]->add_ngb(ngbface, _simplices.size() - 1, index1);
        new_triangles[3]->add_ngb(index1, _simplices.size() - 2, ngbface);
    } else {
        Simplex* triangle[1] = {_simplices[simplex[0]]};
        unsigned int* vorgens = triangle[0]->get_vorgens();
        unsigned int ngbs[3];
        unsigned int ngbfaces[3];
        for(unsigned int i = 0; i < 3; i++) {
            ngbs[i] = triangle[0]->get_ngb(i);
            ngbfaces[i] = triangle[0]->get_ngbface(i);
        }
        remove_tetrahedron(triangle[0]);

        // the points are ordered in positive orientation, so we can just as
        // well use this information to our advantage and order the new
        // triangles correctly right away
        Simplex* new_triangles[3];
        new_triangles[0] = new Simplex(index, vorgens[0], vorgens[1]);
        new_triangles[1] = new Simplex(index, vorgens[2], vorgens[0]);
        new_triangles[2] = new Simplex(index, vorgens[1], vorgens[2]);
        _simplices[simplex[0]] = new_triangles[0];
        _simplices.push_back(new_triangles[1]);
        _simplices.push_back(new_triangles[2]);

        new_triangles[0]->set_previous_check(simplex[0]);
        new_triangles[0]->set_next_check(_simplices.size() - 2);
        new_triangles[1]->set_previous_check(simplex[0]);
        new_triangles[1]->set_next_check(_simplices.size() - 1);
        new_triangles[2]->set_previous_check(_simplices.size() - 2);
        new_triangles[2]->set_next_check(_simplices.size() - 1);

        // neighbours
        if(ngbs[2]) {
            _simplices[ngbs[2]]->change_ngb(ngbfaces[2], simplex[0], 0);
        }
        new_triangles[0]->add_ngb(0, ngbs[2], ngbfaces[2]);

        if(ngbs[1]) {
            _simplices[ngbs[1]]->change_ngb(ngbfaces[1], _simplices.size() - 2,
                                            0);
        }
        new_triangles[1]->add_ngb(0, ngbs[1], ngbfaces[1]);

        if(ngbs[0]) {
            _simplices[ngbs[0]]->change_ngb(ngbfaces[0], _simplices.size() - 1,
                                            0);
        }
        new_triangles[2]->add_ngb(0, ngbs[0], ngbfaces[0]);

        new_triangles[0]->add_ngb(1, _simplices.size() - 1, 2);
        new_triangles[2]->add_ngb(2, simplex[0], 1);

        new_triangles[0]->add_ngb(2, _simplices.size() - 2, 1);
        new_triangles[1]->add_ngb(1, simplex[0], 2);

        new_triangles[1]->add_ngb(2, _simplices.size() - 1, 1);
        new_triangles[2]->add_ngb(1, _simplices.size() - 2, 2);
    }

    unsigned int lastindex = simplex[0];
    while(_simplices[lastindex]->get_next_check()) {
        bool again = check_simplex(lastindex, index);
        while(again) {
            again = check_simplex(lastindex, index);
        }
        _simplices[lastindex]->set_previous_check(0);
        // we need a temporary variable, might as well update _lastchecked at
        // the same time :)
        _lastchecked = _simplices[lastindex]->get_next_check();
        // suppose lastindex contained a self-reference, then this line also
        // ensures that the loop will end :)
        _simplices[lastindex]->set_next_check(0);
        lastindex = _lastchecked;
        // we know for sure that nextindex is the new first element in the list,
        // so we set its previous to a self-reference
        _simplices[lastindex]->set_previous_check(0);
    }
    // garbage collection. If omitted, you're bound to run out of memory some
    // time
    clear_remove_stack();
}
#endif

/**
 * @brief Method to safely remove all memory used by simplices that do no
 * longer exist.
 *
 * The _remove_stack is filled during the testing and insertion of new simplices
 * by the method DelTess::remove_tetrahedron.
 */
void DelTess::clear_remove_stack() {
    for(unsigned int i = 0; i < _remove_stack.size(); i++) {
        delete _remove_stack[i];
    }
    _remove_stack.clear();
}

#if ndim_ == 3
/**
 * @brief The degenerate case where the inserted point lies on a face of the
 * tesselation
 *
 * The two tetrahedra sharing the face are replaced by six new ones.
 * The five original points are divided into two toppoints who each share three
 * tetrahedra and tree other points who share four. The only difficulty is then
 * to find out which subtetrahedron of the "upper" tetrahedron fits which
 * subtetrahedron of the "bottom" one. The rest of the dependencies can be fixed
 * by our knowledge of the indices of the toppoints in the point arrays.
 *
 * @param index The index of the VorGen that is being inserted in the internal
 * _points vector
 * @param simplex Indices of the two Simplices that share the face on which the
 * VorGen lies in the internal _simplices vector
 */
void DelTess::two_to_six_flip(unsigned int index, unsigned int* simplex) {
    vector<Simplex*> tetrahedron(2);
    for(unsigned int i = 0; i < 2; i++) {
        tetrahedron[i] = _simplices[simplex[i]];
    }

    unsigned int* vorgens1 = tetrahedron[0]->get_vorgens();
    unsigned int* vorgens2 = tetrahedron[1]->get_vorgens();
    int toppoints[2];
    toppoints[0] = 0;
    toppoints[1] = 0;
    unsigned int ngbs[2][4];
    unsigned int ngbfaces[2][4];
    for(int i = 0; i < 4; i++) {
        ngbs[0][i] = tetrahedron[0]->get_ngb(i);
        ngbfaces[0][i] = tetrahedron[0]->get_ngbface(i);
        if(vorgens1[i] != vorgens2[0] && vorgens1[i] != vorgens2[1] &&
           vorgens1[i] != vorgens2[2] && vorgens1[i] != vorgens2[3]) {
            toppoints[0] = i;
        }
        ngbs[1][i] = tetrahedron[1]->get_ngb(i);
        ngbfaces[1][i] = tetrahedron[1]->get_ngbface(i);
        if(vorgens2[i] != vorgens1[0] && vorgens2[i] != vorgens1[1] &&
           vorgens2[i] != vorgens1[2] && vorgens2[i] != vorgens1[3]) {
            toppoints[1] = i;
        }
    }
    // toppoints now holds the indices of the toppoints in their respective
    // point arrays. The indices of the other points are
    // (toppoints+1)%4, (toppoints+2)%4 and (toppoints+3)%4

    unsigned int sids[6] = {simplex[0], simplex[1], 0, 0, 0, 0};

    remove_tetrahedron(tetrahedron[0]);
    remove_tetrahedron(tetrahedron[1]);
    Simplex* new_tetrahedra[6];
    new_tetrahedra[0] = new Simplex(index, vorgens1[toppoints[0]],
                                    vorgens1[(toppoints[0] + 1) % 4],
                                    vorgens1[(toppoints[0] + 2) % 4], _points);
    new_tetrahedra[1] = new Simplex(index, vorgens1[toppoints[0]],
                                    vorgens1[(toppoints[0] + 1) % 4],
                                    vorgens1[(toppoints[0] + 3) % 4], _points);
    new_tetrahedra[2] = new Simplex(index, vorgens1[toppoints[0]],
                                    vorgens1[(toppoints[0] + 2) % 4],
                                    vorgens1[(toppoints[0] + 3) % 4], _points);
    new_tetrahedra[3] = new Simplex(index, vorgens2[toppoints[1]],
                                    vorgens2[(toppoints[1] + 1) % 4],
                                    vorgens2[(toppoints[1] + 2) % 4], _points);
    new_tetrahedra[4] = new Simplex(index, vorgens2[toppoints[1]],
                                    vorgens2[(toppoints[1] + 1) % 4],
                                    vorgens2[(toppoints[1] + 3) % 4], _points);
    new_tetrahedra[5] = new Simplex(index, vorgens2[toppoints[1]],
                                    vorgens2[(toppoints[1] + 2) % 4],
                                    vorgens2[(toppoints[1] + 3) % 4], _points);

    _simplices[simplex[0]] = new_tetrahedra[0];
    _simplices[simplex[1]] = new_tetrahedra[1];
    for(unsigned int i = 2; i < 6; i++) {
        if(_free_size) {
            unsigned int free = _free_array[--_free_size];
            _simplices[free] = new_tetrahedra[i];
            sids[i] = free;
        } else {
            _simplices.push_back(new_tetrahedra[i]);
            _lastindex++;
            sids[i] = _lastindex;
        }
    }

    // not a mistake: a self-reference means the beginning of the stack
    _simplices[sids[0]]->set_previous_check(sids[0]);
    _simplices[sids[0]]->set_next_check(sids[1]);
    _simplices[sids[1]]->set_previous_check(sids[0]);
    _simplices[sids[1]]->set_next_check(sids[2]);
    _simplices[sids[2]]->set_previous_check(sids[1]);
    _simplices[sids[2]]->set_next_check(sids[3]);
    _simplices[sids[3]]->set_previous_check(sids[2]);
    _simplices[sids[3]]->set_next_check(sids[4]);
    _simplices[sids[4]]->set_previous_check(sids[3]);
    _simplices[sids[4]]->set_next_check(sids[5]);
    _simplices[sids[5]]->set_previous_check(sids[4]);
    // not a mistake: a self-reference means the end of the stack
    _simplices[sids[5]]->set_next_check(sids[5]);

    // neighbours
    // start with the easy ones: the neighbours opposite to the inserted point
    // (and the corresponding neighbours of these neighbours)
    int index1;
    int ngbface = -1;
    for(int i = 0; i < 3; i++) {
        if(ngbs[0][(toppoints[0] + i + 1) % 4]) {
            // we know this, because we constructed the new tetrahedra in this
            // way
            index1 = 0;
            ngbface = ngbfaces[0][(toppoints[0] + i + 1) % 4];
            _simplices[ngbs[0][(toppoints[0] + i + 1) % 4]]->change_ngb(
                    ngbface, sids[2 - i], index1);
            _simplices[sids[2 - i]]->add_ngb(
                    index1, ngbs[0][(toppoints[0] + i + 1) % 4], ngbface);
        }
        if(ngbs[1][(toppoints[1] + i + 1) % 4]) {
            index1 = 0;
            ngbface = ngbfaces[1][(toppoints[1] + i + 1) % 4];
            _simplices[ngbs[1][(toppoints[1] + i + 1) % 4]]->change_ngb(
                    ngbface, sids[5 - i], index1);
            _simplices[sids[5 - i]]->add_ngb(
                    index1, ngbs[1][(toppoints[1] + i + 1) % 4], ngbface);
        }
    }

    // the other easy part: the new tetrahedra inside the original ones

    // first original tetrahedron
    index1 = _simplices[sids[0]]->get_index(vorgens1[(toppoints[0] + 2) % 4]);
    ngbface = _simplices[sids[1]]->add_ngb_from_vorgen(
            vorgens1[(toppoints[0] + 3) % 4], sids[0], index1);
    _simplices[sids[0]]->add_ngb(index1, sids[1], ngbface);

    index1 = _simplices[sids[0]]->get_index(vorgens1[(toppoints[0] + 1) % 4]);
    ngbface = _simplices[sids[2]]->add_ngb_from_vorgen(
            vorgens1[(toppoints[0] + 3) % 4], sids[0], index1);
    _simplices[sids[0]]->add_ngb(index1, sids[2], ngbface);

    index1 = _simplices[sids[1]]->get_index(vorgens1[(toppoints[0] + 1) % 4]);
    ngbface = _simplices[sids[2]]->add_ngb_from_vorgen(
            vorgens1[(toppoints[0] + 2) % 4], sids[1], index1);
    _simplices[sids[1]]->add_ngb(index1, sids[2], ngbface);
    // second original tetrahedron
    index1 = _simplices[sids[3]]->get_index(vorgens2[(toppoints[1] + 2) % 4]);
    ngbface = _simplices[sids[4]]->add_ngb_from_vorgen(
            vorgens2[(toppoints[1] + 3) % 4], sids[3], index1);
    _simplices[sids[3]]->add_ngb(index1, sids[4], ngbface);

    index1 = _simplices[sids[3]]->get_index(vorgens2[(toppoints[1] + 1) % 4]);
    ngbface = _simplices[sids[5]]->add_ngb_from_vorgen(
            vorgens2[(toppoints[1] + 3) % 4], sids[3], index1);
    _simplices[sids[3]]->add_ngb(index1, sids[5], ngbface);

    index1 = _simplices[sids[4]]->get_index(vorgens2[(toppoints[1] + 1) % 4]);
    ngbface = _simplices[sids[5]]->add_ngb_from_vorgen(
            vorgens2[(toppoints[1] + 2) % 4], sids[4], index1);
    _simplices[sids[4]]->add_ngb(index1, sids[5], ngbface);
    // and now the crappy part: inter-original-tetrahedra stuff. We need to find
    // out which tetrahedron1 fits which tetrahedron2
    // for the moment, we use the awfully heavy function get_common_tetrahedron
    // we can live with this because the 2 to 6 flip is very rarely used
    unsigned int common_tetrahedron[2];
    for(int i = 0; i < 3; i++) {
        get_common_tetrahedron(
                index, vorgens1[(toppoints[0] + i + 2 + (i == 2)) % 4],
                vorgens1[(toppoints[0] + i + 3 + (i == 2 || i == 1)) % 4], sids,
                common_tetrahedron);
        int j = 0;
        if(common_tetrahedron[0] == sids[2 - i]) {
            j = 1;
        }
        index1 = _simplices[sids[2 - i]]->get_index(vorgens1[toppoints[0]]);
        ngbface = _simplices[common_tetrahedron[j]]->add_ngb_from_vorgen(
                vorgens2[toppoints[1]], sids[2 - i], index1);
        _simplices[sids[2 - i]]->add_ngb(index1, common_tetrahedron[j],
                                         ngbface);
    }
}

/**
 * @brief Perform a 4 to 4 flip
 *
 * This degenerate case arises when the line connecting a newly inserted point
 * and the toppoint of a tetrahedron that fails the in sphere test cuts an edge
 * of the current tesselation. Furthermore, this edge has to be a common edge of
 * exactly 4 tetrahedra, otherwise nothing has to be done for this failing
 * tetrahedron.
 * The tetrahedra are ordered as such that
 *  - t1 has neighbours 2 and 3
 *  - t2 has neighbours 1 and 4
 *  - t3 has neighbours 1 and 4
 *  - t4 has neighbours 2 and 3
 *
 * p1 and p2 are the new borders of the line that has to be swapped (p1 is part
 * of t1, p2 of t2)
 *
 * p3 and p4 are the old borders, which are part of all tetrahedra
 *
 * @param v1,v2,v3,v4 Indices of VorGens as explained above
 * @param s1,s2,s3,s4 Indices of Simplices as explained above
 */
void DelTess::four_to_four_flip(unsigned int v1, unsigned int v2,
                                unsigned int v3, unsigned int v4,
                                unsigned int s1, unsigned int s2,
                                unsigned int s3, unsigned int s4) {
    unsigned int simplex[4] = {s1, s2, s3, s4};
    Simplex* t1 = _simplices[s1];
    Simplex* t2 = _simplices[s2];
    Simplex* t3 = _simplices[s3];
    Simplex* t4 = _simplices[s4];

    unsigned int ngbs[4][2];
    unsigned int ngbfaces[4][2];
    unsigned int* vorgens1 = t1->get_vorgens();
    unsigned int* vorgens3 = t3->get_vorgens();
    unsigned int topvorgens[2] = {0, 0};

    ngbs[0][0] = t1->get_ngb_from_vorgen(v3);
    ngbs[1][0] = t2->get_ngb_from_vorgen(v3);
    ngbs[2][0] = t3->get_ngb_from_vorgen(v3);
    ngbs[3][0] = t4->get_ngb_from_vorgen(v3);
    ngbs[0][1] = t1->get_ngb_from_vorgen(v4);
    ngbs[1][1] = t2->get_ngb_from_vorgen(v4);
    ngbs[2][1] = t3->get_ngb_from_vorgen(v4);
    ngbs[3][1] = t4->get_ngb_from_vorgen(v4);

    ngbfaces[0][0] = t1->get_ngbface_from_vorgen(v3);
    ngbfaces[1][0] = t2->get_ngbface_from_vorgen(v3);
    ngbfaces[2][0] = t3->get_ngbface_from_vorgen(v3);
    ngbfaces[3][0] = t4->get_ngbface_from_vorgen(v3);
    ngbfaces[0][1] = t1->get_ngbface_from_vorgen(v4);
    ngbfaces[1][1] = t2->get_ngbface_from_vorgen(v4);
    ngbfaces[2][1] = t3->get_ngbface_from_vorgen(v4);
    ngbfaces[3][1] = t4->get_ngbface_from_vorgen(v4);

    for(int i = 0; i < 4; i++) {
        if(vorgens1[i] != v1 && vorgens1[i] != v3 && vorgens1[i] != v4) {
            topvorgens[0] = vorgens1[i];
        }
        if(vorgens3[i] != v1 && vorgens3[i] != v3 && vorgens3[i] != v4) {
            topvorgens[1] = vorgens3[i];
        }
    }

    Simplex* new_tetrahedra[4];
    new_tetrahedra[0] = new Simplex(v1, v2, topvorgens[0], v3, _points);
    new_tetrahedra[1] = new Simplex(v1, v2, topvorgens[1], v3, _points);
    new_tetrahedra[2] = new Simplex(v1, v2, topvorgens[0], v4, _points);
    new_tetrahedra[3] = new Simplex(v1, v2, topvorgens[1], v4, _points);

    // we have to relink the list before we overwrite things (otherwise, strange
    // things happen when the next_check of t1 is removed)
    if(t3->get_next_check()) {
        // t3 cannot be the first simplex on the stack, since that is already
        // simplex[0] :)
        if(t3->get_next_check() != simplex[2]) {
            _simplices[t3->get_previous_check()]->set_next_check(
                    t3->get_next_check());
            _simplices[t3->get_next_check()]->set_previous_check(
                    t3->get_previous_check());
        } else {
            _simplices[t3->get_previous_check()]->set_next_check(
                    t3->get_previous_check());
        }
    }
    if(t4->get_next_check()) {
        // other_tetrahedron cannot be the first simplex on the stack, since
        // that is already simplex[0] :)
        if(t4->get_next_check() != simplex[3]) {
            _simplices[t4->get_previous_check()]->set_next_check(
                    t4->get_next_check());
            _simplices[t4->get_next_check()]->set_previous_check(
                    t4->get_previous_check());
        } else {
            _simplices[t4->get_previous_check()]->set_next_check(
                    t4->get_previous_check());
        }
    }

    _simplices[simplex[0]] = new_tetrahedra[0];
    _simplices[simplex[1]] = new_tetrahedra[1];
    _simplices[simplex[2]] = new_tetrahedra[2];
    _simplices[simplex[3]] = new_tetrahedra[3];

    _simplices[simplex[0]]->set_previous_check(simplex[0]);
    _simplices[simplex[0]]->set_next_check(simplex[1]);
    _simplices[simplex[1]]->set_previous_check(simplex[0]);
    _simplices[simplex[1]]->set_next_check(simplex[2]);
    _simplices[simplex[2]]->set_previous_check(simplex[1]);
    _simplices[simplex[2]]->set_next_check(simplex[3]);
    _simplices[simplex[3]]->set_previous_check(simplex[2]);
    if(t1->get_next_check() != simplex[0]) {
        _simplices[simplex[3]]->set_next_check(t1->get_next_check());
        _simplices[t1->get_next_check()]->set_previous_check(simplex[3]);
    } else {
        _simplices[simplex[3]]->set_next_check(simplex[3]);
    }

    // neighbours
    // there probably is a more elegant way to write this, but... this will do
    unsigned int index1 = _simplices[simplex[2]]->get_index(v2);
    unsigned int ngbface = -1;
    if(ngbs[0][0]) {
        ngbface = ngbfaces[0][0];
        _simplices[ngbs[0][0]]->change_ngb(ngbface, simplex[2], index1);
        _simplices[simplex[2]]->add_ngb(index1, ngbs[0][0], ngbface);
    }
    index1 = _simplices[simplex[2]]->get_index(v1);
    ngbface = -1;
    if(ngbs[1][0]) {
        ngbface = ngbfaces[1][0];
        _simplices[ngbs[1][0]]->change_ngb(ngbface, simplex[2], index1);
        _simplices[simplex[2]]->add_ngb(index1, ngbs[1][0], ngbface);
    }
    index1 = _simplices[simplex[3]]->get_index(v2);
    ngbface = -1;
    if(ngbs[2][0]) {
        ngbface = ngbfaces[2][0];
        _simplices[ngbs[2][0]]->change_ngb(ngbface, simplex[3], index1);
        _simplices[simplex[3]]->add_ngb(index1, ngbs[2][0], ngbface);
    }
    index1 = _simplices[simplex[3]]->get_index(v1);
    ngbface = -1;
    if(ngbs[3][0]) {
        ngbface = ngbfaces[3][0];
        _simplices[ngbs[3][0]]->change_ngb(ngbface, simplex[3], index1);
        _simplices[simplex[3]]->add_ngb(index1, ngbs[3][0], ngbface);
    }
    index1 = _simplices[simplex[0]]->get_index(v2);
    ngbface = -1;
    if(ngbs[0][1]) {
        ngbface = ngbfaces[0][1];
        _simplices[ngbs[0][1]]->change_ngb(ngbface, simplex[0], index1);
        _simplices[simplex[0]]->add_ngb(index1, ngbs[0][1], ngbface);
    }
    index1 = _simplices[simplex[0]]->get_index(v1);
    ngbface = -1;
    if(ngbs[1][1]) {
        ngbface = ngbfaces[1][1];
        _simplices[ngbs[1][1]]->change_ngb(ngbface, simplex[0], index1);
        _simplices[simplex[0]]->add_ngb(index1, ngbs[1][1], ngbface);
    }
    index1 = _simplices[simplex[1]]->get_index(v2);
    ngbface = -1;
    if(ngbs[2][1]) {
        ngbface = ngbfaces[2][1];
        _simplices[ngbs[2][1]]->change_ngb(ngbface, simplex[1], index1);
        _simplices[simplex[1]]->add_ngb(index1, ngbs[2][1], ngbface);
    }
    index1 = _simplices[simplex[1]]->get_index(v1);
    if(ngbs[3][1]) {
        ngbface = ngbfaces[3][1];
        _simplices[ngbs[3][1]]->change_ngb(ngbface, simplex[1], index1);
        _simplices[simplex[1]]->add_ngb(index1, ngbs[3][1], ngbface);
    }

    index1 = _simplices[simplex[0]]->get_index(v3);
    ngbface =
            _simplices[simplex[2]]->add_ngb_from_vorgen(v4, simplex[0], index1);
    _simplices[simplex[0]]->add_ngb(index1, simplex[2], ngbface);

    index1 = _simplices[simplex[0]]->get_index(topvorgens[0]);
    ngbface = _simplices[simplex[1]]->add_ngb_from_vorgen(topvorgens[1],
                                                          simplex[0], index1);
    _simplices[simplex[0]]->add_ngb(index1, simplex[1], ngbface);

    index1 = _simplices[simplex[1]]->get_index(v3);
    ngbface =
            _simplices[simplex[3]]->add_ngb_from_vorgen(v4, simplex[1], index1);
    _simplices[simplex[1]]->add_ngb(index1, simplex[3], ngbface);

    index1 = _simplices[simplex[2]]->get_index(topvorgens[0]);
    ngbface = _simplices[simplex[3]]->add_ngb_from_vorgen(topvorgens[1],
                                                          simplex[2], index1);
    _simplices[simplex[2]]->add_ngb(index1, simplex[3], ngbface);
}

/**
 * @brief n to 2n flip
 *
 * Degenerate case where the newly inserted point lies on an edge of the
 * current tesselation. In this case, all (n) tetrahedra sharing this edge
 * have to be replaced by 2 new tetrahedra.
 * The tetrahedra in tetrahedron are already ordered in order of occurrence when
 * rotating around the common axis (with no direction of rotation specified).
 * For every tetrahedron, we first determine the indices of the axis points in
 * its point array. This gives us the indices of the other points as well.
 * The neighbours are fixed by using two placeholders for the left and right
 * neighbours (left is in the non-specified rotational direction, right is in
 * the other direction). The placeholders are replaced by the correct tetrahedra
 * in the end, when all new tetrahedra are in place.
 *
 * @param index The index of the VorGen that is being inserted in the internal
 * _points vector
 * @param simplex Array containing the indices of the n Simplices in the
 * internal _simplices vector
 * @param n The number of Simplices passed on
 */
void DelTess::n_to_2n_flip(unsigned int index, unsigned int* simplex,
                           unsigned int n) {
    vector<Simplex*> tetrahedron;
    for(unsigned int i = 0; i < n; i++) {
        tetrahedron.push_back(_simplices[simplex[i]]);
    }

    vector<Simplex*> new_tetrahedra;
    vector<unsigned int> leftvorgens;
    vector<unsigned int> rightvorgens;
    vector<unsigned int> sids;
    // placeholder indices. Only tetrahedron.size() new simplices are added to
    // _simplices, so these indices won't exist for sure
    unsigned int size = _lastindex + n;
    // 42 is the most not random random number I can think of
    unsigned int dummy[2] = {size + 42, size + 43};
    bool has_axis = false;
    unsigned int axis_vorgens[2] = {0, 0};
    for(unsigned int i = 0; i < n; i++) {
        // determine axis and non-axis points (similar to code in
        // find_tetrahedron as it is exactly the same)
        unsigned int* vorgens = tetrahedron[i]->get_vorgens();
        unsigned int ngbs[4];
        unsigned int ngbfaces[4];
        for(int j = 0; j < 4; j++) {
            ngbs[j] = tetrahedron[i]->get_ngb(j);
            ngbfaces[j] = tetrahedron[i]->get_ngbface(j);
        }
        remove_tetrahedron(tetrahedron[i]);
        double tests[6];
        tests[0] = normcross(_points[vorgens[0]], _points[vorgens[1]],
                             _points[index]);
        tests[1] = normcross(_points[vorgens[0]], _points[vorgens[2]],
                             _points[index]);
        tests[2] = normcross(_points[vorgens[0]], _points[vorgens[3]],
                             _points[index]);
        tests[3] = normcross(_points[vorgens[1]], _points[vorgens[2]],
                             _points[index]);
        tests[4] = normcross(_points[vorgens[1]], _points[vorgens[3]],
                             _points[index]);
        tests[5] = normcross(_points[vorgens[2]], _points[vorgens[3]],
                             _points[index]);
        // find minimum
        int min = 0;
        for(int j = 1; j < 6; j++) {
            if(tests[j] < tests[min]) {
                min = j;
            }
        }
        // convert min-value to points on axis and not on axis
        int non_axis[2], axis[2];
        if(!has_axis) {
            if(min < 3) {
                axis[0] = 0;
            } else if(min < 5) {
                axis[0] = 1;
            } else {
                axis[0] = 2;
            }
            if(min == 0) {
                axis[1] = 1;
            } else if(min == 1 || min == 3) {
                axis[1] = 2;
            } else {
                axis[1] = 3;
            }
            axis_vorgens[0] = vorgens[axis[0]];
            axis_vorgens[1] = vorgens[axis[1]];
            has_axis = true;
        } else {
            int j = 0;
            while(vorgens[j] != axis_vorgens[0]) {
                j++;
            }
            axis[0] = j;
            j = 0;
            while(vorgens[j] != axis_vorgens[1]) {
                j++;
            }
            axis[1] = j;
        }
        if(min > 2) {
            non_axis[0] = 0;
        } else if(min > 0) {
            non_axis[0] = 1;
        } else {
            non_axis[0] = 2;
        }
        if(min == 5) {
            non_axis[1] = 1;
        } else if(min == 2 || min == 4) {
            non_axis[1] = 2;
        } else {
            non_axis[1] = 3;
        }
        // nicely divide up the tetrahedron in two new ones
        sids.push_back(simplex[i]);
        new_tetrahedra.push_back(new Simplex(index, vorgens[axis[0]],
                                             vorgens[non_axis[0]],
                                             vorgens[non_axis[1]], _points));
        new_tetrahedra.push_back(new Simplex(index, vorgens[axis[1]],
                                             vorgens[non_axis[0]],
                                             vorgens[non_axis[1]], _points));

        _simplices[simplex[i]] = new_tetrahedra[new_tetrahedra.size() - 2];
        if(_free_size) {
            unsigned int free = _free_array[--_free_size];
            _simplices[free] = new_tetrahedra.back();
            sids.push_back(free);
        } else {
            _simplices.push_back(new_tetrahedra.back());
            _lastindex++;
            sids.push_back(_lastindex);
        }

        // set neighbours
        int index1 = 0;
        int ngbface = -1;
        if(ngbs[axis[1]]) {
            ngbface = ngbfaces[axis[1]];
            _simplices[ngbs[axis[1]]]->change_ngb(
                    ngbface, sids[new_tetrahedra.size() - 2], index1);
        }
        _simplices[sids[new_tetrahedra.size() - 2]]->add_ngb(
                index1, ngbs[axis[1]], ngbface);

        if(ngbs[non_axis[0]] == simplex[(i + 1) % n]) {
            _simplices[sids[new_tetrahedra.size() - 2]]->add_ngb_from_vorgen(
                    vorgens[non_axis[0]], dummy[0]);
            _simplices[sids[new_tetrahedra.size() - 2]]->add_ngb_from_vorgen(
                    vorgens[non_axis[1]], dummy[1]);
            _simplices[sids[new_tetrahedra.size() - 1]]->add_ngb_from_vorgen(
                    vorgens[non_axis[0]], dummy[0]);
            _simplices[sids[new_tetrahedra.size() - 1]]->add_ngb_from_vorgen(
                    vorgens[non_axis[1]], dummy[1]);
            leftvorgens.push_back(vorgens[non_axis[0]]);
            rightvorgens.push_back(vorgens[non_axis[1]]);
            leftvorgens.push_back(vorgens[non_axis[0]]);
            rightvorgens.push_back(vorgens[non_axis[1]]);
        } else {
            _simplices[sids[new_tetrahedra.size() - 2]]->add_ngb_from_vorgen(
                    vorgens[non_axis[1]], dummy[0]);
            _simplices[sids[new_tetrahedra.size() - 2]]->add_ngb_from_vorgen(
                    vorgens[non_axis[0]], dummy[1]);
            _simplices[sids[new_tetrahedra.size() - 1]]->add_ngb_from_vorgen(
                    vorgens[non_axis[1]], dummy[0]);
            _simplices[sids[new_tetrahedra.size() - 1]]->add_ngb_from_vorgen(
                    vorgens[non_axis[0]], dummy[1]);
            leftvorgens.push_back(vorgens[non_axis[1]]);
            rightvorgens.push_back(vorgens[non_axis[0]]);
            leftvorgens.push_back(vorgens[non_axis[1]]);
            rightvorgens.push_back(vorgens[non_axis[0]]);
        }
        index1 = 0;
        ngbface = -1;
        if(ngbs[axis[0]]) {
            ngbface = ngbfaces[axis[0]];
            _simplices[ngbs[axis[0]]]->change_ngb(
                    ngbface, sids[new_tetrahedra.size() - 1], index1);
        }
        _simplices[sids[new_tetrahedra.size() - 1]]->add_ngb(
                index1, ngbs[axis[0]], ngbface);

        index1 = 1;
        ngbface = _simplices[sids[new_tetrahedra.size() - 2]]
                          ->add_ngb_from_vorgen(vorgens[axis[0]],
                                                sids[new_tetrahedra.size() - 1],
                                                index1);
        _simplices[sids[new_tetrahedra.size() - 1]]->add_ngb_from_vorgen(
                vorgens[axis[1]], sids[new_tetrahedra.size() - 2], ngbface);
    }
    for(unsigned int i = 0; i < sids.size(); i++) {
        if(i) {
            _simplices[sids[i]]->set_previous_check(sids[i - 1]);
        } else {
            _simplices[sids[i]]->set_previous_check(sids[i]);
        }
        if(i < sids.size() - 1) {
            _simplices[sids[i]]->set_next_check(sids[i + 1]);
        } else {
            _simplices[sids[i]]->set_next_check(sids[i]);
        }
        // adjust neighbours: dummy1 is placeholder for the next tetrahedron,
        // dummy2 for the previous one
        unsigned int index1, index2;
        index1 = _simplices[sids[(i + 2) % new_tetrahedra.size()]]->get_index(
                rightvorgens[(i + 2) % new_tetrahedra.size()]);
        index2 = _simplices[sids[i]]->get_index(leftvorgens[i]);
        _simplices[sids[i]]->change_ngb(
                index2, sids[(i + 2) % new_tetrahedra.size()], index1);
        index1 = _simplices[sids[(i + new_tetrahedra.size() - 2) %
                                 new_tetrahedra.size()]]
                         ->get_index(
                                 leftvorgens[(i + new_tetrahedra.size() - 2) %
                                             new_tetrahedra.size()]);
        index2 = _simplices[sids[i]]->get_index(rightvorgens[i]);
        _simplices[sids[i]]->change_ngb(
                index2,
                sids[(i + new_tetrahedra.size() - 2) % new_tetrahedra.size()],
                index1);
    }
}
#endif

#if ndim_ == 3
bool DelTess::check_simplex(unsigned int index, unsigned int vindex) {
    // one does not simply use _simplices index instead of tetrahedron,
    // since _simplices[index] changes throughout the flip...
    Simplex* tetrahedron = _simplices[index];

    // tetrahedron cannot be NULL, since we make sure there are now NULL values
    // in the stack...

    unsigned int other_vorgens[3];
    unsigned int* vorgens = tetrahedron->get_vorgens();
    unsigned int current = tetrahedron->get_index(vindex);
    other_vorgens[0] = vorgens[(current + 1) % 4];
    other_vorgens[1] = vorgens[(current + 2) % 4];
    other_vorgens[2] = vorgens[(current + 3) % 4];
    // between 2 and 4 simplices are removed and replaced by between 2 and 4 new
    // simplices
    Simplex* removedra[4] = {NULL};
    unsigned int rindex = 0;
    unsigned int simplex[2] = {index, tetrahedron->get_ngb(current)};

    if(simplex[1]) {
        Simplex* tetrahedron2 = _simplices[simplex[1]];
        unsigned int new_vorgen =
                tetrahedron2->vorgen(tetrahedron->get_ngbface(current));
        VorGen* new_point = _points[new_vorgen];
        if(tetrahedron->in_sphere(new_point, _points)) {
            removedra[rindex++] = tetrahedron;
            removedra[rindex++] = tetrahedron2;
            // test whether the line point-new point cuts the triangle formed by
            // other_points on the inside
            unsigned int intest[3] = {0, 0, 0};
            unsigned int intest_flag =
                    tetrahedron->flip_test(current, new_point, intest, _points);
            if(!intest_flag) {
                // it does: 2 to 3 flip
                unsigned int sids[3] = {simplex[0], simplex[1], 0};

                Simplex* new_tetrahedra[3];
                new_tetrahedra[0] =
                        new Simplex(vindex, new_vorgen, other_vorgens[0],
                                    other_vorgens[1], _points);
                new_tetrahedra[1] =
                        new Simplex(vindex, new_vorgen, other_vorgens[0],
                                    other_vorgens[2], _points);
                new_tetrahedra[2] =
                        new Simplex(vindex, new_vorgen, other_vorgens[1],
                                    other_vorgens[2], _points);

                _simplices[simplex[0]] = new_tetrahedra[0];
                _simplices[simplex[1]] = new_tetrahedra[1];
                if(_free_size) {
                    unsigned int free = _free_array[--_free_size];
                    _simplices[free] = new_tetrahedra[2];
                    sids[2] = free;
                } else {
                    _simplices.push_back(new_tetrahedra[2]);
                    _lastindex++;
                    sids[2] = _lastindex;
                }

                // a self-reference will stay a self-reference in this case.
                // Hurrah! Even better: we know we have a self-reference, since
                // sids[0] is the first simplex on the stack...
                _simplices[sids[0]]->set_previous_check(
                        tetrahedron->get_previous_check());
                _simplices[sids[0]]->set_next_check(sids[1]);
                _simplices[sids[1]]->set_previous_check(sids[0]);
                _simplices[sids[1]]->set_next_check(sids[2]);
                _simplices[sids[2]]->set_previous_check(sids[1]);
                // check if we are on the end of the stack
                if(tetrahedron->get_next_check() != sids[0]) {
                    _simplices[sids[2]]->set_next_check(
                            tetrahedron->get_next_check());
                    _simplices[tetrahedron->get_next_check()]
                            ->set_previous_check(sids[2]);
                } else {
                    _simplices[sids[2]]->set_next_check(sids[2]);
                }

                unsigned int ngbs[2][3];
                unsigned int ngbfaces[2][3];
                for(int k = 0; k < 3; k++) {
                    ngbs[0][k] = tetrahedron->get_ngb((current + k + 1) % 4);
                    ngbs[1][k] =
                            tetrahedron2->get_ngb_from_vorgen(other_vorgens[k]);
                    ngbfaces[0][k] =
                            tetrahedron->get_ngbface((current + k + 1) % 4);
                    ngbfaces[1][k] = tetrahedron2->get_ngbface_from_vorgen(
                            other_vorgens[k]);
                }
                // pn01
                unsigned int index1 = 0;
                int ngbface = -1;
                if(ngbs[1][2]) {
                    ngbface = ngbfaces[1][2];
                    _simplices[ngbs[1][2]]->change_ngb(ngbface, sids[0],
                                                       index1);
                    _simplices[sids[0]]->add_ngb(index1, ngbs[1][2], ngbface);
                }
                index1 = 1;
                ngbface = -1;
                if(ngbs[0][2]) {
                    ngbface = ngbfaces[0][2];
                    _simplices[ngbs[0][2]]->change_ngb(ngbface, sids[0],
                                                       index1);
                    _simplices[sids[0]]->add_ngb(index1, ngbs[0][2], ngbface);
                }
                // pn02
                index1 = 0;
                ngbface = -1;
                if(ngbs[1][1]) {
                    ngbface = ngbfaces[1][1];
                    _simplices[ngbs[1][1]]->change_ngb(ngbface, sids[1],
                                                       index1);
                    _simplices[sids[1]]->add_ngb(index1, ngbs[1][1], ngbface);
                }
                index1 = 1;
                ngbface = -1;
                if(ngbs[0][1]) {
                    ngbface = ngbfaces[0][1];
                    _simplices[ngbs[0][1]]->change_ngb(ngbface, sids[1],
                                                       index1);
                    _simplices[sids[1]]->add_ngb(index1, ngbs[0][1], ngbface);
                }
                // pn12
                index1 = 0;
                ngbface = -1;
                if(ngbs[1][0]) {
                    ngbface = ngbfaces[1][0];
                    _simplices[ngbs[1][0]]->change_ngb(ngbface, sids[2],
                                                       index1);
                    _simplices[sids[2]]->add_ngb(index1, ngbs[1][0], ngbface);
                }
                index1 = 1;
                ngbface = -1;
                if(ngbs[0][0]) {
                    ngbface = ngbfaces[0][0];
                    _simplices[ngbs[0][0]]->change_ngb(ngbface, sids[2],
                                                       index1);
                    _simplices[sids[2]]->add_ngb(index1, ngbs[0][0], ngbface);
                }

                index1 = _simplices[sids[0]]->get_index(other_vorgens[1]);
                ngbface = _simplices[sids[1]]->add_ngb_from_vorgen(
                        other_vorgens[2], sids[0], index1);
                _simplices[sids[0]]->add_ngb(index1, sids[1], ngbface);

                index1 = _simplices[sids[0]]->get_index(other_vorgens[0]);
                ngbface = _simplices[sids[2]]->add_ngb_from_vorgen(
                        other_vorgens[2], sids[0], index1);
                _simplices[sids[0]]->add_ngb(index1, sids[2], ngbface);

                index1 = _simplices[sids[1]]->get_index(other_vorgens[0]);
                ngbface = _simplices[sids[2]]->add_ngb_from_vorgen(
                        other_vorgens[1], sids[1], index1);
                _simplices[sids[1]]->add_ngb(index1, sids[2], ngbface);
            } else {
                // it does not: test whether it lies on an edge of the triangle
                if(intest_flag > 2) {
                    // it does: possible 4 to 4 flip. Find out if we have
                    // exactly 4 tetrahedra who share this edge
                    unsigned int s3 = tetrahedron->get_ngb(intest[2]);
                    unsigned int s4 = tetrahedron2->get_ngb_from_vorgen(
                            tetrahedron->vorgen(intest[2]));
                    if(s3 && s4 && s3 != s4 &&
                       _simplices[s3]->get_ngb_from_vorgen(vindex) == s4) {
                        // yes: 4 to 4 flip
                        // we first remove the simplices, because the flip
                        // changes _simplices!
                        removedra[rindex++] = _simplices[s3];
                        removedra[rindex++] = _simplices[s4];
                        four_to_four_flip(vindex, new_vorgen,
                                          tetrahedron->vorgen(intest[0]),
                                          tetrahedron->vorgen(intest[1]),
                                          simplex[0], simplex[1], s3, s4);
                    } else {
                        // not exactly 4 tetrahedra share the edge: nothing
                        // happens
                        // the failing tetrahedron will be replaced by treating
                        // another failing tetrahedron further on
                        return false;
                    }
                } else {
                    // it does not: possible 3 to 2 flip
                    for(unsigned k = 0; k < 3; k++) {
                        if((current + 1 + k) % 4 != intest[0] &&
                           (current + 1 + k) % 4 != intest[1]) {
                            unsigned int os =
                                    tetrahedron->get_ngb((current + 1 + k) % 4);
                            unsigned int os2 =
                                    tetrahedron2->get_ngb_from_vorgen(
                                            other_vorgens[k]);
                            // special case: the line cuts a face outside
                            // tetrahedron2. Nothing happens
                            // the failing tetrahedron will be replaced by
                            // treating another failing tetrahedron further on
                            if(!os || !os2 || os != os2) {
                                return false;
                            }

                            // 3 to 2 flip
                            Simplex* other_tetrahedron = _simplices[os];
                            removedra[rindex++] = other_tetrahedron;

                            // old version
                            unsigned int simplex2 = os;
                            unsigned int insidevorgen[2];
                            insidevorgen[0] = tetrahedron->vorgen(intest[0]);
                            insidevorgen[1] = tetrahedron->vorgen(intest[1]);

                            if(_free_size < 1000) {
                                _free_array[_free_size++] = simplex2;
                            }
                            Simplex* new_tetrahedra[2];
                            new_tetrahedra[0] = new Simplex(
                                    vindex, new_vorgen, other_vorgens[k],
                                    insidevorgen[0], _points);
                            new_tetrahedra[1] = new Simplex(
                                    vindex, new_vorgen, other_vorgens[k],
                                    insidevorgen[1], _points);

                            // removing simplex2 is totally unrelated to the new
                            // simplices, we just have to adjust the links in
                            // previous and next
                            // if simplex2 was not in the list, then it's
                            // next_check will be 0
                            if(other_tetrahedron->get_next_check()) {
                                // other_tetrahedron cannot be the first simplex
                                // on the stack, since that is already
                                // simplex[0] :)
                                if(other_tetrahedron->get_next_check() !=
                                   simplex2) {
                                    _simplices[other_tetrahedron
                                                       ->get_previous_check()]
                                            ->set_next_check(
                                                    other_tetrahedron
                                                            ->get_next_check());
                                    _simplices[other_tetrahedron
                                                       ->get_next_check()]
                                            ->set_previous_check(
                                                    other_tetrahedron
                                                            ->get_previous_check());
                                } else {
                                    _simplices[other_tetrahedron
                                                       ->get_previous_check()]
                                            ->set_next_check(
                                                    other_tetrahedron
                                                            ->get_previous_check());
                                }
                            }

                            _simplices[simplex[0]] = new_tetrahedra[0];
                            _simplices[simplex[1]] = new_tetrahedra[1];
                            _simplices[simplex2] = NULL;

                            // guaranteed self-reference
                            _simplices[simplex[0]]->set_previous_check(
                                    tetrahedron->get_previous_check());
                            _simplices[simplex[0]]->set_next_check(simplex[1]);
                            _simplices[simplex[1]]->set_previous_check(
                                    simplex[0]);
                            if(tetrahedron->get_next_check() != simplex[0]) {
                                _simplices[simplex[1]]->set_next_check(
                                        tetrahedron->get_next_check());
                                _simplices[tetrahedron->get_next_check()]
                                        ->set_previous_check(simplex[1]);
                            } else {
                                _simplices[simplex[1]]->set_next_check(
                                        simplex[1]);
                            }

                            unsigned int ngbs[6];
                            ngbs[0] = tetrahedron->get_ngb(intest[0]);
                            ngbs[1] = tetrahedron->get_ngb(intest[1]);
                            ngbs[2] = tetrahedron2->get_ngb_from_vorgen(
                                    insidevorgen[0]);
                            ngbs[3] = tetrahedron2->get_ngb_from_vorgen(
                                    insidevorgen[1]);
                            ngbs[4] = other_tetrahedron->get_ngb_from_vorgen(
                                    insidevorgen[0]);
                            ngbs[5] = other_tetrahedron->get_ngb_from_vorgen(
                                    insidevorgen[1]);

                            unsigned int ngbfaces[6];
                            ngbfaces[0] = tetrahedron->get_ngbface(intest[0]);
                            ngbfaces[1] = tetrahedron->get_ngbface(intest[1]);
                            ngbfaces[2] = tetrahedron2->get_ngbface_from_vorgen(
                                    insidevorgen[0]);
                            ngbfaces[3] = tetrahedron2->get_ngbface_from_vorgen(
                                    insidevorgen[1]);
                            ngbfaces[4] =
                                    other_tetrahedron->get_ngbface_from_vorgen(
                                            insidevorgen[0]);
                            ngbfaces[5] =
                                    other_tetrahedron->get_ngbface_from_vorgen(
                                            insidevorgen[1]);

                            // pnk0
                            unsigned int index1 = 0;
                            int ngbface = -1;
                            if(ngbs[3]) {
                                ngbface = ngbfaces[3];
                                _simplices[ngbs[3]]->change_ngb(
                                        ngbface, simplex[0], index1);
                                _simplices[simplex[0]]->add_ngb(index1, ngbs[3],
                                                                ngbface);
                            }
                            index1 = 1;
                            ngbface = -1;
                            if(ngbs[1]) {
                                ngbface = ngbfaces[1];
                                _simplices[ngbs[1]]->change_ngb(
                                        ngbface, simplex[0], index1);
                                _simplices[simplex[0]]->add_ngb(index1, ngbs[1],
                                                                ngbface);
                            }
                            index1 = _simplices[simplex[0]]->get_index(
                                    other_vorgens[k]);
                            ngbface = -1;
                            if(ngbs[5]) {
                                ngbface = ngbfaces[5];
                                _simplices[ngbs[5]]->change_ngb(
                                        ngbface, simplex[0], index1);
                                _simplices[simplex[0]]->add_ngb(index1, ngbs[5],
                                                                ngbface);
                            }
                            // pnk1
                            index1 = 0;
                            ngbface = -1;
                            if(ngbs[2]) {
                                ngbface = ngbfaces[2];
                                _simplices[ngbs[2]]->change_ngb(
                                        ngbface, simplex[1], index1);
                                _simplices[simplex[1]]->add_ngb(index1, ngbs[2],
                                                                ngbface);
                            }
                            index1 = 1;
                            ngbface = -1;
                            if(ngbs[0]) {
                                ngbface = ngbfaces[0];
                                _simplices[ngbs[0]]->change_ngb(
                                        ngbface, simplex[1], index1);
                                _simplices[simplex[1]]->add_ngb(index1, ngbs[0],
                                                                ngbface);
                            }
                            index1 = _simplices[simplex[1]]->get_index(
                                    other_vorgens[k]);
                            if(ngbs[4]) {
                                ngbface = ngbfaces[4];
                                _simplices[ngbs[4]]->change_ngb(
                                        ngbface, simplex[1], index1);
                                _simplices[simplex[1]]->add_ngb(index1, ngbs[4],
                                                                ngbface);
                            }

                            index1 = _simplices[simplex[0]]->get_index(
                                    insidevorgen[0]);
                            ngbface =
                                    _simplices[simplex[1]]->add_ngb_from_vorgen(
                                            insidevorgen[1], simplex[0],
                                            index1);
                            _simplices[simplex[0]]->add_ngb(index1, simplex[1],
                                                            ngbface);
                        }
                    }
                }
            }
            // garbage collection: remove all tetrahedra and unregister them
            // with their points
            for(unsigned k = 0; k < rindex; k++) {
                remove_tetrahedron(removedra[k]);
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
    return true;
}
#else
/**
 * @brief Check if the Simplex with the given index passes the geometrical
 * in-circle or in-sphere test
 *
 * Check a Simplex by subjecting the toppoint of the Simplex sharing the
 * face opposite to VorGen (this sentence is correct, reread it slowly until
 * it's clear) to the in sphere test. If the VorGen fails the test, the Simplex
 * has to be replaced by performing a 2 to 3 or 3 to 2 flip (most general case).
 * In some cases, no flip can be performed, while as a degenerate case, a 4 to 4
 * flip is also possible. In 2D, only a 2 to 2 flip exists.
 *
 * Because of the variety of possibilities, all simplices that have to be
 * removed are added to a local stack and are only removed when the method
 * successfully ends.
 *
 * @param index The index of the Simplex that has to be checked in the internal
 * _simplices vector
 * @param vindex The index of the VorGen that is being added in the internal
 * _points vector
 * @return True if the simplex with the given index has been changed during the
 * check and has to be rechecked.
 */
bool DelTess::check_simplex(unsigned int index, unsigned int vindex) {
    vector<Simplex*> removangles;
    Simplex* triangle = _simplices[index];
    unsigned int simplex[2] = {index, triangle->get_ngb_from_vorgen(vindex)};
    if(!simplex[1]) {
        return false;
    }
    // this can be achieved in a more elegant way, by making use of the
    // ngbfaces...
    Simplex* other_triangle = _simplices[triangle->get_ngb_from_vorgen(vindex)];
    unsigned int* vorgens = triangle->get_vorgens();
    unsigned int j = 0;
    while(vorgens[j] != vindex) {
        j++;
    }
    // j now contains the index of point in points
    unsigned int* other_vorgens = other_triangle->get_vorgens();
    unsigned int i = 0;
    while(other_vorgens[i] == vorgens[0] || other_vorgens[i] == vorgens[1] ||
          other_vorgens[i] == vorgens[2]) {
        i++;
    }
    // i now contains the index of the point of other_triangle that is not part
    // of triangle
    if(triangle->in_sphere(_points[other_vorgens[i]], _points)) {
        // flip side
        removangles.push_back(triangle);
        removangles.push_back(other_triangle);

        Simplex* new_triangles[2];
        // as usual, we might as well order the points in counterclockwise order
        // from the start
        new_triangles[0] = new Simplex(vindex, other_vorgens[i],
                                       other_vorgens[(i + 1) % 3]);
        new_triangles[1] = new Simplex(vindex, other_vorgens[(i + 2) % 3],
                                       other_vorgens[i]);
        for(int k = 0; k < 2; k++) {
            _simplices[simplex[k]] = new_triangles[k];
        }

        new_triangles[0]->set_previous_check(simplex[0]);
        new_triangles[0]->set_next_check(simplex[1]);
        new_triangles[1]->set_previous_check(simplex[0]);
        if(triangle->get_next_check() != simplex[0]) {
            new_triangles[1]->set_next_check(triangle->get_next_check());
        } else {
            new_triangles[1]->set_next_check(simplex[1]);
        }

        unsigned int ngb[2][3];
        unsigned int ngbfaces[2][3];
        for(int k = 0; k < 3; k++) {
            ngb[0][k] = triangle->get_ngb(k);
            ngb[1][k] = other_triangle->get_ngb(k);
            ngbfaces[0][k] = triangle->get_ngbface(k);
            ngbfaces[1][k] = other_triangle->get_ngbface(k);
        }
        // neighbours
        int ngbface = -1;
        int index1 = 0;
        if(ngb[1][(i + 2) % 3]) {
            ngbface = ngbfaces[1][(i + 2) % 3];
            _simplices[ngb[1][(i + 2) % 3]]->change_ngb(ngbface, simplex[0],
                                                        index1);
        }
        new_triangles[0]->add_ngb(index1, ngb[1][(i + 2) % 3], ngbface);

        ngbface = -1;
        index1 = new_triangles[1]->get_index(vindex);
        if(ngb[1][(i + 1) % 3]) {
            ngbface = ngbfaces[1][(i + 1) % 3];
            _simplices[ngb[1][(i + 1) % 3]]->change_ngb(ngbface, simplex[1],
                                                        index1);
        }
        new_triangles[1]->add_ngb(index1, ngb[1][(i + 1) % 3], ngbface);

        ngbface = -1;
        index1 = new_triangles[0]->get_index(other_vorgens[i]);
        if(ngb[0][(j + 1) % 3]) {
            ngbface = ngbfaces[0][(j + 1) % 3];
            _simplices[ngb[0][(j + 1) % 3]]->change_ngb(ngbface, simplex[0],
                                                        index1);
        }
        new_triangles[0]->add_ngb(index1, ngb[0][(j + 1) % 3], ngbface);

        ngbface = -1;
        index1 = new_triangles[1]->get_index(other_vorgens[i]);
        if(ngb[0][(j + 2) % 3]) {
            ngbface = ngbfaces[0][(j + 2) % 3];
            _simplices[ngb[0][(j + 2) % 3]]->change_ngb(ngbface, simplex[1],
                                                        index1);
        }
        new_triangles[1]->add_ngb(index1, ngb[0][(j + 2) % 3], ngbface);

        index1 = new_triangles[1]->get_index(other_vorgens[(i + 2) % 3]);
        ngbface = new_triangles[0]->add_ngb_from_vorgen(
                other_vorgens[(i + 1) % 3], simplex[1], index1);
        new_triangles[1]->add_ngb(index1, simplex[0], ngbface);
    } else {
        return false;
    }

    for(unsigned k = 0; k < removangles.size(); k++) {
        remove_tetrahedron(removangles[k]);
    }
    return true;
}
#endif

#if ndim_ == 3
/**
 * @brief Get the common tetrahedron of the 3 given points
 *
 * Get the 2 tetrahedra that share the face formed by the 3 given points.
 * This function is quite heavy and should not be used very often, as it
 * consists of a triple for-loop with undetermined size (a point has order 10
 * tetrahedra)
 *
 * @param point1,point2,point3 Three points for which we search a common
 * tetrahedron
 * @param sids The indices of the six possible Simplices in the internal
 * _simplices vector
 * @param common_tetrahedron Array in which to store the indices of the two
 * common Simplices
 */
void DelTess::get_common_tetrahedron(unsigned int point1, unsigned int point2,
                                     unsigned int point3, unsigned int* sids,
                                     unsigned int* common_tetrahedron) {
    unsigned int position = 0;
    for(unsigned int i = 6; i--;) {
        if(_simplices[sids[i]] == NULL) {
            continue;
        }
        unsigned int* vorgens = _simplices[sids[i]]->get_vorgens();
        if(vorgens[0] == point1 || vorgens[1] == point1 ||
           vorgens[2] == point1 || vorgens[3] == point1) {
            if(vorgens[0] == point2 || vorgens[1] == point2 ||
               vorgens[2] == point2 || vorgens[3] == point2) {
                if(vorgens[0] == point3 || vorgens[1] == point3 ||
                   vorgens[2] == point3 || vorgens[3] == point3) {
                    common_tetrahedron[position++] = sids[i];
                    if(position > 1) {
                        return;
                    }
                }
            }
        }
    }
}
#endif

#if ndim_ == 3
unsigned int DelTess::find_simplex(VorGen* point, unsigned int* found) {
    unsigned int result = 0;
    unsigned int simplex = _lastindex;
    while(_simplices[simplex] == NULL) {
        simplex--;
    }
    if(!simplex) {
        simplex++;
        while(_simplices[simplex] == NULL) {
            simplex++;
        }
    }
    int test = _simplices[simplex]->inside(point, _points);
    unsigned int previous_simplex = 0;
    unsigned int numit = 0;
    while(test == 0 && numit < 1000) {
        double com[3];
        com[0] = 0;
        com[1] = 0;
        com[2] = 0;
        unsigned int* vorgens = _simplices[simplex]->get_vorgens();
        VorGen* points[4] = {_points[vorgens[0]], _points[vorgens[1]],
                             _points[vorgens[2]], _points[vorgens[3]]};
        for(int i = 0; i < 4; i++) {
            Vec p12 = points[i]->get_p12();
            com[0] += p12.x();
            com[1] += p12.y();
            com[2] += p12.z();
        }
        com[0] *= 0.25;
        com[1] *= 0.25;
        com[2] *= 0.25;
        VorGen center = VorGen(com[0], com[1], com[2]);
        center.set_p12(Vec(com[0], com[1], com[2]));
        Line connection = Line(&center, point);
        unsigned int tmpid = 0;
        double abce = ExactArithmetic::orient3d(
                points[0]->get_p12(), points[1]->get_p12(),
                points[2]->get_p12(), point->get_p12());
        if(abce > 0.) {
            tmpid = _simplices[simplex]->get_ngb(3);
        } else {
            double abde = ExactArithmetic::orient3d(
                    points[0]->get_p12(), points[1]->get_p12(),
                    points[3]->get_p12(), point->get_p12());
            if(abde < 0.) {
                tmpid = _simplices[simplex]->get_ngb(2);
            } else {
                double acde = ExactArithmetic::orient3d(
                        points[0]->get_p12(), points[2]->get_p12(),
                        points[3]->get_p12(), point->get_p12());
                if(acde > 0.) {
                    tmpid = _simplices[simplex]->get_ngb(1);
                } else {
                    double bcde = ExactArithmetic::orient3d(
                            points[1]->get_p12(), points[2]->get_p12(),
                            points[3]->get_p12(), point->get_p12());
                    if(bcde < 0.) {
                        tmpid = _simplices[simplex]->get_ngb(0);
                    } else {
                        cerr << "Point is not inside, but exact tests say "
                                "it is... Confused..."
                             << endl;
                        my_exit();
                    }
                }
            }
        }
        if(tmpid != previous_simplex) {
            previous_simplex = simplex;
            simplex = tmpid;
        } else {
            // special case: we are trapped in an eternal loop
            cerr << "Eternal loop encountered in tetrahedron search. "
                    "Rerouting..."
                 << endl;
            previous_simplex = simplex;
            simplex = _simplices[simplex]->get_ngb(rand() % 4);
        }
        test = _simplices[simplex]->inside(point, _points);
        numit++;
    }
    if(numit == 1000) {
        cerr << "Error! Maximum number of iterations (1000) reached in simplex "
                "search!"
             << endl;
        cerr << "Try increasing the Voronoi tolerance in the .ini file..."
             << endl;
        my_exit();
    }
    found[result++] = simplex;
    // test for degenerate cases
    if(test > 1) {
        if(test > 2) {
            // point is on edge of current tesselation
            // we have to rotate around the (unknown) common axis
            unsigned int* vorgens = _simplices[simplex]->get_vorgens();
            VorGen* points[4] = {_points[vorgens[0]], _points[vorgens[1]],
                                 _points[vorgens[2]], _points[vorgens[3]]};
            double tests[6];
            tests[0] = normcross(points[0], points[1], point);  // ab
            tests[1] = normcross(points[0], points[2], point);  // ac
            tests[2] = normcross(points[0], points[3], point);  // ad
            tests[3] = normcross(points[1], points[2], point);  // bc
            tests[4] = normcross(points[1], points[3], point);  // bd
            tests[5] = normcross(points[2], points[3], point);  // cd
            // find minimum
            int min = 0;
            for(int i = 1; i < 6; i++) {
                if(tests[i] < tests[min]) {
                    min = i;
                }
            }
            // convert min-value to points on axis and not on axis
            int non_axis[2], axis[2];
            if(min > 2) {
                non_axis[0] = 0;
            } else if(min > 0) {
                non_axis[0] = 1;
            } else {
                non_axis[0] = 2;
            }
            if(min == 5) {
                non_axis[1] = 1;
            } else if(min == 2 || min == 4) {
                non_axis[1] = 2;
            } else {
                non_axis[1] = 3;
            }
            if(min < 3) {
                axis[0] = 0;
            } else if(min < 5) {
                axis[0] = 1;
            } else {
                axis[0] = 2;
            }
            if(min == 0) {
                axis[1] = 1;
            } else if(min == 1 || min == 3) {
                axis[1] = 2;
            } else {
                axis[1] = 3;
            }
            // walk around the axis and add tetrahedra on the fly (code similar
            // to voronoicell construction)
            // it could be that this code is too complicated and slow, but as it
            // is seldom processed, this is not of our present concern
            unsigned int vaxis1 = vorgens[axis[0]];
            unsigned int vaxis2 = vorgens[axis[1]];
            unsigned int vstart = vorgens[non_axis[0]];
            unsigned int vother = vorgens[non_axis[1]];
            unsigned int next_simplex =
                    _simplices[simplex]->get_ngb(non_axis[0]);
            unsigned int* next_vorgens =
                    _simplices[next_simplex]->get_vorgens();
            int next_index = 0;
            while(next_index < 4 && (next_vorgens[next_index] == vaxis1 ||
                                     next_vorgens[next_index] == vaxis2 ||
                                     next_vorgens[next_index] == vother)) {
                next_index++;
            }
            unsigned int vnext = next_vorgens[next_index];
            vstart = vother;
            int counter = 0;
            while(next_simplex != found[0] &&
                  counter < FIND_SIMPLEX_BUFFER_SIZE) {
                found[result++] = next_simplex;
                next_simplex =
                        _simplices[next_simplex]->get_ngb_from_vorgen(vstart);
                next_vorgens = _simplices[next_simplex]->get_vorgens();
                next_index = 0;
                while(next_index < 4 && (next_vorgens[next_index] == vaxis1 ||
                                         next_vorgens[next_index] == vaxis2 ||
                                         next_vorgens[next_index] == vnext)) {
                    next_index++;
                }
                vstart = vnext;
                vnext = next_vorgens[next_index];
                counter++;
            }
            if(counter == FIND_SIMPLEX_BUFFER_SIZE) {
                // should not happen
                cout << "encountered an axis with more than"
                     << FIND_SIMPLEX_BUFFER_SIZE
                     << " common "
                        "tetrahedra. This is highly unlikely. We better stop."
                     << endl;
                my_exit();
            }
        } else {
            // point is on face of current tesselation. Find other tetrahedron.
            for(int i = 0; i < 4; i++) {
                unsigned int ngbsimplex = _simplices[simplex]->get_ngb(i);
                if(ngbsimplex &&
                   _simplices[ngbsimplex]->inside(point, _points) > 0) {
                    found[result++] = ngbsimplex;
                }
            }
            if(result < 2) {
                // should not happen
                cout << "Error: point is on face, but only one tetrahedron "
                        "found for face..."
                     << endl;
                my_exit();
            }
        }
    }
    return result;
}
#else
/**
 * @brief Find out which Simplex contains the given VorGen
 *
 * The core function determining the speed of the algorithm: find out in which
 * Simplex the given VorGen resides. To speed up the process, an initial guess
 * can be provided. If the points are added in space filling order and the
 * initial guess is set to the last simplex added, this usually works very
 * fast.
 *
 * If a Simplex fails the inside test, we determine which face of it is cut
 * by the line connecting its centroid with the VorGen tested and move on to its
 * neighbour at that side. This allows us to efficiently locate the correct
 * Simplex. Due to numerical inaccuracies, an endless loop can occur. This is
 * prevented by using exact geometric tests if necessary. In very exceptional
 * cases an endless loop still occurs. We therefore keep track of tested
 * simplices and if an endless loop is detected, we randomly reset the current
 * Simplex to one of its neighbours. This usually fixes the problem.
 *
 * The newly inserted VorGen can also lie on a face or and edge of the current
 * tesselation. This is detected here and we make sure that the resulting vector
 * contains all simplices that fulfill the inside condition.
 *
 * @param point The VorGen that is being added to the tesselation
 * @param simplex Array in which to store the simplex or simplices containing
 * the given VorGen
 * @return The number of simplices found (usually one, two means the VorGen is
 * on a face of the tesselation, three that it is on an edge; in 2D only the
 * edge is possible)
 */
unsigned int DelTess::find_simplex(VorGen* point, unsigned int* simplex) {
    unsigned int result = 0;
    unsigned int triangle = _simplices.size() - 1;
    int test = _simplices[triangle]->inside(point, _points);
    while(!test) {
        unsigned int* vorgens = _simplices[triangle]->get_vorgens();
        VorGen* points[3] = {_points[vorgens[0]], _points[vorgens[1]],
                             _points[vorgens[2]]};
        double ctr[2];
        _simplices[triangle]->get_centroid(ctr, _points);
        Vec centroid1(ctr[0], ctr[1]);
        Vec centroid2 = get_p12(centroid1);
        Vec n = point->get_p12() - centroid2;
        int i = 2;
        bool intertest = false;
        while(i > -1 && !intertest) {
            Vec a = points[(i + 1) % 3]->get_p12();
            Vec b = points[(i + 2) % 3]->get_p12();
            double denom = n[0] * (b[1] - a[1]) + n[1] * (a[0] - b[0]);
            if(denom != 0) {
                denom = 1. / denom;
                double h = (n[0] * centroid2[1] - n[1] * centroid2[0] +
                            n[1] * a[0] - n[0] * a[1]) *
                           denom;
                if(h > _tolerance && h < 1 - _tolerance) {
                    double p = (centroid2[1] * (b[0] - a[0]) +
                                a[1] * (centroid2[0] - b[0]) +
                                b[1] * (a[0] - centroid2[0])) *
                               denom;
                    if(p > _tolerance) {
                        intertest = true;
                    }
                }
            }
            i--;
        }
        if(intertest) {
            // i now contains the index of the point of which the associated
            // neighbour is the next triangle
            triangle = _simplices[triangle]->get_ngb(i + 1);
        } else {
            // need a more exact method
            Vec a = points[0]->get_p12();
            Vec b = points[1]->get_p12();
            Vec c = points[2]->get_p12();
            Vec d = point->get_p12();
            double abd = ExactArithmetic::orient2d(a, b, d);
            if(abd <= 0.) {
                triangle = _simplices[triangle]->get_ngb(2);
            } else {
                double acd = ExactArithmetic::orient2d(a, c, d);
                if(acd >= 0.) {
                    triangle = _simplices[triangle]->get_ngb(1);
                } else {
                    // no need to check the third face, because if it is not
                    // at that side, the point should be in the triangle, which
                    // it is not
                    triangle = _simplices[triangle]->get_ngb(2);
                }
            }
        }
        test = _simplices[triangle]->inside(point, _points);
    }
    simplex[result++] = triangle;
    // test for degenerate case...
    if(test > 1) {
        simplex[result++] = _simplices[triangle]->get_ngb(test / 2 - 1);
    }
    return result;
}
#endif

#if ndim_ == 3
/**
 * @brief Calculate the norm of the cross product of (a-point) and (b-point)
 *
 * This should be 0 if a, b and point are colinear
 *
 * @param a,b,point Points as explained above
 * @return The cross product of (a-point) and (b-point)
 */
double DelTess::normcross(VorGen* a, VorGen* b, VorGen* point) {
    double ap[3], bp[3];
    Vec a12 = a->get_p12();
    Vec b12 = b->get_p12();
    Vec p12 = point->get_p12();
    for(int i = 0; i < 3; i++) {
        ap[i] = a12[i] - p12[i];
        bp[i] = b12[i] - p12[i];
    }
    double v[3];
    v[0] = ap[1] * bp[2] - ap[2] * bp[1];
    v[1] = ap[2] * bp[0] - ap[0] * bp[2];
    v[2] = ap[0] * bp[1] - ap[1] * bp[0];
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}
#endif

/**
 * @brief Remove the given Simplex
 *
 * This boils down to adding the Simplex to the remove stack for later removal.
 *
 * @param tetrahedron The Simplex to remove
 */
void DelTess::remove_tetrahedron(Simplex* tetrahedron) {
    _remove_stack.push_back(tetrahedron);
}

/**
 * @brief Output the tesselation (in ascii) to the given stream
 *
 * The points are outputted in the format \verbatim
p:\t(c1,c2,c3)\n\endverbatim
 * (with the ci's doubles)
 * The tetrahedra are outputted as \verbatim
t:\t(c1,c2,c3)\t(d1,d2,d3)\t(e1,e2,e3)\t(f1,f2,f3)\n\endverbatim
 * This last part does not work for the moment being
 *
 * @param stream An ostream to write to
 */
void DelTess::output_tesselation(ostream& stream) {
    for(vector<VorGen*>::iterator it = _points.begin(); it < _points.end();
        it++) {
        VorGen* point = *it;
        point->print(stream);
        stream << "\n\n";
    }
    stream << flush;
    for(unsigned int i = 1; i < _simplices.size(); i++) {
        if(_simplices[i] != NULL) {
            _simplices[i]->print(stream, _points);
        }
    }
}

/**
 * @brief Get the points of the tesselation.
 *
 * @return A vector containing all the points of the tesselation
 */
vector<VorGen*> DelTess::get_points() { return _points; }

/**
  * @brief  Add mirror particles to the tesselation to impose boundary
  * conditions and to take into account particles residing on other processes.
  *
  * We really need complete cells for all active particles (the ones that are
  * already in the tesselation) and all neighbouring particles of these active
  * particles. We therefore perform a double loop. We set the id of all active
  * points (points with an active particle) to 2 and then loop over all points
  * in the tesselation. We add neighbouring points and give them an id of 1. The
  * id of the original point is set to 0, unless if we did not find enough
  * neighbours to guarantee a complete cell, then it is kept at 2 and the loop
  * is repeated for that point. After the first loop has finished, points have
  * either an id of 0, in which case they correspond to active cells, or an id 1
  * which corresponds to inactive neighbours. To guarantee completeness for the
  * latter, we again perform a loop, this time over all points with id 1. We
  * again add neighbours with an id 1 and set the id of the original points to 0
  * if the cell is guaranteed to be complete. At the end, we are left with
  * points with an id 0 that are complete, and points with id 1 that are
  * incomplete. To distinguish between active and non-active cells, we then set
  * the id of all points corresponding to complete cells to 2. This way, we
  * have 2 possibilities:
  *  - id 1: a ghost/mirror that is incomplete
  *  - id 2: a complete cell
  *
  * There are only 3/4 (depending on number of dimensions) points with id 0: the
  * vertices of the large all encompassing simplex we set up at the start.
  *
  * @param parttree A Tree containing all particles that need to be considered
  */
void DelTess::add_mirrors(Tree& parttree) {
    LOGS("Start adding mirror points");
    int worldsize = MPIGlobal::size;
    int worldrank = MPIGlobal::rank;
    unsigned int counter = 0;
    _exports.resize(worldsize);
    _exportcopies.resize(worldsize);
    _ghosts.resize(worldsize);
    vector<bool> dummyflags(worldsize);

// points may be added during the loop (not deleted), we only want to loop
// over the old points
#if ndim_ == 2
    unsigned int pointssize = _points.size() - 3;
    unsigned int startpoint = 0;
#else
    unsigned int pointssize = _points.size() - 4;
    unsigned int startpoint = 0;
#endif
    // activate all particles, they will get surrounded by complete cells
    unsigned int redo = 0;
    bool initial_radius = false;
    for(unsigned int i = startpoint; i < pointssize; i++) {
        if(_points[i]) {
            _points[i]->set_id(2);
            initial_radius |= !_points[i]->get_particle()->get_max_radius();
            redo++;
        }
    }

    if(initial_radius) {
        // first step, need to initialize the search radius
        for(unsigned int i = 1; i < _simplices.size(); i++) {
            if(_simplices[i] != NULL) {
                VorGen* sp = _simplices[i]->get_special_point(_points);
                unsigned int* vorgens = _simplices[i]->get_vorgens();
#if ndim_ == 3
                VorGen* points[4] = {_points[vorgens[0]], _points[vorgens[1]],
                                     _points[vorgens[2]], _points[vorgens[3]]};
#else
                VorGen* points[3] = {_points[vorgens[0]], _points[vorgens[1]],
                                     _points[vorgens[2]]};
#endif
                double r =
                        (sp->get_position() - points[0]->get_position()).norm();
                for(unsigned int j = ndim_ + 1; j--;) {
                    if(points[j]->get_particle() != NULL &&
                       points[j]->get_id() != 1) {
                        if(points[j]->get_particle()->get_max_radius()) {
                            points[j]->get_particle()->set_max_radius(min(
                                    points[j]->get_particle()->get_max_radius(),
                                    2. * r));
                        } else {
                            points[j]->get_particle()->set_max_radius(2. * r);
                        }
                    }
                }
            }
        }
    }

    vector<unsigned int> other_indices(worldsize, 0);
    unsigned int lastindex = 1;
    if(worldsize > 1) {
        unsigned int redo_glob;
        MyMPI_Allreduce(&redo, &redo_glob, 1, MPI_UNSIGNED, MPI_SUM);
        redo = redo_glob;
    }
    // first add the particles whose cells have to be constructed entirely
    while(redo) {
        vector<vector<Vec> > exportlist(worldsize);
        vector<vector<Vec> > exportlist_outside(worldsize);
        vector<vector<double> > exportradii(worldsize);
        vector<vector<double> > exportradii_outside(worldsize);
        for(unsigned int i = startpoint; i < pointssize; i++) {
            if(!_points[i]) {
                continue;
            }
            if(_points[i]->get_id() == 2) {
                vector<bool> exportflag(worldsize, false);
                double rad = _points[i]->get_particle()->get_max_radius();
                Vec coords = _points[i]->get_position();
                if(_periodic) {
                    vector<VorGen*> ngbs;
                    ngbs.reserve(100);
                    parttree.get_neighbours_outside(coords, rad, ngbs,
                                                    exportflag);
                    for(unsigned int j = ngbs.size(); j--;) {
                        if(ngbs[j]->get_particle()->get_vorgen() == NULL) {
                            add_particle(
                                    ngbs[j]->get_particle(),
                                    ngbs[j]->get_particle()->get_local_id());
                            counter++;
                        }
                        _mirrors.push_back(ngbs[j]);
                        _mirrors.back()->set_id(1);
                        _mirrors.back()->set_original(
                                ngbs[j]->get_particle()->get_local_id());
                        _points.push_back(ngbs[j]);
                        add_point(_points.size() - 1);
                        counter++;
                    }
                } else {
                    vector<VorGen*> ngbs;
                    ngbs.reserve(100);
                    parttree.get_neighbours_outside(coords, rad, ngbs,
                                                    exportflag);
                    for(unsigned int j = ngbs.size(); j--;) {
                        ngbs[j]->set_original(
                                ngbs[j]->get_particle()->get_local_id());
                        ngbs[j]->set_particle(NULL);
                        _mirrors.push_back(ngbs[j]);
                        _mirrors.back()->set_id(1);
                        _points.push_back(ngbs[j]);
                        add_point(_points.size() - 1);
                        counter++;
                    }
                }
                for(unsigned int j = worldsize; j--;) {
                    if(exportflag[j]) {
                        exportlist_outside[j].push_back(coords);
                        exportradii_outside[j].push_back(rad);
                        exportflag[j] = false;
                    }
                }
                NgbSearch ngbsearch(coords, rad, exportflag, 100);
                parttree.walk_tree(ngbsearch);
                vector<GasParticle*> inside_ngbs = ngbsearch.get_ngbs();
                for(unsigned int j = 0; j < inside_ngbs.size(); j++) {
                    if(inside_ngbs[j] != NULL &&
                       inside_ngbs[j]->get_vorgen() == NULL) {
                        add_particle(inside_ngbs[j],
                                     inside_ngbs[j]->get_local_id());
                        counter++;
                    }
                }
                for(unsigned int j = worldsize; j--;) {
                    if(exportflag[j]) {
                        exportlist[j].push_back(coords);
                        exportradii[j].push_back(rad);
                    }
                }
                _points[i]->set_id(0);
            }
        }

        // communicate
        vector<vector<Vec> > inlist(worldsize);
        vector<vector<double> > inradii(worldsize);
        vector<vector<Vec> > inlist_outside(worldsize);
        vector<vector<double> > inradii_outside(worldsize);
        if(MPIGlobal::size > 1) {
            std::vector<MPI_Request> reqs((MPIGlobal::size - 1) * 2,
                                          MPI_REQUEST_NULL);
            vector<unsigned int> buffers((MPIGlobal::size - 1) * 2);
            unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
            for(unsigned int i = 0; i < buffers.size(); i++) {
                buffers[i] = i * bufsize;
            }

            unsigned int maxsize = bufsize / ((ndim_ + 1) * sizeof(double));
            // make sure we have some space left at the end of the buffer to
            // signal more to come
            if(!(bufsize % ((ndim_ + 1) * sizeof(double)))) {
                maxsize--;
            }

            vector<int> allflag((MPIGlobal::size - 1) * 2);
            vector<unsigned int> numsend((MPIGlobal::size - 1) * 2);
            vector<unsigned int> sendsizes(2 * (MPIGlobal::size - 1));

            int allflagsum = 0;
            for(int i = 0; i < MPIGlobal::size - 1; i++) {
                // initialize send sizes
                sendsizes[2 * i] =
                        exportlist[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                .size();
                sendsizes[2 * i + 1] =
                        exportlist_outside[(MPIGlobal::rank + 1 + i) %
                                           MPIGlobal::size]
                                .size();

                int send_pos = 0;
                for(unsigned int si = 0;
                    si < std::min(maxsize, sendsizes[2 * i]); si++) {
                    MyMPI_Pack(&exportlist[(MPIGlobal::rank + 1 + i) %
                                           MPIGlobal::size][si][0],
                               ndim_, MPI_DOUBLE,
                               &MPIGlobal::sendbuffer[buffers[2 * i]], bufsize,
                               &send_pos);
                    MyMPI_Pack(&exportradii[(MPIGlobal::rank + 1 + i) %
                                            MPIGlobal::size][si],
                               1, MPI_DOUBLE,
                               &MPIGlobal::sendbuffer[buffers[2 * i]], bufsize,
                               &send_pos);
                }
                allflag[2 * i] = (sendsizes[2 * i] <= maxsize);
                // if allflagsum equals MPIGlobal::size-1, all data has been
                // sent
                allflagsum += allflag[2 * i];
                numsend[2 * i] = 1;
                // add continuation signal
                MyMPI_Pack(&allflag[2 * i], 1, MPI_INT,
                           &MPIGlobal::sendbuffer[buffers[2 * i]], bufsize,
                           &send_pos);

                MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[2 * i]], send_pos,
                            MPI_PACKED,
                            (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0,
                            &reqs[2 * i]);

                send_pos = 0;
                for(unsigned int si = 0;
                    si < std::min(maxsize, sendsizes[2 * i + 1]); si++) {
                    MyMPI_Pack(&exportlist_outside[(MPIGlobal::rank + 1 + i) %
                                                   MPIGlobal::size][si][0],
                               ndim_, MPI_DOUBLE,
                               &MPIGlobal::sendbuffer[buffers[2 * i + 1]],
                               bufsize, &send_pos);
                    MyMPI_Pack(&exportradii_outside[(MPIGlobal::rank + 1 + i) %
                                                    MPIGlobal::size][si],
                               1, MPI_DOUBLE,
                               &MPIGlobal::sendbuffer[buffers[2 * i + 1]],
                               bufsize, &send_pos);
                }
                allflag[2 * i + 1] = (sendsizes[2 * i + 1] <= maxsize);
                // if allflagsum equals MPIGlobal::size-1, all data has been
                // sent
                allflagsum += allflag[2 * i + 1];
                numsend[2 * i + 1] = 1;
                // add continuation signal
                MyMPI_Pack(&allflag[2 * i + 1], 1, MPI_INT,
                           &MPIGlobal::sendbuffer[buffers[2 * i + 1]], bufsize,
                           &send_pos);

                MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[2 * i + 1]],
                            send_pos, MPI_PACKED,
                            (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 1,
                            &reqs[2 * i + 1]);
            }

            for(int i = 0; i < 2 * (MPIGlobal::size - 1); i++) {
                MPI_Status status;
                for(int j = 0; j < MPIGlobal::size - 1; j++) {
                    if(!allflag[2 * j]) {
                        int flag;
                        MyMPI_Test(&reqs[2 * j], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[2 * j] * maxsize;
                                si < std::min((numsend[2 * j] + 1) * maxsize,
                                              sendsizes[2 * j]);
                                si++) {
                                MyMPI_Pack(
                                        &exportlist[(MPIGlobal::rank + 1 + j) %
                                                    MPIGlobal::size][si][0],
                                        ndim_, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j]],
                                        bufsize, &send_pos);
                                MyMPI_Pack(
                                        &exportradii[(MPIGlobal::rank + 1 + j) %
                                                     MPIGlobal::size][si],
                                        1, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j]],
                                        bufsize, &send_pos);
                            }
                            allflag[2 * j] = (sendsizes[2 * j] <=
                                              (numsend[2 * j] + 1) * maxsize);
                            allflagsum += allflag[2 * j];
                            numsend[2 * j]++;
                            // add continuation signal
                            MyMPI_Pack(&allflag[2 * j], 1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[2 * j]],
                                       bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[2 * j]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    0, &reqs[2 * j]);
                        }
                    }
                    if(!allflag[2 * j + 1]) {
                        int flag;
                        MyMPI_Test(&reqs[2 * j + 1], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[2 * j + 1] * maxsize;
                                si <
                                std::min((numsend[2 * j + 1] + 1) * maxsize,
                                         sendsizes[2 * j + 1]);
                                si++) {
                                MyMPI_Pack(
                                        &exportlist_outside[(MPIGlobal::rank +
                                                             1 + j) %
                                                            MPIGlobal::size][si]
                                                           [0],
                                        ndim_, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j +
                                                                       1]],
                                        bufsize, &send_pos);
                                MyMPI_Pack(
                                        &exportradii_outside[(MPIGlobal::rank +
                                                              1 + j) %
                                                             MPIGlobal::size]
                                                            [si],
                                        1, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j +
                                                                       1]],
                                        bufsize, &send_pos);
                            }
                            allflag[2 * j + 1] =
                                    (sendsizes[2 * j + 1] <=
                                     (numsend[2 * j + 1] + 1) * maxsize);
                            allflagsum += allflag[2 * j + 1];
                            numsend[2 * j + 1]++;
                            // add continuation signal
                            MyMPI_Pack(
                                    &allflag[2 * j + 1], 1, MPI_INT,
                                    &MPIGlobal::sendbuffer[buffers[2 * j + 1]],
                                    bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[2 * j + 1]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    1, &reqs[2 * j + 1]);
                        }
                    }
                }

                int index;
                int tag;
                int recv_pos;

                // wait for a message from any source and with any tag to arrive
                MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
                // query the MPI_Status object to retrieve tag, source and size
                tag = status.MPI_TAG;
                index = status.MPI_SOURCE;
                int nelements;
                MyMPI_Get_count(&status, MPI_PACKED, &nelements);
                // select an index in the buffers to use for receiving and
                // sending
                int freebuffer;
                if(index == MPIGlobal::size - 1) {
                    freebuffer = 2 * MPIGlobal::rank;
                } else {
                    freebuffer = 2 * index;
                }
                if(tag == 0) {
                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                               nelements, MPI_PACKED, index, tag, &status);
                    recv_pos = 0;
                    while(recv_pos < nelements - 4) {
                        Vec inlist_el;
                        double inradii_el;
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &inlist_el[0], ndim_,
                                MPI_DOUBLE);
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &inradii_el, 1,
                                MPI_DOUBLE);
                        inlist[index].push_back(inlist_el);
                        inradii[index].push_back(inradii_el);
                    }
                    int flag;
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &flag, 1, MPI_INT);
                    if(!flag) {
                        i--;
                    }
                }
                if(tag == 1) {
                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                               nelements, MPI_PACKED, index, tag, &status);
                    recv_pos = 0;
                    while(recv_pos < nelements - 4) {
                        Vec inlist_outside_el;
                        double inradii_outside_el;
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &inlist_outside_el[0],
                                ndim_, MPI_DOUBLE);
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &inradii_outside_el, 1,
                                MPI_DOUBLE);
                        inlist_outside[index].push_back(inlist_outside_el);
                        inradii_outside[index].push_back(inradii_outside_el);
                    }
                    if(inlist_outside[index].size() > 0 &&
                       !other_indices[index]) {
                        other_indices[index] = lastindex;
                        lastindex++;
                    }
                    int flag;
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &flag, 1, MPI_INT);
                    if(!flag) {
                        i--;
                    }
                }
            }

            while(allflagsum < 2 * (MPIGlobal::size - 1)) {
                MPI_Status status;
                for(int j = 0; j < MPIGlobal::size - 1; j++) {
                    if(!allflag[2 * j]) {
                        int flag;
                        MyMPI_Test(&reqs[2 * j], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[2 * j] * maxsize;
                                si < std::min((numsend[2 * j] + 1) * maxsize,
                                              sendsizes[2 * j]);
                                si++) {
                                MyMPI_Pack(
                                        &exportlist[(MPIGlobal::rank + 1 + j) %
                                                    MPIGlobal::size][si],
                                        ndim_, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j]],
                                        bufsize, &send_pos);
                                MyMPI_Pack(
                                        &exportradii[(MPIGlobal::rank + 1 + j) %
                                                     MPIGlobal::size][si],
                                        1, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j]],
                                        bufsize, &send_pos);
                            }
                            allflag[2 * j] = (sendsizes[2 * j] <=
                                              (numsend[2 * j] + 1) * maxsize);
                            allflagsum += allflag[2 * j];
                            numsend[2 * j]++;
                            // add continuation signal
                            MyMPI_Pack(&allflag[2 * j], 1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[2 * j]],
                                       bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[2 * j]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    0, &reqs[2 * j]);
                        }
                    }
                    if(!allflag[2 * j + 1]) {
                        int flag;
                        MyMPI_Test(&reqs[2 * j + 1], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[2 * j + 1] * maxsize;
                                si <
                                std::min((numsend[2 * j + 1] + 1) * maxsize,
                                         sendsizes[2 * j + 1]);
                                si++) {
                                MyMPI_Pack(
                                        &exportlist_outside[(MPIGlobal::rank +
                                                             1 + j) %
                                                            MPIGlobal::size]
                                                           [si],
                                        ndim_, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j +
                                                                       1]],
                                        bufsize, &send_pos);
                                MyMPI_Pack(
                                        &exportradii_outside[(MPIGlobal::rank +
                                                              1 + j) %
                                                             MPIGlobal::size]
                                                            [si],
                                        1, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j +
                                                                       1]],
                                        bufsize, &send_pos);
                            }
                            allflag[2 * j + 1] =
                                    (sendsizes[2 * j + 1] <=
                                     (numsend[2 * j + 1] + 1) * maxsize);
                            allflagsum += allflag[2 * j + 1];
                            numsend[2 * j + 1]++;
                            // add continuation signal
                            MyMPI_Pack(
                                    &allflag[2 * j + 1], 1, MPI_INT,
                                    &MPIGlobal::sendbuffer[buffers[2 * j + 1]],
                                    bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[2 * j + 1]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    1, &reqs[2 * j + 1]);
                        }
                    }
                }
            }

            // if we do not put a barrier here, the number of loops above might
            // get screwed up by interfering messages from below...
            vector<MPI_Status> status((MPIGlobal::size - 1) * 2);
            MyMPI_Waitall((MPIGlobal::size - 1) * 2, &reqs[0], &status[0]);
            MyMPI_Barrier();
        }

        // search for imported coordinates
        vector<vector<GasParticle*> > outlist(worldsize);
        for(int i = worldsize; i--;) {
            if(worldrank != i) {
                for(int j = inlist_outside[i].size(); j--;) {
                    vector<VorGen*> import_ngbs;
                    import_ngbs.reserve(100);
                    parttree.get_neighbours_outside(
                            inlist_outside[i][j], inradii_outside[i][j],
                            import_ngbs, dummyflags, other_indices[i]);
                    for(unsigned int k = import_ngbs.size(); k--;) {
                        if(!import_ngbs[k]->get_particle()->is_pseudo()) {
                            if(import_ngbs[k]->get_particle()->get_vorgen() ==
                               NULL) {
                                add_particle(import_ngbs[k]->get_particle(),
                                             import_ngbs[k]
                                                     ->get_particle()
                                                     ->get_local_id());
                                counter++;
                            }
                            _exportcopies[i].push_back(new GasParticle(
                                    *import_ngbs[k]->get_particle()));
                            _exportcopies[i].back()->set_x(import_ngbs[k]->x());
                            _exportcopies[i].back()->set_y(import_ngbs[k]->y());
#if ndim_ == 3
                            _exportcopies[i].back()->set_z(import_ngbs[k]->z());
#endif
                            outlist[i].push_back(_exportcopies[i].back());
                            if(_periodic) {
                                _exports[i].push_back(_exportcopies[i].back());
                            } else {
                                _exportcopies[i].back()->set_h(-1.);
                            }
                        }
                        delete import_ngbs[k];
                    }
                }
                for(unsigned int j = inlist[i].size(); j--;) {
                    NgbSearch ngbsearch(inlist[i][j], inradii[i][j], dummyflags,
                                        100);
                    parttree.walk_tree(ngbsearch);
                    vector<GasParticle*> import_ngbs = ngbsearch.get_ngbs();
                    for(unsigned int k = import_ngbs.size(); k--;) {
                        if(!import_ngbs[k]->is_exported(i) &&
                           !import_ngbs[k]->is_pseudo()) {
                            if(import_ngbs[k]->get_vorgen() == NULL) {
                                add_particle(import_ngbs[k],
                                             import_ngbs[k]->get_local_id());
                                counter++;
                            }
                            outlist[i].push_back(import_ngbs[k]);
                            import_ngbs[k]->do_export(i);
                            _exports[i].push_back(import_ngbs[k]);
                        }
                    }
                }
            }
        }

        // communicate back
        if(MPIGlobal::size > 1) {
            std::vector<MPI_Request> reqs((MPIGlobal::size - 1),
                                          MPI_REQUEST_NULL);
            vector<unsigned int> buffers((MPIGlobal::size - 1));
            unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
            for(unsigned int i = 0; i < buffers.size(); i++) {
                buffers[i] = i * bufsize;
            }

            unsigned int maxsize = bufsize / sizeof(GasParticle);
            if(!(bufsize % sizeof(GasParticle))) {
                maxsize--;
            }

            vector<int> allflag(MPIGlobal::size - 1);
            vector<unsigned int> numsend(MPIGlobal::size - 1);
            vector<unsigned int> sendsizes(MPIGlobal::size - 1);
            int allflagsum = 0;
            for(int i = 0; i < MPIGlobal::size - 1; i++) {
                sendsizes[i] =
                        outlist[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                .size();

                int send_pos = 0;
                for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]);
                    si++) {
                    outlist[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                           [si]->pack_data(&MPIGlobal::sendbuffer[buffers[i]],
                                           bufsize, &send_pos);
                }
                allflag[i] = (sendsizes[i] <= maxsize);
                allflagsum += allflag[i];
                numsend[i] = 1;
                MyMPI_Pack(&allflag[i], 1, MPI_INT,
                           &MPIGlobal::sendbuffer[buffers[i]], bufsize,
                           &send_pos);

                MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos,
                            MPI_PACKED,
                            (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0,
                            &reqs[i]);
            }

            for(int i = 0; i < (MPIGlobal::size - 1); i++) {
                MPI_Status status;
                for(int j = 0; j < MPIGlobal::size - 1; j++) {
                    if(!allflag[j]) {
                        int flag;
                        MyMPI_Test(&reqs[j], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[j] * maxsize;
                                si < std::min((numsend[j] + 1) * maxsize,
                                              sendsizes[j]);
                                si++) {
                                outlist[(MPIGlobal::rank + 1 + j) %
                                        MPIGlobal::size]
                                       [si]->pack_data(&MPIGlobal::sendbuffer
                                                               [buffers[j]],
                                                       bufsize, &send_pos);
                            }
                            allflag[j] = (sendsizes[j] <=
                                          (numsend[j] + 1) * maxsize);
                            allflagsum += allflag[j];
                            numsend[j]++;
                            // add continuation signal
                            MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[j]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    0, &reqs[j]);
                        }
                    }
                }
                int index;
                int tag;
                int recv_pos;

                // wait for a message from any source and with any tag to arrive
                MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
                // query the MPI_Status object to retrieve tag, source and size
                tag = status.MPI_TAG;
                index = status.MPI_SOURCE;
                int nelements;
                MyMPI_Get_count(&status, MPI_PACKED, &nelements);
                // select an index in the buffers to use for receiving and
                // sending
                int freebuffer;
                if(index == MPIGlobal::size - 1) {
                    freebuffer = MPIGlobal::rank;
                } else {
                    freebuffer = index;
                }
                if(tag == 0) {
                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                               nelements, MPI_PACKED, index, tag, &status);
                    recv_pos = 0;
                    while(recv_pos < nelements - 4) {
                        GasParticle part(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos);
                        if(part.h() < 0.) {
                            Vec pos = part.get_position();
#if ndim_ == 3
                            _mirrors.push_back(
                                    new VorGen(pos[0], pos[1], pos[2]));
#else
                            _mirrors.push_back(new VorGen(pos[0], pos[1]));
#endif
                        } else {
                            _ghosts[index].push_back(new GasParticle(part));
                            Vec pos = _ghosts[index].back()->get_position();
#if ndim_ == 3
                            _mirrors.push_back(
                                    new VorGen(pos[0], pos[1], pos[2]));
#else
                            _mirrors.push_back(new VorGen(pos[0], pos[1]));
#endif
                            _mirrors.back()->set_particle(
                                    _ghosts[index].back());
                            _ghosts[index].back()->set_vorgen(_mirrors.back());
                        }
                        _mirrors.back()->set_id(1);
                        _mirrors.back()->set_process(index);
                        _points.push_back(_mirrors.back());
                        add_point(_points.size() - 1);
                        counter++;
                    }
                    int flag;
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &flag, 1, MPI_INT);
                    if(!flag) {
                        i--;
                    }
                }
            }

            while(allflagsum < MPIGlobal::size - 1) {
                MPI_Status status;
                for(int j = 0; j < MPIGlobal::size - 1; j++) {
                    if(!allflag[j]) {
                        int flag;
                        MyMPI_Test(&reqs[j], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[j] * maxsize;
                                si < std::min((numsend[j] + 1) * maxsize,
                                              sendsizes[j]);
                                si++) {
                                outlist[(MPIGlobal::rank + 1 + j) %
                                        MPIGlobal::size]
                                       [si]->pack_data(&MPIGlobal::sendbuffer
                                                               [buffers[j]],
                                                       bufsize, &send_pos);
                            }
                            allflag[j] = (sendsizes[j] <=
                                          (numsend[j] + 1) * maxsize);
                            allflagsum += allflag[j];
                            numsend[j]++;
                            // add continuation signal
                            MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[j]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    0, &reqs[j]);
                        }
                    }
                }
            }
            // if we do not put a barrier here, the number of loops above might
            // get screwed up by interfering messages from below...
            vector<MPI_Status> status((MPIGlobal::size - 1));
            MyMPI_Waitall((MPIGlobal::size - 1), &reqs[0], &status[0]);
            MyMPI_Barrier();
        }

        // check if everything is ok now
        redo = 0;
        for(unsigned int i = 1; i < _simplices.size(); i++) {
            if(_simplices[i] != NULL) {
                VorGen* sp = _simplices[i]->get_special_point(_points);
                unsigned int* vorgens = _simplices[i]->get_vorgens();
#if ndim_ == 3
                VorGen* points[4] = {_points[vorgens[0]], _points[vorgens[1]],
                                     _points[vorgens[2]], _points[vorgens[3]]};
#else
                VorGen* points[3] = {_points[vorgens[0]], _points[vorgens[1]],
                                     _points[vorgens[2]]};
#endif
                double r = (sp->get_position() - points[0]->get_position())
                                   .norm2();
                for(unsigned int j = ndim_ + 1; j--;) {
                    if(points[j]->get_particle() != NULL &&
                       points[j]->get_id() != 1) {
                        double rad =
                                points[j]->get_particle()->get_max_radius();
                        if(4. * r >= rad * rad) {
                            points[j]->set_id(2);
                            redo++;
                        }
                    }
                }
            }
        }
        counter = 0;
        for(unsigned int i = startpoint; i < pointssize; i++) {
            if(_points[i]) {
                if(_points[i]->get_id() == 2) {
                    _points[i]->get_particle()->set_max_radius(
                            1.5 * _points[i]->get_particle()->get_max_radius());
                }
            }
        }
        if(worldsize > 1) {
            unsigned int redo_glob;
            MyMPI_Allreduce(&redo, &redo_glob, 1, MPI_UNSIGNED, MPI_SUM);
            redo = redo_glob;
        }
    }

    redo = 1;
// now add particles to complete the cells of the particles we added above
#if ndim_ == 2
    unsigned int pointssize2 = _points.size() - _mirrors.size() - 3;
    startpoint = 0;
#else
    unsigned int pointssize2 = _points.size() - _mirrors.size() - 4;
    startpoint = 0;
#endif
    while(redo) {
        vector<vector<Vec> > exportlist(worldsize);
        vector<vector<Vec> > exportlist_outside(worldsize);
        vector<vector<double> > exportradii(worldsize);
        vector<vector<double> > exportradii_outside(worldsize);
        // points may be added during the loop, we only want to loop over the
        // old points
        for(unsigned int i = startpoint; i < pointssize2; i++) {
            if(!_points[i]) {
                continue;
            }
            if(_points[i]->get_id()) {
                vector<bool> exportflag(worldsize, false);
                double rad = _points[i]->get_particle()->get_max_radius();
                Vec coords = _points[i]->get_position();
                if(_periodic) {
                    vector<VorGen*> ngbs;
                    ngbs.reserve(100);
                    parttree.get_neighbours_outside(coords, rad, ngbs,
                                                    exportflag);
                    for(unsigned int j = ngbs.size(); j--;) {
                        _mirrors.push_back(ngbs[j]);
                        _mirrors.back()->set_id(1);
                        _mirrors.back()->set_original(
                                ngbs[j]->get_particle()->get_local_id());
                        _points.push_back(ngbs[j]);
                        add_point(_points.size() - 1);
                        counter++;
                    }
                } else {
                    vector<VorGen*> ngbs;
                    ngbs.reserve(100);
                    parttree.get_neighbours_outside(coords, rad, ngbs,
                                                    exportflag);
                    for(unsigned int j = ngbs.size(); j--;) {
                        ngbs[j]->set_original(
                                ngbs[j]->get_particle()->get_local_id());
                        ngbs[j]->set_particle(NULL);
                        _mirrors.push_back(ngbs[j]);
                        _mirrors.back()->set_id(1);
                        _points.push_back(ngbs[j]);
                        add_point(_points.size() - 1);
                        counter++;
                    }
                }
                for(int j = worldsize; j--;) {
                    if(exportflag[j]) {
                        exportlist_outside[j].push_back(coords);
                        exportradii_outside[j].push_back(rad);
                        exportflag[j] = false;
                    }
                }
                NgbSearch ngbsearch(coords, rad, exportflag, 100);
                parttree.walk_tree(ngbsearch);
                vector<GasParticle*> inside_ngbs = ngbsearch.get_ngbs();
                for(unsigned int j = 0; j < inside_ngbs.size(); j++) {
                    if(inside_ngbs[j] != NULL &&
                       inside_ngbs[j]->get_vorgen() == NULL) {
#if ndim_ == 3
                        _mirrors.push_back(new VorGen(inside_ngbs[j]->x(),
                                                      inside_ngbs[j]->y(),
                                                      inside_ngbs[j]->z()));
#else
                        _mirrors.push_back(new VorGen(inside_ngbs[j]->x(),
                                                      inside_ngbs[j]->y()));
#endif
                        inside_ngbs[j]->set_vorgen(_mirrors.back());
                        _mirrors.back()->set_particle(inside_ngbs[j]);
                        _points.push_back(_mirrors.back());
                        add_point(_points.size() - 1);
                        _mirrors.back()->set_id(1);
                        counter++;
                    }
                }
                for(int j = worldsize; j--;) {
                    if(exportflag[j]) {
                        exportlist[j].push_back(coords);
                        exportradii[j].push_back(rad);
                    }
                }
                _points[i]->set_id(0);
            }
        }

        // communicate
        vector<vector<Vec> > inlist(worldsize);
        vector<vector<double> > inradii(worldsize);
        vector<vector<Vec> > inlist_outside(worldsize);
        vector<vector<double> > inradii_outside(worldsize);
        if(MPIGlobal::size > 1) {
            std::vector<MPI_Request> reqs((MPIGlobal::size - 1) * 2,
                                          MPI_REQUEST_NULL);
            vector<unsigned int> buffers((MPIGlobal::size - 1) * 2);
            unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
            for(unsigned int i = 0; i < buffers.size(); i++) {
                buffers[i] = i * bufsize;
            }

            unsigned int maxsize = bufsize / ((ndim_ + 1) * sizeof(double));
            // make sure we have some space left at the end of the buffer to
            // signal more to come
            if(!(bufsize % ((ndim_ + 1) * sizeof(double)))) {
                maxsize--;
            }

            vector<int> allflag((MPIGlobal::size - 1) * 2);
            vector<unsigned int> numsend((MPIGlobal::size - 1) * 2);
            vector<unsigned int> sendsizes(2 * (MPIGlobal::size - 1));

            // for some reason, reusing the same buffer for every send does not
            // work
            // we therefore use a separate send buffer for every process
            int allflagsum = 0;
            for(int i = 0; i < MPIGlobal::size - 1; i++) {
                // initialize send sizes
                sendsizes[2 * i] =
                        exportlist[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                .size();
                sendsizes[2 * i + 1] =
                        exportlist_outside[(MPIGlobal::rank + 1 + i) %
                                           MPIGlobal::size]
                                .size();

                int send_pos = 0;
                for(unsigned int si = 0;
                    si < std::min(maxsize, sendsizes[2 * i]); si++) {
                    MyMPI_Pack(&exportlist[(MPIGlobal::rank + 1 + i) %
                                           MPIGlobal::size][si][0],
                               ndim_, MPI_DOUBLE,
                               &MPIGlobal::sendbuffer[buffers[2 * i]], bufsize,
                               &send_pos);
                    MyMPI_Pack(&exportradii[(MPIGlobal::rank + 1 + i) %
                                            MPIGlobal::size][si],
                               1, MPI_DOUBLE,
                               &MPIGlobal::sendbuffer[buffers[2 * i]], bufsize,
                               &send_pos);
                }
                allflag[2 * i] = (sendsizes[2 * i] <= maxsize);
                // if allflagsum equals MPIGlobal::size-1, all data has been
                // sent
                allflagsum += allflag[2 * i];
                numsend[2 * i] = 1;
                // add continuation signal
                MyMPI_Pack(&allflag[2 * i], 1, MPI_INT,
                           &MPIGlobal::sendbuffer[buffers[2 * i]], bufsize,
                           &send_pos);

                MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[2 * i]], send_pos,
                            MPI_PACKED,
                            (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0,
                            &reqs[2 * i]);

                send_pos = 0;
                for(unsigned int si = 0;
                    si < std::min(maxsize, sendsizes[2 * i + 1]); si++) {
                    MyMPI_Pack(&exportlist_outside[(MPIGlobal::rank + 1 + i) %
                                                   MPIGlobal::size][si][0],
                               ndim_, MPI_DOUBLE,
                               &MPIGlobal::sendbuffer[buffers[2 * i + 1]],
                               bufsize, &send_pos);
                    MyMPI_Pack(&exportradii_outside[(MPIGlobal::rank + 1 + i) %
                                                    MPIGlobal::size][si],
                               1, MPI_DOUBLE,
                               &MPIGlobal::sendbuffer[buffers[2 * i + 1]],
                               bufsize, &send_pos);
                }
                allflag[2 * i + 1] = (sendsizes[2 * i + 1] <= maxsize);
                // if allflagsum equals MPIGlobal::size-1, all data has been
                // sent
                allflagsum += allflag[2 * i + 1];
                numsend[2 * i + 1] = 1;
                // add continuation signal
                MyMPI_Pack(&allflag[2 * i + 1], 1, MPI_INT,
                           &MPIGlobal::sendbuffer[buffers[2 * i + 1]], bufsize,
                           &send_pos);

                MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[2 * i + 1]],
                            send_pos, MPI_PACKED,
                            (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 1,
                            &reqs[2 * i + 1]);
            }

            for(int i = 0; i < 2 * (MPIGlobal::size - 1); i++) {
                MPI_Status status;
                for(int j = 0; j < MPIGlobal::size - 1; j++) {
                    if(!allflag[2 * j]) {
                        int flag;
                        MyMPI_Test(&reqs[2 * j], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[2 * j] * maxsize;
                                si < std::min((numsend[2 * j] + 1) * maxsize,
                                              sendsizes[2 * j]);
                                si++) {
                                MyMPI_Pack(
                                        &exportlist[(MPIGlobal::rank + 1 + j) %
                                                    MPIGlobal::size][si][0],
                                        ndim_, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j]],
                                        bufsize, &send_pos);
                                MyMPI_Pack(
                                        &exportradii[(MPIGlobal::rank + 1 + j) %
                                                     MPIGlobal::size][si],
                                        1, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j]],
                                        bufsize, &send_pos);
                            }
                            allflag[2 * j] = (sendsizes[2 * j] <=
                                              (numsend[2 * j] + 1) * maxsize);
                            allflagsum += allflag[2 * j];
                            numsend[2 * j]++;
                            // add continuation signal
                            MyMPI_Pack(&allflag[2 * j], 1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[2 * j]],
                                       bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[2 * j]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    0, &reqs[2 * j]);
                        }
                    }
                    if(!allflag[2 * j + 1]) {
                        int flag;
                        MyMPI_Test(&reqs[2 * j + 1], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[2 * j + 1] * maxsize;
                                si <
                                std::min((numsend[2 * j + 1] + 1) * maxsize,
                                         sendsizes[2 * j + 1]);
                                si++) {
                                MyMPI_Pack(
                                        &exportlist_outside[(MPIGlobal::rank +
                                                             1 + j) %
                                                            MPIGlobal::size][si]
                                                           [0],
                                        ndim_, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j +
                                                                       1]],
                                        bufsize, &send_pos);
                                MyMPI_Pack(
                                        &exportradii_outside[(MPIGlobal::rank +
                                                              1 + j) %
                                                             MPIGlobal::size]
                                                            [si],
                                        1, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j +
                                                                       1]],
                                        bufsize, &send_pos);
                            }
                            allflag[2 * j + 1] =
                                    (sendsizes[2 * j + 1] <=
                                     (numsend[2 * j + 1] + 1) * maxsize);
                            allflagsum += allflag[2 * j + 1];
                            numsend[2 * j + 1]++;
                            // add continuation signal
                            MyMPI_Pack(
                                    &allflag[2 * j + 1], 1, MPI_INT,
                                    &MPIGlobal::sendbuffer[buffers[2 * j + 1]],
                                    bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[2 * j + 1]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    1, &reqs[2 * j + 1]);
                        }
                    }
                }

                int index;
                int tag;
                int recv_pos;

                // wait for a message from any source and with any tag to arrive
                MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
                // query the MPI_Status object to retrieve tag, source and size
                tag = status.MPI_TAG;
                index = status.MPI_SOURCE;
                int nelements;
                MyMPI_Get_count(&status, MPI_PACKED, &nelements);
                // select an index in the buffers to use for receiving and
                // sending
                int freebuffer;
                if(index == MPIGlobal::size - 1) {
                    freebuffer = 2 * MPIGlobal::rank;
                } else {
                    freebuffer = 2 * index;
                }
                if(tag == 0) {
                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                               nelements, MPI_PACKED, index, tag, &status);
                    recv_pos = 0;
                    while(recv_pos < nelements - 4) {
                        Vec inlist_el;
                        double inradii_el;
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &inlist_el[0], ndim_,
                                MPI_DOUBLE);
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &inradii_el, 1,
                                MPI_DOUBLE);
                        inlist[index].push_back(inlist_el);
                        inradii[index].push_back(inradii_el);
                    }
                    int flag;
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &flag, 1, MPI_INT);
                    if(!flag) {
                        i--;
                    }
                }
                if(tag == 1) {
                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                               nelements, MPI_PACKED, index, tag, &status);
                    recv_pos = 0;
                    while(recv_pos < nelements - 4) {
                        Vec inlist_outside_el;
                        double inradii_outside_el;
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &inlist_outside_el[0],
                                ndim_, MPI_DOUBLE);
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &inradii_outside_el, 1,
                                MPI_DOUBLE);
                        inlist_outside[index].push_back(inlist_outside_el);
                        inradii_outside[index].push_back(inradii_outside_el);
                    }
                    if(inlist_outside[index].size() > 0 &&
                       !other_indices[index]) {
                        other_indices[index] = lastindex;
                        lastindex++;
                    }
                    int flag;
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &flag, 1, MPI_INT);
                    if(!flag) {
                        i--;
                    }
                }
            }

            while(allflagsum < 2 * (MPIGlobal::size - 1)) {
                MPI_Status status;
                for(int j = 0; j < MPIGlobal::size - 1; j++) {
                    if(!allflag[2 * j]) {
                        int flag;
                        MyMPI_Test(&reqs[2 * j], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[2 * j] * maxsize;
                                si < std::min((numsend[2 * j] + 1) * maxsize,
                                              sendsizes[2 * j]);
                                si++) {
                                MyMPI_Pack(
                                        &exportlist[(MPIGlobal::rank + 1 + j) %
                                                    MPIGlobal::size][si],
                                        ndim_, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j]],
                                        bufsize, &send_pos);
                                MyMPI_Pack(
                                        &exportradii[(MPIGlobal::rank + 1 + j) %
                                                     MPIGlobal::size][si],
                                        1, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j]],
                                        bufsize, &send_pos);
                            }
                            allflag[2 * j] = (sendsizes[2 * j] <=
                                              (numsend[2 * j] + 1) * maxsize);
                            allflagsum += allflag[2 * j];
                            numsend[2 * j]++;
                            // add continuation signal
                            MyMPI_Pack(&allflag[2 * j], 1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[2 * j]],
                                       bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[2 * j]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    0, &reqs[2 * j]);
                        }
                    }
                    if(!allflag[2 * j + 1]) {
                        int flag;
                        MyMPI_Test(&reqs[2 * j + 1], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[2 * j + 1] * maxsize;
                                si <
                                std::min((numsend[2 * j + 1] + 1) * maxsize,
                                         sendsizes[2 * j + 1]);
                                si++) {
                                MyMPI_Pack(
                                        &exportlist_outside[(MPIGlobal::rank +
                                                             1 + j) %
                                                            MPIGlobal::size]
                                                           [si],
                                        ndim_, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j +
                                                                       1]],
                                        bufsize, &send_pos);
                                MyMPI_Pack(
                                        &exportradii_outside[(MPIGlobal::rank +
                                                              1 + j) %
                                                             MPIGlobal::size]
                                                            [si],
                                        1, MPI_DOUBLE,
                                        &MPIGlobal::sendbuffer[buffers[2 * j +
                                                                       1]],
                                        bufsize, &send_pos);
                            }
                            allflag[2 * j + 1] =
                                    (sendsizes[2 * j + 1] <=
                                     (numsend[2 * j + 1] + 1) * maxsize);
                            allflagsum += allflag[2 * j + 1];
                            numsend[2 * j + 1]++;
                            // add continuation signal
                            MyMPI_Pack(
                                    &allflag[2 * j + 1], 1, MPI_INT,
                                    &MPIGlobal::sendbuffer[buffers[2 * j + 1]],
                                    bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[2 * j + 1]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    1, &reqs[2 * j + 1]);
                        }
                    }
                }
            }
            vector<MPI_Status> status((MPIGlobal::size - 1) * 2);
            MyMPI_Waitall((MPIGlobal::size - 1) * 2, &reqs[0], &status[0]);
            MyMPI_Barrier();
        }

        // search for imported coordinates
        vector<vector<Vec> > outlist(worldsize);
        for(int i = worldsize; i--;) {
            if(worldrank != i) {
                for(unsigned int j = inlist_outside[i].size(); j--;) {
                    vector<VorGen*> import_ngbs;
                    import_ngbs.reserve(100);
                    parttree.get_neighbours_outside(
                            inlist_outside[i][j], inradii_outside[i][j],
                            import_ngbs, dummyflags, other_indices[i]);
                    for(unsigned int k = import_ngbs.size(); k--;) {
                        if(!import_ngbs[k]->get_particle()->is_pseudo()) {
                            outlist[i].push_back(
                                    import_ngbs[k]->get_position());
                        }
                        delete import_ngbs[k];
                    }
                }
                for(unsigned int j = inlist[i].size(); j--;) {
                    NgbSearch ngbsearch(inlist[i][j], inradii[i][j], dummyflags,
                                        100);
                    parttree.walk_tree(ngbsearch);
                    vector<GasParticle*> import_ngbs = ngbsearch.get_ngbs();
                    for(unsigned int k = import_ngbs.size(); k--;) {
                        if(!import_ngbs[k]->is_exported(i) &&
                           !import_ngbs[k]->is_pseudo()) {
                            outlist[i].push_back(
                                    import_ngbs[k]->get_position());
                            import_ngbs[k]->do_export(i);
                        }
                    }
                }
            }
        }

        // communicate back
        if(MPIGlobal::size > 1) {
            std::vector<MPI_Request> reqs((MPIGlobal::size - 1),
                                          MPI_REQUEST_NULL);
            vector<unsigned int> buffers((MPIGlobal::size - 1));
            unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
            for(unsigned int i = 0; i < buffers.size(); i++) {
                buffers[i] = i * bufsize;
            }

            unsigned int maxsize = bufsize / sizeof(Vec);
            if(!(bufsize % sizeof(Vec))) {
                maxsize--;
            }

            vector<int> allflag(MPIGlobal::size - 1);
            vector<unsigned int> numsend(MPIGlobal::size - 1);
            vector<unsigned int> sendsizes(MPIGlobal::size - 1);
            int allflagsum = 0;
            for(int i = 0; i < MPIGlobal::size - 1; i++) {
                sendsizes[i] =
                        outlist[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                .size();

                int send_pos = 0;
                for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]);
                    si++) {
                    MyMPI_Pack(&outlist[(MPIGlobal::rank + 1 + i) %
                                        MPIGlobal::size][si][0],
                               ndim_, MPI_DOUBLE,
                               &MPIGlobal::sendbuffer[buffers[i]], bufsize,
                               &send_pos);
                }
                allflag[i] = (sendsizes[i] <= maxsize);
                allflagsum += allflag[i];
                numsend[i] = 1;
                MyMPI_Pack(&allflag[i], 1, MPI_INT,
                           &MPIGlobal::sendbuffer[buffers[i]], bufsize,
                           &send_pos);

                MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos,
                            MPI_PACKED,
                            (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0,
                            &reqs[i]);
            }

            for(int i = 0; i < (MPIGlobal::size - 1); i++) {
                MPI_Status status;
                for(int j = 0; j < MPIGlobal::size - 1; j++) {
                    if(!allflag[j]) {
                        int flag;
                        MyMPI_Test(&reqs[j], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[j] * maxsize;
                                si < std::min((numsend[j] + 1) * maxsize,
                                              sendsizes[j]);
                                si++) {
                                MyMPI_Pack(&outlist[(MPIGlobal::rank + 1 + j) %
                                                    MPIGlobal::size][si][0],
                                           ndim_, MPI_DOUBLE,
                                           &MPIGlobal::sendbuffer[buffers[j]],
                                           bufsize, &send_pos);
                            }
                            allflag[j] = (sendsizes[j] <=
                                          (numsend[j] + 1) * maxsize);
                            allflagsum += allflag[j];
                            numsend[j]++;
                            MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[j]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    0, &reqs[j]);
                        }
                    }
                }

                int index;
                int tag;
                int recv_pos;

                // wait for a message from any source and with any tag to arrive
                MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
                // query the MPI_Status object to retrieve tag, source and size
                tag = status.MPI_TAG;
                index = status.MPI_SOURCE;
                int nelements;
                MyMPI_Get_count(&status, MPI_PACKED, &nelements);
                // select an index in the buffers to use for receiving and
                // sending
                int freebuffer;
                if(index == MPIGlobal::size - 1) {
                    freebuffer = MPIGlobal::rank;
                } else {
                    freebuffer = index;
                }
                if(tag == 0) {
                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                               nelements, MPI_PACKED, index, tag, &status);
                    recv_pos = 0;
                    while(recv_pos < nelements - 4) {
                        Vec pos;
                        MyMPI_Unpack(
                                &MPIGlobal::recvbuffer[buffers[freebuffer]],
                                nelements, &recv_pos, &pos[0], ndim_,
                                MPI_DOUBLE);
#if ndim_ == 3
                        _mirrors.push_back(new VorGen(pos[0], pos[1], pos[2]));
#else
                        _mirrors.push_back(new VorGen(pos[0], pos[1]));
#endif
                        _mirrors.back()->set_process(index);
                        _mirrors.back()->set_id(1);
                        _points.push_back(_mirrors.back());
                        add_point(_points.size() - 1);
                        counter++;
                    }
                    int flag;
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &flag, 1, MPI_INT);
                    if(!flag) {
                        i--;
                    }
                }
            }

            while(allflagsum < MPIGlobal::size - 1) {
                MPI_Status status;
                for(int j = 0; j < MPIGlobal::size - 1; j++) {
                    if(!allflag[j]) {
                        int flag;
                        MyMPI_Test(&reqs[j], &flag, &status);
                        if(flag) {
                            int send_pos = 0;
                            for(unsigned int si = numsend[j] * maxsize;
                                si < std::min((numsend[j] + 1) * maxsize,
                                              sendsizes[j]);
                                si++) {
                                MyMPI_Pack(&outlist[(MPIGlobal::rank + 1 + j) %
                                                    MPIGlobal::size][si][0],
                                           ndim_, MPI_DOUBLE,
                                           &MPIGlobal::sendbuffer[buffers[j]],
                                           bufsize, &send_pos);
                            }
                            allflag[j] = (sendsizes[j] <=
                                          (numsend[j] + 1) * maxsize);
                            allflagsum += allflag[j];
                            numsend[j]++;
                            MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);

                            MyMPI_Isend(
                                    &MPIGlobal::sendbuffer[buffers[j]],
                                    send_pos, MPI_PACKED,
                                    (MPIGlobal::rank + 1 + j) % MPIGlobal::size,
                                    0, &reqs[j]);
                        }
                    }
                }
            }

            vector<MPI_Status> status((MPIGlobal::size - 1));
            MyMPI_Waitall((MPIGlobal::size - 1), &reqs[0], &status[0]);
            MyMPI_Barrier();
        }

        // check if everything is ok now
        redo = 0;
        for(unsigned int i = 1; i < _simplices.size(); i++) {
            if(_simplices[i] != NULL) {
                VorGen* sp = _simplices[i]->get_special_point(_points);
                unsigned int* vorgens = _simplices[i]->get_vorgens();
#if ndim_ == 3
                VorGen* points[4] = {_points[vorgens[0]], _points[vorgens[1]],
                                     _points[vorgens[2]], _points[vorgens[3]]};
#else
                VorGen* points[3] = {_points[vorgens[0]], _points[vorgens[1]],
                                     _points[vorgens[2]]};
#endif
                double r = (sp->get_position() - points[0]->get_position())
                                   .norm2();
                for(unsigned int j = ndim_ + 1; j--;) {
                    if(points[j]->get_particle() != NULL &&
                       points[j]->get_id() != 1 &&
                       4. * r >= points[j]->get_particle()->get_max_radius() *
                                         points[j]
                                                 ->get_particle()
                                                 ->get_max_radius()) {
                        points[j]->set_id(1);
                        redo++;
                    }
                }
            }
        }
        counter = 0;
        for(unsigned int i = startpoint; i < pointssize2; i++) {
            if(!_points[i]) {
                continue;
            }
            if(_points[i]->get_id()) {
                _points[i]->get_particle()->set_max_radius(
                        1.5 * _points[i]->get_particle()->get_max_radius());
            }
        }
        if(worldsize > 1) {
            unsigned int redo_glob;
            MyMPI_Allreduce(&redo, &redo_glob, 1, MPI_UNSIGNED, MPI_SUM);
            redo = redo_glob;
        }
    }

#if ndim_ == 2
    pointssize = _points.size() - _mirrors.size() - 3;
    startpoint = 0;
#else
    pointssize = _points.size() - _mirrors.size() - 4;
    startpoint = 0;
#endif
    for(unsigned int i = startpoint; i < pointssize; i++) {
        if(_points[i]) {
            _points[i]->get_particle()->set_max_radius(0.);
            _points[i]->set_id(2);
        }
    }
    for(unsigned int i = 1; i < _simplices.size(); i++) {
        if(_simplices[i] != NULL) {
            VorGen* sp = _simplices[i]->get_special_point(_points);
            unsigned int* vorgens = _simplices[i]->get_vorgens();
#if ndim_ == 3
            VorGen* points[4] = {_points[vorgens[0]], _points[vorgens[1]],
                                 _points[vorgens[2]], _points[vorgens[3]]};
#else
            VorGen* points[3] = {_points[vorgens[0]], _points[vorgens[1]],
                                 _points[vorgens[2]]};
#endif
            double r = (sp->get_position() - points[0]->get_position()).norm();
            for(unsigned int j = ndim_ + 1; j--;) {
                if(points[j]->get_id() == 2) {
                    points[j]->get_particle()->set_max_radius(
                            max(points[j]->get_particle()->get_max_radius(),
                                2. * r));
                }
            }
        }
    }

    LOGS("Done adding mirror points");
}

/**
  * @brief Communicate primitive variables between processes.
  *
  * To calculate the primitive variables, a VorCell is needed. This can only be
  * constructed on the "home" process of a Particle (the process where the
  * original Particle resides). We then communicate the primitive variables to
  * all processes which use a ghost copy of the original Particle, to ensure
  * that these processes use the correct primitive variables.
  */
void DelTess::update_Ws() {
    if(MPIGlobal::size < 2) {
        return;
    }

    vector<MPI_Request> reqs(MPIGlobal::size - 1, MPI_REQUEST_NULL);
    vector<unsigned int> buffers(MPIGlobal::size - 1);
    unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
    for(unsigned int i = 0; i < buffers.size(); i++) {
        buffers[i] = i * bufsize;
    }

    unsigned int maxsize = bufsize / sizeof(StateVector);
    if(!(bufsize % sizeof(StateVector))) {
        maxsize--;
    }

    vector<int> allflag(MPIGlobal::size - 1);
    vector<unsigned int> numsend(MPIGlobal::size - 1);
    vector<unsigned int> sendsizes(MPIGlobal::size - 1);
    int allflagsum = 0;
    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        sendsizes[i] =
                _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size].size();

        int send_pos = 0;
        for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]); si++) {
            StateVector W =
                    _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                            [si]->get_Wvec();
            MyMPI_Pack(&W[0], ndim_ + 2, MPI_DOUBLE,
                       &MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);
        }
        allflag[i] = (sendsizes[i] <= maxsize);
        // if allflagsum equals MPIGlobal::size-1, all data has been sent
        allflagsum += allflag[i];
        numsend[i] = 1;
        // add continuation signal
        MyMPI_Pack(&allflag[i], 1, MPI_INT, &MPIGlobal::sendbuffer[buffers[i]],
                   bufsize, &send_pos);

        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos, MPI_PACKED,
                    (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0, &reqs[i]);
    }

    vector<unsigned int> numreceived(MPIGlobal::size - 1, 0);
    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        StateVector W = _exports[(MPIGlobal::rank + 1 + j) %
                                                 MPIGlobal::size]
                                                [si]->get_Wvec();
                        MyMPI_Pack(&W[0], ndim_ + 2, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }

        int index;
        int tag;
        int recv_pos;

        // wait for a message from any source and with any tag to arrive
        MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
        // query the MPI_Status object to retrieve tag, source and size
        tag = status.MPI_TAG;
        index = status.MPI_SOURCE;
        int nelements;
        MyMPI_Get_count(&status, MPI_PACKED, &nelements);
        // select an index in the buffers to use for receiving and sending
        int freebuffer;
        if(index == MPIGlobal::size - 1) {
            freebuffer = MPIGlobal::rank;
        } else {
            freebuffer = index;
        }

        if(tag == 0) {
            MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                       MPI_PACKED, index, tag, &status);
            recv_pos = 0;
            unsigned int j = numreceived[freebuffer];
            while(recv_pos < nelements - 4) {
                StateVector W;
                MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                             nelements, &recv_pos, &W[0], ndim_ + 2,
                             MPI_DOUBLE);
                _ghosts[index][j]->set_W(W);
                j++;
            }
            numreceived[freebuffer] = j;
            int flag;
            MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                         &recv_pos, &flag, 1, MPI_INT);
            if(!flag) {
                i--;
            }
        }
    }

    while(allflagsum < MPIGlobal::size - 1) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        StateVector W = _exports[(MPIGlobal::rank + 1 + j) %
                                                 MPIGlobal::size]
                                                [si]->get_Wvec();
                        MyMPI_Pack(&W[0], ndim_ + 2, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }
    }

    vector<MPI_Status> status((MPIGlobal::size - 1));
    MyMPI_Waitall((MPIGlobal::size - 1), &reqs[0], &status[0]);
    MyMPI_Barrier();
}

/**
  * @brief Communicate gradients, velocities, total surface areas and centroids
  * between processes.
  *
  * These quantities are VorCell bound and hence can only be calculated on the
  * "home" process (the process where the original Particle resides), since only
  * there a VorCell is calculated for the Particle. The quantities are then
  * communicated to all processes using a ghost copy of the Particle and used to
  * set the quantities for the ghost particles.
  */
void DelTess::update_gradients() {
    if(MPIGlobal::size < 2) {
        return;
    }

    vector<MPI_Request> reqs((MPIGlobal::size - 1), MPI_REQUEST_NULL);
    vector<unsigned int> buffers((MPIGlobal::size - 1));
    unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
    for(unsigned int i = 0; i < buffers.size(); i++) {
        buffers[i] = i * bufsize;
    }

    unsigned int maxsize = bufsize / (ndim_ * sizeof(StateVector) +
                                      sizeof(Vec) + sizeof(double));
    if(!(bufsize %
         (ndim_ * sizeof(StateVector) + sizeof(Vec) + sizeof(double)))) {
        maxsize--;
    }

    vector<int> allflag(MPIGlobal::size - 1);
    vector<unsigned int> numsend(MPIGlobal::size - 1);
    vector<unsigned int> sendsizes(MPIGlobal::size - 1);
    int allflagsum = 0;
    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        sendsizes[i] =
                _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size].size();

        int send_pos = 0;
        for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]); si++) {
            StateVector gradients[ndim_];
            _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                    [si]->get_vorgen()
                            ->get_particle()
                            ->get_gradients(gradients);
            Vec v = _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                            [si]->get_vorgen()
                                    ->get_particle()
                                    ->get_velocity();
            double A = _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                               [si]->get_vorgen()
                                       ->get_particle()
                                       ->get_total_area();
            MyMPI_Pack(&gradients[0][0], ndim_ * (ndim_ + 2), MPI_DOUBLE,
                       &MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);
            MyMPI_Pack(&v[0], ndim_, MPI_DOUBLE,
                       &MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);
            MyMPI_Pack(&A, 1, MPI_DOUBLE, &MPIGlobal::sendbuffer[buffers[i]],
                       bufsize, &send_pos);
        }
        allflag[i] = (sendsizes[i] <= maxsize);
        // if allflagsum equals MPIGlobal::size-1, all data has been sent
        allflagsum += allflag[i];
        numsend[i] = 1;
        // add continuation signal
        MyMPI_Pack(&allflag[i], 1, MPI_INT, &MPIGlobal::sendbuffer[buffers[i]],
                   bufsize, &send_pos);

        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos, MPI_PACKED,
                    (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0, &reqs[i]);
    }

    vector<unsigned int> numreceived(MPIGlobal::size - 1, 0);
    for(int i = 0; i < (MPIGlobal::size - 1); i++) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        StateVector gradients[ndim_];
                        _exports[(MPIGlobal::rank + 1 + j) % MPIGlobal::size]
                                [si]->get_vorgen()
                                        ->get_particle()
                                        ->get_gradients(gradients);
                        Vec v = _exports[(MPIGlobal::rank + 1 + j) %
                                         MPIGlobal::size]
                                        [si]->get_vorgen()
                                                ->get_particle()
                                                ->get_velocity();
                        double A = _exports[(MPIGlobal::rank + 1 + j) %
                                            MPIGlobal::size]
                                           [si]->get_vorgen()
                                                   ->get_particle()
                                                   ->get_total_area();
                        MyMPI_Pack(&gradients[0][0], ndim_ * (ndim_ + 2),
                                   MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                        MyMPI_Pack(&v[0], ndim_, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                        MyMPI_Pack(&A, 1, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }

        int index;
        int tag;
        int recv_pos;

        // wait for a message from any source and with any tag to arrive
        MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
        // query the MPI_Status object to retrieve tag, source and size
        tag = status.MPI_TAG;
        index = status.MPI_SOURCE;
        int nelements;
        MyMPI_Get_count(&status, MPI_PACKED, &nelements);
        // select an index in the buffers to use for receiving and sending
        int freebuffer;
        if(index == MPIGlobal::size - 1) {
            freebuffer = MPIGlobal::rank;
        } else {
            freebuffer = index;
        }

        if(tag == 0) {
            MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                       MPI_PACKED, index, tag, &status);
            recv_pos = 0;
            unsigned int j = numreceived[freebuffer];
            while(recv_pos < nelements - 4) {
                StateVector gradients[ndim_];
                MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                             nelements, &recv_pos, &gradients[0][0],
                             ndim_ * (ndim_ + 2), MPI_DOUBLE);
                _ghosts[index][j]->set_gradients(gradients);
                Vec v;
                MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                             nelements, &recv_pos, &v[0], ndim_, MPI_DOUBLE);
#if ndim_ == 3
                _ghosts[index][j]->set_v(v[0], v[1], v[2]);
#else
                _ghosts[index][j]->set_v(v[0], v[1], 0.);
#endif
                double A;
                MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                             nelements, &recv_pos, &A, 1, MPI_DOUBLE);
                _ghosts[index][j]->set_total_area(A);
                j++;
            }
            numreceived[freebuffer] = j;
            int flag;
            MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                         &recv_pos, &flag, 1, MPI_INT);

            if(!flag) {
                i--;
            }
        }
    }

    while(allflagsum < MPIGlobal::size - 1) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        StateVector gradients[ndim_];
                        _exports[(MPIGlobal::rank + 1 + j) % MPIGlobal::size]
                                [si]->get_vorgen()
                                        ->get_particle()
                                        ->get_gradients(gradients);
                        Vec v = _exports[(MPIGlobal::rank + 1 + j) %
                                         MPIGlobal::size]
                                        [si]->get_vorgen()
                                                ->get_particle()
                                                ->get_velocity();
                        double A = _exports[(MPIGlobal::rank + 1 + j) %
                                            MPIGlobal::size]
                                           [si]->get_vorgen()
                                                   ->get_particle()
                                                   ->get_total_area();
                        MyMPI_Pack(&gradients[0][0], ndim_ * (ndim_ + 2),
                                   MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                        MyMPI_Pack(&v[0], ndim_, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                        MyMPI_Pack(&A, 1, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }
    }

    vector<MPI_Status> status((MPIGlobal::size - 1));
    MyMPI_Waitall((MPIGlobal::size - 1), &reqs[0], &status[0]);
    MyMPI_Barrier();
}

/**
  * @brief Communicate dQs between processes
  *
  * Every ghost particle has registered incoming fluxes and added them to its
  * local dQ. These dQs are now sent back to the corresponding "home" processes
  * of the ghosts (the process where the original particle resides) and are
  * added to the dQ of the original particle.
  */
void DelTess::update_dQs() {
    if(MPIGlobal::size < 2) {
        return;
    }

    vector<MPI_Request> reqs((MPIGlobal::size - 1), MPI_REQUEST_NULL);
    vector<unsigned int> buffers((MPIGlobal::size - 1));
    unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
    for(unsigned int i = 0; i < buffers.size(); i++) {
        buffers[i] = i * bufsize;
    }

    unsigned int maxsize = bufsize / (sizeof(StateVector) + sizeof(Vec));
    if(!(bufsize % (sizeof(StateVector) + sizeof(Vec)))) {
        maxsize--;
    }

    vector<int> allflag(MPIGlobal::size - 1);
    vector<unsigned int> numsend(MPIGlobal::size - 1);
    vector<unsigned int> sendsizes(MPIGlobal::size - 1);
    int allflagsum = 0;
    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        sendsizes[i] =
                _ghosts[(MPIGlobal::rank + 1 + i) % MPIGlobal::size].size();

        int send_pos = 0;
        for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]); si++) {
            StateVector dQ =
                    _ghosts[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                           [si]->get_dQvec();
            Vec dE = _ghosts[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                            [si]->get_delta_E();
            MyMPI_Pack(&dQ[0], ndim_ + 2, MPI_DOUBLE,
                       &MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);
            MyMPI_Pack(&dE[0], ndim_, MPI_DOUBLE,
                       &MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);
        }
        allflag[i] = (sendsizes[i] <= maxsize);
        allflagsum += allflag[i];
        numsend[i] = 1;
        MyMPI_Pack(&allflag[i], 1, MPI_INT, &MPIGlobal::sendbuffer[buffers[i]],
                   bufsize, &send_pos);

        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos, MPI_PACKED,
                    (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0, &reqs[i]);
    }

    vector<unsigned int> numreceived(MPIGlobal::size - 1, 0);
    for(int i = 0; i < (MPIGlobal::size - 1); i++) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        StateVector dQ = _ghosts[(MPIGlobal::rank + 1 + j) %
                                                 MPIGlobal::size]
                                                [si]->get_dQvec();
                        Vec dE = _ghosts[(MPIGlobal::rank + 1 + j) %
                                         MPIGlobal::size]
                                        [si]->get_delta_E();
                        MyMPI_Pack(&dQ[0], ndim_ + 2, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                        MyMPI_Pack(&dE[0], ndim_, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }

        int index;
        int tag;
        int recv_pos;

        // wait for a message from any source and with any tag to arrive
        MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
        // query the MPI_Status object to retrieve tag, source and size
        tag = status.MPI_TAG;
        index = status.MPI_SOURCE;
        int nelements;
        MyMPI_Get_count(&status, MPI_PACKED, &nelements);
        // select an index in the buffers to use for receiving and sending
        int freebuffer;
        if(index == MPIGlobal::size - 1) {
            freebuffer = MPIGlobal::rank;
        } else {
            freebuffer = index;
        }

        if(tag == 0) {
            MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                       MPI_PACKED, index, tag, &status);
            recv_pos = 0;
            unsigned int j = numreceived[freebuffer];
            while(recv_pos < nelements - 4) {
                StateVector dQ;
                Vec dE;
                MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                             nelements, &recv_pos, &dQ[0], ndim_ + 2,
                             MPI_DOUBLE);
                MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                             nelements, &recv_pos, &dE[0], ndim_, MPI_DOUBLE);
                _exports[index][j]->get_vorgen()->get_particle()->increase_dQ(
                        dQ);
                _exports[index]
                        [j]->get_vorgen()
                                ->get_particle()
                                ->increase_delta_E(dE);
                j++;
            }
            numreceived[freebuffer] = j;
            int flag;
            MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                         &recv_pos, &flag, 1, MPI_INT);
            if(!flag) {
                i--;
            }
        }
    }

    while(allflagsum < MPIGlobal::size - 1) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        StateVector dQ = _ghosts[(MPIGlobal::rank + 1 + j) %
                                                 MPIGlobal::size]
                                                [si]->get_dQvec();
                        Vec dE = _ghosts[(MPIGlobal::rank + 1 + j) %
                                         MPIGlobal::size]
                                        [si]->get_delta_E();
                        MyMPI_Pack(&dQ[0], ndim_ + 2, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                        MyMPI_Pack(&dE[0], ndim_, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }
    }

    vector<MPI_Status> status((MPIGlobal::size - 1));
    MyMPI_Waitall((MPIGlobal::size - 1), &reqs[0], &status[0]);
    MyMPI_Barrier();
}

/**
  * @brief Communicate timesteps between processes for active particles.
  *
  * The timestep of a Particle is calculated on the "home" process of this
  * Particle (the process where the original Particle resides). But since the
  * timesteps need to be known on all processes where a ghost of the original
  * Particle is used, we need to communicate these timesteps to the respective
  * processes. Only active particles need to be considered. Inactive particles
  * are disabled by setting their id to 0.
  *
  * @param currentTime The current integer time of the simulation. An active
  * particle has an endtime equal to currentTime
  */
void DelTess::update_dts(unsigned long currentTime) {
    if(MPIGlobal::size < 2) {
        return;
    }

    vector<MPI_Request> reqs(MPIGlobal::size - 1, MPI_REQUEST_NULL);
    vector<unsigned int> buffers(MPIGlobal::size - 1);
    unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
    for(unsigned int i = 0; i < buffers.size(); i++) {
        buffers[i] = i * bufsize;
    }

    unsigned int maxsize = bufsize / sizeof(unsigned long);
    if(!(bufsize % sizeof(unsigned long))) {
        maxsize--;
    }

    vector<int> allflag(MPIGlobal::size - 1);
    vector<unsigned int> numsend(MPIGlobal::size - 1);
    vector<unsigned int> sendsizes(MPIGlobal::size - 1);
    int allflagsum = 0;
    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        sendsizes[i] =
                _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size].size();

        int send_pos = 0;
        for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]); si++) {
            unsigned long dt =
                    _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                            [si]->get_vorgen()
                                    ->get_particle()
                                    ->get_timestep();
            MyMPI_Pack(&dt, 1, MPI_UNSIGNED_LONG,
                       &MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);
        }
        allflag[i] = (sendsizes[i] <= maxsize);
        // if allflagsum equals MPIGlobal::size-1, all data has been sent
        allflagsum += allflag[i];
        numsend[i] = 1;
        // add continuation signal
        MyMPI_Pack(&allflag[i], 1, MPI_INT, &MPIGlobal::sendbuffer[buffers[i]],
                   bufsize, &send_pos);

        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos, MPI_PACKED,
                    (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0, &reqs[i]);
    }

    vector<unsigned int> numreceived(MPIGlobal::size - 1, 0);
    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        unsigned long dt = _exports[(MPIGlobal::rank + 1 + j) %
                                                    MPIGlobal::size]
                                                   [si]->get_vorgen()
                                                           ->get_particle()
                                                           ->get_timestep();
                        MyMPI_Pack(&dt, 1, MPI_UNSIGNED_LONG,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }

        int index;
        int tag;
        int recv_pos;

        // wait for a message from any source and with any tag to arrive
        MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
        // query the MPI_Status object to retrieve tag, source and size
        tag = status.MPI_TAG;
        index = status.MPI_SOURCE;
        int nelements;
        MyMPI_Get_count(&status, MPI_PACKED, &nelements);
        // select an index in the buffers to use for receiving and sending
        int freebuffer;
        if(index == MPIGlobal::size - 1) {
            freebuffer = MPIGlobal::rank;
        } else {
            freebuffer = index;
        }

        if(tag == 0) {
            MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                       MPI_PACKED, index, tag, &status);
            recv_pos = 0;
            unsigned int j = numreceived[freebuffer];
            while(recv_pos < nelements - 4) {
                unsigned long dt;
                MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                             nelements, &recv_pos, &dt, 1, MPI_UNSIGNED_LONG);
                // since the ghost particle was communicated with all of its
                // properties from the original process, it can already have
                // dQ and dE_grav contributions. However, we only want to record
                // the contributions from this process, so we have to reset them
                _ghosts[index][j]->reset_dQ();
                _ghosts[index][j]->reset_dE_grav();
                if(_ghosts[index][j]->get_endtime() == currentTime) {
                    _ghosts[index][j]->set_timestep(dt);
                    _ghosts[index][j]->set_id(2);
                } else {
                    _ghosts[index][j]->get_vorgen()->set_id(0);
                }
                j++;
            }
            numreceived[freebuffer] = j;
            int flag;
            MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                         &recv_pos, &flag, 1, MPI_INT);
            if(!flag) {
                i--;
            }
        }
    }

    while(allflagsum < MPIGlobal::size - 1) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        unsigned long dt = _exports[(MPIGlobal::rank + 1 + j) %
                                                    MPIGlobal::size]
                                                   [si]->get_vorgen()
                                                           ->get_particle()
                                                           ->get_timestep();
                        MyMPI_Pack(&dt, 1, MPI_UNSIGNED_LONG,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }
    }

    vector<MPI_Status> status((MPIGlobal::size - 1));
    MyMPI_Waitall((MPIGlobal::size - 1), &reqs[0], &status[0]);
    MyMPI_Barrier();
}

/**
 * @brief Communicate gravitational correction terms between processes
 *
 * The gravitational correction term for every GasParticle is calculated as a
 * part of the gravity tree walk. However, to add the correction to the
 * gravitational acceleration, we have to calculate a volume integral (converted
 * into a surface integral) over the corresponding Voronoi cell. This means we
 * need the gravitational correction terms for all neighbours, including those
 * that are ghosts imported from other MPI processes. Since the correction terms
 * are only up-to-date after the tree walk, communication needs to be done in a
 * separate routine.
 */
void DelTess::update_gravitational_corrections() {
    if(MPIGlobal::size < 2) {
        return;
    }

    vector<MPI_Request> reqs(MPIGlobal::size - 1, MPI_REQUEST_NULL);
    vector<unsigned int> buffers(MPIGlobal::size - 1);
    unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
    for(unsigned int i = 0; i < buffers.size(); i++) {
        buffers[i] = i * bufsize;
    }

    unsigned int maxsize = bufsize / sizeof(double);
    if(!(bufsize % sizeof(double))) {
        maxsize--;
    }

    vector<int> allflag(MPIGlobal::size - 1);
    vector<unsigned int> numsend(MPIGlobal::size - 1);
    vector<unsigned int> sendsizes(MPIGlobal::size - 1);
    int allflagsum = 0;
    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        sendsizes[i] =
                _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size].size();

        int send_pos = 0;
        for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]); si++) {
            double eta = _exports[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                 [si]->get_vorgen()
                                         ->get_particle()
                                         ->get_eta();
            MyMPI_Pack(&eta, 1, MPI_DOUBLE, &MPIGlobal::sendbuffer[buffers[i]],
                       bufsize, &send_pos);
        }
        allflag[i] = (sendsizes[i] <= maxsize);
        // if allflagsum equals MPIGlobal::size-1, all data has been sent
        allflagsum += allflag[i];
        numsend[i] = 1;
        // add continuation signal
        MyMPI_Pack(&allflag[i], 1, MPI_INT, &MPIGlobal::sendbuffer[buffers[i]],
                   bufsize, &send_pos);

        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos, MPI_PACKED,
                    (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0, &reqs[i]);
    }

    vector<unsigned int> numreceived(MPIGlobal::size - 1, 0);
    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        double eta = _exports[(MPIGlobal::rank + 1 + j) %
                                              MPIGlobal::size]
                                             [si]->get_vorgen()
                                                     ->get_particle()
                                                     ->get_eta();
                        MyMPI_Pack(&eta, 1, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }

        int index;
        int tag;
        int recv_pos;

        // wait for a message from any source and with any tag to arrive
        MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
        // query the MPI_Status object to retrieve tag, source and size
        tag = status.MPI_TAG;
        index = status.MPI_SOURCE;
        int nelements;
        MyMPI_Get_count(&status, MPI_PACKED, &nelements);
        // select an index in the buffers to use for receiving and sending
        int freebuffer;
        if(index == MPIGlobal::size - 1) {
            freebuffer = MPIGlobal::rank;
        } else {
            freebuffer = index;
        }

        if(tag == 0) {
            MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                       MPI_PACKED, index, tag, &status);
            recv_pos = 0;
            unsigned int j = numreceived[freebuffer];
            while(recv_pos < nelements - 4) {
                double eta;
                MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                             nelements, &recv_pos, &eta, 1, MPI_DOUBLE);
                _ghosts[index][j]->set_eta(eta);
                j++;
            }
            numreceived[freebuffer] = j;
            int flag;
            MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]], nelements,
                         &recv_pos, &flag, 1, MPI_INT);
            if(!flag) {
                i--;
            }
        }
    }

    while(allflagsum < MPIGlobal::size - 1) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                        si++) {
                        double eta = _exports[(MPIGlobal::rank + 1 + j) %
                                              MPIGlobal::size]
                                             [si]->get_vorgen()
                                                     ->get_particle()
                                                     ->get_eta();
                        MyMPI_Pack(&eta, 1, MPI_DOUBLE,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);
                    }
                    allflag[j] = (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                    allflagsum += allflag[j];
                    numsend[j]++;
                    // add continuation signal
                    MyMPI_Pack(&allflag[j], 1, MPI_INT,
                               &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                               &send_pos);

                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]], send_pos,
                                MPI_PACKED,
                                (MPIGlobal::rank + 1 + j) % MPIGlobal::size, 0,
                                &reqs[j]);
                }
            }
        }
    }

    vector<MPI_Status> status((MPIGlobal::size - 1));
    MyMPI_Waitall((MPIGlobal::size - 1), &reqs[0], &status[0]);
    MyMPI_Barrier();
}

/**
 * @brief Brute force check of the validity of the tesselation
 *
 * \attention This method has a runtime proportional to the number of simplices
 * times the number of points. It should only be used for small tesselations
 * and for testing purposes!
 *
 * Brute force check if the current tesselation obeys the properties of a
 * Delaunay tesselation, by subjecting all simplices to the in sphere test for
 * all points of the tesselation.
 *
 * Crashes if a point passes the in sphere test.
 */
void DelTess::check_tesselation() {
    for(unsigned int i = 1; i < _simplices.size(); i++) {
        if(_simplices[i] != NULL) {
            for(unsigned int j = 0; j < _points.size(); j++) {
                unsigned int* vorgens = _simplices[i]->get_vorgens();
#if ndim_ == 3
                VorGen* points[4] = {_points[vorgens[0]], _points[vorgens[1]],
                                     _points[vorgens[2]], _points[vorgens[3]]};
                bool check = _points[j] != points[0] &&
                             _points[j] != points[1] &&
                             _points[j] != points[2] && _points[j] != points[3];
#else
                VorGen* points[3] = {_points[vorgens[0]], _points[vorgens[1]],
                                     _points[vorgens[2]]};
                bool check = _points[j] != points[0] &&
                             _points[j] != points[1] && _points[j] != points[2];
#endif
                if(check) {
                    if(_simplices[i]->in_sphere(_points[j], _points)) {
                        cerr << "The Delauny Tesselation algorithm is "
                                "officially broken"
                             << endl;
                        my_exit();
                    }
                }
            }
        }
    }
}

/**
  * @brief Access the Simplex with the given index in the internal simplices
  * vector
  *
  * @param index The index of a Simplex in the internal simplices vector
  * @return A Simplex
  */
Simplex* DelTess::get_simplex(unsigned int index) { return _simplices[index]; }

#if ndim_ == 3
void DelTess::set_relations() {
    for(unsigned int i = 1; i < _simplices.size(); i++) {
        // it is possible that an entry in the list is empty (after a 3 to 2
        // flip, a NULL entry is created and if no new simplices are made
        // afterwards, it stays NULL)
        if(_simplices[i] != NULL) {
            unsigned int* vorgens = _simplices[i]->get_vorgens();
            VorGen* points[4] = {_points[vorgens[0]], _points[vorgens[1]],
                                 _points[vorgens[2]], _points[vorgens[3]]};
            for(unsigned int j = ndim_ + 1; j--;) {
                points[j]->add_tetrahedron(_simplices[i]);
            }
        }
    }
}
#else
/**
  * @brief Add the simplices of the tesselation to their respective vertices
  *
  * This is necessary to be able to access the simplices via the vertices during
  * the VorTess construction.
  */
void DelTess::set_relations() {
    for(unsigned int i = 1; i < _simplices.size(); i++) {
        unsigned int* vorgens = _simplices[i]->get_vorgens();
        VorGen* points[3] = {_points[vorgens[0]], _points[vorgens[1]],
                             _points[vorgens[2]]};
        for(unsigned int j = 0; j < 3; j++) {
            points[j]->add_tetrahedron(_simplices[i]);
        }
    }
}
#endif

/**
 * @brief Get the total number of points in the tesselation
 *
 * @return The total number of points in the tesselation
 */
unsigned int DelTess::get_size() {
#if ndim_ == 2
    return _points.size() - _mirrors.size() - 3;
#else
    return _points.size() - _mirrors.size() - 4;
#endif
}

/**
 * @brief Get the triangles that make up the tesselation
 *
 * @param positions Positions of the vertices of the triangles
 * @param connectivity Connections between the positions that form the actual
 * triangles
 */
void DelTess::get_triangles(vector<float>& positions,
                            vector<int>& connectivity) {
    for(unsigned int i = 0; i < _simplices.size(); i++) {
        bool ghost = false;
        vector<float> to_add;
        for(unsigned int j = 0; j < 3; j++) {
            unsigned int v = _simplices[i]->vorgen(j);
            if(v >= _points.size() - _mirrors.size() - 3) {
                ghost = true;
                break;
            }
            VorGen* point = _points[v];
            to_add.push_back(point->x());
            to_add.push_back(point->y());
            to_add.push_back(0.);
        }
        if(!ghost) {
            for(unsigned int j = 0; j < to_add.size(); j++) {
                positions.push_back(to_add[j]);
            }
            connectivity.push_back(3);
            connectivity.push_back(positions.size() / 3 - 3);
            connectivity.push_back(positions.size() / 3 - 2);
            connectivity.push_back(positions.size() / 3 - 1);
        }
    }
}
