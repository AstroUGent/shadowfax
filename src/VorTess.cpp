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
 * @file VorTess.cpp
 *
 * @brief Voronoi tesselation for the old algorithm: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "VorTess.hpp"
#include "DelTess.hpp"                // for DelTess
#include "ProgramLog.hpp"             // for LOGS
#include "Simplex.hpp"                // for Simplex
#include "VorCell.hpp"                // for VorCell
#include "VorFace.hpp"                // for VorFace
#include "VorGen.hpp"                 // for VorGen
#include "utilities/GasParticle.hpp"  // for GasParticle
#include <algorithm>                  // for sort
#include <cstdlib>                    // for NULL
#include <list>                       // for list, _List_iterator
#include <ostream>                    // for operator<<, stringstream, etc
#include <sstream>
#include <vector>  // for vector
#if ndim_ == 2
#include <string>  // for operator<<
#endif
using namespace std;

/**
 * @brief Constructor
 *
 * @param delcont DelCont specifying the dimensions of the simulation box
 * @param numpart Number of particles that will be added to the tesselation
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective (false)
 * @param tolerance Tolerance value used to distinguish between standard
 * floating point arithmetics and arbitrary precision arithmetics in geometry
 * tests
 */
VorTess::VorTess(DelCont* delcont, unsigned int numpart, bool periodic,
                 double tolerance) {
    _delaunay = new DelTess(delcont, numpart, periodic, tolerance);
    _delaunay->add_voronoi_tesselation(this);

    LOGS("VorTess constructed");
}

/**
 * @brief Destructor
 *
 * Remove all cells and points to ensure proper garbage collection.
 */
VorTess::~VorTess() {
    for(unsigned int i = 0; i < _cells.size(); i++) {
        delete _cells[i];
    }
    _cells.clear();
    for(unsigned int i = _faces.size(); i--;) {
        delete _faces[i];
    }
    _faces.clear();
    delete _delaunay;

    LOGS("VorTess destructed");
}

/**
 * @brief Print the tesselation (in ascii) to the given stream
 *
 * Points of the tesselation are outputted as
 * p:\\t(c1,c2,c3)\n (with the ci doubles)
 * lines are outputted as
 * l:\\t(c1,c2,c3)\\t(d1,d2,d3)\n
 * (for 2D, the respective last coordinates vanish in these expressions)
 *
 * This method calls the print method of the VorCells.
 *
 * @param stream std::ostream to write to
 */
void VorTess::print_tesselation(ostream& stream) {
    //    stream << "Outputting Voronoi tesselation:\n";
    for(unsigned int i = 0; i < _cells.size(); i++) {
        _cells[i]->print(stream);
    }
    stream << flush;
}

/**
 * @brief Print the tesselation to the given stream in a format that can be
 * plotted using gnuplot
 *
 * @param stream std::ostream to write to
 */
void VorTess::print_tesselation_gnuplot(ostream& stream) {
    for(unsigned int i = 0; i < _cells.size(); i++) {
        _cells[i]->print_gnuplot(stream, i);
    }
    stream << flush;
}

/**
 * @brief Access the VorCells of the tesselation
 *
 * This function is not used for the moment.
 *
 * @return std::vector containing the cells of the tesselation
 */
vector<VorCell*> VorTess::get_cells() {
    return _cells;
}

/**
 * @brief Add the given GasParticle with the given index to the tesselation
 *
 * The same index can later be used to retrieve the corresponding cell from the
 * tesselation.
 *
 * @param part GasParticle to add
 * @param index Index of the particle in the internal lists, can be used to
 * retrieve it later
 */
void VorTess::add_point(GasParticle* part, unsigned int index) {
    _delaunay->add_particle(part, index);
}

/**
 * @brief Construct the tesselation
 *
 * Core function of this class.
 *
 * A walk is performed on the points of the DelTess. With each point a VorCell
 * is associated. Then for every axis determined by this point and a
 * neighbouring point (that is part of a neighbouring tetrahedron) a VorFace is
 * constructed (by rotating around the axis), which is then added to the cell.
 */
#if ndim_ == 3
void VorTess::construct() {
    LOGS("Starting VorTess construction");
    _delaunay->set_relations();
    vector<VorGen*> points = _delaunay->get_points();
    for(unsigned int i = 0; i < _delaunay->get_size(); i++) {
        if(!points[i]) {
            continue;
        }
        _cells.push_back(new VorCell(points[i]));
        points[i]->set_cell(_cells.back());
        list<Simplex*> tetrahedra = points[i]->get_tetrahedra();
        Simplex* current_tetrahedron;
        for(list<Simplex*>::iterator it = tetrahedra.begin();
            it != tetrahedra.end(); it++) {
            current_tetrahedron = *it;
            unsigned int* vorgens2 = current_tetrahedron->get_vorgens();
            VorGen* points2[4] = {points[vorgens2[0]], points[vorgens2[1]],
                                  points[vorgens2[2]], points[vorgens2[3]]};
            unsigned int k = current_tetrahedron->get_index(i);
            // k now contains the index of points[i] in points2
            for(int l = 0; l < 3; l++) {
                if(!points2[(k + l + 1) % 4]->flagged()) {
                    points2[(k + l + 1) % 4]->flag();
                    _cells.back()->add_ngb(points2[(k + l + 1) % 4]);
                    _cells.back()->add_ngb_id(vorgens2[(k + l + 1) % 4] - 4);
                    VorFace* face;
                    if(points2[(k + l + 1) % 4]->get_cell() == NULL) {
                        _faces.push_back(new VorFace(i, points));
                        _cells.back()->add_face(_faces.back());
                        face = _faces.back();
                        face->add_ngb(points2[(k + l + 1) % 4]);
                        face->add_ngb_id(
                                points2[(k + l + 1) % 4]->get_particle_id());
                        VorGen* axis = points2[(k + l + 1) % 4];
                        VorGen* start = points2[(k + (l + 1) % 3 + 1) % 4];
                        VorGen* other = points2[(k + (l + 2) % 3 + 1) % 4];
                        unsigned int vstart =
                                vorgens2[(k + (l + 1) % 3 + 1) % 4];
                        unsigned int vother =
                                vorgens2[(k + (l + 2) % 3 + 1) % 4];
                        Simplex* next_tetrahedron = _delaunay->get_simplex(
                                current_tetrahedron->get_ngb(
                                        (k + (l + 1) % 3 + 1) % 4));
                        unsigned int* next_vorgens =
                                next_tetrahedron->get_vorgens();
                        VorGen* next_points[4] = {points[next_vorgens[0]],
                                                  points[next_vorgens[1]],
                                                  points[next_vorgens[2]],
                                                  points[next_vorgens[3]]};
                        unsigned int next_index =
                                current_tetrahedron->get_ngbface(
                                        (k + (l + 1) % 3 + 1) % 4);
                        VorGen* next_point = next_points[next_index];
                        unsigned int vnext = next_vorgens[next_index];
                        start = other;
                        vstart = vother;
                        VorGen* first_special =
                                next_tetrahedron->get_special_point(points);
                        VorGen* prev_special = NULL;
                        unsigned int count = 0;
                        while(next_point != other) {
                            VorGen* next_special =
                                    next_tetrahedron->get_special_point(points);
                            if(prev_special == NULL ||
                               (next_special->distance(*prev_special) >
                                        1.e-13 &&
                                next_special->distance(*first_special) >
                                        1.e-13)) {
                                face->add_vertex(next_special);
                                count++;
                            }
                            face->add_facengb(start);
                            prev_special = next_special;
                            next_tetrahedron = _delaunay->get_simplex(
                                    next_tetrahedron->get_ngb_from_vorgen(
                                            vstart));
                            next_vorgens = next_tetrahedron->get_vorgens();
                            next_points[0] = points[next_vorgens[0]];
                            next_points[1] = points[next_vorgens[1]];
                            next_points[2] = points[next_vorgens[2]];
                            next_points[3] = points[next_vorgens[3]];
                            next_index = 0;
                            while(next_points[next_index] == axis ||
                                  next_points[next_index] == points2[k] ||
                                  next_points[next_index] == next_point) {
                                next_index++;
                            }
                            start = next_point;
                            vstart = vnext;
                            next_point = next_points[next_index];
                            vnext = next_vorgens[next_index];
                        }
                        VorGen* next_special =
                                next_tetrahedron->get_special_point(points);
                        if(next_special->distance(*prev_special) > 1.e-13 &&
                           next_special->distance(*first_special) > 1.e-13) {
                            face->add_vertex(next_special);
                            count++;
                        }
                        face->add_facengb(start);
                        // a face with less than 3 vertices can in principle
                        // exist, but in this case, we don't use it during the
                        // remainder of the calculation. We therefore set its
                        // area to 0
                        if(count < 3) {
                            face->set_area(0.);
                        } else {
                            // the order does matter: midpoint is used in area
                            // calculation!
                            face->calculate_midpoint();
                            face->calculate_area();
                        }
                    } else {
                        face = points2[(k + l + 1) % 4]->get_cell()->get_face(
                                points[i]);
                        _cells.back()->add_face(face);
                    }
                }
            }
        }
        // clear flags
        vector<VorGen*> ngbs = _cells.back()->get_ngbs();
        for(unsigned int j = 0; j < ngbs.size(); j++) {
            ngbs[j]->unflag();
        }

        _cells.back()->calculate_centroid();
        _cells.back()->calculate_volume();
    }
    LOGS("VorTess construction ready");
}
#else
void VorTess::construct() {
    LOGS("Starting VorTess construction");
    _delaunay->set_relations();
    vector<VorGen*> points = _delaunay->get_points();
    for(unsigned int i = 0; i < _delaunay->get_size(); i++) {
        if(!points[i]) {
            continue;
        }
        _cells.push_back(new VorCell(points[i]));
        points[i]->set_cell(_cells.back());
        list<Simplex*> tetrahedra = points[i]->get_tetrahedra();
        // walk around the triangles in counterclockwise order
        Simplex* current = *(tetrahedra.begin());
        // iindex is the index of the central point in the current triangle
        int iindex = current->get_index(i);
        unsigned int* cvorgens = current->get_vorgens();
        VorGen* cpoints[3] = {points[cvorgens[0]], points[cvorgens[1]],
                              points[cvorgens[2]]};
        // start is the first point in counterclockwise order after the central
        // point
        VorGen* start = cpoints[(iindex + 1) % 3];
        // next is the second point in counterclockwise order after the central
        // point
        // it will be the first point after the central point in the next
        // triangle
        VorGen* next = cpoints[(iindex + 2) % 3];
        unsigned int vnext = cvorgens[(iindex + 2) % 3];
        // prevsp is the midpoint of the circumsphere of the current triangle
        VorGen* prevsp = current->get_special_point(points);
        // firstsp is the first midpoint, to be used in the end to connect the
        // final edge
        VorGen* firstsp = prevsp;
        // cursp is the midpoint of the circumsphere of the next triangle, which
        // becomes the current triangle at the start of every iteration
        VorGen* cursp = NULL;
        // iterate until next == start, at this point, we have completed a full
        // walk-around
        while(next != start) {
            // use the information stored in the current triangle to retrieve
            // absolute index information for the next triangle without
            // searching for points
            int nindex = current->get_ngbindex((iindex + 1) % 3);
            // go to the next triangle
            current =
                    _delaunay->get_simplex(current->get_ngb((iindex + 1) % 3));
            cursp = current->get_special_point(points);
            //      double distance = fabs(cursp->distance(*prevsp));
            VorFace* face;
            // every face is only constructed once
            // if the point on the other side does not have a cell
            // associated with it, then it is either a ghost point
            // or it is not yet processed
            // either way, this point (and cell) becomes the owner
            // of the face
            if(next->get_cell() == NULL) {
                _faces.push_back(new VorFace(i, points));
                face = _faces.back();
                _cells.back()->add_face(face);
                face->add_vertex(prevsp);
                face->add_vertex(cursp);
                face->add_ngb(next);
                face->add_ngb_id(next->get_particle_id());
                face->calculate_area();
                // try to solve problems with accuracy
                if(prevsp->distance(*cursp) < 1.e-13) {
                    face->set_area(0.);
                }
                face->calculate_midpoint();
            } else {
                face = next->get_cell()->get_face(points[i]);
                _cells.back()->add_face(face);
            }
            // store the neighbour relations
            _cells.back()->add_ngb(next);
            _cells.back()->add_ngb_id(vnext);
            cvorgens = current->get_vorgens();
            cpoints[0] = points[cvorgens[0]];
            cpoints[1] = points[cvorgens[1]];
            cpoints[2] = points[cvorgens[2]];
            // update loop variables
            next = cpoints[nindex];
            vnext = cvorgens[nindex];
            iindex = (nindex + 1) % 3;
            prevsp = cursp;
        }
        VorFace* face;
        // construct the final face
        if(next->get_cell() == NULL) {
            _faces.push_back(new VorFace(i, points));
            face = _faces.back();
            _cells.back()->add_face(face);
            face->add_vertex(prevsp);
            face->add_vertex(firstsp);
            face->add_ngb(next);
            face->add_ngb_id(next->get_particle_id());
            face->calculate_area();
            // try to solve problems with accuracy
            if(prevsp->distance(*firstsp) < 1.e-13) {
                face->set_area(0.);
            }
            face->calculate_midpoint();
        } else {
            face = next->get_cell()->get_face(points[i]);
            _cells.back()->add_face(face);
        }
        _cells.back()->add_ngb(next);
        _cells.back()->add_ngb_id(vnext);
        _cells.back()->calculate_centroid();
        _cells.back()->calculate_volume();
    }
    LOGS("VorTess construction ready");
}
#endif

/**
 * @brief Do the hydrodynamical integration on the mesh
 *
 * For every face, we first calculate its velocity and then the flux that passes
 * through.
 *
 * @param timeline TimeLine of the simulation
 * @param solver RiemannSolver used to solve the Riemann problem at the faces
 * @param particles ParticleVector containing the particles of the simulation
 */
#ifndef ICMAKER
void VorTess::hydro(TimeLine& timeline, RiemannSolver& solver,
                    ParticleVector& particles) {
    LOGS("Starting flux calculation");
    for(unsigned int i = _faces.size(); i--;) {
        _faces[i]->set_v();
        _faces[i]->calculate_flux(timeline, solver);
    }
    LOGS("Flux calculation done");
}
#endif

/**
 * @brief Advect the hydrodynamical quantities for the given timestep
 *
 * @warning Does not work!
 *
 * @param dt Timestep
 */
void VorTess::advect(double dt) {
    for(unsigned int i = _faces.size(); i--;) {
        _faces[i]->calculate_advection_flux(dt);
    }
    for(unsigned int i = _cells.size(); i--;) {
        _cells[i]->get_particle()->update_Q();
    }
}

/**
 * @brief Calculate and output the total volume of all cells in the tesselation
 * to the given stream
 *
 * @param stream std::ostream to write to
 */
void VorTess::print_cell_statistics(ostream& stream) {
    double totvol = 0.;
    double totsize = 0.;
    for(unsigned int i = 0; i < _cells.size(); i++) {
        totvol += _cells[i]->get_volume();
        totsize += _cells[i]->get_h();
    }
    double avgvolume = totvol / _cells.size();
    double avgsize = totsize / _cells.size();
    double sigmavol = 0.;
    double sigmasize = 0.;
    double diff;
    for(unsigned int i = 0; i < _cells.size(); i++) {
        diff = _cells[i]->get_volume() - avgvolume;
        sigmavol += diff * diff;
        diff = _cells[i]->get_h() - avgsize;
        sigmasize += diff * diff;
    }
    sigmavol = sqrt(sigmavol / _cells.size());
    sigmasize = sqrt(sigmasize / _cells.size());
    stream << "Total volume: " << totvol << endl;
    stream << "Average cell volume: " << avgvolume << " +- " << sigmavol
           << endl;
    stream << "Average cell size: " << avgsize << " +- " << sigmasize << endl;
}

/**
 * @brief Print the Delaunay tesselation corresponding to this tesselation to
 * the given stream
 *
 * @param stream std::ostream to write to
 */
void VorTess::print_delaunay(ostream& stream) {
    _delaunay->output_tesselation(stream);
}

/**
 * @brief Print the tesselation to the given stream in a format that can be read
 * by POVRAY
 *
 * @param stream std::ostream to write to
 */
void VorTess::print_tesselation_pov(ostream& stream) {
    for(unsigned int i = _cells.size(); i--;) {
        _cells[i]->print_pov(stream);
    }
}

#ifndef ICMAKER
/**
 * @brief Print the tesselation to the given stream in a format that can be
 * plotted using leaflet
 *
 * @param colormap ColorMap used to convert values to colors
 * @param maxW Maximal primitive variables Statevector
 * @param minW Minimal primitive variables StateVector
 */
void VorTess::print_tesselation_leaflet(ColorMap* colormap, StateVector maxW,
                                        StateVector minW) {
    double tilesize = 1. / 64.;
    ofstream ofiles[256];
    for(unsigned int xblock = 0; xblock < 4; xblock++) {
        for(unsigned int yblock = 0; yblock < 4; yblock++) {
            for(unsigned int i = xblock * 16; i < (xblock + 1) * 16; i++) {
                for(unsigned int j = yblock * 16; j < (yblock + 1) * 16; j++) {
                    stringstream name;
                    name << "grid/x";
                    name.fill('0');
                    name.width(2);
                    name << i;
                    name << "y";
                    name.fill('0');
                    name.width(2);
                    name << (63 - j);
                    name << ".dat";
                    ofiles[(i - xblock * 16) * 16 + (j - yblock * 16)].open(
                            name.str().c_str(), ios::out | ios::binary);
                }
            }
            for(unsigned int i = _cells.size(); i--;) {
                double box[4] = {xblock * 16 * tilesize, yblock * 16 * tilesize,
                                 (xblock + 1) * 16 * tilesize,
                                 (yblock + 1) * 16 * tilesize};
                if(!_cells[i]->overlap(box)) {
                    continue;
                }
                for(unsigned int ix = xblock * 16; ix < (xblock + 1) * 16;
                    ix++) {
                    for(unsigned int iy = yblock * 16; iy < (yblock + 1) * 16;
                        iy++) {
                        double smallbox[4] = {ix * tilesize, iy * tilesize,
                                              (ix + 1) * tilesize,
                                              (iy + 1) * tilesize};
                        if(_cells[i]->overlap(smallbox)) {
                            _cells[i]->print_leaflet(
                                    ofiles[(ix - xblock * 16) * 16 +
                                           (iy - yblock * 16)],
                                    ix * 256, iy * 256, colormap, maxW, minW);
                        }
                    }
                }
            }
            for(int i = 0; i < 256; i++) {
                ofiles[i].close();
            }
        }
    }
}
#endif

/**
 * @brief Print the tesselation to the given stream in VTK-format
 *
 * @param stream std::ostream to write to
 */
void VorTess::print_tesselation_vtk(ostream& stream) {
#if ndim_ == 2
    stream << "# vtk DataFile Version 3.0\n";
    stream << "vtk output\n";
    stream << "ASCII\n";
    stream << "DATASET POLYDATA\n";
    stringstream points;
    stringstream polygons;
    stringstream density;
    unsigned int numpoints = 0;
    unsigned int numpolygons = 0;
    unsigned int numconnections = 0;
    unsigned int numcells = _cells.size();
    for(unsigned int i = numcells; i--;) {
        _cells[i]->print_vtk(points, numpoints, polygons, numpolygons,
                             numconnections, density);
    }
    stream << "POINTS " << numpoints << " float\n";
    stream << points.str();
    stream << "POLYGONS " << numpolygons << " " << numconnections << "\n";
    stream << polygons.str();
    stream << "CELL_DATA " << numcells << "\n";
    stream << "SCALARS Density float\n";
    stream << "LOOKUP_TABLE default\n";
    stream << density.str();
    stream << endl;
#else
    // this code writes a dummy VTK file containing two entangled tetrahedra
    stream << "# vtk DataFile Version 2.0\n";
    stream << "some label\n";
    stream << "BINARY\n";
    stream << "DATASET UNSTRUCTURED_GRID\n";
    stream << "POINTS 5 float\n";
    float positions[15] = {0., 0., 0., 1., 0., 0., 0., 1.,
                           0., 0., 0., 1., 1., 1., 1.};
    char* bytes = reinterpret_cast<char*>(positions);
    for(unsigned int i = 0; i < 15; i++) {
        char fbytes[4];
        fbytes[0] = bytes[4 * i];
        fbytes[1] = bytes[4 * i + 1];
        fbytes[2] = bytes[4 * i + 2];
        fbytes[3] = bytes[4 * i + 3];
        bytes[4 * i] = fbytes[3];
        bytes[4 * i + 1] = fbytes[2];
        bytes[4 * i + 2] = fbytes[1];
        bytes[4 * i + 3] = fbytes[0];
    }
    stream.write(bytes, 15 * sizeof(float));
    stream << "CELLS 2 10\n";
    int connectivity[10] = {4, 0, 1, 2, 3, 4, 0, 1, 2, 4};
    bytes = reinterpret_cast<char*>(connectivity);
    for(unsigned int i = 0; i < 10; i++) {
        char fbytes[4];
        fbytes[0] = bytes[4 * i];
        fbytes[1] = bytes[4 * i + 1];
        fbytes[2] = bytes[4 * i + 2];
        fbytes[3] = bytes[4 * i + 3];
        bytes[4 * i] = fbytes[3];
        bytes[4 * i + 1] = fbytes[2];
        bytes[4 * i + 2] = fbytes[1];
        bytes[4 * i + 3] = fbytes[0];
    }
    stream.write(bytes, 10 * sizeof(int));
    stream << "CELL_TYPES 2\n";
    int celltype[2] = {10, 10};
    bytes = reinterpret_cast<char*>(celltype);
    for(unsigned int i = 0; i < 2; i++) {
        char fbytes[4];
        fbytes[0] = bytes[4 * i];
        fbytes[1] = bytes[4 * i + 1];
        fbytes[2] = bytes[4 * i + 2];
        fbytes[3] = bytes[4 * i + 3];
        bytes[4 * i] = fbytes[3];
        bytes[4 * i + 1] = fbytes[2];
        bytes[4 * i + 2] = fbytes[1];
        bytes[4 * i + 3] = fbytes[0];
    }
    stream.write(bytes, 2 * sizeof(int));
    stream << "CELL_DATA 2\n";
    stream << "SCALARS name float\n";
    stream << "LOOKUP_TABLE default\n";
    float data[2] = {0., 1.};
    bytes = reinterpret_cast<char*>(data);
    for(unsigned int i = 0; i < 2; i++) {
        char fbytes[4];
        fbytes[0] = bytes[4 * i];
        fbytes[1] = bytes[4 * i + 1];
        fbytes[2] = bytes[4 * i + 2];
        fbytes[3] = bytes[4 * i + 3];
        bytes[4 * i] = fbytes[3];
        bytes[4 * i + 1] = fbytes[2];
        bytes[4 * i + 2] = fbytes[1];
        bytes[4 * i + 3] = fbytes[0];
    }
    stream.write(bytes, 2 * sizeof(float));
    stream << "POINT_DATA 5\n";
#endif
}

/**
 * @brief Calculate the characteristic lengths of all cells
 *
 * The characteristic length of a cell is the radius of a sphere/circle with the
 * same volume/face area as the cell.
 */
void VorTess::set_hs() {
    for(unsigned int i = _cells.size(); i--;) {
        _cells[i]->set_h();
    }
}

/**
 * @brief Finalize the Delaunay tesselation by adding mirror points and MPI
 * ghosts
 *
 * @param parttree Tree for the particle set
 */
void VorTess::complete(Tree& parttree) {
    _delaunay->add_mirrors(parttree);
}

/**
 * @brief Communicate primitive variables between MPI processes
 */
void VorTess::update_Ws() {
    _delaunay->update_Ws();
}

/**
 * @brief Communicate gradients between MPI processes
 */
void VorTess::update_gradients() {
    _delaunay->update_gradients();
}

/**
 * @brief Communicate fluxes between MPI processes
 */
void VorTess::update_dQs() {
    _delaunay->update_dQs();
}

/**
 * @brief Communicate timesteps between MPI processes
 *
 * @param currentTime Current integer time of the simulation
 */
void VorTess::update_dts(unsigned long currentTime) {
    _delaunay->update_dts(currentTime);
}

/**
 * @brief Communicate gravitational correction terms between MPI processes
 */
void VorTess::update_gravitational_corrections() {
    // we first make sure all local factors are added to eta
    for(unsigned int i = 0; i < _cells.size(); i++) {
        _cells[i]->finalize_eta();
    }
    _delaunay->update_gravitational_corrections();
}

/**
 * @brief Get the cell with the given index
 *
 * @param index Index with which the particle corresponding to the cell was
 * added to the tesselation
 * @return VorCell corresponding to the given index
 */
VorCell* VorTess::get_cell(unsigned int index) {
    return _cells[index];
}

/**
 * @brief Check if the Delaunay tesselation fulfils the geometrical criterion
 *
 * @warning This method is computationally very expensive and should only be
 * used for debugging and for small point sets!
 */
void VorTess::check_delaunay() {
    _delaunay->check_tesselation();
}

/**
 * @brief Get the faces of the tesselation
 *
 * @return Reference to the face list
 */
vector<VorFace*>& VorTess::get_faces() {
    return _faces;
}

/**
 * @brief Get the triangles that make up the cells of the tesselation
 *
 * @param positions List of positions to add to
 * @param connectivity List with connectivity information to add to
 * @param data List with particle data to add to
 */
void VorTess::get_triangles(std::vector<float>& positions,
                            std::vector<int>& connectivity,
                            std::vector<StateVector>& data) {
    for(unsigned int i = 0; i < _cells.size(); i++) {
        _cells[i]->get_triangles(positions, connectivity, data);
    }
}
/**
 * @brief Get the triangles that make up the Delaunay tesselation
 *
 * @param positions Positions of the vertices of the triangles
 * @param connectivity Connections between the positions that form the actual
 * triangles
 */
void VorTess::get_delaunay_triangles(vector<float>& positions,
                                     vector<int>& connectivity) {
    _delaunay->get_triangles(positions, connectivity);
}

/**
 * @brief Get the indices of the neighbours of the given cell
 *
 * For debugging purposes, to compare tesselations build by different
 * algorithms.
 *
 * @param cell VorCell for which we want the neighbours
 * @return List of neighbour indices
 */
vector<unsigned long> VorTess::get_ngb_ids(VorCell* cell) {
    vector<unsigned long> ngb_ids;
    vector<VorGen*> ngbs = cell->get_ngbs();
    for(unsigned int i = 0; i < ngbs.size(); i++) {
        if(ngbs[i]->get_particle()) {
            ngb_ids.push_back(ngbs[i]->get_particle()->id());
        } else {
            VorCell* original = _cells[ngbs[i]->get_original()];
            ngb_ids.push_back(original->get_particle()->id());
        }
    }
    sort(ngb_ids.begin(), ngb_ids.end());
    return ngb_ids;
}
