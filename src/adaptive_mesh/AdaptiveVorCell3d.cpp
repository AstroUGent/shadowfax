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
 * @file AdaptiveVorCell3d.cpp
 *
 * @brief 3D mesh evolution Voronoi cell: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdaptiveVorCell3d.hpp"
#include "AdaptiveMeshException.hpp"  // for AdaptiveMeshException
#include "AdaptiveMeshUtils.hpp"      // for get_wallpos, etc
#include "AdaptiveVorFace3d.hpp"      // for AdaptiveVorFace3d
#include "Error.hpp"                  // for my_exit
#include "RestartFile.hpp"            // for RestartFile
#include "StateVector.hpp"            // for StateVector, operator*, etc
#include "predicates.hpp"             // for insphere_old, etc
#include "utilities/Cuboid.hpp"       // for Cuboid
#include "utilities/GasParticle.hpp"  // for GasParticle
#include <cmath>                      // for cbrt, fabs, M_PI
#include <cstdlib>                    // for exit, NULL
#include <ext/alloc_traits.h>
#include <iostream>  // for cerr, cout
#include <sstream>
using namespace std;

//#define VERBOSE

/**
 * @brief Calculate the midpoint of the sphere through points A, B, C and D
 *
 * If no valid midpoint could be calculated, the code will abort.
 *
 * @param a Coordinates of point A
 * @param b Coordinates of point B
 * @param c Coordinates of point C
 * @param d Coordinates of point D
 * @param vert Array to store the result in
 * @return True if a valid midpoint was calculated, false otherwise
 */
bool AdaptiveVorCell3d::get_vertexpoint(double* a, double* b, double* c,
                                        double* d, double* vert) {
    double r1[3], r2[3], r3[3];
    r1[0] = b[0] - a[0];
    r1[1] = b[1] - a[1];
    r1[2] = b[2] - a[2];
    r2[0] = c[0] - a[0];
    r2[1] = c[1] - a[1];
    r2[2] = c[2] - a[2];
    r3[0] = d[0] - a[0];
    r3[1] = d[1] - a[1];
    r3[2] = d[2] - a[2];
    double fac1, fac2, fac3;
    fac1 = r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2];
    fac2 = r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2];
    fac3 = r3[0] * r3[0] + r3[1] * r3[1] + r3[2] * r3[2];
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
    double Vinv = 1. / V;
    vert[0] = a[0] + R[0] * Vinv;
    vert[1] = a[1] + R[1] * Vinv;
    vert[2] = a[2] + R[2] * Vinv;
    if(vert[0] != vert[0] || vert[1] != vert[1] || vert[2] != vert[2]) {
        throw AdaptiveMeshException();
        cerr << "NaN vertex!" << endl;
        cerr << "V: " << V << endl;
        cerr << a[0] << "\t" << a[1] << "\t" << a[2] << "\n";
        cerr << b[0] << "\t" << b[1] << "\t" << b[2] << "\n";
        cerr << c[0] << "\t" << c[1] << "\t" << c[2] << "\n";
        cerr << d[0] << "\t" << d[1] << "\t" << d[2] << "\n";
        vert[0] = 0.25 * (a[0] + b[0] + c[0] + d[0]);
        vert[1] = 0.25 * (a[1] + b[1] + c[1] + d[1]);
        vert[2] = 0.25 * (a[2] + b[2] + c[2] + d[2]);
        cerr << vert[0] << "\t" << vert[1] << "\t" << vert[2] << endl;
        double orient = predicates::orient3d_old(a, b, c, d);
        cerr << "orient = " << orient << endl;
        if(!orient) {
            cerr << "Exactly 0!" << endl;
        }
        my_exit();
        return false;
    }
    return true;
}

// check if point is inside the sphere through a, b, c and d
/**
 * @brief Check if point E is inside the sphere through points A, B, C and D
 *
 * @param a Coordinates of point A
 * @param b Coordinates of point B
 * @param c Coordinates of point C
 * @param d Coordinates of point D
 * @param point Coordinates of point E
 * @return True if E is inside the circumsphere through A, B, C, and D, false
 * otherwise
 */
bool AdaptiveVorCell3d::inside(double* a, double* b, double* c, double* d,
                               double* point) {
    return predicates::insphere_old(a, b, c, d, point) > 0.;
}

/**
 * @brief Restart constructor. Initialize a cell from the given RestartFile
 *
 * @param rfile RestartFile to read from
 * @param particle GasParticle associated with the cell (not used)
 */
AdaptiveVorCell3d::AdaptiveVorCell3d(RestartFile& rfile,
                                     GasParticle* particle) {
    rfile.read(_ngbs);
    rfile.read(_ngbwalls);
    rfile.read(_facengbs);
    rfile.read(_facengbwalls);
    unsigned int vsize;
    rfile.read(vsize);
    _faces.resize(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        _faces[i] = new AdaptiveVorFace3d(rfile);
    }
    rfile.read(_pos, 3);
    rfile.read(_valid_pos, 3);
    rfile.read(_id);
    _particle = particle;
    _flag1 = false;
    _flag2 = false;
    _active = false;
}

/**
 * @brief Dump the cell to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void AdaptiveVorCell3d::dump(RestartFile& rfile) {
    rfile.write(_ngbs);
    rfile.write(_ngbwalls);
    rfile.write(_facengbs);
    rfile.write(_facengbwalls);
    unsigned int vsize = _faces.size();
    rfile.write(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        _faces[i]->dump(rfile);
    }
    rfile.write(_pos, 3);
    rfile.write(_valid_pos, 3);
    rfile.write(_id);
}

/**
 * @brief Constructor.
 *
 * @param pos Coordinates of the generator of the cell
 * @param id Index of the cell in the AdaptiveVorTess3d
 * @param particle GasParticle associated with the cell
 */
AdaptiveVorCell3d::AdaptiveVorCell3d(double* pos, unsigned int id,
                                     GasParticle* particle) {
    _pos[0] = pos[0];
    _pos[1] = pos[1];
    _pos[2] = pos[2];
    _valid_pos[0] = pos[0];
    _valid_pos[1] = pos[1];
    _valid_pos[2] = pos[2];
    _id = id;
    _particle = particle;
    _flag1 = false;
    _flag2 = false;
    _active = false;
    _flip_info[0] = 0;
    _flip_info[1] = 0;
}

/**
 * @brief Free the memory used by the faces
 */
void AdaptiveVorCell3d::cleanup_faces() {
    for(unsigned int i = 0; i < _faces.size(); i++) {
        delete _faces[i];
    }
    _faces.clear();
}

/**
 * @brief Destructor. Clean up the faces
 */
AdaptiveVorCell3d::~AdaptiveVorCell3d() {
    cleanup_faces();
}

/**
 * @brief Add the given neighbour to the end of the list
 *
 * Also creates appropriate fields in the other internal lists to guarantee that
 * all lists always have the same size.
 *
 * @param id Integer index of the new neighbour in the AdaptiveVorTess3d
 * @param wall Boundary flag for the new neighbour
 * @return Unsigned integer index of the neighbour in the internal neighbour
 * list
 */
unsigned int AdaptiveVorCell3d::add_ngb(int id, int wall) {
    _ngbs.push_back(id);
    _ngbwalls.push_back(wall);
    vector<int> facengbs;
    _facengbs.push_back(facengbs);
    vector<int> facengbwalls;
    _facengbwalls.push_back(facengbwalls);
    return _ngbs.size() - 1;
}

/**
 * @brief Add a facengb to the end of the list corresponding to the last added
 * neighbour
 *
 * Also adds an appropriate wall key, to guarantee both lists always have the
 * same size
 *
 * @param id Integer index of the new facengb in the AdaptiveVorTess3d
 * @param wall Boundary flag for the new facengb
 */
void AdaptiveVorCell3d::add_facengb(int id, int wall) {
    _facengbs.back().push_back(id);
    _facengbwalls.back().push_back(wall);
}

/**
 * @brief Debug method: print cell information to the given stream
 *
 * @param stream std::ostream to write to
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost list
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell3d::print_relevant_info(ostream& stream,
                                            vector<AdaptiveVorCell3d*>& cells,
                                            vector<AdaptiveVorCell3d*>& ghosts,
                                            Cuboid& box) {
    stream << "position: " << _pos[0] << "\t" << _pos[1] << "\t" << _pos[2]
           << "\n";
    stream << "_ngbs.size() = " << _ngbs.size() << "\n";
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        double pos[3];
        if(_ngbs[i] < (int)cells.size()) {
            pos[0] = cells[_ngbs[i]]->_pos[0];
            pos[1] = cells[_ngbs[i]]->_pos[1];
            pos[2] = cells[_ngbs[i]]->_pos[2];
            if(_ngbwalls[i] < 0) {
                AdaptiveMeshUtils::get_periodic_position(pos, _ngbwalls[i],
                                                         box);
            }
        } else {
            AdaptiveMeshUtils::get_wall_position(
                    ghosts[_ngbs[i] - cells.size()]->_pos, pos,
                    ghosts[_ngbs[i] - cells.size()]->_ngbs[0], box);
        }
        stream << "_ngbs[" << i << "] = " << _ngbs[i] << " (" << _ngbwalls[i]
               << "): " << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\n";
    }
    stream << "\n";
    unsigned int numfacengb = 0;
    for(unsigned int i = 0; i < _facengbs.size(); i++) {
        stream << "_facengbs[" << i << "]: {";
        for(unsigned int j = 0; j < _facengbs[i].size() - 1; j++) {
            stream << _facengbs[i][j] << " (" << _facengbwalls[i][j] << "), ";
            numfacengb++;
        }
        stream << _facengbs[i].back() << " (" << _facengbwalls[i].back()
               << ")}\n";
        numfacengb++;
    }
    stream << "\n";
}

/**
 * @brief Calculate faces for the cell
 *
 * The vertices of the faces are the midpoints of the circumspheres traced out
 * by the cell generator, one of its neighbours and two consecutive facengbs of
 * this neighbour.
 *
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell3d::complete(vector<AdaptiveVorCell3d*>& cells,
                                 vector<AdaptiveVorCell3d*>& ghosts,
                                 Cuboid& box) {
    // reset the list of faces and make it the same size as the ngb list
    cleanup_faces();
    _faces.resize(_ngbs.size(), NULL);
    // calculate a face for every ngb
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        double ngbpos[3];
        // get the position of the ngb, taking into account the boundaries
        if(_ngbs[i] < (int)cells.size()) {
            AdaptiveVorCell3d* ngb = cells[_ngbs[i]];
            ngbpos[0] = ngb->_pos[0];
            ngbpos[1] = ngb->_pos[1];
            ngbpos[2] = ngb->_pos[2];
            if(_ngbwalls[i] < 0) {
                AdaptiveMeshUtils::get_periodic_position(ngbpos, _ngbwalls[i],
                                                         box);
            }
        } else {
            AdaptiveVorCell3d* ngb = ghosts[_ngbs[i] - cells.size()];
            AdaptiveMeshUtils::get_wall_position(ngb->_pos, ngbpos,
                                                 ngb->_ngbs[0], box);
        }
        // create the face
        _faces[i] = new AdaptiveVorFace3d();
        // walk the facengbs for this ngb and create vertices
        for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
            // get the positions of the current and next facengb,
            // taking into account the boundaries
            double prevpos[3], nextpos[3];
            if(_facengbs[i][j] < (int)cells.size()) {
                AdaptiveVorCell3d* prev = cells[_facengbs[i][j]];
                prevpos[0] = prev->_pos[0];
                prevpos[1] = prev->_pos[1];
                prevpos[2] = prev->_pos[2];
                if(_facengbwalls[i][j] < 0) {
                    AdaptiveMeshUtils::get_periodic_position(
                            prevpos, _facengbwalls[i][j], box);
                }
            } else {
                AdaptiveVorCell3d* prev =
                        ghosts[_facengbs[i][j] - cells.size()];
                AdaptiveMeshUtils::get_wall_position(prev->_pos, prevpos,
                                                     prev->_ngbs[0], box);
            }
            if(_facengbs[i][(j + 1) % _facengbs[i].size()] <
               (int)cells.size()) {
                AdaptiveVorCell3d* next =
                        cells[_facengbs[i][(j + 1) % _facengbs[i].size()]];
                nextpos[0] = next->_pos[0];
                nextpos[1] = next->_pos[1];
                nextpos[2] = next->_pos[2];
                if(_facengbwalls[i][(j + 1) % _facengbwalls[i].size()] < 0) {
                    AdaptiveMeshUtils::get_periodic_position(
                            nextpos,
                            _facengbwalls[i][(j + 1) % _facengbwalls[i].size()],
                            box);
                }
            } else {
                AdaptiveVorCell3d* next =
                        ghosts[_facengbs[i][(j + 1) % _facengbs[i].size()] -
                               cells.size()];
                AdaptiveMeshUtils::get_wall_position(next->_pos, nextpos,
                                                     next->_ngbs[0], box);
            }
            // calculate the vertex position
            double vert[3];
            if(get_vertexpoint(_pos, ngbpos, prevpos, nextpos, vert)) {
                _faces[i]->add_vertexpoint(vert);
            } else {
                // throw an error
                throw AdaptiveMeshException();
                cerr << "Problem in vertex calculation. We better stop."
                     << endl;
                my_exit();
            }
        }
    }
}

/**
 * @brief Print the cell in plottable format to the given stream
 *
 * @param stream std::ostream to write to
 * @param id Integer ID of the cell in the AdaptiveVorTess3d
 */
void AdaptiveVorCell3d::print(ostream& stream, int id) {
    stream << _pos[0] << "\t" << _pos[1] << "\t" << _pos[2];
    if(id >= 0) {
        stream << "\t" << id << endl;
    }
    stream << "\n\n\n";
    for(unsigned int i = 0; i < _faces.size(); i++) {
        if(_faces[i])
            _faces[i]->print(stream);
    }
    stream << endl;
}

/**
 * @brief Print the cell in plottable format to the given stream, applying the
 * given boundary flag to all faces
 *
 * @param stream std::ostream to write to
 * @param id Integer ID of the cell in the AdaptiveVorTess3d
 * @param wall Boundary flag
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell3d::print_copy(ostream& stream, int id, int wall,
                                   Cuboid& box) {
    double pos[3];
    pos[0] = _pos[0];
    pos[1] = _pos[1];
    pos[2] = _pos[2];
    if(wall < 0) {
        AdaptiveMeshUtils::get_periodic_position(pos, wall, box);
    }
    stream << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\t" << id
           << "\n\n\n";
    for(unsigned int i = 0; i < _faces.size(); i++) {
        if(_faces[i])
            _faces[i]->print_copy(stream, wall, box);
    }
    stream << endl;
}

/**
 * @brief Move the generator to the given position
 *
 * @param pos New position for the cell generator
 */
void AdaptiveVorCell3d::move(double* pos) {
    _pos[0] = pos[0];
    _pos[1] = pos[1];
    _pos[2] = pos[2];
}

/**
 * @brief Check if the cell is a valid Voronoi cell
 *
 * @warning This method scales very badly and is only meant for debugging
 * in small simulations
 *
 * @param id Unsigned integer ID of the cell in the AdaptiveVorTess3d cell list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param estream std::ostream to write error messages to
 * @param box Cuboid specifying the dimensions of the simulation box
 * @return True if errors were detected, false otherwise
 */
bool AdaptiveVorCell3d::check(unsigned int id,
                              vector<AdaptiveVorCell3d*>& cells,
                              ostream& estream, Cuboid& box) {
    bool wrong = false;
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        int ngb = _ngbs[i];
        if(ngb < (int)cells.size()) {
            for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
                int next = _facengbs[i][j];
                int prev = _facengbs[i][(j + 1) % _facengbs[i].size()];
                if(next < (int)cells.size() && prev < (int)cells.size()) {
                    AdaptiveVorCell3d* cell[3] = {cells[ngb], cells[next],
                                                  cells[prev]};
                    double cell0pos[3];
                    double cell1pos[3];
                    double cell2pos[3];
                    for(unsigned int k = 0; k < 3; k++) {
                        cell0pos[k] = cell[0]->_pos[k];
                        cell1pos[k] = cell[1]->_pos[k];
                        cell2pos[k] = cell[2]->_pos[k];
                    }
                    if(_ngbwalls[i] < 0) {
                        AdaptiveMeshUtils::get_periodic_position(
                                cell0pos, _ngbwalls[i], box);
                    }
                    if(_facengbwalls[i][j] < 0) {
                        AdaptiveMeshUtils::get_periodic_position(
                                cell1pos, _facengbwalls[i][j], box);
                    }
                    if(_facengbwalls[i][(j + 1) % _facengbwalls[i].size()] <
                       0) {
                        AdaptiveMeshUtils::get_periodic_position(
                                cell2pos,
                                _facengbwalls[i][(j + 1) %
                                                 _facengbwalls[i].size()],
                                box);
                    }
                    //                    if(predicates::orient3d_old(_pos,
                    //                    cell0pos, cell1pos,
                    //                    cell2pos) < 0.){
                    //                        cerr << "Error!" << endl;
                    //                        print(estream, id);
                    //                        estream.flush();
                    //                        exit(1);
                    //                    } else {
                    //                        continue;
                    //                    }
                    for(unsigned int k = 0; k < cells.size(); k++) {
                        if((int)k != ngb && (int)k != next && (int)k != prev) {
                            AdaptiveVorCell3d* cell4 = cells[k];
                            double cell4pos[3];
                            cell4pos[0] = cell4->_pos[0];
                            cell4pos[1] = cell4->_pos[1];
                            cell4pos[2] = cell4->_pos[2];
                            //                                cout <<
                            //                                predicates::orient3d(_pos,
                            // cell0pos, cell1pos, cell2pos) << endl;;
                            if(inside(_pos, cell0pos, cell1pos, cell2pos,
                                      cell4pos)) {
                                double inspherevalue = predicates::insphere_old(
                                        _pos, cell0pos, cell1pos, cell2pos,
                                        cell4pos);
                                cout << inspherevalue << endl;
                                wrong |= true;
                                double vert[3];
                                get_vertexpoint(_pos, cell[0]->_pos,
                                                cell[1]->_pos, cell[2]->_pos,
                                                vert);
                                estream << vert[0] << "\t" << vert[1] << "\t"
                                        << vert[2] << "\n";
                                stringstream filename;
                                filename << "wrongcells.dat";
                                ofstream cfile(filename.str().c_str());
                                this->print(cfile, _id);
                                cell[0]->print_copy(cfile, cell[0]->_id,
                                                    _ngbwalls[i], box);
                                cell[1]->print_copy(cfile, cell[1]->_id,
                                                    _facengbwalls[i][j], box);
                                cell[2]->print_copy(
                                        cfile, cell[2]->_id,
                                        _facengbwalls[i]
                                                     [(j + 1) %
                                                      _facengbwalls[i].size()],
                                        box);
                                cell4->print(cfile, cell4->_id);
                                cfile.close();
                                cerr << _pos[0] << "\t" << _pos[1] << "\t"
                                     << _pos[2] << endl;
                                cerr << cell0pos[0] << "\t" << cell0pos[1]
                                     << "\t" << cell0pos[2] << endl;
                                cerr << cell1pos[0] << "\t" << cell1pos[1]
                                     << "\t" << cell1pos[2] << endl;
                                cerr << cell2pos[0] << "\t" << cell2pos[1]
                                     << "\t" << cell2pos[2] << endl;
                                cerr << cell4pos[0] << "\t" << cell4pos[1]
                                     << "\t" << cell4pos[2] << endl;
                                cerr << id << "\t" << ngb << "\t" << next
                                     << "\t" << prev << endl;
                                exit(1);
                            }
                        }
                    }
                    vector<int> ngbs = cells[_ngbs[i]]->_ngbs;
                    for(unsigned int k = 0; k < ngbs.size(); k++) {
                        int ngbngb = ngbs[k];
                        if(ngbngb != ngb && ngbngb != next && ngbngb != prev &&
                           ngbngb != (int)id) {
                            AdaptiveVorCell3d* cell4 = cells[ngbngb];
                            double cell4pos[3];
                            cell4pos[0] = cell4->_pos[0];
                            cell4pos[1] = cell4->_pos[1];
                            cell4pos[2] = cell4->_pos[2];
                            //                                cout <<
                            // predicates::orient3d(_pos,
                            //                            cell0pos, cell1pos,
                            // cell2pos) << endl;;
                            if(inside(_pos, cell0pos, cell1pos, cell2pos,
                                      cell4pos)) {
                                //                                double
                                // inspherevalue =
                                //                                predicates::insphere(_pos,
                                // cell0pos,
                                //                                cell1pos,
                                // cell2pos, cell4pos);
                                //                                    cout <<
                                // inspherevalue << endl;
                                wrong |= true;
                                double vert[3];
                                get_vertexpoint(_pos, cell[0]->_pos,
                                                cell[1]->_pos, cell[2]->_pos,
                                                vert);
                                estream << vert[0] << "\t" << vert[1] << "\t"
                                        << vert[2] << "\n";
                                stringstream filename;
                                filename << "wrongcells.dat";
                                ofstream cfile(filename.str().c_str());
                                this->print(cfile, _id);
                                cell[0]->print_copy(cfile, cell[0]->_id,
                                                    _ngbwalls[i], box);
                                cell[1]->print_copy(cfile, cell[1]->_id,
                                                    _facengbwalls[i][j], box);
                                cell[2]->print_copy(
                                        cfile, cell[2]->_id,
                                        _facengbwalls[i]
                                                     [(j + 1) %
                                                      _facengbwalls[i].size()],
                                        box);
                                cell4->print(cfile, cell4->_id);
                                cfile.close();
                                cerr << _pos[0] << "\t" << _pos[1] << "\t"
                                     << _pos[2] << endl;
                                cerr << cell0pos[0] << "\t" << cell0pos[1]
                                     << "\t" << cell0pos[2] << endl;
                                cerr << cell1pos[0] << "\t" << cell1pos[1]
                                     << "\t" << cell1pos[2] << endl;
                                cerr << cell2pos[0] << "\t" << cell2pos[1]
                                     << "\t" << cell2pos[2] << endl;
                                cerr << cell4pos[0] << "\t" << cell4pos[1]
                                     << "\t" << cell4pos[2] << endl;
                                cerr << id << "\t" << ngb << "\t" << next
                                     << "\t" << prev << endl;
                                my_exit();
                            }
                        }
                    }
                }
            }
        }
    }
    return wrong;
}

/**
 * @brief Get the position of the cell generator
 *
 * @param pos Array to store the result in
 */
void AdaptiveVorCell3d::get_position(double* pos) {
    pos[0] = _pos[0];
    pos[1] = _pos[1];
    pos[2] = _pos[2];
}

/**
 * @brief Check if the given cell index is in the neighbour list
 *
 * @param id Integer index of a cell in the AdaptiveVorTess3d cell list
 * @return True if the cell is a neighbour of this cell, false otherwise
 */
bool AdaptiveVorCell3d::is_ngb(int id) {
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        if(_ngbs[i] == id) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Remove the neighbour with the given index from the neighbour list
 *
 * If the given index is not found in the neighbour list, the code will abort.
 *
 * If it is found, the corresponding fields in the other internal lists are
 * also removed.
 *
 * @param id Unsigned integer index of a cell in the AdaptiveVorTess3d cell list
 * @return Unsigned integer index of the removed neighbour in the internal list
 */
unsigned int AdaptiveVorCell3d::remove_ngb(unsigned int id) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != (int)id) {
        i++;
    }
    if(i == _ngbs.size()) {
        throw AdaptiveMeshException();
        cerr << "asked to remove ngb that is not ngb!" << endl;
        my_exit();
    }
    _ngbs.erase(_ngbs.begin() + i);
    _ngbwalls.erase(_ngbwalls.begin() + i);
    _facengbs.erase(_facengbs.begin() + i);
    _facengbwalls.erase(_facengbwalls.begin() + i);
    return i;
}

/**
 * @brief Remove the facengb with the given index from the list corresponding to
 * the neighbour with the given index
 *
 * If the given neighbour is not found in the internal neighbour list, or the
 * given facengb is not found in the facengb list corresponding to the
 * neighbour, the code will abort.
 *
 * If both are found, we also remove the corresponding boundary flag.
 *
 * @param ngb Integer index of a cell in the AdaptiveVorTess3d cell list
 * @param facengb Integer index of a cell in the AdaptiveVorTess3d cell list
 */
void AdaptiveVorCell3d::remove_facengb(int ngb, int facengb) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != ngb) {
        i++;
    }
    if(i == _ngbs.size()) {
        throw AdaptiveMeshException();
        cerr << "A terrible internal accident happened. We're sorry..." << endl;
        cerr << "Trying to remove ngb " << ngb << " from " << _id << endl;
        cerr << "... *CRASH* ..." << endl;
        my_exit();
    }
    unsigned int j = 0;
    while(j < _facengbs[i].size() && _facengbs[i][j] != facengb) {
        j++;
    }
    if(j == _facengbs[i].size()) {
        throw AdaptiveMeshException();
        cerr << "A terrible internal accident happened. We're sorry..." << endl;
        cerr << "Trying to remove facengb " << facengb << " from " << _id
             << endl;
        cerr << "... *CRASH* ..." << endl;
        my_exit();
    }
    _facengbs[i].erase(_facengbs[i].begin() + j);
    _facengbwalls[i].erase(_facengbwalls[i].begin() + j);
}

/**
 * @brief Remove the facengb with the given index from the list corresponding to
 * the neighbour with the given index
 *
 * If the given neighbour is not found in the internal neighbour list, or the
 * given facengb is not found in the facengb list corresponding to the
 * neighbour, the code will abort.
 *
 * If both are found, we also remove the corresponding boundary flag.
 *
 * @param edges Array to store information about the edges that are affected by
 * this removal
 * @param ngb Integer index of a cell in the AdaptiveVorTess3d cell list
 * @param facengb Integer index of a cell in the AdaptiveVorTess3d cell list
 */
void AdaptiveVorCell3d::remove_facengb(unsigned int* edges, int ngb,
                                       int facengb) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != ngb) {
        i++;
    }
    if(i == _ngbs.size()) {
        throw AdaptiveMeshException();
        cerr << "A terrible internal accident happened. We're sorry..." << endl;
        cerr << "Trying to remove ngb " << ngb << " from " << _id << endl;
        cerr << "... *CRASH* ..." << endl;
        my_exit();
    }
    unsigned int j = 0;
    while(j < _facengbs[i].size() && _facengbs[i][j] != facengb) {
        j++;
    }
    if(j == _facengbs[i].size()) {
        throw AdaptiveMeshException();
        cerr << "A terrible internal accident happened. We're sorry..." << endl;
        cerr << "Trying to remove facengb " << facengb << " from " << _id
             << endl;
        cerr << "... *CRASH* ..." << endl;
        my_exit();
    }
    _facengbs[i].erase(_facengbs[i].begin() + j);
    _facengbwalls[i].erase(_facengbwalls[i].begin() + j);
    edges[0] = i;
    //    edges[1] = (j+_facengbs[i].size()-1)%_facengbs[i].size();
    //    // it could happen that you just removed the last item in the
    // vector...
    //    edges[2] = j%_facengbs[i].size();
    edges[1] = (j + _facengbs[i].size() - 2) % _facengbs[i].size();
    // it could happen that you just removed the last item in the vector...
    edges[2] = (j + _facengbs[i].size() - 1) % _facengbs[i].size();
}

/**
 * @brief Remove the facengb with the given index from the list corresponding to
 * the neighbour with the given index
 *
 * Unlike AdaptiveVorCell3d::remove_facengb(int ngb, int facengb), this method
 * never aborts, hence its name.
 *
 * If both are found, we also remove the corresponding boundary flag.
 *
 * @param ngb Integer index of a cell in the AdaptiveVorTess3d cell list
 * @param facengb Integer index of a cell in the AdaptiveVorTess3d cell list
 */
void AdaptiveVorCell3d::safely_remove_facengb(int ngb, int facengb) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != ngb) {
        i++;
    }
    if(i == _ngbs.size()) {
        // safely means: no crash here
        return;
    }
    unsigned int j = 0;
    while(j < _facengbs[i].size() && _facengbs[i][j] != facengb) {
        j++;
    }
    if(j == _facengbs[i].size()) {
        // nor here
        return;
    }
    _facengbs[i].erase(_facengbs[i].begin() + j);
    _facengbwalls[i].erase(_facengbwalls[i].begin() + j);
}

/**
 * @brief Add the given facengb in the list for the given neighbour, before the
 * given facengb and with the corresponding boundary flag
 *
 * If the given neighbour is not found, or the given facengb is not found in the
 * list corresponding to this neighbour, the code will abort.
 *
 * @param edges Array to store information about the edges that are affected by
 * this insertion
 * @param ngb Integer index of a cell in the AdaptiveVorTess3d cell list
 * @param facengb Integer index of a cell in the AdaptiveVorTess3d cell list
 * @param before Integer index of a cell in the AdaptiveVorTess3d cell list
 * @param wall Boundary flag for the new facengb
 */
void AdaptiveVorCell3d::add_facengb(unsigned int* edges, int ngb, int facengb,
                                    int before, int wall) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != ngb) {
        i++;
    }
    if(i == _ngbs.size()) {
        throw AdaptiveMeshException();
        cerr << "Not a ngb!" << endl;
        my_exit();
    }
    unsigned int j = 0;
    while(j < _facengbs[i].size() && _facengbs[i][j] != before) {
        j++;
    }
    if(j == _facengbs[i].size()) {
        throw AdaptiveMeshException();
        cerr << "Not a facengb!" << endl;
        my_exit();
    }
    _facengbs[i].insert(_facengbs[i].begin() + j, facengb);
    _facengbwalls[i].insert(_facengbwalls[i].begin() + j, wall);
    edges[0] = i;
    //    edges[1] = (j+_facengbs[i].size()-1)%_facengbs[i].size();
    //    edges[2] = j;
    //    edges[3] = (j+1)%_facengbs[i].size();
    edges[1] = (j + _facengbs[i].size() - 2) % _facengbs[i].size();
    edges[2] = (j + _facengbs[i].size() - 1) % _facengbs[i].size();
    edges[3] = j;
}

/**
 * @brief Add the given facengb in the list for the given neighbour, before the
 * given facengb and with the corresponding boundary flag
 *
 * If the given neighbour is not found, or the given facengb is not found in the
 * list corresponding to this neighbour, the code will abort.
 *
 * @param ngb Integer index of a cell in the AdaptiveVorTess3d cell list
 * @param facengb Integer index of a cell in the AdaptiveVorTess3d cell list
 * @param before Integer index of a cell in the AdaptiveVorTess3d cell list
 * @param wall Boundary flag for the new facengb
 */
void AdaptiveVorCell3d::add_facengb(int ngb, int facengb, int before,
                                    int wall) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != ngb) {
        i++;
    }
    if(i == _ngbs.size()) {
        throw AdaptiveMeshException();
        cerr << "Not a ngb!" << endl;
        my_exit();
    }
    unsigned int j = 0;
    while(j < _facengbs[i].size() && _facengbs[i][j] != before) {
        j++;
    }
    if(j == _facengbs[i].size()) {
        throw AdaptiveMeshException();
        cerr << "Not a facengb!" << endl;
        my_exit();
    }
    _facengbs[i].insert(_facengbs[i].begin() + j, facengb);
    _facengbwalls[i].insert(_facengbwalls[i].begin() + j, wall);
}

/**
 * @brief Create a new face between C and D, with facengbs A, B, and E
 *
 * The method updates the appropriate faces in A, B, and E and removes the
 * flipped edges by removing the appropriate facengbs from A, B, and E.
 *
 * If something goes wrong during the face creation, the code will abort.
 *
 * @param edges_to_check Array to store information about edges that have to be
 * checked after this face creation.
 * @param A Integer index of cell A in the AdaptiveVorTess3d cell list
 * @param B Integer index of cell B in the AdaptiveVorTess3d cell list
 * @param C Integer index of cell C in the AdaptiveVorTess3d cell list
 * @param D Integer index of cell D in the AdaptiveVorTess3d cell list
 * @param E Integer index of cell E in the AdaptiveVorTess3d cell list
 * @param wallA Boundary flag of cell A
 * @param wallB Boundary flag of cell B
 * @param wallC Boundary flag of cell C
 * @param wallD Boundary flag of cell D
 * @param wallE Boundary flag of cell E
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 * @return True if the creation succeeded, false otherwise
 */
bool AdaptiveVorCell3d::create_new_face(unsigned int* edges_to_check, int A,
                                        int B, int C, int D, int E, int wallA,
                                        int wallB, int wallC, int wallD,
                                        int wallE,
                                        vector<AdaptiveVorCell3d*>& cells,
                                        vector<AdaptiveVorCell3d*>& ghosts,
                                        bool periodic) {
#ifdef VERBOSE
    cerr << "A = " << A << ", B = " << B << ", C = " << C << ", D = " << D
         << ", E = " << E << endl;
#endif
    unsigned int currentedge = 0;
    unsigned int edges[4];
    if(C < (int)cells.size()) {
        // it can happen that an insertion is incorrectly flagged due to
        // another insertion that still has to happen and that masks a
        // removal. Long story short: in this case we can safely skip
        // this insertion and continue.
        if(cells[C]->is_ngb(D)) {
            throw AdaptiveMeshException();
            // for closer inspection
            cerr << "Aborted" << endl;
            ofstream afile("abortfile.dat");
            Vec anchor;
            Vec sides(1., 1., 1.);
            Cuboid box(anchor, sides);
            cells[A]->complete(cells, ghosts, box);
            cells[B]->complete(cells, ghosts, box);
            cells[C]->complete(cells, ghosts, box);
            cells[D]->complete(cells, ghosts, box);
            cells[E]->complete(cells, ghosts, box);
            cells[A]->print_copy(afile, 0, wallA, box);
            cells[B]->print_copy(afile, 1, wallB, box);
            cells[C]->print_copy(afile, 2, wallC, box);
            cells[D]->print_copy(afile, 3, wallD, box);
            cells[E]->print_copy(afile, 4, wallE, box);
            afile.close();
            my_exit();
            return false;
        }
        int wallpos = AdaptiveMeshUtils::get_wallpos(wallC, wallD);
        unsigned int ngbi = cells[C]->add_ngb(D, wallpos);
        // order!
        cells[C]->add_facengb(A, AdaptiveMeshUtils::get_wallpos(wallC, wallA));
        cells[C]->add_facengb(E, AdaptiveMeshUtils::get_wallpos(wallC, wallE));
        cells[C]->add_facengb(B, AdaptiveMeshUtils::get_wallpos(wallC, wallB));
        edges_to_check[currentedge++] = 2;  // C
        edges_to_check[currentedge++] = ngbi;
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = 2;  // C
        edges_to_check[currentedge++] = ngbi;
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = 2;  // C
        edges_to_check[currentedge++] = ngbi;
        edges_to_check[currentedge++] = 2;

        // add facengb D before B in the list for ngb A
        cells[C]->add_facengb(edges, A, D, B, wallpos);
        edges_to_check[currentedge++] = 2;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 2;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        edges_to_check[currentedge++] = 2;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[3];
        cells[C]->add_facengb(edges, E, D, A, wallpos);
        edges_to_check[currentedge++] = 2;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 2;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        edges_to_check[currentedge++] = 2;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[3];
        cells[C]->add_facengb(edges, B, D, E, wallpos);
        edges_to_check[currentedge++] = 2;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 2;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        edges_to_check[currentedge++] = 2;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[3];
    }

    if(D < (int)cells.size()) {
        int wallpos = AdaptiveMeshUtils::get_wallpos(wallD, wallC);
        //        unsigned int ngbi = cells[D]->add_ngb(C, wallpos);
        cells[D]->add_ngb(C, wallpos);
        // order!
        cells[D]->add_facengb(A, AdaptiveMeshUtils::get_wallpos(wallD, wallA));
        cells[D]->add_facengb(B, AdaptiveMeshUtils::get_wallpos(wallD, wallB));
        cells[D]->add_facengb(E, AdaptiveMeshUtils::get_wallpos(wallD, wallE));
        // face already added
        //        edges_to_check[currentedge++] = 3; // D
        //        edges_to_check[currentedge++] = ngbi;
        //        edges_to_check[currentedge++] = 0;
        //        edges_to_check[currentedge++] = 3; // D
        //        edges_to_check[currentedge++] = ngbi;
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = 3; // D
        //        edges_to_check[currentedge++] = ngbi;
        //        edges_to_check[currentedge++] = 2;

        cells[D]->add_facengb(edges, A, C, E, wallpos);
        edges_to_check[currentedge++] = 3;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 3;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        edges_to_check[currentedge++] = 3;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[3];
        cells[D]->add_facengb(edges, B, C, A, wallpos);
        edges_to_check[currentedge++] = 3;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 3;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        edges_to_check[currentedge++] = 3;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[3];
        cells[D]->add_facengb(edges, E, C, B, wallpos);
        edges_to_check[currentedge++] = 3;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 3;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        edges_to_check[currentedge++] = 3;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[3];
    }

    if(A < (int)cells.size()) {
        cells[A]->add_facengb(edges, C, D, E,
                              AdaptiveMeshUtils::get_wallpos(wallA, wallD));
        //        edges_to_check[currentedge++] = 0;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 0;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        edges_to_check[currentedge++] = 0;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[3];
        cells[A]->add_facengb(edges, D, C, B,
                              AdaptiveMeshUtils::get_wallpos(wallA, wallC));
        //        edges_to_check[currentedge++] = 0;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 0;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        edges_to_check[currentedge++] = 0;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[3];
        // remove_facengb(<ngb that causes face>, <facengb to remove>)
        cells[A]->remove_facengb(edges, E, B);
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        cells[A]->remove_facengb(edges, B, E);
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
    }

    if(B < (int)cells.size()) {
        cells[B]->add_facengb(edges, C, D, A,
                              AdaptiveMeshUtils::get_wallpos(wallB, wallD));
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[3];
        cells[B]->add_facengb(edges, D, C, E,
                              AdaptiveMeshUtils::get_wallpos(wallB, wallC));
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[3];
        cells[B]->remove_facengb(edges, A, E);
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        cells[B]->remove_facengb(edges, E, A);
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
    }

    if(E < (int)cells.size()) {
        cells[E]->add_facengb(edges, C, D, B,
                              AdaptiveMeshUtils::get_wallpos(wallE, wallD));
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[3];
        cells[E]->add_facengb(edges, D, C, A,
                              AdaptiveMeshUtils::get_wallpos(wallE, wallC));
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[3];
        cells[E]->remove_facengb(edges, A, B);
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        cells[E]->remove_facengb(edges, B, A);
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
    }
    return true;
}

/**
 * @brief Remove the face between C and D, and the appropriate facengbs in A, B,
 * and E
 *
 * The method also creates new facengbs in A, B, and E, and removes the
 * appropriate facengbs in C and D.
 *
 * If something goes wrong during the face removal, the code will abort.
 *
 * @param edges_to_check Array to store information about the edges that should
 * be checked after the face removal
 * @param A Integer index of cell A in the AdaptiveVorTess3d cell list
 * @param B Integer index of cell B in the AdaptiveVorTess3d cell list
 * @param C Integer index of cell C in the AdaptiveVorTess3d cell list
 * @param D Integer index of cell D in the AdaptiveVorTess3d cell list
 * @param E Integer index of cell E in the AdaptiveVorTess3d cell list
 * @param wallA Boundary flag for cell A
 * @param wallB Boundary flag for cell B
 * @param wallC Boundary flag for cell C
 * @param wallD Boundary flag for cell D
 * @param wallE Boundary flag for cell E
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 * @return True if the face removal succeeded, false otherwise
 */
bool AdaptiveVorCell3d::remove_face(unsigned int* edges_to_check, int A, int B,
                                    int C, int D, int E, int wallA, int wallB,
                                    int wallC, int wallD, int wallE,
                                    std::vector<AdaptiveVorCell3d*>& cells,
                                    std::vector<AdaptiveVorCell3d*>& ghosts,
                                    bool periodic) {
#ifdef VERBOSE
    cerr << "A = " << A << ", B = " << B << ", C = " << C << ", D = " << D
         << ", E = " << E << endl;
#endif
    // we first remove the facengbs and then the ngbs, since
    // remove_facengb needs the ngb to be present in the _ngbs list
    // wait a minute... it does not. Anyhow, does not really matter
    unsigned int currentedge = 0;
    unsigned int edges[4];

    if(C < (int)cells.size()) {
        //        unsigned int backup[6];
        cells[C]->remove_facengb(edges, A, D);
        //        edges_to_check[currentedge++] = 2;
        //        backup[0] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 2;
        //        backup[1] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        cells[C]->remove_facengb(edges, B, D);
        //        edges_to_check[currentedge++] = 2;
        //        backup[2] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 2;
        //        backup[3] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        cells[C]->remove_facengb(edges, E, D);
        //        edges_to_check[currentedge++] = 2;
        //        backup[4] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 2;
        //        backup[5] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        unsigned int ngbi = cells[C]->remove_ngb(D);
        cells[C]->remove_ngb(D);
        //        for(unsigned int i = 0; i < 6; i++){
        //            // this value cannot be equal, since we just removed this
        // face
        //            if(edges_to_check[backup[i]] > ngbi){
        //                edges_to_check[backup[i]]--;
        //            }
        //        }
    }

    if(D < (int)cells.size()) {
        //        unsigned int backup[6];
        cells[D]->remove_facengb(edges, A, C);
        //        edges_to_check[currentedge++] = 3;
        //        backup[0] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 3;
        //        backup[1] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        cells[D]->remove_facengb(edges, B, C);
        //        edges_to_check[currentedge++] = 3;
        //        backup[2] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 3;
        //        backup[3] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        cells[D]->remove_facengb(edges, E, C);
        //        edges_to_check[currentedge++] = 3;
        //        backup[4] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 3;
        //        backup[5] = currentedge;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        unsigned int ngbi = cells[D]->remove_ngb(C);
        cells[D]->remove_ngb(C);
        //        for(unsigned int i = 0; i < 6; i++){
        //            // this value cannot be equal, since we just removed this
        // face
        //            if(edges_to_check[backup[i]] > ngbi){
        //                edges_to_check[backup[i]]--;
        //            }
        //        }
    }

    if(A < (int)cells.size()) {
        cells[A]->remove_facengb(edges, C, D);
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        cells[A]->remove_facengb(edges, D, C);
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        cells[A]->add_facengb(edges, E, B, D,
                              AdaptiveMeshUtils::get_wallpos(wallA, wallB));
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[3];
        cells[A]->add_facengb(edges, B, E, C,
                              AdaptiveMeshUtils::get_wallpos(wallA, wallE));
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        edges_to_check[currentedge++] = 0;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[3];
    }

    if(B < (int)cells.size()) {
        cells[B]->remove_facengb(edges, C, D);
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        cells[B]->remove_facengb(edges, D, C);
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        cells[B]->add_facengb(edges, A, E, D,
                              AdaptiveMeshUtils::get_wallpos(wallB, wallE));
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        edges_to_check[currentedge++] = 1;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[3];
        cells[B]->add_facengb(edges, E, A, C,
                              AdaptiveMeshUtils::get_wallpos(wallB, wallA));
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        edges_to_check[currentedge++] = 1;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[3];
    }

    if(E < (int)cells.size()) {
        cells[E]->remove_facengb(edges, C, D);
        edges_to_check[currentedge++] = 4;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 4;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        cells[E]->remove_facengb(edges, D, C);
        edges_to_check[currentedge++] = 4;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[1];
        edges_to_check[currentedge++] = 4;
        edges_to_check[currentedge++] = edges[0];
        edges_to_check[currentedge++] = edges[2];
        cells[E]->add_facengb(edges, A, B, C,
                              AdaptiveMeshUtils::get_wallpos(wallE, wallB));
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[3];
        cells[E]->add_facengb(edges, B, A, D,
                              AdaptiveMeshUtils::get_wallpos(wallE, wallA));
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[1];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[2];
        //        edges_to_check[currentedge++] = 4;
        //        edges_to_check[currentedge++] = edges[0];
        //        edges_to_check[currentedge++] = edges[3];
    }
    return true;
}

/**
 * @brief Check if one of the three faces sharing the edge between cells A, B,
 * and C is a triangle and completely flips
 *
 * The function returns a flag returns one of the following values:
 *  - 0: no triangular face
 *  - 1: triangular face between A and B
 *  - 2: triangular face between A and C
 *  - 3: triangular face between B and C
 *
 * @param A Integer index of cell A in the AdaptiveVorTess3d cell list
 * @param B Integer index of cell B in the AdaptiveVorTess3d cell list
 * @param C Integer index of cell C in the AdaptiveVorTess3d cell list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @return A value in the range 0-3 indicating the result of the query
 */
unsigned int AdaptiveVorCell3d::check_triangular_face(
        int A, int B, int C, vector<AdaptiveVorCell3d*>& cells) {
    //        cerr << "check_triangular_face:" << endl;
    //        cerr << "A: " << A << ", B: " << B << ", C: " << C << endl;
    unsigned int i = 0;
    while(i < cells[A]->_ngbs.size() && cells[A]->_ngbs[i] != B &&
          cells[A]->_ngbs[i] != C) {
        i++;
    }
    if(i == cells[A]->_ngbs.size()) {
        throw AdaptiveMeshException();
        cerr << "Error in ngbsearch for A" << endl;
        my_exit();
    }
    //        cerr << "Asize1: " << cells[A]->_facengbs[i].size() << endl;
    if(cells[A]->_facengbs[i].size() < 4) {
        if(cells[A]->_ngbs[i] == B) {
            return 1;
        } else {
            return 2;
        }
    }
    i++;
    while(i < cells[A]->_ngbs.size() && cells[A]->_ngbs[i] != B &&
          cells[A]->_ngbs[i] != C) {
        i++;
    }
    if(i == cells[A]->_ngbs.size()) {
        throw AdaptiveMeshException();
        cerr << "Error in ngbsearch for A" << endl;
        my_exit();
    }
    //        cerr << "Asize2: " << cells[A]->_facengbs[i].size() << endl;
    if(cells[A]->_facengbs[i].size() < 4) {
        if(cells[A]->_ngbs[i] == B) {
            return 1;
        } else {
            return 2;
        }
    }
    // we only have to check one more face: the one between B and C
    i = 0;
    while(i < cells[B]->_ngbs.size() && cells[B]->_ngbs[i] != C) {
        i++;
    }
    if(i == cells[B]->_ngbs.size()) {
        throw AdaptiveMeshException();
        cerr << "Error in ngbsearch for B" << endl;
        my_exit();
    }
    //        cerr << "Bsize: " << cells[B]->_facengbs[i].size() << endl;
    if(cells[B]->_facengbs[i].size() < 4) {
        return 3;
    }
    return 0;
}

// void AdaptiveVorCell3d::trim_first_position(double *pos, Cuboid &box,
// bool* trims){
//    if(pos[0] > box.get_anchor()[0]+0.5*box.get_sides()[0]){
//        pos[0] -= box.get_sides()[0];
//        trims[0] = true;
//    } else {
//        trims[0] = false;
//    }
//    if(pos[1] > box.get_anchor()[1]+0.5*box.get_sides()[1]){
//        pos[1] -= box.get_sides()[1];
//        trims[1] = true;
//    } else {
//        trims[1] = false;
//    }
//    if(pos[2] > box.get_anchor()[2]+0.5*box.get_sides()[2]){
//        pos[2] -= box.get_sides()[2];
//        trims[2] = true;
//    } else {
//        trims[2] = false;
//    }
//}

/**
 * @brief Convert the given position to its given periodic copy
 *
 * @deprecated This method is no longer used
 *
 * @param pos Position to convert
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param trims Array of bools specifying which coordinates should be converted
 */
void AdaptiveVorCell3d::trim_position(double* pos, Cuboid& box, bool* trims) {
    if(trims[0]) {
        pos[0] -= box.get_sides()[0];
    }
    if(trims[1]) {
        pos[1] -= box.get_sides()[1];
    }
    if(trims[2]) {
        pos[2] -= box.get_sides()[2];
    }
}

/**
 * @brief Check if the generator of cell E is inside the sphere through the
 * generators of A, B, C, and D.
 *
 * @param A Integer index of cell A in the AdaptiveVorTess3d cell list
 * @param B Integer index of cell B in the AdaptiveVorTess3d cell list
 * @param C Integer index of cell C in the AdaptiveVorTess3d cell list
 * @param D Integer index of cell D in the AdaptiveVorTess3d cell list
 * @param E Integer index of cell E in the AdaptiveVorTess3d cell list
 * @param wallA Boundary flag for cell A
 * @param wallB Boundary flag for cell B
 * @param wallC Boundary flag for cell C
 * @param wallD Boundary flag for cell D
 * @param wallE Boundary flag for cell E
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 * @param verbose Bool specifying if we want this method to print information to
 * the stdout
 * @param testdouble Array to optionally store geometrical test values
 * @return True if the generator of E is inside the sphere through the
 * generators of A, B, C, and D, false otherwise
 */
bool AdaptiveVorCell3d::in_sphere_test(int A, int B, int C, int D, int E,
                                       int wallA, int wallB, int wallC,
                                       int wallD, int wallE,
                                       vector<AdaptiveVorCell3d*>& cells,
                                       vector<AdaptiveVorCell3d*>& ghosts,
                                       Cuboid& box, bool periodic, bool verbose,
                                       double* testdouble) {
    double pA[3];
    double pB[3];
    double pC[3];
    double pD[3];
    double pE[3];
    if(A >= (int)cells.size()) {
        pA[0] = ghosts[A - cells.size()]->_pos[0];
        pA[1] = ghosts[A - cells.size()]->_pos[1];
        pA[2] = ghosts[A - cells.size()]->_pos[2];
        AdaptiveMeshUtils::get_wall_position(
                pA, pA, ghosts[A - cells.size()]->_ngbs[0], box);
    } else {
        pA[0] = cells[A]->_pos[0];
        pA[1] = cells[A]->_pos[1];
        pA[2] = cells[A]->_pos[2];
    }
    if(periodic && wallA < 0) {
        AdaptiveMeshUtils::get_periodic_position(pA, wallA, box);
    }
    if(B >= (int)cells.size()) {
        pB[0] = ghosts[B - cells.size()]->_pos[0];
        pB[1] = ghosts[B - cells.size()]->_pos[1];
        pB[2] = ghosts[B - cells.size()]->_pos[2];
        AdaptiveMeshUtils::get_wall_position(
                pB, pB, ghosts[B - cells.size()]->_ngbs[0], box);
    } else {
        pB[0] = cells[B]->_pos[0];
        pB[1] = cells[B]->_pos[1];
        pB[2] = cells[B]->_pos[2];
    }
    if(periodic && wallB < 0) {
        AdaptiveMeshUtils::get_periodic_position(pB, wallB, box);
    }
    if(C >= (int)cells.size()) {
        pC[0] = ghosts[C - cells.size()]->_pos[0];
        pC[1] = ghosts[C - cells.size()]->_pos[1];
        pC[2] = ghosts[C - cells.size()]->_pos[2];
        AdaptiveMeshUtils::get_wall_position(
                pC, pC, ghosts[C - cells.size()]->_ngbs[0], box);
    } else {
        pC[0] = cells[C]->_pos[0];
        pC[1] = cells[C]->_pos[1];
        pC[2] = cells[C]->_pos[2];
    }
    if(periodic && wallC < 0) {
        AdaptiveMeshUtils::get_periodic_position(pC, wallC, box);
    }
    if(D >= (int)cells.size()) {
        pD[0] = ghosts[D - cells.size()]->_pos[0];
        pD[1] = ghosts[D - cells.size()]->_pos[1];
        pD[2] = ghosts[D - cells.size()]->_pos[2];
        AdaptiveMeshUtils::get_wall_position(
                pD, pD, ghosts[D - cells.size()]->_ngbs[0], box);
    } else {
        pD[0] = cells[D]->_pos[0];
        pD[1] = cells[D]->_pos[1];
        pD[2] = cells[D]->_pos[2];
    }
    if(periodic && wallD < 0) {
        AdaptiveMeshUtils::get_periodic_position(pD, wallD, box);
    }
    if(E >= (int)cells.size()) {
        pE[0] = ghosts[E - cells.size()]->_pos[0];
        pE[1] = ghosts[E - cells.size()]->_pos[1];
        pE[2] = ghosts[E - cells.size()]->_pos[2];
        AdaptiveMeshUtils::get_wall_position(
                pE, pE, ghosts[E - cells.size()]->_ngbs[0], box);
    } else {
        pE[0] = cells[E]->_pos[0];
        pE[1] = cells[E]->_pos[1];
        pE[2] = cells[E]->_pos[2];
    }
    if(periodic && wallE < 0) {
        AdaptiveMeshUtils::get_periodic_position(pE, wallE, box);
    }

    //    bool trims[3];
    //    trims[0] = (wallA == -3) || (wallA == -6) || (wallA == -9) ||
    //    (wallA == -12) || (wallA == -14) || (wallA == -17) || (wallA == -20)
    // ||
    //    (wallA == -23) || (wallA == -26);
    //    trims[0] |= (wallB == -3) || (wallB == -6) || (wallB == -9) ||
    //    (wallB == -12) || (wallB == -14) || (wallB == -17) || (wallB == -20)
    // ||
    //    (wallB == -23) || (wallB == -26);
    //    trims[0] |= (wallC == -3) || (wallC == -6) || (wallC == -9) ||
    //    (wallC == -12) || (wallC == -14) || (wallC == -17) || (wallC == -20)
    // ||
    //    (wallC == -23) || (wallC == -26);
    //    trims[0] |= (wallD == -3) || (wallD == -6) || (wallD == -9) ||
    //    (wallD == -12) || (wallD == -14) || (wallD == -17) || (wallD == -20)
    // ||
    //    (wallD == -23) || (wallD == -26);
    //    trims[0] |= (wallE == -3) || (wallE == -6) || (wallE == -9) ||
    //    (wallE == -12) || (wallE == -14) || (wallE == -17) || (wallE == -20)
    // ||
    //    (wallE == -23) || (wallE == -26);
    //    trims[1] = (wallA == -1) || (wallA == -2) || (wallA == -3) ||
    //    (wallA == -10) || (wallA == -11) || (wallA == -12) || (wallA == -18)
    // ||
    //            (wallA == -19) || (wallA == -20);
    //    trims[1] |= (wallB == -1) || (wallB == -2) || (wallB == -3) ||
    //    (wallB == -10) || (wallB == -11) || (wallB == -12) || (wallB == -18)
    // ||
    //            (wallB == -19) || (wallB == -20);
    //    trims[1] |= (wallC == -1) || (wallC == -2) || (wallC == -3) ||
    //    (wallC == -10) || (wallC == -11) || (wallC == -12) || (wallC == -18)
    // ||
    //            (wallC == -19) || (wallC == -20);
    //    trims[1] |= (wallD == -1) || (wallD == -2) || (wallD == -3) ||
    //    (wallD == -10) || (wallD == -11) || (wallD == -12) || (wallD == -18)
    // ||
    //            (wallD == -19) || (wallD == -20);
    //    trims[1] |= (wallE == -1) || (wallE == -2) || (wallE == -3) ||
    //    (wallE == -10) || (wallE == -11) || (wallE == -12) || (wallE == -18)
    // ||
    //            (wallE == -19) || (wallE == -20);
    //    trims[2] = (wallA == -1) || (wallA == -2) || (wallA == -3) ||
    //    (wallA == -4) || (wallA == -5) || (wallA == -6) || (wallA == -7) ||
    //            (wallA == -8) || (wallA == -9);
    //    trims[2] |= (wallB == -1) || (wallB == -2) || (wallB == -3) ||
    //    (wallB == -4) || (wallB == -5) || (wallB == -6) || (wallB == -7) ||
    //            (wallB == -8) || (wallB == -9);
    //    trims[2] |= (wallC == -1) || (wallC == -2) || (wallC == -3) ||
    //    (wallC == -4) || (wallC == -5) || (wallC == -6) || (wallC == -7) ||
    //            (wallC == -8) || (wallC == -9);
    //    trims[2] |= (wallD == -1) || (wallD == -2) || (wallD == -3) ||
    //    (wallD == -4) || (wallD == -5) || (wallD == -6) || (wallD == -7) ||
    //            (wallD == -8) || (wallD == -9);
    //    trims[2] |= (wallE == -1) || (wallE == -2) || (wallE == -3) ||
    //    (wallE == -4) || (wallE == -5) || (wallE == -6) || (wallE == -7) ||
    //            (wallE == -8) || (wallE == -9);
    //    trim_position(pA, box, trims);
    //    trim_position(pB, box, trims);
    //    trim_position(pC, box, trims);
    //    trim_position(pD, box, trims);
    //    trim_position(pE, box, trims);

    double test = predicates::insphere_old(pA, pB, pD, pC, pE);
    if(testdouble) {
        *testdouble = test;
    }
    if(verbose) {
        cout << test << endl;
    }
    //    if(test == 0. || test == -0.){
    //        if(A > (int)cells.size()){
    //            cerr << ghosts[A-cells.size()]->get_id();
    //        } else {
    //            cerr << A;
    //        }
    //        cerr << "\t";
    //        if(B > (int)cells.size()){
    //            cerr << ghosts[B-cells.size()]->get_id();
    //        } else {
    //            cerr << B;
    //        }
    //        cerr << "\t";
    //        if(C > (int)cells.size()){
    //            cerr << ghosts[C-cells.size()]->get_id();
    //        } else {
    //            cerr << C;
    //        }
    //        cerr << "\t";
    //        if(D > (int)cells.size()){
    //            cerr << ghosts[D-cells.size()]->get_id();
    //        } else {
    //            cerr << D;
    //        }
    //        cerr << "\t";
    //        if(E > (int)cells.size()){
    //            cerr << ghosts[E-cells.size()]->get_id();
    //        } else {
    //            cerr << E;
    //        }
    //        cerr << endl;
    //        cerr << pA[0] << "\t" << pA[1] << "\t" << pA[2] << endl;
    //        cerr << pB[0] << "\t" << pB[1] << "\t" << pB[2] << endl;
    //        cerr << pC[0] << "\t" << pC[1] << "\t" << pC[2] << endl;
    //        cerr << pD[0] << "\t" << pD[1] << "\t" << pD[2] << endl;
    //        cerr << pE[0] << "\t" << pE[1] << "\t" << pE[2] << endl;
    //        cerr << test << endl;
    //        cout << "Exactly zero!" << endl;
    //        ofstream ofile("doublengb.dat");
    //        print_relevant_info(ofile, cells, ghosts, box);
    //        ofile.close();
    //        exit(1);
    //    }
    return test > 0.;
}

/**
 * @brief Loop over all faces of the cell and count the edge flips for every
 * face
 *
 * Possible return values:
 *  - 0: no edge flips
 *  - 1: one edge flip
 *  - 2: face flip
 *  - 3: too many flips: invalid cell
 *  - 4: double edge flip in face with four edges
 *
 * @param id Integer index of the cell in the AdaptiveVorTess3d cell list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 * @return Unsigned integer in the range 0-4 specifying the result of the tests
 */
unsigned int AdaptiveVorCell3d::get_flips(
        int id, std::vector<AdaptiveVorCell3d*>& cells,
        std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box, bool periodic) {
    vector<unsigned int> numflips(_facengbs.size(), 0);
    //    _flips.clear();
    //    _flips.resize(_facengbs.size());
    unsigned int flip1[2] = {0, 0};
    unsigned int flip3 = 0;
    unsigned int flip4[2] = {0, 0};
    bool valid4 = false;
    for(unsigned int i = 0; i < _facengbs.size(); i++) {
        //        _flips[i].resize(_facengbs[i].size());
        //        if(id > _ngbs[i]){
        //            continue;
        //        }
        for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
            unsigned int ni = (j + 1) % _facengbs[i].size();
            unsigned int pi = (j + 2) % _facengbs[i].size();
            int cA = id;
            int cB = _ngbs[i];
            int cC = _facengbs[i][ni];
            int cD = _facengbs[i][j];
            int cE = _facengbs[i][pi];
            int cwA = 0;
            int cwB = _ngbwalls[i];
            int cwC = _facengbwalls[i][ni];
            int cwD = _facengbwalls[i][j];
            int cwE = _facengbwalls[i][pi];
            bool test = in_sphere_test(cA, cB, cC, cD, cE, cwA, cwB, cwC, cwD,
                                       cwE, cells, ghosts, box, periodic);
            //            _flips[i][j] = test;
            if(test) {
                numflips[i]++;
                if(numflips[i] == 2) {
                    if(flip1[1] == j - 2) {
                        flip4[0] = i;
                        flip4[1] = j;
                        valid4 = true;
                    }
                }
                if(numflips[i] == 3) {
                    if(!valid4) {
                        flip4[0] = i;
                        flip4[1] = j;
                        valid4 = true;
                    }
                }
                flip1[0] = i;
                flip1[1] = j;
                if(numflips[i] > 3) {
                    return 3;
                }
            }
        }
    }
    // count the total number of single flips and triple flips
    unsigned int num1 = 0;
    unsigned int num4 = 0;
    unsigned int num3 = 0;
    for(unsigned int i = 0; i < numflips.size(); i++) {
        if(numflips[i]) {
            if(numflips[i] == 1) {
                num1++;
            } else {
                if(numflips[i] == 3) {
                    if(_facengbs[i].size() == 3) {
                        num3++;
                        flip3 = i;
                    } else {
                        if(_facengbs[i].size() == 4) {
                            num4++;
                        } else {
                            // a triple flip is only valid if the face has only
                            // 3
                            // vertices
                            // if not, this cell is invalid (return value 3)
                            //                        cerr << "here" << endl;
                            return 3;
                        }
                    }
                } else {
                    if(numflips[i] == 2 && _facengbs[i].size() == 4 && valid4) {
                        num4++;
                    } else {
                        //                        cerr << numflips[i] << endl;
                        //                        cerr << _facengbs[i].size() <<
                        // endl;
                        //                    cerr << _flips[i][0] << "\t" <<
                        // _flips[i][1] << "\t"
                        //                        << _flips[i][2] << "\t" <<
                        // _flips[i][3] << endl;
                        //                    cerr << _facengbs[i][0] << "\t" <<
                        // _facengbs[i][1] << "\t"
                        //                        << _facengbs[i][2] << "\t" <<
                        // _facengbs[i][3] << endl;
                        //                        cerr << "heri" << endl;
                        return 3;
                    }
                }
            }
        }
    }
    //    cerr << "num1: " << num1 << ", num3: " << num3 << ", num4: " << num4
    //    << endl;
    // interprete the results
    //    if(num4){
    //        if(!num3 && num1 < 3){
    //            _flip_info[0] = flip4[0];
    //            _flip_info[1] = flip4[1];
    //            return 4;
    //        } else {
    //            return 3;
    //        }
    //    }
    //    if(num3){
    //        if(num1 < 4){
    //            _flip_info[0] = flip3;
    //            return 2;
    //        } else {
    //            return 3;
    //        }
    //    }
    //    if(num1){
    //        if(num1 < 4){
    //            _flip_info[0] = flip1[0];
    //            _flip_info[1] = flip1[1];
    //            return 1;
    //        } else {
    //            return 3;
    //        }
    //    }
    //    return 0;
    if(num4 && num1 == 2 && num3 == 0) {
        _flip_info[0] = flip4[0];
        _flip_info[1] = flip4[1];
        //        cerr << _ngbs[_flip_info[0]] << "\t"
        //        << _facengbs[_flip_info[0]][_flip_info[1]] << "\t";
        //        if(_flip_info[1] > 1){
        //            cerr << _facengbs[_flip_info[0]][_flip_info[1]-2] << endl;
        //        } else {
        //            cerr << _facengbs[_flip_info[0]][_flip_info[1]+2] << endl;
        //        }
        return 4;
    }
    if(num4 && num1 == 3 && num3 == 0) {
        _flip_info[0] = flip4[0];
        _flip_info[1] = flip4[1];
        return 4;
    }
    if(num1 == 0 && num3 == 0 && num4 == 0) {
        // no flips
        return 0;
    }
    if(num1 == 2 && num3 == 0) {
        // edge flip
        _flip_info[0] = flip1[0];
        _flip_info[1] = flip1[1];
        return 1;
    }
    if(num1 == 3 && num3 == 1) {
        // face flip
        _flip_info[0] = flip3;
        return 2;
    }
    // too many flips: invalid cell
    return 3;
}

/**
 * @brief Loop over all faces of the cell and count the edge flips for every
 * face that is somehow connected to the given neighbour cell
 *
 * Possible return values:
 *  - 0: no edge flips
 *  - 1: one edge flip
 *  - 2: face flip
 *  - 3: too many flips: invalid cell
 *  - 4: double edge flip in face with four edges
 *
 * @param id Integer index of the cell in the AdaptiveVorTess3d cell list
 * @param ngb Integer index of a neighbour cell in the AdaptiveVorTess3d cell
 * list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 * @return Unsigned integer in the range 0-4 specifying the result of the tests
 */
unsigned int AdaptiveVorCell3d::get_flips(
        int id, int ngb, std::vector<AdaptiveVorCell3d*>& cells,
        std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box, bool periodic) {
    vector<unsigned int> numflips(_facengbs.size(), 0);
    unsigned int flip1[2] = {0, 0};
    unsigned int flip3 = 0;
    unsigned int flip4[2] = {0, 0};
    bool valid4 = false;
    for(unsigned int i = 0; i < _facengbs.size(); i++) {
        if(_ngbs[i] != ngb) {
            bool is_ngb = false;
            for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
                if(_facengbs[i][j] == ngb) {
                    is_ngb = true;
                    break;
                }
            }
            if(!is_ngb) {
                continue;
            }
        }
        for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
            unsigned int ni = (j + 1) % _facengbs[i].size();
            unsigned int pi = (j + 2) % _facengbs[i].size();
            int cA = id;
            int cB = _ngbs[i];
            int cC = _facengbs[i][ni];
            int cD = _facengbs[i][j];
            int cE = _facengbs[i][pi];
            int cwA = 0;
            int cwB = _ngbwalls[i];
            int cwC = _facengbwalls[i][ni];
            int cwD = _facengbwalls[i][j];
            int cwE = _facengbwalls[i][pi];
            bool test = in_sphere_test(cA, cB, cC, cD, cE, cwA, cwB, cwC, cwD,
                                       cwE, cells, ghosts, box, periodic);
            if(test) {
                numflips[i]++;
                if(numflips[i] == 2) {
                    if(flip1[1] == j - 2) {
                        flip4[0] = i;
                        flip4[1] = j;
                        valid4 = true;
                    }
                }
                flip1[0] = i;
                flip1[1] = j;
                if(numflips[i] > 3) {
                    return 3;
                }
            }
        }
    }
    // count the total number of single flips and triple flips
    unsigned int num1 = 0;
    unsigned int num4 = 0;
    unsigned int num3 = 0;
    for(unsigned int i = 0; i < numflips.size(); i++) {
        if(numflips[i]) {
            if(numflips[i] == 1) {
                num1++;
            } else {
                if(numflips[i] == 3) {
                    if(_facengbs[i].size() == 3) {
                        num3++;
                        flip3 = i;
                    } else {
                        // a triple flip is only valid if the face has only 3
                        // vertices
                        // if not, this cell is invalid (return value 3)
                        //                        cerr << "here" << endl;
                        return 3;
                    }
                } else {
                    if(numflips[i] == 2 && _facengbs[i].size() == 4 && valid4) {
                        num4++;
                    } else {
                        //                        cerr << numflips[i] << endl;
                        //                        cerr << _facengbs[i].size() <<
                        // endl;
                        //                    cerr << _flips[i][0] << "\t" <<
                        // _flips[i][1] << "\t"
                        //                        << _flips[i][2] << "\t" <<
                        // _flips[i][3] << endl;
                        //                    cerr << _facengbs[i][0] << "\t" <<
                        // _facengbs[i][1] << "\t"
                        //                        << _facengbs[i][2] << "\t" <<
                        // _facengbs[i][3] << endl;
                        //                        cerr << "heri" << endl;
                        return 3;
                    }
                }
            }
        }
    }
    if(num4 && num1 == 2 && num3 == 0) {
        _flip_info[0] = flip4[0];
        _flip_info[1] = flip4[1];
        //        cerr << _ngbs[_flip_info[0]] << "\t"
        //        << _facengbs[_flip_info[0]][_flip_info[1]] << "\t";
        //        if(_flip_info[1] > 1){
        //            cerr << _facengbs[_flip_info[0]][_flip_info[1]-2] << endl;
        //        } else {
        //            cerr << _facengbs[_flip_info[0]][_flip_info[1]+2] << endl;
        //        }
        return 4;
    }
    if(num1 == 0 && num3 == 0 && num4 == 0) {
        // no flips
        return 0;
    }
    if(num1 == 2 && num3 == 0) {
        // edge flip
        _flip_info[0] = flip1[0];
        _flip_info[1] = flip1[1];
        return 1;
    }
    if(num1 == 3 && num3 == 1) {
        // face flip
        _flip_info[0] = flip3;
        return 2;
    }
    // too many flips: invalid cell
    return 3;
}

/**
 * @brief Loop over all faces of the cell and count the edge flips for every
 * face, only considering neighbours with lower indices in the AdaptiveVorTess3d
 * cell list
 *
 * Possible return values:
 *  - 0: no edge flips
 *  - 1: one edge flip
 *  - 2: face flip
 *  - 3: too many flips: invalid cell
 *  - 4: double edge flip in face with four edges
 *
 * Information about the encountered flips is also stored in an array, allowing
 * to count the total number of flips of every type.
 *
 * @param flips Array to store information about the encountered flips
 * @param id Integer index of the cell in the AdapativeVorTess3d cell list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 * @return A value in the range 0-4 indicating the result of the tests
 */
unsigned int AdaptiveVorCell3d::get_flips_minimal(
        unsigned int* flips, int id, std::vector<AdaptiveVorCell3d*>& cells,
        std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box, bool periodic) {
    vector<unsigned int> numflips(_facengbs.size(), 0);
    unsigned int flip1[2] = {0, 0};
    unsigned int flip3 = 0;
    unsigned int flip4[2] = {0, 0};
    bool valid4 = false;
    for(unsigned int i = 0; i < _facengbs.size(); i++) {
        if(id > _ngbs[i]) {
            continue;
        }
        for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
            unsigned int ni = (j + 1) % _facengbs[i].size();
            unsigned int pi = (j + 2) % _facengbs[i].size();
            int cA = id;
            int cB = _ngbs[i];
            int cC = _facengbs[i][ni];
            int cD = _facengbs[i][j];
            int cE = _facengbs[i][pi];
            int cwA = 0;
            int cwB = _ngbwalls[i];
            int cwC = _facengbwalls[i][ni];
            int cwD = _facengbwalls[i][j];
            int cwE = _facengbwalls[i][pi];
            bool test = in_sphere_test(cA, cB, cC, cD, cE, cwA, cwB, cwC, cwD,
                                       cwE, cells, ghosts, box, periodic);
            if(test) {
                numflips[i]++;
                if(numflips[i] == 2) {
                    if(flip1[1] == j - 2) {
                        flip4[0] = i;
                        flip4[1] = j;
                        valid4 = true;
                    }
                }
                flip1[0] = i;
                flip1[1] = j;
                if(numflips[i] > 3) {
                    flips[3]++;
                    //                    cerr << "quit" << endl;
                    //                    cerr << numflips[i] << ", " <<
                    // _facengbs[i].size()
                    //                    << endl;
                    return 3;
                }
            }
        }
    }
    // count the total number of single flips and triple flips
    unsigned int num1 = 0;
    unsigned int num4 = 0;
    unsigned int num3 = 0;
    for(unsigned int i = 0; i < numflips.size(); i++) {
        if(numflips[i]) {
            if(numflips[i] == 1) {
                num1++;
            } else {
                if(numflips[i] == 3) {
                    if(_facengbs[i].size() == 3) {
                        num3++;
                        flip3 = i;
                    } else {
                        flips[3]++;
                        //                        cerr << "quit" << endl;
                        //                        cerr << numflips[i] << endl;
                        //                        cerr << _facengbs[i].size() <<
                        // endl;
                        return 3;
                    }
                } else {
                    if(numflips[i] == 2 && _facengbs[i].size() == 4 && valid4) {
                        num4++;
                    } else {
                        flips[3]++;
                        //                        cerr << "quit" << endl;
                        //                        cerr << numflips[i] << endl;
                        //                        cerr << _facengbs[i].size() <<
                        // endl;
                        return 3;
                    }
                }
            }
        }
    }
    //    cerr << "num1: " << num1 << ", num3: " << num3 << ", num4: " << num4
    //    << endl;
    unsigned int returnval = 0;
    if(num4) {
        _flip_info[0] = flip4[0];
        _flip_info[1] = flip4[1];
        returnval = 4;
    } else {
        if(num3) {
            _flip_info[0] = flip3;
            returnval = 2;
        } else {
            if(num1) {
                _flip_info[0] = flip1[0];
                _flip_info[1] = flip1[1];
                returnval = 1;
            }
        }
    }
    if(num1 || num3 || num4) {
        flips[1] += num1;
        flips[2] += num3;
        flips[4] += num4;
    } else {
        flips[0]++;
    }
    return returnval;
}

/**
 * @brief Loop over all faces of the cell and count the edge flips for every
 * face, only considering neighbours with lower indices in the AdaptiveVorTess3d
 * cell list that are somehow connected to the given neighbour
 *
 * Possible return values:
 *  - 0: no edge flips
 *  - 1: one edge flip
 *  - 2: face flip
 *  - 3: too many flips: invalid cell
 *  - 4: double edge flip in face with four edges
 *
 * Information about the encountered flips is also stored in an array, allowing
 * to count the total number of flips of every type.
 *
 * @param flips Array to store information about the encountered flips
 * @param id Integer index of the cell in the AdaptiveVorTess3d cell list
 * @param ngb Integer index of a neighbouring cell in the AdaptiveVorTess3d cell
 * list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 * @return A value in the range 0-4 indicating the result of the tests
 */
unsigned int AdaptiveVorCell3d::get_flips_minimal(
        unsigned int* flips, int id, int ngb,
        std::vector<AdaptiveVorCell3d*>& cells,
        std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box, bool periodic) {
    vector<unsigned int> numflips(_facengbs.size(), 0);
    unsigned int flip1[2] = {0, 0};
    unsigned int flip3 = 0;
    unsigned int flip4[2] = {0, 0};
    bool valid4 = false;
    for(unsigned int i = 0; i < _facengbs.size(); i++) {
        if(_ngbs[i] != ngb) {
            bool is_ngb = false;
            for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
                if(_facengbs[i][j] == ngb) {
                    is_ngb = true;
                    break;
                }
            }
            if(!is_ngb) {
                continue;
            }
        }
        if(id > _ngbs[i]) {
            continue;
        }
        for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
            unsigned int ni = (j + 1) % _facengbs[i].size();
            unsigned int pi = (j + 2) % _facengbs[i].size();
            int cA = id;
            int cB = _ngbs[i];
            int cC = _facengbs[i][ni];
            int cD = _facengbs[i][j];
            int cE = _facengbs[i][pi];
            int cwA = 0;
            int cwB = _ngbwalls[i];
            int cwC = _facengbwalls[i][ni];
            int cwD = _facengbwalls[i][j];
            int cwE = _facengbwalls[i][pi];
            bool test = in_sphere_test(cA, cB, cC, cD, cE, cwA, cwB, cwC, cwD,
                                       cwE, cells, ghosts, box, periodic);
            if(test) {
                numflips[i]++;
                if(numflips[i] == 2) {
                    if(flip1[1] == j - 2) {
                        flip4[0] = i;
                        flip4[1] = j;
                        valid4 = true;
                    }
                }
                flip1[0] = i;
                flip1[1] = j;
                if(numflips[i] > 3) {
                    flips[3]++;
                    return 3;
                }
            }
        }
    }
    // count the total number of single flips and triple flips
    unsigned int num1 = 0;
    unsigned int num4 = 0;
    unsigned int num3 = 0;
    for(unsigned int i = 0; i < numflips.size(); i++) {
        if(numflips[i]) {
            if(numflips[i] == 1) {
                num1++;
            } else {
                if(numflips[i] == 3) {
                    if(_facengbs[i].size() == 3) {
                        num3++;
                        flip3 = i;
                    } else {
                        flips[3]++;
                        return 3;
                    }
                } else {
                    if(numflips[i] == 2 && _facengbs[i].size() == 4 && valid4) {
                        num4++;
                    } else {
                        flips[3]++;
                        return 3;
                    }
                }
            }
        }
    }
    unsigned int returnval = 0;
    if(num4) {
        _flip_info[0] = flip4[0];
        _flip_info[1] = flip4[1];
        returnval = 4;
    } else {
        if(num3) {
            _flip_info[0] = flip3;
            returnval = 2;
        } else {
            if(num1) {
                _flip_info[0] = flip1[0];
                _flip_info[1] = flip1[1];
                returnval = 1;
            }
        }
    }
    if(num1 || num3 || num4) {
        flips[1] += num1;
        flips[2] += num3;
        flips[4] += num4;
    } else {
        flips[0]++;
    }
    return returnval;
}

/**
 * @brief Old kernel of the mesh evolution algorithm
 *
 * Replaced by explicit flip counting and per cell restoration in
 * AdaptiveVorTess3d
 *
 * @param id Unsigned integer index of the cell in the AdaptiveVorTess3d cell
 * list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost dcell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 * @return True if errors in the mesh were detected and restored, false
 * otherwise
 */
bool AdaptiveVorCell3d::detect_crossovers(unsigned int id,
                                          vector<AdaptiveVorCell3d*>& cells,
                                          vector<AdaptiveVorCell3d*>& ghosts,
                                          Cuboid& box, bool periodic) {
    // for every face: check if the face flipped
    // hmm, we actually already did this...
    for(unsigned int i = 0; i < _facengbs.size(); i++) {
        // not entirely sure this is correct, since a previous operation could
        // have invalidated a valid face
        //        if(!_numflips[i]){
        //            continue;
        //        }
        for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
            // speaking in edge terms, j corresponds to the previous edge
            // j+3 (=ni) is the edge itself (and the index of the
            // corresponding face neighbour in the list) and j+6 (=pi) is
            // the next edge
            unsigned int ni = (j + 1) % _facengbs[i].size();
            unsigned int pi = (j + 2) % _facengbs[i].size();
            // we check if the generator of the next edge lies inside the
            // sphere through this generator, the neighbour associated with
            // this face, the generator of this edge and the generator of
            // the previous edge
            int cA = id;
            int cB = _ngbs[i];
            int cC = _facengbs[i][ni];
            int cD = _facengbs[i][j];
            int cE = _facengbs[i][pi];
            int cwA = 0;
            int cwB = _ngbwalls[i];
            int cwC = _facengbwalls[i][ni];
            int cwD = _facengbwalls[i][j];
            int cwE = _facengbwalls[i][pi];
            bool test = in_sphere_test(cA, cB, cC, cD, cE, cwA, cwB, cwC, cwD,
                                       cwE, cells, ghosts, box, periodic);
            if(test) {
                // we check if we deal with an insertion or a removal
                // a removal occurs when this face or one of its
                // neighbouring faces (in this cell and other cells) has
                // only three vertices
                int cellself = id;
                int cellngb = _ngbs[i];
                int cellother = _facengbs[i][ni];
                unsigned int removal = check_triangular_face(cellself, cellngb,
                                                             cellother, cells);
                int A, B, C, D, E;
                int wallA, wallB, wallC, wallD, wallE;
                bool status;
                if(removal) {
                    // for closer inspection
                    //                        cerr << "Face removal! (case = "
                    // << removal << ")"
                    //                    << endl;
                    // names depend on case at hand...
                    if(removal == 1) {
                        A = _facengbs[i][j];
                        B = _facengbs[i][pi];
                        C = id;
                        D = _ngbs[i];
                        E = _facengbs[i][ni];
                        wallA = _facengbwalls[i][j];
                        wallB = _facengbwalls[i][pi];
                        wallC = 0;
                        wallD = _ngbwalls[i];
                        wallE = _facengbwalls[i][ni];
                    } else {
                        if(removal == 2) {
                            A = _facengbs[i][pi];
                            B = _facengbs[i][j];
                            C = id;
                            D = _facengbs[i][ni];
                            E = _ngbs[i];
                            wallA = _facengbwalls[i][pi];
                            wallB = _facengbwalls[i][j];
                            wallC = 0;
                            wallD = _facengbwalls[i][ni];
                            wallE = _ngbwalls[i];
                        } else {
                            A = _facengbs[i][j];
                            B = _facengbs[i][pi];
                            C = _ngbs[i];
                            D = _facengbs[i][ni];
                            E = id;
                            wallA = _facengbwalls[i][j];
                            wallB = _facengbwalls[i][pi];
                            wallC = _ngbwalls[i];
                            wallD = _facengbwalls[i][ni];
                            wallE = 0;
                        }
                    }
                    unsigned int edges_to_check[126];
                    status = remove_face(edges_to_check, A, B, C, D, E, wallA,
                                         wallB, wallC, wallD, wallE, cells,
                                         ghosts, periodic);
                } else {
                    // for closer inspection
                    //                        cerr << "Face insertion!" << endl;
                    A = id;
                    B = _facengbs[i][ni];
                    C = _facengbs[i][j];
                    D = _facengbs[i][pi];
                    E = _ngbs[i];
                    wallA = 0;
                    wallB = _facengbwalls[i][ni];
                    wallC = _facengbwalls[i][j];
                    wallD = _facengbwalls[i][pi];
                    wallE = _ngbwalls[i];
                    unsigned int edges_to_check[162];
                    status = create_new_face(edges_to_check, A, B, C, D, E,
                                             wallA, wallB, wallC, wallD, wallE,
                                             cells, ghosts, periodic);
                }
                if(status) {
                    cells[A]->detect_crossovers(A, cells, ghosts, box,
                                                periodic);
                    cells[B]->detect_crossovers(B, cells, ghosts, box,
                                                periodic);
                    cells[C]->detect_crossovers(C, cells, ghosts, box,
                                                periodic);
                    cells[D]->detect_crossovers(D, cells, ghosts, box,
                                                periodic);
                    cells[E]->detect_crossovers(E, cells, ghosts, box,
                                                periodic);
                }
                return status;
            }
        }
    }
    return false;
}

/**
 * @brief Print the given neighbours of the cell in plottable format to the
 * given stream
 *
 * Useful for debugging purposes, but crashes if one of the requested neighbours
 * is not a neighbour.
 *
 * @param stream std::ostream to write to
 * @param ngblist List of integer indices of neighbours to print
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell3d::print_copies(ostream& stream, vector<int>& ngblist,
                                     vector<AdaptiveVorCell3d*>& cells,
                                     Cuboid& box) {
    print_copy(stream, 0, 0, box);
    int copy_id = 1;
    for(unsigned int i = 0; i < ngblist.size(); i++) {
        unsigned int id = 0;
        while(id < _ngbs.size() && _ngbs[id] != ngblist[i]) {
            id++;
        }
        if(id == _ngbs.size()) {
            cerr << "Error!" << endl;
            continue;
            //            exit(1);
        }
        cells[_ngbs[id]]->print_copy(stream, copy_id, _ngbwalls[id], box);
        copy_id++;
    }
}

/**
 * @brief Replace boundary flags for the given neighbour taking into account a
 * periodic movement of the neighbour through the given boundary
 *
 * @param id Integer index of the cell in the AdaptiveVorTess3d cell list
 * @param wall Boundary flag for the boundary through which the neighbour moves
 */
void AdaptiveVorCell3d::swap_wallpos(int id, int wall) {
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        if(_ngbs[i] == id) {
            _ngbwalls[i] =
                    AdaptiveMeshUtils::replace_wallpos(_ngbwalls[i], wall);
        }
        for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
            if(_facengbs[i][j] == id) {
                _facengbwalls[i][j] = AdaptiveMeshUtils::replace_wallpos(
                        _facengbwalls[i][j], wall);
            }
        }
    }
}

/**
 * @brief Keep the generator of the cell inside the periodic box
 *
 * If necessary, coordinates and boundary flags are adjusted.
 *
 * @param id Integer index of the cell in the AdaptiveVorTess3d cell list
 * @param cells Reference to the AdapativeVorTess3d cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell3d::keep_inside(int id, vector<AdaptiveVorCell3d*>& cells,
                                    Cuboid& box) {
    double dx[3];
    dx[0] = _pos[0];
    dx[1] = _pos[1];
    dx[2] = _pos[2];
    int wall = AdaptiveMeshUtils::trim_wallpos(_pos, box);
    if(wall < 0) {
        for(unsigned int i = 0; i < _ngbs.size(); i++) {
            cells[_ngbs[i]]->swap_wallpos(id, wall);
            _ngbwalls[i] = AdaptiveMeshUtils::swap_wallpos(_ngbwalls[i], wall);
            for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
                _facengbwalls[i][j] = AdaptiveMeshUtils::swap_wallpos(
                        _facengbwalls[i][j], wall);
            }
        }
        dx[0] -= _pos[0];
        dx[1] -= _pos[1];
        dx[2] -= _pos[2];
        _valid_pos[0] -= dx[0];
        _valid_pos[1] -= dx[1];
        _valid_pos[2] -= dx[2];
        _particle->set_x(_pos[0]);
        _particle->set_y(_pos[1]);
        _particle->set_z(_pos[2]);
        Vec centroid = _particle->get_centroid();
        centroid[0] -= dx[0];
        centroid[1] -= dx[1];
        centroid[2] -= dx[2];
        _particle->set_centroid(centroid);
    }
}

/**
 * @brief Get the indices of the neighbours of this cell
 *
 * @return std::vector of integer indices
 */
vector<int> AdaptiveVorCell3d::get_ngbs() {
    return _ngbs;
}

/**
 * @brief Get the characteristic radius of the cell
 *
 * The characteristic radius is defined as the radius of a sphere with the same
 * volume as the cell.
 *
 * @return The characteristic radius of the cell
 */
double AdaptiveVorCell3d::get_h() {
    double V = get_volume();
    double h = cbrt(V * 3. / 4. / M_PI);
    return h;
}

/**
 * @brief Calculate the volume of the cell
 *
 * The total volume is the sum of the volumes of the tetrahedra formed by the
 * generator of the cell and three consecutive vertices of one of its faces.
 *
 * @return The volume of the cell
 */
double AdaptiveVorCell3d::get_volume() {
    double V = 0;
    for(unsigned int i = 0; i < _faces.size(); i++) {
        vector<double> pos = _faces[i]->get_vertices();
        for(unsigned int j = 3; j < pos.size() - 3; j += 3) {
            double A;
            double r1[3], r2[3], r3[3];
            r1[0] = pos[0] - pos[j];
            r1[1] = pos[1] - pos[j + 1];
            r1[2] = pos[2] - pos[j + 2];
            r2[0] = pos[j] - pos[j + 3];
            r2[1] = pos[j + 1] - pos[j + 4];
            r2[2] = pos[j + 2] - pos[j + 5];
            r3[0] = pos[j + 3] - _pos[0];
            r3[1] = pos[j + 4] - _pos[1];
            r3[2] = pos[j + 5] - _pos[2];
            A = fabs(r1[0] * r2[1] * r3[2] + r1[1] * r2[2] * r3[0] +
                     r1[2] * r2[0] * r3[1] - r1[2] * r2[1] * r3[0] -
                     r2[2] * r3[1] * r1[0] - r3[2] * r1[1] * r2[0]);
            A /= 6;
            V += A;
        }
    }
    return V;
}

/**
 * @brief Calculate the centroid of the cell
 *
 * The centroid of the cell is the volume average of the centroids of the
 * tetrahedra formed by the cell generator and three consecutive vertices of one
 * of its faces.
 *
 * @return The centroid of the cell
 */
Vec AdaptiveVorCell3d::get_centroid() {
    double C_tot[3], A_tot;
    C_tot[0] = 0;
    C_tot[1] = 0;
    C_tot[2] = 0;
    A_tot = 0;
    for(unsigned int i = 0; i < _faces.size(); i++) {
        vector<double> pos = _faces[i]->get_vertices();
        for(unsigned int j = 3; j < pos.size() - 3; j += 3) {
            double C[3], A;
            C[0] = 0.25 * (pos[0] + pos[j] + pos[j + 3] + _pos[0]);
            C[1] = 0.25 * (pos[1] + pos[j + 1] + pos[j + 4] + _pos[1]);
            C[2] = 0.25 * (pos[2] + pos[j + 2] + pos[j + 5] + _pos[2]);
            double r1[3], r2[3], r3[3];
            r1[0] = pos[0] - pos[j];
            r1[1] = pos[1] - pos[j + 1];
            r1[2] = pos[2] - pos[j + 2];
            r2[0] = pos[j] - pos[j + 3];
            r2[1] = pos[j + 1] - pos[j + 4];
            r2[2] = pos[j + 2] - pos[j + 5];
            r3[0] = pos[j + 3] - _pos[0];
            r3[1] = pos[j + 4] - _pos[1];
            r3[2] = pos[j + 5] - _pos[2];
            A = fabs(r1[0] * r2[1] * r3[2] + r1[1] * r2[2] * r3[0] +
                     r1[2] * r2[0] * r3[1] - r1[2] * r2[1] * r3[0] -
                     r2[2] * r3[1] * r1[0] - r3[2] * r1[1] * r2[0]);
            A /= 6;
            C_tot[0] += C[0] * A;
            C_tot[1] += C[1] * A;
            C_tot[2] += C[2] * A;
            A_tot += A;
        }
    }
    C_tot[0] /= A_tot;
    C_tot[1] /= A_tot;
    C_tot[2] /= A_tot;
    Vec centroid(C_tot[0], C_tot[1], C_tot[2]);
    return centroid;
}

/**
 * @brief Get the GasParticle associated with this cell
 *
 * @return The GasParticle associated with this cell
 */
GasParticle* AdaptiveVorCell3d::get_particle() {
    return _particle;
}

/**
 * @brief Get the face with the given index
 *
 * @param i Unsigned integer index of a face in the internal list
 * @return AdaptiveVorFace3d corresponding to the given index
 */
AdaptiveVorFace3d* AdaptiveVorCell3d::get_face(unsigned int i) {
    return _faces[i];
}

/**
 * @brief Get the position of the ghost copy of this cell
 *
 * @warning Not implemented yet!
 *
 * @param pos Array to store the result in
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell3d::get_ghost_position(double* pos, Cuboid& box) {
    // to implement
}

/**
 * @brief Get the position of the periodic copy of the given neighbour of this
 * cell
 *
 * @param i Unsigned integer index of a neighbour in the internal list
 * @param cells Reference to the AdaptiveVorTess3d cel list
 * @param posvec Array to store the result in
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell3d::get_periodic_position(
        unsigned int i, std::vector<AdaptiveVorCell3d*>& cells, double* posvec,
        Cuboid& box) {
    AdaptiveVorCell3d* ngb = cells[_ngbs[i]];
    posvec[0] = ngb->_pos[0];
    posvec[1] = ngb->_pos[1];
    posvec[2] = ngb->_pos[2];
    if(_ngbwalls[i] < 0) {
        AdaptiveMeshUtils::get_periodic_position(posvec, _ngbwalls[i], box);
    }
}

/**
 * @brief Calculate the velocity for this cell
 *
 * @param particle GasParticle associated with this cell
 * @return Velocity for this cell
 */
Vec AdaptiveVorCell3d::get_velocity(GasParticle* particle) {
    StateVector W = particle->get_Wvec();
    // Springel correction
    Vec vel;
    if(particle->get_Wvec().rho()) {
        Vec centroid = get_centroid();
        Vec d = centroid - particle->get_position();
        double Ri = get_volume();
        //    double csnd = W.soundspeed();
        double csnd = particle->get_soundspeed();
        Ri = cbrt(3. * Ri / 4. / M_PI);
        double check = 4. * d.norm() / Ri;
        if(check > 0.9) {
            if(check < 1.1) {
                vel = csnd / d.norm() * (d.norm() - 0.9 * 0.25 * Ri) /
                      (0.2 * 0.25 * Ri) * d;
            } else {
                vel = csnd / d.norm() * d;
            }
        }
    }
    // Vogelsberger correction
    //    Vec vel;
    //    if(_central_point->get_particle()->get_Wvec().rho()){
    //        double maxfaceangle = 0.;
    //        for(unsigned int i = 0; i < _faces.size(); i++){
    //            maxfaceangle = std::max(maxfaceangle,
    //    sqrt(_faces[i]->get_area()/M_PI)/(_central_point->get_position()-
    //    _faces[i]->get_midpoint()).norm());
    //        }
    //        if(maxfaceangle > 1.68){
    //            Vec d = get_centroid() - _central_point->get_position();
    //            double csnd = ((GasParticle*)_central_point
    //    ->get_particle())->get_soundspeed();
    //            vel = csnd/d.norm()*d;
    //        }
    //    }
    //    Vec vel = _central_point->get_particle()->get_mesh_v();
    return Vec(W[1] + vel[0], W[2] + vel[1], W[3] + vel[2]);
}

/**
 * @brief Estimate gradients for the hydrodynamical quantities of this cell
 *
 * @param delta Array to store the results in
 * @param particle GasParticle associated with this cell
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
void AdaptiveVorCell3d::estimate_gradients(
        StateVector* delta, GasParticle* particle,
        std::vector<AdaptiveVorCell3d*>& cells,
        std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box, bool periodic) {
    StateVector W = particle->get_Wvec();
    // gradient estimation
    StateVector Wmaxvec, Wminvec;
    Wmaxvec = W;
    Wminvec = W;
    for(unsigned int i = 0; i < _faces.size(); i++) {
        double Aij = _faces[i]->get_area();
        if(Aij) {
            double midfacevec[3];
            _faces[i]->get_midpoint(midfacevec);
            AdaptiveVorCell3d* ngb;
            double ngbposvec[3];
            if(_ngbs[i] < (int)cells.size()) {
                ngb = cells[_ngbs[i]];
                ngbposvec[0] = ngb->_pos[0];
                ngbposvec[1] = ngb->_pos[1];
                ngbposvec[2] = ngb->_pos[2];
                if(periodic && _ngbwalls[i] < 0) {
                    AdaptiveMeshUtils::get_periodic_position(ngbposvec,
                                                             _ngbwalls[i], box);
                }
            } else {
                ngb = ghosts[_ngbs[i] - cells.size()];
                AdaptiveMeshUtils::get_wall_position(ngb->_pos, ngbposvec,
                                                     ngb->_ngbs[0], box);
            }
            Vec midface(midfacevec[0], midfacevec[1], midfacevec[2]);
            Vec ngbpos(ngbposvec[0], ngbposvec[1], ngbposvec[2]);
            Vec c = midface - 0.5 * (particle->get_position() + ngbpos);
            Vec rij = particle->get_position() - ngbpos;
            double rnorm = rij.norm();
            StateVector Wj;
            // check if we have to mirror this particle
            if(_ngbs[i] >= (int)cells.size()) {
                Wj = particle->get_Wvec();
                _faces[i]->transform(Wj, _pos, ngbposvec);
                Wj[1] = -Wj[1];
                _faces[i]->invtransform(Wj, _pos, ngbposvec);
            } else {
                Wj = ngb->_particle->get_Wvec();
            }
            for(unsigned int l = ndim_; l--;) {
                delta[l] += Aij * ((Wj - W) * c[l] / rnorm -
                                   0.5 * (W + Wj) * rij[l] / rnorm);
            }
            Wmaxvec = max(Wmaxvec, Wj);
            Wminvec = min(Wminvec, Wj);
        }
    }
    for(unsigned int l = ndim_; l--;) {
        delta[l] /= get_volume();
    }
    // slope limiting
    Vec centroid = get_centroid();
    StateVector alphavec(1.);
    for(unsigned int l = _faces.size(); l--;) {
        if(_faces[l]->get_area()) {
            double midfacevec[3];
            _faces[l]->get_midpoint(midfacevec);
            Vec midface(midfacevec[0], midfacevec[1], midfacevec[2]);
            Vec d = midface - centroid;
            StateVector deltap =
                    delta[0] * d[0] + delta[1] * d[1] + delta[2] * d[2];
            alphavec = min(alphavec,
                           max((Wmaxvec - W) / deltap, (Wminvec - W) / deltap));
        }
    }
    for(unsigned int m = ndim_; m--;) {
        delta[m] *= alphavec;
    }

    // due to numerical error, it can happen that the reconstructed values at
    // the faces are still smaller than the minimal value. If this happens, we
    // slightly adjust the gradients to make sure reconstructed values are
    // always positive
    // THIS CODE CAN BE FOUND IN VorCell.cpp
}

/**
 * @brief Flag this cell to be checked during the third detection iteration
 */
void AdaptiveVorCell3d::queue() {
    _flag2 = true;
}
/**
 * @brief Check if this cell has to be tested during the third detection
 * iteration
 *
 * @return True if the cell has to be tested, false otherwise
 */
bool AdaptiveVorCell3d::queued() {
    return _flag2;
}

/**
 * @brief Flag this cell to be tested during the first detection iteration
 */
void AdaptiveVorCell3d::activate() {
    _active = true;
}

/**
 * @brief Check if this cell has to be tested during the first detection
 * iteration
 *
 * @return True if the cell has to be tested, false otherwise
 */
bool AdaptiveVorCell3d::active() {
    return _active;
}

/**
 * @brief Flag this cell to be tested during the second detection iteration
 */
void AdaptiveVorCell3d::semiactivate() {
    _flag1 = true;
}

/**
 * @brief Check whether this cell has to be tested during the second detection
 * iteration
 *
 * @return True if the cell has to tested, false otherwise
 */
bool AdaptiveVorCell3d::semiactive() {
    return _flag1;
}

/**
 * @brief Reset all detection iteration flags
 */
void AdaptiveVorCell3d::reset_flags() {
    _flag1 = false;
    _flag2 = false;
    _active = false;
}

/**
 * @brief Set the last valid position of this cell
 *
 * @param valid_pos Coordinates of the last valid position for this cell
 */
void AdaptiveVorCell3d::set_valid_pos(double* valid_pos) {
    _valid_pos[0] = valid_pos[0];
    _valid_pos[1] = valid_pos[1];
    _valid_pos[2] = valid_pos[2];
}

/**
 * @brief Get the last valid position of this cell
 *
 * @param valid_pos Array to store the result in
 */
void AdaptiveVorCell3d::get_valid_pos(double* valid_pos) {
    valid_pos[0] = _valid_pos[0];
    valid_pos[1] = _valid_pos[1];
    valid_pos[2] = _valid_pos[2];
}

/**
 * @brief Remove the internally stored flipped face
 *
 * @param flips Not used
 * @param id Integer index of this cell in the AdaptiveVorTess3d cell list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
void AdaptiveVorCell3d::remove_face(unsigned int* flips, int id,
                                    std::vector<AdaptiveVorCell3d*>& cells,
                                    std::vector<AdaptiveVorCell3d*>& ghosts,
                                    Cuboid& box, bool periodic) {
    remove_face(flips, id, _flip_info[0], cells, ghosts, box, periodic);
}

/**
 * @brief Remove the face with the given index
 *
 * @param flips Not used
 * @param id Integer index of this cell in the AdaptiveVorTess3d cell list
 * @param face Unsigned integer index of a neighbour and face of this cell
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
void AdaptiveVorCell3d::remove_face(unsigned int* flips, int id,
                                    unsigned int face,
                                    std::vector<AdaptiveVorCell3d*>& cells,
                                    std::vector<AdaptiveVorCell3d*>& ghosts,
                                    Cuboid& box, bool periodic) {
    //    unsigned int i = _flip_info[0];
    unsigned int i = face;
    unsigned int j = 0;
    unsigned int ni = 1;
    unsigned int pi = 2;
    int A = _facengbs[i][j];
    int B = _facengbs[i][pi];
    int C = id;
    int D = _ngbs[i];
    int E = _facengbs[i][ni];
    int wallA = _facengbwalls[i][j];
    int wallB = _facengbwalls[i][pi];
    int wallC = 0;
    int wallD = _ngbwalls[i];
    int wallE = _facengbwalls[i][ni];
    unsigned int edges_to_check[63];
    remove_face(edges_to_check, A, B, C, D, E, wallA, wallB, wallC, wallD,
                wallE, cells, ghosts, periodic);
    int ids[5] = {A, B, C, D, E};
    check_cells(ids, 5, cells, ghosts, box, periodic);
    //    check_cells(ids, 5, edges_to_check, 63, cells, ghosts, box, periodic);
    return;
}

/**
 * @brief Insert a new face to correct for the internally stored flipped edge
 *
 * @param flips Not used
 * @param id Integer index of this cell in the AdaptiveVorTess3d cell list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
void AdaptiveVorCell3d::insert_face(unsigned int* flips, int id,
                                    std::vector<AdaptiveVorCell3d*>& cells,
                                    std::vector<AdaptiveVorCell3d*>& ghosts,
                                    Cuboid& box, bool periodic) {
    insert_face(flips, id, _flip_info[0], _flip_info[1], cells, ghosts, box,
                periodic);
}

/**
 * @brief Insert a new face to correct for the flip of the given edge in the
 * given face
 *
 * The edge is the index of the left facengb that traces out the edge if the
 * three facengbs responsible for the edge are ordered clockwise if you view the
 * face along the vector pointing from the cell generator to the neighbour that
 * shares the face.
 *
 * @warning Check if the description above is correct!
 *
 * @param flips Not used
 * @param id Integer index of this cell in the AdaptiveVorTess3d cell list
 * @param face Unsigned integer index of a neighbour and face of this cell
 * @param edge Unsigned integer index of an edge in the given face
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
void AdaptiveVorCell3d::insert_face(unsigned int* flips, int id,
                                    unsigned int face, unsigned int edge,
                                    std::vector<AdaptiveVorCell3d*>& cells,
                                    std::vector<AdaptiveVorCell3d*>& ghosts,
                                    Cuboid& box, bool periodic) {
    //    unsigned int i = _flip_info[0];
    //    unsigned int j = _flip_info[1];
    unsigned int i = face;
    unsigned int j = edge;
    unsigned int ni = (j + 1) % _facengbs[i].size();
    unsigned int pi = (j + 2) % _facengbs[i].size();
    int A = id;
    int B = _facengbs[i][ni];
    int C = _facengbs[i][j];
    int D = _facengbs[i][pi];
    int E = _ngbs[i];
    //    cerr << C << "\t" << E << endl;
    int wallA = 0;
    int wallB = _facengbwalls[i][ni];
    int wallC = _facengbwalls[i][j];
    int wallD = _facengbwalls[i][pi];
    int wallE = _ngbwalls[i];
    unsigned int edges_to_check[81];
    bool status = create_new_face(edges_to_check, A, B, C, D, E, wallA, wallB,
                                  wallC, wallD, wallE, cells, ghosts, periodic);
    if(status) {
        int ids[5] = {A, B, C, D, E};
        //        cerr << "here" << endl;
        //        check_cells(ids, 5, edges_to_check, 81, cells, ghosts, box,
        // periodic);
        check_cells(ids, 5, cells, ghosts, box, periodic);
    }
    return;
}

/**
 * @brief Flip the internally stored face
 *
 * @param id Integer index of this cell in the AdaptiveVorTess3d cell list
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
void AdaptiveVorCell3d::flip_face(int id,
                                  std::vector<AdaptiveVorCell3d*>& cells,
                                  std::vector<AdaptiveVorCell3d*>& ghosts,
                                  Cuboid& box, bool periodic) {
    flip_face(id, _flip_info[0], _flip_info[1], cells, ghosts, box, periodic);
}

/**
 * @brief Flip the given face to correct for the flip of the given edge in the
 * face
 *
 * The edge is the index of the left facengb that traces out the edge if the
 * three facengbs responsible for the edge are ordered clockwise if you view the
 * face along the vector pointing from the cell generator to the neighbour that
 * shares the face.
 *
 * If a face flips, there are always two edges that flip, which one is passed on
 * is not important for this method to work.
 *
 * @warning Check if the description above is correct!
 *
 * @param id Integer index of this cell in the AdaptiveVorTess3d cell list
 * @param face Unsigned integer index of a neighbour and face of this cell
 * @param edge Unsigned integer index of an edge in the given face
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
void AdaptiveVorCell3d::flip_face(int id, unsigned int face, unsigned int edge,
                                  std::vector<AdaptiveVorCell3d*>& cells,
                                  std::vector<AdaptiveVorCell3d*>& ghosts,
                                  Cuboid& box, bool periodic) {
    //    unsigned int a = _flip_info[1]-2;
    //    unsigned int b = _flip_info[1];
    unsigned int a, b;
    if(edge > 1) {
        a = edge - 2;
        b = edge;
    } else {
        a = edge;
        b = edge + 2;
    }
    unsigned int d, f;
    if(a % 2) {
        d = 0;
        f = 2;
    } else {
        d = 3;
        f = 1;
    }
    //    unsigned int e = _flip_info[0];
    unsigned int e = face;
    int A = _facengbs[e][a];
    int B = _facengbs[e][b];
    int wallA = _facengbwalls[e][a];
    int wallB = _facengbwalls[e][b];
    int E = _ngbs[e];
    int wallE = _ngbwalls[e];
    int D = _facengbs[e][d];
    int F = _facengbs[e][f];
    int wallD = _facengbwalls[e][d];
    int wallF = _facengbwalls[e][f];
    int C = id;
    int wallC = 0;
#ifdef VERBOSE
    cerr << "A: " << A << ", B: " << B << ", C: " << C << ", D: " << D
         << ", E: " << E << ", F: " << F << endl;
#endif
    //    ofstream ofile("face_flip.dat");
    //    cells[A]->complete(cells, ghosts, box);
    //    cells[B]->complete(cells, ghosts, box);
    //    cells[C]->complete(cells, ghosts, box);
    //    cells[D]->complete(cells, ghosts, box);
    //    cells[E]->complete(cells, ghosts, box);
    //    cells[F]->complete(cells, ghosts, box);
    //    cells[A]->print_copy(ofile, 0, wallA, box);
    //    cells[B]->print_copy(ofile, 1, wallB, box);
    //    cells[C]->print_copy(ofile, 2, wallC, box);
    //    cells[D]->print_copy(ofile, 3, wallD, box);
    //    cells[E]->print_copy(ofile, 4, wallE, box);
    //    cells[F]->print_copy(ofile, 5, wallF, box);
    //    ofile.close();
    if(cells[A]->is_ngb(B)) {
        throw AdaptiveMeshException();
        cerr << "Godverdoeme!" << endl;
        ofstream ofile("face_flip.dat");
        cells[A]->complete(cells, ghosts, box);
        cells[B]->complete(cells, ghosts, box);
        cells[C]->complete(cells, ghosts, box);
        cells[D]->complete(cells, ghosts, box);
        cells[E]->complete(cells, ghosts, box);
        cells[F]->complete(cells, ghosts, box);
        cells[A]->print_copy(ofile, 0, wallA, box);
        cells[B]->print_copy(ofile, 1, wallB, box);
        cells[C]->print_copy(ofile, 2, wallC, box);
        cells[D]->print_copy(ofile, 3, wallD, box);
        cells[E]->print_copy(ofile, 4, wallE, box);
        cells[F]->print_copy(ofile, 5, wallF, box);
        ofile.close();
        my_exit();
    }

    //    unsigned int numedge = 0;
    //    unsigned int edges_to_check[10];
    unsigned int edges[4];

    // remove the face between C and E
    cells[C]->remove_facengb(edges, A, E);
    cells[C]->remove_facengb(edges, B, E);
    cells[C]->remove_facengb(edges, D, E);
    cells[C]->remove_facengb(edges, F, E);
    cells[C]->remove_ngb(E);

    cells[E]->remove_facengb(A, C);
    cells[E]->remove_facengb(B, C);
    cells[E]->remove_facengb(D, C);
    cells[E]->remove_facengb(F, C);
    cells[E]->remove_ngb(C);

    // we do not add new facengbs, because they are removed again in the next
    // step
    cells[A]->remove_facengb(C, E);
    cells[A]->remove_facengb(E, C);

    cells[B]->remove_facengb(C, E);
    cells[B]->remove_facengb(E, C);

    cells[D]->remove_facengb(C, E);
    cells[D]->remove_facengb(E, C);

    cells[F]->remove_facengb(C, E);
    cells[F]->remove_facengb(E, C);

    // add a new face between A and B
    int wallpos;
    wallpos = AdaptiveMeshUtils::get_wallpos(wallA, wallB);
    cells[A]->add_ngb(B, wallpos);
    cells[A]->add_facengb(C, AdaptiveMeshUtils::get_wallpos(wallA, wallC));
    cells[A]->add_facengb(D, AdaptiveMeshUtils::get_wallpos(wallA, wallD));
    cells[A]->add_facengb(E, AdaptiveMeshUtils::get_wallpos(wallA, wallE));
    cells[A]->add_facengb(F, AdaptiveMeshUtils::get_wallpos(wallA, wallF));
    cells[A]->add_facengb(C, B, F, wallpos);
    cells[A]->add_facengb(D, B, C, wallpos);
    cells[A]->add_facengb(E, B, D, wallpos);
    cells[A]->add_facengb(F, B, E, wallpos);

    wallpos = AdaptiveMeshUtils::get_wallpos(wallB, wallA);
    cells[B]->add_ngb(A, wallpos);
    cells[B]->add_facengb(F, AdaptiveMeshUtils::get_wallpos(wallB, wallF));
    cells[B]->add_facengb(E, AdaptiveMeshUtils::get_wallpos(wallB, wallE));
    cells[B]->add_facengb(D, AdaptiveMeshUtils::get_wallpos(wallB, wallD));
    cells[B]->add_facengb(C, AdaptiveMeshUtils::get_wallpos(wallB, wallC));
    cells[B]->add_facengb(C, A, D, wallpos);
    cells[B]->add_facengb(D, A, E, wallpos);
    cells[B]->add_facengb(E, A, F, wallpos);
    cells[B]->add_facengb(F, A, C, wallpos);

    cells[C]->add_facengb(A, B, D,
                          AdaptiveMeshUtils::get_wallpos(wallC, wallB));
    cells[C]->add_facengb(B, A, F,
                          AdaptiveMeshUtils::get_wallpos(wallC, wallA));

    cells[D]->add_facengb(A, B, E,
                          AdaptiveMeshUtils::get_wallpos(wallD, wallB));
    cells[D]->add_facengb(B, A, C,
                          AdaptiveMeshUtils::get_wallpos(wallD, wallA));

    cells[E]->add_facengb(A, B, F,
                          AdaptiveMeshUtils::get_wallpos(wallE, wallB));
    cells[E]->add_facengb(B, A, D,
                          AdaptiveMeshUtils::get_wallpos(wallE, wallA));

    cells[F]->add_facengb(A, B, C,
                          AdaptiveMeshUtils::get_wallpos(wallF, wallB));
    cells[F]->add_facengb(B, A, E,
                          AdaptiveMeshUtils::get_wallpos(wallF, wallA));
    int ids[6] = {A, B, C, D, E, F};
    check_cells(ids, 6, cells, ghosts, box, periodic);

    //    ofile.open("face_flip_info.dat");
    //    cells[A]->print_relevant_info(ofile, cells, ghosts, box);
    //    cells[B]->print_relevant_info(ofile, cells, ghosts, box);
    //    cells[C]->print_relevant_info(ofile, cells, ghosts, box);
    //    cells[D]->print_relevant_info(ofile, cells, ghosts, box);
    //    cells[E]->print_relevant_info(ofile, cells, ghosts, box);
    //    cells[F]->print_relevant_info(ofile, cells, ghosts, box);
    //    ofile.close();
    //    exit(1);
    //    ofile.open("face_flip_after.dat");
    //    cells[A]->complete(cells, ghosts, box);
    //    cells[B]->complete(cells, ghosts, box);
    //    cells[C]->complete(cells, ghosts, box);
    //    cells[D]->complete(cells, ghosts, box);
    //    cells[E]->complete(cells, ghosts, box);
    //    cells[F]->complete(cells, ghosts, box);
    //    cells[A]->print_copy(ofile, 0, wallA, box);
    //    cells[B]->print_copy(ofile, 1, wallB, box);
    //    cells[C]->print_copy(ofile, 2, wallC, box);
    //    cells[D]->print_copy(ofile, 3, wallD, box);
    //    cells[E]->print_copy(ofile, 4, wallE, box);
    //    cells[F]->print_copy(ofile, 5, wallF, box);
    //    ofile.close();
}

/**
 * @brief Count the number of edge and face flips in the cells with the given
 * indices and correct for them
 *
 * @param ids Array with integer indices of cells in the AdaptiveVorTess3d cell
 * list
 * @param numid Number of cell indices in ids
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
void AdaptiveVorCell3d::check_cells(int* ids, unsigned int numid,
                                    std::vector<AdaptiveVorCell3d*>& cells,
                                    std::vector<AdaptiveVorCell3d*>& ghosts,
                                    Cuboid& box, bool periodic) {
    //    cerr << "checking" << endl;
    unsigned int flips[5] = {0, 0, 0, 0, 0};
    unsigned int lastflips[3]{0, 0, 0};
    for(unsigned int i = 0; i < numid; i++) {
        unsigned int numflips =
                cells[ids[i]]->get_flips(ids[i], cells, ghosts, box, periodic);
        flips[numflips]++;
        //        unsigned int numflips =
        // cells[ids[i]]->get_flips_minimal(flips,
        //        ids[i], cells, box, periodic);
        if(numflips == 1) {
            lastflips[0] = ids[i];
        }
        if(numflips == 2) {
            lastflips[1] = ids[i];
        }
        if(numflips == 4) {
            lastflips[2] = ids[i];
        }
    }
// we don't check for invalid flips or invalid combinations, since
// all flips are guaranteed to work, independent of the order in which
// they are performed
// this is due to the convergence of the incremental Delaunay construction
//    cerr << "done" << endl;
#ifdef VERBOSE
    cerr << flips[0] << "\t" << flips[1] << "\t" << flips[2] << "\t" << flips[3]
         << "\t" << flips[4] << endl;
#endif
    unsigned loopcount = 0;
    while(flips[1] || flips[2] || flips[3] || flips[4]) {
        if(flips[2]) {
#ifdef VERBOSE
            cerr << "removal" << endl;
#endif
            cells[lastflips[1]]->remove_face(flips, lastflips[1], cells, ghosts,
                                             box, periodic);
        } else {
            if(flips[4]) {
//                if(!flips[1]){
//                    cerr << "This does not compute. Error!" << endl;
//                    cerr << flips[1] << "\t" << flips[2] << "\t" << flips[3]
//                << "\t" << flips[4] << endl;
//                    ofstream itfile("iterations_error.dat");
//                    for(unsigned int i = 0; i < numid; i++){
//                        cells[ids[i]]->complete(cells, ghosts, box);
//                        cells[ids[i]]->print(itfile, i);
//                    }
//                    itfile.close();
//                    exit(1);
//                }
#ifdef VERBOSE
                cerr << "4 to 4 flip" << endl;
#endif
                cells[lastflips[2]]->flip_face(lastflips[2], cells, ghosts, box,
                                               periodic);
            } else {
                if(flips[1]) {
#ifdef VERBOSE
                    cerr << "insertion" << endl;
#endif
                    cells[lastflips[0]]->insert_face(flips, lastflips[0], cells,
                                                     ghosts, box, periodic);
                }
            }
        }
        flips[0] = 0;
        flips[1] = 0;
        flips[2] = 0;
        flips[3] = 0;
        flips[4] = 0;
        for(unsigned int i = 0; i < numid; i++) {
            unsigned int numflips = cells[ids[i]]->get_flips(
                    ids[i], cells, ghosts, box, periodic);
            flips[numflips]++;
            //            unsigned int numflips =
            // cells[ids[i]]->get_flips_minimal(flips,
            //            ids[i], cells, box, periodic);
            if(numflips == 1) {
                lastflips[0] = ids[i];
            }
            if(numflips == 2) {
                lastflips[1] = ids[i];
            }
            if(numflips == 4) {
                lastflips[2] = ids[i];
            }
        }
        loopcount++;
        if(loopcount > 100) {
            throw AdaptiveMeshException();
            cerr << "Too many iterations. We better stop." << endl;
            ofstream itfile("iterations_error.dat");
            for(unsigned int i = 0; i < numid; i++) {
                cells[ids[i]]->complete(cells, ghosts, box);
                cells[ids[i]]->print(itfile, i);
            }
            itfile.close();
            my_exit();
        }
    }
#ifdef VERBOSE
    cerr << flips[0] << "\t" << flips[1] << "\t" << flips[2] << "\t" << flips[3]
         << "\t" << flips[4] << endl;
    cerr << "nothing" << endl;
#endif
}

/**
 * @brief Count the number of edge and face flips in the cells with the given
 * indices and the given affected edges in these cells, and correct for them
 *
 * @param ids Array with integer indices of cells in the AdaptiveVorTess3d cell
 * list
 * @param numid Number of cell indices in ids
 * @param edges Array with unsigned integer indices of edges that are affected
 * and should be checked
 * @param numedge Number of indices in edges
 * @param cells Reference to the AdaptiveVorTess3d cell list
 * @param ghosts Reference to the AdaptiveVorTess3d ghost cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
void AdaptiveVorCell3d::check_cells(int* ids, unsigned int numid,
                                    unsigned int* edges, unsigned int numedge,
                                    std::vector<AdaptiveVorCell3d*>& cells,
                                    std::vector<AdaptiveVorCell3d*>& ghosts,
                                    Cuboid& box, bool periodic) {
    //    cerr << "checking" << endl;
    unsigned int flips[5];
    flips[0] = 0;
    flips[1] = 0;
    flips[2] = 0;
    flips[3] = 0;
    flips[4] = 0;
    unsigned int lastflips[3] = {0, 0, 0};
    unsigned int curcell = 0;
    unsigned int curface = 0;
    unsigned int numflip = 0;
    unsigned int flip1[2] = {0, 0};
    unsigned int realflip1[2] = {0, 0};
    unsigned int flip3 = 0;
    unsigned int flip4[2] = {0, 0};
    for(unsigned int i = 0; i < numedge; i += 3) {
        if(curcell != edges[i] || curface != edges[i + 1]) {
            //            cerr << curcell << ", " << curface << ": " << numflip
            // << endl;
            if(numflip) {
                if(numflip == 1) {
                    flips[1]++;
                    realflip1[0] = flip1[0];
                    realflip1[1] = flip1[1];
                    lastflips[0] = curcell;
                } else {
                    if(numflip == 2) {
                        unsigned int facesize =
                                cells[ids[curcell]]->_facengbs[curface].size();
                        if(facesize == 4) {
                            flips[4]++;
                            lastflips[2] = curcell;
                        } else {
                            if(facesize == 3) {
                                flips[2]++;
                                flip3 = curface;
                                lastflips[1] = curcell;
                            } else {
                                flips[3]++;
                            }
                        }
                    }
                }
            } else {
                flips[0]++;
            }
            numflip = 0;
        }
        curcell = edges[i];
        curface = edges[i + 1];
        AdaptiveVorCell3d* cell = cells[ids[curcell]];
        unsigned int facesize = cell->_facengbs[curface].size();
        unsigned int j = edges[i + 2];
        unsigned int ni = (j + 1) % facesize;
        unsigned int pi = (j + 2) % facesize;
        int cA = ids[curcell];
        int cB = cell->_ngbs[curface];
        int cC = cell->_facengbs[curface][ni];
        int cD = cell->_facengbs[curface][j];
        int cE = cell->_facengbs[curface][pi];
        int cwA = 0;
        int cwB = cell->_ngbwalls[curface];
        int cwC = cell->_facengbwalls[curface][ni];
        int cwD = cell->_facengbwalls[curface][j];
        int cwE = cell->_facengbwalls[curface][pi];
        double testdouble;
        bool test = in_sphere_test(cA, cB, cC, cD, cE, cwA, cwB, cwC, cwD, cwE,
                                   cells, ghosts, box, periodic, false,
                                   &testdouble);
        //        cerr << "Edge: " << cC << "," << cD << " in face " << cA <<
        // "," << cB
        //        << " (" << cE << "): " << testdouble << endl;
        if(test) {
            //            cerr << "yes" << endl;
            numflip++;
            flip1[0] = curface;
            flip1[1] = j;
            if(numflip == 2) {
                flip4[0] = curface;
                flip4[1] = j;
            }
        } else {
            //            cerr << "no" << endl;
        }
    }
    //    cerr << curcell << ", " << curface << ": " << numflip << endl;
    if(numflip) {
        if(numflip == 1) {
            flips[1]++;
            realflip1[0] = flip1[0];
            realflip1[1] = flip1[1];
            lastflips[0] = curcell;
        } else {
            if(numflip == 2) {
                unsigned int facesize =
                        cells[ids[curcell]]->_facengbs[curface].size();
                if(facesize == 4) {
                    flips[4]++;
                    lastflips[2] = curcell;
                } else {
                    if(facesize == 3) {
                        flips[2]++;
                        flip3 = curface;
                        lastflips[1] = curcell;
                    } else {
                        flips[3]++;
                    }
                }
            }
        }
    } else {
        flips[0]++;
    }
//#warning PROBLEM: if the last face is invalid, then flip1 will contain wrong
//    values...
//    if(flips[1] || flips[2] || flips[3] || flips[4]){
//        cerr << flips[0] << "\t" << flips[1] << "\t" << flips[2] << "\t"
//    << flips[3] << "\t" << flips[4] << endl;
//        ofstream ofile("double_insertion.dat");
//        for(unsigned int i = 0; i < numid; i++){
//            cells[ids[i]]->complete(cells, ghosts, box);
//            cells[ids[i]]->print(ofile, i);
//        }
//        ofile.close();
//        cerr << lastflips[0] << endl;
//        cerr << flip1[0] << "\t" << flip1[1] << endl;
//        exit(1);
//    }
// we don't check for invalid flips or invalid combinations, since
// all flips are guaranteed to work, independent of the order in which
// they are performed
// this is due to the convergence of the incremental Delaunay construction
//    cerr << "done" << endl;
#ifdef VERBOSE
    cerr << flips[0] << "\t" << flips[1] << "\t" << flips[2] << "\t" << flips[3]
         << "\t" << flips[4] << endl;
#endif
    if(flips[4]) {
#ifdef VERBOSE
        cerr << "4 to 4 flip" << endl;
#endif
        cells[ids[lastflips[2]]]->flip_face(ids[lastflips[2]], flip4[0],
                                            flip4[1], cells, ghosts, box,
                                            periodic);
        return;
    }
    if(flips[2]) {
#ifdef VERBOSE
        cerr << "removal" << endl;
#endif
        cells[ids[lastflips[1]]]->remove_face(flips, ids[lastflips[1]], flip3,
                                              cells, ghosts, box, periodic);
        return;
    }
    if(flips[1]) {
#ifdef VERBOSE
        cerr << "insertion (" << flips[1] << ", " << flips[3] << ")" << endl;
#endif
        //        cerr << ids[lastflips[0]] << ", " << flip1[0] << ", " <<
        // flip1[1]
        //        << endl;
        //        cerr << cells[ids[lastflips[0]]]->_ngbs[flip1[0]] << "\t"
        //        << cells[ids[lastflips[0]]]->_facengbs[flip1[0]][flip1[1]] <<
        // endl;
        //        cerr << lastflips[0] << "\t" << realflip1[0] << "\t" <<
        // realflip1[1]
        //        << endl;
        cells[ids[lastflips[0]]]->insert_face(flips, ids[lastflips[0]],
                                              realflip1[0], realflip1[1], cells,
                                              ghosts, box, periodic);
        return;
    }
#ifdef VERBOSE
    cerr << "nothing" << endl;
#endif
}

/**
 * @brief Switch the internal neighbour indices with the new values from the
 * AdaptiveVorTess3d list
 *
 * @param new_ids std::vector with new unsigned integer indices
 */
void AdaptiveVorCell3d::update_ngbs(vector<unsigned int>& new_ids) {
    _id = new_ids[_id];
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        if(_ngbs[i] < (int)new_ids.size()) {
            _ngbs[i] = new_ids[_ngbs[i]];
        }
    }
    for(unsigned int i = 0; i < _facengbs.size(); i++) {
        for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
            if(_facengbs[i][j] < (int)new_ids.size()) {
                _facengbs[i][j] = new_ids[_facengbs[i][j]];
            }
        }
    }
}

/**
 * @brief Switch the internal ghost indices with new values from the
 * AdaptiveVorTess3d list
 *
 * @param new_ids std::vector with new unsigned integer indices
 * @param cellsize New size of the AdaptiveVorTess3d cell list
 */
void AdaptiveVorCell3d::update_ghosts(vector<unsigned int>& new_ids,
                                      unsigned int cellsize) {
    _id = new_ids[_id];
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        if(_ngbs[i] >= (int)cellsize) {
            _ngbs[i] = new_ids[_ngbs[i] - cellsize];
        }
    }
    for(unsigned int i = 0; i < _facengbs.size(); i++) {
        for(unsigned int j = 0; j < _facengbs[i].size(); j++) {
            if(_facengbs[i][j] >= (int)cellsize) {
                _facengbs[i][j] = new_ids[_facengbs[i][j] - cellsize];
            }
        }
    }
}

/**
 * @brief Get the index of this cell in the AdaptiveVorTess3d cell list
 *
 * @return Unsigned integer index of this cell
 */
unsigned int AdaptiveVorCell3d::get_id() {
    return _id;
}

/**
 * @brief Check if all facengbs are also neighbours of this cell
 *
 * If not, the code will abort.
 */
void AdaptiveVorCell3d::clean_facengbs() {
    for(unsigned int i = 0; i < _facengbs.size(); i++) {
        for(unsigned int j = _facengbs[i].size(); j--;) {
            if(!is_ngb(_facengbs[i][j])) {
                throw AdaptiveMeshException();
                cerr << "Error" << endl;
                my_exit();
            }
        }
    }
}
