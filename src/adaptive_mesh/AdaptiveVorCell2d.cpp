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
 * @file AdaptiveVorCell2d.cpp
 *
 * @brief 2D mesh evolution Voronoi cell: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdaptiveVorCell2d.hpp"
#include "AdaptiveCellList.hpp"       // for AdaptiveCellList, etc
#include "AdaptiveMeshException.hpp"  // for AdaptiveMeshException
#include "AdaptiveMeshUtils.hpp"      // for get_wallpos, etc
#include "AdaptiveVorFace2d.hpp"      // for AdaptiveVorFace2d
#include "Error.hpp"                  // for my_exit
#include "MPIGlobal.hpp"              // for rank
#include "MPIMethods.hpp"             // for MyMPI_Pack, MyMPI_Unpack
#include "RestartFile.hpp"            // for RestartFile
#include "StateVector.hpp"            // for StateVector, operator*, etc
#include "predicates.hpp"             // for incircle_old
#include "utilities/GasParticle.hpp"  // for GasParticle
#include <algorithm>                  // for max, min
#include <cmath>                      // for fabs, sqrt, M_PI
#include <iostream>                   // for cerr
#include <iterator>                   // for reverse_iterator, etc
#include <set>                        // for set, etc
using namespace std;

/// @cond this statement excludes the block from doxygen documentation
#ifdef DITWOORD
int VorCell::get_common_ngb(unsigned int idnot, VorCell* a, VorCell* b) {
    for(unsigned int i = 0; i < a->_ngbs.size(); i++) {
        for(unsigned int j = 0; j < b->_ngbs.size(); j++) {
            if(a->_ngbs[i] == b->_ngbs[j] && a->_ngbs[i] != idnot) {
                return a->_ngbs[i];
            }
        }
    }
    return -1;
}
#endif
/// @endcond

/**
 * @brief Constructor.
 *
 * Construct a new cell with the given generator and number of neighbours.
 * Optionally, the original ID and a pointer to the GasParticle associated with
 * the cell can also be set.
 *
 * @param pos Position of the generator of the cell
 * @param nngb Number of neighbours of the cell
 * @param oid ID of the particle of which this cell is a ghost
 * @param particle GasParticle associated with the cell
 */
AdaptiveVorCell2d::AdaptiveVorCell2d(Vec pos, unsigned int nngb,
                                     unsigned int oid, GasParticle* particle) {
    _pos[0] = pos[0];
    _pos[1] = pos[1];
    _ngbs.reserve(nngb);
    _orientation.resize(nngb, 0.);
    _exngbs.resize(nngb, -1);
    _faces.resize(nngb, NULL);
    _oid = oid;
    _wall = false;
    _particle = particle;
    _rank = -1;

    if(particle) {
        _particleID = particle->id();
    }
    _newngbs.reserve(nngb);
}

/**
 * @brief Constructor.
 *
 * Construct a new cell with the given generator and number of neighbours.
 * Optionally, the original ID and a pointer to the GasParticle associated with
 * the cell can also be set.
 *
 * @param pos Position of the generator of the cell
 * @param nngb Number of neighbours of the cell
 * @param oid ID of the particle of which this cell is a ghost
 * @param particle GasParticle associated with the cell
 */
AdaptiveVorCell2d::AdaptiveVorCell2d(double* pos, unsigned int nngb,
                                     unsigned int oid, GasParticle* particle) {
    _pos[0] = pos[0];
    _pos[1] = pos[1];
    _ngbs.reserve(nngb);
    _orientation.resize(nngb, 0.);
    _exngbs.resize(nngb, -1);
    _faces.resize(nngb, NULL);
    _oid = oid;
    _wall = false;
    _particle = particle;
    _rank = -1;

    if(particle) {
        _particleID = particle->id();
    }
    _newngbs.reserve(nngb);
}
/**
 * @brief Destructor.
 *
 * Frees the memory used by the faces of the cell.
 */
AdaptiveVorCell2d::~AdaptiveVorCell2d() {
    for(unsigned int i = 0; i < _faces.size(); i++) {
        delete _faces[i];
    }
}

/**
 * @brief Get the ID of the original GasParticle of which this cell is a ghost
 *
 * @return Unsigned integer ID of a GasParticle
 */
unsigned int AdaptiveVorCell2d::get_original() {
    return _oid;
}

/**
 * @brief Get the neighbours of the cell
 *
 * @return std::vector of integer indices of neighbouring cells in the
 * AdaptiveVorTess
 */
vector<int> AdaptiveVorCell2d::get_ngbs() {
    return _ngbs;
}

/**
 * @brief Get boundary flags for the neighbours of the cell
 *
 * Flag 0 signals a normal neighbour. Flags -1 to -8 signal periodic or
 * reflective copies of cells. The precise conversion between flag values and
 * geometrical boundaries is done in AdaptiveMeshUtils.cpp.
 *
 * @return std::vector of integer flags for every neighbour of the cell
 */
vector<int> AdaptiveVorCell2d::get_ngbwalls() {
    return _ngbwalls;
}

/**
 * @brief Change the neighbour with the given index in the neighbour list to the
 * new given value
 *
 * No sanity check on the index value is performed, this method might cause
 * segmentation faults.
 *
 * @param index Index of the neighbour that has to be replaced in the internal
 * neighbour list
 * @param new_ngb Integer index of the new neighbour in the AdaptiveVorTess
 */
void AdaptiveVorCell2d::change_ngb(unsigned int index, int new_ngb) {
    _ngbs[index] = new_ngb;
}

/// @cond
#ifdef DITWOORD
void VorCell::change_ngb_from_value(int value, int new_ngb) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != value) {
        i++;
    }
    if(i == _ngbs.size()) {
        return;
        //        throw ngbexception(2);
    }
    _ngbs[i] = new_ngb;
}
#endif
/// @endcond

/**
 * @brief Change the neighbour with the given value to the new given value
 *
 * If the given neighbour is not found, this function simply returns and nothing
 * happens, hence its name.
 *
 * @param value Neighbour value that should be replaced
 * @param new_ngb New neighbour value
 */
void AdaptiveVorCell2d::safely_change_ngb_from_value(int value, int new_ngb) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != value) {
        i++;
    }
    if(i == _ngbs.size()) {
        return;
    }
    _ngbs[i] = new_ngb;
}

/**
 * @brief Add the given neighbour to the end of the list
 *
 * @param id Integer index of the neighbour in the AdaptiveVorTess
 */
void AdaptiveVorCell2d::add_ngb(int id) {
    _ngbs.push_back(id);
}

/**
 * @brief Add the given neighbour boundary flag to the end of the list
 *
 * @param id Integer boundary flag
 */
void AdaptiveVorCell2d::add_ngbwall(int id) {
    _ngbwalls.push_back(id);
}

/// @cond
#ifdef DITWOORD
void VorCell::set_orientations() {
    for(unsigned int i = 0; i < _faces.size(); i++) {
        double a[2], b[2];
        _faces[i]->get_points(a, b);
        _orientation[i] = predicates::orient2d(a, b, _pos);
    }
}
#endif
/// @endcond

// void AdaptiveVorCell2d::complete(vector<AdaptiveVorCell2d*>& cells,
// vector<AdaptiveVorCell2d*>& ghosts, std::vector<AdaptiveVorCell2d*>& orphans,
// Cuboid& box, bool periodic){
//    for(unsigned int i = 0; i < _faces.size(); i++){
//        delete _faces[i];
//    }
//    _wall = false;
//    for(unsigned int i = 0; i < _ngbs.size(); i++){
//        double ngbpos[2];
//        if(_ngbs[i] < 0 || _ngbs[i] < cells.size()){
//            if(_ngbwalls[i] > 0){
//                AdaptiveVorCell2d* ngb = orphans[-_ngbs[i]];
//                ngbpos[0] = ngb->_pos[0];
//                ngbpos[1] = ngb->_pos[1];
//                if(periodic && _ngbwalls[i]-9 < 0){
//                    AdaptiveMeshUtils::get_periodic_position(ngbpos,
//_ngbwalls[i]-9, box);
//                }
//            } else {
//                AdaptiveVorCell2d* ngb = cells[_ngbs[i]];
//                ngbpos[0] = ngb->_pos[0];
//                ngbpos[1] = ngb->_pos[1];
//                if(periodic && _ngbwalls[i] < 0){
//                    AdaptiveMeshUtils::get_periodic_position(ngbpos,
//_ngbwalls[i], box);
//                }
//            }
//        } else {
//            _wall = true;
//            AdaptiveVorCell2d* ghost = ghosts[_ngbs[i]-cells.size()];
//            AdaptiveMeshUtils::get_wall_position(ghost->_pos, ngbpos,
// ghost->_ngbs[0], box);
//        }
//        double mid[2];
//        mid[0] = 0.5*(_pos[0]+ngbpos[0]);
//        mid[1] = 0.5*(_pos[1]+ngbpos[1]);
//        _faces[i] = new AdaptiveVorFace2d(mid);
//        if(i){
//            double prevpos[2];
//            if(_ngbs[i-1] < 0 || _ngbs[i-1] < cells.size()){
//                if(_ngbwalls[i-1] > 0){
//                    AdaptiveVorCell2d* prev = orphans[-_ngbs[i-1]];
//                    prevpos[0] = prev->_pos[0];
//                    prevpos[1] = prev->_pos[1];
//                    if(periodic && _ngbwalls[i-1]-9 < 0){
//                        AdaptiveMeshUtils::get_periodic_position(prevpos,
//_ngbwalls[i-1]-9, box);
//                    }
//                } else {
//                    AdaptiveVorCell2d* prev;
//                    if(_faces[i-1]->active()){
//                        prev = cells[_ngbs[i-1]];
//                    } else {
//                        prev = cells[_ngbs[i-2]];
//                    }
//                    prevpos[0] = prev->_pos[0];
//                    prevpos[1] = prev->_pos[1];
//                    if(periodic && _ngbwalls[i-1] < 0){
//                        AdaptiveMeshUtils::get_periodic_position(prevpos,
//_ngbwalls[i-1], box);
//                    }
//                }
//            } else {
//                AdaptiveVorCell2d* ghost = ghosts[_ngbs[i-1]-cells.size()];
//                AdaptiveMeshUtils::get_wall_position(ghost->_pos, prevpos,
// ghost->_ngbs[0], box);
//            }
//            double vert[2];
//            if(prevpos[0] == ngbpos[0] && prevpos[1] == ngbpos[1]){
//                cerr << "ERROR" << endl;
//                exit(1);
//            }
//            _faces[i]->add_Lvertex(prevpos, ngbpos, _pos, vert);
//            _faces[i-1]->set_Rvertex(vert);
//        }
//    }
//    double ngbpos[2], prevpos[2];
//    if(_ngbs[0] < 0 || _ngbs[0] < cells.size()){
//        if(_ngbwalls[0] > 0){
//            AdaptiveVorCell2d* ngb = orphans[-_ngbs[0]];
//            ngbpos[0] = ngb->_pos[0];
//            ngbpos[1] = ngb->_pos[1];
//            if(periodic && _ngbwalls[0] < 9){
//                AdaptiveMeshUtils::get_periodic_position(ngbpos,
//_ngbwalls[0]-9, box);
//            }
//        } else {
//            AdaptiveVorCell2d* ngb = cells[_ngbs[0]];
//            ngbpos[0] = ngb->_pos[0];
//            ngbpos[1] = ngb->_pos[1];
//            if(periodic && _ngbwalls[0] < 0){
//                AdaptiveMeshUtils::get_periodic_position(ngbpos, _ngbwalls[0],
// box);
//            }
//        }
//    } else {
//        AdaptiveVorCell2d* ghost = ghosts[_ngbs[0]-cells.size()];
//        AdaptiveMeshUtils::get_wall_position(ghost->_pos, ngbpos,
// ghost->_ngbs[0], box);
//    }
//    if(_ngbs.back() < 0 || _ngbs.back() < cells.size()){
//        if(_ngbwalls.back() > 0){
//            AdaptiveVorCell2d* prev = orphans[-_ngbs.back()];
//            prevpos[0] = prev->_pos[0];
//            prevpos[1] = prev->_pos[1];
//            if(periodic && _ngbwalls.back() < 9){
//                AdaptiveMeshUtils::get_periodic_position(prevpos,
//_ngbwalls.back()-9, box);
//            }
//        } else {
//            AdaptiveVorCell2d* prev;
//            if(_faces.back()->active()){
//                prev = cells[_ngbs.back()];
//            } else {
//                prev = cells[_ngbs[_ngbs.size()-2]];
//            }
//            prevpos[0] = prev->_pos[0];
//            prevpos[1] = prev->_pos[1];
//            if(periodic && _ngbwalls.back() < 0){
//                AdaptiveMeshUtils::get_periodic_position(prevpos,
//_ngbwalls.back(), box);
//            }
//        }
//    } else {
//        AdaptiveVorCell2d* ghost = ghosts[_ngbs.back()-cells.size()];
//        AdaptiveMeshUtils::get_wall_position(ghost->_pos, prevpos,
// ghost->_ngbs[0], box);
//    }
//    if(_faces[0]->active()){
//        double vert[2];
//        if(prevpos[0] == ngbpos[0] && prevpos[1] == ngbpos[1]){
//            cerr << "ERROR" << endl;
//            cerr << _ngbs[0] << "\t" << _ngbs.back() << endl;
//            exit(1);
//        }
//        _faces[0]->add_Lvertex(prevpos, ngbpos, _pos, vert);
//        _faces.back()->set_Rvertex(vert);
//    }

//    for(unsigned int i = 0; i < _faces.size(); i++){
//        if(!_faces[i]->active()){
//            double L[2], R[2];
//            _faces[(i+1)%_faces.size()]->get_points(L, R);
//            _faces[(i+_faces.size()-1)%_faces.size()]->set_Rvertex(L);
//            _faces[i]->set_points(L, L);
//        }
//    }
//}

/**
 * @brief Calculate the faces of this cell to make it geometrically complete
 *
 * The list of faces is erased and new faces are created by looping over the
 * list of neighbours. For every neighbour, the vertices of the corresponding
 * face are the midpoints of the circumcircles through respectively the cell,
 * the neighbour and the previous neighbour in the list, and the cell, the
 * neighbour and the next neighbour in the list. If a neighbour is a ghost, the
 * boundary flag for that neighbour is used to convert its coordinates to usable
 * coordinates.
 *
 * @param cells List of all the cells in the AdaptiveVorTess, used to retrieve
 * neighbour information
 * @param box Cuboid specifying the dimensions of the simulation box, used to
 * convert ghost positions
 * @param periodic Flag signaling if a periodic (true) or reflective (false)
 * simulation box is used.
 */
void AdaptiveVorCell2d::complete(AdaptiveCellList<AdaptiveVorCell2d>& cells,
                                 Cuboid& box, bool periodic) {
    for(unsigned int i = 0; i < _faces.size(); i++) {
        delete _faces[i];
    }
    _wall = false;
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        double ngbpos[2];
        if(!cells.is_ghost(_ngbs[i])) {
            if(cells.is_orphan(_ngbs[i])) {
                AdaptiveVorCell2d* ngb = cells.get_orphan(_ngbs[i]);
                ngbpos[0] = ngb->_pos[0];
                ngbpos[1] = ngb->_pos[1];
                if(periodic && _ngbwalls[i] < 0) {
                    AdaptiveMeshUtils::get_periodic_position(ngbpos,
                                                             _ngbwalls[i], box);
                }
            } else {
                AdaptiveVorCell2d* ngb = cells.get_normal(_ngbs[i]);
                ngbpos[0] = ngb->_pos[0];
                ngbpos[1] = ngb->_pos[1];
                if(periodic && _ngbwalls[i] < 0) {
                    AdaptiveMeshUtils::get_periodic_position(ngbpos,
                                                             _ngbwalls[i], box);
                }
            }
        } else {
            _wall = true;
            AdaptiveVorCell2d* ghost = cells.get_ghost(_ngbs[i]);
            AdaptiveMeshUtils::get_wall_position(ghost->_pos, ngbpos,
                                                 ghost->_ngbs[0], box);
        }
        double mid[2];
        mid[0] = 0.5 * (_pos[0] + ngbpos[0]);
        mid[1] = 0.5 * (_pos[1] + ngbpos[1]);
        _faces[i] = new AdaptiveVorFace2d(mid);
        if(i) {
            double prevpos[2];
            if(!cells.is_ghost(_ngbs[i - 1])) {
                if(cells.is_orphan(_ngbs[i - 1])) {
                    AdaptiveVorCell2d* prev = cells.get_orphan(_ngbs[i - 1]);
                    prevpos[0] = prev->_pos[0];
                    prevpos[1] = prev->_pos[1];
                    if(periodic && _ngbwalls[i - 1] < 0) {
                        AdaptiveMeshUtils::get_periodic_position(
                                prevpos, _ngbwalls[i - 1], box);
                    }
                } else {
                    AdaptiveVorCell2d* prev = cells.get_normal(_ngbs[i - 1]);
                    prevpos[0] = prev->_pos[0];
                    prevpos[1] = prev->_pos[1];
                    if(periodic && _ngbwalls[i - 1] < 0) {
                        AdaptiveMeshUtils::get_periodic_position(
                                prevpos, _ngbwalls[i - 1], box);
                    }
                }
            } else {
                AdaptiveVorCell2d* ghost = cells.get_ghost(_ngbs[i - 1]);
                AdaptiveMeshUtils::get_wall_position(ghost->_pos, prevpos,
                                                     ghost->_ngbs[0], box);
            }
            double vert[2];
            if(prevpos[0] == ngbpos[0] && prevpos[1] == ngbpos[1]) {
                cerr << "ERROR" << endl;
                my_exit();
            }
            _faces[i]->add_Lvertex(prevpos, ngbpos, _pos, vert);
            _faces[i - 1]->set_Rvertex(vert);
        }
    }
    double ngbpos[2], prevpos[2];
    if(!cells.is_ghost(_ngbs[0])) {
        if(cells.is_orphan(_ngbs[0])) {
            AdaptiveVorCell2d* ngb = cells.get_orphan(_ngbs[0]);
            ngbpos[0] = ngb->_pos[0];
            ngbpos[1] = ngb->_pos[1];
            if(periodic && _ngbwalls[0] < 0) {
                AdaptiveMeshUtils::get_periodic_position(ngbpos, _ngbwalls[0],
                                                         box);
            }
        } else {
            AdaptiveVorCell2d* ngb = cells.get_normal(_ngbs[0]);
            ngbpos[0] = ngb->_pos[0];
            ngbpos[1] = ngb->_pos[1];
            if(periodic && _ngbwalls[0] < 0) {
                AdaptiveMeshUtils::get_periodic_position(ngbpos, _ngbwalls[0],
                                                         box);
            }
        }
    } else {
        AdaptiveVorCell2d* ghost = cells.get_ghost(_ngbs[0]);
        AdaptiveMeshUtils::get_wall_position(ghost->_pos, ngbpos,
                                             ghost->_ngbs[0], box);
    }
    if(!cells.is_ghost(_ngbs.back())) {
        if(cells.is_orphan(_ngbs.back())) {
            AdaptiveVorCell2d* prev = cells.get_orphan(_ngbs.back());
            prevpos[0] = prev->_pos[0];
            prevpos[1] = prev->_pos[1];
            if(periodic && _ngbwalls.back() < 0) {
                AdaptiveMeshUtils::get_periodic_position(prevpos,
                                                         _ngbwalls.back(), box);
            }
        } else {
            AdaptiveVorCell2d* prev = cells.get_normal(_ngbs.back());
            prevpos[0] = prev->_pos[0];
            prevpos[1] = prev->_pos[1];
            if(periodic && _ngbwalls.back() < 0) {
                AdaptiveMeshUtils::get_periodic_position(prevpos,
                                                         _ngbwalls.back(), box);
            }
        }
    } else {
        AdaptiveVorCell2d* ghost = cells.get_ghost(_ngbs.back());
        AdaptiveMeshUtils::get_wall_position(ghost->_pos, prevpos,
                                             ghost->_ngbs[0], box);
    }

    double vert[2];
    if(prevpos[0] == ngbpos[0] && prevpos[1] == ngbpos[1]) {
        cerr << "ERROR" << endl;
        cerr << _ngbs[0] << "\t" << _ngbs.back() << endl;
        my_exit();
    }
    _faces[0]->add_Lvertex(prevpos, ngbpos, _pos, vert);
    _faces.back()->set_Rvertex(vert);
}

/**
 * @brief Calculate the faces of the cell
 *
 * @param list NewCellList holding accurate positions for all IDs
 */
void AdaptiveVorCell2d::complete(NewCellList& list) {
    for(unsigned int i = 0; i < _faces.size(); i++) {
        delete _faces[i];
    }
    Vec pos = list.get_position(_particleID);
    _pos[0] = pos.x();
    _pos[1] = pos.y();
    _wall = false;
    for(unsigned int i = 0; i < _newngbs.size(); i++) {
        double ngbpos[2];
        pos = list.get_position(_newngbs[i]);
        ngbpos[0] = pos.x();
        ngbpos[1] = pos.y();
        double mid[2];
        mid[0] = 0.5 * (_pos[0] + ngbpos[0]);
        mid[1] = 0.5 * (_pos[1] + ngbpos[1]);
        _faces[i] = new AdaptiveVorFace2d(mid);
        if(i) {
            double prevpos[2];
            pos = list.get_position(_newngbs[i - 1]);
            prevpos[0] = pos.x();
            prevpos[1] = pos.y();
            double vert[2];
            if(prevpos[0] == ngbpos[0] && prevpos[1] == ngbpos[1]) {
                cerr << "ERROR" << endl;
                cerr << prevpos[0] << "\t" << prevpos[1] << endl;
                cerr << ngbpos[0] << "\t" << ngbpos[1] << endl;
                my_exit();
            }
            _faces[i]->add_Lvertex(prevpos, ngbpos, _pos, vert);
            _faces[i - 1]->set_Rvertex(vert);
        }
    }
    double ngbpos[2], prevpos[2];
    pos = list.get_position(_newngbs[0]);
    ngbpos[0] = pos.x();
    ngbpos[1] = pos.y();
    pos = list.get_position(_newngbs.back());
    prevpos[0] = pos.x();
    prevpos[1] = pos.y();

    double vert[2];
    if(prevpos[0] == ngbpos[0] && prevpos[1] == ngbpos[1]) {
        cerr << "ERROR" << endl;
        cerr << _ngbs[0] << "\t" << _ngbs.back() << endl;
        cerr << prevpos[0] << "\t" << prevpos[1] << endl;
        cerr << ngbpos[0] << "\t" << ngbpos[1] << endl;
        my_exit();
    }
    _faces[0]->add_Lvertex(prevpos, ngbpos, _pos, vert);
    _faces.back()->set_Rvertex(vert);
}

/**
 * @brief Debug method: print cell information to the given output stream
 *
 * @param stream std::ostream to print to
 * @param id ID of the cell
 */
void AdaptiveVorCell2d::print(ostream& stream, int id = -1) {
    stream << _pos[0] << "\t" << _pos[1];
    if(id >= 0) {
        stream << "\t" << _particleID << endl;
    }
    stream << "\n\n";
    for(unsigned int i = 0; i < _faces.size(); i++) {
        _faces[i]->print(stream);
    }
    stream << endl;
}

/// @cond
#ifdef DITWOORD
void VorCell::fill_exngbs(unsigned int id, vector<VorCell*>& cells) {
    for(unsigned int i = 0; i < _ngbs.size() - 1; i++) {
        if(_ngbs[i] >= 0 && _ngbs[i + 1] >= 0) {
            _exngbs[i] =
                    get_common_ngb(id, cells[_ngbs[i]], cells[_ngbs[i + 1]]);
        }
    }
    if(_ngbs[0] >= 0 && _ngbs.back() >= 0) {
        _exngbs[_exngbs.size() - 1] =
                get_common_ngb(id, cells[_ngbs.back()], cells[_ngbs[0]]);
    }
}

void VorCell::print_exngbs(ostream& stream, vector<VorCell*>& cells) {
    for(unsigned int i = 0; i < _exngbs.size(); i++) {
        if(_exngbs[i] >= 0) {
            cells[_exngbs[i]]->print(stream);
        }
    }
}
#endif
/// @endcond

/**
 * @brief Calculate the centroid of the cell
 *
 * The centroid is the weighted mean of the centroids of the triangles
 * constituted of the generator of the cell and the two vertices of its faces.
 *
 * @param centroid Array to store the result in
 */
void AdaptiveVorCell2d::get_centroid(double* centroid) {
    centroid[0] = 0.;
    centroid[1] = 0.;
    double A_tot = 0.;
    for(unsigned int i = 0; i < _faces.size(); i++) {
        double L[2], R[2];
        _faces[i]->get_points(L, R);
        double C[2];
        C[0] = (L[0] + R[0] + _pos[0]) / 3.;
        C[1] = (L[1] + R[1] + _pos[1]) / 3.;
        double a[2], b[2];
        a[0] = L[0] - _pos[0];
        a[1] = L[1] - _pos[1];
        b[0] = R[0] - _pos[0];
        b[1] = R[1] - _pos[1];
        double A = 0.5 * fabs(a[0] * b[1] - a[1] * b[0]);
        centroid[0] += A * C[0];
        centroid[1] += A * C[1];
        A_tot += A;
    }
    centroid[0] /= A_tot;
    centroid[1] /= A_tot;
}

/**
 * @brief Move the generator of the cell to the given position
 *
 * @param pos New position for the generator of the cell
 */
void AdaptiveVorCell2d::move(double* pos) {
    _pos[0] = pos[0];
    _pos[1] = pos[1];
}

/// @cond
#ifdef DITWOORD
bool VorCell::check(unsigned int id, vector<VorCell*>& cells,
                    ostream& estream) {
    bool wrong = false;
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        int ngb = _ngbs[i];
        int prev = _ngbs[(i + _ngbs.size() - 1) % _ngbs.size()];
        if(ngb < cells.size() && prev < cells.size()) {
            VorCell* ngbcell = cells[ngb];
            VorCell* prevcell = cells[prev];
            for(unsigned int j = 0; j < cells.size(); j++) {
                if(j != ngb && j != prev && j != id) {
                    if(inside(ngbcell->_pos, prevcell->_pos, _pos,
                              cells[j]->_pos)) {
                        wrong |= true;
                        double vert[2];
                        get_vertexpoint(ngbcell->_pos, prevcell->_pos,
                                        cells[j]->_pos, vert);
                        estream << vert[0] << "\t" << vert[1] << "\n";
                    }
                }
            }
        }
    }
    return wrong;
}
#endif
/// @endcond

/**
 * @brief Retrieve the position of the cell generator
 *
 * @param pos Array to store the result in
 */
void AdaptiveVorCell2d::get_position(double* pos) {
    pos[0] = _pos[0];
    pos[1] = _pos[1];
}

/**
 * @brief Remove the neighbour with the given index from the list
 *
 * If the neighbour is not found, the code will abort.
 *
 * @param id Unsigned integer index of the neighbour that has to be removed
 * @param periodic Flag specifying if the simulation box is periodic (true) or
 * reflective (false).
 */
void AdaptiveVorCell2d::remove_ngb(unsigned int id, bool periodic) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != (int)id) {
        i++;
    }
    if(i == _ngbs.size()) {
        cerr << "Ngbexception(1)" << endl;
        cerr << "oid: " << _oid << endl;
        cerr << "id: " << id << endl;
        cerr << "ngbs:" << endl;
        for(unsigned int j = 0; j < _ngbs.size(); j++) {
            cerr << j << ": " << _ngbs[j] << endl;
        }
        my_exit();
    }
    _ngbs.erase(_ngbs.begin() + i);
    // don't do this for orphans
    if(_faces.size()) {
        _faces.erase(_faces.begin() + i);
        _orientation.erase(_orientation.begin() + i);
        _exngbs.erase(_exngbs.end() - 1);
    }
    if(periodic) {
        _ngbwalls.erase(_ngbwalls.begin() + i);
    }
}

/**
 * @brief Remove the neighbour with the given index from the list
 *
 * If the neighbour is not found, the code will abort.
 *
 * @param id ID of the neighbour that has to be removed
 * @param periodic Flag specifying if the simulation box is periodic (true) or
 * reflective (false).
 */
void AdaptiveVorCell2d::remove_newngb(unsigned long id, bool periodic) {
    unsigned int i = 0;
    while(i < _newngbs.size() && _newngbs[i] != id) {
        i++;
    }
    if(i == _newngbs.size()) {
        cerr << "Newngb to remove not found!" << endl;
        cerr << id << endl;
        cerr << _particleID << endl;
        cerr << "Ngbs:" << endl;
        for(unsigned int j = 0; j < _newngbs.size(); j++) {
            cerr << _newngbs[j] << endl;
        }
        my_exit();
    }
    _newngbs.erase(_newngbs.begin() + i);
    _faces.erase(_faces.begin() + i);
}

/**
 * @brief Remove the neighbour with the given index from the list
 *
 * If the neighbour is not found, nothing will happen, hence the name.
 *
 * @param id Unsigned integer index of the neighbour that has to be removed
 * @param periodic Flag specifying if the simulation box is periodic (true) or
 * reflective (false).
 */
void AdaptiveVorCell2d::safely_remove_ngb(unsigned int id, bool periodic) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != (int)id) {
        i++;
    }
    if(i == _ngbs.size()) {
        // safely = do not crash on this one
        return;
    }
    _ngbs.erase(_ngbs.begin() + i);
    _faces.erase(_faces.begin() + i);
    _orientation.erase(_orientation.begin() + i);
    _exngbs.erase(_exngbs.end() - 1);
    if(periodic) {
        _ngbwalls.erase(_ngbwalls.begin() + i);
    }
}

/**
 * @brief Add the given neighbour to the list, in between the neighbours with
 * the given indices
 *
 * If none of the given neighbour indices is found in the list, or only one is
 * found, the neighbour is added to the end of the list. If the given neighbour
 * is already a neighbour of this cell, the code will abort.
 *
 * @param id1 Index of an existing neighbour of the cell
 * @param id2 Index of another existing neighbour of the cell
 * @param ngb Index of the new neighbour that should be added to the cell
 */
void AdaptiveVorCell2d::add_ngb(unsigned int id1, unsigned int id2, int ngb) {
    if(is_ngb(ngb)) {
        cerr << "Already a ngb(1)!" << endl;
        my_exit();
    }
    _exngbs.push_back(-1);
    unsigned int idx1 = 0, idx2 = 0;
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        if(_ngbs[i] == (int)id1) {
            idx1 = i;
        }
        if(_ngbs[i] == (int)id2) {
            idx2 = i;
        }
    }
    unsigned int diff = max(idx1, idx2) - min(idx1, idx2);
    if(diff == 1) {
        _ngbs.insert(_ngbs.begin() + max(idx1, idx2), ngb);
        _orientation.insert(_orientation.begin() + max(idx1, idx2), 0.);
        _faces.insert(_faces.begin() + max(idx1, idx2), NULL);
    } else {
        _ngbs.push_back(ngb);
        _orientation.push_back(0.);
        _faces.push_back(NULL);
    }
}

/**
 * @brief Add the given neighbour to the list, in between the neighbours with
 * the given IDs
 *
 * If none of the given neighbour indices is found in the list, or only one is
 * found, the neighbour is added to the end of the list. If the given neighbour
 * is already a neighbour of this cell, the code will abort.
 *
 * @param id1 ID of an existing neighbour of the cell
 * @param id2 ID of another existing neighbour of the cell
 * @param ngb ID of the new neighbour that should be added to the cell
 */
void AdaptiveVorCell2d::add_newngb(unsigned long id1, unsigned long id2,
                                   unsigned long ngb) {
    unsigned int idx1 = 0, idx2 = 0;
    for(unsigned int i = 0; i < _newngbs.size(); i++) {
        if(_newngbs[i] == id1) {
            idx1 = i;
        }
        if(_newngbs[i] == id2) {
            idx2 = i;
        }
    }
    unsigned int diff = max(idx1, idx2) - min(idx1, idx2);
    if(diff == 1) {
        _newngbs.insert(_newngbs.begin() + max(idx1, idx2), ngb);
        _faces.insert(_faces.begin() + max(idx1, idx2), NULL);
    } else {
        _newngbs.push_back(ngb);
        _faces.push_back(NULL);
    }
}

/**
 * @brief Add the given neighbour and associated boundary flag to the list, in
 * between the neighbours with the given indices
 *
 * If none of the given neighbour indices is found in the list, or only one is
 * found, the neighbour is added to the end of the list. If the given neighbour
 * is already a neighbour of this cell, the code will abort.
 *
 * @param id1 Index of an existing neighbour of the cell
 * @param id2 Index of another existing neighbour of the cell
 * @param ngb Index of the new neighbour that should be added to the cell
 * @param wallpos Boundary flag for the new neighbour
 */
void AdaptiveVorCell2d::add_ngb(unsigned int id1, unsigned int id2, int ngb,
                                int wallpos) {
    if(is_ngb(ngb)) {
        cerr << "Already a ngb(2)!" << endl;
        my_exit();
    }
    _exngbs.push_back(-1);
    unsigned int idx1 = 0, idx2 = 0;
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        if(_ngbs[i] == (int)id1) {
            idx1 = i;
        }
        if(_ngbs[i] == (int)id2) {
            idx2 = i;
        }
    }
    unsigned int diff = max(idx1, idx2) - min(idx1, idx2);
    if(diff == 1) {
        _ngbs.insert(_ngbs.begin() + max(idx1, idx2), ngb);
        _orientation.insert(_orientation.begin() + max(idx1, idx2), 0.);
        _faces.insert(_faces.begin() + max(idx1, idx2), NULL);
        _ngbwalls.insert(_ngbwalls.begin() + max(idx1, idx2), wallpos);
    } else {
        _ngbs.push_back(ngb);
        _orientation.push_back(0.);
        _faces.push_back(NULL);
        _ngbwalls.push_back(wallpos);
    }
}

/**
 * @brief Get the two neighbours of a given neighbour in the neighbour list
 *
 * Suppose A, B, and C are three consecutive neighbours in the neighbourlist. If
 * we pass on the ID of B to this method, it will return C and A, in that order.
 *
 * If id is not a neighbour, this method might get stuck in an endless loop or
 * crash with a segmentation fault.
 *
 * @param id ID of the neighbour of which we want the two neighbours
 * @param ngbs Array to store the result in
 */
void AdaptiveVorCell2d::get_ngbs(unsigned int id, int* ngbs) {
    unsigned int i = 0;
    while(_ngbs[i] != (int)id) {
        i++;
    }
    ngbs[0] = _ngbs[(i + 1) % _ngbs.size()];
    ngbs[1] = _ngbs[(i + _ngbs.size() - 1) % _ngbs.size()];
}

/**
 * @brief Add the given new neighbour before the given neighbour in the list
 *
 * If the given new neighbour is already a neighbour or if the given old
 * neighbour is not found in the list, the code will abort.
 *
 * @param id Unsigned integer index of a neighbour that is already in the list
 * @param ngb Integer index of a new neighbour to add to the list
 */
void AdaptiveVorCell2d::add_ngb_before(unsigned int id, int ngb) {
    if(is_ngb(ngb)) {
        cerr << "Already a ngb(3)!" << endl;
        my_exit();
    }
    _exngbs.push_back(-1);
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != (int)id) {
        i++;
    }
    if(i == _ngbs.size()) {
        //        throw ngbexception(3);
        cerr << "Ngbexception(3)" << endl;
        my_exit();
    }
    _ngbs.insert(_ngbs.begin() + i, ngb);
    _orientation.insert(_orientation.begin() + i, 0.);
    _faces.insert(_faces.begin() + i, NULL);
}

/**
 * @brief Add the given new neighbour and the associated boundary flag before
 * the given neighbour in the list
 *
 * If the given new neighbour is already a neighbour or if the given old
 * neighbour is not found in the list, the code will abort.
 *
 * @param id Unsigned integer index of a neighbour that is already in the list
 * @param ngb Integer index of a new neighbour to add to the list
 * @param wallpos Boundary flag for the new neighbour
 */
void AdaptiveVorCell2d::add_ngb_before(unsigned int id, int ngb, int wallpos) {
    if(is_ngb(ngb)) {
        cerr << "Already a ngb(4)!" << endl;
        my_exit();
    }
    _exngbs.push_back(-1);
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != (int)id) {
        i++;
    }
    if(i == _ngbs.size()) {
        //        throw ngbexception(3);
        cerr << "Ngbexception(3)" << endl;
        my_exit();
    }
    _ngbs.insert(_ngbs.begin() + i, ngb);
    _orientation.insert(_orientation.begin() + i, 0.);
    _faces.insert(_faces.begin() + i, NULL);
    _ngbwalls.insert(_ngbwalls.begin() + i, wallpos);
}

/**
 * @brief Add the given new neighbour after the given neighbour in the list
 *
 * If the given new neighbour is already a neighbour or if the given old
 * neighbour is not found in the list, the code will abort.
 *
 * @param id Unsigned integer index of a neighbour that is already in the list
 * @param ngb Integer index of a new neighbour to add to the list
 */
void AdaptiveVorCell2d::add_ngb_after(unsigned int id, int ngb) {
    if(is_ngb(ngb)) {
        cerr << "Already a ngb(5)!" << endl;
        my_exit();
    }
    _exngbs.push_back(-1);
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != (int)id) {
        i++;
    }
    if(i == _ngbs.size()) {
        //        throw ngbexception(3);
        cerr << "Ngbexception(3)" << endl;
        my_exit();
    }
    _ngbs.insert(_ngbs.begin() + i + 1, ngb);
    _orientation.insert(_orientation.begin() + i + 1, 0.);
    _faces.insert(_faces.begin() + i + 1, NULL);
}

/**
 * @brief Add the given new neighbour and the associated boundary flag after the
 * given neighbour in the list
 *
 * If the given new neighbour is already a neighbour or if the given old
 * neighbour is not found in the list, the code will abort.
 *
 * @param id Unsigned integer index of a neighbour that is already in the list
 * @param ngb Integer index of a new neighbour to add to the list
 * @param wallpos Boundary flag for the new neighbour
 */
void AdaptiveVorCell2d::add_ngb_after(unsigned int id, int ngb, int wallpos) {
    if(is_ngb(ngb)) {
        cerr << "Already a ngb(6)!" << endl;
        my_exit();
    }
    _exngbs.push_back(-1);
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != (int)id) {
        i++;
    }
    if(i == _ngbs.size()) {
        //        throw ngbexception(3);
        cerr << "Ngbexception(3)" << endl;
        my_exit();
    }
    _ngbs.insert(_ngbs.begin() + i + 1, ngb);
    _orientation.insert(_orientation.begin() + i + 1, 0.);
    _faces.insert(_faces.begin() + i + 1, NULL);
    _ngbwalls.insert(_ngbwalls.begin() + i + 1, wallpos);
}

/**
 * @brief Test if the generator of the cell is inside the polygon formed by its
 * faces
 *
 * It could happen that this is no longer the case after a movement of the
 * generator. In this case, our evolution algorithm will not work and we have to
 * abort.
 *
 * The test is performed by calculating the so-called winding number of the
 * polygon around the point. If the point is inside the polygon, this winding
 * number will always be zero, irrespective of the exact form of the polygon. If
 * it is outside, the winding number will have a non-zero value.
 *
 * The winding number is ill-defined if the point is part of the polygon, but
 * this cannot happen in this specific case due to the way the polygon is
 * constructed.
 */
void AdaptiveVorCell2d::compute_winding_number() {
    // source: http://geomalgorithms.com/a03-_inclusion.html
    int wn = 0;
    vector<int> wnv(_faces.size(), 0);
    vector<double> isL(_faces.size(), 0.);
    for(unsigned int i = 0; i < _faces.size(); i++) {
        double a[2], b[2];
        _faces[i]->get_points(a, b);
        if(a[1] <= _pos[1]) {
            if(b[1] > _pos[1]) {
                double isLeft = (b[0] - a[0]) * (_pos[1] - a[1]) -
                                (_pos[0] - a[0]) * (b[1] - a[1]);
                isL[i] = isLeft;
                if(isLeft > 0.) {
                    wn++;
                    wnv[i] = 1;
                }
            }
        } else {
            if(b[1] <= _pos[1]) {
                double isLeft = (b[0] - a[0]) * (_pos[1] - a[1]) -
                                (_pos[0] - a[0]) * (b[1] - a[1]);
                isL[i] = isLeft;
                if(isLeft < 0.) {
                    wn--;
                    wnv[i] = -1;
                }
            }
        }
    }
    if(!wn) {
        cerr << "Point outside cell!" << endl;
//#define CRASHONERROR

#ifdef CRASHONERROR
        cerr << _pos[0] << "\t" << _pos[1] << endl;
        for(unsigned int i = 0; i < _faces.size(); i++) {
            double a[2], b[2];
            _faces[i]->get_points(a, b);
            cerr << wnv[i] << " -- " << isL[i] << " [(" << a[0] << "," << a[1]
                 << "), (" << b[0] << "," << b[1] << ")]" << endl;
        }
        ofstream failfile("failcell.dat");
        print(failfile);
        failfile.close();
        my_exit();
#else
        throw AdaptiveMeshException();
#endif
    }
}

/**
 * @brief Kernel of the mesh evolution algorithm: detect errors in the evolved
 * mesh and correct them
 *
 * @param id Index of the cell in the AdaptiveVorTess cell list
 * @param newlist Reference to the AdaptiveVorTess cell list
 * @param insertions Temporary buffer for parallel algorithm
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param periodic Flag discriminating between a periodic (true) and a
 * reflective (false) simulation box
 * @return True if errors where corrected, false otherwise
 */
bool AdaptiveVorCell2d::detect_crossovers(
        unsigned int id, AdaptiveCellList<AdaptiveVorCell2d>& newlist,
        std::vector<int>& insertions, Cuboid& box, bool periodic) {
    if(_faces.size() != _ngbs.size() || _faces.size() != _exngbs.size()) {
        cerr << "Sizes of internal vectors in VorCell are different!" << endl;
        cerr << _faces.size() << "\t" << _ngbs.size() << "\t" << _exngbs.size()
             << endl;
        my_exit();
    }
    if(periodic) {
        if(_ngbs.size() != _ngbwalls.size()) {
            cerr << "Sizes of internal vectors in VorCell are different!"
                 << endl;
            my_exit();
        }
    }
    complete(newlist, box, periodic);
    compute_winding_number();
    set<unsigned int> to_remove;
    vector<unsigned int> to_check;
    for(unsigned int i = 0; i < _faces.size(); i++) {
        if(!_faces[i]) {
            continue;
        }
        unsigned int previ = (i + _faces.size() - 1) % _faces.size();
        // it can happen that the previous or next face was already removed
        // earlier on
        while(!_faces[previ]) {
            previ += _faces.size() - 1;
            previ %= _faces.size();
        }
        unsigned int nexti = (i + 1) % _faces.size();
        while(!_faces[nexti]) {
            nexti += 1;
            nexti %= _faces.size();
        }
        //        double a[2], b[2];
        //        _faces[i]->get_points(a, b);
        //        double orientation = predicates::orient2d(a, b, _pos);
        double tngbpos[2], tnextpos[2], tprevpos[2];
        if(newlist.is_normal(_ngbs[i])) {
            AdaptiveVorCell2d* ngb = newlist.get_normal(_ngbs[i]);
            tngbpos[0] = ngb->_pos[0];
            tngbpos[1] = ngb->_pos[1];
            if(periodic && _ngbwalls[i] < 0) {
                AdaptiveMeshUtils::get_periodic_position(tngbpos, _ngbwalls[i],
                                                         box);
            }
        } else {
            if(newlist.is_ghost(_ngbs[i])) {
                AdaptiveVorCell2d* ghost = newlist.get_ghost(_ngbs[i]);
                AdaptiveMeshUtils::get_wall_position(ghost->_pos, tngbpos,
                                                     ghost->_ngbs[0], box);
            } else {
                AdaptiveVorCell2d* ngb = newlist.get_orphan(_ngbs[i]);
                tngbpos[0] = ngb->_pos[0];
                tngbpos[1] = ngb->_pos[1];
                if(periodic && _ngbwalls[i] < 0) {
                    AdaptiveMeshUtils::get_periodic_position(tngbpos,
                                                             _ngbwalls[i], box);
                }
            }
        }
        if(newlist.is_normal(_ngbs[nexti])) {
            AdaptiveVorCell2d* ngb = newlist.get_normal(_ngbs[nexti]);
            tnextpos[0] = ngb->_pos[0];
            tnextpos[1] = ngb->_pos[1];
            if(periodic && _ngbwalls[nexti] < 0) {
                AdaptiveMeshUtils::get_periodic_position(tnextpos,
                                                         _ngbwalls[nexti], box);
            }
        } else {
            if(newlist.is_ghost(_ngbs[nexti])) {
                AdaptiveVorCell2d* ghost = newlist.get_ghost(_ngbs[nexti]);
                AdaptiveMeshUtils::get_wall_position(ghost->_pos, tnextpos,
                                                     ghost->_ngbs[0], box);
            } else {
                AdaptiveVorCell2d* ngb = newlist.get_orphan(_ngbs[nexti]);
                tnextpos[0] = ngb->_pos[0];
                tnextpos[1] = ngb->_pos[1];
                if(periodic && _ngbwalls[nexti] < 0) {
                    AdaptiveMeshUtils::get_periodic_position(
                            tnextpos, _ngbwalls[nexti], box);
                }
            }
        }
        if(newlist.is_normal(_ngbs[previ])) {
            AdaptiveVorCell2d* ngb = newlist.get_normal(_ngbs[previ]);
            tprevpos[0] = ngb->_pos[0];
            tprevpos[1] = ngb->_pos[1];
            if(periodic && _ngbwalls[previ] < 0) {
                AdaptiveMeshUtils::get_periodic_position(tprevpos,
                                                         _ngbwalls[previ], box);
            }
        } else {
            if(newlist.is_ghost(_ngbs[previ])) {
                AdaptiveVorCell2d* ghost = newlist.get_ghost(_ngbs[previ]);
                AdaptiveMeshUtils::get_wall_position(ghost->_pos, tprevpos,
                                                     ghost->_ngbs[0], box);
            } else {
                AdaptiveVorCell2d* ngb = newlist.get_orphan(_ngbs[previ]);
                tprevpos[0] = ngb->_pos[0];
                tprevpos[1] = ngb->_pos[1];
                if(periodic && _ngbwalls[previ] < 0) {
                    AdaptiveMeshUtils::get_periodic_position(
                            tprevpos, _ngbwalls[previ], box);
                }
            }
        }
        double orientation =
                predicates::incircle_old(_pos, tngbpos, tprevpos, tnextpos);
        // note that for periodic boundaries and vertices of order four (or very
        // close) this test might give false results for cells at the
        // boundaries. The reason is that all generators are translated to the
        // frame of reference of the current cell, so that cells are treated
        // differently. The test is arbitrary exact, but the translation of
        // generators is not. So roundoff during translation can cause false
        // positives (or rather double positives/negatives which are mutually
        // exclusive). This could be solved by always translating all
        // coordinates to one of the periodic domains, but this would be
        // cumbersome. We choose to allow for wrong fourth order vertices, since
        // the associated faces have area zero anyway. The algorithm will not
        // crash on these and they will not affect the hydro either.
        //
        // Note however that this can cause errors to be detected when
        // explicitly checking the validity of the mesh (with
        // AdaptiveVorTess::check_mesh).
        if(orientation < 0.) {
            // enable this if you only want to test if all errors are detected
            //                continue;
            to_remove.insert(i);
            // make sure we can check that this ngb was removed
            delete _faces[i];
            _faces[i] = NULL;
            if(!newlist.is_ghost(_ngbs[i])) {
                // normal case: remove this cell from the neighbour that causes
                // the error
                if(newlist.is_normal(_ngbs[i])) {
                    newlist.get_normal(_ngbs[i])->remove_ngb(id, periodic);
                } else {
                    // if we do not do this, orphans keep a link to this cell,
                    // even if they are no longer a ngb. This is dangerous
                    // If an orphan is no longer linked, it should be removed
                    // however, we don't do this for now. So we have to make
                    // sure the ngbs of the orphan do not contain actual
                    // information
                    // this line causes trouble... Apparently, there are BAD
                    // orphans :(
                    newlist.get_orphan(_ngbs[i])->remove_ngb(id, periodic);
                }
            } else {
                // THIS CASE CANNOT HAPPEN WITH PERIODIC==TRUE!
                // special case: the current cell is moving away from the border
                // we have to treat this case separately, since we have to make
                // sure that all ghosts are removed (also those that are not a
                // mirror of this cell)
                if(newlist.get_ghost(_ngbs[i])->get_original() == id) {
                    //                    unsigned int previplu1 =
                    // (previ+1)%_ngbs.size();
                    //                    unsigned int nextimin1 =
                    // (nexti+_ngbs.size()-1)
                    //                    %_ngbs.size();
                    // search the first neighbour to the right that is not a
                    // mirror for this wall
                    while(newlist.is_ghost(_ngbs[previ]) &&
                          newlist.get_ghost(_ngbs[previ])->_ngbs[0] ==
                                  newlist.get_ghost(_ngbs[i])->_ngbs[0]) {
                        // fix for the case where previ is already removed by a
                        // previous detection
                        // in this case, previ+1 will not be the first ngb of
                        // previ anymore (after id), in which case nexti may be
                        // inserted in the wrong place in the ngb list for previ
                        if(!to_remove.count(previ)) {
                            //                            previplu1 = previ;
                            to_remove.insert(previ);
                        }
                        // important: make sure the face is not processed later
                        // in the loop over i by deleting it
                        delete _faces[previ];
                        _faces[previ] = NULL;
                        previ += _ngbs.size() - 1;
                        previ %= _ngbs.size();
                    }
                    // same as above, but for the left side
                    while(newlist.is_ghost(_ngbs[nexti]) &&
                          newlist.get_ghost(_ngbs[nexti])->_ngbs[0] ==
                                  newlist.get_ghost(_ngbs[i])->_ngbs[0]) {
                        if(!to_remove.count(nexti)) {
                            //                            nextimin1 = nexti;
                            to_remove.insert(nexti);
                        }
                        delete _faces[nexti];
                        _faces[nexti] = NULL;
                        nexti += 1;
                        nexti %= _ngbs.size();
                    }
                    // the ngb can also be another wall; in this case we don't
                    // do anything for this ngb
                    if(newlist.is_normal(_ngbs[previ])) {
                        //                        cells[_ngbs[previ]]->add_ngb(_ngbs[previplu1],
                        // id,
                        //                        _ngbs[nexti]);
                        newlist.get_normal(_ngbs[previ])
                                ->add_ngb_before(id, _ngbs[nexti]);
                        // it could happen that _ngbs[i] is not a ngb of
                        // _ngbs[previ], so we use the safe method
                        newlist.get_normal(_ngbs[previ])
                                ->safely_remove_ngb(_ngbs[i], periodic);
                        newlist.get_normal(_ngbs[previ])
                                ->complete(newlist, box, periodic);
                    }
                    // same as above
                    if(newlist.is_normal(_ngbs[nexti])) {
                        //                        cells[_ngbs[nexti]]->add_ngb(_ngbs[nextimin1],
                        // id,
                        //                        _ngbs[previ]);
                        newlist.get_normal(_ngbs[nexti])
                                ->add_ngb_after(id, _ngbs[previ]);
                        newlist.get_normal(_ngbs[nexti])
                                ->safely_remove_ngb(_ngbs[i], periodic);
                        newlist.get_normal(_ngbs[nexti])
                                ->complete(newlist, box, periodic);
                    }
                    // skip the rest of the operations
                    continue;
                }
            }
            if(newlist.is_normal(_ngbs[previ])) {
                int to_change = _ngbs[nexti];
                // special case where a cell touches the border: we have to add
                // a new ghost for this cell
                if(newlist.is_ghost(_ngbs[nexti]) &&
                   !newlist.get_normal(_ngbs[previ])->_wall) {
                    // THIS CASE CANNOT HAPPEN WHEN PERIODIC==TRUE!
                    to_change = newlist.add_ghost_cell(new AdaptiveVorCell2d(
                            newlist.get_normal(_ngbs[previ])->_pos, 0,
                            _ngbs[previ], NULL));
                    newlist.ghostback()->add_ngb(
                            newlist.get_ghost(_ngbs[nexti])->_ngbs[0]);
                    // at least one of the cells (id or _ngbs[i]) contains the
                    // ghost of the other as a ngb. We change this ghost to the
                    // newly created ghost
                    if(id == newlist.get_ghost(_ngbs[nexti])->get_original()) {
                        newlist.get_normal(_ngbs[i])
                                ->safely_change_ngb_from_value(_ngbs[nexti],
                                                               to_change);
                    } else {
                        // if nexti was already signaled to be removed by a
                        // previous switch, make sure it is not removed by
                        // erasing it from the remove queue
                        if(to_remove.find(nexti) != to_remove.end()) {
                            to_remove.erase(nexti);
                            // it can happen that _ngbs[i] contains the mirror
                            // of id (which is then nexti+1) due to a previous
                            // switch by safely removing it, we are sure that it
                            // is removed
                            newlist.get_normal(_ngbs[i])->safely_remove_ngb(
                                    _ngbs[(nexti + 1) % _ngbs.size()]);
                        }
                        _ngbs[nexti] = to_change;
                        // this ngb changed, so it is unwise to still keep its
                        // associated face...
                        delete _faces[nexti];
                        _faces[nexti] = NULL;
                    }
                }
                if(periodic) {
                    int ngbwall = _ngbwalls[nexti];
                    ngbwall = AdaptiveMeshUtils::get_wallpos(_ngbwalls[previ],
                                                             ngbwall);
                    // for cells with fourth order vertices at the periodic
                    // boundary, neighbour relations can be a bit tricky at the
                    // start. This condition ensures we never add a neighbour
                    // twice
                    if(!newlist.get_normal(_ngbs[previ])->is_ngb(to_change)) {
                        newlist.get_normal(_ngbs[previ])
                                ->add_ngb(_ngbs[i], id, to_change, ngbwall);
                    }
                } else {
                    newlist.get_normal(_ngbs[previ])
                            ->add_ngb(_ngbs[i], id, to_change);
                }
                newlist.get_normal(_ngbs[previ])
                        ->complete(newlist, box, periodic);
            } else {
                if(!newlist.is_ghost(_ngbs[previ]) &&
                   newlist.is_normal(_ngbs[nexti]) &&
                   newlist.is_normal(_ngbs[i])) {
                    // add to insertions:
                    //  - the rank to send to
                    //  - the original id of the cell on rank
                    //  - the local ids (on this process) of the ngbs
                    //  - the id of this cell
                    //  - the ngbwall
                    AdaptiveVorCell2d* orphan =
                            newlist.get_orphan(_ngbs[previ]);
                    insertions.push_back(orphan->get_rank());
                    insertions.push_back(orphan->get_original());
                    insertions.push_back(id);
                    insertions.push_back(_ngbs[i]);
                    insertions.push_back(_ngbs[nexti]);
                    int ngbwall = _ngbwalls[nexti];
                    ngbwall = AdaptiveMeshUtils::get_wallpos(_ngbwalls[previ],
                                                             ngbwall);
                    insertions.push_back(ngbwall);
                }
            }
            if(newlist.is_normal(_ngbs[nexti])) {
                int to_change = _ngbs[previ];
                // see above
                if(newlist.is_ghost(_ngbs[previ]) &&
                   !newlist.get_normal(_ngbs[nexti])->_wall) {
                    // THIS CASE CANNOT HAPPEN WHEN PERIODIC==TRUE!
                    to_change = newlist.add_ghost_cell(new AdaptiveVorCell2d(
                            newlist.get_normal(_ngbs[nexti])->_pos, 0,
                            _ngbs[nexti], NULL));
                    newlist.ghostback()->add_ngb(
                            newlist.get_ghost(_ngbs[previ])->_ngbs[0]);
                    if(id == newlist.get_ghost(_ngbs[previ])->get_original()) {
                        newlist.get_normal(_ngbs[i])
                                ->safely_change_ngb_from_value(_ngbs[previ],
                                                               to_change);
                    } else {
                        if(to_remove.find(previ) != to_remove.end()) {
                            to_remove.erase(previ);
                            newlist.get_normal(_ngbs[i])->safely_remove_ngb(
                                    _ngbs[(previ + _ngbs.size() - 1) %
                                          _ngbs.size()]);
                        }
                        _ngbs[previ] = to_change;
                        delete _faces[previ];
                        _faces[previ] = NULL;
                    }
                }
                if(periodic) {
                    int ngbwall = _ngbwalls[previ];
                    ngbwall = AdaptiveMeshUtils::get_wallpos(_ngbwalls[nexti],
                                                             ngbwall);
                    // see comment above
                    if(!newlist.get_normal(_ngbs[nexti])->is_ngb(to_change)) {
                        newlist.get_normal(_ngbs[nexti])
                                ->add_ngb(_ngbs[i], id, to_change, ngbwall);
                    }
                } else {
                    newlist.get_normal(_ngbs[nexti])
                            ->add_ngb(_ngbs[i], id, to_change);
                }
                newlist.get_normal(_ngbs[nexti])
                        ->complete(newlist, box, periodic);
            } else {
                if(!newlist.is_ghost(_ngbs[nexti]) &&
                   newlist.is_normal(_ngbs[previ]) &&
                   newlist.is_normal(_ngbs[i])) {
                    // add to insertions:
                    //  - the rank to send to
                    //  - the original id of the cell on rank
                    //  - the local ids (on this process) of the ngbs
                    //  - the id of this cell
                    //  - the ngbwall
                    AdaptiveVorCell2d* orphan =
                            newlist.get_orphan(_ngbs[nexti]);
                    insertions.push_back(orphan->get_rank());
                    insertions.push_back(orphan->get_original());
                    insertions.push_back(id);
                    insertions.push_back(_ngbs[i]);
                    insertions.push_back(_ngbs[previ]);
                    int ngbwall = _ngbwalls[previ];
                    ngbwall = AdaptiveMeshUtils::get_wallpos(_ngbwalls[nexti],
                                                             ngbwall);
                    insertions.push_back(ngbwall);
                }
            }
            if(newlist.is_normal(_ngbs[i])) {
                newlist.get_normal(_ngbs[i])->complete(newlist, box, periodic);
                // we cannot check this ngb now, since this can cause spurious
                // changes in the ngb list of the current cell
                to_check.push_back(_ngbs[i]);
            }
        }
    }
    // order is important! (removing the first element first will shift
    // the other elements we have to remove...)
    //    sort(to_remove.begin(), to_remove.end());
    //    for(unsigned int i = to_remove.size(); i--;){
    for(set<unsigned int>::reverse_iterator it = to_remove.rbegin();
        it != to_remove.rend(); it++) {
        _ngbs.erase(_ngbs.begin() + *it);
        _faces.erase(_faces.begin() + *it);
        _exngbs.erase(_exngbs.end() - 1);
        if(periodic) {
            _ngbwalls.erase(_ngbwalls.begin() + *it);
        }
    }
    //    for(unsigned int i = 0; i < to_check.size(); i++){
    //        newlist.get_normal(to_check[i])->detect_crossovers(to_check[i],
    //    newlist, insertions, box, periodic);
    //    }
    return to_remove.size() > 0;
}

/**
 * @brief Kernel of the mesh evolution algorithm: detect errors in the evolved
 * mesh and correct them
 *
 * @param id Index of the cell in the AdaptiveVorTess cell list
 * @param cells Reference to the AdaptiveVorTess cell list
 * @param positionlist NewCellList holding accurate particle positions
 * @param box Cuboid specifying the dimensions of the simulation box
 * @param insertions List with undetected insertions on other processes
 * @param periodic Flag discriminating between a periodic (true) and a
 * reflective (false) simulation box
 * @return True if errors where corrected, false otherwise
 */
bool AdaptiveVorCell2d::detect_crossovers(
        unsigned int id, std::vector<AdaptiveVorCell2d*>& cells,
        NewCellList& positionlist, Cuboid& box,
        std::vector<unsigned long>& insertions, bool periodic) {
    complete(positionlist);
    compute_winding_number();
    set<unsigned int> to_remove;
    //    vector<unsigned int> to_check;
    for(unsigned int i = 0; i < _faces.size(); i++) {
        if(!_faces[i]) {
            continue;
        }
        unsigned int previ = (i + _faces.size() - 1) % _faces.size();
        // it can happen that the previous or next face was already removed
        // earlier on
        while(!_faces[previ]) {
            previ += _faces.size() - 1;
            previ %= _faces.size();
        }
        unsigned int nexti = (i + 1) % _faces.size();
        while(!_faces[nexti]) {
            nexti += 1;
            nexti %= _faces.size();
        }
        //        double a[2], b[2];
        //        _faces[i]->get_points(a, b);
        //        double orientation = predicates::orient2d(a, b, _pos);
        double tngbpos[2], tnextpos[2], tprevpos[2];
        //        if(newlist.is_normal(_ngbs[i])){
        //            AdaptiveVorCell2d* ngb = newlist.get_normal(_ngbs[i]);
        //            tngbpos[0] = ngb->_pos[0];
        //            tngbpos[1] = ngb->_pos[1];
        //            if(periodic && _ngbwalls[i] < 0){
        //                AdaptiveMeshUtils::get_periodic_position(tngbpos,
        //                _ngbwalls[i],
        //                                                         box);
        //            }
        //        } else {
        //            if(newlist.is_ghost(_ngbs[i])){
        //                AdaptiveVorCell2d* ghost =
        //                newlist.get_ghost(_ngbs[i]);
        //                AdaptiveMeshUtils::get_wall_position(ghost->_pos,
        //                tngbpos,
        //                                                     ghost->_ngbs[0],
        //                                                     box);
        //            } else {
        //                AdaptiveVorCell2d* ngb = newlist.get_orphan(_ngbs[i]);
        //                tngbpos[0] = ngb->_pos[0];
        //                tngbpos[1] = ngb->_pos[1];
        //                if(periodic && _ngbwalls[i] < 0){
        //                    AdaptiveMeshUtils::get_periodic_position(tngbpos,
        //                                                             _ngbwalls[i],
        //                                                             box);
        //                }
        //            }
        //        }
        Vec pos = positionlist.get_position(_newngbs[i]);
        tngbpos[0] = pos.x();
        tngbpos[1] = pos.y();
        //        if(newlist.is_normal(_ngbs[nexti])){
        //            AdaptiveVorCell2d* ngb = newlist.get_normal(_ngbs[nexti]);
        //            tnextpos[0] = ngb->_pos[0];
        //            tnextpos[1] = ngb->_pos[1];
        //            if(periodic && _ngbwalls[nexti] < 0){
        //                AdaptiveMeshUtils::get_periodic_position(tnextpos,
        //                                                         _ngbwalls[nexti],
        //                                                         box);
        //            }
        //        } else {
        //            if(newlist.is_ghost(_ngbs[nexti])){
        //                AdaptiveVorCell2d* ghost =
        //                newlist.get_ghost(_ngbs[nexti]);
        //                AdaptiveMeshUtils::get_wall_position(ghost->_pos,
        //                tnextpos,
        //                                                     ghost->_ngbs[0],
        //                                                     box);
        //            } else {
        //                AdaptiveVorCell2d* ngb =
        //                newlist.get_orphan(_ngbs[nexti]);
        //                tnextpos[0] = ngb->_pos[0];
        //                tnextpos[1] = ngb->_pos[1];
        //                if(periodic && _ngbwalls[nexti] < 0){
        //                    AdaptiveMeshUtils::get_periodic_position(tnextpos,
        //                                                             _ngbwalls[nexti],
        //                                                             box);
        //                }
        //            }
        //        }
        pos = positionlist.get_position(_newngbs[nexti]);
        tnextpos[0] = pos.x();
        tnextpos[1] = pos.y();
        //        if(newlist.is_normal(_ngbs[previ])){
        //            AdaptiveVorCell2d* ngb = newlist.get_normal(_ngbs[previ]);
        //            tprevpos[0] = ngb->_pos[0];
        //            tprevpos[1] = ngb->_pos[1];
        //            if(periodic && _ngbwalls[previ] < 0){
        //                AdaptiveMeshUtils::get_periodic_position(tprevpos,
        //                                                         _ngbwalls[previ],
        //                                                         box);
        //            }
        //        } else {
        //            if(newlist.is_ghost(_ngbs[previ])){
        //                AdaptiveVorCell2d* ghost =
        //                newlist.get_ghost(_ngbs[previ]);
        //                AdaptiveMeshUtils::get_wall_position(ghost->_pos,
        //                tprevpos,
        //                                                     ghost->_ngbs[0],
        //                                                     box);
        //            } else {
        //                AdaptiveVorCell2d* ngb =
        //                newlist.get_orphan(_ngbs[previ]);
        //                tprevpos[0] = ngb->_pos[0];
        //                tprevpos[1] = ngb->_pos[1];
        //                if(periodic && _ngbwalls[previ] < 0){
        //                    AdaptiveMeshUtils::get_periodic_position(tprevpos,
        //                                                             _ngbwalls[previ],
        //                                                             box);
        //                }
        //            }
        //        }
        pos = positionlist.get_position(_newngbs[previ]);
        tprevpos[0] = pos.x();
        tprevpos[1] = pos.y();
        double orientation =
                predicates::incircle_old(_pos, tngbpos, tprevpos, tnextpos);
        // note that for periodic boundaries and vertices of order four (or very
        // close) this test might give false results for cells at the
        // boundaries. The reason is that all generators are translated to the
        // frame of reference of the current cell, so that cells are treated
        // differently. The test is arbitrary exact, but the translation of
        // generators is not. So roundoff during translation can cause false
        // positives (or rather double positives/negatives which are mutually
        // exclusive). This could be solved by always translating all
        // coordinates to one of the periodic domains, but this would be
        // cumbersome. We choose to allow for wrong fourth order vertices, since
        // the associated faces have area zero anyway. The algorithm will not
        // crash on these and they will not affect the hydro either.
        //
        // Note however that this can cause errors to be detected when
        // explicitly checking the validity of the mesh (with
        // AdaptiveVorTess::check_mesh).
        if(orientation < 0.) {
            // enable this if you only want to test if all errors are detected
            //                continue;
            to_remove.insert(i);
            // make sure we can check that this ngb was removed
            delete _faces[i];
            _faces[i] = NULL;
            // only periodic boundaries at the moment
            if(positionlist.is_local(_newngbs[i])) {
                unsigned int index = positionlist.get_locid(_newngbs[i]);
                int wall = positionlist.get_wall(_newngbs[i]);
                //                cerr << _newngbs[i] << "\t" << wall << endl;
                unsigned long to_remove = _particleID;
                if(wall) {
                    // convert the ID to the appropriate copy
                    wall = AdaptiveMeshUtils::get_wallpos(wall, 0);
                    // _particleID will always be a real particle, we do not
                    // have to convert the ID
                    to_remove = positionlist.get_ghost(_particleID, wall);
                }
                // if _newngbs[i] is on another process, the flip will be
                // detected there as well and the removal will be done there
                cells[index]->remove_newngb(to_remove, periodic);
            }
            //            if(!newlist.is_ghost(_ngbs[i])){
            //                // normal case: remove this cell from the
            //                neighbour that causes
            //                // the error
            //                if(newlist.is_normal(_ngbs[i])){
            //                    newlist.get_normal(_ngbs[i])->remove_ngb(id,
            //                    periodic);
            //                } else {
            //                    // if we do not do this, orphans keep a link
            //                    to this cell,
            //                    // even if they are no longer a ngb. This is
            //                    dangerous
            //                    // If an orphan is no longer linked, it should
            //                    be removed
            //                    // however, we don't do this for now. So we
            //                    have to make
            //                    // sure the ngbs of the orphan do not contain
            //                    actual
            //                    // information
            //                    // this line causes trouble... Apparently,
            //                    there are BAD
            //                    // orphans :(
            //                    newlist.get_orphan(_ngbs[i])->remove_ngb(id,
            //                    periodic);
            //                }
            //            } else {
            //                // THIS CASE CANNOT HAPPEN WITH PERIODIC==TRUE!
            //                // special case: the current cell is moving away
            //                from the border
            //                // we have to treat this case separately, since we
            //                have to make
            //                // sure that all ghosts are removed (also those
            //                that are not a
            //                // mirror of this cell)
            //                if(newlist.get_ghost(_ngbs[i])->get_original() ==
            //                id){
            ////                    unsigned int previplu1 =
            ///(previ+1)%_ngbs.size();
            ////                    unsigned int nextimin1 =
            ///(nexti+_ngbs.size()-1)
            ////                    %_ngbs.size();
            //                    // search the first neighbour to the right
            //                    that is not a
            //                    // mirror for this wall
            //                    while(newlist.is_ghost(_ngbs[previ]) &&
            //                          newlist.get_ghost(_ngbs[previ])->_ngbs[0]
            //                          ==
            //                          newlist.get_ghost(_ngbs[i])->_ngbs[0]){
            //                        // fix for the case where previ is already
            //                        removed by a
            //                        // previous detection
            //                        // in this case, previ+1 will not be the
            //                        first ngb of
            //                        // previ anymore (after id), in which case
            //                        nexti may be
            //                        // inserted in the wrong place in the ngb
            //                        list for previ
            //                        if(!to_remove.count(previ)){
            ////                            previplu1 = previ;
            //                            to_remove.insert(previ);
            //                        }
            //                        // important: make sure the face is not
            //                        processed later
            //                        // in the loop over i by deleting it
            //                        delete _faces[previ];
            //                        _faces[previ] = NULL;
            //                        previ += _ngbs.size()-1;
            //                        previ %= _ngbs.size();
            //                    }
            //                    // same as above, but for the left side
            //                    while(newlist.is_ghost(_ngbs[nexti]) &&
            //                          newlist.get_ghost(_ngbs[nexti])->_ngbs[0]
            //                          ==
            //                          newlist.get_ghost(_ngbs[i])->_ngbs[0]){
            //                        if(!to_remove.count(nexti)){
            ////                            nextimin1 = nexti;
            //                            to_remove.insert(nexti);
            //                        }
            //                        delete _faces[nexti];
            //                        _faces[nexti] = NULL;
            //                        nexti += 1;
            //                        nexti %= _ngbs.size();
            //                    }
            //                    // the ngb can also be another wall; in this
            //                    case we don't
            //                    // do anything for this ngb
            //                    if(newlist.is_normal(_ngbs[previ])){
            //// cells[_ngbs[previ]]->add_ngb(_ngbs[previplu1], id,
            ////                        _ngbs[nexti]);
            //                        newlist.get_normal(_ngbs[previ])
            //                                ->add_ngb_before(id,
            //                                _ngbs[nexti]);
            //                        // it could happen that _ngbs[i] is not a
            //                        ngb of
            //                        // _ngbs[previ], so we use the safe method
            //                        newlist.get_normal(_ngbs[previ])
            //                                ->safely_remove_ngb(_ngbs[i],
            //                                periodic);
            //                        newlist.get_normal(_ngbs[previ])
            //                                ->complete(newlist, box,
            //                                periodic);
            //                    }
            //                    // same as above
            //                    if(newlist.is_normal(_ngbs[nexti])){
            //// cells[_ngbs[nexti]]->add_ngb(_ngbs[nextimin1], id,
            ////                        _ngbs[previ]);
            //                        newlist.get_normal(_ngbs[nexti])
            //                                ->add_ngb_after(id, _ngbs[previ]);
            //                        newlist.get_normal(_ngbs[nexti])
            //                                ->safely_remove_ngb(_ngbs[i],
            //                                periodic);
            //                        newlist.get_normal(_ngbs[nexti])
            //                                ->complete(newlist, box,
            //                                periodic);
            //                    }
            //                    // skip the rest of the operations
            //                    continue;
            //                }
            //            }
            if(positionlist.is_local(_newngbs[previ])) {
                unsigned long to_change = _newngbs[nexti];
                unsigned int index = positionlist.get_locid(_newngbs[previ]);
                // none of the IDs needs to be correct when periodic boundaries
                // are used...
                int leftwall = positionlist.get_wall(_newngbs[previ]);
                int ngbwall = positionlist.get_wall(_newngbs[i]);
                int rightwall = positionlist.get_wall(_newngbs[nexti]);
                // these IDs might be copies, we have to convert them to real
                // particle IDs
                unsigned long id1 = positionlist.get_id(_newngbs[i]);
                unsigned long id2 = positionlist.get_id(_particleID);
                unsigned long newid = positionlist.get_id(to_change);
                // we need the correct periodic copy of id1 = _newngbs[i]
                int id1wall = AdaptiveMeshUtils::get_wallpos(leftwall, ngbwall);
                if(id1wall) {
                    id1 = positionlist.get_ghost(id1, id1wall);
                }
                // we need the correct periodic copy of id2 = _particleID
                int id2wall = AdaptiveMeshUtils::get_wallpos(leftwall, 0);
                if(id2wall) {
                    id2 = positionlist.get_ghost(id2, id2wall);
                }
                // we need the correct periodic copy of newid = to_change
                int newwall =
                        AdaptiveMeshUtils::get_wallpos(leftwall, rightwall);
                if(newwall) {
                    newid = positionlist.get_ghost(newid, newwall);
                }
                cells[index]->add_newngb(id1, id2, newid);
                cells[index]->complete(positionlist);
            } else {
                // it is possible that the insertion in _newngbs[previ] is not
                // detected on its home process...
                // this happens if _newngbs[i] is not on the same process as
                // _newngbs[previ]
                // however, if _newngbs[i] is on another process than this cell,
                // the insertion might also be communicated from there and we
                // only communicate from the cell with the lowest id
                if(positionlist.get_rank(_newngbs[i]) !=
                   positionlist.get_rank(_newngbs[previ])) {
                    unsigned long to_change = _newngbs[nexti];
                    // these IDs might be copies, we have to convert them to
                    // real
                    // particle IDs
                    unsigned long id1 = positionlist.get_id(_newngbs[i]);
                    unsigned long id2 = positionlist.get_id(_particleID);
                    if(positionlist.get_rank(_newngbs[i]) == MPIGlobal::rank ||
                       id1 < id2) {
                        unsigned long newid = positionlist.get_id(to_change);
                        // none of the IDs needs to be correct when periodic
                        // boundaries
                        // are used...
                        int leftwall = positionlist.get_wall(_newngbs[previ]);
                        int ngbwall = positionlist.get_wall(_newngbs[i]);
                        int rightwall = positionlist.get_wall(_newngbs[nexti]);
                        // we need the correct periodic copy of id1 =
                        // _newngbs[i]
                        int id1wall = AdaptiveMeshUtils::get_wallpos(leftwall,
                                                                     ngbwall);
                        if(id1wall) {
                            id1 = positionlist.get_ghost(id1, id1wall);
                        }
                        // we need the correct periodic copy of id2 =
                        // _particleID
                        int id2wall =
                                AdaptiveMeshUtils::get_wallpos(leftwall, 0);
                        if(id2wall) {
                            id2 = positionlist.get_ghost(id2, id2wall);
                        }
                        // we need the correct periodic copy of newid =
                        // to_change
                        int newwall = AdaptiveMeshUtils::get_wallpos(leftwall,
                                                                     rightwall);
                        if(newwall) {
                            newid = positionlist.get_ghost(newid, newwall);
                        }
                        insertions.push_back(_newngbs[previ]);
                        insertions.push_back(id1);
                        insertions.push_back(id2);
                        insertions.push_back(newid);
                    }
                }
            }
            //            if(newlist.is_normal(_ngbs[previ])){
            //                int to_change = _ngbs[nexti];
            //                // special case where a cell touches the border:
            //                we have to add
            //                // a new ghost for this cell
            //                if(newlist.is_ghost(_ngbs[nexti]) &&
            //                        !newlist.get_normal(_ngbs[previ])->_wall){
            //                    // THIS CASE CANNOT HAPPEN WHEN
            //                    PERIODIC==TRUE!
            //                    to_change = newlist.add_ghost_cell(
            //                                new AdaptiveVorCell2d(
            //                                    newlist.get_normal(_ngbs[previ])->_pos,
            //                                    0,
            //                                    _ngbs[previ], NULL)
            //                                );
            //                    newlist.ghostback()->add_ngb(
            //                                newlist.get_ghost(_ngbs[nexti])->_ngbs[0]
            //                                );
            //                    // at least one of the cells (id or _ngbs[i])
            //                    contains the
            //                    // ghost of the other as a ngb. We change this
            //                    ghost to the
            //                    // newly created ghost
            //                    if(id ==
            //                    newlist.get_ghost(_ngbs[nexti])->get_original()){
            //                        newlist.get_normal(_ngbs[i])
            //                                ->safely_change_ngb_from_value(_ngbs[nexti],
            //                                                               to_change);
            //                    } else {
            //                        // if nexti was already signaled to be
            //                        removed by a
            //                        // previous switch, make sure it is not
            //                        removed by
            //                        // erasing it from the remove queue
            //                        if(to_remove.find(nexti) !=
            //                        to_remove.end()){
            //                            to_remove.erase(nexti);
            //                            // it can happen that _ngbs[i]
            //                            contains the mirror
            //                            // of id (which is then nexti+1) due
            //                            to a previous
            //                            // switch by safely removing it, we
            //                            are sure that it
            //                            // is removed
            //                            newlist.get_normal(_ngbs[i])
            //                                    ->safely_remove_ngb(_ngbs[(nexti+1)
            //                                                        %_ngbs.size()]);
            //                        }
            //                        _ngbs[nexti] = to_change;
            //                        // this ngb changed, so it is unwise to
            //                        still keep its
            //                        // associated face...
            //                        delete _faces[nexti];
            //                        _faces[nexti] = NULL;
            //                    }
            //                }
            //                if(periodic){
            //                    int ngbwall = _ngbwalls[nexti];
            //                    ngbwall =
            //                    AdaptiveMeshUtils::get_wallpos(_ngbwalls[previ],
            //                                                             ngbwall);
            //                    // for cells with fourth order vertices at the
            //                    periodic
            //                    // boundary, neighbour relations can be a bit
            //                    tricky at the
            //                    //start. This condition ensures we never add a
            //                    neighbour
            //                    // twice
            //                    if(!newlist.get_normal(_ngbs[previ])->is_ngb(to_change)){
            //                        newlist.get_normal(_ngbs[previ])
            //                                ->add_ngb(_ngbs[i], id, to_change,
            //                                ngbwall);
            //                    }
            //                } else {
            //                    newlist.get_normal(_ngbs[previ])
            //                            ->add_ngb(_ngbs[i], id, to_change);
            //                }
            //                newlist.get_normal(_ngbs[previ])->complete(newlist,
            //                box,
            //                                                           periodic);
            //            } else {
            //                if(!newlist.is_ghost(_ngbs[previ]) &&
            //                        newlist.is_normal(_ngbs[nexti]) &&
            //                        newlist.is_normal(_ngbs[i])){
            //                    // add to insertions:
            //                    //  - the rank to send to
            //                    //  - the original id of the cell on rank
            //                    //  - the local ids (on this process) of the
            //                    ngbs
            //                    //  - the id of this cell
            //                    //  - the ngbwall
            //                    AdaptiveVorCell2d* orphan =
            //                    newlist.get_orphan(_ngbs[previ]);
            //                    insertions.push_back(orphan->get_rank());
            //                    insertions.push_back(orphan->get_original());
            //                    insertions.push_back(id);
            //                    insertions.push_back(_ngbs[i]);
            //                    insertions.push_back(_ngbs[nexti]);
            //                    int ngbwall = _ngbwalls[nexti];
            //                    ngbwall =
            //                    AdaptiveMeshUtils::get_wallpos(_ngbwalls[previ],
            //                                                             ngbwall);
            //                    insertions.push_back(ngbwall);
            //                }
            //            }
            if(positionlist.is_local(_newngbs[nexti])) {
                unsigned long to_change = _newngbs[previ];
                unsigned int index = positionlist.get_locid(_newngbs[nexti]);
                // none of the IDs needs to be correct when periodic boundaries
                // are used...
                int leftwall = positionlist.get_wall(_newngbs[previ]);
                int ngbwall = positionlist.get_wall(_newngbs[i]);
                int rightwall = positionlist.get_wall(_newngbs[nexti]);
                // these IDs might be copies, we have to convert them to real
                // particle IDs
                unsigned long id1 = positionlist.get_id(_newngbs[i]);
                unsigned long id2 = positionlist.get_id(_particleID);
                unsigned long newid = positionlist.get_id(to_change);
                // we need the correct periodic copy of id1 = _newngbs[i]
                int id1wall =
                        AdaptiveMeshUtils::get_wallpos(rightwall, ngbwall);
                if(id1wall) {
                    id1 = positionlist.get_ghost(id1, id1wall);
                }
                // we need the correct periodic copy of id2 = _particleID
                int id2wall = AdaptiveMeshUtils::get_wallpos(rightwall, 0);
                if(id2wall) {
                    id2 = positionlist.get_ghost(id2, id2wall);
                }
                // we need the correct periodic copy of newid = to_change
                int newwall =
                        AdaptiveMeshUtils::get_wallpos(rightwall, leftwall);
                if(newwall) {
                    newid = positionlist.get_ghost(newid, newwall);
                }
                cells[index]->add_newngb(id1, id2, newid);
                cells[index]->complete(positionlist);
            } else {
                // it is possible that the insertion in _newngbs[nexti] is not
                // detected on its home process...
                // this happens if _newngbs[i] is not on the same process as
                // _newngbs[previ]
                // however, if _newngbs[i] is on another process than this cell,
                // the insertion might also be communicated from there and we
                // only communicate from the cell with the lowest id
                if(positionlist.get_rank(_newngbs[i]) !=
                   positionlist.get_rank(_newngbs[nexti])) {
                    unsigned long to_change = _newngbs[previ];
                    // these IDs might be copies, we have to convert them to
                    // real
                    // particle IDs
                    unsigned long id1 = positionlist.get_id(_newngbs[i]);
                    unsigned long id2 = positionlist.get_id(_particleID);
                    if(positionlist.get_rank(_newngbs[i]) == MPIGlobal::rank ||
                       id1 < id2) {
                        unsigned long newid = positionlist.get_id(to_change);
                        // none of the IDs needs to be correct when periodic
                        // boundaries
                        // are used...
                        int leftwall = positionlist.get_wall(_newngbs[previ]);
                        int ngbwall = positionlist.get_wall(_newngbs[i]);
                        int rightwall = positionlist.get_wall(_newngbs[nexti]);
                        // we need the correct periodic copy of id1 =
                        // _newngbs[i]
                        int id1wall = AdaptiveMeshUtils::get_wallpos(rightwall,
                                                                     ngbwall);
                        if(id1wall) {
                            id1 = positionlist.get_ghost(id1, id1wall);
                        }
                        // we need the correct periodic copy of id2 =
                        // _particleID
                        int id2wall =
                                AdaptiveMeshUtils::get_wallpos(rightwall, 0);
                        if(id2wall) {
                            id2 = positionlist.get_ghost(id2, id2wall);
                        }
                        // we need the correct periodic copy of newid =
                        // to_change
                        int newwall = AdaptiveMeshUtils::get_wallpos(rightwall,
                                                                     leftwall);
                        if(newwall) {
                            newid = positionlist.get_ghost(newid, newwall);
                        }
                        insertions.push_back(_newngbs[nexti]);
                        insertions.push_back(id1);
                        insertions.push_back(id2);
                        insertions.push_back(newid);
                    }
                }
            }
            //            if(newlist.is_normal(_ngbs[nexti])){
            //                int to_change = _ngbs[previ];
            //                // see above
            //                if(newlist.is_ghost(_ngbs[previ]) &&
            //                        !newlist.get_normal(_ngbs[nexti])->_wall){
            //                    // THIS CASE CANNOT HAPPEN WHEN
            //                    PERIODIC==TRUE!
            //                    to_change = newlist.add_ghost_cell(
            //                                new AdaptiveVorCell2d(
            //                                    newlist.get_normal(_ngbs[nexti])->_pos,
            //                                    0,
            //                                    _ngbs[nexti], NULL)
            //                                );
            //                    newlist.ghostback()->add_ngb(
            //                                newlist.get_ghost(_ngbs[previ])->_ngbs[0]
            //                                );
            //                    if(id ==
            //                    newlist.get_ghost(_ngbs[previ])->get_original()){
            //                        newlist.get_normal(_ngbs[i])
            //                                ->safely_change_ngb_from_value(_ngbs[previ],
            //                                                               to_change);
            //                    } else {
            //                        if(to_remove.find(previ) !=
            //                        to_remove.end()){
            //                            to_remove.erase(previ);
            //                            newlist.get_normal(_ngbs[i])->safely_remove_ngb(
            //                                        _ngbs[(previ+_ngbs.size()-1)
            //                                        %_ngbs.size()]);
            //                        }
            //                        _ngbs[previ] = to_change;
            //                        delete _faces[previ];
            //                        _faces[previ] = NULL;
            //                    }
            //                }
            //                if(periodic){
            //                    int ngbwall = _ngbwalls[previ];
            //                    ngbwall =
            //                    AdaptiveMeshUtils::get_wallpos(_ngbwalls[nexti],
            //                                                             ngbwall);
            //                    // see comment above
            //                    if(!newlist.get_normal(_ngbs[nexti])->is_ngb(to_change)){
            //                        newlist.get_normal(_ngbs[nexti])->add_ngb(_ngbs[i],
            //                        id,
            //                                                                  to_change,
            //                                                                  ngbwall);
            //                    }
            //                } else {
            //                    newlist.get_normal(_ngbs[nexti])->add_ngb(_ngbs[i],
            //                    id,
            //                                                              to_change);
            //                }
            //                newlist.get_normal(_ngbs[nexti])->complete(newlist,
            //                box,
            //                                                           periodic);
            //            } else {
            //                if(!newlist.is_ghost(_ngbs[nexti]) &&
            //                        newlist.is_normal(_ngbs[previ]) &&
            //                        newlist.is_normal(_ngbs[i])){
            //                    // add to insertions:
            //                    //  - the rank to send to
            //                    //  - the original id of the cell on rank
            //                    //  - the local ids (on this process) of the
            //                    ngbs
            //                    //  - the id of this cell
            //                    //  - the ngbwall
            //                    AdaptiveVorCell2d* orphan =
            //                    newlist.get_orphan(_ngbs[nexti]);
            //                    insertions.push_back(orphan->get_rank());
            //                    insertions.push_back(orphan->get_original());
            //                    insertions.push_back(id);
            //                    insertions.push_back(_ngbs[i]);
            //                    insertions.push_back(_ngbs[previ]);
            //                    int ngbwall = _ngbwalls[previ];
            //                    ngbwall =
            //                    AdaptiveMeshUtils::get_wallpos(_ngbwalls[nexti],
            //                                                             ngbwall);
            //                    insertions.push_back(ngbwall);
            //                }
            //            }
            if(positionlist.is_local(_newngbs[i])) {
                unsigned int index = positionlist.get_locid(_newngbs[i]);
                cells[index]->complete(positionlist);
            }
            //            if(newlist.is_normal(_ngbs[i])){
            //                newlist.get_normal(_ngbs[i])->complete(newlist,
            //                box, periodic);
            //                // we cannot check this ngb now, since this can
            //                cause spurious
            //                // changes in the ngb list of the current cell
            //                to_check.push_back(_ngbs[i]);
            //            }
        }
    }
    // order is important! (removing the first element first will shift
    // the other elements we have to remove...)
    //    sort(to_remove.begin(), to_remove.end());
    //    for(unsigned int i = to_remove.size(); i--;){
    for(set<unsigned int>::reverse_iterator it = to_remove.rbegin();
        it != to_remove.rend(); it++) {
        //        _ngbs.erase(_ngbs.begin()+*it);
        //        _faces.erase(_faces.begin()+*it);
        //        _exngbs.erase(_exngbs.end()-1);
        //        if(periodic){
        //            _ngbwalls.erase(_ngbwalls.begin()+*it);
        //        }
        _newngbs.erase(_newngbs.begin() + *it);
        _faces.erase(_faces.begin() + *it);
    }
    //    for(unsigned int i = 0; i < to_check.size(); i++){
    //        newlist.get_normal(to_check[i])->detect_crossovers(to_check[i],
    //    newlist, insertions, box, periodic);
    //    }
    return to_remove.size() > 0;
}

// void AdaptiveVorCell2d::save_restart(ostream& stream){
//    unsigned int vsize = _ngbs.size();
//    stream.write((char*) &vsize, sizeof(unsigned int));
//    stream.write((char*) &_ngbs[0], vsize*sizeof(int));
//    vsize = _orientation.size();
//    stream.write((char*) &vsize, sizeof(unsigned int));
//    stream.write((char*) &_orientation[0], vsize*sizeof(double));
//    vsize = _exngbs.size();
//    stream.write((char*) &vsize, sizeof(unsigned int));
//    stream.write((char*) &_exngbs[0], vsize*sizeof(int));
//    vsize = _faces.size();
//    stream.write((char*) &vsize, sizeof(unsigned int));
//    for(unsigned int i = 0; i < vsize; i++){
//        if(_faces[i]){
//            char flag = 1;
//            stream.write(&flag, 1);
//            _faces[i]->save_restart(stream);
//        } else {
//            char flag = 0;
//            stream.write(&flag, 1);
//        }
//    }
//    stream.write((char*) _pos, 2*sizeof(double));
//    stream.write((char*) &_oid, sizeof(unsigned int));
//    char wall = _wall;
//    stream.write(&wall, 1);
//}

/// @cond
#ifdef DITWOORD
VorCell::VorCell(istream& stream) {
    unsigned int vsize;
    stream.read((char*)&vsize, sizeof(unsigned int));
    _ngbs.resize(vsize);
    stream.read((char*)&_ngbs[0], vsize * sizeof(int));
    stream.read((char*)&vsize, sizeof(unsigned int));
    _orientation.resize(vsize);
    stream.read((char*)&_orientation[0], vsize * sizeof(double));
    stream.read((char*)&vsize, sizeof(unsigned int));
    _exngbs.resize(vsize);
    stream.read((char*)&_exngbs[0], vsize * sizeof(unsigned int));
    stream.read((char*)&vsize, sizeof(unsigned int));
    _faces.resize(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        char flag;
        stream.read(&flag, 1);
        if(flag > 0) {
            _faces[i] = new VorFace(stream);
        } else {
            _faces[i] = NULL;
        }
    }
    stream.read((char*)_pos, 2 * sizeof(double));
    stream.read((char*)&_oid, sizeof(unsigned int));
    char wall;
    stream.read(&wall, 1);
    _wall = wall > 0;
}
#endif
/// @endcond

/**
 * @brief Check if the given index is in the neighbour list
 *
 * @param id Unsigned integer index of a cell in the AdaptiveVorTess
 * @return True if id is a neighbour of the cell, false otherwise
 */
bool AdaptiveVorCell2d::is_ngb(unsigned int id) {
    unsigned int i = 0;
    while(i < _ngbs.size() && _ngbs[i] != (int)id) {
        i++;
    }
    return i < _ngbs.size();
}

/// @cond
#ifdef DITWOORD
bool VorCell::check_ngbs(unsigned int id, vector<VorCell*>& cells) {
    bool check = true;
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        if(_ngbs[i] >= cells.size()) {
            continue;
        }
        if(!cells[_ngbs[i]]->is_ngb(id)) {
            cerr << _ngbs[i] << " does not contain " << id << endl;
            check = false;
        }
    }
    return check;
}
#endif
/// @endcond

/**
 * @brief Calculate the 2D volume (area) of the cell
 *
 * The total area is the sum of the areas of the triangles constituted by the
 * cell generator and the two vertices of its faces.
 *
 * @return The area of the cell
 */
double AdaptiveVorCell2d::get_volume() {
    double V = 0.;
    for(unsigned int i = 0; i < _faces.size(); i++) {
        double L[2], R[2];
        _faces[i]->get_points(L, R);
        L[0] -= _pos[0];
        L[1] -= _pos[1];
        R[0] -= _pos[0];
        R[1] -= _pos[1];
        V += 0.5 * fabs(L[0] * R[1] - L[1] * R[0]);
    }
    return V;
}

/**
 * @brief Get the characteristic radius of the cell
 *
 * The characterstic radius is the radius of a disc with the same area as the
 * cell.
 *
 * @return The characteristic radius of the cell
 */
double AdaptiveVorCell2d::get_h() {
    return sqrt(get_volume() / M_PI);
}

/**
 * @brief Get the velocity of the cell, based on the primitive quantities of its
 * associated GasParticle
 *
 * @param particle GasParticle associated with this cell
 * @return Velocity of the cell in a co-moving mesh scheme
 */
Vec AdaptiveVorCell2d::get_velocity(GasParticle* particle) {
    StateVector W = particle->get_Wvec();
    // Springel correction
    Vec vel;
    if(particle->get_Wvec().rho()) {
        double centroidvec[2];
        get_centroid(centroidvec);
#if ndim_ == 2
        Vec centroid(centroidvec[0], centroidvec[1]);
#else
        Vec centroid;
#endif
        Vec d = centroid - particle->get_position();
        double Ri = get_volume();
        double csnd = particle->get_soundspeed();
        Ri = sqrt(Ri / M_PI);
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
    //    sqrt(_faces[i]->get_area()/M_PI)/
    //    (_central_point->get_position()-_faces[i]->get_midpoint()).norm());
    //        }
    //        if(maxfaceangle > 1.68){
    //            Vec d = get_centroid() - _central_point->get_position();
    //            double csnd = ((GasParticle*)_central_point->get_particle())
    //    ->get_soundspeed();
    //            vel = csnd/d.norm()*d;
    //        }
    //    }
    //    Vec vel = _central_point->get_particle()->get_mesh_v();
    double cellv[2];
    cellv[0] = W[1] + vel[0];
    cellv[1] = W[2] + vel[1];
    if(fabs(cellv[0]) < 1.e-10) {
        cellv[0] = 0.;
    }
    if(fabs(cellv[1]) < 1.e-10) {
        cellv[1] = 0.;
    }
    return Vec(cellv[0], cellv[1]);
}

/**
 * @brief Estimate gradients for the primitive quantities of this cell based on
 * the positions of its faces and the primitive quantities of its neighbours
 *
 * @param delta Array to store the results in
 * @param particle GasParticle associated with this cell
 * @param cells Reference to the cell list of the AdaptiveVorTess
 * @param box Cuboid specifying dimensions of the simulation box
 * @param periodic Flag signaling if we deal with a periodic (true) or a
 * reflective (false) simulation box
 */
void AdaptiveVorCell2d::estimate_gradients(
        StateVector* delta, GasParticle* particle,
        AdaptiveCellList<AdaptiveVorCell2d>& cells, Cuboid& box,
        bool periodic) {
    StateVector W = particle->get_Wvec();
    // gradient estimation
    StateVector Wmaxvec, Wminvec;
    Wmaxvec = W;
    Wminvec = W;
    for(unsigned int i = 0; i < _faces.size(); i++) {
        double Aij = _faces[i]->get_area();
        if(Aij) {
            double midfacevec[2];
            _faces[i]->get_midpoint(midfacevec);
            AdaptiveVorCell2d* ngb;
            double ngbposvec[2];
            if(cells.is_normal(_ngbs[i])) {
                ngb = cells.get_normal(_ngbs[i]);
                ngbposvec[0] = ngb->_pos[0];
                ngbposvec[1] = ngb->_pos[1];
                if(periodic && _ngbwalls[i] < 0) {
                    AdaptiveMeshUtils::get_periodic_position(ngbposvec,
                                                             _ngbwalls[i], box);
                }
            } else {
                ngb = cells.get_ghost(_ngbs[i]);
                AdaptiveMeshUtils::get_wall_position(ngb->_pos, ngbposvec,
                                                     ngb->_ngbs[0], box);
            }
#if ndim_ == 2
            Vec midface(midfacevec[0], midfacevec[1]);
            Vec ngbpos(ngbposvec[0], ngbposvec[1]);
#else
            Vec midface;
            Vec ngbpos;
#endif
            Vec c = midface - 0.5 * (particle->get_position() + ngbpos);
            Vec rij = particle->get_position() - ngbpos;
            double rnorm = rij.norm();
            StateVector Wj;
            // check if we have to mirror this particle
            if(cells.is_ghost(_ngbs[i])) {
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
    double centroidvec[2];
    get_centroid(centroidvec);
#if ndim_ == 2
    Vec centroid(centroidvec[0], centroidvec[1]);
#else
    Vec centroid;
#endif
    StateVector alphavec(1.);
    for(unsigned int l = _faces.size(); l--;) {
        if(_faces[l]->get_area()) {
            double midfacevec[2];
            _faces[l]->get_midpoint(midfacevec);
            Vec midface(midfacevec[0], midfacevec[1]);
            Vec d = midface - centroid;
            StateVector deltap = delta[0] * d[0] + delta[1] * d[1];
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
    //    bool check = true;
    //    unsigned int numloop = 0;
    //    while(check && numloop < 100){
    //        numloop ++;
    //        check = false;
    //        for(unsigned int l = _faces.size(); l--;){
    //            if(_faces[l]->get_area()){
    //                double midfacevec[2];
    //                _faces[l]->get_midpoint(midfacevec);
    //                Vec midface(midfacevec[0], midfacevec[1]);
    //                Vec d = midface - centroid;
    //                StateVector deltap = W + delta[0]*d[0] + delta[1]*d[1];
    //                if(deltap[0]-Wminvec[0] < 0. ||
    //    deltap.p() - Wminvec.p() < 0.){
    //                    delta[0] *= 0.99;
    //                    delta[1] *= 0.99;
    //                    check = true;
    //                    break;
    //                }
    //            }
    //        }
    //    }
    //    if(numloop == 100){
    //        cerr << "Error! Maximum number of iterations (100) reached in"
    //    "gradient estimation!" << endl;
    //        exit(32);
    //    }
}

/**
 * @brief Estimate gradients for the primitive quantities of this cell based on
 * the positions of its faces and the primitive quantities of its neighbours
 *
 * @param delta Array to store the results in
 * @param particle GasParticle associated with this cell
 * @param positionlist Reference to the cell list of the AdaptiveVorTess
 * @param periodic Flag signaling if we deal with a periodic (true) or a
 * reflective (false) simulation box
 */
void AdaptiveVorCell2d::estimate_gradients(StateVector* delta,
                                           GasParticle* particle,
                                           NewCellList& positionlist,
                                           bool periodic) {
    StateVector W = particle->get_Wvec();
    // gradient estimation
    StateVector Wmaxvec, Wminvec;
    Wmaxvec = W;
    Wminvec = W;
    for(unsigned int i = 0; i < _faces.size(); i++) {
        double Aij = _faces[i]->get_area();
        if(Aij) {
            double midfacevec[2];
            _faces[i]->get_midpoint(midfacevec);
            //            AdaptiveVorCell2d* ngb;
            //            double ngbposvec[2];
            //            if(cells.is_normal(_ngbs[i])){
            //                ngb = cells.get_normal(_ngbs[i]);
            //                ngbposvec[0] = ngb->_pos[0];
            //                ngbposvec[1] = ngb->_pos[1];
            //                if(periodic && _ngbwalls[i] < 0){
            //                    AdaptiveMeshUtils::get_periodic_position(ngbposvec,
            //                                                             _ngbwalls[i],
            //                                                             box);
            //                }
            //            } else {
            //                ngb = cells.get_ghost(_ngbs[i]);
            //                AdaptiveMeshUtils::get_wall_position(ngb->_pos,
            //                ngbposvec,
            //                                                     ngb->_ngbs[0],
            //                                                     box);
            //            }
            Vec ngbpos = positionlist.get_position(_ngbs[i]);
#if ndim_ == 2
            Vec midface(midfacevec[0], midfacevec[1]);
//            Vec ngbpos(ngbposvec[0], ngbposvec[1]);
#else
            Vec midface;
            Vec ngbpos;
#endif
            Vec c = midface - 0.5 * (particle->get_position() + ngbpos);
            Vec rij = particle->get_position() - ngbpos;
            double rnorm = rij.norm();
            StateVector Wj;
            // check if we have to mirror this particle
            // only periodic boundaries for now!
            //            if(cells.is_ghost(_ngbs[i])){
            //                Wj = particle->get_Wvec();
            //                _faces[i]->transform(Wj, _pos, ngbposvec);
            //                Wj[1] = -Wj[1];
            //                _faces[i]->invtransform(Wj, _pos, ngbposvec);
            //            } else {
            //                Wj = ngb->_particle->get_Wvec();
            Wj = positionlist.get_W(_ngbs[i]);
            //            }
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
    double centroidvec[2];
    get_centroid(centroidvec);
#if ndim_ == 2
    Vec centroid(centroidvec[0], centroidvec[1]);
#else
    Vec centroid;
#endif
    StateVector alphavec(1.);
    for(unsigned int l = _faces.size(); l--;) {
        if(_faces[l]->get_area()) {
            double midfacevec[2];
            _faces[l]->get_midpoint(midfacevec);
            Vec midface(midfacevec[0], midfacevec[1]);
            Vec d = midface - centroid;
            StateVector deltap = delta[0] * d[0] + delta[1] * d[1];
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
    //    bool check = true;
    //    unsigned int numloop = 0;
    //    while(check && numloop < 100){
    //        numloop ++;
    //        check = false;
    //        for(unsigned int l = _faces.size(); l--;){
    //            if(_faces[l]->get_area()){
    //                double midfacevec[2];
    //                _faces[l]->get_midpoint(midfacevec);
    //                Vec midface(midfacevec[0], midfacevec[1]);
    //                Vec d = midface - centroid;
    //                StateVector deltap = W + delta[0]*d[0] + delta[1]*d[1];
    //                if(deltap[0]-Wminvec[0] < 0. ||
    //    deltap.p() - Wminvec.p() < 0.){
    //                    delta[0] *= 0.99;
    //                    delta[1] *= 0.99;
    //                    check = true;
    //                    break;
    //                }
    //            }
    //        }
    //    }
    //    if(numloop == 100){
    //        cerr << "Error! Maximum number of iterations (100) reached in"
    //    "gradient estimation!" << endl;
    //        exit(32);
    //    }
}

/**
 * @brief Get the GasParticle associated with this cell
 *
 * @return The GasParticle associated with this cell
 */
GasParticle* AdaptiveVorCell2d::get_particle() {
    return _particle;
}

/**
 * @brief Get the AdaptiveVorFace2d with the given index in the internal list
 *
 * @param i Unsigned integer index of the AdaptiveVorFace2d in the internal list
 * @return AdaptiveVorFace2d with the given index
 */
AdaptiveVorFace2d* AdaptiveVorCell2d::get_face(unsigned int i) {
    return _faces[i];
}

/**
 * @brief Get the position of the ghost represented by this cell
 *
 * We always store the actual position of the cell generator, so we need this
 * method to retrieve the position of its ghost if we need it.
 *
 * @param ngbposvec Array to store the result in
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell2d::get_ghost_position(double* ngbposvec, Cuboid& box) {
    AdaptiveMeshUtils::get_wall_position(_pos, ngbposvec, _ngbs[0], box);
}

/**
 * @brief Check if the generator of this cell still lies inside the periodic
 * simulation box and correct its position and its internal state if not
 *
 * The correction consists of
 *  -# moving the generator to a position within the box by adding or
 *     subtracting the box length to one or multiple coordinates
 *  -# changing the boundary flags of all its neighbours
 *  -# changing the boundary flag of the cell in all of its neighbours
 *  -# correcting the cell centroid
 *
 * @param id Index of the cell in the AdaptiveVorTess cell list
 * @param cells Reference to the AdaptiveVorTess cell list
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell2d::apply_periodic_boundaries(
        unsigned int id, AdaptiveCellList<AdaptiveVorCell2d>& cells,
        Cuboid& box) {
    int wallpos = AdaptiveMeshUtils::trim_wallpos(_pos, box);
    if(wallpos < 0) {
        // the cell was translated, we have to update the neighbours
        for(unsigned int i = 0; i < _ngbs.size(); i++) {
            if(!cells.is_orphan(_ngbs[i])) {
                if(wallpos == -4) {
                    vector<int> oldwallpos =
                            cells.get_normal(_ngbs[i])->get_ngbwalls();
                    cerr << "oldwallpos:" << endl;
                    for(unsigned int j = 0; j < oldwallpos.size(); j++) {
                        cerr << j << ": " << oldwallpos[j] << endl;
                    }
                }
                cells.get_normal(_ngbs[i])->swap_wallpos(id, wallpos);
                if(wallpos == -4) {
                    vector<int> oldwallpos =
                            cells.get_normal(_ngbs[i])->get_ngbwalls();
                    cerr << "newwallpos:" << endl;
                    for(unsigned int j = 0; j < oldwallpos.size(); j++) {
                        cerr << j << ": " << oldwallpos[j] << endl;
                    }
                }
                _ngbwalls[i] =
                        AdaptiveMeshUtils::swap_wallpos(_ngbwalls[i], wallpos);
            } else {
                int ngbwall = _ngbwalls[i];
                ngbwall = AdaptiveMeshUtils::swap_wallpos(ngbwall, wallpos);
                _ngbwalls[i] = ngbwall;
            }
        }
        // we also trim the centroid, in case this cell is inactive
        // otherwise, it could happen that the wrong centroid is used in
        // AdaptiveFace.
        if(cells.is_normal(id)) {
            double d[2];
            d[0] = _particle->x() - _pos[0];
            d[1] = _particle->y() - _pos[1];
            _particle->set_x(_pos[0]);
            _particle->set_y(_pos[1]);
            Vec centroid = _particle->get_centroid();
            centroid[0] -= d[0];
            centroid[1] -= d[1];
            _particle->set_centroid(centroid);
        }
    }
}

/**
 * @brief Change the boundary flag for the given neighbour due to the periodic
 * movement specified by the given boundary flag
 *
 * If the given neighbour is not found in the list, the code will abort.
 *
 * @param id Unsigned integer index of a valid neighbour in the list
 * @param wallpos Boundary flag specifying the periodic movement of the
 * neighbour
 */
void AdaptiveVorCell2d::swap_wallpos(unsigned int id, int wallpos) {
    vector<int>::iterator it = _ngbs.begin();
    while(it != _ngbs.end() && *it != (int)id) {
        it++;
    }
    if(it == _ngbs.end()) {
        cerr << "NgbException(inf)" << endl;
        my_exit();
    }
    unsigned int i = it - _ngbs.begin();
    _ngbwalls[i] = AdaptiveMeshUtils::replace_wallpos(_ngbwalls[i], wallpos);
}

/**
 * @brief Get the position of the periodic copy of the given neighbour
 *
 * @param index Unsigned integer index of a valid neighbour in the list
 * @param cells Reference to the AdaptiveVorTess cell list
 * @param pos Array to store the result in
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveVorCell2d::get_periodic_position(
        unsigned int index, AdaptiveCellList<AdaptiveVorCell2d>& cells,
        double* pos, Cuboid& box) {
    AdaptiveVorCell2d* ngb = cells.get_normal(_ngbs[index]);
    pos[0] = ngb->_pos[0];
    pos[1] = ngb->_pos[1];
    if(_ngbwalls[index] < 0) {
        AdaptiveMeshUtils::get_periodic_position(pos, _ngbwalls[index], box);
    }
}

/**
 * @brief Flag this cell to be checked during the third detection iteration
 */
void AdaptiveVorCell2d::queue() {
    _flag2 = true;
}

/**
 * @brief Check if this cell has to be tested during the third detection
 * iteration
 *
 * @return True if the cell has to be tested, false otherwise
 */
bool AdaptiveVorCell2d::queued() {
    return _flag2;
}

/**
 * @brief Flag this cell to be tested during the first detection iteration
 */
void AdaptiveVorCell2d::activate() {
    _active = true;
}

/**
 * @brief Undo the action of AdaptiveVorCell2d::activate().
 */
void AdaptiveVorCell2d::deactivate() {
    _active = false;
}

/**
 * @brief Check if this cell has to be tested during the first detection
 * iteration
 *
 * @return True if the cell has to be tested, false otherwise
 */
bool AdaptiveVorCell2d::active() {
    return _active;
}

/**
 * @brief Flag this cell to be tested during the second detection iteration
 */
void AdaptiveVorCell2d::semiactivate() {
    _flag1 = true;
}

/**
 * @brief Check whether this cell has to be tested during the second detection
 * iteration
 *
 * @return True if the cell has to tested, false otherwise
 */
bool AdaptiveVorCell2d::semiactive() {
    return _flag1;
}

/**
 * @brief Reset all detection iteration flags
 */
void AdaptiveVorCell2d::reset_flags() {
    _flag1 = false;
    _flag2 = false;
    _active = false;
}

/**
 * @brief Update the internal neighbour list with the given new indices of the
 * cells in the AdaptiveVorTess
 *
 * @param new_ids std::vector of unsigned integer new indices
 */
void AdaptiveVorCell2d::update_ngbs(vector<unsigned int>& new_ids) {
    //    _oid = new_ids[_oid];
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        if(_ngbs[i] < (int)new_ids.size() && _ngbs[i] >= 0) {
            _ngbs[i] = new_ids[_ngbs[i]];
        }
    }
}

/**
 * @brief Update the original ID of the cell this cell is an orphan of using the
 * given new indices
 *
 * @param new_ids std::vector of unsigned integer new indices
 */
void AdaptiveVorCell2d::update_local_id(vector<unsigned int>& new_ids) {
    _oid = new_ids[_oid];
}

/**
 * @brief Set the original ID of the cell this cell is an orphan of to the new
 * given value
 *
 * @param new_id New unsigned integer index
 */
void AdaptiveVorCell2d::update_local_id(unsigned int new_id) {
    _oid = new_id;
}

/**
 * @brief Update the internal neighbour list of orphans with the given new
 * indices
 *
 * @param new_ids std::vector of unsigned integer new indices
 * @param cells Reference to the AdaptiveVorTess cell list
 */
void AdaptiveVorCell2d::update_orphans(
        vector<unsigned int>& new_ids,
        AdaptiveCellList<AdaptiveVorCell2d>& cells) {
    for(unsigned int i = 0; i < _ngbs.size(); i++) {
        if(cells.is_orphan(_ngbs[i])) {
            _ngbs[i] = new_ids[-_ngbs[i]];
        }
    }
}

/**
 * @brief Get the first neighbour in the internal list
 *
 * This is the only neighbour if the cell is a ghost.
 *
 * @return The first neighbour in the list
 */
int AdaptiveVorCell2d::get_first_ngb() {
    return _ngbs[0];
}

/**
 * @brief Set the ID of the cell
 *
 * @param id New value for the ID of the cell
 */
void AdaptiveVorCell2d::set_id(unsigned long id) {
    _id = id;
}

/**
 * @brief Get the ID of the cell
 *
 * @return The ID of the cell
 */
unsigned long AdaptiveVorCell2d::get_id() {
    return _id;
}

/**
 * @brief Set the rank of the cell
 *
 * The rank is equal to the rank of the current process for all normal and ghost
 * cells and is equal to the rank of the process that holds the original cell
 * for orphans.
 *
 * @param rank Integer rank of the MPI process on which the original cell
 * resides
 */
void AdaptiveVorCell2d::set_rank(int rank) {
    _rank = rank;
}

/**
 * @brief Get the rank of the cell
 *
 * The rank is equal to the rank of the current process for all normal and ghost
 * cells and is equal to the rank of the process that holds the original cell
 * for orphans.
 *
 * @return Integer rank of the MPI process on which the original cell resides
 */
int AdaptiveVorCell2d::get_rank() {
    return _rank;
}

/**
 * @brief Set the ID of the particle corresponding to this cell
 *
 * @param particleID ID of the particle that corresponds to this cell
 */
void AdaptiveVorCell2d::set_particleID(unsigned long particleID) {
    _particleID = particleID;
}

/**
 * @brief Get the ID of the particle corresponding to this cell
 *
 * @return particleID ID of the particle that corresponds to this cell
 */
unsigned long AdaptiveVorCell2d::get_particleID() {
    return _particleID;
}

/**
 * @brief Get a reference to the new ngb list
 *
 * @return Reference to the new ngb list
 */
vector<unsigned long>& AdaptiveVorCell2d::get_newngbs() {
    return _newngbs;
}

/**
 * @brief Add a neighbour by its unique ID
 *
 * @param particleID Unique ID of a particle or copy of a particle
 */
void AdaptiveVorCell2d::add_ngb_particleID(unsigned long particleID) {
    _newngbs.push_back(particleID);
}

/**
 * @brief Dump the cell to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void AdaptiveVorCell2d::dump(RestartFile& rfile) {
    rfile.write(_ngbs);
    rfile.write(_ngbwalls);
    rfile.write(_orientation);
    rfile.write(_exngbs);
    unsigned int vsize = _faces.size();
    rfile.write(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        _faces[i]->dump(rfile);
    }
    rfile.write(_pos, 2);
    rfile.write(_oid);
    rfile.write(_wall);
    rfile.write(_flag1);
    rfile.write(_flag2);
    rfile.write(_active);
}

/**
 * @brief Restart constructor. Initialize the cell from the given RestartFile.
 *
 * @param rfile RestartFile to read from
 * @param particle GasParticle associated with this cell
 */
AdaptiveVorCell2d::AdaptiveVorCell2d(RestartFile& rfile,
                                     GasParticle* particle) {
    rfile.read(_ngbs);
    rfile.read(_ngbwalls);
    rfile.read(_orientation);
    rfile.read(_exngbs);
    unsigned int vsize;
    rfile.read(vsize);
    _faces.resize(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        _faces[i] = new AdaptiveVorFace2d(rfile);
    }
    rfile.read(_pos, 2);
    rfile.read(_oid);
    rfile.read(_wall);
    rfile.read(_flag1);
    rfile.read(_flag2);
    rfile.read(_active);
    _particle = particle;
}

/**
 * @brief MPI packer
 *
 * @param buffer Buffer to write to
 * @param bufsize Size of the buffer
 * @param position Current position in the buffer (is updated)
 * @param cells Reference to the cell list
 */
void AdaptiveVorCell2d::pack_data(void* buffer, int bufsize, int* position,
                                  AdaptiveCellList<AdaptiveVorCell2d>& cells) {
    Hilbert_Object::pack_data(buffer, bufsize, position);
    unsigned int vsize = _ngbs.size();
    MyMPI_Pack(&vsize, 1, MPI_UNSIGNED, buffer, bufsize, position);
    for(unsigned int i = 0; i < vsize; i++) {
        int rank = MPIGlobal::rank;
        if(cells.is_orphan(_ngbs[i])) {
            rank = cells.get_orphan(_ngbs[i])->get_rank();
        }
        unsigned int oid = cells.get_cell(_ngbs[i])->get_original();
        int wall = _ngbwalls[i];
        MyMPI_Pack(&rank, 1, MPI_INT, buffer, bufsize, position);
        MyMPI_Pack(&oid, 1, MPI_UNSIGNED, buffer, bufsize, position);
        MyMPI_Pack(&wall, 1, MPI_INT, buffer, bufsize, position);
    }

    MyMPI_Pack(_pos, 2, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_oid, 1, MPI_UNSIGNED, buffer, bufsize, position);
    char boolean = _wall;
    MyMPI_Pack(&boolean, 1, MPI_BYTE, buffer, bufsize, position);
    vsize = _newngbs.size();
    MyMPI_Pack(&vsize, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_newngbs[0], vsize, MPI_UNSIGNED_LONG, buffer, bufsize,
               position);
}

/**
 * @brief MPI constructor
 *
 * @param buffer Buffer to read from
 * @param bufsize Size of the buffer
 * @param position Current position in the buffer (is updated)
 * @param cells Reference to the cell list
 */
AdaptiveVorCell2d::AdaptiveVorCell2d(void* buffer, int bufsize, int* position,
                                     AdaptiveCellList<AdaptiveVorCell2d>& cells)
        : Hilbert_Object(buffer, bufsize, position) {
    unsigned int vsize;
    MyMPI_Unpack(buffer, bufsize, position, &vsize, 1, MPI_UNSIGNED);
    _ngbs.resize(vsize);
    _ngbwalls.resize(vsize);
    _faces.resize(vsize, NULL);
    for(unsigned int i = 0; i < vsize; i++) {
        int rank;
        unsigned int oid;
        int wall;
        MyMPI_Unpack(buffer, bufsize, position, &rank, 1, MPI_INT);
        MyMPI_Unpack(buffer, bufsize, position, &oid, 1, MPI_UNSIGNED);
        MyMPI_Unpack(buffer, bufsize, position, &wall, 1, MPI_INT);
        if(rank == MPIGlobal::rank) {
            _ngbs[i] = oid;
            _ngbwalls[i] = wall;
        } else {
            int orphanid;
            if((orphanid = cells.get_orphan(rank, oid))) {
                _ngbs[i] = orphanid;
                _ngbwalls[i] = wall;
            } else {
                cerr << "Creating a new orphan is not yet supported!" << endl;
                my_exit();
            }
        }
    }

    MyMPI_Unpack(buffer, bufsize, position, _pos, 2, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_oid, 1, MPI_UNSIGNED);
    char boolean;
    MyMPI_Unpack(buffer, bufsize, position, &boolean, 1, MPI_BYTE);
    _wall = boolean > 0;

    // should be set to the correct value!
    _particle = NULL;
    _id = 0;

    // correct?
    _flag1 = false;
    _flag2 = false;
    _active = false;

    _rank = MPIGlobal::rank;  // or -1?

    MyMPI_Unpack(buffer, bufsize, position, &vsize, 1, MPI_UNSIGNED);
    _newngbs.resize(vsize);
    MyMPI_Unpack(buffer, bufsize, position, &_newngbs[0], vsize,
                 MPI_UNSIGNED_LONG);
}

/**
 * @brief MPI packer
 *
 * @param buffer Buffer to write to
 * @param bufsize Size of the buffer
 * @param position Current position in the buffer (is updated)
 */
void AdaptiveVorCell2d::pack_data(void* buffer, int bufsize, int* position) {
    Hilbert_Object::pack_data(buffer, bufsize, position);
    MyMPI_Pack(&_particleID, 1, MPI_UNSIGNED_LONG, buffer, bufsize, position);
    unsigned int vsize = _newngbs.size();
    MyMPI_Pack(&vsize, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_newngbs[0], vsize, MPI_UNSIGNED_LONG, buffer, bufsize,
               position);
}

/**
 * @brief MPI constructor
 *
 * @param buffer Buffer to read from
 * @param bufsize Size of the buffer
 * @param position Current position in the buffer (is updated)
 */
AdaptiveVorCell2d::AdaptiveVorCell2d(void* buffer, int bufsize, int* position)
        : Hilbert_Object(buffer, bufsize, position) {
    MyMPI_Unpack(buffer, bufsize, position, &_particleID, 1, MPI_UNSIGNED_LONG);
    unsigned int vsize;
    MyMPI_Unpack(buffer, bufsize, position, &vsize, 1, MPI_UNSIGNED);
    _newngbs.resize(vsize);
    MyMPI_Unpack(buffer, bufsize, position, &_newngbs[0], vsize,
                 MPI_UNSIGNED_LONG);
    _faces.resize(vsize, NULL);
}
