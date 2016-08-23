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
 * @file AdaptiveVorTess3d.cpp
 *
 * @brief 3D evolving mesh: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdaptiveVorTess3d.hpp"
#include "AdaptiveFace3d.hpp"            // for AdaptiveFace3d
#include "AdaptiveMeshException.hpp"     // for AdaptiveMeshException
#include "AdaptiveMeshUtils.hpp"         // for get_wall
#include "AdaptiveVorCell3d.hpp"         // for AdaptiveVorCell3d
#include "AdaptiveVorFace3d.hpp"         // for AdaptiveVorFace3d
#include "Error.hpp"                     // for my_exit
#include "RestartFile.hpp"               // for RestartFile
#include "VorCell.hpp"                   // for VorCell
#include "VorFace.hpp"                   // for VorFace
#include "VorGen.hpp"                    // for VorGen
#include "VorTess.hpp"                   // for VorTess
#include "predicates.hpp"                // for orient3d_old
#include "utilities/GasParticle.hpp"     // for GasParticle
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include <cstddef>                       // for NULL
#include <ext/alloc_traits.h>
#include <iostream>  // for cerr, cout
using namespace std;

//#define FULLCHECK
//#define VERBOSE
//#define CHECK_EVERY_STEP 0

// unsigned int loop = 0;

/**
 * @brief Constructor
 *
 * Initialize an evolving mesh based on the given valid VorTess.
 *
 * @param tesselation Valid VorTess constructed during another algorithm
 * @param cuboid Cuboid specifying the dimensions of the simulation box
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
AdaptiveVorTess3d::AdaptiveVorTess3d(VorTess* tesselation, Cuboid cuboid,
                                     bool periodic)
        : _periodic(periodic), _cuboid(cuboid) {
    vector<VorCell*> cells = tesselation->get_cells();
    unsigned int cellsize = cells.size();
    _cells.resize(cellsize);
    for(unsigned int i = 0; i < cellsize; i++) {
        double pos[3];
        Vec position = cells[i]->get_particle()->get_position();
        pos[0] = position[0];
        pos[1] = position[1];
        pos[2] = position[2];
        unsigned int id = cells[i]->get_particle()->get_local_id();
        _cells[id] = new AdaptiveVorCell3d(pos, id, cells[i]->get_particle());
        vector<VorGen*> ngbs = cells[i]->get_ngbs();
        vector<VorFace*> faces = cells[i]->get_faces();
        for(unsigned int j = 0; j < ngbs.size(); j++) {
            //            if(faces[j]->get_area()){
            if(ngbs[j]->get_particle()) {
                VorGen* ngb = ngbs[j];
                int wall = 0;
                if((ngb->get_position() - ngb->get_particle()->get_position())
                           .norm2()) {
                    double ngbpos[3];
                    ngbpos[0] = ngb->x();
                    ngbpos[1] = ngb->y();
                    ngbpos[2] = ngb->z();
                    double opos[3];
                    opos[0] = ngb->get_particle()->x();
                    opos[1] = ngb->get_particle()->y();
                    opos[2] = ngb->get_particle()->z();
                    wall = AdaptiveMeshUtils::get_wall(ngbpos, opos);
                }
                _cells[id]->add_ngb(ngb->get_particle()->get_local_id(), wall);
                double b[3], c[3], d[3];
                vector<VorGen*> ngbfaces = faces[j]->get_facengbs();
                b[0] = ngb->x();
                b[1] = ngb->y();
                b[2] = ngb->z();
                c[0] = ngbfaces[0]->x();
                c[1] = ngbfaces[0]->y();
                c[2] = ngbfaces[0]->z();
                d[0] = ngbfaces[1]->x();
                d[1] = ngbfaces[1]->y();
                d[2] = ngbfaces[1]->z();
                bool reverse = predicates::orient3d_old(pos, b, c, d) < 0.;
                for(unsigned int k = 0; k < ngbfaces.size(); k++) {
                    unsigned int index = k;
                    if(reverse) {
                        index = ngbfaces.size() - index - 1;
                    }
                    if(ngbfaces[index]->get_particle()) {
                        wall = 0;
                        if((ngbfaces[index]->get_position() -
                            ngbfaces[index]->get_particle()->get_position())
                                   .norm2()) {
                            double ngbpos[3];
                            ngbpos[0] = ngbfaces[index]->x();
                            ngbpos[1] = ngbfaces[index]->y();
                            ngbpos[2] = ngbfaces[index]->z();
                            double opos[3];
                            opos[0] = ngbfaces[index]->get_particle()->x();
                            opos[1] = ngbfaces[index]->get_particle()->y();
                            opos[2] = ngbfaces[index]->get_particle()->z();
                            wall = AdaptiveMeshUtils::get_wall(ngbpos, opos);
                        }
                        _cells[id]->add_facengb(
                                ngbfaces[index]->get_particle()->get_local_id(),
                                wall);
                    } else {
                        double ngbfacepos[3];
                        ngbfacepos[0] = ngbfaces[index]->x();
                        ngbfacepos[1] = ngbfaces[index]->y();
                        ngbfacepos[2] = ngbfaces[index]->z();
                        _ghosts.push_back(new AdaptiveVorCell3d(
                                ngbfacepos, ngbfaces[index]->get_original(),
                                NULL));
                        _cells[id]->add_facengb(_ghosts.size() - 1 + cellsize,
                                                0);
                    }
                }
            } else {
                double b[3];
                b[0] = ngbs[j]->get_position().x();
                b[1] = ngbs[j]->get_position().y();
                b[2] = ngbs[j]->get_position().z();
                _ghosts.push_back(new AdaptiveVorCell3d(
                        b, ngbs[j]->get_original(), NULL));
                _cells[id]->add_ngb(_ghosts.size() - 1 + cellsize, 0);
                double c[3], d[3];
                vector<VorGen*> ngbfaces = faces[j]->get_facengbs();
                c[0] = ngbfaces[0]->x();
                c[1] = ngbfaces[0]->y();
                c[2] = ngbfaces[0]->z();
                d[0] = ngbfaces[1]->x();
                d[1] = ngbfaces[1]->y();
                d[2] = ngbfaces[1]->z();
                bool reverse = predicates::orient3d_old(pos, b, c, d) < 0.;
                for(unsigned int k = 0; k < ngbfaces.size(); k++) {
                    unsigned int index = k;
                    if(reverse) {
                        index = ngbfaces.size() - index - 1;
                    }
                    if(ngbfaces[index]->get_particle()) {
                        _cells[id]->add_facengb(
                                ngbfaces[index]->get_particle()->get_local_id(),
                                0);
                    } else {
                        double ngbfacepos[3];
                        ngbfacepos[0] = ngbfaces[index]->x();
                        ngbfacepos[1] = ngbfaces[index]->y();
                        ngbfacepos[2] = ngbfaces[index]->z();
                        _ghosts.push_back(new AdaptiveVorCell3d(
                                ngbfacepos, ngbfaces[index]->get_original(),
                                NULL));
                        _cells[id]->add_facengb(_ghosts.size() - 1 + cellsize,
                                                0);
                    }
                }
            }
            //            }
        }
        //        _cells[id]->clean_facengbs();
    }

    if(!_periodic) {
        vector<vector<unsigned int> > ghostids(cellsize);
        for(unsigned int i = 0; i < cellsize; i++) {
            ghostids[i].resize(26, 0);
        }

        vector<unsigned int> repids(_ghosts.size());
        for(unsigned int i = 0; i < _ghosts.size(); i++) {
            double pos[3];
            _ghosts[i]->get_position(pos);
            double opos[3];
            _cells[_ghosts[i]->get_id()]->get_position(opos);
            int id = AdaptiveMeshUtils::get_wall(pos, opos);
            _ghosts[i]->add_ngb(id, 0);
            _ghosts[i]->move(opos);
            if(!ghostids[_ghosts[i]->get_id()][id + 26]) {
                ghostids[_ghosts[i]->get_id()][id + 26] = i + cellsize;
            }
            repids[i] = ghostids[_ghosts[i]->get_id()][id + 26];
        }

        for(unsigned int i = 0; i < cellsize; i++) {
            _cells[i]->update_ghosts(repids, cellsize);
        }
    }

    for(unsigned int i = 0; i < _cells.size(); i++) {
        _cells[i]->complete(_cells, _ghosts, _cuboid);
    }

    _tolerance = 0.01;
}

/**
 * @brief Destructor
 *
 * Clean up cells and ghosts.
 */
AdaptiveVorTess3d::~AdaptiveVorTess3d() {
    for(unsigned int i = 0; i < _cells.size(); i++) {
        delete _cells[i];
    }
    for(unsigned int i = 0; i < _ghosts.size(); i++) {
        delete _ghosts[i];
    }
}

/**
 * @brief Print out the tesselation to the given stream in a format that can be
 * plotted using gnuplot
 *
 * Use the following gnuplot command:
\verbatim
splot "<FILENAME>" with linespoints
\endverbatim
 *
 * @param stream std::ostream to write to
 */
void AdaptiveVorTess3d::print_tesselation_gnuplot(ostream& stream) {
    for(unsigned int i = 0; i < _cells.size(); i++) {
        _cells[i]->complete(_cells, _ghosts, _cuboid);
        _cells[i]->print(stream, i);
    }
}

/**
 * @brief Check if the given combination of edge and face flips is valid
 *
 * @param flips Array with counts for different edge and face flips
 * @return True if the combination of flips corresponds to a single mesh
 * deviation and can be corrected using a single face insertion or removal,
 * false otherwise
 */
bool AdaptiveVorTess3d::valid_flips(unsigned int* flips) {
    if(flips[3]) {
        return false;
    }
    if(flips[1] > 3) {
        return false;
    }
    if(flips[2] > 2) {
        return false;
    }
    if(flips[4] > 2) {
        return false;
    }
    if(flips[2] && flips[4]) {
        return false;
    }
    //    if(flips[1] && flips[1] != 3){
    //        return false;
    //    }
    if(flips[2] && flips[2] != 2) {
        return false;
    }
    if(flips[4] && flips[4] != 2) {
        return false;
    }
    // if an edge flip is the only thing that happens, it is detected in
    // exactly 2 cells. However, it is possible that the second face is not
    // processed and hence only one edge flip is detected.
    //    if(!flips[4] && !flips[2] && flips[1] && flips[1] > 2){
    //        return false;
    //    }
    //    if(flips[4] && flips[1] > 2){
    //        return false;
    //    }
    return true;
}

/**
 * @brief Check if the given combination of edge and face flips is valid
 *
 * @param flips Array with counts for different edge and face flips
 * @return True if the combination of flips corresponds to a single mesh
 * deviation and can be corrected using a single face insertion or removal,
 * false otherwise
 */
bool AdaptiveVorTess3d::valid_flips_minimal(unsigned int* flips) {
    if(flips[3]) {
        return false;
    }
    if(flips[1] > 6) {
        return false;
    }
    if(flips[2] > 1) {
        return false;
    }
    if(flips[4] > 1) {
        return false;
    }
    if(flips[2] && flips[4]) {
        return false;
    }
    if(!flips[4] && !flips[2] && flips[1] && flips[1] > 3) {
        return false;
    }
    if(flips[4] && flips[1] > 4) {
        return false;
    }
    return true;
}

/**
 * @brief Update the tesselation so that is a valid Voronoi tesselation again
 *
 * Only the active particles and their direct neighbours are evolved, using the
 * provided value of currentTime to check for active particles.
 *
 * @param particles std::vector of GasParticle instances associated with the
 * cells of the tesselation
 * @param currentTime unsigned long integer current simulation time
 * @return The number of active cells that will have valid associated Voronoi
 * cells
 */
unsigned int AdaptiveVorTess3d::update(ParticleVector& particles,
                                       unsigned long currentTime) {
    vector<AdaptiveVorCell3d*> backup(_cells.size());
    vector<Vec> new_positions_backup(_new_positions.size());
    vector<unsigned int> new_ids(_cells.size());
    for(unsigned int i = 0; i < _cells.size(); i++) {
        backup[i] = _cells[i];
        new_positions_backup[i] = _new_positions[i];
        new_ids[i] = _cells[i]->get_particle()->get_local_id();
    }

    for(unsigned int i = 0; i < _cells.size(); i++) {
        _cells[new_ids[i]] = backup[i];
        _new_positions[new_ids[i]] = new_positions_backup[i];
        _cells[new_ids[i]]->update_ngbs(new_ids);
    }

    // flag and count all active cells
    // also flag their neighbours, based on the neighbours from the previous
    // step
    // NOTE: we should also activate new neighbours, this could happen later on
    // so maybe we should only flag the active particles and then flag all
    // their neighbours in a second step...
    unsigned int count = 0;
    for(unsigned int i = 0; i < _cells.size(); i++) {
        if(particles.gas(i)->get_endtime() == currentTime) {
            unsigned int index = particles.gas(i)->get_local_id();
            _cells[index]->activate();
            count++;
            //            vector<int> ngbs = _cells[index]->get_ngbs();
            //            for(unsigned int j = 0; j < ngbs.size(); j++){
            //                if(ngbs[j] < _cells.size()){
            //                    _cells[ngbs[j]]->queue();
            //                    vector<int> ngbs2 =
            //                    _cells[ngbs[j]]->get_ngbs();
            //                    for(unsigned int k = 0; k < ngbs2.size();
            //                    k++){
            //                        if(ngbs2[k] < _cells.size()){
            //                            _cells[ngbs2[k]]->queue();
            //                        }
            //                    }
            //                }
            //            }
        }
    }
    // for debugging purposes: activate all cells
    //    for(unsigned int i = 0; i < _cells.size(); i++){
    //        _cells[i]->queue();
    //    }

    // evolve the active cells one-by-one
    for(unsigned int mode = 0; mode < 3; mode++) {
        for(unsigned int i = 0; i < _cells.size(); i++) {
            // if the cell is not active, skip it
            if(mode) {
                if(mode == 1) {
                    if(!_cells[i]->semiactive() || _cells[i]->active()) {
                        continue;
                    }
                } else {
                    if(!_cells[i]->queued() || _cells[i]->active() ||
                       _cells[i]->semiactive()) {
                        continue;
                    }
                }
            } else {
                if(!_cells[i]->active()) {
                    continue;
                }
            }
//        if(loop == 176 && i == 11){
//            exit(1);
//        }
#ifdef VERBOSE
            cerr << "generator " << i << endl;
#endif
            // backup the last valid position
            double oldpos[3];
            _cells[i]->get_valid_pos(oldpos);
            double newpos[3];
            newpos[0] = _new_positions[i].x();
            newpos[1] = _new_positions[i].y();
            newpos[2] = _new_positions[i].z();
            // calculate the magnitude of the displacement
            // this is used to determine when a displacement can be considered
            // to
            // be "small"
            double dx[4];
            dx[0] = oldpos[0] - newpos[0];
            dx[1] = oldpos[1] - newpos[1];
            dx[2] = oldpos[2] - newpos[2];
            dx[3] = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            // generic size of the cell, defined as the radius of a sphere with
            // the same volume as the cell. Value is based on the last active
            // version of the cell
            double h = _cells[i]->get_particle()->h();
// move the generator
#ifdef VERBOSE
            cerr << oldpos[0] << "\t" << oldpos[1] << "\t" << oldpos[2] << endl;
            cerr << newpos[0] << "\t" << newpos[1] << "\t" << newpos[2] << endl;
#endif
            _cells[i]->move(newpos);
            //        ofstream lfile("last_mesh.dat");
            //        _cells[i]->complete(_cells, _ghosts, _cuboid);
            //        _cells[i]->print(lfile, i);
            //        lfile.close();

            // now calculate the number of edge flips caused by moving the
            // generator
            // edge flips can only occur in faces that are shared by the current
            // cell or that have the current cell as a face neighbour
            // we also keep track of which cell reported the last flip of every
            // type
            unsigned int lastflips[3];
            lastflips[0] = i;
            lastflips[1] = i;
            lastflips[2] = i;
            unsigned int flips[5];
            flips[0] = 0;
            flips[1] = 0;
            flips[2] = 0;
            flips[3] = 0;
            flips[4] = 0;
#ifdef FULLCHECK
            unsigned int numflips =
                    _cells[i]->get_flips(i, _cells, _cuboid, _periodic);
            flips[numflips]++;
#else
            unsigned int numflips = _cells[i]->get_flips_minimal(
                    flips, i, _cells, _ghosts, _cuboid, _periodic);
#endif
            vector<int> ngbs = _cells[i]->get_ngbs();
            for(unsigned int j = 0; j < ngbs.size(); j++) {
                if(ngbs[j] < (int)_cells.size()) {
#ifdef FULLCHECK
                    numflips = _cells[ngbs[j]]->get_flips(ngbs[j], i, _cells,
                                                          _cuboid, _periodic);
                    flips[numflips]++;
#else
                    numflips = _cells[ngbs[j]]->get_flips_minimal(
                            flips, ngbs[j], i, _cells, _ghosts, _cuboid,
                            _periodic);
#endif
                    if(numflips == 1) {
                        lastflips[0] = ngbs[j];
                    }
                    if(numflips == 2) {
                        lastflips[1] = ngbs[j];
                    }
                    if(numflips == 4) {
                        lastflips[2] = ngbs[j];
                    }
                }
            }
#ifdef VERBOSE
            cerr << "A" << endl;
            cerr << flips[0] << "\t" << flips[1] << "\t" << flips[2] << "\t"
                 << flips[3] << "\t" << flips[4] << endl;
#endif
            // the total number of flips tells us something about the deviations
            // that happened
            // we accept three possible deviations:
            //  * a face insertion: signalled by 3 cells where 2 edges flipped
            //  * a face removal: signalled by 2 cells where a face and 3 edges
            //                    flipped + 3 cells where 2 edges flipped
            //  * a face flip: signalled by 2 cells where a face of 4 vertices
            //                 and 2 edges flipped + 2 cells where 2 edges
            //                 flipped
            // A cell with more than the described flips is invalid, as are all
            // cases with more flips.
            // In these cases, we have to split up the cell movement in smaller
            // movements until we are left with a single deviation.
            if(flips[1] || flips[2] || flips[3] || flips[4]) {
// not ok. We need to restore some things
// if we have more than one deviation, we split up the movement
#ifdef FULLCHECK
                bool test = valid_flips(flips);
#else
                bool test = valid_flips_minimal(flips);
#endif
                unsigned int numloop = 0;
                while(!test || dx[3] > _tolerance * h * h) {
                    double temppos[3];
                    double factor = 0.5;
                    // we start by halving the displacement and further divide
                    // the
                    // displacement by factors of 2 until a single or no
                    // deviation
                    // is left and the total displacement is below a certain
                    // fraction of the generic cell size
                    while(!test || dx[3] > _tolerance * h * h) {
                        // calculate a new position in between newpos and oldpos
                        temppos[0] =
                                (1. - factor) * oldpos[0] + factor * newpos[0];
                        temppos[1] =
                                (1. - factor) * oldpos[1] + factor * newpos[1];
                        temppos[2] =
                                (1. - factor) * oldpos[2] + factor * newpos[2];
                        // move the cell
                        _cells[i]->move(temppos);
                        // calculate the magnitude of the new displacement
                        // NOTE: can't we use the old value and factor to do
                        // this
                        // more efficiently?
                        dx[0] = oldpos[0] - temppos[0];
                        dx[1] = oldpos[1] - temppos[1];
                        dx[2] = oldpos[2] - temppos[2];
                        dx[3] = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                        // count the flips and store the last flips of every
                        // type
                        lastflips[0] = i;
                        lastflips[1] = i;
                        lastflips[2] = i;
                        flips[0] = 0;
                        flips[1] = 0;
                        flips[2] = 0;
                        flips[3] = 0;
                        flips[4] = 0;
#ifdef FULLCHECK
                        numflips = _cells[i]->get_flips(i, _cells, _cuboid,
                                                        _periodic);
                        flips[numflips]++;
#else
                        numflips = _cells[i]->get_flips_minimal(
                                flips, i, _cells, _ghosts, _cuboid, _periodic);
#endif
                        // ngbs has not changed since the last time it was set
                        for(unsigned int j = 0; j < ngbs.size(); j++) {
                            if(ngbs[j] < (int)_cells.size()) {
#ifdef FULLCHECK
                                numflips = _cells[ngbs[j]]->get_flips(
                                        ngbs[j], i, _cells, _cuboid, _periodic);
                                flips[numflips]++;
#else
                                numflips = _cells[ngbs[j]]->get_flips_minimal(
                                        flips, ngbs[j], i, _cells, _ghosts,
                                        _cuboid, _periodic);
#endif
                                if(numflips == 1) {
                                    lastflips[0] = ngbs[j];
                                }
                                if(numflips == 2) {
                                    lastflips[1] = ngbs[j];
                                }
                                if(numflips == 4) {
                                    lastflips[2] = ngbs[j];
                                }
                            }
                        }
#ifdef VERBOSE
                        cerr << "B" << endl;
                        cerr << flips[0] << "\t" << flips[1] << "\t" << flips[2]
                             << "\t" << flips[3] << "\t" << flips[4] << endl;
#endif
                        // divide the factor by two
                        factor *= 0.5;
                        // if factor goes below 1.e-10, we will never solve the
                        // deviation and we better crash
                        if(factor < 1.e-10) {
                            cerr << "Problem could not be resolved. Crashing!"
                                 << endl;
                            throw AdaptiveMeshException();
                            ofstream ofile("wrong_cell.dat");
                            _cells[i]->complete(_cells, _ghosts, _cuboid);
                            _cells[i]->print(ofile, 0);
                            ofile.close();
                            cerr << flips[0] << "\t" << flips[1] << "\t"
                                 << flips[2] << "\t" << flips[3] << "\t"
                                 << flips[4] << endl;
                            cerr << dx[3] << endl;
                            my_exit();
                        }
#ifdef FULLCHECK
                        test = valid_flips(flips);
#else
                        test = valid_flips_minimal(flips);
#endif
                    }
                    //                ofstream lfile("last_mesh.dat");
                    //                _cells[i]->complete(_cells, _ghosts,
                    //                _cuboid);
                    //                _cells[i]->print(lfile, i);
                    //                lfile.close();
                    // we are left with a single deviation. Solve it
                    // all cells have stored the face and facengb corresponding
                    // to the last deviation of a certain type internally
                    // every operation can cause subsequent flips in the newly
                    // created edges and faces, which are checked for internally
                    // when the program exits the functions below, we should
                    // have
                    // a valid mesh for the current value of the displacement
                    if(flips[4]) {
#ifdef VERBOSE
                        cerr << "4 to 4 flip, haha" << endl;
#endif
                        _cells[lastflips[2]]->flip_face(lastflips[2], _cells,
                                                        _ghosts, _cuboid,
                                                        _periodic);
                    } else {
                        if(flips[2]) {
// removal
#ifdef VERBOSE
                            cerr << "removing face" << endl;
#endif
                            _cells[lastflips[1]]->remove_face(
                                    flips, lastflips[1], _cells, _ghosts,
                                    _cuboid, _periodic);
                        } else {
                            if(flips[1]) {
// insertion
#ifdef VERBOSE
                                cerr << "inserting face" << endl;
#endif
                                _cells[lastflips[0]]->insert_face(
                                        flips, lastflips[0], _cells, _ghosts,
                                        _cuboid, _periodic);
                            }
                        }
                    }
//                if(i > 382){
//                    ofstream vfile("wrongcell.dat");
////                    _cells[i]->complete(_cells, _ghosts, _cuboid);
////                    _cells[i]->print(vfile, 0);
////                    ngbs = _cells[i]->get_ngbs();
////                    for(unsigned int j = 0; j < ngbs.size(); j++){
////                        _cells[ngbs[j]]->complete(_cells, _ghosts, _cuboid);
////                        _cells[ngbs[j]]->print(vfile, j+1);
////                    }
//                    vector<int> ngbs(4);
//                    ngbs[0] = 630;
//                    ngbs[1] = 927;
//                    ngbs[2] = 383;
//                    ngbs[3] = 122;
//                    _cells[67]->print_copies(vfile, ngbs, _cells, _cuboid);
//                }
#ifdef CHECK_EVERY_STEP
                    if(i >= CHECK_EVERY_STEP) {
                        ofstream ofile("lastcheck.dat");
                        check_mesh(ofile);
                    }
#endif
                    // temppos is now a valid position and is the new oldpos
                    oldpos[0] = temppos[0];
                    oldpos[1] = temppos[1];
                    oldpos[2] = temppos[2];
                    // repeat the movement to the new position
                    _cells[i]->move(newpos);
                    // recalculate the magnitude of the displacement
                    dx[0] = oldpos[0] - newpos[0];
                    dx[1] = oldpos[1] - newpos[1];
                    dx[2] = oldpos[2] - newpos[2];
                    dx[3] = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                    // ngbs might have changed!
                    ngbs = _cells[i]->get_ngbs();
                    // count flips etc.
                    lastflips[0] = i;
                    lastflips[1] = i;
                    lastflips[2] = i;
                    flips[0] = 0;
                    flips[1] = 0;
                    flips[2] = 0;
                    flips[3] = 0;
                    flips[4] = 0;
#ifdef FULLCHECK
                    numflips =
                            _cells[i]->get_flips(i, _cells, _cuboid, _periodic);
                    flips[numflips]++;
#else
                    numflips = _cells[i]->get_flips_minimal(
                            flips, i, _cells, _ghosts, _cuboid, _periodic);
#endif
                    for(unsigned int j = 0; j < ngbs.size(); j++) {
                        if(ngbs[j] < (int)_cells.size()) {
#ifdef FULLCHECK
                            numflips = _cells[ngbs[j]]->get_flips(
                                    ngbs[j], i, _cells, _cuboid, _periodic);
                            flips[numflips]++;
#else
                            numflips = _cells[ngbs[j]]->get_flips_minimal(
                                    flips, ngbs[j], i, _cells, _ghosts, _cuboid,
                                    _periodic);
#endif
                            if(numflips == 1) {
                                lastflips[0] = ngbs[j];
                            }
                            if(numflips == 2) {
                                lastflips[1] = ngbs[j];
                            }
                            if(numflips == 4) {
                                lastflips[2] = ngbs[j];
                            }
                        }
                    }
#ifdef VERBOSE
                    cerr << "C" << endl;
                    cerr << flips[0] << "\t" << flips[1] << "\t" << flips[2]
                         << "\t" << flips[3] << "\t" << flips[4] << endl;
#endif
#ifdef FULLCHECK
                    test = valid_flips(flips);
#else
                    test = valid_flips_minimal(flips);
#endif
                    numloop++;
                    if(numloop > 100) {
                        cerr << "Exceeded maximum number of loops. Stopping."
                             << endl;
                        throw AdaptiveMeshException();
                        ofstream ofile("errorcells.dat");
                        _cells[i]->complete(_cells, _ghosts, _cuboid);
                        _cells[i]->print(ofile, i);
                        //                    for(unsigned int j = 0; j <
                        //                    ngbs.size(); j++){
                        //                        _cells[ngbs[j]]->complete(_cells,
                        //                        _ghosts, _cuboid);
                        //                        _cells[ngbs[j]]->print(ofile,
                        //                        ngbs[j]);
                        //                    }
                        _cells[i]->print_copies(ofile, ngbs, _cells, _cuboid);
                        ofile.close();
                        my_exit();
                    }
                }
                // we are left with a single deviation. Solve it
                if(flips[4]) {
#ifdef VERBOSE
                    cerr << "4 to 4 flip, hihi" << endl;
#endif
                    _cells[lastflips[2]]->flip_face(
                            lastflips[2], _cells, _ghosts, _cuboid, _periodic);
                } else {
                    if(flips[2]) {
// removal
#ifdef VERBOSE
                        cerr << "removal" << endl;
#endif
                        _cells[lastflips[1]]->remove_face(flips, lastflips[1],
                                                          _cells, _ghosts,
                                                          _cuboid, _periodic);
                    } else {
                        if(flips[1]) {
// insertion
#ifdef VERBOSE
                            cerr << "insertion" << endl;
#endif
                            _cells[lastflips[0]]->insert_face(
                                    flips, lastflips[0], _cells, _ghosts,
                                    _cuboid, _periodic);
                        }
                    }
                }
            }
            if(!mode) {
                // ngbs might have changed
                ngbs = _cells[i]->get_ngbs();
                // flag all neighbours to be moved
                for(unsigned int j = 0; j < ngbs.size(); j++) {
                    if(ngbs[j] < (int)_cells.size()) {
                        _cells[ngbs[j]]->semiactivate();
                    }
                    //                vector<int> ngbs2 =
                    //                _cells[ngbs[j]]->get_ngbs();
                    //                for(unsigned int k = 0; k < ngbs2.size();
                    //                k++){
                    //                    _cells[ngbs2[k]]->queue();
                    //                }
                }
            } else {
                if(mode == 1) {
                    ngbs = _cells[i]->get_ngbs();
                    // flag all neighbours to be moved
                    for(unsigned int j = 0; j < ngbs.size(); j++) {
                        if(ngbs[j] < (int)_cells.size()) {
                            _cells[ngbs[j]]->queue();
                        }
                    }
                }
            }
            // finalize
            _cells[i]->set_valid_pos(newpos);
            // keep the generator inside the periodic box if periodic boundaries
            // are used
            if(_periodic) {
                _cells[i]->keep_inside(i, _cells, _cuboid);
            }
#ifdef CHECK_EVERY_STEP
            if(i >= CHECK_EVERY_STEP) {
                ofstream ofile("lastcheck.dat");
                check_mesh(ofile);
            }
#endif
        }
    }
    //    ofstream ofile("lastcheck.dat");
    //    check_mesh(ofile);

    // calculate new vertex positions for all active cells
    for(unsigned int i = 0; i < _cells.size(); i++) {
        if(_cells[i]->queued() || _cells[i]->semiactive() ||
           _cells[i]->active()) {
            _cells[i]->complete(_cells, _ghosts, _cuboid);
        }
    }
    //    loop++;
    return count;
}

/**
 * @brief Get the characteristic radius of the cell with the given index
 *
 * The characteristic radius if the radius of a sphere with the same volume as
 * the cell.
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @return The characteristic radius of the cell
 */
double AdaptiveVorTess3d::get_h(unsigned int i) {
    return _cells[i]->get_h();
}

/**
 * @brief Get a list of AdaptiveFace3d faces to use during the hydrodynamical
 * flux exchange calculation
 *
 * Only faces that are active, i.e. have at least one active neighbouring cell,
 * are returned.
 *
 * @param currentTime unsigned long integer current simulation time
 * @return std::vector of AdaptiveFace3d instances that should be used during
 * the hydrodynamical flux exchange
 */
vector<AdaptiveFace3d*> AdaptiveVorTess3d::get_faces(
        unsigned long currentTime) {
    vector<AdaptiveFace3d*> faces;
    for(unsigned int i = 0; i < _cells.size(); i++) {
        vector<int> ngbs = _cells[i]->get_ngbs();
        bool active =
                (_cells[i]->get_particle()->get_starttime() == currentTime);
        for(unsigned int j = 0; j < ngbs.size(); j++) {
            // only construct the faces once
            if((int)i < ngbs[j]) {
                if(ngbs[j] < (int)_cells.size()) {
                    if(!active &&
                       _cells[ngbs[j]]->get_particle()->get_starttime() !=
                               currentTime) {
                        // both cells are inactive; we don't add the face
                        continue;
                    }
                    double midface[3];
                    double area = _cells[i]->get_face(j)->calculate_quantities(
                            midface);
                    // if the area is 0, we neglect this face
                    if(area) {
                        Vec pos;
                        if(_periodic) {
                            double posvec[3];
                            _cells[i]->get_periodic_position(j, _cells, posvec,
                                                             _cuboid);
                            pos.set(posvec[0], posvec[1], posvec[2]);
                        } else {
                            pos = _cells[ngbs[j]]
                                          ->get_particle()
                                          ->get_position();
                        }
                        faces.push_back(new AdaptiveFace3d(
                                _cells[i]->get_particle(),
                                _cells[ngbs[j]]->get_particle(), pos, midface,
                                area));
                    }
                } else {
                    if(!active) {
                        continue;
                    }
                    if(_periodic) {
                        if(ngbs[j] - _cells.size() < i) {
                            continue;
                        }
                    }
                    double midface[3];
                    _cells[i]->get_face(j)->get_midpoint(midface);
                    double area = _cells[i]->get_face(j)->get_area();
                    AdaptiveVorCell3d* ngb = _ghosts[ngbs[j] - _cells.size()];
                    GasParticle* right = NULL;
                    if(_periodic) {
                        right = ngb->get_particle();
                    }
                    double ngbposvec[3];
                    ngb->get_ghost_position(ngbposvec, _cuboid);
                    Vec ngbpos(ngbposvec[0], ngbposvec[1], ngbposvec[2]);
                    faces.push_back(
                            new AdaptiveFace3d(_cells[i]->get_particle(), right,
                                               ngbpos, midface, area));
                }
            }
        }
    }
    return faces;
}

/**
 * @brief Calculate the velocity of the cell with the given index, based on the
 * hydrodynamical variables of the given GasParticle
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @param particle GasParticle associated with the given cell
 * @return Velocity for the generator of the cell
 */
Vec AdaptiveVorTess3d::get_velocity(unsigned int i, GasParticle* particle) {
    return _cells[i]->get_velocity(particle);
}

/**
 * @brief Estimate gradients in the cell with the given index, based on the
 * hydrodynamical variables of the given GasParticle
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @param delta Array to store the results in
 * @param particle GasParticle associated with the cell
 */
void AdaptiveVorTess3d::estimate_gradients(unsigned int i, StateVector* delta,
                                           GasParticle* particle) {
    _cells[i]->estimate_gradients(delta, particle, _cells, _ghosts, _cuboid,
                                  _periodic);
}

/**
 * @brief Calculate the volume of the cell with the given index
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @return The volume of the cell
 */
double AdaptiveVorTess3d::get_volume(unsigned int i) {
    return _cells[i]->get_volume();
}

/**
 * @brief Calculate the centroid of the cell with the given index
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @return Coordinates of the centroid of the cell
 */
Vec AdaptiveVorTess3d::get_centroid(unsigned int i) {
    return _cells[i]->get_centroid();
}

/**
 * @brief Update the positions of the mesh generators using the given new
 * positions
 *
 * @param particles std::vector of GasParticle instances which contain the new
 * positions
 */
void AdaptiveVorTess3d::update_positions(ParticleVector& particles) {
    _new_positions.resize(particles.gassize());
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        _new_positions[particles.gas(i)->get_local_id()] =
                particles.gas(i)->get_position();
        _cells[particles.gas(i)->get_local_id()]->reset_flags();
    }
}

/**
 * @brief Check the validity of the Voronoi tesselation
 *
 * The validity is checked using the geometrical in sphere criterion.
 *
 * @param stream std::ostream to write error information to
 */
void AdaptiveVorTess3d::check_mesh(ostream& stream) {
#ifdef VERBOSE
    cerr << "checking mesh" << endl;
#endif
    for(unsigned int i = 0; i < _cells.size(); i++) {
        AdaptiveVorCell3d* cell = _cells[i];
        cell->complete(_cells, _ghosts, _cuboid);
    }
    for(unsigned int i = 0; i < _cells.size(); i++) {
        AdaptiveVorCell3d* cell = _cells[i];
        if(!(i % 1000)) {
            cout << "checking cell " << i << endl;
        }
        cell->check(i, _cells, stream, _cuboid);
    }
#ifdef VERBOSE
    cerr << "done" << endl;
#endif
}

/**
 * @brief Get the indices of neighbouring cell of the cell with the given index
 *
 * @note This function is a dummy and is not functional (yet).
 *
 * @param index unsigned integer index of a cell in the tesselation
 * @return std::vector of unsigned long indices of neighbouring cells
 */
vector<unsigned long> AdaptiveVorTess3d::get_ngb_ids(unsigned int index) {
    vector<unsigned long> ngb_ids;
    return ngb_ids;
}

/**
 * @brief Get the new positions
 *
 * We need this function when the evolution algorithm crashes, since we then
 * have to update all positions before calling the old algorithm.
 *
 * @return The new positions for all mesh generators
 */
vector<Vec> AdaptiveVorTess3d::get_new_positions() {
    return _new_positions;
}

/**
 * @brief Dump the tesselation to the given RestartFile
 *
 * @param rfile RestartFile to write to
 * @param particles std::vector of GasParticle instances used for sorting the
 * output
 */
void AdaptiveVorTess3d::dump(RestartFile& rfile,
                             std::vector<GasParticle*>& particles) {
    _cuboid.dump(rfile);
    unsigned int vsize = _cells.size();
    rfile.write(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        _cells[particles[i]->get_local_id()]->dump(rfile);
    }
    vsize = _ghosts.size();
    rfile.write(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        _ghosts[i]->dump(rfile);
    }
    vsize = _new_positions.size();
    rfile.write(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        Vec new_pos = _new_positions[i];
        rfile.write(new_pos);
    }
    rfile.write(_periodic);
    rfile.write(_tolerance);
}

/**
 * @brief Restart constructor. Initialize the tesselation based on the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 * @param particles std::vector of GasParticle instances used for sorting the
 * cells
 */
AdaptiveVorTess3d::AdaptiveVorTess3d(RestartFile& rfile,
                                     std::vector<GasParticle*>& particles)
        : _cuboid(rfile) {
    unsigned int vsize;
    rfile.read(vsize);
    _cells.resize(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        _cells[particles[i]->get_local_id()] =
                new AdaptiveVorCell3d(rfile, particles[i]);
    }
    rfile.read(vsize);
    _ghosts.resize(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        _ghosts[i] = new AdaptiveVorCell3d(rfile, NULL);
    }
    rfile.read(vsize);
    _new_positions.resize(vsize);
    for(unsigned int i = 0; i < vsize; i++) {
        Vec new_pos;
        rfile.read(new_pos);
        _new_positions[i] = new_pos;
    }
    rfile.read(_periodic);
    rfile.read(_tolerance);
}
