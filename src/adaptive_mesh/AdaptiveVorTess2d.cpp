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
 * @file AdaptiveVorTess2d.cpp
 *
 * @brief 2D evolving mesh: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdaptiveVorTess2d.hpp"
#include "AdaptiveCellParallelSorter.hpp"  // for AdaptiveCellParallelSorter
#include "AdaptiveFace2d.hpp"              // for AdaptiveFace2d
#include "AdaptiveMeshUtils.hpp"           // for get_wallpos, get_wall, etc
#include "AdaptiveVorCell2d.hpp"           // for AdaptiveVorCell2d
#include "AdaptiveVorFace2d.hpp"           // for AdaptiveVorFace2d
#include "Error.hpp"                       // for my_exit
#include "MPIGlobal.hpp"                   // for size, rank, sendbuffer, etc
#include "MPIMethods.hpp"                  // for MyMPI_Pack, MyMPI_Isend, etc
#include "RestartFile.hpp"                 // for RestartFile
#include "VorCell.hpp"                     // for VorCell
#include "VorGen.hpp"                      // for VorGen
#include "VorTess.hpp"                     // for VorTess
#include "mpi.h"                           // for MPI_Status, MPI_INT, etc
#include "predicates.hpp"                  // for incircle_old
#include "utilities/Cuboid.hpp"            // for Cuboid
#include "utilities/GasParticle.hpp"       // for GasParticle
#include "utilities/ParticleVector.hpp"    // for ParticleVector
#include <algorithm>                       // for min, max, sort
#include <cstddef>                         // for NULL
#include <ext/alloc_traits.h>
#include <iostream>  // for cerr, cout
using namespace std;

/**
 * @brief Constructor
 *
 * Initialize a new evolving mesh bashed on the given already constructed
 * VorTess.
 *
 * @param tesselation Valid Voronoi tesselation
 * @param cuboid Cuboid specifying the dimensions of the simulation box
 * @param maxid Maximum ID of any particle that will ever be added to the
 * tesselation
 * @param particles ParticleVector holding the local particles of the simulation
 * @param periodic Bool specifying if we deal with a periodic (true) or
 * reflective (false) simulation box
 */
AdaptiveVorTess2d::AdaptiveVorTess2d(VorTess* tesselation, Cuboid cuboid,
                                     unsigned long maxid,
                                     ParticleVector& particles, bool periodic)
        : _newlist(tesselation->get_cells().size()), _newerlist(maxid, cuboid),
          _periodic(periodic), _cuboid(cuboid) {

    double minanchor = std::min(cuboid.get_anchor()[0], cuboid.get_anchor()[1]);
    double maxsides = std::max(cuboid.get_sides()[0], cuboid.get_sides()[1]);
#if ndim_ == 3
    minanchor = std::min(minanchor, cuboid.get_anchor()[2]);
    maxsides = std::max(maxsides, cuboid.get_sides()[2]);
    Vec anchor(2. * minanchor, 2. * minanchor, 2. * minanchor);
    Vec sides(3. * maxsides, 3. * maxsides, 3. * maxsides);
#else
    Vec anchor(2. * minanchor, 2. * minanchor);
    Vec sides(3. * maxsides, 3. * maxsides);
#endif
    _box12 = Cuboid(anchor, sides);

    vector<VorCell*> cells = tesselation->get_cells();
    unsigned int cellsize = cells.size();
    _cells.resize(cellsize, NULL);
    // orphans are ghost cells for cells that are located on other processes
    // since we already have a list of ghosts for local ghosts, we have to
    // somehow distinguish between:
    //  (a) normal cells
    //  (b) local ghosts
    //  (c) domain ghosts
    // and guarantee that all of these have distinct indices in ranges that
    // can be tested.
    // We opt for the following conventions:
    //  (a) have indices in the range [0,cellsize[
    //  (b) have indices >= cellsize
    //  (c) have indices strictly smaller than 0
    // To easily handle (c), we fill the 0th element in the list. This means
    // we only have to make sure never to use the 0th element. If we need
    // an orphan, we can simply use -index to access it...
    //    _orphans.push_back(NULL);
    //    cerr << "Creating cells" << endl;
    for(unsigned int i = 0; i < cellsize; i++) {
        VorCell* cell = cells[i];
        Vec position = cell->get_particle()->get_position();
        vector<VorGen*> ngbs = cell->get_ngbs();
        AdaptiveVorCell2d* newcell = new AdaptiveVorCell2d(
                position, ngbs.size(), 0, cell->get_particle());
        _newerlist.obtain_position(cell->get_particle()->id(), MPIGlobal::rank);
        newcell->set_id(cell->get_particle()->id());
        _cells[cell->get_particle()->get_local_id()] = newcell;
        _newlist.add_normal_cell(cell->get_particle()->get_local_id(), newcell);
        for(unsigned int j = 0; j < ngbs.size(); j++) {
            // there are two types of ghosts: ghosts with a particle and ghosts
            // without
            // ghosts without particle are reflective ghosts
            // ghosts with particle can be further subdivided in periodic ghosts
            // and domain ghosts
            // periodic ghosts are ghosts with a particle that resides on the
            // local process and where the ghost position differs from the
            // particle position
            // domain ghosts are ghosts that reside on another process and that
            // have the same position as their associated (ghost) particle
            if(ngbs[j]->get_particle()) {
                // ngbwall is only used with periodic boundaries
                //
                // the conditions are quite vaguely structured. Here is the
                // idea:
                //  - we always check if we have a ghost
                //  - if we have a ghost, we discriminate between periodic and
                //    domain ghosts
                if((ngbs[j]->get_position() -
                    ngbs[j]->get_particle()->get_position())
                           .norm2()) {

                    newcell->add_ngb(ngbs[j]->get_particle()->get_local_id());
                    double pos[2];
                    pos[0] = ngbs[j]->x();
                    pos[1] = ngbs[j]->y();
                    double opos[2];
                    opos[0] = ngbs[j]->get_particle()->x();
                    opos[1] = ngbs[j]->get_particle()->y();
                    int wall = AdaptiveMeshUtils::get_wall(pos, opos);
                    unsigned long id =
                            _newerlist.add_ghost(ngbs[j]->get_particle()->id(),
                                                 wall, MPIGlobal::rank);
                    newcell->add_ngb_particleID(id);
                    if(_periodic) {
                        double pos[2];
                        pos[0] = ngbs[j]->x();
                        pos[1] = ngbs[j]->y();
                        double opos[2];
                        opos[0] = ngbs[j]->get_particle()->x();
                        opos[1] = ngbs[j]->get_particle()->y();
                        int id = AdaptiveMeshUtils::get_wall(pos, opos);
                        newcell->add_ngbwall(id);
                        // the line below queues this cell at the boundary for
                        // a consistency check after all cells have been added
                        newcell->queue();
                    } else {
                        cerr << "You should never end up here. Massive error!"
                             << endl;
                        my_exit();
                    }
                } else {
                    if(ngbs[j]->get_id() == 1) {
                        double pos[2];
                        pos[0] = ngbs[j]->get_position().x();
                        pos[1] = ngbs[j]->get_position().y();
                        int wall =
                                AdaptiveMeshUtils::trim_wallpos(pos, _cuboid);
                        Vec orphanpos(pos[0], pos[1]);
                        AdaptiveVorCell2d* orphan = new AdaptiveVorCell2d(
                                orphanpos, 0,
                                ngbs[j]->get_particle()->get_local_id());
                        orphan->set_id(ngbs[j]->get_particle()->id());
                        //                        _orphans.push_back(orphan);
                        orphan->set_rank(ngbs[j]->get_process());
                        orphan->add_ngb(cell->get_particle()->get_local_id());
                        int index = _newlist.add_orphan_cell(orphan);
                        newcell->add_ngb(index);
                        newcell->add_ngbwall(wall);
                        unsigned long id;
                        if(wall) {
                            id = _newerlist.add_ghost(
                                    ngbs[j]->get_particle()->id(), wall,
                                    ngbs[j]->get_process());
                        } else {
                            id = ngbs[j]->get_particle()->id();
                            _newerlist.obtain_position(id,
                                                       ngbs[j]->get_process());
                        }
                        newcell->add_ngb_particleID(id);
                    } else {
                        newcell->add_ngb(
                                ngbs[j]->get_particle()->get_local_id());
                        newcell->add_ngbwall(0);
                        newcell->add_ngb_particleID(
                                ngbs[j]->get_particle()->id());
                        _newerlist.obtain_position(
                                ngbs[j]->get_particle()->id(), MPIGlobal::rank);
                    }
                }
            } else {
                AdaptiveVorCell2d* ghost = new AdaptiveVorCell2d(
                        ngbs[j]->get_position(), 0, ngbs[j]->get_original());
                //                _ghosts.push_back(ghost);
                int index = _newlist.add_ghost_cell(ghost);
                newcell->add_ngb(index);
                newcell->add_ngbwall(0);
                // wall 0 seems incorrect...
                unsigned long id =
                        _newerlist.add_ghost(ngbs[j]->get_original(), 0);
                newcell->add_ngb_particleID(id);
            }
        }
    }
    //    cerr << "Done." << endl;

    _newerlist.calculate_positions(particles, _box12);
    for(unsigned int i = 0; i < _cells.size(); i++) {
        _cells[i]->complete(_newerlist);
    }

    //    stringstream oname;
    //    oname << "initial_grid_" << MPIGlobal::rank << ".txt";
    //    ofstream ofile(oname.str());
    //    for(auto it = _newlist.normalbegin(); it != _newlist.normalend();
    //    ++it){
    //        it->print(ofile, it.index());
    //    }
    //    ofile.close();
    //    MyMPI_Barrier();
    //    my_exit();
    return;

    //    for(auto it = _newlist.normalbegin(); it != _newlist.normalend();
    //    ++it){
    //        cerr << "Cell: " << it->get_particleID() << endl;
    //        vector<unsigned long> &ngbs = it->get_newngbs();
    //        for(unsigned int j = 0; j < ngbs.size(); j++){
    //            cerr << ngbs[j];
    //            if(ngbs[j] > 10){
    //                unsigned long oid = (ngbs[j]-10)>>3;
    //                int wall = 10+8*oid-1-ngbs[j];
    //                cerr << " (" << oid << "," << wall << ")";
    //            }
    //            cerr << " (" << _newerlist.is_local(ngbs[j]) << ")";
    //            cerr << endl;
    //        }
    //    }

    //    stringstream oname;
    //    oname << "initial_grid_" << MPIGlobal::rank << ".txt";
    //    ofstream ofile(oname.str());
    //    for(auto it = _newlist.normalbegin(); it != _newlist.normalend();
    //    ++it){
    //        it->complete(_newerlist);
    //        it->print(ofile, it.index());
    //    }
    //    ofile.close();
    //    MyMPI_Barrier();
    //    my_exit();

    if(!_periodic) {
        vector<vector<unsigned int> > ghostids(cellsize);
        for(unsigned int i = 0; i < cellsize; i++) {
            ghostids[i].resize(8, 0);
        }

        vector<unsigned int> replist;
        vector<unsigned int> repids;
        for(AdaptiveCellList<AdaptiveVorCell2d>::ghostiterator it =
                    _newlist.ghostbegin();
            it != _newlist.ghostend(); ++it) {
            double pos[2];
            it->get_position(pos);
            double opos[2];
            _newlist.get_normal(it->get_original())->get_position(opos);
            int id = AdaptiveMeshUtils::get_wall(pos, opos);
            it->add_ngb(id);
            it->move(opos);
            if(ghostids[it->get_original()][id + 8]) {
                replist.push_back(it.index());
                repids.push_back(ghostids[it->get_original()][id + 8]);
            } else {
                ghostids[it->get_original()][id + 8] = it.index();
            }
        }

        for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                    _newlist.normalbegin();
            it != _newlist.normalend(); ++it) {
            vector<int> ngbs = it->get_ngbs();
            for(unsigned int j = 0; j < ngbs.size(); j++) {
                unsigned int k = 0;
                while(k < replist.size() && ngbs[j] != (int)replist[k]) {
                    k++;
                }
                if(k < replist.size()) {
                    it->change_ngb(j, repids[k]);
                }
            }
        }
    } else {
        cerr << "Checking mesh consistency..." << endl;
        // since the periodic ghosts are added on a cell by cell basis, it
        // can happen that neighbour relations between cells at the boundary
        // are not mutual (especially if the vertices are actually fourth
        // order), since the periodic positions are different depending on
        // the cell you use as a reference point
        // long story short: if a ngb of cell i does not have i as a ngb, we
        // add it manually
        for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                    _newlist.normalbegin();
            it != _newlist.normalend(); ++it) {
            if(it->queued()) {
                vector<int> ngbs = it->get_ngbs();
                vector<int> ngbwalls = it->get_ngbwalls();
                for(unsigned int j = 0; j < ngbs.size(); j++) {
                    if(_newlist.is_normal(ngbs[j])) {
                        if(!_newlist.get_normal(ngbs[j])->is_ngb(it.index())) {
                            _newlist.get_normal(ngbs[j])->add_ngb(
                                    ngbs[(j + ngbs.size() - 1) % ngbs.size()],
                                    ngbs[(j + 1) % ngbs.size()], it.index(),
                                    AdaptiveMeshUtils::get_wallpos(ngbwalls[j],
                                                                   0));
                        }
                    }
                }
            }
        }
        cerr << "Done." << endl;
    }

    cerr << "Cleaning up duplicate orphans..." << endl;
    // we need to make sure every orphan is only represented once:
    // we want only one orphan copy of cell x on process y on this process
    // (if unclear: cell x is on process y, now re-read the line above)
    // the situation is similar to the one with the ghosts, where we want
    // to make sure we only have one copy for every wall
    // problem here is that we do not know the number of cells on other
    // processes
    // this makes it difficult to initialize the vector...
    vector<vector<unsigned int> > orphanids(cellsize);
    for(unsigned int i = 0; i < cellsize; i++) {
        orphanids[i].resize(MPIGlobal::size, 0);
    }

    vector<unsigned int> repids(_newlist.orphansize());
    vector<int> to_keep;
    for(AdaptiveCellList<AdaptiveVorCell2d>::orphaniterator it =
                _newlist.orphanbegin();
        it != _newlist.orphanend(); ++it) {
        while(it->get_original() >= orphanids.size()) {
            orphanids.push_back(vector<unsigned int>(MPIGlobal::size, 0));
        }
        if(!orphanids[it->get_original()][it->get_rank()]) {
            to_keep.push_back(it.index());
            orphanids[it->get_original()][it->get_rank()] = -to_keep.size();
        }
        repids[-it.index()] = orphanids[it->get_original()][it->get_rank()];
    }

    cerr << "To keep:" << endl;
    for(unsigned int i = 0; i < to_keep.size(); i++) {
        cerr << to_keep[i] << endl;
    }

    // now we have a list of orphans and for every orphan also the index
    // of the orphan to use for that orphan. Most of the orphans will be
    // duplicates and will never be used thanks to this list with indices.
    // we want to remove them.
    _newlist.update_orphans(to_keep);

    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                _newlist.normalbegin();
        it != _newlist.normalend(); ++it) {
        it->update_orphans(repids, _newlist);
    }
    cerr << "Done." << endl;

    // last part: we need to make sure the orphans have correct ngbs
    // to do this, we can
    // communicate the ngbs from the process that holds the original cell
    //  or
    // use the ngbs from the local cells that have this orphan as a ngb
    //
    // the first version has the advantage of not requiring an awful lot of
    // sorting. However, we
    // (1) communicate too much, because we also include ngbs that have
    //     no orphans on this process
    // (2) will experience trouble trying to link local orphan ngbs to imported
    //     cells. The orphans hold info about their original cell, but not the
    //     other way around...
    //
    // the second version means we have to loop over all orphans, pick a ngb
    // (we can store one ngb when making the orphan. Since we only retain one
    // copy, which one will be quite random) and then walk around the orphan
    // using ngbs and the ordering in these ngbs. This will pick up other
    // orphan ngbs (but we do have to walk in two directions then...)
    cerr << "Setting orphan ngbs..." << endl;
    for(AdaptiveCellList<AdaptiveVorCell2d>::orphaniterator it =
                _newlist.orphanbegin();
        it != _newlist.orphanend(); ++it) {
        int firstngb = it->get_first_ngb();
        AdaptiveVorCell2d* cell = _newlist.get_normal(firstngb);
        vector<int> ngbs = cell->get_ngbs();
        unsigned int index = 0;
        while(index < ngbs.size() && ngbs[index] != it.index()) {
            index++;
        }
        if(index == ngbs.size()) {
            cerr << "Orphan is not ngb of its ngb!" << endl;
            my_exit();
        }
        // index now holds the index of the current orphan in the ngb list of
        // the current cell
        // we want to find the first ngb of orphan (which will be an orphan), so
        // we need to go back in the orphan ngb list. This means going forward
        // in the current ngb list
        int nextcell = ngbs[(index + 1) % ngbs.size()];
        while(!_newlist.is_orphan(nextcell)) {
            cell = _newlist.get_normal(nextcell);
            ngbs = cell->get_ngbs();
            index = 0;
            while(index < ngbs.size() && ngbs[index] != it.index()) {
                index++;
            }
            if(index == ngbs.size()) {
                cerr << "Orphan is not a ngb of this cell!" << endl;
                my_exit();
            }
            nextcell = ngbs[(index + 1) % ngbs.size()];
        }
        // nextcell now contains the first ngb of orphan, which is an orphan
        // itself
        // change the first ngb to orphan
        it->change_ngb(0, nextcell);
        // we need to add the ngbwall for this orphan
        // this is the ngbwall of the orphan in the current cell, but also
        // taking into account the possible ngbwall of it in cell
        vector<int> ngbwalls = cell->get_ngbwalls();
        int itwall = ngbwalls[index];
        int orphanwall = ngbwalls[(index + 1) % ngbs.size()];
        it->add_ngbwall(AdaptiveMeshUtils::get_wallpos(itwall, orphanwall));
        // cell still contains the last valid normal cell. We again loop until
        // we find an orphan, but now in the opposite direction. And we add
        // new ngbs to orphan on the fly
        it->add_ngb(cell->get_particle()->get_local_id());
        it->add_ngbwall(AdaptiveMeshUtils::get_wallpos(itwall, 0));
        nextcell = ngbs[(index + ngbs.size() - 1) % ngbs.size()];
        while(!_newlist.is_orphan(nextcell)) {
            cell = _newlist.get_normal(nextcell);
            ngbs = cell->get_ngbs();
            index = 0;
            while(index < ngbs.size() && ngbs[index] != it.index()) {
                index++;
            }
            if(index == ngbs.size()) {
                cerr << "Orphan is not a ngb of this cell!" << endl;
                my_exit();
            }
            ngbwalls = cell->get_ngbwalls();
            itwall = ngbwalls[index];
            it->add_ngb(nextcell);
            it->add_ngbwall(AdaptiveMeshUtils::get_wallpos(itwall, 0));
            nextcell = ngbs[(index + ngbs.size() - 1) % ngbs.size()];
        }
        // add the orphan as a last ngb
        ngbwalls = cell->get_ngbwalls();
        itwall = ngbwalls[index];
        orphanwall = ngbwalls[(index + ngbs.size() - 1) % ngbs.size()];
        it->add_ngb(nextcell);
        it->add_ngbwall(AdaptiveMeshUtils::get_wallpos(itwall, orphanwall));
    }
    // for very small grids, it can happen that an orphan has more neighbours
    // than were added above, because the algorithm above only adds neighbours
    // in between other neighbouring orphans. If a neihgbouring orphan is
    // enclosed by two local cells (which can happen in some cases), one of
    // the normal neighbours might not be included. We explicitly check for this
    // eventuality below
    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                _newlist.normalbegin();
        it != _newlist.normalend(); ++it) {
        vector<int> ngbs = it->get_ngbs();
        vector<int> thiswalls = it->get_ngbwalls();
        for(unsigned int i = 0; i < ngbs.size(); i++) {
            if(_newlist.is_orphan(ngbs[i])) {
                AdaptiveVorCell2d* orphan = _newlist.get_orphan(ngbs[i]);
                if(!orphan->is_ngb(it.index())) {
                    cerr << "Not a neighbour!" << endl;
                    cerr << ngbs[i] << endl;
                    cerr << it.index() << endl;
                    cerr << "ngbs:" << endl;
                    vector<int> ongbs = orphan->get_ngbs();
                    for(unsigned int j = 0; j < ongbs.size(); j++) {
                        cerr << j << ": " << ongbs[j] << endl;
                    }
                    cerr << "ngbs[i-1]: "
                         << ngbs[(i + ngbs.size() - 1) % ngbs.size()] << endl;
                    cerr << "ngbs[i+1]: " << ngbs[(i + 1) % ngbs.size()]
                         << endl;
                    // ngbs and walls to add
                    vector<int> to_add;
                    vector<int> addwalls;

                    // loop forward until we find an orphan
                    int nextngb = ngbs[(i + 1) % ngbs.size()];
                    int nextwall = thiswalls[(i + 1) % ngbs.size()];
                    while(!_newlist.is_orphan(nextngb)) {
                        to_add.insert(to_add.begin(), nextngb);
                        AdaptiveVorCell2d* cell = _newlist.get_normal(nextngb);
                        vector<int> nngbs = cell->get_ngbs();
                        unsigned int index = 0;
                        while(index < nngbs.size() && nngbs[index] != ngbs[i]) {
                            index++;
                        }
                        if(index == nngbs.size()) {
                            cerr << "Orphan is not a ngb of this cell!" << endl;
                            my_exit();
                        }
                        vector<int> ngbwalls = cell->get_ngbwalls();
                        int itwall = ngbwalls[index];
                        addwalls.insert(addwalls.begin(), itwall);
                        nextngb = nngbs[(index + 1) % nngbs.size()];
                        nextwall = ngbwalls[(index + 1) % nngbs.size()];
                    }
                    // add the orphan if it was not already added
                    if(!orphan->is_ngb(nextngb)) {
                        to_add.insert(to_add.begin(), nextngb);
                        addwalls.insert(addwalls.begin(), nextwall);
                    }

                    // add this ngb
                    to_add.push_back(it.index());
                    addwalls.push_back(0);

                    // loop backward until we find an orphan
                    nextngb = ngbs[(i + ngbs.size() - 1) % ngbs.size()];
                    nextwall = thiswalls[(i + ngbs.size() - 1) % ngbs.size()];
                    while(!_newlist.is_orphan(nextngb)) {
                        to_add.push_back(nextngb);
                        AdaptiveVorCell2d* cell = _newlist.get_normal(nextngb);
                        vector<int> nngbs = cell->get_ngbs();
                        unsigned int index = 0;
                        while(index < nngbs.size() && nngbs[index] != ngbs[i]) {
                            index++;
                        }
                        if(index == nngbs.size()) {
                            cerr << "Orphan is not a ngb of this cell!" << endl;
                            my_exit();
                        }
                        vector<int> ngbwalls = cell->get_ngbwalls();
                        int itwall = ngbwalls[index];
                        addwalls.push_back(itwall);
                        nextngb = nngbs[(index + nngbs.size() - 1) %
                                        nngbs.size()];
                        nextwall = ngbwalls[(index + nngbs.size() - 1) %
                                            nngbs.size()];
                    }
                    // add the orphan if it was not already added
                    if(!orphan->is_ngb(nextngb)) {
                        to_add.push_back(nextngb);
                        addwalls.push_back(nextwall);
                    }

                    // now actually add the new neighbours
                    for(unsigned int k = 0; k < to_add.size(); k++) {
                        orphan->add_ngb(to_add[k]);
                        if(_newlist.is_orphan(to_add[k])) {
                            orphan->add_ngbwall(AdaptiveMeshUtils::get_wallpos(
                                    thiswalls[i], addwalls[k]));
                        } else {
                            orphan->add_ngbwall(AdaptiveMeshUtils::get_wallpos(
                                    thiswalls[i], 0));
                        }
                    }
                }
            }
        }
    }
    cerr << "Done." << endl;

    cerr << "Outputting..." << endl;
    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                _newlist.normalbegin();
        it != _newlist.normalend(); ++it) {
        AdaptiveVorCell2d* cell = *it;
        cerr << "Cell " << cell->get_id() << endl;
        double pos[2];
        cell->get_position(pos);
        cerr << pos[0] << "\t" << pos[1] << endl;
        vector<int> ngbs = cell->get_ngbs();
        vector<int> ngbwalls = cell->get_ngbwalls();
        for(unsigned int i = 0; i < ngbs.size(); i++) {
            cerr << ngbs[i] << "\t" << ngbwalls[i] << endl;
        }
        cerr << endl;
    }
    cerr << "Done." << endl;

    cerr << "Outputting orphans..." << endl;
    for(AdaptiveCellList<AdaptiveVorCell2d>::orphaniterator it =
                _newlist.orphanbegin();
        it != _newlist.orphanend(); ++it) {
        AdaptiveVorCell2d* orphan = *it;
        cerr << "Orphan " << it.index() << " (" << orphan->get_rank() << ": "
             << orphan->get_original() << ")" << endl;
        double pos[2];
        orphan->get_position(pos);
        cerr << pos[0] << "\t" << pos[1] << endl;
        vector<int> ngbs = orphan->get_ngbs();
        vector<int> ngbwalls = orphan->get_ngbwalls();
        for(unsigned int i = 0; i < ngbs.size(); i++) {
            cerr << ngbs[i] << "\t" << ngbwalls[i] << endl;
        }
        cerr << endl;
    }

    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                _newlist.normalbegin();
        it != _newlist.normalend(); ++it) {
        //        it->complete(_cells, _ghosts, _orphans, _cuboid, _periodic);
        it->complete(_newlist, _cuboid, _periodic);
    }

    //    save_restart(0);
}

/**
 * @brief Destructor
 *
 * Since we now use a custom list to store the cells, nothing has to be done
 * here anymore.
 */
AdaptiveVorTess2d::~AdaptiveVorTess2d() {
    //    for(unsigned int i = 0; i < _cells.size(); i++){
    //        delete _cells[i];
    //    }
    //    for(unsigned int i = 0; i < _ghosts.size(); i++){
    //        delete _ghosts[i];
    //    }
}

/**
 * @brief Output the tesselation to the given stream in a format that can be
 * plotted using gnuplot
 *
 * This can be done using the command
\verbatim
plot "<FILENAME>" with linespoints
\endverbatim
 *
 * @param stream std::ostream to write to
 */
void AdaptiveVorTess2d::print_tesselation_gnuplot(ostream& stream) {
    //    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
    //        _newlist.normalbegin(); it != _newlist.normalend(); ++it){
    ////        it->complete(_newlist, _cuboid, _periodic);
    //        it->complete(_newerlist);
    //        it->print(stream, it.index());
    //    }
    for(unsigned int i = 0; i < _cells.size(); i++) {
        _cells[i]->complete(_newerlist);
        _cells[i]->print(stream, i);
    }
}

/**
 * @brief Kernel of the mesh evolution algorithm: make the Voronoi grid valid
 * again
 *
 * Here we should put some more info on the algorithm...
 *
 * @param particles GasParticle list corresponding to the cells of the grid
 * @param currentTime unsigned long integer simulation time
 * @return The number of active cells in the grid for which the GasParticle time
 * equals the current simulation time
 */
unsigned int AdaptiveVorTess2d::update(ParticleVector& particles,
                                       unsigned long currentTime) {
    //    for(unsigned int i = 0; i < _cells.size(); i++){
    //        _cells[i]->complete(_cells, _ghosts, _cuboid, _periodic);
    //    }
    //    ofstream ofile("wrongfile.dat");
    //    print_tesselation_gnuplot(ofile);
    //    ofile.close();

    // flag the cells that have to be updated: the active cells and their
    // neighbours + their neighbours. The first two are needed for the hydro
    // and have to be fully reconstructed. The last ones need to be checked
    // to ensure the correctness of the neighbours, since a wrong face can
    // only be detected in 2 of the 4 involved cells.
    //    unsigned int count = 0;
    //    for(unsigned int i = 0; i < _cells.size(); i++){
    //        if(particles[i]->get_endtime() == currentTime){
    //            unsigned int index = particles[i]->id();
    //            _cells[index]->queue();
    //            count++;
    //            vector<int> ngbs = _cells[index]->get_ngbs();
    //            for(unsigned int j = 0; j < ngbs.size(); j++){
    //                if(ngbs[j] < _cells.size()){
    //                    _cells[ngbs[j]]->queue();
    //                    vector<int> ngbs2 = _cells[ngbs[j]]->get_ngbs();
    //                    for(unsigned int k = 0; k < ngbs2.size(); k++){
    //                        if(ngbs2[k] < _cells.size()){
    //                            _cells[ngbs2[k]]->queue();
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }

    //    vector<AdaptiveVorCell2d*> backup(_newlist.normalsize());
    //    vector<unsigned int> new_ids(_newlist.normalsize());
    //    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
    //        _newlist.normalbegin(); it != _newlist.normalend(); ++it){
    //        backup[it.index()] = *it;
    //        new_ids[it.index()] = it->get_particle()->get_local_id();
    //    }

    ////    for(unsigned int i = 0; i < _cells.size(); i++){
    ////        _cells[new_ids[i]] = backup[i];
    ////    }

    //    _newlist.reorder(new_ids);
    //    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
    //        _newlist.normalbegin(); it != _newlist.normalend(); ++it){
    //        it->update_ngbs(new_ids);
    //    }

    //    // ghosts store the id of the particle they are a ghost of
    //    // the ghost id itself does not change, but the particle id does
    //    for(AdaptiveCellList<AdaptiveVorCell2d>::ghostiterator it =
    //        _newlist.ghostbegin(); it != _newlist.ghostend(); ++it){
    //        it->update_local_id(new_ids);
    //    }

    //    // update the local ids for the orphans.
    //    // we have to communicate for this
    //    if(MPIGlobal::size > 1){
    //        // (1) fill buffers
    //        vector< vector<int> > request_ids(MPIGlobal::size);
    //        vector< vector<unsigned int> > requests(MPIGlobal::size);
    //        for(AdaptiveCellList<AdaptiveVorCell2d>::orphaniterator it =
    //            _newlist.orphanbegin(); it != _newlist.orphanend(); ++it){
    //            request_ids[it->get_rank()].push_back(it->get_original());
    //            requests[it->get_rank()].push_back(it.index());
    //        }

    //        std::vector<MPI_Request> reqs(MPIGlobal::size*2,
    //        MPI_REQUEST_NULL);
    //        // we have to store the positions we send in a buffer that is not
    //        // deleted when an iteration ends, since we use a non-blocking
    //        send
    //        // when receiving, we can use a buffer with iteration scoop, since
    //        // the data is only used inside the iteration
    //        vector< vector<unsigned int> > sendids(MPIGlobal::size);
    //        for(int i = 0; i < MPIGlobal::size; i++){
    //            if(i != MPIGlobal::rank){
    //                unsigned int vsize = request_ids[i].size();
    //                MyMPI_Isend(&request_ids[i][0], vsize, MPI_INT, i, 0,
    //                            &reqs[2*i]);
    //            }
    //        }

    //        for(int i = 0; i < 2*(MPIGlobal::size-1); i++){
    //            MPI_Status status;
    //            MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
    //            int tag = status.MPI_TAG;
    //            int source = status.MPI_SOURCE;
    //            int nelements;
    //            if(tag == 0){
    //                MyMPI_Get_count(&status, MPI_INT, &nelements);
    //                vector<int> ids(nelements);
    //                MyMPI_Recv(&ids[0], nelements, MPI_INT, source, 0,
    //                &status);
    ////                cerr << MPIGlobal::rank << ": received " << nelements
    ////                << " ids from " << source << endl;
    //                sendids[source].resize(nelements);
    //                for(unsigned int j = 0; j < ids.size(); j++){
    //                    sendids[source][j] = new_ids[ids[j]];
    //                }
    //                MyMPI_Isend(&sendids[source][0], sendids[source].size(),
    //                            MPI_INT, source, 1, &reqs[2*source+1]);
    ////                cerr << MPIGlobal::rank << ": sent "
    ////                << sendpositions[source].size() << " positions to " <<
    /// source
    ////                << endl;
    //            }
    //            if(tag == 1){
    //                MyMPI_Get_count(&status, MPI_INT, &nelements);
    //                vector<unsigned int> recvids(nelements);
    //                MyMPI_Recv(&recvids[0], nelements, MPI_INT, source, 1,
    //                         &status);
    ////                cerr << MPIGlobal::rank << ": received " << nelements
    ////                << " positions from " << source << endl;
    //                for(unsigned int j = 0; j < recvids.size(); j++){
    //                    _newlist.get_orphan(requests[source][j])
    //                            ->update_local_id(recvids[j]);
    //                }
    //            }
    //        }
    //        vector<MPI_Status> status((MPIGlobal::size-1)*2);
    //        MyMPI_Waitall((MPIGlobal::size-1)*2, &reqs[0], &status[0]);
    //        MyMPI_Barrier();
    //    }

    //    AdaptiveCellParallelSorter sorter;
    //    sorter.sort(_indices);

    //    vector<unsigned int> new_indices(_indices.size());
    //    vector<int> new_ranks(_indices.size());
    //    vector< vector<unsigned int> > sendvec(MPIGlobal::size);
    //    for(unsigned int i = 0; i < _indices.size(); i++){
    //        ArgHilbertObject* index = _indices[i];
    //        if(index->get_rank() == MPIGlobal::rank){
    //            new_indices[index->get_index()] = i;
    //            new_ranks[index->get_index()] = MPIGlobal::rank;
    //        } else {
    //            sendvec[index->get_rank()].push_back(i);
    //            sendvec[index->get_rank()].push_back(index->get_index());
    //        }
    //    }

    //    // communicate:
    //    // we need to send the sendvecs to the respective processes and add
    //    the
    //    // information to the local new_indices and new_ranks vectors
    //    if(MPIGlobal::size > 1){
    //        vector<MPI_Request> reqs((MPIGlobal::size-1), MPI_REQUEST_NULL);
    //        vector<unsigned int> buffers((MPIGlobal::size-1));
    //        unsigned int bufsize = MPIGlobal::sendsize/buffers.size();
    //        for(unsigned int i = 0; i < buffers.size(); i++){
    //            buffers[i] = i*bufsize;
    //        }

    //        unsigned int typesize = 2*sizeof(unsigned int);
    //        unsigned int maxsize = bufsize/typesize;
    //        if(!(bufsize%typesize)){
    //            maxsize--;
    //        }

    //        vector<int> allflag(MPIGlobal::size-1);
    //        vector<unsigned int> numsend(MPIGlobal::size-1);
    //        vector<unsigned int> sendsizes(MPIGlobal::size-1);
    //        int allflagsum = 0;
    //        for(int i = 0; i < MPIGlobal::size-1; i++){
    //            sendsizes[i] =
    //                    sendvec[(MPIGlobal::rank+1+i)%MPIGlobal::size].size()/2;

    //            int send_pos = 0;
    //            for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]);
    //                si++){
    //                MyMPI_Pack(&sendvec[(MPIGlobal::rank+1+i)%MPIGlobal::size][2*si],
    //                           2, MPI_UNSIGNED,
    //                           &MPIGlobal::sendbuffer[buffers[i]],
    //                           bufsize, &send_pos);
    //            }
    //            allflag[i] = (sendsizes[i] <= maxsize);
    //            // if allflagsum equals MPIGlobal::size-1, all data has been
    //            sent
    //            allflagsum += allflag[i];
    //            numsend[i] = 1;
    //            // add continuation signal
    //            MyMPI_Pack(&allflag[i], 1, MPI_INT,
    //                       &MPIGlobal::sendbuffer[buffers[i]], bufsize,
    //                       &send_pos);

    //            MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos,
    //                        MPI_PACKED, (MPIGlobal::rank+1+i)%MPIGlobal::size,
    //                        0,
    //                        &reqs[i]);
    //        }

    //        vector<unsigned int> numreceived(MPIGlobal::size-1, 0);
    //        for(int i = 0; i < (MPIGlobal::size-1); i++){
    //            MPI_Status status;
    //            for(int j = 0; j < MPIGlobal::size-1; j++){
    //                if(!allflag[j]){
    //                    int flag;
    //                    MyMPI_Test(&reqs[j], &flag, &status);
    //                    if(flag){
    //                        int send_pos = 0;
    //                        for(unsigned int si = numsend[j]*maxsize;
    //                            si < std::min((numsend[j]+1)*maxsize,
    //                            sendsizes[j]);
    //                            si++){
    //                            MyMPI_Pack(&sendvec[(MPIGlobal::rank+1+j)
    //                                       %MPIGlobal::size][2*si], 2,
    //                                       MPI_UNSIGNED,
    //                                       &MPIGlobal::sendbuffer[buffers[j]],
    //                                       bufsize, &send_pos);
    //                        }
    //                        allflag[j] = (sendsizes[j] <=
    //                        (numsend[j]+1)*maxsize);
    //                        allflagsum += allflag[j];
    //                        numsend[j]++;
    //                        // add continuation signal
    //                        MyMPI_Pack(&allflag[j], 1, MPI_INT,
    //                                   &MPIGlobal::sendbuffer[buffers[j]],
    //                                   bufsize,
    //                                   &send_pos);

    //                        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]],
    //                                    send_pos, MPI_PACKED,
    //                                    (MPIGlobal::rank+1+j)%MPIGlobal::size,
    //                                    0,
    //                                    &reqs[j]);
    //                    }
    //                }
    //            }

    //            int index;
    //            int tag;
    //            int recv_pos;

    //            // wait for a message from any source and with any tag to
    //            arrive
    //            MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
    //            // query the MPI_Status object to retrieve tag, source and
    //            size
    //            tag = status.MPI_TAG;
    //            index = status.MPI_SOURCE;
    //            int nelements;
    //            MyMPI_Get_count(&status, MPI_PACKED, &nelements);
    //            // select an index in the buffers to use for receiving and
    //            sending
    //            int freebuffer;
    //            if(index == MPIGlobal::size-1){
    //                freebuffer = MPIGlobal::rank;
    //            } else {
    //                freebuffer = index;
    //            }

    //            if(tag == 0){
    //                MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
    //                           nelements, MPI_PACKED, index, tag, &status);
    //                recv_pos = 0;
    //                unsigned int j = numreceived[freebuffer];
    //                while(recv_pos < nelements-4){
    //                    unsigned int newindex, oldindex;
    //                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
    //                                 nelements, &recv_pos, &newindex, 1,
    //                                 MPI_UNSIGNED);
    //                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
    //                                 nelements, &recv_pos, &oldindex, 1,
    //                                 MPI_UNSIGNED);
    //                    new_indices[oldindex] = newindex;
    //                    new_ranks[oldindex] = index;
    //                    j++;
    //                }
    //                numreceived[freebuffer] = j;
    //                int flag;
    //                MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
    //                             nelements, &recv_pos, &flag, 1, MPI_INT);

    //                if(!flag){
    //                    i--;
    //                }
    //            }
    //        }

    //        while(allflagsum < MPIGlobal::size-1){
    //            MPI_Status status;
    //            for(int j = 0; j < MPIGlobal::size-1; j++){
    //                if(!allflag[j]){
    //                    int flag;
    //                    MyMPI_Test(&reqs[j], &flag, &status);
    //                    if(flag){
    //                        int send_pos = 0;
    //                        for(unsigned int si = numsend[j]*maxsize;
    //                            si < std::min((numsend[j]+1)*maxsize,
    //                            sendsizes[j]);
    //                            si++){
    //                            MyMPI_Pack(&sendvec[(MPIGlobal::rank+1+j)
    //                                       %MPIGlobal::size][2*si], 2,
    //                                       MPI_UNSIGNED,
    //                                       &MPIGlobal::sendbuffer[buffers[j]],
    //                                       bufsize, &send_pos);
    //                        }
    //                        allflag[j] = (sendsizes[j] <=
    //                        (numsend[j]+1)*maxsize);
    //                        allflagsum += allflag[j];
    //                        numsend[j]++;
    //                        // add continuation signal
    //                        MyMPI_Pack(&allflag[j], 1, MPI_INT,
    //                                   &MPIGlobal::sendbuffer[buffers[j]],
    //                                   bufsize,
    //                                   &send_pos);

    //                        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]],
    //                                    send_pos, MPI_PACKED,
    //                                    (MPIGlobal::rank+1+j)%MPIGlobal::size,
    //                                    0,
    //                                    &reqs[j]);
    //                    }
    //                }
    //            }
    //        }

    //        vector<MPI_Status> status((MPIGlobal::size-1));
    //        MyMPI_Waitall((MPIGlobal::size-1), &reqs[0], &status[0]);
    //        MyMPI_Barrier();
    //    }

    //    cerr << "New indices:" << endl;
    //    for(unsigned int i = 0; i < new_indices.size(); i++){
    //        cerr << i << ": " << new_indices[i] << " (" << new_ranks[i] << ")"
    //             << endl;
    //    }

    //    // _indices now tells us where we should get the new data: cell i can
    //    be
    //    // found on the process with rank _indices[i]->get_rank() and has
    //    local id
    //    // _indices[i]->get_index()
    //    // Question: how can we use this information to update ngb information
    //    and
    //    // thus replace the old index of a ngb by its new value?
    //    // Answer: what we actually need is the reverse list, _rindices, where
    //    e.g.
    //    // _rindices[i]->get_new_index() gives us the new index of old index i
    //    // we could in principle make this list by sorting again on some kind
    //    of
    //    // generalized local index...

    //    // update the ngb information
    //    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
    //        _newlist.normalbegin(); it != _newlist.normalend(); ++it){
    //        it->update_local_id(new_indices);
    //        it->update_ngbs(new_indices);
    //    }
    //    // reorder the list
    //    _newlist.reorder(new_indices);

    //    if(MPIGlobal::size > 1){
    //        // update orphan information: the orphans do not change, but the
    //        ids of the
    //        // original particles do
    //        // Problem: this information is on the original processes: we have
    //        to
    //        // communicate it in some way. Since the original process does not
    //        know
    //        // which information to communicate, we have to request it
    //        first...
    //        vector< vector<unsigned int> > ngbrequest(MPIGlobal::size);
    //        vector< vector<unsigned int> > ngbanswer(MPIGlobal::size);
    //        vector< vector<double> > ngbpos(MPIGlobal::size);
    //        for(AdaptiveCellList<AdaptiveVorCell2d>::orphaniterator it =
    //            _newlist.orphanbegin(); it != _newlist.orphanend(); ++it){
    //            ngbrequest[it->get_rank()].push_back(it->get_original());
    //        }

    //        std::vector<MPI_Request> reqs((MPIGlobal::size-1)*2,
    //                                      MPI_REQUEST_NULL);
    //        std::vector<unsigned int> buffers((MPIGlobal::size-1)*2);
    //        unsigned int bufsize = MPIGlobal::sendsize/buffers.size();
    //        for(unsigned int i = 0; i < buffers.size(); i++){
    //            buffers[i] = i*bufsize;
    //        }

    //        unsigned int typesize = sizeof(unsigned int)+2*sizeof(double);
    //        unsigned int maxsize = bufsize/typesize;
    //        if(!(bufsize%typesize)){
    //            maxsize--;
    //        }

    //        // pack the local indices and send them to their respective
    //        processes
    //        vector<int> allflag(MPIGlobal::size-1);
    //        vector<unsigned int> numsend(MPIGlobal::size-1);
    //        vector<unsigned int> sendsizes(MPIGlobal::size-1);
    //        // we keep track of the total number of messages that should be
    //        sent
    //        // and received. We send and receive at least one message for
    //        every
    //        // process. For every message that is sent, an answer has to be
    //        // received. For incomplete sends, we need to send additional
    //        // messages. For incomplete receives, we have to receive
    //        additional
    //        // messages.
    //        int numtoreceive = MPIGlobal::size-1;
    //        int numtosend = MPIGlobal::size-1;
    //        int numrecv = 0;
    //        int numsent = 0;
    //        for(int i = 0; i < MPIGlobal::size-1; i++){
    //            sendsizes[i] =
    //                    ngbrequest[(MPIGlobal::rank+1+i)%MPIGlobal::size].size();

    //            int send_pos = 0;
    //            for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]);
    //                si++){
    //                MyMPI_Pack(&ngbrequest[(MPIGlobal::rank+1+i)%MPIGlobal::size][si],
    //                        1, MPI_UNSIGNED,
    //                        &MPIGlobal::sendbuffer[buffers[2*i]],
    //                        bufsize, &send_pos);
    //            }
    //            allflag[i] = (sendsizes[i] <= maxsize);
    //            if(!allflag[i]){
    //                numtosend++;
    //            }
    //            numsend[i] = 1;
    //            // add continuation signal
    //            MyMPI_Pack(&allflag[i], 1, MPI_INT,
    //                       &MPIGlobal::sendbuffer[buffers[2*i]], bufsize,
    //                       &send_pos);

    //            MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[2*i]], send_pos,
    //            MPI_PACKED,
    //                    (MPIGlobal::rank+1+i)%MPIGlobal::size, 0, &reqs[2*i]);
    //            numsent++;
    //            numtoreceive++;
    //        }

    //        // receive loop: receive 2 messages from every external process: a
    //        // buffer with indices and a buffer with new indices
    //        // We use MPI_Probe to detect arriving messages, independent of
    //        // their origin or tag. We then receive and process these messages
    //        // according to their tag.
    //        vector<unsigned int> numreceived(MPIGlobal::size-1, 0);
    //        while(numrecv < numtoreceive || numsent < numtosend){
    //            MPI_Status status;
    //            if(numsent < numtosend){
    //                for(int j = 0; j < MPIGlobal::size-1; j++){
    //                    if(!allflag[j]){
    //                        int flag;
    //                        MyMPI_Test(&reqs[j], &flag, &status);
    //                        if(flag){
    //                            int send_pos = 0;
    //                            for(unsigned int si = numsend[j]*maxsize;
    //                                si < std::min((numsend[j]+1)*maxsize,
    //                                              sendsizes[j]);
    //                                si++){
    //                                MyMPI_Pack(&ngbrequest[(MPIGlobal::rank+1+j)
    //                                        %MPIGlobal::size][si],
    //                                        1, MPI_UNSIGNED,
    //                                        &MPIGlobal::sendbuffer[buffers[2*j]],
    //                                        bufsize, &send_pos);
    //                            }
    //                            allflag[j] = (sendsizes[j] <=
    //                                          (numsend[j]+1)*maxsize);
    //                            if(!allflag[j]){
    //                                numtosend++;
    //                            }
    //                            numsend[j]++;
    //                            // add continuation signal
    //                            MyMPI_Pack(&allflag[j], 1, MPI_INT,
    //                                       &MPIGlobal::sendbuffer[buffers[2*j]],
    //                                    bufsize, &send_pos);

    //                            MyMPI_Isend(&MPIGlobal::sendbuffer
    //                                        [buffers[2*j]], send_pos,
    //                                    MPI_PACKED, (MPIGlobal::rank+1+j)%
    //                                    MPIGlobal::size, 0, &reqs[2*j]);
    //                            numsent++;
    //                            numtoreceive++;
    //                        }
    //                    }
    //                }
    //            }

    //            if(numrecv < numtoreceive){
    //                int index;
    //                int tag;
    //                int recv_pos;
    //                int send_pos;

    //                // wait for a message from any source and with any tag to
    //                // arrive
    //                MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
    //                numrecv++;
    //                // query the MPI_Status object to retrieve tag, source and
    //                // size
    //                tag = status.MPI_TAG;
    //                index = status.MPI_SOURCE;
    //                int nelements;
    //                MyMPI_Get_count(&status, MPI_PACKED, &nelements);
    //                // select an index in the buffers to use for receiving and
    //                // sending
    //                int freebuffer;
    //                if(index == MPIGlobal::size-1){
    //                    freebuffer = 2*MPIGlobal::rank;
    //                } else {
    //                    freebuffer = 2*index;
    //                }
    //                if(tag == 0){
    //                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
    //                            nelements, MPI_PACKED, index, 0, &status);
    //                    recv_pos = 0;
    //                    send_pos = 0;
    //                    while(recv_pos < nelements-4){
    //                        unsigned int oldindex;
    //                        MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
    //                                nelements, &recv_pos, &oldindex, 1,
    //                                MPI_UNSIGNED);
    //                        unsigned int newindex = new_indices[oldindex];
    //                        MyMPI_Pack(&newindex, 1, MPI_UNSIGNED,
    //                                   &MPIGlobal::sendbuffer[buffers[freebuffer+1]],
    //                                bufsize, &send_pos);
    //                        double newpos[2];
    //                        _newlist.get_normal(newindex)->get_position(newpos);
    //                        MyMPI_Pack(newpos, 2, MPI_DOUBLE,
    //                                   &MPIGlobal::sendbuffer[buffers[freebuffer+1]],
    //                                bufsize, &send_pos)
    //                    }
    //                    int flag;
    //                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
    //                            nelements, &recv_pos, &flag, 1, MPI_INT);
    //                    if(!flag){
    //                        numtoreceive++;
    //                    }

    //                    MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[freebuffer+1]],
    //                            send_pos, MPI_PACKED, index, 1,
    //                            &reqs[freebuffer+1]);
    //                }

    //                if(tag == 1){
    //                    MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer+1]],
    //                            nelements, MPI_PACKED, index, 1, &status);
    //                    unsigned int j = numreceived[freebuffer/2];
    //                    recv_pos = 0;
    //                    while(recv_pos < nelements){
    //                        unsigned int newindex;
    //                        MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer+1]],
    //                                nelements, &recv_pos, &newindex, 1,
    //                                MPI_UNSIGNED);
    //                        double newpos[2];
    //                        MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer+1]],
    //                                nelements, &recv_pos, newpos, 2,
    //                                MPI_DOUBLE);
    //                        ngbanswer[index].push_back(newindex);
    //                        ngbpos[index].push_back(newpos[0]);
    //                        ngbpos[index].push_back(newpos[1]);
    //                        j++;
    //                    }
    //                    numreceived[freebuffer/2] = j;
    //                }
    //            }
    //        }
    //        // We do not have to do an extra MPI_Test loop, since we know that
    //        // every original send corresponds to a receive above
    //        vector<MPI_Status> status((MPIGlobal::size-1)*2);
    //        MyMPI_Waitall((MPIGlobal::size-1)*2, &reqs[0], &status[0]);
    //        MyMPI_Barrier();

    //        cerr << "updating ghosts..." << endl;
    //        vector<unsigned int> ngbi(MPIGlobal::size, 0);
    //        for(AdaptiveCellList<AdaptiveVorCell2d>::orphaniterator it =
    //            _newlist.orphanbegin(); it != _newlist.orphanend(); ++it){
    //            cerr << it->get_original() << " --> ";
    //            it->update_local_id(ngbanswer[it->get_rank()][ngbi[it->get_rank()]]);
    //            it->update_ngbs(new_indices);
    //            //        double new_pos[2];
    //            //        new_pos[0] =
    //            ngbpos[it->get_rank()][2*ngbi[it->get_rank()]];
    //            //        new_pos[1] =
    //            ngbpos[it->get_rank()][2*ngbi[it->get_rank()]+1];
    //            //        it->move(new_pos);
    //            ngbi[it->get_rank()]++;
    //            cerr << it->get_original() << endl;
    //        }
    //    }

    ////    cerr << "Outputting..." << endl;
    ////    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
    ////        _newlist.normalbegin(); it != _newlist.normalend(); ++it){
    ////        AdaptiveVorCell2d* cell = *it;
    ////        cerr << "Cell " << cell->get_id() << endl;
    ////        double pos[2];
    ////        cell->get_position(pos);
    ////        cerr << pos[0] << "\t" << pos[1] << endl;
    ////        vector<int> ngbs = cell->get_ngbs();
    ////        vector<int> ngbwalls = cell->get_ngbwalls();
    ////        for(unsigned int i = 0; i < ngbs.size(); i++){
    ////            cerr << ngbs[i] << "\t" << ngbwalls[i] << endl;
    ////        }
    ////        cerr << endl;
    ////    }
    ////    cerr << "Done." << endl;

    //    stringstream fname;
    //    fname << "orphans_" << MPIGlobal::rank << ".dat";
    //    ofstream ofile(fname.str());
    //    print_tesselation_gnuplot(ofile);
    //    ofile.close();
    //    MyMPI_Barrier();

    //    cerr << "Sorted indices" << endl;
    //    for(unsigned int i = 0; i < _indices.size(); i++){
    //        cerr << i << ": " << _indices[i]->get_index() << " ("
    //             << _indices[i]->get_rank() << ")" << endl;
    //    }

    //    // delete indices to free memory
    //    for(unsigned int i = 0; i < _indices.size(); i++){
    //        delete _indices[i];
    //    }

    // sort the cells
    AdaptiveCellParallelSorter sorter;
    sorter.sort(_cells);

    // flag the particles for which we want accurate positions
    for(unsigned int i = 0; i < _cells.size(); i++) {
        vector<unsigned long> ngbs = _cells[i]->get_newngbs();
        for(unsigned int j = 0; j < ngbs.size(); j++) {
            // we do not specify a rank
            _newerlist.obtain_position(ngbs[j]);
        }
    }

    // update the position links
    _newerlist.reset_positions();
    _newerlist.calculate_positions(particles, _box12);

    //    MyMPI_Barrier();
    //    my_exit();

    unsigned int count = 0;
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        if(particles.gas(i)->get_endtime() == currentTime) {
            unsigned int index = particles.gas(i)->get_local_id();
            //            _newlist.get_normal(index)->activate();
            _cells[index]->activate();
            count++;
        }
    }
    // uncomment this to update the complete mesh
    //    for(unsigned int i = 0; i < _cells.size(); i++){
    //        _cells[i]->queue();
    //    }

    // we have 4 types of cells:
    //  - active cells with an active particle that have to be correct
    //  - ngbs of active cells, that have to be correct
    //  - ngbs of ngbs of active cells, that need to be processed to ensure the
    //    correctness of the ngbs of active cells
    //  - other cells, that can be safely left as they are
    //
    // Since a non-correct cell can lose or gain extra ngbs, we are only sure
    // which ngbs belong to category 2 after we processed category 1
    // The same goes for category 3.
    // Category 4 is not processed, that's the easy one.
    //
    // what we actually need is:
    // we activate all cells with an active particle
    // while(cells active){
    //   for(every cell){
    //     deactivate cell
    //     detect_crossovers:
    //       perform flips
    //     (re)activate all ngbs with an active particle involved
    //     semiactivate all ngbs with a passive particle (semi, because we do
    //                                                    not want to activate
    //                                                    their ngbs in turn)
    //   }
    //   perform flips for orphans -> communication
    //   (semi)activate orphan ngbs (involved)
    //   count active cells
    // }
    // do the same for semiactive cells
    // do the same for nonactive cells that are involved ngbs of semiactive
    // cells
    //
    // if we work with flags, then we need:
    //  - a flag to signal if a cell's particle is active
    //  - an active flag
    //  - a semiactive flag
    //  - a semisemiactive flag
    // we then need 3 mode iterations and for each iteration a while loop
    // with communication at the end of the iteration
    //
    // last problem we need to solve: orphans can only be handled once
    // the easiest (probably not the best) way to do this is by only handling
    // the orphans on the process with the smallest rank. But then we need
    // to determine the lowest rank for every flip.
    // Which means we need to restructure the flipping method to a cell
    // based method
    unsigned int numactive = count;
    while(numactive) {
        numactive = 0;
        vector<unsigned long> insertions;
        //        for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
        //            _newlist.normalbegin(); it != _newlist.normalend(); ++it){
        //            if(it->active()){
        //                if(it->detect_crossovers(it.index(), _newlist,
        //                insertions,
        //                                         _cuboid, _periodic)){
        //                    numactive++;
        //                } else {
        //                    it->deactivate();
        //                }
        //            }
        //        }
        for(unsigned int i = 0; i < _cells.size(); i++) {
            AdaptiveVorCell2d* cell = _cells[i];
            if(cell->active()) {
                if(cell->detect_crossovers(i, _cells, _newerlist, _cuboid,
                                           insertions, _periodic)) {
                    numactive++;
                } else {
                    cell->deactivate();
                }
            }
        }
        //        cerr << "Insertions size: " << insertions.size() << endl;
        //        for(unsigned int i = 0; i < insertions.size(); i+=4){
        //            cerr << insertions[i] << ": " << insertions[i+1] << ","
        //                 << insertions[i+2] << "," << insertions[i+3] << endl;
        //        }

        if(MPIGlobal::size == 1) {
            continue;
        }

        // communicate
        // we first have to find out to which processes we have to send
        vector<vector<unsigned long> > sendvec(MPIGlobal::size);
        for(unsigned int i = 0; i < insertions.size(); i += 4) {
            int rank = _newerlist.get_rank(insertions[i]);
            sendvec[rank].push_back(insertions[i]);
            sendvec[rank].push_back(insertions[i + 1]);
            sendvec[rank].push_back(insertions[i + 2]);
            sendvec[rank].push_back(insertions[i + 3]);
        }

        vector<MPI_Request> reqs(MPIGlobal::size - 1, MPI_REQUEST_NULL);
        vector<unsigned int> buffers(MPIGlobal::size - 1);
        unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
        for(unsigned int i = 0; i < buffers.size(); i++) {
            buffers[i] = i * bufsize;
        }

        unsigned int typesize = 4 * sizeof(unsigned long);
        unsigned int maxsize = bufsize / typesize;
        if(!(bufsize % typesize)) {
            maxsize--;
        }

        vector<int> allflag(MPIGlobal::size - 1);
        vector<unsigned int> numsend(MPIGlobal::size - 1);
        vector<unsigned int> sendsizes(MPIGlobal::size - 1);
        int allflagsum = 0;
        for(int i = 0; i < MPIGlobal::size - 1; i++) {
            sendsizes[i] =
                    sendvec[(MPIGlobal::rank + 1 + i) % MPIGlobal::size].size();

            int send_pos = 0;
            for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]);
                si += 4) {
                MyMPI_Pack(&sendvec[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                   [si],
                           4, MPI_UNSIGNED_LONG,
                           &MPIGlobal::sendbuffer[buffers[i]], bufsize,
                           &send_pos);
            }
            allflag[i] = (sendsizes[i] <= maxsize);
            // if allflagsum equals MPIGlobal::size-1, all data has been sent
            allflagsum += allflag[i];
            numsend[i] = 1;
            // add continuation signal
            MyMPI_Pack(&allflag[i], 1, MPI_INT,
                       &MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);

            MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos,
                        MPI_PACKED, (MPIGlobal::rank + 1 + i) % MPIGlobal::size,
                        0, &reqs[i]);
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
                            si <
                            std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                            si += 4) {
                            MyMPI_Pack(&sendvec[(MPIGlobal::rank + 1 + j) %
                                                MPIGlobal::size][si],
                                       4, MPI_UNSIGNED_LONG,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);
                        }
                        allflag[j] =
                                (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                        allflagsum += allflag[j];
                        numsend[j]++;
                        // add continuation signal
                        MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);

                        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]],
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
            // select an index in the buffers to use for receiving and sending
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
                unsigned int j = numreceived[freebuffer];
                //                cerr << "received:" << endl;
                while(recv_pos < nelements - 4) {
                    unsigned long ids[4];
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, ids, 4,
                                 MPI_UNSIGNED_LONG);
                    //                    cerr << ids[0] << ": " << ids[1] <<
                    //                    "," << ids[2] << ","
                    //                                   << ids[3] << endl;
                    unsigned int locid = _newerlist.get_locid(ids[0]);
                    AdaptiveVorCell2d* cell = _cells[locid];
                    cell->add_newngb(ids[1], ids[2], ids[3]);
                    //                    cell->complete(_newerlist);
                    // make sure we have the position of the particle
                    _newerlist.obtain_position(ids[3]);
                    cell->activate();
                    numactive++;
                    j++;
                }
                numreceived[freebuffer] = j;
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
                            si <
                            std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                            si += 4) {
                            MyMPI_Pack(&sendvec[(MPIGlobal::rank + 1 + j) %
                                                MPIGlobal::size][si],
                                       4, MPI_UNSIGNED_LONG,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);
                        }
                        allflag[j] =
                                (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                        allflagsum += allflag[j];
                        numsend[j]++;
                        // add continuation signal
                        MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);

                        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]],
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
        unsigned int lnumactive = numactive;
        MyMPI_Allreduce(&lnumactive, &numactive, 1, MPI_INT, MPI_SUM);
        // make sure the positions of newly added copies are in the list
        _newerlist.calculate_positions(particles, _box12);

        // communication:
        // if an orphan was involved in a face flip, it might or might not have
        // changed on the original process, depending on whether or not another
        // ngb that was involved is also on the original process
        // more clearly stated: there are 4 cells involved in the face flip and
        // the flip can be detected in two of them. If two cells are on process
        // A and two on B, then the face flip will be detected and performed
        // on both processes. If three cells are on A and one B, it might not
        // if the flip cannot be detected in the cell on B.
        // The same holds if two cells are on A, one on B and one on C. Then
        // on either B or C the face flip will not be detected for sure. If
        // one cell is on A, B, C and D, then one two of these processes no
        // flip is detected...
        //
        // In the following step, we have to somehow communicate information
        // between the processes to fix this.
        // Two options:
        //  - either we somehow accumulate the necessary information locally
        //   or
        //  - we communicate the new state for all orphans
        //
        // One sounds tricky, but two sounds completely impossible. So let's
        // think about one.
        //
        // We only need to communicate if we can detect something that will
        // not be detected on the other process. This can be checked, since
        // we know which of the 4 cells can detect the flip and we also know
        // on which process each orphan resides.
        // This also means that we are sure that what we will communicate will
        // be new for this process and not a duplicate of something that already
        // happened locally. An exception is when we have four cells on four
        // different processes, since then two processes can send information
        // to two other processes. It seems reasonable to only send information
        // from the process with the lowest/highest rank in this case.
        // So we can guarantee that we only send necessary information and that
        // this information will be received only once and will be useful.
        //
        // What do we need to send?
        // Since the face flip is detected in the cells where the face has to
        // be removed, we will never have to communicate ngb removals. We hence
        // have to communicate ngb insertions, which require 5 pieces of
        // information:
        // (1) the id of the cell that gets the ngb
        // (2)+(3) the ids of the ngbs that surround the new ngb
        // (4) the id of the new ngb
        // (5) the wallpos of the new ngb
        // Note that (2) and (3) can be replaced by a single id if we impose a
        // strict ordering. (4) WILL be a new orphan that has to be created on
        // the receiving process. (5) is easier than you think.
        // (2)/(3) are difficult, since we only have local ngb ids, and not the
        // ids of the orphans that will be stored on the local process. So this
        // requires some backtracing: we send the local ids and compare these
        // to the local ids stored in the orphans on the other process. This
        // is slightly more expensive than just knowing what you look for, but
        // since we only have to search the ngb list of a single cell, it is not
        // that hard a problem...
        //        if(MPIGlobal::size > 1){
        //            // (1) fill buffers
        //            vector< vector<int> > sendinfo(MPIGlobal::size);
        //            vector< vector<double> > sendpos(MPIGlobal::size);
        //            for(unsigned int i = 0; i < insertions.size()/6; i++){
        //                int sendto = insertions[6*i];
        //                // the local id of the original cell
        //                sendinfo[sendto].push_back(insertions[6*i+1]);
        //                // the ngb ids
        //                sendinfo[sendto].push_back(insertions[6*i+2]);
        //                sendinfo[sendto].push_back(insertions[6*i+3]);
        //                // the ngbwall of the new ngb
        //                sendinfo[sendto].push_back(insertions[6*i+5]);
        //                // the information needed to make a new orphan:
        //                //  - original id
        //                //  - position
        //                sendinfo[sendto].push_back(insertions[6*i+4]);
        //                double pos[2];
        //                _newlist.get_normal(insertions[6*i+4])->get_position(pos);
        //                sendpos[sendto].push_back(pos[0]);
        //                sendpos[sendto].push_back(pos[1]);
        //            }

        //            std::vector<MPI_Request> reqs(MPIGlobal::size*2,
        //            MPI_REQUEST_NULL);
        //            for(int i = 0; i < MPIGlobal::size; i++){
        //                if(i != MPIGlobal::rank){
        //                    MyMPI_Isend(&sendinfo[i][0], sendinfo[i].size(),
        //                    MPI_INT, i,
        //                              0, &reqs[2*i]);
        //                    MyMPI_Isend(&sendpos[i][0], sendpos[i].size(),
        //                    MPI_DOUBLE,
        //                                i, 1, &reqs[2*i+1]);
        //                }
        //            }

        //            for(int i = 0; i < MPIGlobal::size-1; i++){
        //                MPI_Status status;
        //                MyMPI_Probe(MPI_ANY_SOURCE, 0, &status);
        //                int source = status.MPI_SOURCE;
        //                int nelements;
        //                MyMPI_Get_count(&status, MPI_INT, &nelements);
        //                vector<int> recvinfo(nelements);
        //                vector<double> recvpos(nelements*2/5);
        //                MyMPI_Recv(&recvinfo[0], nelements, MPI_INT, source,
        //                0,
        //                         &status);
        //                MyMPI_Recv(&recvpos[0], nelements*2/5, MPI_DOUBLE,
        //                source, 1,
        //                         &status);
        //                for(int j = 0; j < nelements/5; j++){
        //                    int in_between[2] = {-1, -1};
        //                    vector<int> ngbs =
        //                            _newlist.get_normal(recvinfo[5*j])->get_ngbs();
        //                    vector<int> ngbwalls =
        //                            _newlist.get_normal(recvinfo[5*j])->get_ngbwalls();
        //                    for(unsigned int k = 0; k < ngbs.size(); k++){
        //                        if(_newlist.is_orphan(ngbs[k])){
        //                            if((int)_newlist.get_orphan(ngbs[k])->get_original()
        //                                    == recvinfo[5*j+1]){
        //                                in_between[0] = k;
        //                            }
        //                            if((int)_newlist.get_orphan(ngbs[k])->get_original()
        //                                    == recvinfo[5*j+2]){
        //                                in_between[1] = k;
        //                            }
        //                        }
        //                    }
        //                    // in the category 'filthy tricks I hope no one
        //                    will ever
        //                    // notice':
        //                    // if for some reason one of the ngbs was already
        //                    removed
        //                    // locally (because there is an overlap of changes
        //                    for a
        //                    // cell), then we can skip this insertion
        //                    // this seems to work, so no questions asked...
        //                    if(in_between[0] == -1 || in_between[1] == -1){
        //                        continue;
        //                    }
        //                    Vec pos(recvpos[2*j], recvpos[2*j+1]);
        //                    AdaptiveVorCell2d* orphan =
        //                            new AdaptiveVorCell2d(pos, 0,
        //                            recvinfo[5*j+4]);
        //                    // this orphan needs ngbs!
        //                    // this is quite easy, since the ngbs are
        //                    in_between [0] and
        //                    // [1] and the cell that receives the orphan
        //                    // walls are a different problem, but I assume it
        //                    should be
        //                    // possible to deduct them from the received wall
        //                    and the
        //                    // ngbwalls of the in_betweens in the cell that
        //                    receives
        //                    orphan->set_rank(source);
        //                    int orphanid = _newlist.add_orphan_cell(orphan);
        //                    // add the ngbs to the orphan
        //                    if(in_between[0] ==
        //                    (int)((in_between[1]+1)%ngbs.size())){
        //                        orphan->add_ngb(ngbs[in_between[0]]);
        //                        orphan->add_ngbwall(
        //                                    AdaptiveMeshUtils::get_wallpos(
        //                                        recvinfo[5*j+3],
        //                                        ngbwalls[in_between[0]]
        //                                        )
        //                                    );
        //                        orphan->add_ngb(recvinfo[5*j]);
        //                        orphan->add_ngbwall(
        //                                    AdaptiveMeshUtils::get_wallpos(
        //                                        recvinfo[5*j+3], 0
        //                                        )
        //                                    );
        //                        orphan->add_ngb(ngbs[in_between[1]]);
        //                        orphan->add_ngbwall(
        //                                    AdaptiveMeshUtils::get_wallpos(
        //                                        recvinfo[5*j+3],
        //                                        ngbwalls[in_between[1]]
        //                                        )
        //                                    );
        //                    } else {
        //                        orphan->add_ngb(ngbs[in_between[1]]);
        //                        orphan->add_ngbwall(
        //                                    AdaptiveMeshUtils::get_wallpos(
        //                                        recvinfo[5*j+3],
        //                                        ngbwalls[in_between[1]]
        //                                        )
        //                                    );
        //                        orphan->add_ngb(recvinfo[5*j]);
        //                        orphan->add_ngbwall(
        //                                    AdaptiveMeshUtils::get_wallpos(
        //                                        recvinfo[5*j+3], 0
        //                                        )
        //                                    );
        //                        orphan->add_ngb(ngbs[in_between[0]]);
        //                        orphan->add_ngbwall(
        //                                    AdaptiveMeshUtils::get_wallpos(
        //                                        recvinfo[5*j+3],
        //                                        ngbwalls[in_between[0]]
        //                                        )
        //                                    );
        //                    }
        //                    _newlist.get_normal(recvinfo[5*j])
        //                            ->add_ngb(ngbs[in_between[0]],
        //                            ngbs[in_between[1]],
        //                                      orphanid, recvinfo[5*j+3]);
        //                    _newlist.get_normal(recvinfo[5*j])->complete(_newlist,
        //                                                                 _cuboid,
        //                                                                 _periodic);
        //                }
        //            }
        //            vector<MPI_Status> status((MPIGlobal::size-1)*2);
        //            MyMPI_Waitall((MPIGlobal::size-1)*2, &reqs[0],
        //            &status[0]);
        //        }
    }
// implement the code below to support inactive particles
#ifdef EENCONSTANTE
    numactive = 0;
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        if(particles.gas(i)->get_endtime() == currentTime) {
            unsigned int index = particles.gas(i)->get_local_id();
            vector<int> ngbs = _newlist.get_normal(index)->get_ngbs();
            for(unsigned int j = 0; j < ngbs.size(); j++) {
                if(_newlist.is_normal(ngbs[j]) &&
                   particles.gas(ngbs[j])->get_endtime() != currentTime) {
                    if(!_newlist.get_normal(ngbs[j])->active()) {
                        _newlist.get_normal(ngbs[j])->activate();
                        _newlist.get_normal(ngbs[j])->semiactivate();
                        numactive++;
                    }
                }
            }
        }
    }
    while(numactive) {
        numactive = 0;
        vector<int> insertions;
        for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                    _newlist.normalbegin();
            it != _newlist.normalend(); ++it) {
            if(it->active()) {
                if(it->detect_crossovers(it.index(), _newlist, insertions,
                                         _cuboid, _periodic)) {
                    numactive++;
                } else {
                    it->deactivate();
                }
            }
        }
    }
    numactive = 0;
    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                _newlist.normalbegin();
        it != _newlist.normalend(); ++it) {
        if(it->semiactive()) {
            vector<int> ngbs = it->get_ngbs();
            for(unsigned int j = 0; j < ngbs.size(); j++) {
                if(_newlist.is_normal(ngbs[j]) &&
                   particles.gas(ngbs[j])->get_endtime() != currentTime) {
                    if(!_newlist.get_normal(ngbs[j])->active() &&
                       !_newlist.get_normal(ngbs[j])->semiactive()) {
                        _newlist.get_normal(ngbs[j])->activate();
                        numactive++;
                    }
                }
            }
        }
    }
    while(numactive) {
        numactive = 0;
        vector<int> insertions;
        for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                    _newlist.normalbegin();
            it != _newlist.normalend(); ++it) {
            if(it->active()) {
                if(it->detect_crossovers(it.index(), _newlist, insertions,
                                         _cuboid, _periodic)) {
                    numactive++;
                } else {
                    it->deactivate();
                }
            }
        }
    }
#endif

    //    for(unsigned int mode = 0; mode < 3; mode++){
    //        unsigned int numactive = 1;
    //        while(numactive){
    //            numactive = 0;
    //        for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
    //    _newlist.normalbegin(); it != _newlist.normalend(); ++it){
    //            unsigned int index = it.index();
    //            if(mode){
    //                if(mode == 1){
    //                    if(!it->semiactive() || it->active()){
    //                        continue;
    //                    }
    //                } else {
    //                    if(!it->queued() || it->active() || it->semiactive()){
    //                        continue;
    //                    }
    //                }
    //            } else {
    //                if(!it->active()){
    //                    continue;
    //                }
    //            }
    //            it->complete(_newlist, _cuboid, _periodic);
    //            if(it->detect_crossovers(index, _newlist, true,
    //    _cuboid, _periodic)){
    ////                it->complete(_newlist, _cuboid, _periodic);
    //                numactive++;
    //            }
    //            if(!mode){
    //                vector<int> ngbs = it->get_ngbs();
    //                // flag all neighbours to be moved
    //                for(unsigned int j = 0; j < ngbs.size(); j++){
    //                    if(_newlist.is_normal(ngbs[j])){
    //                        _newlist.get_normal(ngbs[j])->semiactivate();
    //                    }
    //                }
    //            } else {
    //                if(mode == 1){
    //                    vector<int> ngbs = it->get_ngbs();
    //                    // flag all neighbours to be moved
    //                    for(unsigned int j = 0; j < ngbs.size(); j++){
    //                        if(_newlist.is_normal(ngbs[j])){
    //                            _newlist.get_normal(ngbs[j])->queue();
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //        }
    //    }

    //#warning Turn this off for production simulations!!
    //    if(count == _cells.size()){
    //        cout << "Total volume: " << get_total_volume() << endl;
    //    }

    //    if(currentTime == 348606758656147456){
    //        ofstream tfile("cell910_0.dat");
    //        _cells[910]->print(tfile, 910);
    //    }

    //    if(currentTime == 348747496144502784){
    //        ofstream tfile("cell910_1.dat");
    //        _cells[910]->print(tfile, 910);
    //    }

    //    unsigned int j = 0;
    //    while(_newlist[j]->get_particle()->id() != 4202){
    //        j++;
    //    }
    //    ofstream tfile("cell4202.dat");
    //    _newlist[j]->print(tfile, 4202);
    //    tfile.close();
    //    exit(1);

    return count;
}

// void AdaptiveVorTess2d::save_restart(unsigned int nr){
//    stringstream name;
//    name << "restart";
//    name.fill('0');
//    name.width(3);
//    name << nr;
//    name << ".dat";
//    ofstream rfile(name.str().c_str(), ios::binary | ios::out);
//    unsigned int vsize = _cells.size();
//    rfile.write((char*) &vsize, sizeof(unsigned int));
//    for(unsigned int i = 0; i < vsize; i++){
//        _cells[i]->save_restart(rfile);
//    }
//    vsize = _ghosts.size();
//    rfile.write((char*) &vsize, sizeof(unsigned int));
//    for(unsigned int i = 0; i < vsize; i++){
//        _ghosts[i]->save_restart(rfile);
//    }
//}

/**
 * @brief Get the volume of the cell with the given index
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @return The volume of the cell
 */
double AdaptiveVorTess2d::get_volume(unsigned int i) {
    //    return _newlist.get_normal(i)->get_volume();
    return _cells[i]->get_volume();
}

/**
 * @brief Get the characteristic radius of the cell with the given index
 *
 * The characteristic radius is the radius of a disc with the same area as the
 * cell.
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @return The characteristic radius of the cell
 */
double AdaptiveVorTess2d::get_h(unsigned int i) {
    //    return _newlist.get_normal(i)->get_h();
    return _cells[i]->get_h();
}

/**
 * @brief Calculate a velocity for the cell with the given index, based on the
 * hydrodynamical variables of the given GasParticle
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @param particle GasParticle associated with the given cell
 * @return The velocity for the mesh generator of the cell
 */
Vec AdaptiveVorTess2d::get_velocity(unsigned int i, GasParticle* particle) {
    //    return _newlist.get_normal(i)->get_velocity(particle);
    return _cells[i]->get_velocity(particle);
}

/**
 * @brief Calculate the centroid of the cell with the given index
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @return Coordinates of the centroid of the cell
 */
Vec AdaptiveVorTess2d::get_centroid(unsigned int i) {
#if ndim_ == 2
    double centroid[2];
    //    _newlist.get_normal(i)->get_centroid(centroid);
    _cells[i]->get_centroid(centroid);
    return Vec(centroid[0], centroid[1]);
#else
    return Vec(0., 0., 0.);
#endif
}

/**
 * @brief Estimate gradients in the cell with the given index, based on the
 * hydrodynamical variables of the given GasParticle
 *
 * @param i unsigned integer index of a cell in the tesselation
 * @param delta Array to store the results in
 * @param particle GasParticle associated with the cell
 */
void AdaptiveVorTess2d::estimate_gradients(unsigned int i, StateVector* delta,
                                           GasParticle* particle) {
    //    _newlist.get_normal(i)->estimate_gradients(delta, particle, _newlist,
    //                                               _cuboid, _periodic);
    _cells[i]->estimate_gradients(delta, particle, _newerlist, _periodic);
}

/**
 * @brief Get a list of AdaptiveFace2d faces for which fluxes should be
 * calculated
 *
 * Only faces that are active, i.e. have at least one active neighbouring cell,
 * are returned.
 *
 * @param currentTime unsigned long integer current simulation time
 * @return std::vector of AdaptiveFace2d instances that should be used to do the
 * flux exchange for the hydrodynamical integration
 */
vector<AdaptiveFace2d*> AdaptiveVorTess2d::get_faces(
        unsigned long currentTime) {
#ifdef OUD
    vector<AdaptiveFace2d*> faces;
    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                _newlist.normalbegin();
        it != _newlist.normalend(); ++it) {
        vector<int> ngbs = it->get_ngbs();
        bool active = (it->get_particle()->get_starttime() == currentTime);
        for(unsigned int j = 0; j < ngbs.size(); j++) {
            // only construct the faces once
            if(it.index() < ngbs[j]) {
                if(_newlist.is_normal(ngbs[j])) {
                    if(!active &&
                       _newlist.get_normal(ngbs[j])
                                       ->get_particle()
                                       ->get_starttime() != currentTime) {
                        // both cells are inactive; we don't add the face
                        continue;
                    }
                    double midface[2];
                    it->get_face(j)->get_midpoint(midface);
                    double area = it->get_face(j)->get_area();
                    Vec pos;
                    if(_periodic) {
                        double posvec[2];
                        it->get_periodic_position(j, _newlist, posvec, _cuboid);
                        pos.set(posvec[0], posvec[1]);
                    } else {
                        pos = _newlist.get_normal(ngbs[j])
                                      ->get_particle()
                                      ->get_position();
                    }
                    faces.push_back(new AdaptiveFace2d(
                            it->get_particle(),
                            _newlist.get_normal(ngbs[j])->get_particle(), pos,
                            midface, area));
                } else {
                    if(!active) {
                        continue;
                    }
                    if(_periodic) {
                        // this does never happen. Right?
                        //                        if(ngbs[j]-_cells.size() < i){
                        //                            continue;
                        //                        }
                    }
                    double midface[2];
                    it->get_face(j)->get_midpoint(midface);
                    double area = it->get_face(j)->get_area();
                    AdaptiveVorCell2d* ngb = _newlist.get_ghost(ngbs[j]);
                    GasParticle* right = NULL;
                    if(_periodic) {
                        // this cannot happen...
                        right = ngb->get_particle();
                    }
                    double ngbposvec[2];
                    ngb->get_ghost_position(ngbposvec, _cuboid);
#if ndim_ == 2
                    Vec ngbpos(ngbposvec[0], ngbposvec[1]);
#else
                    Vec ngbpos;
#endif
                    faces.push_back(new AdaptiveFace2d(
                            it->get_particle(), right, ngbpos, midface, area));
                }
            }
        }
    }
    return faces;
#else
    vector<AdaptiveFace2d*> faces;
    for(unsigned int i = 0; i < _cells.size(); i++) {
        vector<unsigned long> ngbs = _cells[i]->get_newngbs();
        //        bool active = (it->get_particle()->get_starttime() ==
        //        currentTime);
        bool active = _newerlist.get_starttime(_cells[i]->get_particleID());
        for(unsigned int j = 0; j < ngbs.size(); j++) {
            // only construct the faces once
            if(_cells[i]->get_particleID() < ngbs[j]) {
                if(!active &&
                   _newerlist.get_starttime(ngbs[j]) != currentTime) {
                    // both cells are inactive; we don't add the face
                    continue;
                }
                double midface[2];
                _cells[i]->get_face(j)->get_midpoint(midface);
                double area = _cells[i]->get_face(j)->get_area();
                Vec pos = _newerlist.get_position(ngbs[j]);
                // we need to fill in the neighbours!
                faces.push_back(
                        new AdaptiveFace2d(NULL, NULL, pos, midface, area));
            }
        }
    }
    return faces;
#endif
}

/**
 * @brief Calculate the total volume of all cells in the tesselation
 *
 * This can be used as a check for the validity of the total tesselation. For a
 * valid tesselation, the total volume equals the volume of the simulation box.
 *
 * @return The total volume of all cells in the tesselation
 */
double AdaptiveVorTess2d::get_total_volume() {
    double totvol = 0.;
    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                _newlist.normalbegin();
        it != _newlist.normalend(); ++it) {
        totvol += it->get_volume();
    }
    return totvol;
}

/**
 * @brief Update the positions of the generators of the cells with the given new
 * values
 *
 * @param particles std::vector of GasParticle instances whose positions
 * correspond to the new positions of the mesh generators
 */
void AdaptiveVorTess2d::update_positions(ParticleVector& particles) {
    //    cerr << "Local particles (real):" << endl;
    //    for(unsigned int i = 0; i < particles.gassize(); i++){
    //        cerr << particles.gas(i)->id() << endl;
    //    }
    //    cerr << "Local particles (list):" << endl;
    //    for(unsigned long id = 0; id < 10; id++){
    //        if(_newerlist.is_local(id)){
    //            cerr << id << endl;
    //        }
    //    }
    // keep particles inside the periodic box
    vector<vector<unsigned long> > commids(MPIGlobal::size);
    vector<vector<unsigned long> > commngbs(MPIGlobal::size);
    vector<vector<int> > commwalls(MPIGlobal::size);
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        // some magic
        double pos[2];
        pos[0] = particles.gas(i)->get_position().x();
        pos[1] = particles.gas(i)->get_position().y();
        int wallpos = AdaptiveMeshUtils::trim_wallpos(pos, _cuboid);
        particles.gas(i)->set_x(pos[0]);
        particles.gas(i)->set_y(pos[1]);
        //        cerr << "particle:" << endl;
        //        cerr << particles.gas(i)->id() << " (" <<
        //        _newerlist.is_local(particles.gas(i)->id()) << ")" << endl;
        if(wallpos) {
            // the particle has moved trough a periodic wall, we have to update
            // its neighbours
            // also the neighbours on other processes...
            //            AdaptiveVorCell2d *cell = _newlist.get_normal(i);
            AdaptiveVorCell2d* cell = _cells[i];
            vector<unsigned long>& newngbs = cell->get_newngbs();
            for(unsigned int j = 0; j < newngbs.size(); j++) {
                //                cerr << newngbs[j] << " - " <<
                //                _newerlist.get_id(newngbs[j]) << " (" <<
                //                _newerlist.is_local(newngbs[j]) << ")" <<
                //                endl;
                if(_newerlist.is_local(newngbs[j])) {
                    unsigned int ingb = _newerlist.get_locid(newngbs[j]);
                    //                    AdaptiveVorCell2d *ngb =
                    //                    _newlist.get_normal(ingb);
                    AdaptiveVorCell2d* ngb = _cells[ingb];
                    vector<unsigned long>& ngbngbs = ngb->get_newngbs();
                    for(unsigned int k = 0; k < ngbngbs.size(); k++) {
                        if(_newerlist.is_copy(ngbngbs[k],
                                              particles.gas(i)->id())) {
                            int oldwall = _newerlist.get_wall(ngbngbs[k]);
                            //                            cerr << ngbngbs[k] <<
                            //                            "\t" << newngbs[j] <<
                            //                            endl;
                            int newwall = AdaptiveMeshUtils::replace_wallpos(
                                    oldwall, wallpos);
                            unsigned long newid;
                            if(newwall) {
                                newid = _newerlist.add_ghost(
                                        particles.gas(i)->id(), newwall,
                                        MPIGlobal::rank);
                            } else {
                                newid = particles.gas(i)->id();
                            }
                            ngbngbs[k] = newid;
                        }
                    }
                } else {
                    int rank = _newerlist.get_rank(newngbs[j]);
                    if(rank < 0) {
                        cerr << "Negative rank. Oops..." << endl;
                        cerr << newngbs[j] << endl;
                        cerr << rank << endl;
                        my_exit();
                    }
                    commids[rank].push_back(newngbs[j]);
                    commngbs[rank].push_back(particles.gas(i)->id());
                    commwalls[rank].push_back(wallpos);
                }
                // change the ngb
                int oldwall = _newerlist.get_wall(newngbs[j]);
                int newwall = AdaptiveMeshUtils::swap_wallpos(oldwall, wallpos);
                unsigned long newid;
                if(newwall) {
                    newid = _newerlist.add_ghost(_newerlist.get_id(newngbs[j]),
                                                 newwall, MPIGlobal::rank);
                } else {
                    newid = _newerlist.get_id(newngbs[j]);
                }
                newngbs[j] = newid;
            }
        }
    }

    if(MPIGlobal::size > 1) {
        // communicate (we only need to send and receive once)
        vector<MPI_Request> reqs(MPIGlobal::size - 1, MPI_REQUEST_NULL);
        vector<unsigned int> buffers(MPIGlobal::size - 1);
        unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
        for(unsigned int i = 0; i < buffers.size(); i++) {
            buffers[i] = i * bufsize;
        }

        unsigned int typesize = 2 * sizeof(unsigned long) + sizeof(int);
        unsigned int maxsize = bufsize / typesize;
        if(!(bufsize % typesize)) {
            maxsize--;
        }

        vector<int> allflag(MPIGlobal::size - 1);
        vector<unsigned int> numsend(MPIGlobal::size - 1);
        vector<unsigned int> sendsizes(MPIGlobal::size - 1);
        int allflagsum = 0;
        for(int i = 0; i < MPIGlobal::size - 1; i++) {
            sendsizes[i] =
                    commids[(MPIGlobal::rank + 1 + i) % MPIGlobal::size].size();

            int send_pos = 0;
            for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]);
                si++) {
                MyMPI_Pack(&commids[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                   [si],
                           1, MPI_UNSIGNED_LONG,
                           &MPIGlobal::sendbuffer[buffers[i]], bufsize,
                           &send_pos);
                MyMPI_Pack(
                        &commngbs[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                 [si],
                        1, MPI_UNSIGNED_LONG,
                        &MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);
                MyMPI_Pack(
                        &commwalls[(MPIGlobal::rank + 1 + i) % MPIGlobal::size]
                                  [si],
                        1, MPI_INT, &MPIGlobal::sendbuffer[buffers[i]], bufsize,
                        &send_pos);
            }
            allflag[i] = (sendsizes[i] <= maxsize);
            // if allflagsum equals MPIGlobal::size-1, all data has been sent
            allflagsum += allflag[i];
            numsend[i] = 1;
            // add continuation signal
            MyMPI_Pack(&allflag[i], 1, MPI_INT,
                       &MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);

            MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos,
                        MPI_PACKED, (MPIGlobal::rank + 1 + i) % MPIGlobal::size,
                        0, &reqs[i]);
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
                            si <
                            std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                            si++) {
                            MyMPI_Pack(&commids[(MPIGlobal::rank + 1 + j) %
                                                MPIGlobal::size][si],
                                       1, MPI_UNSIGNED_LONG,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);
                            MyMPI_Pack(&commngbs[(MPIGlobal::rank + 1 + j) %
                                                 MPIGlobal::size][si],
                                       1, MPI_UNSIGNED_LONG,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);
                            MyMPI_Pack(&commwalls[(MPIGlobal::rank + 1 + j) %
                                                  MPIGlobal::size][si],
                                       1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);
                        }
                        allflag[j] =
                                (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                        allflagsum += allflag[j];
                        numsend[j]++;
                        // add continuation signal
                        MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);

                        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]],
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
            // select an index in the buffers to use for receiving and sending
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
                unsigned int j = numreceived[freebuffer];
                while(recv_pos < nelements - 4) {
                    unsigned long id;
                    unsigned long ngbid;
                    int wall;
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &id, 1,
                                 MPI_UNSIGNED_LONG);
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &ngbid, 1,
                                 MPI_UNSIGNED_LONG);
                    MyMPI_Unpack(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                 nelements, &recv_pos, &wall, 1, MPI_INT);
                    unsigned int ingb = _newerlist.get_locid(id);
                    //                AdaptiveVorCell2d *ngb =
                    //                _newlist.get_normal(ingb);
                    AdaptiveVorCell2d* ngb = _cells[ingb];
                    vector<unsigned long>& ngbngbs = ngb->get_newngbs();
                    for(unsigned int k = 0; k < ngbngbs.size(); k++) {
                        if(_newerlist.is_copy(ngbngbs[k], ngbid)) {
                            int oldwall = _newerlist.get_wall(ngbngbs[k]);
                            int newwall = AdaptiveMeshUtils::replace_wallpos(
                                    oldwall, wall);
                            unsigned long newid;
                            if(newwall) {
                                newid = _newerlist.add_ghost(ngbid, newwall,
                                                             index);
                            } else {
                                newid = ngbid;
                            }
                            ngbngbs[k] = newid;
                        }
                    }
                }
                numreceived[freebuffer] = j;
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
                            si <
                            std::min((numsend[j] + 1) * maxsize, sendsizes[j]);
                            si++) {
                            MyMPI_Pack(&commids[(MPIGlobal::rank + 1 + j) %
                                                MPIGlobal::size][si],
                                       1, MPI_UNSIGNED_LONG,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);
                            MyMPI_Pack(&commngbs[(MPIGlobal::rank + 1 + j) %
                                                 MPIGlobal::size][si],
                                       1, MPI_UNSIGNED_LONG,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);
                            MyMPI_Pack(&commwalls[(MPIGlobal::rank + 1 + j) %
                                                  MPIGlobal::size][si],
                                       1, MPI_INT,
                                       &MPIGlobal::sendbuffer[buffers[j]],
                                       bufsize, &send_pos);
                        }
                        allflag[j] =
                                (sendsizes[j] <= (numsend[j] + 1) * maxsize);
                        allflagsum += allflag[j];
                        numsend[j]++;
                        // add continuation signal
                        MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                   &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                   &send_pos);

                        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[j]],
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

    _newerlist.reset_positions();
    _newerlist.calculate_positions(particles, _box12);

    //    stringstream oname;
    //    oname << "after_update_" << MPIGlobal::rank << ".txt";
    //    ofstream ofile(oname.str());
    //    for(auto it = _newlist.normalbegin(); it != _newlist.normalend();
    //    ++it){
    //        it->complete(_newerlist);
    //        it->print(ofile, it.index());
    //    }
    //    ofile.close();
    //    MyMPI_Barrier();
    ////    my_exit();

    //    for(unsigned int i = 0; i < particles.gassize(); i++){
    //        double pos[2];
    //        pos[0] = particles.gas(i)->x();
    //        pos[1] = particles.gas(i)->y();
    //        AdaptiveVorCell2d* cell =
    //                _newlist.get_normal(particles.gas(i)->get_local_id());
    //        cell->move(pos);
    //        cell->reset_flags();
    //    }

    //    if(MPIGlobal::size > 1){
    //        // update the positions of the orphans
    //        // we do this in tree steps:
    //        // first, we loop over all orphans and obtain the process and
    //        local_id
    //        // on that process of the real particle. We store these values in
    //        // separate buffers for every process.
    //        // then, we send a request to all other processes containing the
    //        // local_ids we found. Every other process then stores the
    //        positions
    //        // we requested and sends them back to the requesting process
    //        // at last, we use the resulting imported values to update the
    //        local
    //        // positions

    //        // we also send the comm data that are needed to change wall
    //        positions
    //        // comm holds: (local index, id of ngb, wallpos)

    //        // (1) fill buffers
    //        vector< vector<int> > request_ids(MPIGlobal::size);
    //        vector< vector<unsigned int> > requests(MPIGlobal::size);
    //        for(AdaptiveCellList<AdaptiveVorCell2d>::orphaniterator it =
    //            _newlist.orphanbegin(); it != _newlist.orphanend(); ++it){
    //            request_ids[it->get_rank()].push_back(it->get_original());
    //            requests[it->get_rank()].push_back(it.index());
    //        }

    //        // (2) communicate
    //        // as always, we need a communication strategy
    //        // we will first send a message containing the buffer size and
    //        buffer
    //        // to the respective processes
    //        // then, we will receive any buffer size message from any source
    //        and
    //        // initialize a non blocking receive for that buffer
    //        // if a buffer is received, we do work for that buffer and send
    //        the
    //        // result back
    //        // we first also issue a non-blocking receive for the answers
    //        std::vector<MPI_Request> reqs(MPIGlobal::size*2,
    //        MPI_REQUEST_NULL);
    //        // we have to store the positions we send in a buffer that is not
    //        // deleted when an iteration ends, since we use a non-blocking
    //        send
    //        // when receiving, we can use a buffer with iteration scoop, since
    //        // the data is only used inside the iteration
    //        vector< vector<double> > sendpositions(MPIGlobal::size);
    //        for(int i = 0; i < MPIGlobal::size; i++){
    //            if(i != MPIGlobal::rank){
    //                unsigned int vsize = request_ids[i].size();
    //                MyMPI_Isend(&request_ids[i][0], vsize, MPI_INT, i, 0,
    //                            &reqs[2*i]);
    ////                cerr << MPIGlobal::rank << ": sent " << vsize << " ids
    /// to "
    ////                << i << endl;
    //            }
    //        }

    //        for(int i = 0; i < 2*(MPIGlobal::size-1); i++){
    //            MPI_Status status;
    //            MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
    //            int tag = status.MPI_TAG;
    //            int source = status.MPI_SOURCE;
    //            int nelements;
    //            if(tag == 0){
    //                MyMPI_Get_count(&status, MPI_INT, &nelements);
    //                vector<int> ids(nelements);
    //                MyMPI_Recv(&ids[0], nelements, MPI_INT, source, 0,
    //                &status);
    ////                cerr << MPIGlobal::rank << ": received " << nelements
    ////                << " ids from " << source << endl;
    //                sendpositions[source].resize(2*nelements);
    //                for(unsigned int j = 0; j < ids.size(); j++){
    //                    double pos[2];
    //                    _newlist.get_normal(ids[j])->get_position(pos);
    //                    sendpositions[source][2*j] = pos[0];
    //                    sendpositions[source][2*j+1] = pos[1];
    //                }
    //                MyMPI_Isend(&sendpositions[source][0],
    //                          sendpositions[source].size(), MPI_DOUBLE,
    //                          source, 1,
    //                          &reqs[2*source+1]);
    ////                cerr << MPIGlobal::rank << ": sent "
    ////                << sendpositions[source].size() << " positions to " <<
    /// source
    ////                << endl;
    //            }
    //            if(tag == 1){
    //                MyMPI_Get_count(&status, MPI_DOUBLE, &nelements);
    //                vector<double> recvpositions(nelements);
    //                MyMPI_Recv(&recvpositions[0], nelements, MPI_DOUBLE,
    //                source, 1,
    //                        &status);
    ////                cerr << MPIGlobal::rank << ": received " << nelements
    ////                << " positions from " << source << endl;
    //                for(unsigned int j = 0; j < recvpositions.size()/2; j++){
    //                    double pos[2];
    //                    pos[0] = recvpositions[2*j];
    //                    pos[1] = recvpositions[2*j+1];
    //                    _newlist.get_orphan(requests[source][j])->move(pos);
    //                }
    //            }
    //        }
    //        vector<MPI_Status> status((MPIGlobal::size-1)*2);
    //        MyMPI_Waitall((MPIGlobal::size-1)*2, &reqs[0], &status[0]);
    //        MyMPI_Barrier();
    //    }

    //    if(!_periodic){
    //        for(AdaptiveCellList<AdaptiveVorCell2d>::ghostiterator it =
    //            _newlist.ghostbegin(); it != _newlist.ghostend(); ++it){
    //            double pos[2];
    //            _newlist.get_normal(it->get_original())->get_position(pos);
    //            it->move(pos);
    //        }
    //    } else {
    //        for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
    //            _newlist.normalbegin(); it != _newlist.normalend(); ++it){
    //            it->apply_periodic_boundaries(it.index(), _newlist, _cuboid);
    //        }
    //        for(AdaptiveCellList<AdaptiveVorCell2d>::orphaniterator it =
    //            _newlist.orphanbegin(); it != _newlist.orphanend(); ++it){
    //            if(it->get_original() == 4){
    //                vector<int> ngbs = it->get_ngbs();
    //                cerr << "ngbs:" << endl;
    //                for(unsigned int i = 0; i < ngbs.size(); i++){
    //                    cerr << i << ": " << ngbs[i] << endl;
    //                }
    //            }
    //            it->apply_periodic_boundaries(it.index(), _newlist, _cuboid);
    //        }
    //    }
    //    // since the periodic boundaries can change the particle keys, we
    //    cannot
    //    // update the indices list before this point
    //    _indices.clear();
    //    _indices.resize(particles.gassize());
    //    for(unsigned int i = 0; i < particles.gassize(); i++){
    //        double pos[2];
    //        pos[0] = particles.gas(i)->x();
    //        pos[1] = particles.gas(i)->y();
    //        unsigned long key = AdaptiveMeshUtils::get_key(pos, _cuboid);
    //        _indices[particles.gas(i)->get_local_id()] =
    //                new ArgHilbertObject(key,
    //                particles.gas(i)->get_local_id(),
    //                                     MPIGlobal::rank);
    //    }
    for(unsigned int i = 0; i < _cells.size(); i++) {
        double pos[2];
        pos[0] = particles.gas(i)->x();
        pos[1] = particles.gas(i)->y();
        unsigned long key = AdaptiveMeshUtils::get_key(pos, _cuboid);
        _cells[i]->set_key(key);
    }
}

/**
 * @brief Check the validity of the Voronoi tesselation by explicitly checking
 * the geometrical in circle criterion for every cell
 *
 * @warning Should only be used for debugging in small simulations!
 *
 * @param stream std::ostream to write error information to
 */
void AdaptiveVorTess2d::check_mesh(ostream& stream) {
#ifdef OUD
    for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it =
                _newlist.normalbegin();
        it != _newlist.normalend(); ++it) {
        double pos[2];
        it->get_position(pos);
        vector<int> ngbs = it->get_ngbs();
        vector<int> ngbwalls = it->get_ngbwalls();
        for(unsigned int j = 0; j < ngbs.size(); j++) {
            double ngbpos[2], nextngb[2];
            _newlist[ngbs[j]]->get_position(ngbpos);
            _newlist[ngbs[(j + 1) % ngbs.size()]]->get_position(nextngb);
            if(ngbwalls[j]) {
                AdaptiveMeshUtils::get_periodic_position(ngbpos, ngbwalls[j],
                                                         _cuboid);
            }
            if(ngbwalls[(j + 1) % ngbs.size()]) {
                AdaptiveMeshUtils::get_periodic_position(
                        nextngb, ngbwalls[(j + 1) % ngbs.size()], _cuboid);
            }
            for(AdaptiveCellList<AdaptiveVorCell2d>::normaliterator it2 =
                        _newlist.normalbegin();
                it2 != _newlist.normalend(); ++it2) {
                if(it2.index() != it.index() && it2.index() != ngbs[j] &&
                   it2.index() != ngbs[(j + 1) % ngbs.size()]) {
                    double cell4[2];
                    it2->get_position(cell4);
                    double test = predicates::incircle_old(pos, ngbpos, nextngb,
                                                           cell4);
                    if(test > 0.) {
                        cerr << test << endl;
                        //                        cerr << pos[0] << "\t" <<
                        // pos[1] << endl;
                        //                        cerr << ngbpos[0] << "\t" <<
                        // ngbpos[1] << endl;
                        //                        cerr << nextngb[0] << "\t" <<
                        // nextngb[1] << endl;
                        //                        cerr << cell4[0] << "\t" <<
                        // cell4[1] << endl;
                        cerr << it.index() << "\t" << ngbs[j] << "\t"
                             << ngbs[(j + 1) % ngbs.size()] << "\t"
                             << it2.index() << endl;
                        cerr << "ERROR!" << endl;
                        //                        int plist[4] = {i, ngbs[j],
                        // ngbs[(j+1)%ngbs.size()]
                        //                        , k};
                        //                        for(unsigned int n = 0; n < 4;
                        // n++){
                        //                            stringstream filename;
                        //                            filename << "cell" <<
                        // plist[n] << ".dat";
                        //                            ofstream
                        // ofile(filename.str().c_str());
                        //                            _cells[plist[n]]->print(ofile,
                        // plist[n]);
                        //                        }
                        //                        exit(1);
                    }
                }
            }
        }
    }
#else
    for(unsigned int i = 0; i < _cells.size(); i++) {
        if(i % 100 == 0) {
            cout << "Cell " << i << "/" << _cells.size() << endl;
        }
        double pos[2];
        _cells[i]->get_position(pos);
        vector<unsigned long> ngbs = _cells[i]->get_newngbs();
        for(unsigned int j = 0; j < ngbs.size(); j++) {
            double ngbpos[2], nextngb[2];
            Vec newpos = _newerlist.get_position(ngbs[j]);
            ngbpos[0] = newpos.x();
            ngbpos[1] = newpos.y();
            newpos = _newerlist.get_position(ngbs[(j + 1) % ngbs.size()]);
            nextngb[0] = newpos.x();
            nextngb[1] = newpos.y();
            for(unsigned int k = 0; k < _cells.size(); k++) {
                if(i != k && _cells[k]->get_particleID() != ngbs[j] &&
                   _cells[k]->get_particleID() != ngbs[(j + 1) % ngbs.size()]) {
                    double cell4[2];
                    _cells[k]->get_position(cell4);
                    double test = predicates::incircle_old(pos, ngbpos, nextngb,
                                                           cell4);
                    if(test > 0.) {
                        cerr << test << endl;
                        //                        cerr << pos[0] << "\t" <<
                        //                        pos[1] << endl;
                        //                        cerr << ngbpos[0] << "\t" <<
                        //                        ngbpos[1] << endl;
                        //                        cerr << nextngb[0] << "\t" <<
                        //                        nextngb[1] << endl;
                        //                        cerr << cell4[0] << "\t" <<
                        //                        cell4[1] << endl;
                        cerr << i << "\t" << ngbs[j] << "\t"
                             << ngbs[(j + 1) % ngbs.size()] << "\t" << k
                             << endl;
                        cerr << "ERROR!" << endl;
                        my_exit();
                        //                        int plist[4] = {i, ngbs[j],
                        //                        ngbs[(j+1)%ngbs.size()]
                        //                        , k};
                        //                        for(unsigned int n = 0; n < 4;
                        //                        n++){
                        //                            stringstream filename;
                        //                            filename << "cell" <<
                        //                            plist[n] << ".dat";
                        //                            ofstream
                        //                            ofile(filename.str().c_str());
                        //                            _cells[plist[n]]->print(ofile,
                        //                            plist[n]);
                        //                        }
                        //                        exit(1);
                    }
                }
            }
        }
    }
#endif
}

/**
 * @brief Get the indices of the neighbours of the cell with the given index
 *
 * The values are sorted so that the tesselation structure can be directly
 * compared with the structure returned by e.g. the old Voronoi construction
 * algorithm.
 *
 * @param index unsigned integer index of a cell in the tesselation
 * @return std::vector with unsigned long indices of neighbouring cells
 */
vector<unsigned long> AdaptiveVorTess2d::get_ngb_ids(unsigned int index) {
    vector<unsigned long> ngb_ids;
    AdaptiveVorCell2d* cell = _newlist.get_normal(index);
    vector<int> ngbs = cell->get_ngbs();
    for(unsigned int i = 0; i < ngbs.size(); i++) {
        if(_newlist.is_normal(ngbs[i])) {
            ngb_ids.push_back(
                    _newlist.get_normal(ngbs[i])->get_particle()->id());
        } else {
            if(_newlist.is_ghost(ngbs[i])) {
                AdaptiveVorCell2d* original = _newlist.get_normal(
                        _newlist.get_ghost(ngbs[i])->get_original());
                ngb_ids.push_back(original->get_particle()->id());
            } else {
                ngb_ids.push_back(_newlist.get_orphan(ngbs[i])->get_id());
            }
        }
    }
    sort(ngb_ids.begin(), ngb_ids.end());
    return ngb_ids;
}

/**
 * @brief Dump the tesselation to the given RestartFile
 *
 * @param rfile RestartFile to write to
 * @param particles List of GasParticle instances corresponding to the cells of
 * the tesselation, used to sort the output
 */
void AdaptiveVorTess2d::dump(RestartFile& rfile,
                             vector<GasParticle*>& particles) {
    _newlist.dump(rfile, particles);
    _newerlist.dump(rfile);
    _cuboid.dump(rfile);
    _box12.dump(rfile);
    rfile.write(_periodic);
}

/**
 * @brief Restart constructor. Initialize the tesselation based on the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 * @param particles List of GasParticle instances corresponding to the cells of
 * the tesselation, used to sort the cells
 */
AdaptiveVorTess2d::AdaptiveVorTess2d(RestartFile& rfile,
                                     std::vector<GasParticle*>& particles)
        : _newlist(rfile, particles), _newerlist(rfile), _cuboid(rfile),
          _box12(rfile) {
    rfile.read(_periodic);
}
