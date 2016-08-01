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
 * @file DelTess_mpi.cpp
 *
 * @brief Delaunay tesselation: part of implementation that depends on MPI
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
#include "NgbSearch.hpp"
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
            W.pack_data(&MPIGlobal::sendbuffer[buffers[i]], bufsize, &send_pos);
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
                        W.pack_data(&MPIGlobal::sendbuffer[buffers[j]], bufsize,
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
                StateVector W(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                              nelements, &recv_pos);
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
                        W.pack_data(&MPIGlobal::sendbuffer[buffers[j]], bufsize,
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
            for(unsigned int j = 0; j < ndim_; j++) {
                gradients[j].pack_data(&MPIGlobal::sendbuffer[buffers[i]],
                                       bufsize, &send_pos);
            }
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
                        for(unsigned int k = 0; k < ndim_; k++) {
                            gradients[k].pack_data(
                                    &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                    &send_pos);
                        }
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
                for(unsigned int k = 0; k < ndim_; k++) {
                    gradients[k] = StateVector(
                            &MPIGlobal::recvbuffer[buffers[freebuffer]],
                            nelements, &recv_pos);
                }
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
                        for(unsigned int k = 0; k < ndim_; k++) {
                            gradients[k].pack_data(
                                    &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                    &send_pos);
                        }
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
