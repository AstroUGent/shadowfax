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
 * @file Tree.cpp
 *
 * @brief Octree used for force calculations, neighbour searches...:
 * implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Tree.hpp"
#include "../io/PhysicalConstant.hpp"
#include "../io/Unit.hpp"
#include "../io/UnitSet.hpp"
#include "GasParticle.hpp"
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include "Particle.hpp"
#include "StateVector.hpp"
#include "TreeWalker.hpp"
#include "VorGen.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

// Tree

/**
 * @brief Constructor
 *
 * @param box Cuboid specifying the dimensions of the box in which the tree is
 * embedded
 * @param periodic Flag indicating if the box is periodic or reflective
 * @param do_ewald Flag indicating if the Ewald tables should be loaded, which
 * are used to calculated periodic force corrections for gravity
 * @param alpha Ewald alpha parameter
 * @param size Size of the Ewald table
 */
Tree::Tree(Cuboid box, bool periodic, bool do_ewald, double alpha,
           unsigned int size)
        : _box(box) {
    _periodic = periodic;

    if(do_ewald) {
        _ewald_table = new EwaldTable(alpha, size);
    } else {
        _ewald_table = NULL;
    }

    Vec center;
    double side = 0.;
    for(unsigned int i = ndim_; i--;) {
        center[i] = box.get_anchor()[i] + 0.5 * box.get_sides()[i];
        side = std::max(side, box.get_sides()[i]);
    }
    _side = side;
    Box treebox(center, side);
    _root = new TreeNode(treebox);
    if(_ewald_table && side) {
        _ewald_table->set_boxsize(side);
    }
}

/**
 * @brief Destructor
 *
 * Clean up the root of the tree (which recursively cleans up all children), the
 * exportnodes and the ewald table (if loaded).
 */
Tree::~Tree() {
    delete _root;
    // _root contains all implementations of =Node
    // ExportNode is no implementation of Node and has to be deleted
    // seperately
    for(unsigned int i = _exportnodes.size(); i--;) {
        delete _exportnodes[i];
    }
    if(_ewald_table) {
        delete _ewald_table;
    }
}

/**
 * @brief Set the periodicity flag of the tree
 *
 * @param periodic True if the tree is embedded in a periodic box, false if it
 * is embedded in a reflective box
 */
void Tree::set_periodic(bool periodic) { _periodic = periodic; }

/**
 * @brief Reset the tree
 *
 * The internal state of the tree is reset as if the tree was deleted and
 * constructed again, but without having to reload Ewald tables or without
 * having to specify the kind of boundary conditions again.
 *
 * @param box Cuboid specifying the dimensions of the box in which the tree is
 * embedded
 */
void Tree::reset(Cuboid box) {
    delete _root;
    for(unsigned int i = _exportnodes.size(); i--;) {
        delete _exportnodes[i];
    }
    _exportnodes.clear();
    _pseudonodes.clear();
    _box = box;
    Vec center;
    double side = 0.;
    for(unsigned int i = ndim_; i--;) {
        center[i] = box.get_anchor()[i] + 0.5 * box.get_sides()[i];
        side = std::max(side, box.get_sides()[i]);
    }
    // if for some reason the boxsize changed: rescale the ewald table
    if(_ewald_table && side) {
        if(_side) {
            if(_side != side) {
                double ratio = side / _side;
                _ewald_table->set_boxsize(ratio);
            }
        } else {
            _ewald_table->set_boxsize(side);
        }
    }
    _side = side;
    Box treebox(center, side);
    _root = new TreeNode(treebox);
}

/**
 * @brief Add the given Particle to the tree
 *
 * @param p Particle to add
 */
void Tree::add_particle(Particle* p) { ((TreeNode*)_root)->add_particle(p); }

/**
 * @brief Signal that no more particles will be added to the tree and
 * restructure the tree so that walking it is more efficient
 *
 * @return Number of nodes in the tree
 */
unsigned int Tree::finalize() {
    return ((TreeNode*)_root)->finalize(NULL, false);
}

/**
 * @brief Exchange node information between MPI-processes
 *
 * Every MPI-process holds a complete copy of the tree, with nodes that are not
 * present locally represented by a PseudoNode. In this method, we update the
 * information stored in the local pseudonodes with accurate information from
 * the MPI-process that contains the actual TreeNode represented by the
 * pseudonode.
 */
void Tree::exchange_pseudonodes() {
    // sort _exportnodes
    sort(_exportnodes.begin(), _exportnodes.end(), HB::sortfunc);
    // sort _pseudonodes
    sort(_pseudonodes.begin(), _pseudonodes.end(), HB::sortfunc);
    // fill export vector with FancyNodeInfo
    vector<NodeInfo> exportvector(_exportnodes.size());
    // order is quite the key of this whole thing, so let's say it's important
    for(unsigned int i = 0; i < _exportnodes.size(); i++) {
        Node* node = _exportnodes[i]->get_node();
        Vec com;
        if(node) {
            for(unsigned int j = ndim_; j--;) {
                com[j] = node->get_center_of_mass(j);
            }
            exportvector[i] = NodeInfo(_exportnodes[i]->get_key(),
                                       node->get_cmax(), node->get_vmax(),
                                       node->get_mass(), node->get_hmax(), com);
        } else {
            exportvector[i] =
                    NodeInfo(_exportnodes[i]->get_key(), 0., 0., 0., 0., com);
        }
    }

    // communicate
    vector<vector<NodeInfo> > importvectors(MPIGlobal::size);
    vector<MPI_Request> reqs(MPIGlobal::size - 1, MPI_REQUEST_NULL);
    vector<unsigned int> buffers(MPIGlobal::size - 1);
    unsigned int bufsize = MPIGlobal::sendsize / buffers.size();
    for(unsigned int i = 0; i < buffers.size(); i++) {
        buffers[i] = i * bufsize;
    }

    unsigned int maxsize = bufsize / sizeof(NodeInfo);
    // make sure we have some space left at the end of the buffer to signal
    // more to come
    if(!(bufsize % sizeof(NodeInfo))) {
        maxsize--;
    }

    // how does it work for small buffers:
    // we start by sending as many elements as possible. At the end of the
    // buffer, we add an extra flag that signals if all elements were send or
    // not. When receiving the message, we read this flag and this way know if
    // we have to do an extra loop to wait for an extra message.
    // At the start of every receive loop, we check if a local send request has
    // finished. If so, we can reuse the associated buffer to again send as many
    // elements as possible, with a flag. We then continue to wait for requests.
    // If all requests have been received, we check if there are still messages
    // that need to be sent and then wait for the other processes to finish.

    vector<int> allflag(MPIGlobal::size - 1);
    vector<unsigned int> numsend(MPIGlobal::size - 1);
    unsigned int sendsize = exportvector.size();
    // for some reason, reusing the same buffer for every send does not work
    // we therefore use a separate send buffer for every process
    int allflagsum = 0;
    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        int send_pos = 0;
        for(unsigned int si = 0; si < std::min(maxsize, sendsize); si++) {
            exportvector[si].pack_data(&MPIGlobal::sendbuffer[buffers[i]],
                                       bufsize, &send_pos);
        }
        allflag[i] = (sendsize <= maxsize);
        // if allflagsum equals MPIGlobal::size-1, all data has been sent
        allflagsum += allflag[i];
        numsend[i] = 1;
        // add continuation signal
        MyMPI_Pack(&allflag[i], 1, MPI_INT, &MPIGlobal::sendbuffer[buffers[i]],
                   bufsize, &send_pos);

        MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[i]], send_pos, MPI_PACKED,
                    (MPIGlobal::rank + 1 + i) % MPIGlobal::size, 0, &reqs[i]);
    }

    for(int i = 0; i < MPIGlobal::size - 1; i++) {
        MPI_Status status;
        for(int j = 0; j < MPIGlobal::size - 1; j++) {
            if(!allflag[j]) {
                int flag;
                MyMPI_Test(&reqs[j], &flag, &status);
                if(flag) {
                    int send_pos = 0;
                    for(unsigned int si = numsend[j] * maxsize;
                        si < std::min((numsend[j] + 1) * maxsize, sendsize);
                        si++) {
                        exportvector[si].pack_data(
                                &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                &send_pos);
                    }
                    allflag[j] = (sendsize <= (numsend[j] + 1) * maxsize);
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
            MyMPI_Recv(&MPIGlobal::recvbuffer[freebuffer], nelements,
                       MPI_PACKED, index, tag, &status);
            recv_pos = 0;
            while(recv_pos < nelements - 4) {
                importvectors[index].push_back(
                        NodeInfo(&MPIGlobal::recvbuffer[freebuffer], nelements,
                                 &recv_pos));
            }
            int flag;
            MyMPI_Unpack(&MPIGlobal::recvbuffer[freebuffer], nelements,
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
                        si < std::min((numsend[j] + 1) * maxsize, sendsize);
                        si++) {
                        exportvector[si].pack_data(
                                &MPIGlobal::sendbuffer[buffers[j]], bufsize,
                                &send_pos);
                    }
                    allflag[j] = (sendsize <= (numsend[j] + 1) * maxsize);
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

    vector<MPI_Status> status(MPIGlobal::size - 1);
    MyMPI_Waitall(MPIGlobal::size - 1, &reqs[0], &status[0]);
    MyMPI_Barrier();

    // fill local _pseudonodes
    vector<NodeInfo> importvector;
    // again, order is important!
    for(int i = 0; i < MPIGlobal::size; i++) {
        if(i != MPIGlobal::rank) {
            for(unsigned int j = 0; j < importvectors[i].size(); j++) {
                importvector.push_back(importvectors[i][j]);
            }
        }
    }

    // we made sure index i in _pseudonodes corresponds to index i in
    // importvector, so now order is not important anymore
    for(unsigned int i = importvector.size(); i--;) {
        PseudoNode* node = _pseudonodes[i];
        node->set_center_of_mass(importvector[i].get_center_of_mass());
        node->set_cmax(importvector[i].get_cmax());
        node->set_vmax(importvector[i].get_vmax());
        node->set_mass(importvector[i].get_mass());
        node->set_hmax(importvector[i].get_hmax());
    }
    // update nodes that are affected by these new data
    ((TreeNode*)_root)->update_quantities();
}

/**
 * @brief Add the pseudoparticles for the given MPI-process in the given
 * Hilbert key range to the tree
 *
 * @param keylow Lower bound of the Hilbert key range
 * @param keyhigh Upper bound of the Hilbert key range
 * @param pindex Rank of the MPI-process which holds the actual treenodes
 * represented by the newly added pseudonodes
 */
void Tree::add_pseudoparticles(unsigned long keylow, unsigned long keyhigh,
                               unsigned int pindex) {
    ((TreeNode*)_root)
            ->add_pseudoparticles(keylow, keyhigh, pindex, 0, 1, _pseudonodes);
}

/**
 * @brief Set the Hilbert key range of the tree on this MPI-process
 *
 * @param keylow Lower bound of the Hilbert key range
 * @param keyhigh Upper bound of the Hilbert key range
 */
void Tree::set_keyrange(unsigned long keylow, unsigned long keyhigh) {
    _keyrange[0] = keylow;
    _keyrange[1] = keyhigh;
    ((TreeNode*)_root)->get_exportnodes(keylow, keyhigh, 0, 1, _exportnodes);
}

/**
 * @brief Print the tree to the given stream
 * @param out std::ostream to write to
 */
void Tree::print(ostream& out) { _root->print(out); }

/**
 * @brief Print global properties of the tree to the given stream
 *
 * @param out std::ostream to write to
 */
void Tree::print_top_quantities(ostream& out) {
    out << "Total mass: " << _root->get_mass() << endl;
    out << "Mass center: (" << _root->get_center_of_mass(0) << ","
        << _root->get_center_of_mass(1) << ")" << endl;
    out << "Cmax: " << _root->get_cmax() << endl;
    out << "Vmax: " << _root->get_vmax() << endl;
}

/**
 * @brief Set velocities and related quantities for the nodes of the tree
 *
 * If Particle information changes after the tree construction, information in
 * the nodes of the tree can get outdated. We update this information by
 * recursively updating the nodes from bottom to top.
 */
void Tree::set_velocities() { ((TreeNode*)_root)->update_quantities(true); }

/**
 * @brief Update the mesh masses
 *
 * @deprecated This method is not used
 */
void Tree::set_mesh_masses() { ((TreeNode*)_root)->update_mesh_masses(); }

/**
 * @brief Get the particles inside the sphere with given origin and radius
 *
 * @param coords Origin of the search sphere
 * @param radius Radius of the search sphere
 * @param ngblist Reference to the list to fill
 * @param exportlist List used for MPI communication
 */
void Tree::get_neighbours(Vec& coords, double radius,
                          vector<GasParticle*>& ngblist,
                          vector<bool>& exportlist) {
    Node* stack = ((TreeNode*)_root)->get_child();
    double r2 = radius * radius;
    while(stack) {
        if(stack->is_leaf()) {
            if(stack->is_pseudo()) {
                if(((PseudoNode*)stack)->get_distance(coords) <= r2) {
                    exportlist[((PseudoNode*)stack)->get_source()] = true;
                }
                stack = ((PseudoNode*)stack)->get_sibling();
            } else {
                if(((Leaf*)stack)->get_particle()->type() == PARTTYPE_GAS) {
                    if((coords - ((Leaf*)stack)->get_particle()->get_position())
                               .norm2() <= r2) {
                        ngblist.push_back(
                                (GasParticle*)((Leaf*)stack)->get_particle());
                    }
                }
                stack = ((Leaf*)stack)->get_sibling();
            }
        } else {
            if(((TreeNode*)stack)->get_distance(coords) <= r2) {
                stack = ((TreeNode*)stack)->get_child();
            } else {
                stack = ((TreeNode*)stack)->get_sibling();
            }
        }
    }
}

/**
 * @brief Get the neighbours inside the sphere with given origin and radius that
 * lie inside on of the periodic or reflective copies of the box in which the
 * tree is embedded, but not in the box itself
 *
 * @param coords Origin of the search sphere
 * @param radius Radius of the search sphere
 * @param ngblist Reference to the list to fill
 * @param exportlist List used for MPI communication
 * @param index Rank of the MPI process on which the search is performed,
 * possibly always 0
 */
void Tree::get_neighbours_outside(Vec& coords, double radius,
                                  vector<VorGen*>& ngblist,
                                  vector<bool>& exportlist,
                                  unsigned int index) {
#if ndim_ == 3
    if(_periodic) {
        for(unsigned int i = 3; i--;) {
            for(unsigned int j = 3; j--;) {
                for(unsigned int k = 3; k--;) {
                    if(!(i == 1 && j == 1 && k == 1)) {
                        double xfac = 1. - i;
                        double yfac = 1. - j;
                        double zfac = 1. - k;
                        Vec coords_copy(coords[0] - xfac * _box.get_sides()[0],
                                        coords[1] - yfac * _box.get_sides()[1],
                                        coords[2] - zfac * _box.get_sides()[2]);
                        vector<GasParticle*> copy_ngbs;
                        get_neighbours(coords_copy, radius, copy_ngbs,
                                       exportlist);
                        for(unsigned int l = copy_ngbs.size(); l--;) {
                            if(!copy_ngbs[l]->is_copied(i * 9 + j * 3 + k + 1,
                                                        index)) {
                                ngblist.push_back(new VorGen(
                                        copy_ngbs[l]->x() +
                                                xfac * _box.get_sides()[0],
                                        copy_ngbs[l]->y() +
                                                yfac * _box.get_sides()[1],
                                        copy_ngbs[l]->z() +
                                                zfac * _box.get_sides()[2]));
                                copy_ngbs[l]->make_copy(i * 9 + j * 3 + k + 1,
                                                        index);
                                ngblist.back()->set_particle(copy_ngbs[l]);
                            }
                        }
                    }
                }
            }
        }
    } else {
        for(unsigned int i = 3; i--;) {
            for(unsigned int j = 3; j--;) {
                for(unsigned int k = 3; k--;) {
                    if(!(i == 1 && j == 1 && k == 1)) {
                        Vec coords_copy(coords[0], coords[1], coords[2]);
                        if(i) {
                            if(i > 1) {
                                coords_copy[0] = 2. * (_box.get_anchor()[0] +
                                                       _box.get_sides()[0]) -
                                                 coords_copy[0];
                            }
                        } else {
                            coords_copy[0] =
                                    2. * _box.get_anchor()[0] - coords_copy[0];
                        }
                        if(j) {
                            if(j > 1) {
                                coords_copy[1] = 2. * (_box.get_anchor()[1] +
                                                       _box.get_sides()[1]) -
                                                 coords_copy[1];
                            }
                        } else {
                            coords_copy[1] =
                                    2. * _box.get_anchor()[1] - coords_copy[1];
                        }
                        if(k) {
                            if(k > 1) {
                                coords_copy[2] = 2. * (_box.get_anchor()[2] +
                                                       _box.get_sides()[2]) -
                                                 coords_copy[2];
                            }
                        } else {
                            coords_copy[2] =
                                    2. * _box.get_anchor()[2] - coords_copy[2];
                        }
                        vector<GasParticle*> copy_ngbs;
                        get_neighbours(coords_copy, radius, copy_ngbs,
                                       exportlist);
                        for(unsigned int l = copy_ngbs.size(); l--;) {
                            if(!copy_ngbs[l]->is_copied(i * 9 + j * 3 + k +
                                                        1)) {
                                double x[3] = {copy_ngbs[l]->x(),
                                               copy_ngbs[l]->y(),
                                               copy_ngbs[l]->z()};
                                if(i) {
                                    if(i > 1) {
                                        x[0] = 2. * (_box.get_anchor()[0] +
                                                     _box.get_sides()[0]) -
                                               x[0];
                                    }
                                } else {
                                    x[0] = 2. * _box.get_anchor()[0] - x[0];
                                }
                                if(j) {
                                    if(j > 1) {
                                        x[1] = 2. * (_box.get_anchor()[1] +
                                                     _box.get_sides()[1]) -
                                               x[1];
                                    }
                                } else {
                                    x[1] = 2. * _box.get_anchor()[1] - x[1];
                                }
                                if(k) {
                                    if(k > 1) {
                                        x[2] = 2. * (_box.get_anchor()[2] +
                                                     _box.get_sides()[2]) -
                                               x[2];
                                    }
                                } else {
                                    x[2] = 2. * _box.get_anchor()[2] - x[2];
                                }
                                ngblist.push_back(new VorGen(x[0], x[1], x[2]));
                                copy_ngbs[l]->make_copy(i * 9 + j * 3 + k + 1);
                                ngblist.back()->set_particle(copy_ngbs[l]);
                            }
                        }
                    }
                }
            }
        }
    }
#else
    if(_periodic) {
        for(unsigned int i = 3; i--;) {
            for(unsigned int j = 3; j--;) {
                if(!(i == 1 && j == 1)) {
                    double xfac = 1. - i;
                    double yfac = 1. - j;
                    Vec coords_copy(coords[0] - xfac * _box.get_sides()[0],
                                    coords[1] - yfac * _box.get_sides()[1]);
                    vector<GasParticle*> copy_ngbs;
                    get_neighbours(coords_copy, radius, copy_ngbs, exportlist);
                    for(unsigned int k = copy_ngbs.size(); k--;) {
                        if(!copy_ngbs[k]->is_copied(i * 3 + j + 1, index)) {
                            ngblist.push_back(new VorGen(
                                    copy_ngbs[k]->x() +
                                            xfac * _box.get_sides()[0],
                                    copy_ngbs[k]->y() +
                                            yfac * _box.get_sides()[1]));
                            copy_ngbs[k]->make_copy(i * 3 + j + 1, index);
                            ngblist.back()->set_particle(copy_ngbs[k]);
                        }
                    }
                }
            }
        }
    } else {
        for(unsigned int i = 3; i--;) {
            for(unsigned int j = 3; j--;) {
                if(!(i == 1 && j == 1)) {
                    Vec coords_copy(coords[0], coords[1]);
                    if(i) {
                        if(i > 1) {
                            coords_copy[0] = 2. * (_box.get_anchor()[0] +
                                                   _box.get_sides()[0]) -
                                             coords_copy[0];
                        }
                    } else {
                        coords_copy[0] =
                                2. * _box.get_anchor()[0] - coords_copy[0];
                    }
                    if(j) {
                        if(j > 1) {
                            coords_copy[1] = 2. * (_box.get_anchor()[1] +
                                                   _box.get_sides()[1]) -
                                             coords_copy[1];
                        }
                    } else {
                        coords_copy[1] =
                                2. * _box.get_anchor()[1] - coords_copy[1];
                    }
                    vector<GasParticle*> copy_ngbs;
                    get_neighbours(coords_copy, radius, copy_ngbs, exportlist);
                    for(unsigned int k = copy_ngbs.size(); k--;) {
                        if(!copy_ngbs[k]->is_copied(i * 3 + j + 1)) {
                            double x[2] = {copy_ngbs[k]->x(),
                                           copy_ngbs[k]->y()};
                            if(i) {
                                if(i > 1) {
                                    x[0] = 2. * (_box.get_anchor()[0] +
                                                 _box.get_sides()[0]) -
                                           x[0];
                                }
                            } else {
                                x[0] = 2. * _box.get_anchor()[0] - x[0];
                            }
                            if(j) {
                                if(j > 1) {
                                    x[1] = 2. * (_box.get_anchor()[1] +
                                                 _box.get_sides()[1]) -
                                           x[1];
                                }
                            } else {
                                x[1] = 2. * _box.get_anchor()[1] - x[1];
                            }
                            ngblist.push_back(new VorGen(x[0], x[1]));
                            copy_ngbs[k]->make_copy(i * 3 + j + 1);
                            ngblist.back()->set_particle(copy_ngbs[k]);
                        }
                    }
                }
            }
        }
    }
#endif
}

/**
 * @brief Get the closest particle to the given coordinates, limiting the
 * initial search to a sphere with given radius around the coordinates
 *
 * If after the initial search no Particle was found, the radius of the sphere
 * is doubled and the search is repeated. This is repeated until a Particle is
 * found.
 *
 * @param coords Coordinates for which the search is done
 * @param r Maximum radius of the search sphere around the coordinates
 * @return Closest Particle to the given coordinates
 */
Particle* Tree::get_closest(Vec& coords, double r) {
    Particle* closest = NULL;
    double r2 = r * r;
    while(!closest) {
        Node* current = ((TreeNode*)_root)->get_child();
        while(current) {
            if(current->is_leaf()) {
                double check =
                        (((Leaf*)current)->get_particle()->get_position() -
                         coords).norm2();
                if(check < r2) {
                    closest = ((Leaf*)current)->get_particle();
                    r2 = check;
                }
                current = ((Leaf*)current)->get_sibling();
            } else {
                if(((TreeNode*)current)->get_distance(coords) < r2) {
                    current = ((TreeNode*)current)->get_child();
                } else {
                    current = ((TreeNode*)current)->get_sibling();
                }
            }
        }
        r2 *= 4.;
    }
    return closest;
}

/**
 * @brief Get the closest particle to the given coordinates, taking into account
 * the periodic boundaries
 *
 * @deprecated This function just links to Tree::get_closest()
 *
 * @param coords Coordinates for which the search is done
 * @param r Maximum radius around the coordinates used for the first search
 * @return Closest Particle to the given coordinates
 */
Particle* Tree::get_periodic_closest(Vec& coords, double r) {
    Particle* closest = get_closest(coords, r);
    return closest;
}

// TreeNode

/**
 * @brief Constructor
 *
 * @param box Box specifying the dimensions of the embedding box of this node
 */
TreeNode::TreeNode(Box box) : Node(false), _box(box) {
    for(unsigned int i = numnode_; i--;) {
        _nodes[i] = NULL;
    }
}

/**
 * @brief Destructor
 *
 * Recursively clean up all children of the node and if it exists, the sibling
 * node on the same level.
 */
TreeNode::~TreeNode() {
    delete _child;
    // every node deletes its own children, so if the nextnode is child of a
    // different parent, we don't delete it here
    if(!is_last()) {
        delete _sibling;
    }
}

/**
 * @brief Get the total mass in the node
 *
 * This version is the implementation of the corresponding method in Node. If
 * you know the Node is a TreeNode, it is faster to use TreeNode::get_mass_node,
 * which saves a function call.
 *
 * @return Total mass in the node
 */
double TreeNode::get_mass() { return _total_mass; }

/**
 * @brief Get the total mass in the node
 *
 * @return Total mass in the node
 */
double TreeNode::get_mass_node() { return _total_mass; }

/**
 * @brief Get the maximal softening length in the node
 *
 * @return Maximal softening length in the node
 */
double TreeNode::get_hmax() { return _hmax; }

/**
 * @brief Get the given coordinate of the center of mass of the node
 *
 * This version is the implementation of the corresponding method in Node. If
 * you know the Node is a TreeNode, it is faster to use
 * TreeNode::get_center_of_mass_node, which saves a function call.
 *
 * @param index Requested coordinate
 * @return Coordinate of the center of mass
 */
double TreeNode::get_center_of_mass(unsigned int index) {
    return _center_of_mass[index];
}

/**
 * @brief Get the given coordinate of the center of mass of the node
 *
 * @param index Requested coordinate
 * @return Coordinate of the center of mass
 */
double TreeNode::get_center_of_mass_node(unsigned int index) {
    return _center_of_mass[index];
}

/**
 * @brief Get the center of mass of the node
 *
 * @return The center of mass of the node
 */
Vec TreeNode::get_center_of_mass_node() {
#if ndim_ == 3
    return Vec(_center_of_mass[0], _center_of_mass[1], _center_of_mass[2]);
#else
    return Vec(_center_of_mass[0], _center_of_mass[1]);
#endif
}

/**
 * @brief Get the maximal soundspeed in the node
 *
 * @return Maximal soundspeed in the node
 */
double TreeNode::get_cmax() { return _cmax; }

/**
 * @brief Get the maximal fluid velocity in the node
 *
 * @return Maximal fluid velocity in the node
 */
double TreeNode::get_vmax() { return _vmax; }

/**
 * @brief Set the mesh mass of the node
 *
 * @deprecated Not used
 *
 * @param mesh_m Mesh mass of the node
 */
void TreeNode::set_mesh_m(double mesh_m) { _mesh_m = mesh_m; }

/**
 * @brief Get the mesh mass of the node
 *
 * @deprecated Not used
 *
 * @return Mesh mass of the node
 */
double TreeNode::get_mesh_m() { return _mesh_m; }

/**
 * @brief Get the distance between the given coordinates and the box of this
 * node
 *
 * @param coords Vec with coordinates
 * @return Distance between the node and the coordinates
 */
double TreeNode::get_distance(Vec& coords) {
    double distance = 0.;
    double x;
    for(unsigned int i = ndim_; i--;) {
        x = fabs(_box.get_anchor()[i] - coords[i]);
        if(x < 0.5 * _box.get_side()) {
            x = 0.;
        } else {
            x -= 0.5 * _box.get_side();
        }
        distance += x * x;
    }
    return distance;
}

/**
 * @brief Get the first child of the node
 *
 * @return First child of the node
 */
Node* TreeNode::get_child() { return _child; }

/**
 * @brief Get the Node on the same level that is next in the tree traversal
 * algorithm
 *
 * @return The next Node in the tree traversal algorithm
 */
Node* TreeNode::get_sibling() { return _sibling; }

/**
 * @brief Get the width of the node
 *
 * @return Width of the node
 */
double TreeNode::get_width() { return _box.get_side(); }

/**
 * @brief Get the center of the node
 *
 * @return Center of the node
 */
Vec TreeNode::get_center() { return _box.get_anchor(); }

/**
 * @brief Add the given particle on the given level to the node
 *
 * This method is the kernel of the tree construction algorithm. We check in
 * which child node of the node the particle lies by using its Hilbert key.
 * If this node is a TreeNode, we call TreeNode::add_particle on this node. If
 * it is a Leaf, we delete it and create a new TreeNode containing the Particle
 * stored in the old Leaf and the new Particle.
 *
 * @param p Particle to add
 * @param level Level of the current node in the tree
 */
void TreeNode::add_particle(Particle* p, unsigned int level) {
    unsigned long key = p->get_key();
    unsigned long nodekey = key >> (60 - ndim_ * (level + 1));
// we isolate the two or three last bits (3 = ..0011; 7 = ..000111)
#if ndim_ == 3
    unsigned int index = nodekey & 7;
#else
    unsigned int index = nodekey & 3;
#endif
    if(_nodes[index] == NULL) {
        _nodes[index] = new Leaf();
        ((Leaf*)_nodes[index])->add_particle(p);
    } else {
        if(_nodes[index]->is_leaf()) {
            Particle* p2 = ((Leaf*)_nodes[index])->get_particle();
            delete _nodes[index];
#if ndim_ == 3
            Vec anchor = _box.get_anchor();
            unsigned long coords[3] = {0};
            HB::get_coords(nodekey << (60 - ndim_ * (level + 1)), 60, coords);
            if(!(coords[0] >> (20 - (level + 1)) & 1)) {
                anchor[0] -= 0.25 * _box.get_side();
            } else {
                anchor[0] += 0.25 * _box.get_side();
            }
            if(!(coords[1] >> (20 - (level + 1)) & 1)) {
                anchor[1] -= 0.25 * _box.get_side();
            } else {
                anchor[1] += 0.25 * _box.get_side();
            }
            if(!(coords[2] >> (20 - (level + 1)) & 1)) {
                anchor[2] -= 0.25 * _box.get_side();
            } else {
                anchor[2] += 0.25 * _box.get_side();
            }
#else
            Vec anchor = _box.get_anchor();
            unsigned long coords[2] = {0};
            HB::get_coords(nodekey << (60 - ndim_ * (level + 1)), 60, coords);
            if(!(coords[0] >> (30 - (level + 1)) & 1)) {
                anchor[0] -= 0.25 * _box.get_side();
            } else {
                anchor[0] += 0.25 * _box.get_side();
            }
            if(!(coords[1] >> (30 - (level + 1)) & 1)) {
                anchor[1] -= 0.25 * _box.get_side();
            } else {
                anchor[1] += 0.25 * _box.get_side();
            }
#endif
            Box box(anchor, 0.5 * _box.get_side());
            _nodes[index] = new TreeNode(box);
            ((TreeNode*)_nodes[index])->add_particle(p2, level + 1);
            ((TreeNode*)_nodes[index])->add_particle(p, level + 1);
        } else {
            ((TreeNode*)_nodes[index])->add_particle(p, level + 1);
        }
    }
}

/**
 * @brief Restructure the treenode so that it can be traversed more efficiently
 * and recursively calculate all node quantities
 *
 * During construction, every TreeNode stores pointers to its 8 (4 in 2D) child
 * nodes. We replace these by two pointers: a pointer to the first child and a
 * pointer to the next node on the same level. If the tree is walked, two things
 * can happen:
 *  - the node is opened: we continue the tree walk with the first child
 *  - the node is treated as a whole: we continue the walk with the next node
 *    on the same level
 *
 * @param sibling Node on the same level that is next in the efficient tree walk
 * algorithm
 * @param last Flag indicating whether this node is the last child of its parent
 * (in which case it does not have a sibling)
 * @return Number of subnodes in this node
 */
unsigned int TreeNode::finalize(Node* sibling, bool last) {
    unsigned int count = 1;
    // _total_mass and _center_of_mass are in a union with _nodes, so we have to
    // calculate them in different variables
    double mass = 0.;
    double center_of_mass[ndim_] = {0.};
    double cmax = 0.;
    double vmax = 0.;
    double hmax = 0.;
    bool only_local = true;
    // call finalize on the children (order is important!)
    for(unsigned int i = 0; i < numnode_; i++) {
        if(_nodes[i]) {
            unsigned int j = i + 1;
            while(j < numnode_ && _nodes[j] == NULL) {
                j++;
            }
            Node* next = sibling;
            bool lastchild = true;
            if(j < numnode_) {
                next = _nodes[j];
                lastchild = false;
            }
            // we only want the properties of the local particles
            count += _nodes[i]->finalize(next, lastchild);
            if(!_nodes[i]->is_pseudo()) {
                mass += _nodes[i]->get_mass();
                for(unsigned int k = ndim_; k--;) {
                    center_of_mass[k] += _nodes[i]->get_mass() *
                                         _nodes[i]->get_center_of_mass(k);
                }
                cmax = std::max(cmax, _nodes[i]->get_cmax());
                vmax = std::max(vmax, _nodes[i]->get_vmax());
                hmax = std::max(hmax, _nodes[i]->get_hmax());
            }
            only_local &= _nodes[i]->is_local();
        }
    }
    // set sibling and child
    unsigned int index = 0;
    while(_nodes[index] == NULL) {
        index++;
    }
    // index will always have a value < numnode_, because otherwise this
    // TreeNode would be a Leaf...
    // except of course if the tree is completely empty...
    Node* child = NULL;
    if(index < numnode_) {
        child = _nodes[index];
    }
    // we don't need _nodes anymore
    _total_mass = mass;
    if(mass) {
        for(unsigned int k = ndim_; k--;) {
            _center_of_mass[k] = center_of_mass[k] / mass;
        }
    }
    _flag |= (last << 1);

    // if the node contains pseudonodes on a lower level, set the appropriate
    // flag
    _flag |= ((!only_local) << 3);

    _sibling = sibling;
    _child = child;
    _cmax = cmax;
    _vmax = vmax;
    _hmax = hmax;
    return count;
}

/**
 * @brief Print the node to the given stream
 *
 * @param out std::ostream to write to
 */
void TreeNode::print(ostream& out) {
    for(unsigned int i = 0; i < ndim_; i++) {
        out << _box.get_anchor()[i] << "\t";
    }
    out << _box.get_side() << "\n";
    for(unsigned int i = 0; i < numnode_; i++) {
        if(_nodes[i] != NULL) {
            _nodes[i]->print(out);
        }
    }
    out << "\n";
}

/**
  * @brief Function to check if the intervals [a1, a2] and [b1, b2] (both with
  * inclusive borders) overlap
  *
  * We have to split the problem depending on the relative sizes of the
  * intervals: if the smallest interval is entirely inside the largest, then the
  * endpoints of the largest will be outside the smallest interval. Hence we
  * always check if the smallest is inside the largest interval.
  *
  * @param a1 Lower boundary of the first interval
  * @param a2 Upper boundary of the first interval
  * @param b1 Lower boundary of the second interval
  * @param b2 Upper boundary of the second interval
  * @return True if the intervals partially or completely overlap
  */
bool TreeNode::interval_overlap(unsigned long a1, unsigned long a2,
                                unsigned long b1, unsigned long b2) {
    bool overlap;
    if((a2 - a1) >= (b2 - b1)) {
        overlap = ((a1 <= b1 && b1 <= a2) || (a1 <= b2 && b2 <= a2));
    } else {
        overlap = ((b1 <= a1 && a1 <= b2) || (b1 <= a2 && a2 <= b2));
    }
    return overlap;
}

/**
 * @brief Add pseudoparticles to the children of the node, corresponding to
 * parts of the tree that are located on other MPI-processes
 *
 * This method is very simular to TreeNode::add_particle, with the difference
 * that the particles that are added are pseudonodes.
 *
 * @param keylow Lower bound of the Hilbert key range for the MPI process
 * @param keyhigh Upper bound of the Hilbert key range for the MPI process
 * @param pindex Rank of the MPI process
 * @param level Level of this node in the tree
 * @param key Hilbert key of a position inside this node
 * @param pseudonodes Reference to the PseudoNode list of the Tree
 */
void TreeNode::add_pseudoparticles(unsigned long keylow, unsigned long keyhigh,
                                   unsigned int pindex, unsigned int level,
                                   unsigned long key,
                                   std::vector<PseudoNode*>& pseudonodes) {
    unsigned long nodekey = key << ndim_;
    unsigned long keysize = 1;
    keysize <<= (60 - ndim_ * (level + 1));
    keysize -= 1;
    for(unsigned int i = numnode_; i--;) {
        unsigned long subnodekey = nodekey + i;
        subnodekey <<= (60 - ndim_ * (level + 1));
        if(interval_overlap(keylow, keyhigh, subnodekey,
                            subnodekey + keysize)) {
            Vec anchor = _box.get_anchor();
            unsigned long coords[ndim_] = {0};
            HB::get_coords(subnodekey, 60, coords);
            for(unsigned int j = ndim_; j--;) {
                coords[j] >>= (60 / ndim_ - (level + 1));
                if(!(coords[j] & 1)) {
                    anchor[j] -= 0.25 * _box.get_side();
                } else {
                    anchor[j] += 0.25 * _box.get_side();
                }
            }
            Box box(anchor, 0.5 * _box.get_side());
            if(_nodes[i] == NULL) {
                if(keylow <= subnodekey && subnodekey + keysize <= keyhigh) {
                    _nodes[i] = new PseudoNode(box, pindex, subnodekey);
                    pseudonodes.push_back(((PseudoNode*)_nodes[i]));
                } else {
                    _nodes[i] = new TreeNode(box);
                    ((TreeNode*)_nodes[i])
                            ->add_pseudoparticles(keylow, keyhigh, pindex,
                                                  level + 1, (key << ndim_) + i,
                                                  pseudonodes);
                }
            } else {
                if(_nodes[i]->is_leaf()) {
                    // this is possible: if the borderkey is a subkey of the
                    // current node, then the split is somewhat lower in the
                    // hierarchy
                    Particle* part2 = ((Leaf*)_nodes[i])->get_particle();
                    delete _nodes[i];
                    _nodes[i] = new TreeNode(box);
                    ((TreeNode*)_nodes[i])->add_particle(part2, level + 1);
                    ((TreeNode*)_nodes[i])
                            ->add_pseudoparticles(keylow, keyhigh, pindex,
                                                  level + 1, (key << ndim_) + i,
                                                  pseudonodes);
                } else {
                    ((TreeNode*)_nodes[i])
                            ->add_pseudoparticles(keylow, keyhigh, pindex,
                                                  level + 1, (key << ndim_) + i,
                                                  pseudonodes);
                }
            }
        }
    }
}

/**
 * @brief Get ExportNode information for the children of this node
 *
 * @param keylow Lower bound of the Hilbert key range of the local tree
 * @param keyhigh Upper bound of the Hilbert key range of the local tree
 * @param level Level of the node in the tree
 * @param key Key for a position inside this node
 * @param exportnodes Reference to the list to fill
 */
void TreeNode::get_exportnodes(unsigned long keylow, unsigned long keyhigh,
                               unsigned int level, unsigned long key,
                               vector<ExportNode*>& exportnodes) {
    unsigned long nodekey = key << ndim_;
    unsigned long keysize = 1;
    keysize <<= (60 - ndim_ * (level + 1));
    keysize -= 1;
    for(unsigned int i = numnode_; i--;) {
        unsigned long subnodekey = nodekey + i;
        subnodekey <<= (60 - ndim_ * (level + 1));
        if(interval_overlap(keylow, keyhigh, subnodekey,
                            subnodekey + keysize)) {
            if(_nodes[i] == NULL) {
                // this can only happen if the node is completely inside the
                // range
                // otherwise, there would be a TreeNode containing PseudoLeaves
                exportnodes.push_back(new ExportNode(NULL, subnodekey));
            } else {
                if(keylow <= subnodekey && subnodekey + keysize <= keyhigh) {
                    exportnodes.push_back(
                            new ExportNode(_nodes[i], subnodekey));
                } else {
                    // we do not encounter individual leaves that are not
                    // entirely inside the interval, since FancyPseudoNodes were
                    // added to represent the blocks on the lowest level
                    // (think about it and you will get it)
                    if(!_nodes[i]->is_leaf()) {
                        ((TreeNode*)_nodes[i])
                                ->get_exportnodes(keylow, keyhigh, level + 1,
                                                  (key << ndim_) + i,
                                                  exportnodes);
                    }
                }
            }
        }
    }
}

/**
 * @brief Recursively update the node properties
 *
 * There are two reasons to do this:
 *  - local particle information changed: all nodes have to be updated
 *  - new information was imported from another MPI process: only nodes with
 *    pseudonode children have to be updated
 *
 * @param all Flag indicating whether all nodes should be updated, or only the
 * nodes that have pseudonode children
 */
void TreeNode::update_quantities(bool all) {
    _total_mass = 0.;
    for(unsigned int i = ndim_; i--;) {
        _center_of_mass[i] = 0.;
    }
    _cmax = 0.;
    _vmax = 0.;
    _hmax = 0.;
    Node* current = _child;
    while(current) {
        Node* next;
        if(!current->is_leaf()) {
            // if the TreeNode is entirely local, the recursive algorithm stops
            // here
            if(all || !current->is_local()) {
                ((TreeNode*)current)->update_quantities(all);
            }
            next = ((TreeNode*)current)->get_sibling();
        } else {
            if(current->is_pseudo()) {
                next = ((PseudoNode*)current)->get_sibling();
            } else {
                next = ((Leaf*)current)->get_sibling();
            }
        }
        _total_mass += current->get_mass();
        for(unsigned int i = ndim_; i--;) {
            _center_of_mass[i] +=
                    current->get_mass() * current->get_center_of_mass(i);
        }
        _cmax = std::max(_cmax, current->get_cmax());
        _vmax = std::max(_vmax, current->get_vmax());
        _hmax = std::max(_hmax, current->get_hmax());
        if(current->is_last()) {
            current = NULL;
        } else {
            current = next;
        }
    }
    for(unsigned int i = ndim_; i--;) {
        _center_of_mass[i] /= _total_mass;
    }
}

/**
 * @brief Update the mesh masses of the node
 *
 * @deprecated Not used
 */
void TreeNode::update_mesh_masses() {
    _mesh_m = 0.;
    Node* current = _child;
    while(current) {
        Node* next;
        if(!current->is_leaf()) {
            ((TreeNode*)current)->update_mesh_masses();
            _mesh_m += ((TreeNode*)current)->get_mesh_m();
            next = ((TreeNode*)current)->get_sibling();
        } else {
            if(current->is_pseudo()) {
                next = ((PseudoNode*)current)->get_sibling();
            } else {
                _mesh_m += ((GasParticle*)((Leaf*)current)->get_particle())
                                   ->get_mesh_m();
                next = ((Leaf*)current)->get_sibling();
            }
        }
        if(current->is_last()) {
            current = NULL;
        } else {
            current = next;
        }
    }
}

/**
 * @brief Walk this node using the given TreeWalker
 *
 * @deprecated Replaced by a non-recursive method in Tree
 *
 * @param walker TreeWalker used to walk the tree
 */
void TreeNode::walk(TreeWalker& walker) {
    if(walker.splitnode(this)) {
        _child->walk(walker);
    } else {
        if(_sibling) {
            _sibling->walk(walker);
        }
    }
}

// Leaf

/**
 * @brief Constructor
 */
Leaf::Leaf() : Node(true) { _sibling = NULL; }

/**
 * @brief Destructor
 *
 * Clean up the next Node on the same level if necessary.
 */
Leaf::~Leaf() {
    if(!is_last()) {
        delete _sibling;
    }
}

/**
 * @brief Get the mass of the Particle in the leaf
 *
 * @return Mass of the Particle
 */
double Leaf::get_mass() { return _particle->get_mass(); }

/**
 * @brief Get the softening length of the Particle in the leaf
 *
 * @return Softening length of the Particle
 */
double Leaf::get_hmax() { return _particle->get_hsoft(); }

/**
 * @brief Get the given coordinate of the  Particle in the leaf
 *
 * @param index Requested coordinate
 * @return Requested coordinate of the Particle
 */
double Leaf::get_center_of_mass(unsigned int index) {
    return _particle->get_position()[index];
}

/**
 * @brief Get the soundspeed of the Particle in the leaf
 *
 * This of course only makes sense when the Particle is a GasParticle. If not,
 * zero is returned.
 *
 * @return The soundspeed of the Particle if it is a GasParticle
 */
double Leaf::get_cmax() {
    if(_particle->type() == PARTTYPE_GAS) {
        return ((GasParticle*)_particle)->get_soundspeed();
    } else {
        return 0.;
    }
}

/**
 * @brief Get the fluid velocity of the Particle
 *
 * This of course only makes sense when the Particle is a GasParticle. If not,
 * zero is returned.
 *
 * @return The fluid velocity of the Particle if it is a GasParticle
 */
double Leaf::get_vmax() {
    if(_particle->type() == PARTTYPE_GAS) {
        StateVector W = ((GasParticle*)_particle)->get_Wvec();
#if ndim_ == 3
        Vec v(W[1], W[2], W[3]);
#else
        Vec v(W[1], W[2]);
#endif
        return v.norm();
    } else {
        return 0.;
    }
}

/**
 * @brief Get the next Node in the efficient tree traversal algorithm
 *
 * @return Next Node on this level
 */
Node* Leaf::get_sibling() { return _sibling; }

/**
 * @brief Add the given Particle to the leaf
 *
 * @param p Particle to add
 */
void Leaf::add_particle(Particle* p) { _particle = p; }

/**
 * @brief Get the Particle in this leaf
 *
 * @return Particle in the leaf
 */
Particle* Leaf::get_particle() { return _particle; }

/**
 * @brief Restructure the leaf for the efficient tree traversal algorithm
 *
 * @param sibling Next Node on the same level in the tree traversal algorithm
 * @param last Flag indicating whether this leaf is the last child of its parent
 * (if so, there will be no next node)
 * @return 0, because we count the nodes and not the leaves
 */
unsigned int Leaf::finalize(Node* sibling, bool last) {
    _sibling = sibling;
    _flag |= (last << 1);
    return 0;
}

/**
 * @brief Print the leaf to the given stream
 *
 * @param out std::stream to write to
 */
void Leaf::print(ostream& out) {
    Vec position = _particle->get_position();
    for(unsigned int i = 0; i < ndim_; i++) {
        out << position[i] << "\t";
    }
    out << "-1\n";
}

/**
 * @brief Recursively walk the tree using the given TreeWalker
 *
 * @deprecated Replaced by a non-recursive method in Tree
 *
 * @param walker TreeWalker used to walk the tree
 */
void Leaf::walk(TreeWalker& walker) {
    walker.leafaction(this);
    if(_sibling) {
        _sibling->walk(walker);
    }
}

// PseudoNode

/**
 * @brief Constructor
 *
 * @param box Box specifying the dimensions of the embedded box
 * @param src Rank of the MPI process holding the corresponding TreeNode
 * @param key Hilbert key of this node
 */
PseudoNode::PseudoNode(Box box, unsigned int src, unsigned long key)
        : Node(true, true), Hilbert_Object(key), _box(box) {
    _sibling = NULL;
    _src = src;
    _vmax = 0.;
    _cmax = 0.;
    _total_mass = 0.;
    // flag this node as not local
    _flag |= 8;
}

/**
 * @brief Destructor
 *
 * Clean up the next Node on the same level if necessary.
 */
PseudoNode::~PseudoNode() {
    if(!is_last()) {
        delete _sibling;
    }
}

/**
 * @brief Restructure the node for the efficient tree traversal algorithm
 *
 * @param sibling Next Node on the same level
 * @param last Flag indicating whether this node is the last child of its parent
 * (if so, there will be no next node)
 * @return 0, because we only count the nodes of the tree
 */
unsigned int PseudoNode::finalize(Node* sibling, bool last) {
    _sibling = sibling;
    _flag |= (last << 1);
    return 0;
}

/**
 * @brief Get the mass of the node
 *
 * @return Mass of the node
 */
double PseudoNode::get_mass() { return _total_mass; }

/**
 * @brief Get the maximal softening length in the node
 *
 * @return Maximal softening length in the node
 */
double PseudoNode::get_hmax() { return _hmax; }

/**
 * @brief Get the given coordinate of the center of mass of the node
 *
 * @param index Requested coordinate
 * @return Requested coordinate of the center of the mass
 */
double PseudoNode::get_center_of_mass(unsigned int index) {
    return _center_of_mass[index];
}

/**
 * @brief Get the center of mass of the node
 *
 * @return Center of mass of the node
 */
Vec PseudoNode::get_center_of_mass_node() { return _center_of_mass; }

/**
 * @brief Get the next Node in the efficient tree traversal algorithm
 *
 * @return Next Node on the same level
 */
Node* PseudoNode::get_sibling() { return _sibling; }

/**
 * @brief Get the distance between the given coordinates and the box of the node
 *
 * @param coords Vec with coordinates
 * @return Distance between the node and the given coordinates
 */
double PseudoNode::get_distance(Vec& coords) {
    double distance = 0.;
    double x;
    for(unsigned int i = ndim_; i--;) {
        x = fabs(_box.get_anchor()[i] - coords[i]);
        if(x < 0.5 * _box.get_side()) {
            x = 0.;
        } else {
            x -= 0.5 * _box.get_side();
        }
        distance += x * x;
    }
    return distance;
}

/**
 * @brief Get the maximal fluid velocity in the node
 *
 * @return Maximal fluid velocity in the node
 */
double PseudoNode::get_vmax() { return _vmax; }

/**
 * @brief Get the maximal soundspeed in the node
 *
 * @return Maximal soundspeed in the node
 */
double PseudoNode::get_cmax() { return _cmax; }

/**
 * @brief Get the rank of the MPI process which holds the TreeNode corresponding
 * to this pseudonode
 *
 * @return Rank of the original MPI process of the node
 */
unsigned int PseudoNode::get_source() { return _src; }

/**
 * @brief Set the center of mass of the node
 *
 * @param com New center of mass for the node
 */
void PseudoNode::set_center_of_mass(Vec com) { _center_of_mass = com; }

/**
 * @brief Set the maximal soundspeed in the node
 *
 * @param cmax New value for the maximal soundspeed in the node
 */
void PseudoNode::set_cmax(double cmax) { _cmax = cmax; }

/**
 * @brief Set the maximal fluid velocity in the node
 * @param vmax New value for the maximal fluid velocity in the node
 */
void PseudoNode::set_vmax(double vmax) { _vmax = vmax; }

/**
 * @brief Get the width of the node
 *
 * @return The width of the node
 */
double PseudoNode::get_width() { return _box.get_side(); }

/**
 * @brief Get the center of the node
 *
 * @return Center of the node
 */
Vec PseudoNode::get_center() { return _box.get_anchor(); }

/**
 * @brief Set the total mass in the node
 *
 * @param mass New value for the total mass in the node
 */
void PseudoNode::set_mass(double mass) { _total_mass = mass; }

/**
 * @brief Set the maximal softening length in the node
 *
 * @param hmax New value for the maximal softening length in the node
 */
void PseudoNode::set_hmax(double hmax) { _hmax = hmax; }

/**
 * @brief Print the pseudonode to the given stream
 *
 * @param out std::stream to write to
 */
void PseudoNode::print(ostream& out) {
    out << _box.get_anchor()[0] << "\t" << _box.get_anchor()[1];
#if ndim_ == 3
    out << "\t" << _box.get_anchor()[2];
#endif
    out << "\t-2\n";
}

/**
 * @brief Walk the node using the given TreeWalker
 *
 * @deprecated Replaced by a non-recursive method in Tree
 *
 * @param walker TreeWalker used to walk the tree
 */
void PseudoNode::walk(TreeWalker& walker) {
    walker.pseudonodeaction(this);
    if(_sibling) {
        _sibling->walk(walker);
    }
}
