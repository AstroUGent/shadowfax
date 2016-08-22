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
 * @file FixedGrid.cpp
 *
 * @brief A fixed cartesian grid: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "FixedGrid.hpp"
#include "Error.hpp"    // for my_exit
#include "Vec.hpp"      // for Vec
#include "VorCell.hpp"  // for VorCell
#include "VorFace.hpp"  // for VorFace
#include "VorGen.hpp"   // for VorGen
#include "utilities/GasParticle.hpp"
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include <cmath>                         // for sqrt, M_PI
#include <iostream>                      // for operator<<, basic_ostream, etc
using namespace std;

/**
 * @brief Constructor
 *
 * Set up the fixed cartesian mesh based on the given ParticleVector. For now,
 * we assume a cubic grid; the size in one dimension is set by the square or
 * cubic root of the number of gasparticles in the particlevector.
 *
 * Grid cells are associated to the GasParticle that has coordinates inside that
 * particular cell.
 *
 * If the total number of gasparticles does not fit in a cubic grid or there are
 * cells with no or more than one associated GasParticle, the initialization
 * will fail.
 *
 * @param particles ParticleVector containing GasParticles with valid positions
 * @param periodic Flag to indicate if the simulation box is periodic or not
 */
FixedGrid::FixedGrid(ParticleVector& particles, bool periodic) {
// let's assume a square/cubic grid
#if ndim_ == 3
    _size[0] = int(cbrt(particles.gassize()));
    _size[1] = _size[0];
    _size[2] = _size[0];
    vector<unsigned int> igrid(_size[0] * _size[1] * _size[2],
                               _size[0] * _size[1] * _size[2]);
    _idx.resize(_size[0] * _size[1] * _size[2]);
#else
    _size[0] = int(sqrt(particles.gassize()));
    _size[1] = _size[0];
    vector<unsigned int> igrid(_size[0] * _size[1], _size[0] * _size[1]);
    _idx.resize(_size[0] * _size[1]);
#endif

    // now we need to store the pointers in the correct grid cell
    // find the cell size in every direction
    double box[ndim_ + ndim_];
    particles.get_container().get_bounding_box(box);
    double cellsize[ndim_];
    for(int i = 0; i < ndim_; i++) {
        cellsize[i] = box[i + ndim_] / _size[i];
    }

// calculate grid quantities
#if ndim_ == 3
    _cellvolume = cellsize[0] * cellsize[1] * cellsize[2];
    _cellh = cbrt(_cellvolume * 3. / 4. / M_PI);
    _cellarea = 2. * (cellsize[0] * cellsize[1] + cellsize[0] * cellsize[2] +
                      cellsize[1] * cellsize[2]);
#else
    _cellvolume = cellsize[0] * cellsize[1];
    _cellh = sqrt(_cellvolume / M_PI);
    _cellarea = 2. * (cellsize[0] + cellsize[1]);
#endif

    _vorgens.resize(particles.gassize());
    _cells.resize(particles.gassize());
    // loop over the particles and find the correct cell
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        // make a vorgen
        Vec x = particles.gas(i)->get_position();

        // store the particle in the grid
        unsigned int ijk[ndim_];
        for(int k = 0; k < ndim_; k++) {
            ijk[k] = 0;
            while(ijk[k] < _size[k] && x[k] > box[k] + ijk[k] * cellsize[k]) {
                ijk[k]++;
            }
            if(ijk[k] == _size[k] && x[k] > box[k] + ijk[k] * cellsize[k]) {
                cerr << "Particle outside grid!" << endl;
                my_exit();
            }
            // we check the upper bound of the cell, but store the lower bound
            ijk[k]--;
        }
#if ndim_ == 3
        unsigned int idx =
                ijk[0] * _size[1] * _size[2] + ijk[1] * _size[2] + ijk[2];
#else
        unsigned int idx = ijk[0] * _size[1] + ijk[1];
#endif
        _idx[i] = idx;
        igrid[idx] = i;
#if ndim_ == 3
        _vorgens[idx] = new VorGen(x.x(), x.y(), x.z());
#else
        _vorgens[idx] = new VorGen(x.x(), x.y());
#endif
        _cells[idx] = new VorCell(_vorgens[idx]);
        _cells[idx]->set_volume(_cellvolume);
        _cells[idx]->set_centroid(particles.gas(i)->get_position());
        // finalize vorgen
        _vorgens[idx]->set_particle(particles.gas(igrid[idx]));
        _vorgens[idx]->set_id(2);
    }
#if ndim_ == 3
    for(unsigned int i = 0; i < _size[0] * _size[1] * _size[2]; i++) {
        if(igrid[i] == _size[0] * _size[1] * _size[2]) {
            cerr << "Empty cell!" << endl;
            cerr << i << endl;
            my_exit();
        }
    }
#else
    for(unsigned int i = 0; i < _size[0] * _size[1]; i++) {
        if(igrid[i] == _size[0] * _size[1]) {
            cerr << "Empty cell!" << endl;
        }
    }
#endif

// set connectivity

#if ndim_ == 3
    _ghosts.resize(2 * _size[1] * _size[2] + 2 * _size[0] * _size[2] +
                   2 * _size[0] * _size[1]);

    // x direction
    for(unsigned int i = 0; i < _size[1]; i++) {
        for(unsigned int j = 0; j < _size[2]; j++) {
            Vec left = particles.gas(igrid[i * _size[2] + j])->get_position();
            Vec right =
                    particles
                            .gas(igrid[(_size[0] - 1) * _size[1] * _size[2] +
                                       i * _size[2] + j])
                            ->get_position();
            _ghosts[2 * (i * _size[2] + j)] =
                    new VorGen(2. * box[0] - left.x(), left.y(), left.z());
            _ghosts[2 * (i * _size[2] + j) + 1] = new VorGen(
                    2. * (box[0] + box[3]) - right.x(), right.y(), right.z());
            if(periodic) {
                _ghosts[2 * (i * _size[2] + j)]->set_particle(particles.gas(
                        igrid[(_size[0] - 1) * _size[1] * _size[2] +
                              i * _size[2] + j]));
                _ghosts[2 * (i * _size[2] + j) + 1]->set_particle(
                        particles.gas(igrid[i * _size[2] + j]));
            }
            _ghosts[2 * (i * _size[2] + j)]->set_id(1);
            _ghosts[2 * (i * _size[2] + j) + 1]->set_id(1);
        }
    }

    // y direction
    for(unsigned int i = 0; i < _size[0]; i++) {
        for(unsigned int j = 0; j < _size[2]; j++) {
            Vec front = particles.gas(igrid[i * _size[1] * _size[2] + j])
                                ->get_position();
            Vec back = particles
                               .gas(igrid[i * _size[1] * _size[2] +
                                          (_size[1] - 1) * _size[2] + j])
                               ->get_position();
            _ghosts[2 * _size[1] * _size[2] + 2 * (i * _size[2] + j)] =
                    new VorGen(front.x(), 2. * box[1] - front.y(), front.z());
            _ghosts[2 * _size[1] * _size[2] + 2 * (i * _size[2] + j) + 1] =
                    new VorGen(back.x(), 2. * (box[1] + box[4]) - back.y(),
                               back.z());
            if(periodic) {
                _ghosts[2 * _size[1] * _size[2] + 2 * (i * _size[2] + j)]
                        ->set_particle(particles.gas(
                                igrid[i * _size[1] * _size[2] +
                                      (_size[1] - 1) * _size[2] + j]));
                _ghosts[2 * _size[1] * _size[2] + 2 * (i * _size[2] + j) + 1]
                        ->set_particle(particles.gas(
                                igrid[i * _size[1] * _size[2] + j]));
            }
            _ghosts[2 * _size[1] * _size[2] + 2 * (i * _size[2] + j)]->set_id(
                    1);
            _ghosts[2 * _size[1] * _size[2] + 2 * (i * _size[2] + j) + 1]
                    ->set_id(1);
        }
    }

    // z direction
    for(unsigned int i = 0; i < _size[0]; i++) {
        for(unsigned int j = 0; j < _size[1]; j++) {
            Vec bottom =
                    particles
                            .gas(igrid[i * _size[1] * _size[2] + j * _size[2]])
                            ->get_position();
            Vec top = particles
                              .gas(igrid[i * _size[1] * _size[2] +
                                         j * _size[2] + _size[2] - 1])
                              ->get_position();
            _ghosts[2 * _size[1] * _size[2] + 2 * _size[0] * _size[2] +
                    2 * (i * _size[1] + j)] =
                    new VorGen(bottom.x(), bottom.y(),
                               2. * box[2] - bottom.z());
            _ghosts[2 * _size[1] * _size[2] + 2 * _size[0] * _size[2] +
                    2 * (i * _size[1] + j) + 1] =
                    new VorGen(top.x(), top.y(),
                               2. * (box[2] + box[5]) - top.z());
            if(periodic) {
                _ghosts[2 * _size[1] * _size[2] + 2 * _size[0] * _size[2] +
                        2 * (i * _size[1] + j)]
                        ->set_particle(particles.gas(
                                igrid[i * _size[1] * _size[2] + j * _size[2] +
                                      _size[2] - 1]));
                _ghosts[2 * _size[1] * _size[2] + 2 * _size[0] * _size[2] +
                        2 * (i * _size[1] + j) + 1]
                        ->set_particle(particles.gas(
                                igrid[i * _size[1] * _size[2] + j * _size[2]]));
            }
            _ghosts[2 * _size[1] * _size[2] + 2 * _size[0] * _size[2] +
                    2 * (i * _size[1] + j)]
                    ->set_id(1);
            _ghosts[2 * _size[1] * _size[2] + 2 * _size[0] * _size[2] +
                    2 * (i * _size[1] + j) + 1]
                    ->set_id(1);
        }
    }

    // create faces...
    _faces.resize((_size[0] + 1) * _size[1] * _size[2] +
                  _size[0] * (_size[1] + 1) * _size[2] +
                  _size[0] * _size[1] * (_size[2] + 1));
    // x direction
    for(unsigned int i = 0; i < _size[0] + 1; i++) {
        for(unsigned int j = 0; j < _size[1]; j++) {
            for(unsigned int k = 0; k < _size[2]; k++) {
                VorGen* left;
                VorGen* right;
                unsigned int lefti;
                unsigned int righti;
                if(!i) {
                    left = _ghosts[2 * (j * _size[2] + k)];
                    lefti = _size[0] * _size[1] * _size[2];
                } else {
                    left = _vorgens[(i - 1) * _size[1] * _size[2] +
                                    j * _size[2] + k];
                    lefti = (i - 1) * _size[1] * _size[2] + j * _size[2] + k;
                }
                if(i == _size[0]) {
                    right = _ghosts[2 * (j * _size[2] + k) + 1];
                    righti = _size[0] * _size[1] * _size[2];
                } else {
                    right = _vorgens[i * _size[1] * _size[2] + j * _size[2] +
                                     k];
                    righti = i * _size[1] * _size[2] + j * _size[2] + k;
                }
                if(lefti < _size[0] * _size[1] * _size[2]) {
                    _faces[i * _size[1] * _size[2] + j * _size[2] + k] =
                            new VorFace(lefti, _vorgens);
                    _faces[i * _size[1] * _size[2] + j * _size[2] + k]->add_ngb(
                            right);
                    _faces[i * _size[1] * _size[2] + j * _size[2] + k]
                            ->add_ngb_id(righti);
                } else {
                    _faces[i * _size[1] * _size[2] + j * _size[2] + k] =
                            new VorFace(righti, _vorgens);
                    _faces[i * _size[1] * _size[2] + j * _size[2] + k]->add_ngb(
                            left);
                    _faces[i * _size[1] * _size[2] + j * _size[2] + k]
                            ->add_ngb_id(lefti);
                }
                _faces[i * _size[1] * _size[2] + j * _size[2] + k]->set_area(
                        cellsize[1] * cellsize[2]);
                double midpoint[3] = {box[0] + i * cellsize[0],
                                      box[1] + (j + 0.5) * cellsize[1],
                                      box[2] + (k + 0.5) * cellsize[2]};
                _faces[i * _size[1] * _size[2] + j * _size[2] + k]
                        ->set_midpoint(midpoint);
                if(lefti < _size[0] * _size[1] * _size[2]) {
                    _cells[lefti]->add_ngb(right);
                    _cells[lefti]->add_face(
                            _faces[i * _size[1] * _size[2] + j * _size[2] + k]);
                }
                if(righti < _size[0] * _size[1] * _size[2]) {
                    _cells[righti]->add_ngb(left);
                    _cells[righti]->add_face(
                            _faces[i * _size[1] * _size[2] + j * _size[2] + k]);
                }
            }
        }
    }

    // y direction
    for(unsigned int i = 0; i < _size[0]; i++) {
        for(unsigned int j = 0; j < _size[1] + 1; j++) {
            for(unsigned int k = 0; k < _size[2]; k++) {
                VorGen* front;
                VorGen* back;
                unsigned int fronti;
                unsigned int backi;
                if(!j) {
                    front = _ghosts[2 * _size[1] * _size[2] +
                                    2 * (i * _size[2] + k)];
                    fronti = _size[0] * _size[1] * _size[2];
                } else {
                    front = _vorgens[i * _size[1] * _size[2] +
                                     (j - 1) * _size[2] + k];
                    fronti = i * _size[1] * _size[2] + (j - 1) * _size[2] + k;
                }
                if(j == _size[1]) {
                    back = _ghosts[2 * _size[1] * _size[2] +
                                   2 * (i * _size[2] + k) + 1];
                    backi = _size[0] * _size[1] * _size[2];
                } else {
                    back = _vorgens[i * _size[1] * _size[2] + j * _size[2] + k];
                    backi = i * _size[1] * _size[2] + j * _size[2] + k;
                }
                if(fronti < _size[0] * _size[1] * _size[2]) {
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           j * _size[0] * _size[2] + i * _size[2] + k] =
                            new VorFace(fronti, _vorgens);
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           j * _size[0] * _size[2] + i * _size[2] + k]
                            ->add_ngb(back);
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           j * _size[0] * _size[2] + i * _size[2] + k]
                            ->add_ngb_id(backi);
                } else {
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           j * _size[0] * _size[2] + i * _size[2] + k] =
                            new VorFace(backi, _vorgens);
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           j * _size[0] * _size[2] + i * _size[2] + k]
                            ->add_ngb(front);
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           j * _size[0] * _size[2] + i * _size[2] + k]
                            ->add_ngb_id(fronti);
                }
                _faces[(_size[0] + 1) * _size[1] * _size[2] +
                       j * _size[0] * _size[2] + i * _size[2] + k]
                        ->set_area(cellsize[0] * cellsize[2]);
                double midpoint[3] = {box[0] + (i + 0.5) * cellsize[0],
                                      box[1] + j * cellsize[1],
                                      box[2] + (k + 0.5) * cellsize[2]};
                _faces[(_size[0] + 1) * _size[1] * _size[2] +
                       j * _size[0] * _size[2] + i * _size[2] + k]
                        ->set_midpoint(midpoint);
                if(fronti < _size[0] * _size[1] * _size[2]) {
                    _cells[fronti]->add_ngb(back);
                    _cells[fronti]->add_face(
                            _faces[(_size[0] + 1) * _size[1] * _size[2] +
                                   j * _size[0] * _size[2] + i * _size[2] + k]);
                }
                if(backi < _size[0] * _size[1] * _size[2]) {
                    _cells[backi]->add_ngb(front);
                    _cells[backi]->add_face(
                            _faces[(_size[0] + 1) * _size[1] * _size[2] +
                                   j * _size[0] * _size[2] + i * _size[2] + k]);
                }
            }
        }
    }

    // z direction
    for(unsigned int i = 0; i < _size[0]; i++) {
        for(unsigned int j = 0; j < _size[1]; j++) {
            for(unsigned int k = 0; k < _size[2] + 1; k++) {
                VorGen* bottom;
                VorGen* top;
                unsigned int bottomi;
                unsigned int topi;
                if(!k) {
                    bottom = _ghosts[2 * _size[1] * _size[2] +
                                     2 * _size[0] * _size[2] +
                                     2 * (i * _size[1] + j)];
                    bottomi = _size[0] * _size[1] * _size[2];
                } else {
                    bottom = _vorgens[i * _size[1] * _size[2] + j * _size[2] +
                                      k - 1];
                    bottomi = i * _size[1] * _size[2] + j * _size[2] + k - 1;
                }
                if(k == _size[2]) {
                    top = _ghosts[2 * _size[1] * _size[2] +
                                  2 * _size[0] * _size[2] +
                                  2 * (i * _size[1] + j) + 1];
                    topi = _size[0] * _size[1] * _size[2];
                } else {
                    top = _vorgens[i * _size[1] * _size[2] + j * _size[2] + k];
                    topi = i * _size[1] * _size[2] + j * _size[2] + k;
                }
                if(bottomi < _size[0] * _size[1] * _size[2]) {
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           _size[0] * (_size[1] + 1) * _size[2] +
                           k * _size[0] * _size[1] + i * _size[1] + j] =
                            new VorFace(bottomi, _vorgens);
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           _size[0] * (_size[1] + 1) * _size[2] +
                           k * _size[0] * _size[1] + i * _size[1] + j]
                            ->add_ngb(top);
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           _size[0] * (_size[1] + 1) * _size[2] +
                           k * _size[0] * _size[1] + i * _size[1] + j]
                            ->add_ngb_id(topi);
                } else {
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           _size[0] * (_size[1] + 1) * _size[2] +
                           k * _size[0] * _size[1] + i * _size[1] + j] =
                            new VorFace(topi, _vorgens);
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           _size[0] * (_size[1] + 1) * _size[2] +
                           k * _size[0] * _size[1] + i * _size[1] + j]
                            ->add_ngb(bottom);
                    _faces[(_size[0] + 1) * _size[1] * _size[2] +
                           _size[0] * (_size[1] + 1) * _size[2] +
                           k * _size[0] * _size[1] + i * _size[1] + j]
                            ->add_ngb_id(bottomi);
                }
                _faces[(_size[0] + 1) * _size[1] * _size[2] +
                       _size[0] * (_size[1] + 1) * _size[2] +
                       k * _size[0] * _size[1] + i * _size[1] + j]
                        ->set_area(cellsize[0] * cellsize[1]);
                double midpoint[3] = {box[0] + (i + 0.5) * cellsize[0],
                                      box[1] + (j + 0.5) * cellsize[1],
                                      box[2] + k * cellsize[2]};
                _faces[(_size[0] + 1) * _size[1] * _size[2] +
                       _size[0] * (_size[1] + 1) * _size[2] +
                       k * _size[0] * _size[1] + i * _size[1] + j]
                        ->set_midpoint(midpoint);
                if(bottomi < _size[0] * _size[1] * _size[2]) {
                    _cells[bottomi]->add_ngb(top);
                    _cells[bottomi]->add_face(
                            _faces[(_size[0] + 1) * _size[1] * _size[2] +
                                   _size[0] * (_size[1] + 1) * _size[2] +
                                   k * _size[0] * _size[1] + i * _size[1] + j]);
                }
                if(topi < _size[0] * _size[1] * _size[2]) {
                    _cells[topi]->add_ngb(bottom);
                    _cells[topi]->add_face(
                            _faces[(_size[0] + 1) * _size[1] * _size[2] +
                                   _size[0] * (_size[1] + 1) * _size[2] +
                                   k * _size[0] * _size[1] + i * _size[1] + j]);
                }
            }
        }
    }
#else
    _ghosts.resize(2 * _size[0] + 2 * _size[1]);

    for(unsigned int i = 0; i < _size[0]; i++) {
        Vec bottom = particles.gas(igrid[i * _size[1]])->get_position();
        Vec top = particles.gas(igrid[(i + 1) * _size[1] - 1])->get_position();
        _ghosts[2 * i] = new VorGen(bottom.x(), 2. * box[1] - bottom.y());
        _ghosts[2 * i + 1] =
                new VorGen(top.x(), 2. * (box[1] + box[3]) - top.y());
        // we only set the particle if the box is periodic
        // if we do not set a particle, the boundary is treated as being
        // reflective
        if(periodic) {
            _ghosts[2 * i]->set_particle(
                    particles.gas(igrid[(i + 1) * _size[1] - 1]));
            _ghosts[2 * i + 1]->set_particle(
                    particles.gas(igrid[i * _size[1]]));
        }
        _ghosts[2 * i]->set_id(1);
        _ghosts[2 * i + 1]->set_id(1);
    }

    for(unsigned int i = 0; i < _size[1]; i++) {
        Vec left = particles.gas(igrid[i])->get_position();
        Vec right = particles.gas(igrid[(_size[0] - 1) * _size[1] + i])
                            ->get_position();
        _ghosts[2 * _size[0] + 2 * i] =
                new VorGen(2. * box[0] - left.x(), left.y());
        _ghosts[2 * _size[0] + 2 * i + 1] =
                new VorGen(2. * (box[0] + box[2]) - right.x(), right.y());
        if(periodic) {
            _ghosts[2 * _size[0] + 2 * i]->set_particle(
                    particles.gas(igrid[(_size[0] - 1) * _size[1] + i]));
            _ghosts[2 * _size[0] + 2 * i + 1]->set_particle(
                    particles.gas(igrid[i]));
        }
        _ghosts[2 * _size[0] + 2 * i]->set_id(1);
        _ghosts[2 * _size[0] + 2 * i + 1]->set_id(1);
    }

    // create faces
    _faces.resize((_size[0] + 1) * _size[1] + _size[0] * (_size[1] + 1));
    for(unsigned int k = 0; k < _size[1]; k++) {
        for(unsigned int i = 0; i < _size[0] + 1; i++) {
            VorGen* left;
            VorGen* right;
            unsigned int lefti;
            unsigned int righti;
            if(!i) {
                left = _ghosts[2 * _size[0] + 2 * k];
                lefti = _size[0] * _size[1];
            } else {
                left = _vorgens[(i - 1) * _size[1] + k];
                lefti = (i - 1) * _size[1] + k;
            }
            if(i == _size[0]) {
                right = _ghosts[2 * _size[0] + 2 * k + 1];
                righti = _size[0] * _size[1];
            } else {
                right = _vorgens[i * _size[1] + k];
                righti = i * _size[1] + k;
            }
            if(lefti < _size[0] * _size[1]) {
                _faces[i * _size[1] + k] = new VorFace(lefti, _vorgens);
                _faces[i * _size[1] + k]->add_ngb(right);
                _faces[i * _size[1] + k]->add_ngb_id(righti);
            } else {
                _faces[i * _size[1] + k] = new VorFace(righti, _vorgens);
                _faces[i * _size[1] + k]->add_ngb(left);
                _faces[i * _size[1] + k]->add_ngb_id(lefti);
            }
            _faces[i * _size[1] + k]->set_area(cellsize[1]);
            double midpoint[2] = {box[0] + i * cellsize[0],
                                  box[1] + (k + 0.5) * cellsize[1]};
            _faces[i * _size[1] + k]->set_midpoint(midpoint);
            if(lefti < _size[0] * _size[1]) {
                _cells[lefti]->add_ngb(right);
                _cells[lefti]->add_face(_faces[i * _size[1] + k]);
            }
            if(righti < _size[0] * _size[1]) {
                _cells[righti]->add_ngb(left);
                _cells[righti]->add_face(_faces[i * _size[1] + k]);
            }
        }
    }
    for(unsigned int k = 0; k < _size[0]; k++) {
        for(unsigned int i = 0; i < _size[1] + 1; i++) {
            VorGen* bottom;
            VorGen* top;
            unsigned int bottomi;
            unsigned int topi;
            if(!i) {
                bottom = _ghosts[2 * k];
                bottomi = _size[0] * _size[1];
            } else {
                bottom = _vorgens[k * _size[1] + (i - 1)];
                bottomi = k * _size[1] + (i - 1);
            }
            if(i == _size[1]) {
                top = _ghosts[2 * k + 1];
                topi = _size[0] * _size[1];
            } else {
                top = _vorgens[k * _size[1] + i];
                topi = k * _size[1] + i;
            }
            if(bottomi < _size[0] * _size[1]) {
                _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k] =
                        new VorFace(bottomi, _vorgens);
                _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k]->add_ngb(
                        top);
                _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k]
                        ->add_ngb_id(topi);
            } else {
                _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k] =
                        new VorFace(topi, _vorgens);
                _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k]->add_ngb(
                        bottom);
                _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k]
                        ->add_ngb_id(bottomi);
            }
            _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k]->set_area(
                    cellsize[0]);
            double midpoint[2] = {box[0] + (k + 0.5) * cellsize[0],
                                  box[1] + i * cellsize[1]};
            _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k]->set_midpoint(
                    midpoint);
            if(bottomi < _size[0] * _size[1]) {
                _cells[bottomi]->add_ngb(top);
                _cells[bottomi]->add_face(
                        _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k]);
            }
            if(topi < _size[0] * _size[1]) {
                _cells[topi]->add_ngb(bottom);
                _cells[topi]->add_face(
                        _faces[(_size[0] + 1) * _size[1] + i * _size[0] + k]);
            }
        }
    }
#endif

    cout << "Fixed grid constructed!" << endl;
}

/**
 * @brief Destructor
 *
 * Clean up generators, ghosts, faces and cells.
 */
FixedGrid::~FixedGrid() {
    for(unsigned int i = 0; i < _vorgens.size(); i++) {
        delete _vorgens[i];
    }
    for(unsigned int i = 0; i < _ghosts.size(); i++) {
        delete _ghosts[i];
    }
    for(unsigned int i = 0; i < _faces.size(); i++) {
        delete _faces[i];
    }
    for(unsigned int i = 0; i < _cells.size(); i++) {
        delete _cells[i];
    }
}

/**
 * @brief Get the characteristic length of any cell in the grid
 *
 * The characteristic length is the radius of a sphere with the same volume as
 * a cell. Since this radius is constant and the same for all cells, it is
 * precomputed once and then just returned.
 *
 * @return The characteristic length of any cell in the grid
 */
double FixedGrid::get_h() {
    return _cellh;
}

/**
 * @brief Get the volume of any cell in the grid
 *
 * Since the volume is constant and the same for all cells, it is precomputed
 * once and then just returned.
 *
 * @return The volume of any cell in the grid
 */
double FixedGrid::get_volume() {
    return _cellvolume;
}

/**
 * @brief Get the total surface area of any cell in the grid
 *
 * Since the surface area is constant and the same for all cells, it is
 * precomputed and then just returned
 *
 * @return The total surface area of any cell in the grid
 */
double FixedGrid::get_area() {
    return _cellarea;
}

/**
 * @brief Calculate hydrodynamical fluxes for all faces in the grid
 *
 * @param timeline TimeLine of the simulation, used to convert integer time to
 * real time
 * @param solver RiemannSolver used to solve the Riemann problem at the faces
 */
void FixedGrid::hydro(TimeLine& timeline, RiemannSolver& solver) {
    for(unsigned int i = 0; i < _faces.size(); i++) {
        VorFace* face = _faces[i];
        face->calculate_flux(timeline, solver);
    }
}

/**
 * @brief Estimate the gradients for the cell with the given index
 *
 * The index is actually that of the associated GasParticle, to obtain the index
 * of this cell in the internal grid, we use an internal lookup table.
 *
 * @param index Index of the GasParticle associated with the cell for which we
 * want to compute gradients
 * @param delta Array to store the resulting gradients in
 */
void FixedGrid::get_gradients(unsigned int index, StateVector* delta) {
    _cells[_idx[index]]->estimate_gradient(delta);
}
