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
 * @file VorTessManager.cpp
 *
 * @brief General interface for Voronoi tesselations: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "VorTessManager.hpp"
#include "ProgramLog.hpp"                // for LOGS
#include "RestartFile.hpp"               // for RestartFile
#include "VorCell.hpp"                   // for VorCell
#include "VorFace.hpp"                   // for VorFace
#include "VorTess.hpp"                   // for VorTess
#include "utilities/GasParticle.hpp"     // for GasParticle
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include "utilities/Timer.hpp"           // for Timer
#include <iostream>                      // for operator<<, ostream, etc
#include <vector>                        // for vector
using namespace std;

//#define FIXED_GRID

#ifndef FIXED_GRID

/**
 * @brief Constructor
 *
 * @param particles ParticleVector for which the tesselation is constructed
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective (false)
 * @param tolerance Tolerance value used to distinguish between standard
 * floating point arithmetics and arbitrary precision arithmetics in geometry
 * tests
 * @param evolve Flag indicating if the grid should be evolved rather than
 * reconstructed for every timestep (currently ignored; grid is always
 * reconstructed)
 */
VorTessManager::VorTessManager(ParticleVector& particles, bool periodic,
                               double tolerance, bool evolve)
        : _particles(particles), _periodic(periodic), _tolerance(tolerance) {
    _gridtimer.start();
    _tesselation = new VorTess(&particles.get_container(), particles.gassize(),
                               periodic, tolerance);
    LOGS("Start adding points to VorTess");
    for(unsigned int i = 0; i < _particles.gassize(); i++) {
        _tesselation->add_point(_particles.gas(i),
                                _particles.gas(i)->get_local_id());
        _particles.gas(i)->reset_copies();
        _particles.gas(i)->reset_export();
    }
    LOGS("Done adding points.");
    _tesselation->complete(_particles.get_tree());
    _tesselation->construct();
    LOGS("VorTess ready");

    _gridtimer.stop();
}

/**
 * @brief Destructor
 *
 * Clean up. Print the timers to the stdout.
 */
VorTessManager::~VorTessManager() {
    delete _tesselation;
    cout << "Spent " << _gridtimer.value() << "s in grid construction" << endl;
    cout << "Spent " << _hydrotimer.value() << "s in hydro related functions"
         << endl;
}

/**
 * @brief Reset the tesselation
 *
 * For the old algorithm, this deletes the old tesselation and initializes a new
 * one.
 *
 * @param cont DelCont specifying the dimensions of the simulation box
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective (false)
 * @param tolerance Tolerance value used to distinguish between standard
 * floating point arithmetics and arbitrary precision arithmetics in geometry
 * tests
 */
void VorTessManager::reset(DelCont* cont, bool periodic, double tolerance) {
    delete _tesselation;
    _tesselation = new VorTess(cont, _particles.gassize(), periodic, tolerance);
}

/**
 * @brief Update the tesselation to the given current system time
 *
 * For the old algorithm, this adds all points to the new tesselation and
 * constructs the new grid.
 *
 * @param currentTime Current integer system time
 * @return Number of active cells in the tesselation
 */
unsigned int VorTessManager::update(unsigned long currentTime) {
    _gridtimer.start();
    unsigned int count = 0;
    for(unsigned int i = 0; i < _particles.gassize(); i++) {
        if(_particles.gas(i)->get_endtime() == currentTime) {
            _tesselation->add_point(_particles.gas(i),
                                    _particles.gas(i)->get_local_id());
            count++;
        }
        _particles.gas(i)->reset_copies();
        _particles.gas(i)->reset_export();
    }
    _tesselation->complete(_particles.get_tree());
    _tesselation->construct();
    _gridtimer.stop();
    return count;
}

/**
 * @brief Communicate primitive variables between MPI processes
 *
 * Currently only done for the old algorithm.
 */
void VorTessManager::update_Ws() {
    _tesselation->update_Ws();
}

/**
 * @brief Calculate characteristic lengths for the active cells of the
 * tesselation
 *
 * The characteristic length is the radius of a sphere/circle with the same
 * volume/face area as the cell
 *
 * @param currentTime Current integer system time
 */
void VorTessManager::set_hs(unsigned long currentTime) {
    for(unsigned int i = 0; i < _particles.gassize(); i++) {
        if(_particles.gas(i)->get_endtime() == currentTime) {
            double h = _particles.gas(i)->get_cell()->get_h();
            _particles.gas(i)->set_h(h);
        }
    }
}

/**
 * @brief Communicate gradients between MPI processes
 *
 * Currently only done for the old algorithm.
 */
void VorTessManager::update_gradients() {
    _tesselation->update_gradients();
}

/**
 * @brief Communicate timesteps between MPI processes
 *
 * Currently only done for the old algorithm.
 *
 * @param currentTime Current integer system time
 */
void VorTessManager::update_dts(unsigned long currentTime) {
    _tesselation->update_dts(currentTime);
}

/**
 * @brief Communicate gravitational correction terms between MPI processes
 *
 * Currently only done for the old algorithm.
 */
void VorTessManager::update_gravitational_corrections() {
    _tesselation->update_gravitational_corrections();
}

/**
 * @brief Do the hydrodynamical integration by calculating the flux exchanges
 * between cells through the faces of the tesselation
 *
 * @param timeline Simulation TimeLine
 * @param solver RiemannSolver used to solve the Riemann problem at the faces
 * @param particles ParticleVector containing the particles of the simulation
 */
void VorTessManager::hydro(TimeLine& timeline, RiemannSolver& solver,
                           ParticleVector& particles) {
    _gridtimer.start();
    vector<VorFace*>& faces = _tesselation->get_faces();
    cout << "Faces size: " << faces.size() << endl;
    _gridtimer.stop();
    _hydrotimer.start();
    for(unsigned int i = 0; i < faces.size(); i++) {
        faces[i]->set_v(_particles);
#ifndef ICMAKER
        faces[i]->calculate_flux(timeline, solver);
#endif
    }
    _hydrotimer.stop();
}

/**
 * @brief Communicate fluxes between MPI processes
 *
 * Currently only done for the old algorithm.
 */
void VorTessManager::update_dQs() {
    _tesselation->update_dQs();
}

/**
 * @brief Print the tesselation to the given stream in a format that can be
 * plotted using gnuplot
 *
 * @param stream std::ostream to write to
 */
void VorTessManager::print_tesselation_gnuplot(ostream& stream) {
    _tesselation->print_tesselation_gnuplot(stream);
}

/**
 * @brief Dump the tesselation to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void VorTessManager::dump(RestartFile& rfile) {
    _gridtimer.dump(rfile);
    _hydrotimer.dump(rfile);
    rfile.write(_periodic);
    rfile.write(_tolerance);
}

/**
 * @brief Restart constructor. Initialize the tesselation from the given
 * RestartFile
 *
 * Fot the old algorithm, we reconstruct the grid, since this is easier than
 * storing it and reading it in again.
 *
 * @param rfile RestartFile to read from
 * @param particles ParticleVector containing the particles of the simulation
 */
VorTessManager::VorTessManager(RestartFile& rfile, ParticleVector& particles)
        : _particles(particles), _gridtimer(rfile), _hydrotimer(rfile) {
    rfile.read(_periodic);
    rfile.read(_tolerance);
    _tesselation = new VorTess(&particles.get_container(), particles.gassize(),
                               _periodic, _tolerance);
    for(unsigned int i = 0; i < _particles.gassize(); i++) {
        _tesselation->add_point(_particles.gas(i),
                                _particles.gas(i)->get_local_id());
        _particles.gas(i)->reset_copies();
        _particles.gas(i)->reset_export();
    }
    _tesselation->complete(_particles.get_tree());
    _tesselation->construct();
}

/**
 * @brief Get the velocity of the cell corresponding to the particle with the
 * given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return The velocity of the cell corresponding to the given particle
 */
Vec VorTessManager::get_velocity(unsigned int index) {
    VorCell* cell = _particles.gas(index)->get_cell();
    return cell->get_velocity();
}

/**
 * @brief Estimate slope-limited gradients for the cell corresponding to the
 * particle with the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @param delta Array to store the resulting gradients in
 */
void VorTessManager::estimate_gradients(unsigned int index,
                                        StateVector* delta) {
    _hydrotimer.start();
    VorCell* cell = _particles.gas(index)->get_cell();
    cell->estimate_gradient(delta);
    _hydrotimer.stop();
}

/**
 * @brief Calculate the volume of the cell corresponding to the particle with
 * the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Volume of the cell corresponding to the given particle
 */
double VorTessManager::get_volume(unsigned int index) {
    VorCell* cell = _particles.gas(index)->get_cell();
    return cell->get_volume();
}

/**
 * @brief Get the total surface area of the Voronoi cell with the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Total surface area of the Voronoi cell corresponding to the given
 * particle
 */
double VorTessManager::get_total_area(unsigned int index) {
    VorCell* cell = _particles.gas(index)->get_cell();
    return cell->get_total_area();
}

/**
 * @brief Calculate the centroid of the cell corresponding to the particle with
 * the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Centroid of the cell corresponding to the given particle
 */
Vec VorTessManager::get_centroid(unsigned int index) {
    VorCell* cell = _particles.gas(index)->get_cell();
    return cell->get_centroid();
}

/**
 * @brief Update the positions of the generators of the tesselation after a
 * movement of the simulation particles
 *
 * For the old algorithm, we only make sure all generators stay inside the
 * simulation box (by trimming the coordinates for a periodic box).
 *
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective
 */
void VorTessManager::update_positions(bool periodic) {
    _gridtimer.start();
    if(periodic) {
        for(unsigned int i = 0; i < _particles.gassize(); i++) {
            // make sure the particle stays inside the periodic box
            _particles.get_container().keep_inside(_particles.gas(i));
        }
    }
    _gridtimer.stop();
}

/**
 * @brief Get the correction term to the gravitational acceleration due to the
 * use of variable softening lengths for the GasParticle with the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Correction to the acceleration
 */
Vec VorTessManager::get_gravitational_correction(unsigned int index) {
    VorCell* cell = _particles.gas(index)->get_cell();
    return cell->get_gravitational_correction();
}

/**
 * @brief Write connectivity information to the given stream
 *
 * This information can be used to compare tesselations constructed with
 * different algorithms.
 *
 * @param stream std::ostream to write to
 * @param currentTime Current integer system time
 */
void VorTessManager::dump_connectivity(ostream& stream,
                                       unsigned long currentTime) {
    for(unsigned int i = 0; i < 3 * _particles.gassize(); i++) {
        unsigned int j = 0;
        while(j < _particles.gassize() && _particles.gas(j)->id() != i) {
            j++;
        }
        if(j == _particles.gassize()) {
            continue;
        }
        if(_particles.gas(j)->get_endtime() == currentTime) {
            unsigned long id = _particles.gas(j)->id();
            char c = 'c';
            stream.write(&c, 1);
            stream.write(reinterpret_cast<char*>(&id), sizeof(unsigned long));
            vector<unsigned long> ids =
                    _tesselation->get_ngb_ids(_particles.gas(j)->get_cell());
            for(unsigned int k = 0; k < ids.size(); k++) {
                stream.write(reinterpret_cast<char*>(&ids[k]),
                             sizeof(unsigned long));
            }
        }
    }
}

/**
 * @brief Print cell statistics for the Voronoi mesh
 */
void VorTessManager::print_statistics() {
    _tesselation->print_cell_statistics(cout);
}

#else  // FIXED_GRID

/**
 * @brief Constructor
 *
 * @param particles ParticleVector for which the tesselation is constructed
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective (false)
 * @param tolerance Tolerance value used to distinguish between standard
 * floating point arithmetics and arbitrary precision arithmetics in geometry
 * tests
 * @param evolve Flag indicating if the grid should be evolved rather than
 * reconstructed for every timestep (currently ignored; grid is always
 * reconstructed)
 */
VorTessManager::VorTessManager(ParticleVector& particles, bool periodic,
                               double tolerance, bool evolve)
        : _particles(particles) {
    _fixedgrid = new FixedGrid(particles, periodic);
}

/**
 * @brief Destructor
 *
 * Clean up the FixedGrid.
 */
VorTessManager::~VorTessManager() {
    delete _fixedgrid;
}

/**
 * @brief Reset the tesselation
 *
 * Nothing is done, since the grid is fixed.
 *
 * @param cont DelCont specifying the dimensions of the simulation box
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective (false)
 * @param tolerance Tolerance value used to distinguish between standard
 * floating point arithmetics and arbitrary precision arithmetics in geometry
 * tests
 */
void VorTessManager::reset(DelCont* cont, bool periodic, double tolerance) {
    // nothing needs to happen, because the grid is fixed
}

/**
 * @brief Update the tesselation
 *
 * Nothing is done, since the grid is fixed.
 *
 * @param currentTime Current integer system time
 * @return 0, because nothing changes
 */
unsigned int VorTessManager::update(unsigned long currentTime) {
    // nothing happens, the grid is fixed
    return 0;
}

/**
 * @brief Communicate primitive variables between processes
 *
 * @warning Not implemented for fixed grid!
 */
void VorTessManager::update_Ws() {
    // do communication
}

/**
 * @brief Calculate characteristic lengths for the active cells of the
 * tesselation
 *
 * The characteristic length is the radius of a sphere/circle with the same
 * volume/face area as the cell
 *
 * @param currentTime Current integer system time
 */
void VorTessManager::set_hs(unsigned long currentTime) {
    for(unsigned int i = 0; i < _particles.gassize(); i++) {
        if(_particles.gas(i)->get_endtime() == currentTime) {
            double h = _fixedgrid->get_h();
            _particles.gas(i)->set_h(h);
        }
    }
}

/**
 * @brief Communicate gradients between MPI processes
 *
 * @warning Not implemented for fixed grid!
 */
void VorTessManager::update_gradients() {
    // do communication
}

/**
 * @brief Communicate timesteps between MPI processes
 *
 * @warning Not implemented for fixed grid!
 *
 * @param currentTime Current integer system time
 */
void VorTessManager::update_dts(unsigned long currentTime) {
    // do communication
}

/**
 * @brief Do the hydrodynamical integration by calculating the flux exchanges
 * between cells through the faces of the tesselation
 *
 * @param timeline Simulation TimeLine
 * @param solver RiemannSolver used to solve the Riemann problem at the faces
 * @param particles ParticleVector containing the particles of the simulation
 */
void VorTessManager::hydro(TimeLine& timeline, RiemannSolver& solver,
                           ParticleVector& particles) {
    _fixedgrid->hydro(timeline, solver);
}

/**
 * @brief Communicate fluxes between MPI processes
 *
 * @warning Not implemented for fixed grid!
 */
void VorTessManager::update_dQs() {
    // do communication
}

/**
 * @brief Print the tesselation to the given stream in a format that can be
 * plotted using gnuplot
 *
 * @warning Not implemented for fixed grid!
 *
 * @param stream std::ostream to write to
 */
void VorTessManager::print_tesselation_gnuplot(std::ostream& stream) {
    // do printing
}

/**
 * @brief Dump the tesselation to the given RestartFile
 *
 * @warning Not implemented for fixed grid!
 *
 * @param rfile RestartFile to write to
 */
void VorTessManager::dump(RestartFile& rfile) {
    // do dumping
}

/**
 * @brief Restart constructor. Initialize the tesselation from the given
 * RestartFile
 *
 * @warning Not implemented for fixed grid!
 *
 * @param rfile RestartFile to write to
 * @param particles ParticleVector containing the particles of the simulation
 */
VorTessManager::VorTessManager(RestartFile& rfile, ParticleVector& particles)
        : _particles(particles) {
    // do restart initialization
}

/**
 * @brief Get the velocity of the cell corresponding to the particle with the
 * given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return A zero velocity, because the grid is fixed
 */
Vec VorTessManager::get_velocity(unsigned int index) {
    // return zero velocity; the grid is fixed
    return Vec();
}

/**
 * @brief Estimate slope-limited gradients for the cell corresponding to the
 * particle with the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @param delta Array to store the resulting gradients in
 */
void VorTessManager::estimate_gradients(unsigned int index,
                                        StateVector* delta) {
    _fixedgrid->get_gradients(index, delta);
}

/**
 * @brief Calculate the volume of the cell corresponding to the particle with
 * the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Volume of the cell corresponding to the given particle
 */
double VorTessManager::get_volume(unsigned int index) {
    return _fixedgrid->get_volume();
}

/**
 * @brief Calculate the centroid of the cell corresponding to the particle with
 * the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Centroid of the cell corresponding to the given particle
 */
Vec VorTessManager::get_centroid(unsigned int index) {
    // the centroid is the actual position of the particle
    Vec centroid = _particles.gas(index)->get_position();
    return centroid;
}

/**
 * @brief Update the positions of the generators of the tesselation after a
 * movement of the simulation particles
 *
 * Does nothing, since the grid is fixed.
 *
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective
 */
void VorTessManager::update_positions(bool periodic) {
    // do nothing
}

/**
 * @brief Get the correction term to the gravitational acceleration due to the
 * use of variable softening lengths for the GasParticle with the given index
 *
 * @warning Not implemented for fixed grid!
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Correction to the acceleration
 */
Vec VorTessManager::get_gravitational_correction(unsigned int index) {
    // currently not used
    return Vec();
}

/**
 * @brief Write connectivity information to the given stream
 *
 * @warning Not implemented for fixed grid!
 *
 * @param stream std::ostream to write to
 * @param currentTime Current integer system time
 */
void VorTessManager::dump_connectivity(std::ostream& stream,
                                       unsigned long currentTime) {
    // not used
}

#endif  // FIXED_GRID
