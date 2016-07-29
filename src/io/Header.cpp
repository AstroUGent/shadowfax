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
 * @file Header.cpp
 *
 * @brief Snapshot header: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Header.hpp"
#include "Error.hpp"
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include "RestartFile.hpp"
#include "StateVector.hpp"
#include <iostream>
using namespace std;

/**
 * @brief Constructor
 *
 * Initialize variables based on code compilation flags.
 * Time and particle number are set to 0, global timestep is disabled.
 * The simulation box is set to a point in the origin of the coordinate system
 * (so only 0's).
 */
Header::Header() {
    _npart = 0;
    for(unsigned int i = PARTTYPE_COUNTER; i--;) {
        _npartspec[i] = 0;
    }
    _ndim = ndim_;
    for(unsigned int i = ndim_ + ndim_; i--;) {
        _box[i] = 0.;
    }
    _time = 0.;
    _periodic = false;
#ifdef NOMUSCL
    _second_order = false;
#else
    _second_order = true;
#endif
#ifdef STATIC
    _static = true;
#else
    _static = false;
#endif
    _global_timestep = false;
    _gamma = GAMMA;
    _gravity = false;
    _hsoft = 0.;
}

/**
  * @brief Set the number of gas particles in the simulation
  *
  * @param npart The number of gas particles in the simulation
  */
void Header::set_ngaspart(unsigned int npart) {
    _npartspec[PARTTYPE_GAS] = npart;
    _npart += npart;
}

/**
  * @brief Set the number of DM particles in the simulation
  *
  * @param npart The number of DM particles in the simulation
  */
void Header::set_ndmpart(unsigned int npart) {
    _npartspec[PARTTYPE_DM] = npart;
    _npart += npart;
}

/**
  * @brief Set the number of star particles in the simulation
  *
  * @param npart The number of star particles in the simulation
  */
void Header::set_nstarpart(unsigned int npart) {
    _npartspec[PARTTYPE_STAR] = npart;
    _npart += npart;
}

/**
  * @brief Set the simulation box
  *
  * The first 3 (or 2) coordinates are the origin of the box, the last
  * coordinate is the side of the (cubic) box.
  *
  * @param box A 3- or 4-element array specifying the origin and side of the
  * simulation box
  */
void Header::set_box(double* box) {
    for(unsigned int i = ndim_ + ndim_; i--;) {
        _box[i] = box[i];
    }
}

/**
  * @brief Set the simulation time
  *
  * @param time The current time of the simulation
  */
void Header::set_time(double time) {
    _time = time;
}

/**
  * @brief Set the global timestep flag
  *
  * @param global_timestep Boolean specifying if individual timesteps (0) or a
  * global timestep (1) is used
  */
void Header::set_global_timestep(bool global_timestep) {
    _global_timestep = global_timestep;
}

/**
 * @brief Set the box periodicity flag
 *
 * @param periodic Boolean specifying if the simulation box is periodic (true)
 * or reflective (false)
 */
void Header::set_periodic(bool periodic) {
    _periodic = periodic;
}

/**
 * @brief Set the gravity flag
 *
 * @param gravity Boolean indicating whether we use gravity or not
 */
void Header::set_gravity(bool gravity) {
    _gravity = gravity;
}

/**
 * @brief Set the global softening length
 *
 * @param hsoft Softening length used if no individual softening lenght is
 * specified for a Particle
 */
void Header::set_hsoft(double hsoft) {
    _hsoft = hsoft;
}

/**
  * @brief Get a pointer to buffer containing the number of particles per type,
  * as used by C-style methods to write to files.
  *
  * @return A pointer to the memory location holding the number of particles per
  * type
  */
void* Header::get_npartspec() {
    return &_npartspec;
}

/**
  * @brief Get a pointer to buffer containing the number of dimensions, as used
  * by C-style methods to write to files.
  *
  * @return A pointer to the memory location holding the number of dimensions
  */
void* Header::get_ndim() {
    return &_ndim;
}

/**
  * @brief Get a pointer to buffer containing the simulation box, as used by
  * C-style methods to write to files.
  *
  * @return A pointer to the memory location holding the simulation box
  */
void* Header::get_box() {
    return _box;
}

/**
  * @brief Get a pointer to buffer containing the simulation time, as used by
  * C-style methods to write to files.
  *
  * @return A pointer to the memory location holding the simulation time
  */
void* Header::get_time() {
    return &_time;
}

/**
  * @brief Get a pointer to buffer containing the periodic boundary flag, as
  * used by C-style methods to write to files.
  *
  * @return A pointer to the memory location holding the periodic boundary flag
  */
void* Header::get_periodic() {
    return &_periodic;
}

/**
  * @brief Get a pointer to buffer containing the second order flag, as used by
  * C-style methods to write to files.
  *
  * @return A pointer to the memory location holding the second order flag
  */
void* Header::get_second_order() {
    return &_second_order;
}

/**
  * @brief Get a pointer to buffer containing the static flag, as used by
  * C-style methods to write to files.
  *
  * @return A pointer to the memory location holding the static flag
  */
void* Header::get_static() {
    return &_static;
}

/**
  * @brief Get a pointer to buffer containing the global timestep flag, as used
  * by C-style methods to write to files.
  *
  * @return A pointer to the memory location holding the global timestep flag
  */
void* Header::get_global_timestep() {
    return &_global_timestep;
}

/**
  * @brief Get a pointer to buffer containing the adiabatic index, as used by
  * C-style methods to write to files.
  *
  * @return A pointer to the memory location holding the adiabatic index
  */
void* Header::get_gamma() {
    return &_gamma;
}

/**
 * @brief Get a pointer to buffer containing the gravity flag, as used by
 * C-style methods to write to files
 *
 * @return A pointer to the memory location holding the gravity flag
 */
void* Header::get_gravity() {
    return &_gravity;
}

/**
 * @brief Get a pointer to buffer containing the softening length, as used by
 * C-style methods to write to files
 *
 * @return A pointer to the memory location holding the softening length
 */
void* Header::get_hsoft() {
    return &_hsoft;
}

/**
  * @brief Check if the compile flag options and the Header variables are
  * compatible
  *
  * If not, issue warnings or errors and in some cases crash.
  */
void Header::check_makeflags() {
    if(_ndim != ndim_) {
        cerr << "Error: dimensions of code and snapshot don't match (code was "
                "compiled with ndim_="
             << ndim_ << ", while snapshot contains "
                         "ndim="
             << _ndim << ")!" << endl;
        my_exit();
    }
#ifdef NOMUSCL
    if(_second_order) {
        cerr << "Warning: code compiled with first order accuracy, snapshot "
                "has second order accuracy"
             << endl;
    }
#else
    if(!_second_order) {
        cerr << "Warning: code compiled with second order accuracy, snapshot "
                "has only first order accuracy"
             << endl;
    }
#endif
#ifdef STATIC
    if(!_static) {
        cerr << "Warning: code compiled with static mesh, snapshot contains "
                "moving mesh"
             << endl;
    }
#else
    if(_static) {
        cerr << "Warning: code compiled with moving mesh, snapshot contains "
                "static mesh"
             << endl;
    }
#endif
    if(_gamma != GAMMA) {
        cerr << "Error: code compiled with gamma=" << GAMMA
             << ", while the "
                "snapshot contains gamma="
             << _gamma << "!" << endl;
        my_exit();
    }
}

/**
  * @brief Get the number of particles
  *
  * @return The number of particles in the simulation
  */
unsigned int Header::npart() {
    if(!_npart) {
        _npart = _npartspec[0] + _npartspec[1];
    }
    return _npart;
}

/**
  * @brief Get the number of gas particles
  *
  * @return The number of gas particles in the simulation
  */
unsigned int Header::ngaspart() {
    return _npartspec[0];
}

/**
  * @brief Get the number of DM particles
  *
  * @return The number of DM particles in the simulation
  */
unsigned int Header::ndmpart() {
    return _npartspec[1];
}

/**
  * @brief Get the number of star particles
  *
  * @return The number of star particles in the simulation
  */
unsigned int Header::nstarpart() {
    return _npartspec[PARTTYPE_STAR];
}

/**
  * @brief Get the simulation box
  *
  * @param box A 4- or 3-element array to fill
  */
void Header::box(double* box) {
    for(unsigned int i = ndim_ + ndim_; i--;) {
        box[i] = _box[i];
    }
}

/**
  * @brief Get the simulation time
  *
  * @return The simulation time
  */
double Header::time() {
    return _time;
}

/**
  * @brief Get the global timestep flag
  *
  * @return true if a global timestep is used, false if individual timesteps are
  * desired
  */
bool Header::global_timestep() {
    return _global_timestep & 1;
}

/**
 * @brief Get the periodic simulation box flag
 *
 * @return true if the simulation box is periodic, false if it is reflective
 */
bool Header::periodic() {
    return _periodic & 1;
}

/**
 * @brief Get the gravity flag
 *
 * @return true if gravity is used, false otherwise
 */
bool Header::gravity() {
    return _gravity;
}

/**
 * @brief Get the general softening length
 *
 * @return The general softening length used for particles that do not have an
 * individual softening length
 */
double Header::hsoft() {
    return _hsoft;
}

/**
 * @brief Pack data to MPI buffer for communication
 *
 * @param buffer MPI buffer
 * @param bufsize Buffer size
 * @param position Position in the buffer (is updated)
 */
void Header::pack_data(void* buffer, int bufsize, int* position) {
    MyMPI_Pack(&_npart, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_npartspec[0], PARTTYPE_COUNTER, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_ndim, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_box[0], ndim_ + ndim_, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_time, 1, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_periodic, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_second_order, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_static, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_global_timestep, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_gamma, 1, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_gravity, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_hsoft, 1, MPI_DOUBLE, buffer, bufsize, position);
}

/**
 * @brief MPI constructor. Initialize the header based on the given MPI buffer
 *
 * @param buffer MPI buffer
 * @param bufsize Buffer size
 * @param position Position in buffer (is updated)
 */
Header::Header(void* buffer, int bufsize, int* position) {
    MyMPI_Unpack(buffer, bufsize, position, &_npart, 1, MPI_UNSIGNED);
    MyMPI_Unpack(buffer, bufsize, position, &_npartspec[0], PARTTYPE_COUNTER, MPI_UNSIGNED);
    MyMPI_Unpack(buffer, bufsize, position, &_ndim, 1, MPI_UNSIGNED);
    MyMPI_Unpack(buffer, bufsize, position, &_box[0], ndim_ + ndim_,
                 MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_time, 1, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_periodic, 1, MPI_UNSIGNED);
    MyMPI_Unpack(buffer, bufsize, position, &_second_order, 1, MPI_UNSIGNED);
    MyMPI_Unpack(buffer, bufsize, position, &_static, 1, MPI_UNSIGNED);
    MyMPI_Unpack(buffer, bufsize, position, &_global_timestep, 1, MPI_UNSIGNED);
    MyMPI_Unpack(buffer, bufsize, position, &_gamma, 1, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_gravity, 1, MPI_UNSIGNED);
    MyMPI_Unpack(buffer, bufsize, position, &_hsoft, 1, MPI_DOUBLE);
}

/**
 * @brief Dump header information to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void Header::dump(RestartFile& rfile) {
    rfile.write(_npart);
    rfile.write(_npartspec, PARTTYPE_COUNTER);
    rfile.write(_ndim);
    rfile.write(_box, ndim_ + ndim_);
    rfile.write(_time);
    rfile.write(_periodic);
    rfile.write(_second_order);
    rfile.write(_static);
    rfile.write(_global_timestep);
    rfile.write(_gamma);
    rfile.write(_gravity);
    rfile.write(_hsoft);
}

/**
 * @brief Restart constructor. Initialize the header using the given RestartFile
 *
 * @param rfile RestartFile to read from
 */
Header::Header(RestartFile& rfile) {
    rfile.read(_npart);
    rfile.read(_npartspec, PARTTYPE_COUNTER);
    rfile.read(_ndim);
    rfile.read(_box, ndim_ + ndim_);
    rfile.read(_time);
    rfile.read(_periodic);
    rfile.read(_second_order);
    rfile.read(_static);
    rfile.read(_global_timestep);
    rfile.read(_gamma);
    rfile.read(_gravity);
    rfile.read(_hsoft);
}
