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
 * @file Particle.cpp
 *
 * @brief Common properties of all particles: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Particle.hpp"
#include "RestartFile.hpp"
#include <iostream>  // for ostream
using namespace std;

/**
 * \brief Construct a default empty Particle with id, position, velocity... set
 * to zero
 */
Particle::Particle() : Hilbert_Object() {
    _old_a = 0.;
    _id = 0;
    _starttime = 0;
    _endtime = 0;
    _hsoft = 0.;
    _epot = 0.;
    reset_comp_cost();
}

/**
 * \brief Construct a Particle with the given position and all other properties
 * set to zero
 *
 * @param pos Vec containing the position of the particle
 */
Particle::Particle(Vec pos) : Hilbert_Object() {
    _x = pos;
    _old_a = 0.;
    _id = 0;
    _starttime = 0;
    _endtime = 0;
    _hsoft = 0.;
    _epot = 0.;
    reset_comp_cost();
}

/**
 * \brief Get the x-coordinate of the position of the particle
 *
 * @return The x-coordinate of the position of the particle
 */
double Particle::x() {
    return _x[0];
}

/**
 * \brief Get the y-coordinate of the position of the particle
 *
 * @return The y-coordinate of the position of the particle
 */
double Particle::y() {
    return _x[1];
}

/**
 * \brief Get the z-coordinate of the position of the particle
 *
 * This function only works properly when the code is compiled for 3D usage,
 * otherwise, it just returns 0.
 *
 * @return The z-coordinate of the position of the particle
 */
double Particle::z() {
#if ndim_ == 3
    return _x[2];
#else
    return 0.;
#endif
}

/**
 * \brief Get the specified coordinate of the position of the particle
 *
 * Possible parameter values are dependent on the number of dimensions. Usage of
 * any other values will compile, but will also lead to strange bugs in your
 * program.
 *
 * @param index The index of the coordinate you want to access (0 = x, 1 = y,
 * 2 = z)
 * @return The coordinate of the particle position corresponding to the given
 * index
 */
double Particle::pos(int index) {
    return _x[index];
}

/**
  * \brief Get the velocity of the particle
  *
  * @return Vec containing the velocity of the particle
  */
Vec Particle::get_velocity() {
    return _v;
}

/**
 * \brief Get the x-component of the velocity of the particle
 *
 * @return The x-component of the velocity of the particle
 */
double Particle::vx() {
    return _v[0];
}

/**
 * \brief Get the y-component of the velocity of the particle
 *
 * @return The y-component of the velocity of the particle
 */
double Particle::vy() {
    return _v[1];
}

/**
 * \brief Get the z-component of the velocity of the particle
 *
 * This function only works properly when the code is compiled for 3D usage,
 * otherwise, it just returns 0.
 *
 * @return The z-component of the velocity of the particle
 */
double Particle::vz() {
#if ndim_ == 3
    return _v[2];
#else
    return 0.;
#endif
}

/**
 * \brief Get the specified coordinate of the velocity of the particle
 *
 * Possible parameter values are dependent on the number of dimensions. Usage of
 * any other values will compile, but will also lead to strange bugs in your
 * program.
 *
 * @param index The index of the coordinate you want to access (0 = x, 1 = y,
 * 2 = z)
 * @return The coordinate of the particle velocity corresponding to the given
 * index
 */
double Particle::vel(int index) {
    return _v[index];
}

/**
 * \brief Set the x-coordinate of the position of the particle
 *
 * @param x The new x-coordinate for the position of the particle
 */
void Particle::set_x(double x) {
    _x[0] = x;
}

/**
 * \brief Set the y-coordinate of the position of the particle
 *
 * @param y The new y-coordinate for the position of the particle
 */
void Particle::set_y(double y) {
    _x[1] = y;
}

/**
 * \brief Set the z-coordinate of the position of the particle
 *
 * This function only works properly when the code is compiled for 3D usage.
 * Other uses will not result in compilation errors, but will lead to strange
 * results.
 *
 * @param z The new z-coordinate for the position of the particle
 */
void Particle::set_z(double z) {
    _x[2] = z;
}

/**
 * \brief Print the Particle to the given ostream
 *
 * This method currently does nothing
 *
 * @param stream The ostream to write to
 */
void Particle::print(ostream& stream) {
    // do nothing
}

/**
  * \brief Print the position and id of the particle to a line in the given
  * stream
  *
  * The coordinates of the position are separated by tabs, there are always 3
  * coordinates outputted, even in the 2D case where the third coordinate will
  * be undefined.
  *
  * @param stream std::ostream to write to
  */
void Particle::print_gen(ostream& stream) {
    stream << _x[0] << "\t" << _x[1] << "\t" << _x[2] << "\t" << _id << "\n";
}

/**
  * \brief Set the velocity of the particle
  *
  * We always require 3 coordinates (even in the 2D case), but the third
  * coordinate is only used in the 3D case.
  *
  * @param vx The new x-coordinate of the velocity of the particle
  * @param vy The new y-coordinate of the velocity of the particle
  * @param vz The new z-coordinate of the velocity of the particle
  */
void Particle::set_v(double vx, double vy, double vz) {
    _v[0] = vx;
    _v[1] = vy;
#if ndim_ == 3
    _v[2] = vz;
#endif
}

/**
 * @brief Set the velocity of the particle to the given value
 *
 * @param v New velocity for the particle
 */
void Particle::set_velocity(Vec& v) {
    _v[0] = v[0];
    _v[1] = v[1];
#if ndim_ == 3
    _v[2] = v[2];
#endif
}

/**
 * \brief Get the unique identifier for the Particle
 *
 * @return The id of the Particle
 */
unsigned long Particle::id() {
    return _id;
}

/**
 * Set the unique identifier for the Particle
 * @param id The new ID for the Particle
 */
void Particle::set_id(unsigned long id) {
    _id = id;
}

/**
 * @brief Drift the particle for the given time with a constant velocity equal
 * to the particle velocity
 *
 * @param dt double precision floating-point time interval during which to drift
 */
void Particle::drift(double dt) {
    _x += dt * _v;
}

/**
 * @brief Get the integer timestep for this particle
 *
 * @return unsigned long integer timestep for this particle on the integer
 * simulation timeline
 */
unsigned long Particle::get_timestep() {
    return _endtime - _starttime;
}

/**
 * @brief Set the integer timestep for this particle to the given value
 *
 * The location of the particle on the integer simulation timeline is also
 * updated by setting the starting time of the particle to the previous end
 * time.
 *
 * @param timestep unsigned long integer timestep for this particle
 */
void Particle::set_timestep(unsigned long timestep) {
    _starttime = _endtime;
    _endtime = _starttime + timestep;
}

/**
 * @brief Reset the particle timestep to the given value
 *
 * This method is almost the same as Particle::set_timestep(), but it does not
 * set the particle starting time to the previous end time.
 *
 * @param timestep unsigned long integer timestep for this particle
 */
void Particle::reset_timestep(unsigned long timestep) {
    _endtime = _starttime + timestep;
}

/**
 * @brief Get the end time of the particle on the integer simulation timeline
 *
 * @return The end time of the current particle timestep on the integer
 * simulation timeline
 */
unsigned long Particle::get_endtime() {
    return _endtime;
}

/**
 * @brief Set the end time of the current particle timestep on the integer
 * simulation timeline to the given value
 *
 * @param endtime unsigned long integer end time for this particle
 */
void Particle::set_endtime(unsigned long endtime) {
    _endtime = endtime;
}

/**
 * @brief Set the starting time of the current particle timestep on the integer
 * simulation timeline to the given value
 *
 * @param integertime unsigned long integer starting time for this particle
 */
void Particle::set_starttime(unsigned long integertime) {
    _starttime = integertime;
    _endtime = integertime;
}

/**
 * @brief Get the starting time of the particle on the integer simulation
 * timeline
 *
 * @return The starting time of the current particle timestep on the integer
 * simulation timeline
 */
unsigned long Particle::get_starttime() {
    return _starttime;
}

/**
 * @brief Get the gravitational acceleration for this particle
 *
 * @return The gravitational acceleration for this particle
 */
Vec Particle::get_gravitational_acceleration() {
    return _a_grav_new;
}

/**
 * @brief Set the magnitude of the gravitational acceleration during the
 * previous timestep
 *
 * This value is used for the relative opening criterion during the
 * gravitational treewalk.
 *
 * @param a_grav The magnitude of the gravitational acceleration during the
 * previous timestep this particle was active
 */
void Particle::set_old_acceleration(double a_grav) {
    _old_a = a_grav;
}

/**
 * @brief Add the given velocity to the particle velocity
 *
 * @param dv Acceleration multiplied with a real time interval
 */
void Particle::accelerate(Vec dv) {
    _v += dv;
}

/**
 * @brief Identical to Particle::drift()
 *
 * @param dt double precision floating-point time interval during which to drift
 */
void Particle::move(double dt) {
    _x += _v * dt;
}

/**
 * @brief Set the gravitational potential for this particle
 *
 * @param epot New value for the gravitational potential
 */
void Particle::set_gravitational_potential(double epot) {
    _epot = epot;
}

/**
 * @brief Get the gravitational potential for this particle
 *
 * @return Gravitational potential of this particle
 */
double Particle::get_gravitational_potential() {
    return _epot;
}

/**
 * @brief Set the gravitational softening length for this particle
 *
 * @param hsoft New value for the gravitational softening length
 */
void Particle::set_hsoft(double hsoft) {
    _hsoft = hsoft;
}

/**
 * @brief Get the gravitational softening length for this particle
 *
 * @return The gravitational softening length of this particle
 */
double Particle::get_hsoft() {
    return _hsoft;
}

/**
 * @brief MPI constructor. Initialize the particle using the given MPI buffer
 *
 * @param buffer Buffer to read from
 * @param bufsize Size of the buffer
 * @param position Position of the buffer (is updated)
 */
Particle::Particle(void* buffer, int bufsize, int* position)
        : Hilbert_Object(buffer, bufsize, position) {
    MyMPI_Unpack(buffer, bufsize, position, &_x[0], ndim_, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_v[0], ndim_, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_a_grav_new[0], ndim_, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_old_a, 1, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_id, 1, MPI_UNSIGNED_LONG);
    MyMPI_Unpack(buffer, bufsize, position, &_starttime, 1, MPI_UNSIGNED_LONG);
    MyMPI_Unpack(buffer, bufsize, position, &_endtime, 1, MPI_UNSIGNED_LONG);
    MyMPI_Unpack(buffer, bufsize, position, &_comp_cost, 1, MPI_UNSIGNED);
    MyMPI_Unpack(buffer, bufsize, position, &_hsoft, 1, MPI_DOUBLE);
    _epot = 0.;
}

/**
 * @brief Write particle data to the given MPI buffer for communication
 *
 * @param buffer Buffer to write to
 * @param bufsize Size of the buffer
 * @param position Position of the buffer (is updated)
 */
void Particle::pack_data(void* buffer, int bufsize, int* position) {
    Hilbert_Object::pack_data(buffer, bufsize, position);
    MyMPI_Pack(&_x[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_v[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_a_grav_new[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_old_a, 1, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_id, 1, MPI_UNSIGNED_LONG, buffer, bufsize, position);
    MyMPI_Pack(&_starttime, 1, MPI_UNSIGNED_LONG, buffer, bufsize, position);
    MyMPI_Pack(&_endtime, 1, MPI_UNSIGNED_LONG, buffer, bufsize, position);
    MyMPI_Pack(&_comp_cost, 1, MPI_UNSIGNED, buffer, bufsize, position);
    MyMPI_Pack(&_hsoft, 1, MPI_DOUBLE, buffer, bufsize, position);
}

/**
 * @brief Dump the particle to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void Particle::dump(RestartFile& rfile) {
    Hilbert_Object::dump(rfile);
    rfile.write(_x);
    rfile.write(_v);
    rfile.write(_a_grav_new);
    rfile.write(_old_a);
    rfile.write(_id);
    rfile.write(_starttime);
    rfile.write(_endtime);
    rfile.write(_comp_cost);
    rfile.write(_epot);
    rfile.write(_hsoft);
}

/**
 * @brief Restart constructor. Initialize the particle using the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 */
Particle::Particle(RestartFile& rfile) : Hilbert_Object(rfile) {
    rfile.read(_x);
    rfile.read(_v);
    rfile.read(_a_grav_new);
    rfile.read(_old_a);
    rfile.read(_id);
    rfile.read(_starttime);
    rfile.read(_endtime);
    rfile.read(_comp_cost);
    rfile.read(_epot);
    rfile.read(_hsoft);
}
