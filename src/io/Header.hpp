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
 * @file Header.hpp
 *
 * @brief Snapshot header: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEADER_HPP
#define HEADER_HPP

#include "utilities/ParticleTypes.hpp"
#include <ostream>

class RestartFile;

/**
  * \brief Header information about a ParticleVector
  *
  * A Header contains important information about a ParticleVector that can be
  * directly written to a snapshot, like the total number of particle, the
  * number of dimensions and all sorts of integration options.
  *
  * It also contains methods to check if header information is compatible with
  * the current version of the code.
  */
class Header {
  private:
    /*! \brief The total number of particles in the simulation */
    unsigned int _npart;

    /*! \brief The total number of gas and dm particles in the simulation */
    unsigned int _npartspec[PARTTYPE_COUNTER];

    /*! \brief The number of dimensions of the simulation (2 or 3) */
    unsigned int _ndim;

    /*! \brief Origin and side of the simulation box */
    double _box[ndim_ + ndim_];

    /*! \brief The current time of the simulation */
    double _time;

    /*! \brief Boolean specifying if the code is run with periodic or reflective
     *  boundaries (0 = reflective, 1 = periodic) */
    unsigned int _periodic;

    /*! \brief Boolean specifying if the code is run with a second or first
     *  order integration scheme (0 = first order, 1 = second order) */
    unsigned int _second_order;

    /*! \brief Boolean specifying if the code is run with a static or a moving
     *  mesh (0 = static, 1 = moving) */
    unsigned int _static;

    /*! \brief Boolean specifying if a global timestep is used (0 = individual
     *  timesteps, 1 = global timestep) */
    unsigned int _global_timestep;

    /*! \brief Value used for the adiabatic index of the gas (not sure if this
     *  works correctly) */
    double _gamma;

    /*! \brief Boolean specifying if the code is run with gravity or not (0 = no
     *  gravity, 1 = gravity) */
    unsigned int _gravity;

    /*! \brief Value of the softening length for collisionless particles */
    double _hsoft;

  public:
    Header();
    ~Header() {}

    void set_ngaspart(unsigned int npart);
    void set_ndmpart(unsigned int npart);
    void set_nstarpart(unsigned int npart);
    void set_box(double* box);
    void set_time(double time);
    void set_global_timestep(bool global_timestep);
    void set_periodic(bool periodic);
    void set_gravity(bool gravity);
    void set_hsoft(double hsoft);

    void* get_npartspec();
    void* get_ndim();
    void* get_box();
    void* get_time();
    void* get_periodic();
    void* get_second_order();
    void* get_static();
    void* get_global_timestep();
    void* get_gamma();
    void* get_gravity();
    void* get_hsoft();

    unsigned int npart();
    unsigned int ngaspart();
    unsigned int ndmpart();
    unsigned int nstarpart();
    void box(double* box);
    double time();
    bool global_timestep();
    bool periodic();
    bool gravity();
    double hsoft();

    void check_makeflags();

    void pack_data(void* buffer, int bufsize, int* position);
    Header(void* buffer, int bufsize, int* position);

    void dump(RestartFile& rfile);
    Header(RestartFile& rfile);
};

#endif  // HEADER_HPP
