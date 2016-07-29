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
 * @file ParticleVector.hpp
 *
 * @brief Custom Particle vector: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PARTICLEVECTOR_HPP
#define PARTICLEVECTOR_HPP

#include "../src/io/Header.hpp"
#include "DelCont.hpp"
#include "ParticleTypes.hpp"
#include "Tree.hpp"
#include "utilities/Timer.hpp"
#include <vector>

class RectangularBox;
class GasParticle;
class DMParticle;
class ParameterFile;
class ParticleConverter;
class StarParticle;

#define PARTICLEVECTOR_DEFAULT_PERIODICFLAG false
#define PARTICLEVECTOR_DEFAULT_EWALDFLAG false
#define PARTICLEVECTOR_DEFAULT_ALPHA 2.
#define PARTICLEVECTOR_DEFAULT_SIZE 64

/**
 * @brief Specialized container to store the particles of the simulation
 *
 * The container contains an internal vector storing the particles, and utility
 * arrays for handling the different particle types: an array with the number
 * of particles of each type (locally) and the total number of particles, and an
 * array with the offsets of the different particle types in the internal
 * vector. We make sure that the particles are always grouped together per type,
 * in the order imposed by the order in ParticleType.
 *
 * The container also contains a Tree structure for the particles, a
 * RectangularBox that encompasses them and a Header that can be used to write
 * the particles to a snapshot.
 *
 * It implements methods to sort the particles (in parallel) and to build and
 * update the Tree.
 */
class ParticleVector {
  private:
    /*! \brief List of particles */
    std::vector<Particle*> _particles;

    /*! \brief Offsets of the different particle types in the internal particle
     *  vector */
    unsigned int _offsets[PARTTYPE_COUNTER];

    /*! \brief Sizes of the different particle types on the local process */
    unsigned int _sizes[PARTTYPE_COUNTER + 1];

    /*! \brief The number of active particles at the current system time */
    unsigned int _numactive;

    /*! \brief Octtree used for neighbour searches and tree walks */
    Tree _tree;

    /*! \brief RectangularBox specifying the dimensions of the enclosing box
     *  that contains all particles */
    RectangularBox _container;

    /*! \brief Header containing global information on all particles over all
     *  MPI processes, used to write Header information to snapshot files */
    Header _header;

    /*! \brief Timer used to quantify the time spent in sorting */
    Timer _sorttimer;
    /*! \brief Timer used to quantify the time spent in tree building */
    Timer _treetimer;

    /**
     * @brief Comparison function used to sort particles on type
     *
     * @param a Particle A
     * @param b Particle B
     * @return True if the ParticleType of B is larger than the ParticleType of
     * A
     */
    static bool compare_type(Particle* a, Particle* b) {
        return a->type() < b->type();
    }

  public:
    ParticleVector(bool global_timestep, RectangularBox container,
                   bool periodic = PARTICLEVECTOR_DEFAULT_PERIODICFLAG,
                   bool do_ewald = PARTICLEVECTOR_DEFAULT_EWALDFLAG,
                   double alpha = PARTICLEVECTOR_DEFAULT_ALPHA,
                   unsigned int size = PARTICLEVECTOR_DEFAULT_SIZE);
    ParticleVector(ParameterFile* parameters, RectangularBox container,
                   bool periodic = PARTICLEVECTOR_DEFAULT_PERIODICFLAG,
                   bool do_ewald = PARTICLEVECTOR_DEFAULT_EWALDFLAG);
    ~ParticleVector();

    void add_gas_particle(GasParticle* particle);
    void add_DM_particle(DMParticle* particle);
    void add_star_particle(StarParticle* particle);
    void construct_tree();
    void sort();

    void set_container(RectangularBox container);
    void set_periodic(bool periodic);

    void finalize();

    /**
     * @brief Get the size of the gas particle list
     *
     * @return Size of the local gas particle list
     */
    unsigned int gassize() {
        return _sizes[PARTTYPE_GAS];
    }

    /**
     * @brief Get the size of the dark matter particle list
     *
     * @return Size of the local dark matter particle list
     */
    unsigned int dmsize() {
        return _sizes[PARTTYPE_DM];
    }

    /**
     * @brief Get the size of the star particle list
     *
     * @return Size of the local star particle list
     */
    unsigned int starsize() {
        return _sizes[PARTTYPE_STAR];
    }

    /**
     * @brief Resize the gas particle list
     *
     * @param size New size for the gas particle list
     */
    void resizegas(unsigned int size) {
        _sizes[PARTTYPE_GAS] = size;
        _sizes[PARTTYPE_COUNTER] = 0;
        for(int type = 0; type < PARTTYPE_COUNTER; type++) {
            _sizes[PARTTYPE_COUNTER] += _sizes[type];
        }
        _particles.resize(_sizes[PARTTYPE_COUNTER]);
        for(int type = PARTTYPE_GAS + 1; type < PARTTYPE_COUNTER; type++) {
            _offsets[type] = _offsets[type - 1] + _sizes[type - 1];
        }
    }

    /**
     * @brief Resize the dark matter particle list
     *
     * @param size New size for the dark matter particle list
     */
    void resizedm(unsigned int size) {
        _sizes[PARTTYPE_DM] = size;
        _sizes[PARTTYPE_COUNTER] = 0;
        for(int type = 0; type < PARTTYPE_COUNTER; type++) {
            _sizes[PARTTYPE_COUNTER] += _sizes[type];
        }
        _particles.resize(_sizes[PARTTYPE_COUNTER]);
        for(int type = PARTTYPE_DM + 1; type < PARTTYPE_COUNTER; type++) {
            _offsets[type] = _offsets[type - 1] + _sizes[type - 1];
        }
    }

    /**
     * @brief Resize the star particle list
     *
     * @param size New size for the star particle list
     */
    void resizestar(unsigned int size) {
        _sizes[PARTTYPE_STAR] = size;
        _sizes[PARTTYPE_COUNTER] = 0;
        for(int type = 0; type < PARTTYPE_COUNTER; type++) {
            _sizes[PARTTYPE_COUNTER] += _sizes[type];
        }
        _particles.resize(_sizes[PARTTYPE_COUNTER]);
        for(int type = PARTTYPE_STAR + 1; type < PARTTYPE_COUNTER; type++) {
            _offsets[type] = _offsets[type - 1] + _sizes[type - 1];
        }
    }

    /**
     * @brief Access the ith element of the gas particle list
     *
     * @param i unsigned integer index of the requested gas particle in the list
     * @return Reference to the ith gas particle in the list
     */
    GasParticle*& gas(unsigned int i) {
        return (GasParticle*&)_particles[_offsets[PARTTYPE_GAS] + i];
    }

    /**
     * @brief Access the ith element of the dark matter particle list
     *
     * @param i unsigned integer index of the requested dark matter particle in
     * the list
     * @return Reference to the ith dark matter particle in the list
     */
    DMParticle*& dm(unsigned int i) {
        return (DMParticle*&)_particles[_offsets[PARTTYPE_DM] + i];
    }

    /**
     * @brief Access the ith element of the star particle list
     *
     * @param i unsigned integer index of the requested star particle in the
     * list
     * @return Reference to the ith star particle in the list
     */
    StarParticle*& star(unsigned int i) {
        return (StarParticle*&)_particles[_offsets[PARTTYPE_STAR] + i];
    }

    /**
     * @brief Get a reference to the particle octtree
     *
     * @return Refernce to the Tree
     */
    Tree& get_tree() {
        return _tree;
    }

    /**
     * @brief Get a reference to the last gas particle in the list
     *
     * @return Reference to the last gas particle in the list
     */
    GasParticle*& gasback() {
        return (GasParticle*&)
                _particles[_offsets[PARTTYPE_GAS] + _sizes[PARTTYPE_GAS] - 1];
    }

    /**
     * @brief Get a reference to the last dark matter particle in the list
     *
     * @return Reference to the last dark matter particle in the list
     */
    DMParticle*& dmback() {
        return (DMParticle*&)
                _particles[_offsets[PARTTYPE_DM] + _sizes[PARTTYPE_DM] - 1];
    }

    /**
     * @brief Get a reference to the last dark matter particle in the list
     *
     * @return Reference to the last dark matter particle in the list
     */
    StarParticle*& starback() {
        return (StarParticle*&)
                _particles[_offsets[PARTTYPE_STAR] + _sizes[PARTTYPE_STAR] - 1];
    }

    /**
     * @brief Get a reference to the container containing all particles
     *
     * @return DelCont used for the Voronoi grid construction
     */
    DelCont& get_container() {
        return _container;
    }

    /**
     * @brief Check if we use a global timestep or individual timesteps
     *
     * @return True if we use a global timestep, false otherwise
     */
    bool global_timestep() {
        return _header.global_timestep();
    }

    void print_local_particles(std::string filename);

    Header& get_header();

    /**
     * @brief Get a reference to the Header
     *
     * @return Reference to the Header
     */
    Header& get_local_header() {
        return _header;
    }

    void set_numactive(unsigned int numactive);
    unsigned int get_numactive();

    void convert(ParticleConverter& converter, unsigned long current_time);

    void dump(RestartFile& rfile);
    ParticleVector(RestartFile& rfile, RectangularBox& box,
                   bool periodic = PARTICLEVECTOR_DEFAULT_PERIODICFLAG,
                   bool do_ewald = PARTICLEVECTOR_DEFAULT_EWALDFLAG,
                   double alpha = PARTICLEVECTOR_DEFAULT_ALPHA,
                   unsigned int size = PARTICLEVECTOR_DEFAULT_SIZE);
    ParticleVector(RestartFile& rfile, ParameterFile* parameters,
                   RectangularBox& box,
                   bool periodic = PARTICLEVECTOR_DEFAULT_PERIODICFLAG,
                   bool do_ewald = PARTICLEVECTOR_DEFAULT_EWALDFLAG);
};

#endif  // PARTICLEVECTOR_HPP
