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
 * @file Particle.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief Common properties of all particles: header
 *
 * contains class Particle
 *
 * A Particle stores information on a particle: its coordinates, mass, key and
 * the corresponding point in the Delaunay/VoronoiTesselation.
 */
#ifndef HEAD_PARTICLE
#define HEAD_PARTICLE

#include "Hilbert.hpp"        // for Hilbert_Object
#include "ParticleTypes.hpp"  // for ParticleType
#include "Vec.hpp"
#include <iostream>  // for ostream

class RestartFile;

/**
 * @brief Representation of a particle
 *
 * A Particle has properties like a position, a velocity and an acceleration and
 * offers an interface to some mass variable.
 * It also holds a timestep for a numerical integration scheme.
 */
class Particle : public Hilbert_Object {

  protected:
    /*! @brief Position of the particle */
    Vec _x;

    /*! @brief Velocity of the particle */
    Vec _v;

    /*! @brief Gravitational acceleration of the particle */
    Vec _a_grav_new;

    /*! @brief Norm of the gravitational acceleration during the previous
     *  integration timestep */
    double _old_a;

    /*! @brief Unique identifier for this Particle */
    unsigned long _id;

    /*! @brief Time on the integer timeline of the beginning of the current
     *   timestep */
    unsigned long _starttime;

    /*! @brief Time on the integer timeline of the end of the current
     *  timestep */
    unsigned long _endtime;

    /*! @brief Computational cost of the gravitational treewalk during the
     *  previous timestep, used for load-balancing */
    unsigned int _comp_cost;

    /*! @brief Gravitational potential energy of the particle */
    double _epot;

    /*! @brief Gravitational softening length of the particle */
    double _hsoft;

  public:
    Particle();
    Particle(Vec pos);
    Particle(void* buffer, int bufsize, int* position);
    virtual ~Particle() {}

    /**
      * @brief Convenience function to distinguish between different
      * implementations of Particle
      */
    virtual ParticleType type() = 0;

    double x();
    double y();
    double z();
    double vx();
    double vy();
    double vz();
    double pos(int index);

    /**
      * @brief Get the position of the particle
      *
      * @return Vec containing the position of the particle
      */
    Vec& get_position() {
        return _x;
    }

    Vec get_velocity();
    double vel(int index);
    void set_x(double x);
    void set_y(double y);
    void set_z(double z);
    void set_v(double vx, double vy, double vz);
    void set_velocity(Vec& v);

    void drift(double dt);

    void print(std::ostream& stream);
    void print_gen(std::ostream& stream);

    unsigned long id();
    void set_id(unsigned long id);

    unsigned long get_timestep();
    void set_timestep(unsigned long timestep);
    void reset_timestep(unsigned long timestep);

    unsigned long get_endtime();
    unsigned long get_starttime();

    void set_endtime(unsigned long endtime);
    void set_starttime(unsigned long integertime);

    /**
      * @brief Set the gravitational acceleration of the particle
      *
      * @param a_grav New gravitational acceleration for the particle
      */
    virtual inline void set_gravitational_acceleration(Vec a_grav) {
        _a_grav_new = a_grav;
    }

    Vec get_gravitational_acceleration();

    void set_old_acceleration(double a_grav);

    /**
      * @brief Get the norm of the gravitational acceleration during the
      * previous integration timestep
      *
      * @return The norm of the gravitational acceleration during the previous
      * integration timestep
      */
    inline double get_old_acceleration() {
        return _old_a;
    }

    /**
     * @brief Add the given number of force calculations to the computational
     * cost for this particle
     *
     * @param comp_cost unsigned integer number of force calculations
     */
    inline void add_comp_cost(unsigned int comp_cost) {
        _comp_cost += comp_cost;
    }

    /**
     * @brief Reset the computational cost for the particle to zero
     */
    inline void reset_comp_cost() {
        _comp_cost = 0;
    }

    /**
     * @brief Get the computational cost of the gravity calculation for this
     * particle during the previous timestep
     *
     * @return The computational cost for this particle
     */
    inline unsigned int get_comp_cost() {
        return _comp_cost;
    }

    /**
     * @brief Set the mass of the particle
     *
     * @param mass New mass of the particle
     */
    virtual void set_mass(double mass) = 0;

    /**
     * @brief Get the mass of the particle
     *
     * @return Mass of the particle
     */
    virtual double get_mass() = 0;

    void accelerate(Vec dv);
    void move(double dt);

    void set_gravitational_potential(double epot);
    double get_gravitational_potential();

    void set_hsoft(double hsoft);
    double get_hsoft();

    virtual void pack_data(void* buffer, int bufsize, int* position);

    virtual void dump(RestartFile& rfile);
    Particle(RestartFile& rfile);
};

#endif
