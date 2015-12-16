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
 * @file GasParticle.hpp
 *
 * @brief Gas particle: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef GASPARTICLE_HPP
#define GASPARTICLE_HPP

#include <vector>

#include "Particle.hpp"
#include "StateVector.hpp"

class VorGen;
class VorCell;
class Block;

/**
  * \brief Representation of a gas particle (in SPH style thinking)
  *
  * The particle is in fact a voronoi grid cell, its position is the position of
  * the mesh generator. For most purposes however, it is easier to think of the
  * data being contained in a particle (e.g. for visualization).
  */
class GasParticle : public Particle{
private:
    /*! \brief Gravitational acceleration from the previous timestep */
    Vec _a_grav_old;

    /*! \brief Soundspeed of the gas */
    double _csnd;

    /*! \brief Pointer to the VorGen that generates the current VorCell
     *  associated to the gas particle */
    VorGen* _vorgen;

    /*! \brief Index of the VorGen associated with this particle in the
     *  DelTess VorGen list */
    unsigned int _vorgenid;

    /*! \brief Conserved quantities of the previous timestep */
    StateVector _old_Q;

    /*! \brief Size of the voronoi cell (assuming a sphere with the same volume
     *  as the cell) */
    double _h;

    /*! \brief Vector holding an integer for every MPI process specifying which
     *  periodic/reflective copies of the particle already exist */
    std::vector<unsigned int> _copies;

    /*! \brief Booleans specifying to which MPI processes this particle was
     *  already exported */
    bool *_exports;

    /*! \brief Primitive variables (density, velocity and pressure) for the
     *  gas */
    StateVector _W;

    /*! \brief Conserved quantities (mass, momentum and energy) for the gas */
    StateVector _Q;

    /*! \brief Integrated fluxes for the conserved quantities during the current
     *  timestep */
    StateVector _dQ;

    /*! \brief Gravitational energy flux through the faces during the current
     *  timestep */
    Vec _delta_E_grav;

    /*! \brief Gradients for the primitive variables in the cell */
    StateVector _gradientvecs[ndim_];

    /*! \brief Radius used to search for neighbours to complete the VorCell
     *  associated to this particle */
    double _max_radius;

    /*! \brief Flag marking this particle as a pseudo particle */
    bool _is_pseudo;

    /*! \brief Velocity for a Springel like regularization algorithm */
    Vec _mesh_v;

    /*! \brief Mass for a Springel like regularization algorithm */
    double _mesh_m;

    /*! \brief Maximum Mach number occuring for all shockwaves in flux
     *  calculations for faces during the current timestep */
    double _max_mach;

    /*! \brief Coordinates of the geometrical centroid of the latest Voronoi
     *  cell that was associated with this particle */
    Vec _centroid;

    /*! \brief Index of the particle in the local ParticleVector */
    unsigned int _local_id;

    /*! \brief Factor for gravitational softening correction term */
    double _eta;

    /*! \brief Velocity of the particle during the previous timestep */
    Vec _old_v;

    /*! \brief Real timestep of the particle (as opposed to the integer value
     *  stored indirectly in Particle */
    double _real_dt;

    /*! \brief Total surface area of the Voronoi cell associated with this
     *  particle */
    double _total_area;

public:
    GasParticle();
    GasParticle(Vec pos);
    GasParticle(GasParticle& p);
    GasParticle(void* buffer, int bufsize, int* position);
    virtual ~GasParticle();

    /**
      * \brief Mark this Particle as a GasParticle
      *
      * @return PARTTYPE_GAS
      */
    ParticleType type(){
        return PARTTYPE_GAS;
    }

    void set_soundspeed(double csnd);
    double get_soundspeed();

    void set_vorgen(VorGen* point);
    VorGen* get_vorgen();
    VorCell* get_cell();

    void set_vorgenid(unsigned int vorgenid);
    unsigned int get_vorgenid();

    void set_W(StateVector W);
    const StateVector& get_Wvec();

    void set_Q(StateVector Q);
    StateVector& get_Qvec();

    double h();
    void set_h(double h);

    double get_total_area();
    void set_total_area(double total_area);

    void apply_gravity(double dt);

    void save_Q();

    void set_mass(double mass);

    /**
      * \brief Set the mass of the particle during the previous timestep
      *
      * @param old_mass Value to store in the conserved quantities StateVector
      */
    void set_old_mass(double old_mass){
        _old_Q[0] = old_mass;
    }

    /**
      * \brief Get the mass of the particle
      *
      * The mass is the first component of the conserved quantities StateVector
      *
      * @return The mass of the gas inside the Voronoi grid cell
      */
    inline double get_mass(){
        return _Q[0];
    }

    /**
      * \brief Get the mass of the particle during the previous timestep
      *
      * @return The mass of the gas inside the Voronoi grid cell during the
      * previous timestep
      */
    inline double get_old_mass(){
        return _old_Q[0];
    }

    void set_mesh_v(Vec& mesh_v);
    Vec& get_mesh_v();

    void set_mesh_m(double mesh_m);
    double get_mesh_m();

    void set_pseudo();
    bool is_pseudo();

    void update_W(double gamma);

    double get_gradient(unsigned int Windex, unsigned int dimindex);

    void set_gradients(StateVector* gradients);

    void get_gradients(StateVector* gradients);

    void increase_dQ(StateVector dQ);
    StateVector& get_dQvec();
    void reset_dQ();

    void increase_delta_E(Vec dE);
    void add_dE_grav_cell(double dt);
    void reset_dE_grav();
    Vec get_delta_E();

    void update_Q();

    void make_copy(unsigned int id, unsigned int index=0);
    bool is_copied(unsigned int id, unsigned int index=0);
    void reset_copies();

    void do_export(unsigned int id);
    bool is_exported(unsigned int id);
    void reset_export();

    void set_max_radius(double max_radius);
    double get_max_radius();

    void set_paq(unsigned int index, double paq);

    void save_props(Block& block);
    void save_grid(Block& block);

    void save_props_gadget(Block& block);

    /**
      * \brief Set the gravitational acceleration of this particle
      *
      * The current gravitational acceleration is saved to another variable
      * before it is overwritten.
      *
      * @param a_grav New gravitational acceleration for the particle
      */
    inline void set_gravitational_acceleration(Vec a_grav){
        _a_grav_old = _a_grav_new;
        _a_grav_new = a_grav;
    }

    void set_max_mach(double mach);
    double get_max_mach();

    virtual void pack_data(void* buffer, int bufsize, int* position);

    virtual void dump(RestartFile &rfile);
    GasParticle(RestartFile &rfile);

    void dump_ascii(std::ostream& stream);

    void set_centroid(Vec centroid);
    Vec& get_centroid();

    void set_local_id(unsigned int local_id);
    unsigned int get_local_id();

    void set_eta(double eta);
    double get_eta();

    void set_real_timestep(double real_dt);
    double get_real_timestep();
};

/**
 * @brief Exception thrown when something goes wrong during the update of the
 * conserved quantities
 */
class UpdateQException : public std::exception{
private:
    /*! \brief GasParticle throwing the exception */
    GasParticle& _p;

public:
    UpdateQException(GasParticle& p);

    virtual const char* what() const throw();
};

#endif // GASPARTICLE_HPP
