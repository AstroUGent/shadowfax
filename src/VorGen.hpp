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
 * @file VorGen.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief contains class VorGen
 *
 * A VorGen represents a point. It holds the coordinates of the point and a list
 * of triangles/tetrahedra of which it is a vertex. It also stores a pointer to
 * the VorCell of which it is the center.
 */
#ifndef HEADER_VORGEN
#define HEADER_VORGEN

#include "Vec.hpp"  // for Vec, operator-
#include <list>     // for list
#include <ostream>  // for ostream

class GasParticle;
class Simplex;
class VorCell;

/**
 * @brief Voronoi generator
 *
 * When a GasParticle is added to a VorTess, a corresponding VorGen with the
 * same coordinates is created and added instead. A VorGen is also created for
 * periodic, reflective and ghost copies of particles.
 *
 * During the Delaunay tesselation construction, the VorGen holds a list of
 * simplices of which it is part.
 */
class VorGen {
  private:
    /*! \brief Position of the point */
    Vec _p;
    /*! \brief Position of the point in the interval [1,2], used for exact
     *  geometrical tests */
    Vec _p12;
    /*! \brief VorCell of which this point is the generator */
    VorCell* _cell;
    /*! \brief GasParticle associated with this generator */
    GasParticle* _particle;
    /*! \brief Index of the particle in the local ParticleVector */
    unsigned int _particle_index;
    /*! \brief Simplices that have this point as one of their vertices */
    std::list<Simplex*> _tetrahedra;
    /*! \brief Search radius for Delaunay tesselation completion */
    double _sr;
    /*! \brief Variable that is apparently no longer used... */
    VorGen* _mirror;
    /*! \brief Indices of the simplices that have this point as one of their
     *  vertices in the DelTess simplex list */
    std::list<unsigned int> _simplices;
    /*! \brief Integer keeping track of which reflective copies of this point
     *  exist */
    unsigned int _mirrored;
    /*! \brief ID of the GasParticle corresponding to this generator */
    unsigned long _original;
    /*! \brief Variable holding the process number where this point resides (for
     *  ghost points) */
    unsigned int _process;

    /*! \brief Flag used to flag this VorGen */
    bool _flag;

  public:
#if ndim_ == 3
    VorGen(double x, double y, double z);
#else
    VorGen(double x, double y);
#endif
    ~VorGen();
    void add_tetrahedron(Simplex* tetrahedron);
    void add_simplex(unsigned int index);
    void remove_tetrahedron(Simplex* tetrahedron);
    void remove_simplex(unsigned int index);

    /**
     * @brief Calculate the distance between this point and the given point
     *
     * @param point Second point
     * @return Distance between this point and the given point
     */
    inline double distance(VorGen point) {
        return (_p - point._p).norm();
    }

    /**
     * @brief Get the x-component of the position of this point
     *
     * @return x-coordinate of the point
     */
    inline double x() {
        return _p[0];
    }

    /**
     * @brief Get the y-component of the position of this point
     *
     * @return y-coordinate of the point
     */
    inline double y() {
        return _p[1];
    }

#if ndim_ == 3
    /**
     * @brief Get the z-component of the position of this point
     *
     * @return z-coordinate of the point
     */
    inline double z() {
        return _p[2];
    }
#endif

    /**
     * @brief Get the component with the given index of the position of this
     * point
     *
     * @param index Index of a coordinate of this point
     * @return Value of the requested coordinate
     */
    inline double pos(int index) {
        return _p[index];
    }

    void move(double x, double y, double z);
    void set_cell(VorCell* cell);
    VorCell* get_cell();
    void set_mirror(VorGen* mirror);
    VorGen* get_mirror();
    void set_particle(GasParticle* particle);
    GasParticle* get_particle();
    void set_particle_id(unsigned int index);
    unsigned int get_particle_id();
    void print(std::ostream& stream);
    std::list<Simplex*> get_tetrahedra();
    std::list<unsigned int> get_simplices();
    void set_search_radius(double sr);
    double get_search_radius();
    void mirror(unsigned int casus);
    bool mirrored(unsigned int casus);

    Vec get_position();

    void set_id(unsigned int id);
    unsigned int get_id();
    void reset_simplices();

    void set_original(unsigned long id);
    unsigned long get_original();

    void set_process(unsigned int process);
    unsigned int get_process();

    void set_p12(Vec p12);
    Vec& get_p12();

    void flag();
    void unflag();
    bool flagged();
};

#endif
