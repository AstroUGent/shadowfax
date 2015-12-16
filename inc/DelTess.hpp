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
 * \file DelTess.hpp
 *
 * @author Bert Vandenbroucke
 *
 * @brief Representation of a Delaunay tesselation
 *
 * contains class DelaunayTesselation
 *
 * The DelaunayTesselation builds and stores a DelaunayTesselation for an area/
 * volume specified by the given Container. Particles are added one by one by
 * calls to the method add_particle. For some problems, this class needs an
 * external octree, which is provided by the Tree argument.
 */

#ifndef HEAD_DELTESS
#define HEAD_DELTESS

#include <vector>
#include <ostream>
#include <string>
#include <set>
#include <map>
#include <cmath>

#include "utilities/Cuboid.hpp"

class VorTess;
class DelCont;
class GasParticle;
class VorGen;
class Tree;
class Simplex;
class VorCell;

/**
 * \brief The Delaunay triangulation
 *
 * Base class to construct a Delaunay triangulation, by adding points one by
 * one.
 * The triangulation itself is stored in Point instances by linking these
 * through Tetrahedron instances.
 */
class DelTess{
private:
    /*! \brief Tolerance value used to discriminate between exact and inexact
     *   arithmetics */
    double _tolerance;

    /*! \brief vector containing the points of the tesselation */
    std::vector<VorGen*> _points;

    /*! \brief vector containg the mirror points of the tesselation */
    std::vector<VorGen*> _mirrors;

    /*! \brief DelCont specifying the volume that contains the points of the
     *  tesselation and the boundaries */
    DelCont* _container;

    /*! \brief VorTess associated with this DelTess */
    VorTess* _voronoi;

    /*! \brief Stack of Simplex pointers that have to be deleted after
     *  restoration of Delaunayhood */
    std::vector<Simplex*> _remove_stack;

    /*! \brief Index of the most recently pushed back Simplex. Used as an
     *  independent measure of the number of simplices */
    unsigned int _lastindex;

    /*! \brief Index of the latest index that was checked for Delaunayhood. Used
     *  as a starting point in the simplex search during the next point
     *  insertion */
    unsigned int _lastchecked;

    /*! \brief Index of the last VorGen that was added to the tesselation */
    unsigned int _lastvorgen;

    /*! \brief A vector containing the simplices that make up the tesselation */
    std::vector<Simplex*> _simplices;

    /*! \brief Cuboid specifying a box that contains all points that could
     *  possibly be added to the tesselation, used to rescale point coordinates
     *  to the range [1,2] for exact geometrical tests */
    Cuboid _box12;

#if ndim_==3
    /** \brief An array holding free entries in the simplices vector.
     *
     * Used in 3D where not every Simplex that is removed is replaced by a new
     * one, so that empty entries may occur. The maximal number of entries in
     * the array is roughly determined by the maximal number of 3 to 2 flips
     * that occurs during a point insertion (since only this flip creates empty
     * entries). For some cases, this can get quite large. The value of 1000
     * seemed to work until now, but for very large cartesian grids, more could
     * be necessary...
     */
    unsigned int _free_array[1000];

    /*! \brief The number of valid elements in the _free_array. The latest valid
     *  element has index _free_size-1. */
    unsigned int _free_size;
#endif

    /*! \brief A vector of vectors of ghost particles (one vector for every
     *  process) */
    std::vector< std::vector<GasParticle*> > _ghosts;

    /*! \brief A vector of vectors of exported particles (one vector for every
     *  process) */
    std::vector< std::vector<GasParticle*> > _exports;

    /*! \brief A vector of vectors of exported particle copies (one vector for
     *  every process) */
    std::vector< std::vector<GasParticle*> > _exportcopies;

    /*! \brief A vector of ghost cells (cells associated with particles on
     *  another process) */
    std::vector<VorCell*> _ghostcells;

    /*! \brief Flag indicating whether the simulation box is periodic (true) or
     *  reflective (false) */
    bool _periodic;

    void add_ghosts();
#if ndim_==3
    void two_to_six_flip(unsigned int index, unsigned int* indices);
    void four_to_four_flip(unsigned int v1, unsigned int v2, unsigned int v3,
                           unsigned int v4, unsigned int s1, unsigned int s2,
                           unsigned int s3, unsigned int s4);
    void n_to_2n_flip(unsigned int index, unsigned int* simplex,
                      unsigned int n);
    void get_common_tetrahedron(unsigned int point1, unsigned int point2,
                                unsigned int point3, unsigned int* sids,
                                unsigned int* common_tetrahedron);
    double normcross(VorGen* a, VorGen* b, VorGen* point);
#endif
    bool check_simplex(unsigned int index, unsigned int vindex);
    unsigned int find_simplex(VorGen* point, unsigned int* found);
    void remove_tetrahedron(Simplex* tetrahedron);
    void clear_remove_stack();

    Vec get_p12(Vec pos);

public:
    DelTess(DelCont* container, unsigned int numpart, bool periodic = false,
            double tolerance = 1.e-9);
    ~DelTess();

    void add_particle(GasParticle* particle, unsigned int index);
    void add_point(unsigned int index);
    void add_voronoi_tesselation(VorTess* voronoi);
    void output_tesselation(std::ostream &stream);
    DelCont* get_container();
    std::vector<VorGen*> get_points();
    void add_mirrors(Tree& parttree);
    void check_tesselation();
    Simplex* get_simplex(unsigned int index);
    void set_relations();
    unsigned int get_size();

    void check_methods();

    void update_Ws();
    void update_gradients();
    void update_dQs();
    void update_dts(unsigned long currentTime);
    void update_gravitational_corrections();

    void get_triangles(std::vector<float>& positions,
                       std::vector<int>& connectivity);
};

#endif
