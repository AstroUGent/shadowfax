/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file AdaptiveVorTess2d.hpp
 *
 * @brief 2D evolving mesh: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_ADAPTIVEVORTESS2D
#define HEAD_ADAPTIVEVORTESS2D

#include "AdaptiveCellList.hpp"  // for AdaptiveCellList, etc
#include "Vec.hpp"               // for Vec
#include "utilities/Cuboid.hpp"  // for Cuboid
#include <ostream>               // for ostream
#include <vector>                // for vector

class AdaptiveFace2d;
class AdaptiveVorCell2d;
class ArgHilbertObject;
class GasParticle;
class ParticleVector;
class RestartFile;
class StateVector;
class VorTess;

/**
 * @brief 2D adaptive Voronoi tesselation
 *
 * Used to store the Voronoi tesselation for the mesh evolution algorithm
 */
class AdaptiveVorTess2d {
  private:
    /*! @brief List of cells */
    std::vector<AdaptiveVorCell2d*> _cells;
    //    std::vector<AdaptiveVorCell2d*> _ghosts;
    //    // _orphans are cells which do not have a local particle
    //    std::vector<AdaptiveVorCell2d*> _orphans;
    /*! @brief List of cells */
    AdaptiveCellList<AdaptiveVorCell2d> _newlist;

    /*! @brief List linking Particle information and cells */
    NewCellList _newerlist;

    /*! @brief Flag to signal if we reside in a periodic (true) or reflective
     *  (false) simulation box */
    bool _periodic;

    /*! @brief Cuboid specifying the dimensions of the simulation box */
    Cuboid _cuboid;

    /*! @brief Cuboid specifying the maximal box in which particles can be
     * added, needed for conversion of coordinates to the interval [1,2] */
    Cuboid _box12;

    /*! @brief List with ArgHilbertObjects, used to find out where particles
     *  moved to during the parallel sort */
    std::vector<ArgHilbertObject*> _indices;

    //    void save_restart(unsigned int nr = 0);

  public:
    AdaptiveVorTess2d(VorTess* tesselation, Cuboid cuboid, unsigned long maxid,
                      ParticleVector& particles, bool periodic = false);
    ~AdaptiveVorTess2d();

    unsigned int update(ParticleVector& particles, unsigned long currentTime);

    //    void print_tesselation(std::ostream &stream);
    void print_tesselation_gnuplot(std::ostream& stream);
//    void print_delaunay(std::ostream &stream);
//    void print_tesselation_pov(std::ostream &stream);
//    void print_tesselation_vtk(std::ostream &stream);
#ifndef ICMAKER
//    void print_tesselation_leaflet(ColorMap* colormap, StateVector maxW,
//                                   StateVector minW);
#endif
//    std::vector<VorCell*> get_cells();
//    std::vector<VorGen*> get_points();
//    void add_point(GasParticle* part);
#ifndef ICMAKER
//    void hydro(TimeLine& timeline, Solver& solver);
#endif
    //    void advect(double dt);
    //    void construct();
    //    void print_cell_statistics(std::ostream &stream);
    //    void set_hs();
    //    void complete(Tree& parttree);
    //    void update_Ws();
    //    void update_gradients();
    //    void update_dQs();
    //    void update_dts(unsigned long currentTime);
    //    VorCell* get_cell(unsigned int index);

    //    void check_delaunay();

    double get_volume(unsigned int i);
    double get_h(unsigned int i);
    Vec get_velocity(unsigned int i, GasParticle* particle);
    Vec get_centroid(unsigned int i);
    void estimate_gradients(unsigned int i, StateVector* delta,
                            GasParticle* particle);
    std::vector<AdaptiveFace2d*> get_faces(unsigned long currentTime);

    double get_total_volume();

    void update_positions(ParticleVector& particles);

    void check_mesh(ostream& stream);

    std::vector<unsigned long> get_ngb_ids(unsigned int index);

    void dump(RestartFile& rfile, std::vector<GasParticle*>& particles);
    AdaptiveVorTess2d(RestartFile& rfile, std::vector<GasParticle*>& particles);
};

#endif  // HEAD_ADAPTIVEVORTESS2D
