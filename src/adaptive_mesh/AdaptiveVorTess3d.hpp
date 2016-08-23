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
 * @file AdaptiveVorTess3d.hpp
 *
 * @brief 3D evolving mesh: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_ADAPTIVEVORTESS3D
#define HEAD_ADAPTIVEVORTESS3D

#include "Vec.hpp"               // for Vec
#include "utilities/Cuboid.hpp"  // for Cuboid
#include <ostream>               // for ostream
#include <vector>                // for vector

class AdaptiveFace3d;
class AdaptiveVorCell3d;
class GasParticle;
class ParticleVector;
class RestartFile;
class StateVector;
class VorTess;

/**
 * @brief 3D adaptive Voronoi tesselation
 *
 * The 3D Voronoi mesh that is used in the mesh evolution algorithm.
 */
class AdaptiveVorTess3d {
  private:
    /*! @brief List of normal cells */
    std::vector<AdaptiveVorCell3d*> _cells;

    /*! @brief List of ghost cells */
    std::vector<AdaptiveVorCell3d*> _ghosts;

    /*! @brief List of new positions for the generators */
    std::vector<Vec> _new_positions;

    /*! @brief Flag specifying if we deal with a periodic (true) or reflective
     *  (false) simulation box */
    bool _periodic;

    /*! @brief Cuboid specifying the dimensions of the simulation box */
    Cuboid _cuboid;

    /*! @brief Tolerance parameter used internally */
    double _tolerance;

    //    void save_restart(unsigned int nr = 0);

  public:
    AdaptiveVorTess3d(VorTess* tesselation, Cuboid cuboid,
                      bool periodic = false);
    ~AdaptiveVorTess3d();

    bool valid_flips(unsigned int* flips);
    bool valid_flips_minimal(unsigned int* flips);
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
    std::vector<AdaptiveFace3d*> get_faces(unsigned long currentTime);

    //    double get_total_volume();

    void update_positions(ParticleVector& particles);

    void check_mesh(ostream& stream);

    std::vector<unsigned long> get_ngb_ids(unsigned int index);

    vector<Vec> get_new_positions();

    void dump(RestartFile& rfile, std::vector<GasParticle*>& particles);
    AdaptiveVorTess3d(RestartFile& rfile, std::vector<GasParticle*>& particles);
};

#endif
