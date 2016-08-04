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
 * @file VorTess.hpp
 *
 * @brief Voronoi tesselation for the old algorithm: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_VORTESS
#define HEAD_VORTESS

#include "StateVector.hpp"  // for StateVector
#include <ostream>          // for ostream
#include <vector>           // for vector

class ColorMap;
class DelCont;
class DelTess;
class GasParticle;
class ParticleVector;
class RiemannSolver;
class TimeLine;
class Tree;
class VorCell;
class VorFace;

/**
 * @brief 2D or 3D Voronoi tesselation
 *
 * Used for the old algorithm.
 */
class VorTess {
  private:
    /*! \brief Delaunay tesselation of the point set */
    DelTess* _delaunay;

    /*! \brief Cells of the tesselation */
    std::vector<VorCell*> _cells;
    /*! \brief Faces of the tesselation */
    std::vector<VorFace*> _faces;

  public:
    VorTess(DelCont* cont, unsigned int numpart, bool periodic = false,
            double tolerance = 1.e-9);
    ~VorTess();

    void print_tesselation(std::ostream& stream);
    void print_tesselation_gnuplot(std::ostream& stream);
    void print_delaunay(std::ostream& stream);
    void print_tesselation_pov(std::ostream& stream);
    void print_tesselation_vtk(std::ostream& stream);
    void print_tesselation_fastvor(std::ostream& stream, bool periodic,
                                   bool binary = true);
#ifndef ICMAKER
    void print_tesselation_leaflet(ColorMap* colormap, StateVector maxW,
                                   StateVector minW);
#endif
    std::vector<VorCell*> get_cells();
    void add_point(GasParticle* part, unsigned int index);
#ifndef ICMAKER
    void hydro(TimeLine& timeline, RiemannSolver& solver,
               ParticleVector& particles);
#endif
    void advect(double dt);
    void construct();
    void print_cell_statistics(std::ostream& stream);
    void set_hs();
    void complete(Tree& parttree);
    void update_Ws();
    void update_gradients();
    void update_dQs();
    void update_dts(unsigned long currentTime);
    void update_gravitational_corrections();
    VorCell* get_cell(unsigned int index);

    void check_delaunay();

    std::vector<VorFace*>& get_faces();

    void get_triangles(std::vector<float>& positions,
                       std::vector<int>& connectivity,
                       std::vector<StateVector>& data);

    void get_delaunay_triangles(std::vector<float>& positions,
                                std::vector<int>& connectivity);

    std::vector<unsigned long> get_ngb_ids(VorCell* cell);
};

#endif
