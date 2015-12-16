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
 * @file VorCell.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief contains class VorCell
 *
 * A VorCell holds all information concering a voronoi cell: the central point
 * and the faces (3D) or edges (2D). It allows for the calculation of
 * cell-specific properties like the centroid.
 */
#ifndef HEAD_VORCELL
#define HEAD_VORCELL

#include "Vec.hpp"
#include <vector>
#include <ostream>
#include <map>
#include "Image.hpp"
#include "StateVector.hpp"

class Simplex;
class VorFace;
class VorGen;
class GasParticle;

/**
 * @brief 2D or 3D Voronoi cell
 *
 * Stores the list of VorFaces for this cell and performs the volume and
 * centroid calculations as well as the gradient estimation and generator
 * velocity computation.
 */
class VorCell{
private:
    /*! \brief Faces that border the cell */
    std::vector<VorFace*> _faces;
    /*! \brief Neighbouring points of this cell */
    std::vector<VorGen*> _ngbs;
    /*! \brief Generator of this cell */
    VorGen* _central_point;

    /*! \brief Indices of the faces in the VorTess face list */
    std::vector<unsigned int> _face_ids;
    /*! \brief Indices of the neighbours in the DelTess point list */
    std::vector<unsigned int> _ngb_ids;

    /*! \brief Centroid of the cell */
    Vec _centroid;
    /*! \brief Volume of the cell */
    double _volume;
    /*! \brief Total area of the faces of the cell */
    double _total_area;

public:
    VorCell(VorGen* point);
    ~VorCell();

    void print(std::ostream &stream);
    void print_gnuplot(std::ostream &stream, unsigned int id = 0);
    void print_pov(std::ostream &stream);
#ifndef ICMAKER
    void print_leaflet(std::ostream &vstream, int ox, int oy,
                       ColorMap* colormap, StateVector maxW, StateVector minW);
#endif
    bool overlap(double* box);
    void print_vtk(std::ostream &vstream, unsigned int &numv,
                   std::ostream &pstream, unsigned int &nump,
                   unsigned int &numc, std::ostream &dstream);
    void print_fastvor(std::ostream& stream, unsigned int &numghost,
                       bool periodic, bool binary = true);

    Vec& get_centroid();
    void set_centroid(double* centroid);
    void set_centroid(Vec centroid);
    void calculate_centroid();

    double get_volume();
    void set_volume(double volume);
    void calculate_volume();

    double get_total_area();
    void calculate_total_area();

    void add_ngb(VorGen* ngb);
    std::vector<VorGen*> get_ngbs();

    void add_face(VorFace* face);
    int number_of_faces();
    std::vector<VorFace*> get_faces();
    VorFace* get_face(VorGen* point);

    void add_ngb_id(unsigned int ngb);
    std::vector<unsigned int> get_ngb_ids();

    void add_face_id(unsigned int face);
    std::vector<unsigned int> get_face_ids();

    void estimate_gradient(StateVector* delta);
    void set_h();
    double get_h();
    Vec get_velocity();

    GasParticle* get_particle();

    double overlap_volume(double* corner, double side);
    double periodic_overlap_volume(double* corner, double side);

    void finalize_eta();
    Vec get_gravitational_correction();

    void get_triangles(std::vector<float>& positions,
                       std::vector<int>& connectivity,
                       std::vector<StateVector>& data);
};

#endif
