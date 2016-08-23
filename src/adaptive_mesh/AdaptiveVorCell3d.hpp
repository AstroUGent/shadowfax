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
 * @file AdaptiveVorCell3d.hpp
 *
 * @brief 3D mesh evolution Voronoi cell: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_ADAPTIVEVORCELL3D
#define HEAD_ADAPTIVEVORCELL3D

#include "Vec.hpp"  // for Vec
#include <cstddef>  // for NULL
#include <ostream>  // for ostream
#include <vector>   // for vector

class AdaptiveVorFace3d;
class Cuboid;
class GasParticle;
class RestartFile;
class StateVector;

/**
 * @brief 3D adaptive voronoi cell
 *
 * Used in the 3D version of the grid evolution algorithm to represent a single
 * Voronoi cell.
 */
class AdaptiveVorCell3d {
  private:
    /*! @brief List of indices of neighbouring cells in the AdaptiveVorTess3d */
    std::vector<int> _ngbs;

    /*! @brief Boundary flags for the neighbours */
    std::vector<int> _ngbwalls;

    // facengbs: for every ngb, we store a list of facengbs. This list is
    // ordered w.r.t. the vector pointing from _pos to the position of the
    // corresponding neighbour generator
    // for every facengb, we also store the associated wall key
    /*! @brief List of facengbs for every neighbour */
    std::vector<std::vector<int> > _facengbs;

    /*! @brief List of boundary flags for every facengb */
    std::vector<std::vector<int> > _facengbwalls;

    // vector of faces of this cell
    // when explicitly filled, _faces[i] holds the face for _ngbs[i]
    // however, this vector is not always up-to-date!
    /*! @brief Face of the cell */
    std::vector<AdaptiveVorFace3d*> _faces;

    /*! @brief Array containing information about the last flip that was
     *  performed in this cell */
    unsigned int _flip_info[2];

    /*! @brief Actual position of the generator of the cell */
    double _pos[3];

    /*! @brief Last valid position of the generator of this cell */
    double _valid_pos[3];

    /*! @brief id of the cell in the cell list */
    unsigned int _id;

    /*! @brief GasParticle associated with the cell */
    GasParticle* _particle;

    /*! @brief Flag signaling if this cell should be considered during the
     *  second detection loop */
    bool _flag1;

    /*! @brief Flag signaling if this cell should be considered during the third
     *  and last detection loop */
    bool _flag2;

    /*! @brief Flag signaling if this cell should be considered during the first
     *  detection loop */
    bool _active;

    bool get_vertexpoint(double* a, double* b, double* c, double* d,
                         double* vert);
    bool inside(double* a, double* b, double* c, double* d, double* point);

  public:
    AdaptiveVorCell3d(RestartFile& rfile, GasParticle* particle);
    void dump(RestartFile& rfile);

    AdaptiveVorCell3d(double* pos, unsigned int id, GasParticle* particle);
    void cleanup_faces();
    ~AdaptiveVorCell3d();

    unsigned int add_ngb(int id, int wall = 0);
    void add_facengb(int id, int wall = 0);
    void print_relevant_info(std::ostream& stream,
                             std::vector<AdaptiveVorCell3d*>& cells,
                             std::vector<AdaptiveVorCell3d*>& ghosts,
                             Cuboid& box);
    void complete(std::vector<AdaptiveVorCell3d*>& cells,
                  std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box);
    void print(std::ostream& stream, int id = -1);
    void print_copy(std::ostream& stream, int id, int wall, Cuboid& box);
    void move(double* pos);
    bool check(unsigned int id, std::vector<AdaptiveVorCell3d*>& cells,
               std::ostream& estream, Cuboid& box);
    void get_position(double* pos);
    bool is_ngb(int id);
    unsigned int remove_ngb(unsigned int id);
    void remove_facengb(int ngb, int facengb);
    void remove_facengb(unsigned int* edges, int ngb, int facengb);
    void safely_remove_facengb(int ngb, int facengb);
    void add_facengb(unsigned int* edges, int ngb, int facengb, int before,
                     int wall = 0);
    void add_facengb(int ngb, int facengb, int before, int wall = 0);
    bool create_new_face(unsigned int* edges_to_check, int A, int B, int C,
                         int D, int E, int wallA, int wallB, int wallC,
                         int wallD, int wallE,
                         std::vector<AdaptiveVorCell3d*>& cells,
                         std::vector<AdaptiveVorCell3d*>& ghosts,
                         bool periodic = false);
    bool remove_face(unsigned int* edges_to_check, int A, int B, int C, int D,
                     int E, int wallA, int wallB, int wallC, int wallD,
                     int wallE, std::vector<AdaptiveVorCell3d*>& cells,
                     std::vector<AdaptiveVorCell3d*>& ghosts,
                     bool periodic = false);
    unsigned int check_triangular_face(int A, int B, int C,
                                       std::vector<AdaptiveVorCell3d*>& cells);
    //    void trim_first_position(double* pos, Cuboid& box, bool* trims);
    void trim_position(double* pos, Cuboid& box, bool* trims);
    bool in_sphere_test(int A, int B, int C, int D, int E, int wallA, int wallB,
                        int wallC, int wallD, int wallE,
                        std::vector<AdaptiveVorCell3d*>& cells,
                        std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                        bool periodic = false, bool verbose = false,
                        double* testdouble = NULL);
    unsigned int get_flips(int id, std::vector<AdaptiveVorCell3d*>& cells,
                           std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                           bool periodic = false);
    unsigned int get_flips(int id, int ngb,
                           std::vector<AdaptiveVorCell3d*>& cells,
                           std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                           bool periodic = false);
    unsigned int get_flips_minimal(unsigned int* flips, int id,
                                   std::vector<AdaptiveVorCell3d*>& cells,
                                   std::vector<AdaptiveVorCell3d*>& ghosts,
                                   Cuboid& box, bool periodic = false);
    unsigned int get_flips_minimal(unsigned int* flips, int id, int ngb,
                                   std::vector<AdaptiveVorCell3d*>& cells,
                                   std::vector<AdaptiveVorCell3d*>& ghosts,
                                   Cuboid& box, bool periodic = false);
    bool detect_crossovers(unsigned int id,
                           std::vector<AdaptiveVorCell3d*>& cells,
                           std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                           bool periodic = false);
    void print_copies(std::ostream& stream, std::vector<int>& ngblist,
                      std::vector<AdaptiveVorCell3d*>& cells, Cuboid& box);
    void swap_wallpos(int id, int wall);
    void keep_inside(int id, std::vector<AdaptiveVorCell3d*>& cells,
                     Cuboid& box);

    std::vector<int> get_ngbs();
    double get_h();
    double get_volume();
    Vec get_centroid();
    GasParticle* get_particle();
    AdaptiveVorFace3d* get_face(unsigned int i);
    void get_ghost_position(double* pos, Cuboid& box);
    void get_periodic_position(unsigned int i,
                               std::vector<AdaptiveVorCell3d*>& cells,
                               double* posvec, Cuboid& box);
    Vec get_velocity(GasParticle* particle);
    void estimate_gradients(StateVector* delta, GasParticle* particle,
                            std::vector<AdaptiveVorCell3d*>& cells,
                            std::vector<AdaptiveVorCell3d*>& ghosts,
                            Cuboid& box, bool periodic = false);

    void queue();
    bool queued();
    void activate();
    bool active();
    void semiactivate();
    bool semiactive();
    void reset_flags();

    void set_valid_pos(double* valid_pos);
    void get_valid_pos(double* valid_pos);

    void remove_face(unsigned int* flips, int id,
                     std::vector<AdaptiveVorCell3d*>& cells,
                     std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                     bool periodic = false);
    void insert_face(unsigned int* flips, int id,
                     std::vector<AdaptiveVorCell3d*>& cells,
                     std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                     bool periodic = false);
    void flip_face(int id, std::vector<AdaptiveVorCell3d*>& cells,
                   std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                   bool periodic = false);

    void remove_face(unsigned int* flips, int id, unsigned int face,
                     std::vector<AdaptiveVorCell3d*>& cells,
                     std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                     bool periodic = false);
    void insert_face(unsigned int* flips, int id, unsigned int face,
                     unsigned int edge, std::vector<AdaptiveVorCell3d*>& cells,
                     std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                     bool periodic = false);
    void flip_face(int id, unsigned int face, unsigned int edge,
                   std::vector<AdaptiveVorCell3d*>& cells,
                   std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                   bool periodic = false);

    void check_cells(int* ids, unsigned int numid,
                     std::vector<AdaptiveVorCell3d*>& cells,
                     std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                     bool periodic = false);
    void check_cells(int* ids, unsigned int numid, unsigned int* edges,
                     unsigned int numedge,
                     std::vector<AdaptiveVorCell3d*>& cells,
                     std::vector<AdaptiveVorCell3d*>& ghosts, Cuboid& box,
                     bool periodic = false);

    void update_ngbs(std::vector<unsigned int>& new_ids);
    void update_ghosts(std::vector<unsigned int>& new_ids,
                       unsigned int cellsize);

    unsigned int get_id();

    void clean_facengbs();
};

#endif
