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
 * @file AdaptiveVorCell2d.hpp
 *
 * @brief 2D mesh evolution Voronoi cell: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_ADAPTIVEVORCELL2D
#define HEAD_ADAPTIVEVORCELL2D

#include "Vec.hpp"                // for Vec
#include "utilities/Hilbert.hpp"  // for Hilbert_Object
#include <cstddef>                // for NULL
#include <ostream>                // for ostream
#include <vector>                 // for vector

class AdaptiveVorFace2d;
class Cuboid;
class GasParticle;
class NewCellList;
class RestartFile;
class StateVector;
template <class Cell> class AdaptiveCellList;

/**
 * @brief 2D adaptive Voronoi cell
 *
 * Used to represent a cell of the Voronoi tesselation for the mesh evolution
 * algorithm.
 */
class AdaptiveVorCell2d : public Hilbert_Object {
  private:
    /*! @brief List of neighbours of this cell */
    std::vector<unsigned long> _newngbs;

    /*! @brief ID of the particle corresponding to this cell */
    unsigned long _particleID;

    /*! @brief List of neighbours of this cell */
    std::vector<int> _ngbs;

    /*! @brief List of neighbour wall flags */
    std::vector<int> _ngbwalls;

    /*! @brief Deprecated list with unclear purpose */
    std::vector<double> _orientation;

    /*! @brief Another deprecated list with unclear purpose */
    std::vector<int> _exngbs;

    /*! @brief List of faces of this cell */
    std::vector<AdaptiveVorFace2d*> _faces;

    /*! @brief Position of the generator of this cell */
    double _pos[2];

    /*! @brief Index of the original cell if this cell represents a ghost */
    unsigned int _oid;

    /*! @brief True if this cell is adjacent to a wall of the simulation box */
    bool _wall;

    /*! @brief GasParticle associated with this cell */
    GasParticle* _particle;

    /*! @brief Flag set on the first pass of the mesh update algorithm, when
     *  all neighbours of active cells are activated */
    bool _flag1;

    /*! @brief Flag set on the second pass of the mesh update algorithm, when
     *  all neighbours of semiactive cells are activated */
    bool _flag2;

    /*! @brief Flag set when a cell is active */
    bool _active;

    /*! @brief ID of the cell */
    unsigned long _id;

    /*! @brief The rank of the process this cell belongs to */
    // only interesting for orphan cells
    // initialized as -1, which means this rank corresponds to the local process
    int _rank;

    //    int get_common_ngb(unsigned int idnot, AdaptiveVorCell2d* a,
    //    AdaptiveVorCell2d* b);

  public:
    AdaptiveVorCell2d(Vec pos, unsigned int nngb, unsigned int oid = 0,
                      GasParticle* particle = NULL);
    AdaptiveVorCell2d(double* pos, unsigned int nngb, unsigned int oid = 0,
                      GasParticle* particle = NULL);
    ~AdaptiveVorCell2d();

    unsigned int get_original();
    std::vector<int> get_ngbs();
    std::vector<int> get_ngbwalls();
    void change_ngb(unsigned int index, int new_ngb);
    //    void change_ngb_from_value(int value, int new_ngb);
    void safely_change_ngb_from_value(int value, int new_ngb);
    void add_ngb(int id);
    void add_ngbwall(int id);
    //    void set_orientations();
    //    void complete(std::vector<AdaptiveVorCell2d*>& cells,
    //    std::vector<AdaptiveVorCell2d*>& ghosts,
    //    std::vector<AdaptiveVorCell2d*>& orphans,
    //    Cuboid &box, bool periodic = false);
    void complete(AdaptiveCellList<AdaptiveVorCell2d>& cells, Cuboid& box,
                  bool periodic = false);
    void complete(NewCellList& list);
    void print(std::ostream& stream, int id);
    //    void fill_exngbs(unsigned int id, std::vector<AdaptiveVorCell2d*>&
    // cells);
    //    void print_exngbs(std::ostream& stream,
    //    std::vector<AdaptiveVorCell2d*>& cells);
    void get_centroid(double* centroid);
    void move(double* pos);
    //    bool check(unsigned int id, std::vector<AdaptiveVorCell2d*>& cells,
    //    std::ostream& estream);
    void get_position(double* pos);
    void remove_ngb(unsigned int id, bool periodic = false);
    void remove_newngb(unsigned long id, bool periodic = false);
    void safely_remove_ngb(unsigned int id, bool periodic = false);
    void add_ngb(unsigned int id1, unsigned int id2, int ngb);
    void add_newngb(unsigned long id1, unsigned long id2, unsigned long ngb);
    void add_ngb(unsigned int id1, unsigned int id2, int ngb, int wallpos);
    void get_ngbs(unsigned int id, int* ngbs);
    void compute_winding_number();
    bool detect_crossovers(unsigned int id,
                           AdaptiveCellList<AdaptiveVorCell2d>& newlist,
                           std::vector<int>& insertions, Cuboid& box,
                           bool periodic = false);

    bool detect_crossovers(unsigned int id,
                           std::vector<AdaptiveVorCell2d*>& cells,
                           NewCellList& positionlist, Cuboid& box,
                           std::vector<unsigned long>& insertions,
                           bool periodic = false);

    bool is_ngb(unsigned int id);
    //    bool check_ngbs(unsigned int id, std::vector<AdaptiveVorCell2d*>&
    // cells);

    void add_ngb_before(unsigned int id, int ngb);
    void add_ngb_before(unsigned int id, int ngb, int wallpos);
    void add_ngb_after(unsigned int id, int ngb);
    void add_ngb_after(unsigned int id, int ngb, int wallpos);

    //    void save_restart(std::ostream& stream);
    //    AdaptiveVorCell2d(std::istream& stream);

    double get_volume();
    double get_h();
    Vec get_velocity(GasParticle* particle);
    void estimate_gradients(StateVector* delta, GasParticle* particle,
                            AdaptiveCellList<AdaptiveVorCell2d>& cells,
                            Cuboid& box, bool periodic = false);
    void estimate_gradients(StateVector* delta, GasParticle* particle,
                            NewCellList& positionlist, bool periodic = false);
    GasParticle* get_particle();
    AdaptiveVorFace2d* get_face(unsigned int i);

    void get_ghost_position(double* ngbposvec, Cuboid& box);
    void apply_periodic_boundaries(unsigned int id,
                                   AdaptiveCellList<AdaptiveVorCell2d>& cells,
                                   Cuboid& box);

    void swap_wallpos(unsigned int id, int wallpos);

    void get_periodic_position(unsigned int index,
                               AdaptiveCellList<AdaptiveVorCell2d>& cells,
                               double* pos, Cuboid& box);

    void queue();
    bool queued();
    void activate();
    void deactivate();
    bool active();
    void semiactivate();
    bool semiactive();
    void reset_flags();

    void update_ngbs(vector<unsigned int>& new_ids);
    void update_local_id(vector<unsigned int>& new_ids);
    void update_local_id(unsigned int new_id);
    void update_orphans(vector<unsigned int>& new_ids,
                        AdaptiveCellList<AdaptiveVorCell2d>& cells);

    int get_first_ngb();

    void set_id(unsigned long id);
    unsigned long get_id();

    void set_rank(int rank);
    int get_rank();

    void set_particleID(unsigned long particleID);
    unsigned long get_particleID();

    std::vector<unsigned long>& get_newngbs();
    void add_ngb_particleID(unsigned long particleID);

    void dump(RestartFile& rfile);
    AdaptiveVorCell2d(RestartFile& rfile, GasParticle* particle = NULL);

    void pack_data(void* buffer, int bufsize, int* position,
                   AdaptiveCellList<AdaptiveVorCell2d>& cells);
    AdaptiveVorCell2d(void* buffer, int bufsize, int* position,
                      AdaptiveCellList<AdaptiveVorCell2d>& cells);

    void pack_data(void* buffer, int bufsize, int* position);
    AdaptiveVorCell2d(void* buffer, int bufsize, int* position);
};

#endif  // HEAD_ADAPTIVEVORCELL2D
