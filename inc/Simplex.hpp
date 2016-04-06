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
 * @file Simplex.hpp
 *
 * @brief A general simplex in X dimensions. If X=2, this is a triangle, if X=3,
 * we call it a tetrahedron: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_TETRAHEDRON
#define HEAD_TETRAHEDRON

#include "Error.hpp"
#include "ExArith.h"
#include "VorGen.hpp"
#include <iostream>
#include <map>
#include <vector>

/**
 * @brief 2D triangle or 3D tetrahedron
 *
 * The simplex is the basic building block of the DelTess. It contains indices
 * of 3 or 4 (depending on the dimension) VorGen elements of the DelTess and
 * links them together. It is also used to calculate the vertices of the
 * Voronoi tesselation, which are the midpoints of the circumspheres of the
 * simplices that make up the Delaunay tesselation.
 */
class Simplex {
  private:
    /*! \brief Indices of the vertices that make up the simplex in the DelTess
     *  VorGen list */
    unsigned int _vorgens[ndim_ + 1];

    /*! \brief Indices of the neighbouring simplices in the DelTess Simplex
     *  list */
    unsigned int _ngbs[ndim_ + 1];

    /*! \brief Indices in the neighbouring simplices of this simplex */
    int _ngbface[ndim_ + 1];

    /*! \brief Midpoint of the circumsphere of this simplex. This will be a
     *  vertex of the final Voronoi grid, if this Simplex is retained */
    VorGen* _midpoint_circumsphere;

    /*! \brief Linked list for DelTess Simplex queue: next simplex to check */
    unsigned int _next_check;

    /*! \brief Linked list for DelTess Simplex queue: previous simplex to
     *  check */
    unsigned int _previous_check;

  public:
#if ndim_ == 3
    /**
     * @brief Empty constructor
     *
     * Construct an empty Simplex
     */
    inline Simplex() {
        _vorgens[0] = 0;
        _vorgens[1] = 0;
        _vorgens[2] = 0;
        _vorgens[3] = 0;
        _ngbs[0] = 0;
        _ngbs[1] = 0;
        _ngbs[2] = 0;
        _ngbs[3] = 0;
        _ngbface[0] = -1;
        _ngbface[1] = -1;
        _ngbface[2] = -1;
        _ngbface[3] = -1;
        _midpoint_circumsphere = NULL;
    }

    /**
     * @brief Construct a tetrahedron in 3D with the four given points as
     * vertices
     *
     * @param id1 Index of the first point in the DelTess VorGen list
     * @param id2 Index of the second point in the DelTess VorGen list
     * @param id3 Index of the third point in the DelTess VorGen list
     * @param id4 Index of the fourth point in the DelTess VorGen list
     * @param points Reference to the DelTess VorGen list
     */
    inline Simplex(unsigned int id1, unsigned int id2, unsigned int id3,
                   unsigned int id4, std::vector<VorGen*>& points) {
        VorGen* point1 = points[id1];
        VorGen* point2 = points[id2];
        VorGen* point3 = points[id3];
        VorGen* point4 = points[id4];
        double det_or;
        det_or =
                ExactArithmetic::orient3d(point1->get_p12(), point2->get_p12(),
                                          point3->get_p12(), point4->get_p12());

        // to make things more interesting, orient3d uses the opposite notion of
        // positively oriented. We stick however to Springel's definition.
        if(det_or <= 0) {
            _vorgens[0] = id1;
            _vorgens[1] = id2;
            _vorgens[2] = id3;
            _vorgens[3] = id4;
        } else {
            _vorgens[0] = id1;
            _vorgens[1] = id2;
            _vorgens[2] = id4;
            _vorgens[3] = id3;
        }
        // initializing variables is always a good idea
        _ngbs[0] = 0;
        _ngbs[1] = 0;
        _ngbs[2] = 0;
        _ngbs[3] = 0;
        _ngbface[0] = -1;
        _ngbface[1] = -1;
        _ngbface[2] = -1;
        _ngbface[3] = -1;
        // _midpoint_circumsphere is empty for the time being
        _midpoint_circumsphere = NULL;
    }
#else
    Simplex(unsigned int id1, unsigned int id2, unsigned int id3);
#endif

    /**
     * @brief Destructor
     *
     * Clean up the midpoint, if it exists
     */
    ~Simplex() {
        // deleting a NULL pointer does not give rise to errors
        delete _midpoint_circumsphere;
    }

    int inside(VorGen* point, std::vector<VorGen*>& points);

#if ndim_ == 3
    /**
     * @brief Check if the given testpoint lies inside the face between an
     * invalid tetrahedron and its neighbour opposite to the newly added point,
     * to discriminate between a 2-to-3 and a 3-to-2 flip
     *
     * @param tid Index of the newly added point in the DelTess VorGen list
     * @param testpoint VorGen to be tested
     * @param indices Indices in the DelTess VorGen list of the vertices of the
     * line of the triangle closest to the testpoint (or all vertices of the
     * triangle in case testpoint lies on the boundary of the triangle)
     * @param points Reference to the DelTess VorGen list
     * @return 0, 2 or 3, depending on the number of indices stored in the array
     */
    inline unsigned int flip_test(unsigned int tid, VorGen* testpoint,
                                  unsigned int* indices,
                                  std::vector<VorGen*>& points) {
        // tid contains the index of toppoint in _points (we call it i for the
        // rest of this function)
        unsigned int result = 0;
        int j = 0;
        double test;
        while(!result && j < 3) {
            // for even i: have to check (i,i+1,i+2); (i,i+3,i+1); (i,i+2,i+3)
            // (always oriented counterclockwise with the upward direction
            // inside the tetrahedron)
            // for odd i: have to check (i,i+2,i+1); (i,i+1,i+3); (i,i+3,i+2)
            if(tid % 2 == 0) {
                test = ExactArithmetic::orient3d(
                        points[_vorgens[tid]]->get_p12(),
                        points[_vorgens[(tid + j + 1) % 4]]->get_p12(),
                        points[_vorgens[(tid + j + 2 + (j == 2)) % 4]]
                                ->get_p12(),
                        testpoint->get_p12());
            } else {
                test = ExactArithmetic::orient3d(
                        points[_vorgens[tid]]->get_p12(),
                        points[_vorgens[(tid + j + 2 + (j == 2)) % 4]]
                                ->get_p12(),
                        points[_vorgens[(tid + j + 1) % 4]]->get_p12(),
                        testpoint->get_p12());
            }
            if(test >= 0) {
                indices[0] = (tid + j + 1) % 4;
                indices[1] = (tid + j + 2 + (j == 2)) % 4;
                result = 2;
                if(test == 0) {
                    // we have to somehow signal this case; the most efficient
                    // way is by adding an element to result, so its size is >2
                    // we need the point of the triangle that is not on the side
                    // anyway, so may as well use that one
                    indices[2] = (tid + j + 3 + (j == 2 || j == 1)) % 4;
                    result = 3;
                }
            }
            j++;
        }
        return result;
    }
#endif

#if ndim_ == 3
    /**
     * @brief Test whether the given point is inside the circumsphere of the
     * tetrahedron
     *
     * To prevent roundoff error from skewing the result, we use arbitrary
     * precision arithmetics.
     *
     * @param point VorGen to test
     * @param points Reference to the DelTess VorGen list
     * @return True if the point lies inside the sphere, false otherwise
     */
    inline bool in_sphere(VorGen* point, std::vector<VorGen*>& points) {
        return ExactArithmetic::insphere(points[_vorgens[0]]->get_p12(),
                                         points[_vorgens[1]]->get_p12(),
                                         points[_vorgens[2]]->get_p12(),
                                         points[_vorgens[3]]->get_p12(),
                                         point->get_p12()) < 0;
    }
#else
    /**
     * @brief Test whether the given point is inside the circumcircle of the
     * triangle
     *
     * To prevent roundoff error from skewing the result, we use arbitrary
     * precision arithmetics.
     *
     * @param point VorGen to test
     * @param points Reference to the DelTess VorGen list
     * @return True if the point liest inside the circle, false otherwise
     */
    inline bool in_sphere(VorGen* point, std::vector<VorGen*>& points) {
        return ExactArithmetic::incircle(points[_vorgens[0]]->get_p12(),
                                         points[_vorgens[1]]->get_p12(),
                                         points[_vorgens[2]]->get_p12(),
                                         point->get_p12()) >= 0.;
    }
#endif

    /**
     * @brief Set the next Simplex to check in a linked list
     *
     * @param next_check Next Simplex to check
     */
    void set_next_check(unsigned int next_check) { _next_check = next_check; }

    /**
     * @brief Get the next Simplex to check in a linked list
     *
     * @return Next Simplex to check
     */
    unsigned int get_next_check() { return _next_check; }

    /**
     * @brief Set the previous Simplex to check in a linked list
     *
     * @param previous_check Previous Simplex to check
     */
    void set_previous_check(unsigned int previous_check) {
        _previous_check = previous_check;
    }

    /**
     * @brief Get the previous Simplex to check in a linked list
     *
     * @return Previous Simplex to check
     */
    unsigned int get_previous_check() { return _previous_check; }

    VorGen* get_special_point(std::vector<VorGen*>& points);

    /**
     * @brief Get a pointer to the internal vertex list
     *
     * @return A pointer to the internal vertex list
     */
    unsigned int* get_vorgens() { return _vorgens; }

    /**
     * @brief Get the vertex with the given index
     *
     * @param index Index of the vertex
     * @return Index of the requested vertex in the DelTess VorGen list
     */
    inline unsigned int vorgen(int index) { return _vorgens[index]; }

    void print(std::ostream& stream, std::vector<VorGen*>& points);

    /**
     * @brief Add a neighbour to the list using the given point for index lookup
     *
     * @param point Index of a vertex of this simplex in the DelTess VorGen list
     * @param simplex Index of a Simplex in the DelTess Simplex list
     * @param ngbindex Index of this simplex in the neighbouring Simplex
     * @return Index of the neighbour in the internal neighbour list
     */
    inline int add_ngb_from_vorgen(unsigned int point, unsigned int simplex,
                                   int ngbindex = 0) {
        int i = 0;
        while(_vorgens[i] != point) {
            i++;
        }
        _ngbs[i] = simplex;
        _ngbface[i] = ngbindex;
        return i;
    }

    /**
     * @brief Get the neighbour corresponding to the given point
     *
     * The neighbour corresponding to a point is the neighbour that shares the
     * face opposite to the point with this simplex.
     *
     * @param point Index of a vertex in the DelTess VorGen list
     * @return Index of the neighbouring Simplex in the DelTess Simplex list
     */
    inline unsigned int get_ngb_from_vorgen(unsigned int point) {
        unsigned int i = 0;
        while(i < ndim_ + 1 && _vorgens[i] != point) {
            i++;
        }
        if(i == ndim_ + 1) {
            // should not happen
            std::cerr << "Asked neighbour of point that is not part of "
                         "tetrahedron. Error!"
                      << std::endl;
            my_exit();
        }
        return _ngbs[i];
    }

    /**
     * @brief Add the given neighbour at the given position to the neighbour
     * list
     *
     * @param index Index of the new neighbour in the internal list
     * @param simplex Index of the neighbour in the DelTess Simplex list
     * @param ngbindex Index of this simplex in the neighbouring simplex
     */
    inline void add_ngb(int index, unsigned int simplex, int ngbindex = 0) {
        _ngbs[index] = simplex;
        _ngbface[index] = ngbindex;
    }

    /**
     * @brief Change the neighbour at the given index to the given value
     *
     * @param index Index of the neighbour in the internal list
     * @param new_value Index of the neighbour in the DelTess Simplex list
     * @param ngbindex Index of this simplex in the neighbouring simplex
     */
    inline void change_ngb(unsigned int index, unsigned int new_value,
                           int ngbindex = 0) {
        _ngbs[index] = new_value;
        _ngbface[index] = ngbindex;
    }

    /**
     * @brief Get the neighbour at the given index
     *
     * @param index Index of the neighbour in the internal list
     * @return Index of the neighbour in the DelTess Simplex list
     */
    inline unsigned int get_ngb(int index) { return _ngbs[index]; }

    /**
     * @brief Get the index of the given point in the internal vertex list
     *
     * @param point Index of the point in the DelTess VorGen list
     * @return Index of the point in the internal vertex list
     */
    inline int get_index(unsigned int point) {
        int i = ndim_;
        while(i >= 0 && _vorgens[i] != point) {
            i--;
        }
        if(i < 0) {
            std::cerr << "Asked index of point that is not inside simplex. "
                         "Error!"
                      << std::endl;
            my_exit();
        }
        return i;
    }

    /**
     * @brief Get the index of this simplex in the neighbour with the given
     * index in the neighbour list
     *
     * @param ngb Index of a simplex in the internal neighbour list
     * @return Index of this simplex in the neighbour list of the given
     * neighbour
     */
    inline int get_ngbindex(unsigned int ngb) { return _ngbface[ngb]; }

    void get_centroid(double* centroid, std::vector<VorGen*>& points);

    unsigned int get_idsum(std::vector<VorGen*>& points);

    /**
     * @brief Get the index of this simplex in the neighbour with the given
     * index in the neighbour list
     *
     * @param index Index of a simplex in the internal neighbour list
     * @return Index of this simplex in the neighbour list of the given
     * neighbour
     */
    inline unsigned int get_ngbface(unsigned int index) {
        return _ngbface[index];
    }

    /**
     * @brief Get the index of this simplex in the neighbour corresponding to
     * the given point
     *
     * The neighbour corresponding to the point is the neighbour that shares the
     * face opposite to point with this simplex.
     *
     * @param point Index of a vertex in the DelTess VorGen list
     * @return Index of this simplex in the neighbour corresponding to the
     * given point
     */
    inline unsigned int get_ngbface_from_vorgen(unsigned int point) {
        unsigned int i = 0;
        while(_vorgens[i] != point) {
            i++;
        }
        return _ngbface[i];
    }
};

#endif
