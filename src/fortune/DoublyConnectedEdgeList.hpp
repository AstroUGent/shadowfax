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
 * @file DoublyConnectedEdgeList.hpp
 *
 * @brief Geometrical structure used to store Fortune's Voronoi mesh
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef DOUBLYCONNECTEDEDGELIST_HPP
#define DOUBLYCONNECTEDEDGELIST_HPP

#if ndim_ == 2

#include "Vec.hpp"  // for Vec
#include <cmath>    // for atan2
#include <cstddef>  // for NULL
#include <ostream>  // for ostream
#include <vector>   // for vector

class DelCont;
class Face;

/**
 * @brief Vertex of the Voronoi diagram in Fortune's algorithm
 */
class Vertex {
  private:
    /*! @brief Position of the vertex */
    Vec _position;

    /*! @brief Number of references to this vertex */
    unsigned int _refs;

  public:
    /**
     * @brief Constructor
     *
     * @param position Position of the vertex
     */
    Vertex(Vec position) : _position(position), _refs(0) {}

    ~Vertex() {}

    /**
     * @brief Get the position of this vertex
     *
     * @return Position of the vertex
     */
    Vec get_position() {
        return _position;
    }

    /**
     * @brief Change the position of the vertex to the given value
     *
     * @param position New position for the vertex
     */
    void reset(Vec position) {
        _position = position;
    }

    /**
     * @brief Reference this vertex
     *
     * Every object that stores a pointer to this particular vertex should call
     * this method exactly once.
     */
    void reference() {
        _refs++;
    }

    /**
     * @brief Dereference this vertex
     *
     * Every object that stores a pointer to this particular vertex should call
     * this method exactly once. If the return value is positive, it should
     * delete the vertex.
     *
     * @return True if the vertex is no longer referenced by any other object
     * and should be deleted, false otherwise
     */
    bool dereference() {
        _refs--;
        return !_refs;
    }
};

/**
 * @brief Edge of the Voronoi diagram in Fortune's algorithm
 */
class Edge {
  private:
    /*! @brief Vertex that acts as origin for this edge */
    Vertex* _origin;

    /*! @brief Twin edge that has the other vertex of the edge as origin */
    Edge* _twin;

    /*! @brief Next edge that has the other vertex of the edge as origin but
     *  is no twin of this edge (the origin of this edge is no vertex of next
     *  edge) */
    Edge* _next_edge;

    /*! @brief Previous edge that has the origin of this edge as a vertex (but
     *  not as origin) and is no twin of this edge (the origin of previous edge
     *  is not a vertex of this edge) */
    Edge* _prev_edge;

    /*! @brief Cell corresponding to this edge */
    Face* _face;

  public:
    /**
     * @brief Constructor
     *
     * Initializes an empty Edge,
     */
    Edge()
            : _origin(NULL), _twin(NULL), _next_edge(NULL), _prev_edge(NULL),
              _face(NULL) {}

    /**
     * @brief Destructor
     *
     * Dereference the origin and delete it if necessary. We have to be
     * careful about deleting the origin, since it is also the origin of the
     * twin of the previous edge.
     */
    ~Edge() {
        if(_origin && _origin->dereference()) {
            delete _origin;
        }
    }

    /**
     * @brief Set the origin of this edge, if it does not exist yet
     *
     * The origin Vertex is referenced to ensure proper memory handling.
     *
     * @param origin Vertex that is the origin of this edge
     */
    void set_origin(Vertex* origin) {
        if(!_origin) {
            _origin = origin;
            origin->reference();
        }
    }

    /**
     * @brief Get the position of the origin of this edge
     *
     * @return The position of the origin, a Vec
     */
    Vec get_origin() {
        return _origin->get_position();
    }

    /**
     * @brief Set the twin edge for this edge
     *
     * @param twin Twin edge for this edge
     */
    void set_twin(Edge* twin) {
        _twin = twin;
    }

    /**
     * @brief Get the twin edge for this edge
     *
     * @return Twin edge for this edge
     */
    Edge* get_twin() {
        return _twin;
    }

    void plot(std::ostream& stream);
    void plot_gnuplot(std::ostream& stream);

    bool trim_edge(DelCont* container, std::vector<Vertex*>& vertices,
                   std::vector<Edge*>& edges);

    /**
     * @brief Set the next edge of this edge
     *
     * @param next_edge Next edge of this edge
     */
    void set_next_edge(Edge* next_edge) {
        _next_edge = next_edge;
    }

    /**
     * @brief Get the next edge of this edge
     *
     * @return The next edge of this edge
     */
    Edge* get_next_edge() {
        return _next_edge;
    }

    /**
     * @brief Set the previous edge of this edge
     *
     * @param prev_edge Previous edge of this edge
     */
    void set_prev_edge(Edge* prev_edge) {
        _prev_edge = prev_edge;
    }

    /**
     * @brief Get the previous edge of this edge
     *
     * @return The previous edge of this edge
     */
    Edge* get_prev_edge() {
        return _prev_edge;
    }

    bool outside(DelCont* container);

    /**
     * @brief Set the cell of this edge
     *
     * @param face Cell of this edge
     */
    void set_face(Face* face) {
        _face = face;
    }

    /**
     * @brief Get the cell of this edge
     *
     * @return The cell of this edge
     */
    Face* get_face() {
        return _face;
    }

    /**
     * @brief Check if this edge has a valid origin
     *
     * @return True if the origin is set, false otherwise
     */
    bool has_origin() {
        return _origin != NULL;
    }
};

/**
 * @brief Cell of the Voronoi diagram in Fortune's algorithm
 *
 * A cell has an area and a centroid.
 */
class Face {
  private:
    /*! @brief One of the edges associated with this cell */
    Edge* _edge;

    /*! @brief Geometrical area of this cell */
    double _area;

    /*! @brief Geometrical centroid of this cell */
    Vec _centroid;

  public:
    /**
     * @brief Constructor
     *
     * Initializes a cell with given reference Edge and zero area and centroid.
     *
     * @param edge Edge associated with this cell
     */
    Face(Edge* edge) : _edge(edge) {
        _area = 0.;
    }

    ~Face() {}

    void calculate_quantities();

    /**
     * @brief Get the area of this cell
     *
     * @warning The area has to be calculated first by calling
     * Face::calculate_quantities().
     *
     * @return The area of the cell
     */
    double get_area() {
        return _area;
    }

    /**
     * @brief Get the centroid of this cell
     *
     * @warning The centroid has to be calculated first by calling
     * Face::calculate_quantities().
     *
     * @return The centroid of the cell, a Vec
     */
    Vec get_centroid() {
        return _centroid;
    }
};

/**
 * @brief Class used to sort vertices on their angle around the origin
 */
class VertexSorter {
  private:
    /*! @brief List of corresponding angles on which to sort */
    std::vector<double> _angles;

  public:
    /**
     * @brief Constructor
     *
     * Initialize the reference to the list to sort (why?) and fill the list
     * with corresponding angles around the origin (why?).
     *
     * @param vertices List of vertices to sort on angle.
     */
    VertexSorter(std::vector<Vertex*>& vertices) {
        _angles.resize(vertices.size(), 0.);
        for(unsigned int i = 0; i < vertices.size(); i++) {
            _angles[i] = atan2(vertices[i]->get_position().y() - 0.5,
                               vertices[i]->get_position().x() - 0.5);
        }
    }

    ~VertexSorter() {}

    /**
     * @brief Functor used for sorting on angle
     *
     * @param i unsigned integer index of a vertex in the list
     * @param j unsigned integer index of a vertex in the list
     * @return True if the angle corresponding to the first vertex is smaller
     * than that corresponding to the second vertex, false otherwise
     */
    bool operator()(unsigned int i, unsigned int j) {
        return _angles[i] < _angles[j];
    }
};

/**
 * @brief Voronoi diagram abstraction used for Fortune's algorithm
 */
class DoublyConnectedEdgeList {
  private:
    /*! @brief List of edges in the Voronoi tesselation */
    std::vector<Edge*> _edges;

    /*! @brief List of cells in the Voronoi tesselation */
    std::vector<Face*> _faces;

  public:
    DoublyConnectedEdgeList() {}
    ~DoublyConnectedEdgeList();

    void add_edge(Edge* edge);

    void plot(std::ostream& stream);
    void plot_gnuplot(std::ostream& stream);

    void trim_edges(DelCont* container);

    void print_cell(std::ostream& stream, Edge* edge);

    void test_connections(std::ostream& stream);

    void set_faces();
};

#endif

#endif  // DOUBLYCONNECTEDEDGELIST_HPP
