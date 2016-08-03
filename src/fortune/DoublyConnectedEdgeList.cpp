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
 * @file DoublyConnectedEdgeList.cpp
 *
 * @brief Geometrical structure used to store Fortune's Voronoi mesh
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#if ndim_ == 2

#include "DoublyConnectedEdgeList.hpp"
#include "DelCont.hpp"
#include <algorithm>
#include <iostream>
using namespace std;

/**
 * @brief Write the edge in plottable format to the given stream
 *
 * @param stream std::ostream to write to
 */
void Edge::plot(ostream& stream) {
    if(_origin && _twin && _twin->_origin) {
        stream << "sb-\t" << _origin->get_position()[0] << "\t"
               << _twin->_origin->get_position()[0] << "\t"
               << _origin->get_position()[1] << "\t"
               << _twin->_origin->get_position()[1] << "\n";
    } else {
        if(_origin) {
            stream << "pbo\t" << _origin->get_position()[0] << "\t"
                   << _origin->get_position()[1] << "\n";
        }
    }
}

/**
 * @brief Write the edge in gnuplot plottable format to the given stream
 *
 * @param stream std::ostream to write to
 */
void Edge::plot_gnuplot(ostream& stream) {
    if(_origin && _twin && _twin->_origin) {
        stream << _origin->get_position()[0] << "\t"
               << _origin->get_position()[1] << "\n";
        stream << _twin->_origin->get_position()[0] << "\t"
               << _twin->_origin->get_position()[1] << "\n\n";
    }
}

/**
 * @brief Cut off the edge at the border of the given box
 *
 * @warning This function does not work (yet?)
 *
 * @param container DelCont specifying the dimensions of the simulation box
 * @param vertices Reference to the list of vertices of the tesselation
 * @param edges Reference to the list of edges of the tesselation
 * @return True if both the edge and its twin have an origin within the
 * simulation box, false otherwise, although not completely sure...
 */
bool Edge::trim_edge(DelCont* container, vector<Vertex*>& vertices,
                     vector<Edge*>& edges) {
    if(_origin && _twin && _twin->_origin) {
        bool this_in = container->inside(_origin->get_position());
        bool twin_in = container->inside(_twin->_origin->get_position());
        if(this_in && twin_in) {
            return true;
        } else {
            if(!this_in && !twin_in) {
                return false;
            }
            if(!this_in) {
                double box[3];
                container->get_bounding_box(box);
                double fraction = 0.;
                if(_origin->get_position().x() < box[0]) {
                    fraction = std::max(
                            fraction,
                            fabs(_origin->get_position().x() - box[0]) /
                                    fabs(_origin->get_position().x() -
                                         _twin->get_origin().x()));
                }
                if(_origin->get_position().x() > box[0] + box[2]) {
                    fraction = std::max(
                            fraction, fabs(_origin->get_position().x() -
                                           box[0] - box[2]) /
                                              fabs(_origin->get_position().x() -
                                                   _twin->get_origin().x()));
                }
                if(_origin->get_position().y() < box[1]) {
                    fraction = std::max(
                            fraction,
                            fabs(_origin->get_position().y() - box[1]) /
                                    fabs(_origin->get_position().y() -
                                         _twin->get_origin().y()));
                }
                if(_origin->get_position().y() > box[1] + box[2]) {
                    fraction = std::max(
                            fraction, fabs(_origin->get_position().y() -
                                           box[1] - box[2]) /
                                              fabs(_origin->get_position().y() -
                                                   _twin->get_origin().y()));
                }
                fraction = 1. - fraction;
                Vec origin(_twin->_origin->get_position()[0] +
                                   (_origin->get_position()[0] -
                                    _twin->_origin->get_position()[0]) *
                                           fraction,
                           _twin->_origin->get_position()[1] +
                                   (_origin->get_position()[1] -
                                    _twin->_origin->get_position()[1]) *
                                           fraction);
                if(_origin->dereference()) {
                    delete _origin;
                }
                _origin = new Vertex(origin);
                _origin->reference();
                vertices.push_back(_origin);
                edges.push_back(_twin);
            } else {
                return _twin->trim_edge(container, vertices, edges);
            }
            return true;
        }
    } else {
        return false;
    }
}

/**
 * @brief Check if this edge lies completely outside the given simulation box
 *
 * @param container DelCont specifying the dimensions of the simulation box
 * @return True if the origin of the edge and of its twin are outside the
 * simulation box
 */
bool Edge::outside(DelCont* container) {
    bool this_in = container->inside(_origin->get_position());
    bool twin_in = container->inside(_twin->_origin->get_position());
    return !this_in && !twin_in;
}

/**
 * @brief Calculate the area and the centroid of the cell
 */
void Face::calculate_quantities() {
    Edge* current = _edge;
    Edge* next = _edge->get_next_edge();
    while(next != _edge) {
        Vec pi = current->get_origin();
        Vec pj = next->get_origin();
        double fac = pi.x() * pj.y() - pj.x() * pi.y();
        _area += fac;
        _centroid += fac * (pi + pj);
        current = next;
        next = current->get_next_edge();
    }
    Vec pi = current->get_origin();
    Vec pj = next->get_origin();
    double fac = pi.x() * pj.y() - pj.x() * pi.y();
    _area += fac;
    _centroid += fac * (pi + pj);
    _area *= 0.5;
    _centroid /= (6. * _area);
}

/**
 * @brief Destructor
 *
 * Clean up the edges and the cells of the Voronoi tesselation. The vertices are
 * deleted by the edges using a custom smart pointer method.
 */
DoublyConnectedEdgeList::~DoublyConnectedEdgeList() {
    for(unsigned int i = 0; i < _edges.size(); i++) {
        delete _edges[i];
    }
    for(unsigned int i = 0; i < _faces.size(); i++) {
        delete _faces[i];
    }
}

/**
 * @brief Write the tesselation in plottable format to the given stream
 *
 * @param stream std::ostream to write to
 */
void DoublyConnectedEdgeList::plot(ostream& stream) {
    for(unsigned int i = 0; i < _edges.size(); i++) {
        _edges[i]->plot(stream);
    }
}

/**
 * @brief Write the tesselation in gnuplot plottable format to the given stream
 *
 * @param stream std::ostream to write to
 */
void DoublyConnectedEdgeList::plot_gnuplot(ostream& stream) {
    for(unsigned int i = 0; i < _edges.size(); i++) {
        _edges[i]->plot_gnuplot(stream);
    }
}

/**
 * @brief Add the given Edge to the tesselation
 *
 * @param edge Edge to add to the tesselation
 */
void DoublyConnectedEdgeList::add_edge(Edge* edge) {
    _edges.push_back(edge);
}

/**
 * @brief Trim the edges of the tesselation to complete cells at the border
 * of the simulation box
 *
 * @warning This part of the algorithm does not work (yet?)
 *
 * @param container DelCont specifying the dimensions of the simulation box
 */
void DoublyConnectedEdgeList::trim_edges(DelCont* container) {
    vector<Edge*> remaining_edges;
    vector<Vertex*> border;
    vector<Edge*> connections;
    cout << _edges.size() << endl;
    for(unsigned int i = 0; i < _edges.size(); i++) {
        bool trim = _edges[i]->trim_edge(container, border, connections);
        if(trim) {
            remaining_edges.push_back(_edges[i]);
        } else {
            if(_edges[i]->get_twin()) {
                _edges[i]->get_twin()->set_twin(NULL);
            }
            delete _edges[i];
        }
    }
    _edges = remaining_edges;
    //    cout << border.size() << endl;

    double box[3];
    container->get_bounding_box(box);
    Vec origin(box[0], box[1]);
    Vertex* origin_vertex = new Vertex(origin);
    border.push_back(origin_vertex);
    origin.set(box[0], box[1] + box[2]);
    origin_vertex = new Vertex(origin);
    border.push_back(origin_vertex);
    origin.set(box[0] + box[2], box[1]);
    origin_vertex = new Vertex(origin);
    border.push_back(origin_vertex);
    origin.set(box[0] + box[2], box[1] + box[2]);
    origin_vertex = new Vertex(origin);
    border.push_back(origin_vertex);

    VertexSorter sorter(border);
    vector<unsigned int> indices(border.size(), 0);
    for(unsigned int i = 0; i < border.size(); i++) {
        indices[i] = i;
    }
    sort(indices.begin(), indices.end(), sorter);

    Edge* prev_added = NULL;
    cout << indices.size() << endl;
    for(unsigned int i = 0; i < indices.size(); i++) {
        unsigned int index1 = indices[i];
        unsigned int index2 = indices[0];
        if(i < indices.size() - 1) {
            index2 = indices[i + 1];
        }
        Edge* edge1 = new Edge();
        Edge* edge2 = new Edge();
        edge1->set_origin(border[index1]);
        edge2->set_origin(border[index2]);
        edge1->set_twin(edge2);
        edge2->set_twin(edge1);
        _edges.push_back(edge1);
        _edges.push_back(edge2);

        // if the index represents a corner of the box, there is no associated
        // edge connection
        if(index1 < connections.size()) {
            connections[index1]->set_next_edge(edge1);
            edge1->set_prev_edge(connections[index1]);
            // connections[index1]->get_twin()->set_prev_edge is the previously
            // added edge
            // BUT connections[index1]->get_twin was
            // connections[index2]->get_twin in the previous round
        } else {
            edge1->set_prev_edge(prev_added);
            if(prev_added) {
                prev_added->set_next_edge(edge1);
            }
        }
        if(index2 < connections.size()) {
            connections[index2]->get_twin()->set_prev_edge(edge1);
            edge1->set_next_edge(connections[index2]->get_twin());
            // connections[index2]->set_next_edge is the next added edge
            // BUT connections[index2] in the previous round is now
            // connections[index1] :)

        } else {
            // edge1->set_next_edge is the next edge that is added
            // this is already handled by setting the next edge for prev_added
        }
        // edge2->set_prev_edge is the next edge that is added
        // edge2->set_next_edge is the previously added edge
        // BUT maybe we don't want to set this, because the outer area is not
        // really a cell (is it?)
        //        if(prev_added){
        //            edge2->set_next_edge(prev_added->get_twin());
        //            prev_added->get_twin()->set_prev_edge(edge2);
        //        }

        prev_added = edge1;
    }
}

/**
 * @brief Print the cell corresponding to the given Edge in gnuplot plottable
 * format to the given stream
 *
 * Different edges correspond to the same cell. The cell is reconstructed by
 * iteratively plotting the next edge until the next edge is equal to the
 * current edge.
 *
 * @param stream std::ostream to write to
 * @param edge Edge used to identify the cell
 */
void DoublyConnectedEdgeList::print_cell(ostream& stream, Edge* edge) {
    edge->plot_gnuplot(stream);
    Edge* next = edge->get_next_edge();
    while(next && next != edge) {
        next->plot_gnuplot(stream);
        next = next->get_next_edge();
    }
    //    if(!next){
    //        next = edge->get_prev_edge();
    //        while(next){
    //            next->plot_gnuplot(stream);
    //            next = next->get_prev_edge();
    //        }
    //    }
}

/**
 * @brief Print some cells as a test of the algorithm
 *
 * @param stream std::ostream to write to
 */
void DoublyConnectedEdgeList::test_connections(ostream& stream) {
    print_cell(stream, _edges[30]);
    print_cell(stream, _edges[4]);
    print_cell(stream, _edges[50]);
}

/**
 * @brief Set the cells of the tesselation and calculate their properties
 *
 * We iterate over all edges and associate a new cell to every edge that does
 * not yet have a cell. The current edge is then used as reference edge for the
 * cell. All other edges that correspond to the same cell receive a pointer to
 * the cell.
 *
 * When the cell is complete, its geometrical properties are calculated.
 *
 * As a check of the algorithm, this method prints the total area of all cells
 * to the stdout. If everything works, this area should be equal to the area of
 * the simulation box.
 */
void DoublyConnectedEdgeList::set_faces() {
    double area = 0.;
    for(unsigned int i = 0; i < _edges.size(); i++) {
        if(!_edges[i]->get_face()) {
            if(!_edges[i]->get_next_edge()) {
                continue;
            }
            Face* face = new Face(_edges[i]);
            _edges[i]->set_face(face);
            Edge* current = _edges[i]->get_next_edge();
            while(current != _edges[i]) {
                current->set_face(face);
                current = current->get_next_edge();
            }
            _faces.push_back(face);
            face->calculate_quantities();
            area += face->get_area();
        }
    }
    cout << _faces.size() << " faces, total area: " << area << endl;
}

#endif
