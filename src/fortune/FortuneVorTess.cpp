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
 * @file FortuneVorTess.cpp
 *
 * @brief Voronoi tesselation created using Fortune's algorithm: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#if ndim_ == 2

#include "FortuneVorTess.hpp"
#include "../utilities/GasParticle.hpp"  // for GasParticle
#include "BinarySearchTree.hpp"          // for FortuneSite, BinaryLeaf, etc
#include "PriorityQueue.hpp"             // for CircleEvent, PriorityQueue, etc
#include "Vec.hpp"                       // for Vec
#include <cmath>                         // for fabs
#include <cstddef>                       // for NULL
#include <ostream>                       // for operator<<, ostream, etc
using namespace std;

/**
 * @brief Constructor
 *
 * @param container DelCont specifying the dimensions of the simulation box
 * @param npart Number of generators that will be added to the tesselation
 */
FortuneVorTess::FortuneVorTess(DelCont* container, unsigned int npart) {
    _container = container;
    if(npart) {
        _sites.reserve(npart);
    }
}

/**
 * @brief Destructor
 *
 * Clean up the generators
 */
FortuneVorTess::~FortuneVorTess() {
    for(unsigned int i = 0; i < _sites.size(); i++) {
        delete _sites[i];
    }
}

/**
 * @brief Add the given GasParticle to the tesselation
 *
 * @param particle GasParticle to add to the tesselation
 */
void FortuneVorTess::add_point(GasParticle* particle) {
    _sites.push_back(new FortuneSite(particle->get_position()));
}

/**
 * @brief Construct the Voronoi tesselation
 */
void FortuneVorTess::construct() {
    BinarySearchTree stree;
    PriorityQueue queue;

    for(unsigned int i = 0; i < _sites.size(); i++) {
        queue.add_event(
                new SiteEvent(_sites[i]->get_position().y(), _sites[i]));
    }

    Event* event = queue.get_event();
    while(event) {
        //        if(event->get_priority() < 0.){
        //            break;
        //        }
        if(event->isCircleEvent()) {
            // do circle event stuff
            //            cout << "CircleEvent: "
            //            << ((CircleEvent*)event)->get_arc()
            //            ->get_site()->get_position().x() << "\t"
            //            << ((CircleEvent*)event)
            //            ->get_arc()->get_site()->get_position().y() << endl;
            CircleEvent* circle_event = (CircleEvent*)event;
            CircleEvent* events[2];
            Edge* edges[3];
            stree.remove(circle_event->get_arc(), circle_event->get_priority(),
                         events, edges);
            if(events[0]) {
                queue.add_event(events[0]);
            }
            if(events[1]) {
                queue.add_event(events[1]);
            }

            Vertex* eventpoint = new Vertex(circle_event->get_midpoint());
            for(int i = 0; i < 3; i++) {
                if(edges[i]) {
                    edges[i]->set_origin(eventpoint);
                    //                    _vorlist.add_edge(edges[i]);
                }
            }
            //            if(edges[0]){
            //                _vorlist.add_edge(edges[0]);
            //            }
            if(edges[1]) {
                _vorlist.add_edge(edges[1]);
            }
            if(edges[2]) {
                _vorlist.add_edge(edges[2]);
            }

            //            handleCircleEvent((CircleEvent*)event, stree, queue);
        } else {
            //            cout << "SiteEvent: " << ((SiteEvent*)event)
            //            ->get_site()->get_position().x() << "\t" <<
            // ((SiteEvent*)event)
            //            ->get_site()->get_position().y() << endl;
            //            CircleEvent* events[2];
            //            FortuneSite* site = ((SiteEvent*)event)->get_site();
            //            stree.add(site, site->get_position().y(), events);
            //            if(events[0]){
            //                queue.add_event(events[0]);
            //            }
            //            if(events[1]){
            //                queue.add_event(events[1]);
            //            }

            handleSiteEvent((SiteEvent*)event, stree, queue);
        }
        delete event;
        event = queue.get_event();
    }
    stree.get_halfedges(_vorlist, _container);
    _vorlist.trim_edges(_container);

    _vorlist.set_faces();
}

/**
 * @brief Check if the edges in between the given BinaryLeafs converge and give
 * rise to a new CircleEvent
 *
 * @param left Left BinaryLeaf
 * @param middle Middle BinaryLeaf
 * @param right Right BinaryLeaf
 * @param ybeach y-coordinate of the current beachline
 * @return Pointer to a newly created CircleEvent (NULL if the edges diverge)
 */
CircleEvent* FortuneVorTess::check_triple(BinaryLeaf* left, BinaryLeaf* middle,
                                          BinaryLeaf* right, double ybeach) {
    if(!(left->get_site()->get_position().x() ==
                 middle->get_site()->get_position().x() &&
         middle->get_site()->get_position().x() ==
                 right->get_site()->get_position().x())) {
        Vec midpoint;
        double radius = BinaryInternalNode::get_circle(
                left->get_site()->get_position(),
                middle->get_site()->get_position(),
                right->get_site()->get_position(), midpoint);
        double distance = ybeach - midpoint[1] + radius;
        if(radius && distance >= 0.) {
            double d1 =
                    fabs(BinaryInternalNode::intersection(
                                 left->get_site()->get_position(),
                                 middle->get_site()->get_position(), ybeach) -
                         midpoint[0]);
            double d2 = fabs(BinaryInternalNode::intersection(
                                     left->get_site()->get_position(),
                                     middle->get_site()->get_position(),
                                     ybeach - 0.1 * distance) -
                             midpoint[0]);
            double d3 =
                    fabs(BinaryInternalNode::intersection(
                                 middle->get_site()->get_position(),
                                 right->get_site()->get_position(), ybeach) -
                         midpoint[0]);
            double d4 = fabs(BinaryInternalNode::intersection(
                                     middle->get_site()->get_position(),
                                     right->get_site()->get_position(),
                                     ybeach - 0.1 * distance) -
                             midpoint[0]);
            if(d2 < d1 && d4 < d3) {
                CircleEvent* event =
                        new CircleEvent(midpoint[1] - radius, midpoint, middle);
                middle->set_event(event);
                return event;
            }
        }
    }
    return NULL;
}

/**
 * @brief Handle a SiteEvent by adding the new site to the tesselation and by
 * creating a new arc for the site that is added to the beachline
 *
 * @param event SiteEvent that is being handled
 * @param stree BinarySearchTree that represents the beachline
 * @param queue PriorityQueue to add newly created events to
 */
void FortuneVorTess::handleSiteEvent(SiteEvent* event, BinarySearchTree& stree,
                                     PriorityQueue& queue) {
    CircleEvent* events[2];
    BinaryLeaf* new_arc = stree.add(
            event->get_site(), event->get_site()->get_position().y(), events);
    BinaryLeaf* old_arc_left = new_arc->get_prev_leaf();
    if(old_arc_left) {
        BinaryLeaf* arc_left = old_arc_left->get_prev_leaf();
        if(arc_left) {
            // check event
            CircleEvent* circle_event =
                    check_triple(arc_left, old_arc_left, new_arc,
                                 event->get_site()->get_position().y());
            if(circle_event) {
                queue.add_event(circle_event);
            }
        }
    }
    BinaryLeaf* old_arc_right = new_arc->get_next_leaf();
    if(old_arc_right) {
        BinaryLeaf* arc_right = old_arc_right->get_next_leaf();
        if(arc_right) {
            // check event
            CircleEvent* circle_event =
                    check_triple(new_arc, old_arc_right, arc_right,
                                 event->get_site()->get_position().y());
            if(circle_event) {
                queue.add_event(circle_event);
            }
        }
    }
    // check for events and add them to the queue
    //    if(events[0]){
    //        queue.add_event(events[0]);
    //    }
    //    if(events[1]){
    //        queue.add_event(events[1]);
    //    }

    //    BinaryInternalNode* parent_left = new_arc->get_left_parent();
    //    BinaryInternalNode* parent_right = new_arc->get_right_parent();
    //    // do edges
}

/**
 * @brief Handle a CircleEvent by creating a new vertex in the tesselation
 *
 * @deprecated This method does not do anything...
 *
 * @param event CircleEvent that is being handled
 * @param stree BinarySearchTree that represents the beachline
 * @param queue PriorityQueue to add newly created events to
 */
void FortuneVorTess::handleCircleEvent(CircleEvent* event,
                                       BinarySearchTree& stree,
                                       PriorityQueue& queue) {
    //    BinaryLeaf* middle_arc = event->get_arc();
    //    BinaryLeaf* left_middle_arc = middle_arc->get_prev_leaf();
    //    BinaryLeaf* left_arc = left_middle_arc->get_prev_leaf();
    //    BinaryLeaf* right_middle_arc = middle_arc->get_next_leaf();
    //    BinaryLeaf* right_arc = right_middle_arc->get_next_leaf();
    //    // check for events and add them to the queue

    //    BinaryInternalNode* parent_left = middle_arc->get_left_parent();
    //    BinaryInternalNode* parent_right = middle_arc->get_right_parent();
    //    // do edge stuff

    //    // now we can remove the arc

    //    // we have to find out which parent was removed and which one gets a
    // new
    //    // edge...
    //    // or we can return it...
}

/**
 * @brief Write the tesselation in gnuplot-plottable format to the given stream
 *
 * @param stream std::ostream to write to
 */
void FortuneVorTess::print_tesselation_gnuplot(ostream& stream) {
    for(unsigned int i = 0; i < _sites.size(); i++) {
        stream << _sites[i]->get_position().x() << "\t"
               << _sites[i]->get_position().y() << "\n\n";
    }
    _vorlist.plot_gnuplot(stream);
}

/**
 * @brief Write the tesselation in plottable format to the given stream
 *
 * @param stream std::ostream to write to
 */
void FortuneVorTess::print_tesselation(ostream& stream) {
    for(unsigned int i = 0; i < _sites.size(); i++) {
        stream << "pro\t" << _sites[i]->get_position().x() << "\t"
               << _sites[i]->get_position().y() << "\n";
    }
    _vorlist.plot(stream);
}

/**
 * @brief Test the connections of the tesselation
 *
 * @param stream std::ostream to write output to
 */
void FortuneVorTess::test_connections(ostream& stream) {
    _vorlist.test_connections(stream);
}

#endif
