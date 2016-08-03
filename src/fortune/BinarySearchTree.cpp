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
 * @file BinarySearchTree.cpp
 *
 * @brief Binary search tree for Fortune's algorithm: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#if ndim_ == 2

#include "BinarySearchTree.hpp"
#include "DelCont.hpp"                  // for DelCont
#include "DoublyConnectedEdgeList.hpp"  // for Edge, Vertex, etc
#include "PriorityQueue.hpp"            // for CircleEvent
#include <algorithm>                    // for max
#include <cstdlib>                      // for NULL, abs
using namespace std;

// BinaryInternalNode

/**
 * @brief Constructor
 *
 * @param pleft Left generator of the Edge
 * @param pright Right generator of the Edge
 * @param parent Parent node
 */
BinaryInternalNode::BinaryInternalNode(FortuneSite* pleft, FortuneSite* pright,
                                       BinaryInternalNode* parent) {
    _pleft = pleft;
    _pright = pright;
    _parent = parent;
    _leaf = NULL;
    _children[0] = NULL;
    _children[1] = NULL;
    _edge = NULL;
}

/**
 * @brief Destructor. Clean up children.
 */
BinaryInternalNode::~BinaryInternalNode() {
    // if one of these is NULL, deleting it does not cause problems
    delete _leaf;
    delete _children[0];
    delete _children[1];
}

/**
 * @brief Erase the information stored in this node by setting all children to
 * NULL
 */
void BinaryInternalNode::erase() {
    _leaf = NULL;
    _children[0] = NULL;
    _children[1] = NULL;
}

/**
 * @brief Check if the tree is balanced and rebalance it if not
 *
 * @param n Maximum allowed depth difference between children of this node
 */
void BinaryInternalNode::check_balance(int n) {
    if(abs(_children[0]->depth() - _children[1]->depth()) > n) {
        rebalance();
    } else {
        if(_parent) {
            _parent->check_balance(1);
        }
    }
}

/**
 * @brief Calculated the depth of the node
 *
 * The depth is the maximal number of nodes between a child of this node and a
 * leaf of the tree.
 *
 * A leaf has depth 0, the parent of that leaf has depth 1 and so on until we
 * are at the root of the tree.
 *
 * @return The depth of this node
 */
int BinaryInternalNode::depth() {
    if(_leaf) {
        return 0;
    } else {
        return std::max(_children[0]->depth(), _children[1]->depth()) + 1;
    }
}

/**
 * @brief Restructure the nodes of this node so that the depth difference
 * between its children is at most 1
 *
 * There are four different scenarios for an unbalanced node, which are
 * distinguished by the depth of the children and their children. Once we know
 * which scenario applies, rebalancing is just a matter of swapping some
 * internal pointers.
 */
void BinaryInternalNode::rebalance() {
    BinaryInternalNode* T[4] = {NULL};
    BinaryInternalNode* subnodes[2] = {NULL};
    FortuneSite* points[6] = {NULL};
    Edge* edges[3] = {NULL};
    BinaryInternalNode* children[2] = {NULL};
    if(_children[0]->depth() < _children[1]->depth()) {
        T[0] = _children[0];
        points[0] = _pleft;
        points[1] = _pright;
        edges[0] = _edge;
        subnodes[0] = _children[1];
        children[0] = _children[1]->get_child(0);
        children[1] = _children[1]->get_child(1);
        if(children[0]->depth() < children[1]->depth()) {
            T[1] = children[0];
            T[2] = children[1]->get_child(0);
            T[3] = children[1]->get_child(1);
            points[2] = _children[1]->get_pleft();
            points[3] = _children[1]->get_pright();
            points[4] = children[1]->get_pleft();
            points[5] = children[1]->get_pright();
            edges[1] = _children[1]->get_edge();
            edges[2] = children[1]->get_edge();
            subnodes[1] = children[1];
        } else {
            T[1] = children[0]->get_child(0);
            T[2] = children[0]->get_child(1);
            T[3] = children[1];
            points[2] = children[0]->get_pleft();
            points[3] = children[0]->get_pright();
            points[4] = _children[1]->get_pleft();
            points[5] = _children[1]->get_pright();
            edges[1] = children[0]->get_edge();
            edges[2] = _children[1]->get_edge();
            subnodes[1] = children[0];
        }
    } else {
        T[3] = _children[1];
        points[4] = _pleft;
        points[5] = _pright;
        edges[2] = _edge;
        subnodes[0] = _children[0];
        children[0] = _children[0]->get_child(0);
        children[1] = _children[0]->get_child(1);
        if(children[0]->depth() < children[1]->depth()) {
            T[0] = children[0];
            T[1] = children[1]->get_child(0);
            T[2] = children[1]->get_child(1);
            points[0] = _children[0]->get_pleft();
            points[1] = _children[0]->get_pright();
            points[2] = children[1]->get_pleft();
            points[3] = children[1]->get_pright();
            edges[0] = _children[0]->get_edge();
            edges[1] = children[1]->get_edge();
            subnodes[1] = children[1];
        } else {
            T[0] = children[0]->get_child(0);
            T[1] = children[0]->get_child(1);
            T[2] = children[1];
            points[0] = children[0]->get_pleft();
            points[1] = children[0]->get_pright();
            points[2] = _children[0]->get_pleft();
            points[3] = _children[0]->get_pright();
            edges[0] = children[0]->get_edge();
            edges[1] = _children[0]->get_edge();
            subnodes[1] = children[0];
        }
    }

    _pleft = points[2];
    _pright = points[3];
    _children[0] = subnodes[0];
    _children[1] = subnodes[1];
    _edge = edges[1];

    subnodes[0]->set_pleft(points[0]);
    subnodes[0]->set_pright(points[1]);
    subnodes[0]->set_child(T[0], 0);
    subnodes[0]->set_child(T[1], 1);
    subnodes[0]->set_parent(this);
    subnodes[0]->set_edge(edges[0]);
    T[0]->set_parent(subnodes[0]);
    T[1]->set_parent(subnodes[0]);

    subnodes[1]->set_pleft(points[4]);
    subnodes[1]->set_pright(points[5]);
    subnodes[1]->set_child(T[2], 0);
    subnodes[1]->set_child(T[3], 1);
    subnodes[1]->set_parent(this);
    subnodes[1]->set_edge(edges[2]);
    T[2]->set_parent(subnodes[1]);
    T[3]->set_parent(subnodes[1]);
}

/**
 * @brief Get the child node of this node with the given index
 *
 * @param index Index of the requested child (0 or 1)
 * @return Child node of this node
 */
BinaryInternalNode* BinaryInternalNode::get_child(unsigned int index) {
    return _children[index];
}

/**
 * @brief Get the left generator of the Edge associated with this node
 *
 * @return Left generator
 */
FortuneSite* BinaryInternalNode::get_pleft() {
    return _pleft;
}

/**
 * @brief Get the right generator of the Edge associated with this node
 *
 * @return Right generator
 */
FortuneSite* BinaryInternalNode::get_pright() {
    return _pright;
}

/**
 * @brief Get the Edge associated with this node
 *
 * @return Edge associated with this node
 */
Edge* BinaryInternalNode::get_edge() {
    return _edge;
}

/**
 * @brief Set the left generator of the Edge associated with this node
 *
 * @param pleft Left generator
 */
void BinaryInternalNode::set_pleft(FortuneSite* pleft) {
    _pleft = pleft;
}

/**
 * @brief Set the right generator of the Edge associated with this node
 *
 * @param pright Right generator
 */
void BinaryInternalNode::set_pright(FortuneSite* pright) {
    _pright = pright;
}

/**
 * @brief Set the child node of this node with the given index
 *
 * @param child Child node for this node
 * @param index Index of the child node (0 or 1)
 */
void BinaryInternalNode::set_child(BinaryInternalNode* child,
                                   unsigned int index) {
    _children[index] = child;
}

/**
 * @brief Set the parent node of this node
 *
 * @param parent Parent node for this node
 */
void BinaryInternalNode::set_parent(BinaryInternalNode* parent) {
    _parent = parent;
    if(_leaf) {
        _leaf->set_parent(parent);
    }
}

/**
 * @brief Set the Edge associated with this node
 *
 * @param edge Edge associated with this node
 */
void BinaryInternalNode::set_edge(Edge* edge) {
    _edge = edge;
}

/**
 * @brief Recursively check the balance of the tree by first checking the
 * balance at this level and then continuing on the level of the parent
 */
void BinaryInternalNode::check_balance_recursive() {
    check_balance(1);
    if(_parent) {
        _parent->check_balance_recursive();
    }
}

/**
 * @brief Add the given FortuneSite to the tree using the given beachline
 *
 * If this node is an internal node, we check in which of its two children we
 * should add the FortuneSite and recursively go down the tree.
 *
 * If we encounter a leaf, we replace it by a new internal node with the old
 * leaf and a new leaf containing the FortuneSite as its children. Once this is
 * done, we recursively check the balance of the tree and rebalance if
 * necessary.
 *
 * @param site FortuneSite to add to the tree
 * @param ybeach y-coordinate of the current beachline
 * @param events Array to store newly generated CircleEvent instances in
 * @return New leaf holding the FortuneSite
 */
BinaryLeaf* BinaryInternalNode::add(FortuneSite* site, double ybeach,
                                    CircleEvent** events) {
    // do a lot of stuff and so
    double xbreak = intersection(_pleft->get_position(),
                                 _pright->get_position(), ybeach);

    if(_pleft->get_position().y() == _pright->get_position().y() &&
       _pright->get_position().y() == ybeach) {
        if(_pleft->get_position().x() > _pright->get_position().x()) {
            // make sure we end up on the left
            xbreak = 1.e99;
        } else {
            xbreak = _pright->get_position().x();
        }
    }

    unsigned int index = 1;
    if(site->get_position().x() < xbreak) {
        index = 0;
    }

    if(_children[index]->get_leaf()) {
        // we found the arc!

        BinaryLeaf* arc = _children[index]->get_leaf();
        BinaryLeaf* ngbleft = arc->get_prev_leaf();
        BinaryLeaf* ngbright = arc->get_next_leaf();

        arc->deactivate();

        if(_pleft->get_position().y() == _pright->get_position().y() &&
           _pright->get_position().y() == ybeach) {
            // degenerate case at the start of the grid construction
            FortuneSite* site2 = arc->get_site();
            _children[index]->erase();
            delete _children[index];
            delete arc;
            BinaryLeaf* leaf_left;
            BinaryLeaf* leaf_right;
            BinaryLeaf* new_leaf;
            if(site->get_position().x() < site2->get_position().x()) {
                _children[index] = new BinaryInternalNode(site, site2, this);
                leaf_left = new BinaryLeaf(site, _children[index]);
                leaf_right = new BinaryLeaf(site2, _children[index]);
                new_leaf = leaf_left;
                if(index) {
                    _pright = site;
                } else {
                    // find first left parent
                    BinaryInternalNode* current = this;
                    BinaryInternalNode* parent = current->get_parent();
                    while(parent && parent->get_child(0) == current) {
                        current = parent;
                        parent = current->get_parent();
                    }
                    if(parent) {
                        parent->set_pright(site);
                    }
                }
            } else {
                _children[index] = new BinaryInternalNode(site2, site, this);
                leaf_left = new BinaryLeaf(site2, _children[index]);
                leaf_right = new BinaryLeaf(site, _children[index]);
                new_leaf = leaf_right;
                if(index) {
                    BinaryInternalNode* current = this;
                    BinaryInternalNode* parent = current->get_parent();
                    while(parent && parent->get_child(1) == current) {
                        current = parent;
                        parent = current->get_parent();
                    }
                    if(parent) {
                        parent->set_pleft(site);
                    }
                } else {
                    _pleft = site;
                }
            }
            _children[index]->add_leaf(leaf_left, 0);
            _children[index]->add_leaf(leaf_right, 1);

            _children[index]->set_edge(new Edge());
            _children[index]->get_edge()->set_twin(new Edge());
            _children[index]->get_edge()->get_twin()->set_twin(
                    _children[index]->get_edge());
            Vec origin(0.5 * (site2->get_position().x() +
                              site->get_position().x()),
                       100.);
            Vertex* originpoint = new Vertex(origin);
            _children[index]->get_edge()->get_twin()->set_origin(originpoint);
            // since the twin edge is not connected to any arc, we have created
            // a memory leak here (the Edge will never be added to the DCEL)
            // we should solve it, but for the moment, we rather choose to
            // ignore it

            if(ngbleft) {
                ngbleft->set_next_leaf(leaf_left);
            }
            leaf_left->set_next_leaf(leaf_right);
            if(ngbright) {
                leaf_right->set_next_leaf(ngbright);
            }

            check_balance(2);
            check_balance(1);
            return new_leaf;
        }

        FortuneSite* site2 = arc->get_site();
        _children[index]->erase();
        delete _children[index];
        delete arc;
        BinaryLeaf* leaves_left[2];
        BinaryLeaf* leaf_right;
        _children[index] = new BinaryInternalNode(site, site2, this);
        leaf_right = new BinaryLeaf(site2, _children[index]);
        _children[index]->add_subnode(site2, site, 0, leaves_left);
        _children[index]->add_leaf(leaf_right, 1);

        _children[index]->set_edge(new Edge());
        _children[index]->get_edge()->set_twin(
                _children[index]->get_child(0)->get_edge());
        _children[index]->get_child(0)->get_edge()->set_twin(
                _children[index]->get_edge());

        if(ngbleft) {
            ngbleft->set_next_leaf(leaves_left[0]);
        }
        leaves_left[0]->set_next_leaf(leaves_left[1]);
        leaves_left[1]->set_next_leaf(leaf_right);
        if(ngbright) {
            leaf_right->set_next_leaf(ngbright);
        }

        if(ngbleft &&
           !(ngbleft->get_site()->get_position().x() ==
                     site->get_position().x() &&
             site->get_position().x() == site2->get_position().x())) {
            Vec midpoint;
            double radius = BinaryInternalNode::get_circle(
                    ngbleft->get_site()->get_position(), site2->get_position(),
                    site->get_position(), midpoint);
            double distance = ybeach - midpoint[1] + radius;
            if(radius && distance >= 0.) {
                double d1 =
                        fabs(intersection(ngbleft->get_site()->get_position(),
                                          site2->get_position(), ybeach) -
                             midpoint[0]);
                double d2 =
                        fabs(intersection(ngbleft->get_site()->get_position(),
                                          site2->get_position(),
                                          ybeach - 0.1 * distance) -
                             midpoint[0]);
                double d3 = fabs(site->get_position().x() - midpoint[0]);
                double d4 = fabs(intersection(site2->get_position(),
                                              site->get_position(),
                                              ybeach - 0.1 * distance) -
                                 midpoint[0]);
                if(d2 < d1 && d4 < d3) {
                    events[0] = new CircleEvent(midpoint[1] - radius, midpoint,
                                                leaves_left[0]);
                    leaves_left[0]->set_event(events[0]);
                }
            }
        }
        if(ngbright &&
           !(ngbright->get_site()->get_position().x() ==
                     site->get_position().x() &&
             site->get_position().x() == site2->get_position().x())) {
            Vec midpoint;
            double radius = BinaryInternalNode::get_circle(
                    site->get_position(), site2->get_position(),
                    ngbright->get_site()->get_position(), midpoint);
            double distance = ybeach - midpoint[1] + radius;
            if(radius && distance >= 0.) {
                double d1 = fabs(intersection(site->get_position(),
                                              site2->get_position(), ybeach) -
                                 midpoint[0]);
                double d2 = fabs(intersection(site->get_position(),
                                              site2->get_position(),
                                              ybeach - 0.1 * distance) -
                                 midpoint[0]);
                double d3 = fabs(site->get_position().x() - midpoint[0]);
                double d4 =
                        fabs(intersection(site2->get_position(),
                                          ngbright->get_site()->get_position(),
                                          ybeach - 0.1 * distance) -
                             midpoint[0]);
                if(d2 < d1 && d4 < d3) {
                    events[1] = new CircleEvent(midpoint[1] - radius, midpoint,
                                                leaf_right);
                    leaf_right->set_event(events[1]);
                }
            }
        }

        check_balance(2);
        check_balance(1);
        return leaves_left[1];
    } else {
        return _children[index]->add(site, ybeach, events);
    }
}

/**
 * @brief Make a leaf of this node of the tree
 *
 * @param leaf BinaryLeaf that is the only child of this node
 */
void BinaryInternalNode::make_leaf(BinaryLeaf* leaf) {
    _leaf = leaf;
}

/**
 * @brief Replace the leaf child at the given index by a new internal node
 *
 * @param pleft Left generator for the new internal node
 * @param pright Right generator for the new internal node
 * @param index Index of the child that should be replaced (0 or 1)
 * @param new_leaves Array to store newly created BinaryLeaf instances in
 */
void BinaryInternalNode::add_subnode(FortuneSite* pleft, FortuneSite* pright,
                                     unsigned int index,
                                     BinaryLeaf** new_leaves) {
    _children[index] = new BinaryInternalNode(pleft, pright, this);
    new_leaves[0] = new BinaryLeaf(pleft, _children[index]);
    _children[index]->add_leaf(new_leaves[0], 0);
    new_leaves[1] = new BinaryLeaf(pright, _children[index]);
    _children[index]->add_leaf(new_leaves[1], 1);
    _children[index]->set_edge(new Edge());
}

/**
 * @brief Add the given leaf as a child at the given index to this internal node
 *
 * @param leaf BinaryLeaf to add
 * @param index Index of the new child (0 or 1)
 */
void BinaryInternalNode::add_leaf(BinaryLeaf* leaf, unsigned int index) {
    _children[index] = new BinaryInternalNode(NULL, NULL, this);
    _children[index]->make_leaf(leaf);
}

/**
 * @brief Remove the given BinaryLeaf from the tree using the given beachline
 *
 * @param leaf BinaryLeaf to remove
 * @param ybeach y-coordinate of the current beachline
 * @param events Array to store newly created CircleEvent instances in
 * @param edges Array to store newly created Edge instances in
 */
void BinaryInternalNode::remove(BinaryLeaf* leaf, double ybeach,
                                CircleEvent** events, Edge** edges) {
    if(_children[0]->get_leaf()) {
        if(_children[0]->get_leaf() == leaf) {
            _parent->convert_child(this, _children[1], leaf, ybeach, events,
                                   edges);
            return;
        }
    }
    // let's suppose everything is correct... Then in this case, _children[0]
    // has to be a leaf and has to be the given leaf
    _parent->convert_child(this, _children[0], leaf, ybeach, events, edges);
}

/**
  * @brief Effectively removes a leaf in the BinarySearchTree by converting a
  * grandchild BinaryInternalNode of the current BinaryInternalNode to a leaf.
  *
  * This function is the main handler of a CircleEvent. Two existing edges A and
  * B converge to the midpoint of the CircleEvent (which is not involved here).
  * The disappearing arc gives rise to two new triples of consecutive arcs that
  *have to be checked for convergence and that can generate a new CircleEvent.
  * A new Edge originates from the midpoint of the CircleEvent and is connected
  *to the correct BinaryInternalNode.
  *
  * The edges A and B are determined as follows (A is called edges[1] and B
  * edges[2]):
  *  - If the disappearing arc is part of the left child:
  *    - If the disappearing arc is the left child of the left child:
  *        - A is the Edge connected to the first ancestor of the current
  *          BinaryInternalNode that is to the left of the current
  *          BinaryInternalNode
  *        - B is the Edge connected to the BinaryInternalNode that is converted
  *    - If the disappearing arc is the right child of the left child:
  *        - A is the Edge connected to the BinaryInternalNode that is converted
  *        - B is the Edge connected to the current BinaryInternalNode
  *  - If the disappearing arc is part of the right child:
  *    - If the disappearing arc is the left child of the right child:
  *        - A is the Edge connected to the current BinaryInternalNode
  *        - B is the Edge connected to the BinaryInternalNode that is converted
  *    - If the disappearing arc is the right child of the right child:
  *        - A is the Edge connected to the BinaryInternalNode that is converted
  *        - B is the Edge connected to the first ancestor of the current
  *          BinaryInternalNode that is to the right of the current
  *          BinaryInternalNode
  *
  * We can check what the first ancestor to the left or to the right of a
  * BinaryInternalNode is by recursively walking up the nodes: call the current
  * BinaryInternalNode arc1 and suppose that arc2 is its parent. As long as arc1
  * is the left/right child of arc2, arc2 becomes the new arc1 and the
  * parent of arc2 becomes the new arc2. If arc2 is NULL (and hence arc1 is the
  * root of the tree), there is no ancestor to the left/right (but this cannot
  * happen geometrically).
  *
  * A converging Edge may originate from
  *  - A SiteEvent, in which case it already has a twin Edge
  *  - A previous CircleEvent, in which case there is no twin yet
  * If there is no twin, a twin Edge is created and connected to the
  * corresponding Edge.
  *
  * All newly created edges get the midpoint of the CircleEvent as origin. This
  * means:
  *  - All edges following from a SiteEvent that converge on this event (and not
  *    their twins)
  *  - The Edge originating from this CircleEvent
  *  - All twin edges of edges following from a previous CircleEvent
  *
  * The new triples of arcs are:
  *  - The triple leaf->_prev_leaf->_prev_leaf, leaf->_prev_leaf,
  *    leaf->_next_leaf
  *  - The triple leaf->_prev_leaf, leaf->_next_leaf,
  *    leaf->_next_leaf->_next_leaf
  *
  * The newly created Edge is always connected to the BinaryInternalNode that
  * has a converging Edge and that does not disappear (so A or B, where the one
  * that is not the BinaryInternalNode that is converted is chosen).
  *
  * @param child Child BinaryInternalNode of the current BinaryInternalNode that
  * has to be converted to a leaf of the tree
  * @param grandchild Child of the child that is converted that is not changed
  * by the conversion (but rises a level in the tree)
  * @param leaf BinaryLeaf that is removed from the tree. This is a child of the
  * given child of the current BinaryInternalNode
  * @param ybeach Current position of the beachline (needed for convergence
  * check)
  * @param events Empty array with storage for two possible CircleEvent events
  * that can arise due to the two new triples of arcs
  * @param edges Empty array with storage for up to three edges that are created
  * or affected by the event
  */
void BinaryInternalNode::convert_child(BinaryInternalNode* child,
                                       BinaryInternalNode* grandchild,
                                       BinaryLeaf* leaf, double ybeach,
                                       CircleEvent** events, Edge** edges) {
    for(int i = 0; i < 3; i++) {
        edges[i] = NULL;
    }

    // maybe you should check if both ngbs exist. But no, wait a minute, this
    // does not make sense geometrically...
    leaf->get_prev_leaf()->set_next_leaf(leaf->get_next_leaf());

    BinaryLeaf* prev_leaf = leaf->get_prev_leaf()->get_prev_leaf();
    if(prev_leaf && prev_leaf->get_site() != leaf->get_site()) {
        Vec midpoint;
        Vec pleft = prev_leaf->get_site()->get_position();
        Vec pmid = leaf->get_prev_leaf()->get_site()->get_position();
        Vec pright = leaf->get_next_leaf()->get_site()->get_position();
        double radius = get_circle(pleft, pmid, pright, midpoint);
        double distance = ybeach - midpoint[1] + radius;
        if(radius && distance >= 1.e-13) {
            distance = fabs(distance);
            double d1 = fabs(intersection(pleft, pmid, ybeach) - midpoint[0]);
            double d2 =
                    fabs(intersection(pleft, pmid, ybeach - 0.1 * distance) -
                         midpoint[0]);
            double d3 = fabs(intersection(pmid, pright, ybeach) - midpoint[0]);
            double d4 =
                    fabs(intersection(pmid, pright, ybeach - 0.1 * distance) -
                         midpoint[0]);
            if(d2 - 1.e-13 < d1 && d4 - 1.e-13 < d3) {
                events[0] = new CircleEvent(midpoint[1] - radius, midpoint,
                                            leaf->get_prev_leaf());
                leaf->get_prev_leaf()->set_event(events[0]);
            }
        }
    }
    BinaryLeaf* next_leaf = leaf->get_next_leaf()->get_next_leaf();
    if(next_leaf && next_leaf->get_site() != leaf->get_site()) {
        Vec midpoint;
        Vec pleft = leaf->get_prev_leaf()->get_site()->get_position();
        Vec pmid = leaf->get_next_leaf()->get_site()->get_position();
        Vec pright = next_leaf->get_site()->get_position();
        double radius = get_circle(pleft, pmid, pright, midpoint);
        double distance = ybeach - midpoint[1] + radius;
        if(radius && distance >= -1.e-13) {
            distance = fabs(distance);
            double d1 = fabs(intersection(pleft, pmid, ybeach) - midpoint[0]);
            double d2 =
                    fabs(intersection(pleft, pmid, ybeach - 0.1 * distance) -
                         midpoint[0]);
            double d3 = fabs(intersection(pmid, pright, ybeach) - midpoint[0]);
            double d4 =
                    fabs(intersection(pmid, pright, ybeach - 0.1 * distance) -
                         midpoint[0]);
            if(d2 - 1.e-13 < d1 && d4 - 1.e-13 < d3) {
                events[1] = new CircleEvent(midpoint[1] - radius, midpoint,
                                            leaf->get_next_leaf());
                leaf->get_next_leaf()->set_event(events[1]);
            }
        }
    }

    if(_children[0] == child) {
        _children[0] = grandchild;
        if(child->get_child(0)->get_leaf() &&
           child->get_child(0)->get_leaf() == leaf) {
            child->get_child(0)->erase();
            delete child->get_child(0);
            // we have to find the parent that holds the wrong _pright and
            // change it
            BinaryInternalNode* current = this;
            BinaryInternalNode* parent = current->get_parent();
            while(parent && parent->get_child(0) == current) {
                current = parent;
                parent = current->get_parent();
            }
            parent->set_pright(leaf->get_next_leaf()->get_site());
            edges[0] = new Edge();
            edges[1] = parent->get_edge();
            edges[2] = child->get_edge();
            parent->set_edge(edges[0]);
        } else {
            child->get_child(1)->erase();
            delete child->get_child(1);
            _pleft = leaf->get_prev_leaf()->get_site();
            edges[0] = new Edge();
            edges[1] = child->get_edge();
            edges[2] = _edge;
            _edge = edges[0];
        }
    } else {
        _children[1] = grandchild;
        if(child->get_child(0)->get_leaf() &&
           child->get_child(0)->get_leaf() == leaf) {
            child->get_child(0)->erase();
            delete child->get_child(0);
            _pright = leaf->get_next_leaf()->get_site();
            edges[0] = new Edge();
            edges[1] = _edge;
            edges[2] = child->get_edge();
            _edge = edges[0];
        } else {
            child->get_child(1)->erase();
            delete child->get_child(1);
            BinaryInternalNode* current = this;
            BinaryInternalNode* parent = current->get_parent();
            while(parent && parent->get_child(1) == current) {
                current = parent;
                parent = current->get_parent();
            }
            parent->set_pleft(leaf->get_prev_leaf()->get_site());
            edges[0] = new Edge();
            edges[1] = child->get_edge();
            edges[2] = parent->get_edge();
            parent->set_edge(edges[0]);
        }
    }

    edges[0]->set_twin(new Edge());
    edges[0]->get_twin()->set_twin(edges[0]);
    edges[0] = edges[0]->get_twin();

    // edge connections
    edges[0]->set_prev_edge(edges[2]->get_twin());
    edges[1]->set_prev_edge(edges[0]->get_twin());
    edges[2]->set_prev_edge(edges[1]->get_twin());

    edges[0]->get_twin()->set_next_edge(edges[1]);
    edges[1]->get_twin()->set_next_edge(edges[2]);
    edges[2]->get_twin()->set_next_edge(edges[0]);

    child->erase();
    delete child;
    grandchild->set_parent(this);
    check_balance_recursive();
}

/**
 * @brief Get the parent of this node
 *
 * @return Parent of this node
 */
BinaryInternalNode* BinaryInternalNode::get_parent() {
    return _parent;
}

/**
 * @brief Write the node to the given stream in a format that can be plotted by
 * a specialized Python-script
 *
 * @param stream std::ostream to write to
 * @param origin Origin of the 1D plot
 * @param width Width of the node in the plot
 * @param level Level of the node
 */
void BinaryInternalNode::plot(ostream& stream, double origin, double width,
                              int level) {
    stream << "pro\t" << (origin + 0.5 * width) << "\t" << (-level) << "\n";
    int labels[2];
    labels[0] = 100 * _pleft->get_position().x();
    labels[1] = 100 * _pright->get_position().x();
    float roundlabels[2];
    roundlabels[0] = labels[0] / 100.;
    roundlabels[1] = labels[1] / 100.;
    stream << "t\t" << (origin + 0.5 * width - 0.06) << "\t" << (-level) << "\t"
           << roundlabels[0] << "\n";
    stream << "t\t" << (origin + 0.5 * width + 0.01) << "\t" << (-level) << "\t"
           << roundlabels[1] << "\n";
    //    stream << "t\t" << (origin+0.5*width+0.01) << "\t" << (-level) << "\t"
    //    << depth() << "\n";
    stream << "sr-\t" << (origin + 0.5 * width) << "\t"
           << (origin + 0.25 * width) << "\t" << (-level) << "\t"
           << (-level - 1) << "\n";
    if(_children[0]->get_leaf()) {
        stream << "pbo\t" << (origin + 0.25 * width) << "\t" << (-level - 1)
               << "\n";
    } else {
        _children[0]->plot(stream, origin, 0.5 * width, level + 1);
    }
    stream << "sr-\t" << (origin + 0.5 * width) << "\t"
           << (origin + 0.75 * width) << "\t" << (-level) << "\t"
           << (-level - 1) << "\n";
    if(_children[1]->get_leaf()) {
        stream << "pbo\t" << (origin + 0.75 * width) << "\t" << (-level - 1)
               << "\n";
    } else {
        _children[1]->plot(stream, origin + 0.5 * width, 0.5 * width,
                           level + 1);
    }
}

/**
 * @brief Cut off the edge corresponding to this node to limit it to the
 * simulation box
 */
void BinaryInternalNode::trim_edges() {
    if(!_leaf) {
        Vec origin = _edge->get_origin();
        origin.set(origin.x(), origin.y() + 1.);
        Vertex* originpoint = new Vertex(origin);
        if(_edge->get_twin()) {
            _edge->get_twin()->set_origin(originpoint);
        } else {
            _edge->set_twin(new Edge());
            _edge->get_twin()->set_origin(originpoint);
        }

        _children[0]->trim_edges();
        _children[1]->trim_edges();
    }
}

/**
 * @brief Add the edges corresponding to this node to the Voronoi grid
 *
 * @param vorlist Reference to the Voronoi grid
 * @param container DelCont specifying the dimensions of the box
 */
void BinaryInternalNode::get_halfedges(DoublyConnectedEdgeList& vorlist,
                                       DelCont* container) {
    if(!_leaf) {
        if(_edge->get_twin()) {
            Vec origin = _edge->get_twin()->get_origin();
            if(container->inside(origin)) {
                Vec rico = _pleft->get_position() - _pright->get_position();
                if(rico[1]) {
                    double perprico = -rico[0] / rico[1];
                    double dx = 1.;
                    if(origin.x() < 0.5) {
                        dx = -1.;
                    }
                    origin.set(origin.x() + dx, origin.y() + perprico * dx);
                } else {
                    double dy = 2.;
                    if(origin.y() < 0.5) {
                        dy = -2.;
                    }
                    origin.set(origin.x(), origin.y() + dy);
                }
                Vertex* originpoint = new Vertex(origin);
                _edge->set_origin(originpoint);
                vorlist.add_edge(_edge);
            }
        } else {
            // this code is never executed, since every edge gets a twin when it
            // is created
            Vec origin = _edge->get_origin();
            if(container->inside(origin)) {
                Vec rico = _pleft->get_position() - _pright->get_position();
                if(rico[1]) {
                    double perprico = -rico[0] / rico[1];
                    double dx = 1.;
                    if(origin.x() < 0.5) {
                        dx = -1.;
                    }
                    origin.set(origin.x() + dx, origin.y() + perprico * dx);
                } else {
                    double dy = 2.;
                    if(origin.y() < 0.5) {
                        dy = -2.;
                    }
                    origin.set(origin.x(), origin.y() + dy);
                }
                Vertex* originpoint = new Vertex(origin);
                _edge->set_twin(new Edge());
                _edge->get_twin()->set_twin(_edge);
                _edge->get_twin()->set_origin(originpoint);
                vorlist.add_edge(_edge);
                vorlist.add_edge(_edge->get_twin());
            }
        }
        _children[0]->get_halfedges(vorlist, container);
        _children[1]->get_halfedges(vorlist, container);
    }
}

// BinaryLeaf

/**
 * @brief Constructor
 *
 * @param site FortuneSite corresponding to this leaf
 * @param parent Parent node
 */
BinaryLeaf::BinaryLeaf(FortuneSite* site, BinaryInternalNode* parent) {
    _site = site;
    _event = NULL;
    _prev_leaf = NULL;
    _next_leaf = NULL;
    _parent = parent;
}

/**
 * @brief Set a pointer to the next leaf in the tree
 *
 * @param next_leaf Next BinaryLeaf in the tree
 */
void BinaryLeaf::set_next_leaf(BinaryLeaf* next_leaf) {
    _next_leaf = next_leaf;
    next_leaf->set_prev_leaf(this);
}

/**
 * @brief Set a pointer to the previous leaf in the tree
 *
 * @param prev_leaf Previous BinaryLeaf in the tree
 */
void BinaryLeaf::set_prev_leaf(BinaryLeaf* prev_leaf) {
    _prev_leaf = prev_leaf;
}

/**
 * @brief Signal the CircleEvent of this BinaryLeaf to be a false alarm
 *
 * This means another event that affects this leaf will occur before the
 * CircleEvent of this leaf is encountered.
 */
void BinaryLeaf::deactivate() {
    if(_event) {
        _event->false_alarm();
    }
}

/**
 * @brief Get the y-coordinate of the arc corresponding to the FortuneSite of
 * this leaf for the given beachline and x-coordinate
 *
 * @param x x-coordinate at which to evaluate the arc
 * @param ybeach y-coordinate of the current beachline
 * @return y-coordinate of the arc
 */
double BinaryLeaf::get_value(double x, double ybeach) {
    double a = 2. * (_site->get_position().y() - ybeach);
    if(!a) {
        return ybeach;
    }
    a = 1. / a;
    double b = a * _site->get_position().x();
    double c = 0.5 * (_site->get_position().y() + ybeach) +
               b * _site->get_position().x();
    b *= -2.;
    return a * x * x + b * x + c;
}

// BinarySearchTree

/**
 * @brief Constructor
 *
 * Initialize a tree with a single root BinaryInternalNode which is a leaf of
 * the tree.
 */
BinarySearchTree::BinarySearchTree() {
    _root = new BinaryInternalNode(NULL, NULL, NULL);
    BinaryLeaf* leaf = new BinaryLeaf(NULL, NULL);
    _root->make_leaf(leaf);
    //    _first_leaf = NULL;
}

/**
 * @brief Destructor. Recursively free all nodes by deleting the root node
 */
BinarySearchTree::~BinarySearchTree() {
    delete _root;
}

/**
 * @brief Add the given FortuneSite to the tree, using the given beachline
 *
 * @param site FortuneSite to add to the tree
 * @param ybeach y-coordinate of the current beachline
 * @param events Array to add newly created CircleEvent instances to
 * @return BinaryLeaf that was added to the tree
 */
BinaryLeaf* BinarySearchTree::add(FortuneSite* site, double ybeach,
                                  CircleEvent** events) {
    events[0] = NULL;
    events[1] = NULL;
    if(_root->get_leaf()) {
        BinaryLeaf* leaf = _root->get_leaf();
        if(leaf->get_site()) {
            // make node
            FortuneSite* site2 = leaf->get_site();
            _root->erase();
            delete _root;
            delete leaf;
            if(site2->get_position().y() == site->get_position().y()) {
                BinaryLeaf* leaf_right;
                BinaryLeaf* leaf_left;
                BinaryLeaf* new_leaf;
                if(site2->get_position().x() < site->get_position().x()) {
                    _root = new BinaryInternalNode(site2, site, NULL);
                    leaf_left = new BinaryLeaf(site2, _root);
                    leaf_right = new BinaryLeaf(site, _root);
                    new_leaf = leaf_right;
                } else {
                    _root = new BinaryInternalNode(site, site2, NULL);
                    leaf_left = new BinaryLeaf(site, _root);
                    leaf_right = new BinaryLeaf(site2, _root);
                    new_leaf = leaf_left;
                }
                _root->add_leaf(leaf_left, 0);
                _root->add_leaf(leaf_right, 1);
                //                _first_leaf = leaf_left;
                _root->set_edge(new Edge());
                _root->get_edge()->set_twin(new Edge());
                _root->get_edge()->get_twin()->set_twin(_root->get_edge());
                Vec origin(0.5 * (site2->get_position().x() +
                                  site->get_position().x()),
                           100.);
                Vertex* originpoint = new Vertex(origin);
                _root->get_edge()->get_twin()->set_origin(originpoint);
                leaf_left->set_next_leaf(leaf_right);
                return new_leaf;
            }
            BinaryLeaf* leaves_left[2];
            BinaryLeaf* leaf_right;
            _root = new BinaryInternalNode(site, site2, NULL);
            _root->set_edge(new Edge());
            leaf_right = new BinaryLeaf(site2, _root);
            _root->add_subnode(site2, site, 0, leaves_left);
            _root->add_leaf(leaf_right, 1);
            //            _first_leaf = leaves_left[0];

            _root->get_edge()->set_twin(_root->get_child(0)->get_edge());
            _root->get_child(0)->get_edge()->set_twin(_root->get_edge());

            leaves_left[0]->set_next_leaf(leaves_left[1]);
            leaves_left[1]->set_next_leaf(leaf_right);
            return leaves_left[1];
        } else {
            delete leaf;
            leaf = new BinaryLeaf(site, NULL);
            _root->make_leaf(leaf);
            return leaf;
            //            _first_leaf = _root->get_leaf();
        }
    } else {
        return _root->add(site, ybeach, events);
        //        if(_first_leaf->get_prev_leaf()){
        //            _first_leaf = _first_leaf->get_prev_leaf();
        //        }
    }
}

/**
 * @brief Remove the given BinaryLeaf from the tree
 *
 * @param leaf BinaryLeaf to remove
 * @param ybeach y-coordinate of the current beachline
 * @param events Array to store newly created CircleEvent instances in
 * @param edges Array to store newly created Edge instances in
 */
void BinarySearchTree::remove(BinaryLeaf* leaf, double ybeach,
                              CircleEvent** events, Edge** edges) {
    events[0] = NULL;
    events[1] = NULL;
    if(leaf->get_prev_leaf()) {
        leaf->get_prev_leaf()->deactivate();
    }
    if(leaf->get_next_leaf()) {
        leaf->get_next_leaf()->deactivate();
    }
    //    if(leaf == _first_leaf){
    //        _first_leaf = leaf->get_next_leaf();
    //    }
    leaf->get_parent()->remove(leaf, ybeach, events, edges);
    delete leaf;
}

/**
 * @brief Write the tree to a given stream in a format that can be plotted using
 * a specialized Python-script
 *
 * @param stream std::ostream to write to
 */
void BinarySearchTree::plot(ostream& stream) {
    if(_root->get_leaf()) {
        if(_root->get_leaf()->get_site()) {
            stream << "pro\t0.5\t-1\n";
        }
    } else {
        _root->plot(stream, 0., 1., 1);
    }
}

/**
 * @brief Write the beachline to the given stream in a format that can be
 * plotted using a specialized Python-script
 *
 * @deprecated This method no longer works
 *
 * @param stream std::ostream to write to
 * @param ybeach y-coordinate of the current beachline
 */
void BinarySearchTree::plot_beachline(ostream& stream, double ybeach) {
    //    if(_first_leaf){
    //        vector<double> breakpoints;
    //        vector<BinaryLeaf*> arcs;
    //        BinaryLeaf* arc1 = _first_leaf;
    //        BinaryLeaf* arc2 = arc1->get_next_leaf();
    //        while(arc2){
    //            double xbreak = BinaryInternalNode::intersection(
    //    arc1->get_site()->get_position(), arc2->get_site()->get_position(),
    //    ybeach);
    ////            if(arc1->get_site()->get_position().y() ==
    ////    arc2->get_site()->get_position().y()){
    ////                if(arc1->get_site()->get_position().x() >
    ////    arc2->get_site()->get_position().x()){
    ////                    xbreak = 10000;
    ////                }
    ////            }
    //            breakpoints.push_back(xbreak);
    //            arcs.push_back(arc1);
    //            arc1 = arc2;
    //            arc2 = arc1->get_next_leaf();
    //        }
    //        arcs.push_back(arc1);
    //        for(int i = 0; i < 1000; i++){
    //            unsigned int j = 0;
    //            double x = i*0.001;
    //            while(j < breakpoints.size() && x > breakpoints[j]){
    //                j++;
    //            }
    //            stream << x << "\t" << arcs[j]->get_value(x, ybeach) << "\n";
    //        }
    //    }
}

/**
 * @brief Cut off trailing edges at the boundaries of the simulation box
 */
void BinarySearchTree::trim_edges() {
    _root->trim_edges();
}

/**
 * @brief Add the edges of the internal nodes to the Voronoi grid
 *
 * @param vorlist Reference to the Voronoi tesselation
 * @param container DelCont specifying the dimensions of the simulation box
 */
void BinarySearchTree::get_halfedges(DoublyConnectedEdgeList& vorlist,
                                     DelCont* container) {
    _root->get_halfedges(vorlist, container);
}

/**
 * @brief Test the BinarySearchTree::remove() method
 */
void BinarySearchTree::removetest() {
    CircleEvent* events[2];
    Edge* edges[5];
    // ybeach does not matter, leaf is completely random
    remove(_root->get_child(0)
                   ->get_child(1)
                   ->get_child(1)
                   ->get_child(1)
                   ->get_leaf(),
           0., events, edges);
}

/**
 * @brief Print out the tree to the stdout
 *
 * @deprecated Method does no longer work
 */
void BinarySearchTree::print() {
    //    unsigned int numarcs = 0;
    //    BinaryLeaf* arc1 = _first_leaf;
    //    BinaryLeaf* arc2 = arc1->get_next_leaf();
    //    while(arc2){
    //        numarcs++;
    //        cout << arc1->get_site()->get_position().x() << "\t"
    //    << arc1->get_site()->get_position().y() << endl;
    //        arc1 = arc2;
    //        arc2 = arc1->get_next_leaf();
    //    }
    //    cout << arc1->get_site()->get_position().x() << "\t"
    //    << arc1->get_site()->get_position().y() << endl;
    //    cout << numarcs << " arcs" << endl;
}

#endif
