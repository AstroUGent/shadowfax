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
 * @file NgbSearch.cpp
 *
 * @brief Neighbour search treewalks
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "utilities/Particle.hpp"
#include "utilities/Tree.hpp"
#include "utilities/TreeWalker.hpp"
using namespace std;

/**
 * @brief Constructor
 *
 * @param center Center for the neighbour search
 * @param radius Radius of the neighbour search
 * @param exportlist Flags for MPI-communication
 * @param reserve Expected number of neighbours
 */
NgbSearch::NgbSearch(Vec center, double radius, vector<bool>& exportlist,
                     unsigned int reserve)
        : _center(center), _exportlist(exportlist) {
    if(reserve) {
        _ngbs.reserve(reserve);
    }
    _radius2 = radius * radius;
}

/**
 * @brief Get the resulting list of neighbours
 *
 * @return List of neighbours
 */
vector<GasParticle*> NgbSearch::get_ngbs() {
    return _ngbs;
}

/**
 * @brief Decide whether to open the given node
 *
 * The node should be opened if the search sphere overlaps with it.
 *
 * @param node TreeNode on which to operate
 * @return True if the node should be opened, false otherwise
 */
bool NgbSearch::splitnode(TreeNode* node) {
    return node->get_distance(_center) <= _radius2;
}

/**
 * @brief Action to perform on a leaf of the tree
 *
 * If the corresponding Particle falls inside the search sphere, we add it to
 * the list.
 *
 * @param leaf Leaf on which to operate
 */
void NgbSearch::leafaction(Leaf* leaf) {
    Particle* p = leaf->get_particle();
    if(p->type() == PARTTYPE_GAS) {
        if((p->get_position() - _center).norm2() <= _radius2) {
            _ngbs.push_back((GasParticle*)p);
        }
    }
}

/**
 * @brief Action to perform when a PseudoNode is encountered
 *
 * @param pseudonode PseudoNode on which to operate
 */
void NgbSearch::pseudonodeaction(PseudoNode* pseudonode) {
    if(pseudonode->get_distance(_center) <= _radius2) {
        _exportlist[pseudonode->get_source()] = true;
    }
}

// ClosestNgbSearch

/**
 * @brief Constructor
 *
 * @param center Center of the closest neighbour search
 * @param radius Radius of the closest neighbour search
 */
ClosestNgbSearch::ClosestNgbSearch(Vec center, double radius)
        : _center(center) {
    _closest = NULL;
    _radius = radius * radius;
}

/**
 * @brief Set the position for the closest neighbour search
 *
 * @param position Position for the closest neighbour search
 */
void ClosestNgbSearch::set_position(Vec position) {
    _center = position;
}

/**
 * @brief Get the position for the closest neighbour search
 *
 * @return Position for the closest neighbour search
 */
Vec ClosestNgbSearch::get_position() {
    return _center;
}

/**
 * @brief Get the closest Particle
 *
 * @return Result of the treewalk
 */
GasParticle* ClosestNgbSearch::get_closest() {
    return _closest;
}

/**
 * @brief Increase the search radius with a factor 2
 */
void ClosestNgbSearch::increase_radius() {
    _radius *= 4.;
}

/**
 * @brief Decide whether to open the node
 *
 * The node is opened if the current search sphere overlaps with it.
 *
 * @param node TreeNode on which we operate
 * @return True if the node should be opened, false otherwise
 */
bool ClosestNgbSearch::splitnode(TreeNode* node) {
    return node->get_distance(_center) < _radius;
}

/**
 * @brief Action to perform on a single leaf of the tree
 *
 * We check if the corresponding Particle is inside the current search sphere.
 * If it is, we replace the current value of the closest Particle with the
 * Particle and replace the current value of the search radius squared by the
 * distance between the center of the search and the Particle.
 *
 * @param leaf Leaf on which to operate
 */
void ClosestNgbSearch::leafaction(Leaf* leaf) {
    Particle* p = leaf->get_particle();
    if(p->type() == PARTTYPE_GAS) {
        double r2 = (p->get_position() - _center).norm2();
        if(r2 < _radius) {
            _closest = (GasParticle*)p;
            _radius = r2;
        }
    }
}

/**
 * @brief Action to perform when a PseudoNode is encountered
 *
 * @warning Not implemented!
 *
 * @param pseudonode PseudoNode on which to operate
 */
void ClosestNgbSearch::pseudonodeaction(PseudoNode* pseudonode) {
    // do something
}
