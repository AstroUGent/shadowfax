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
 * @file ClosestNgbSearch.cpp
 *
 * @brief Neighbour search treewalks
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ClosestNgbSearch.hpp"
#include "Vec.hpp"                      // for Vec, operator-
#include "utilities/GasParticle.hpp"    // for GasParticle
#include "utilities/Particle.hpp"       // for Particle
#include "utilities/ParticleTypes.hpp"  // for ParticleType::PARTTYPE_GAS
#include "utilities/Tree.hpp"           // for Leaf, PseudoNode, TreeNode
#include <cfloat>
#include <cstddef>  // for NULL
#include <vector>   // for vector, vector<>::reference
using namespace std;

/**
 * @brief Constructor
 *
 * @param star StarParticle for which the closest neighbour search is performed
 * @param center Center of the closest neighbour search
 * @param radius Radius of the closest neighbour search
 */
ClosestNgbSearch::ClosestNgbSearch(StarParticle* star, Vec center,
                                   double radius = DBL_MAX)
        : _center(center) {
    _star = star;
    _closest = NULL;
    _radius = radius * radius;
}

/**
 * @brief Import constructor
 *
 * @param import Import containing data communicated from another MPI process
 */
ClosestNgbSearch::ClosestNgbSearch(Import& import)
        : _center(import.get_center()) {
    _star = NULL;
    _closest = NULL;
    _radius = import.get_radius2();
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
    Vec r = _center - node->get_center();
    if(_boxsize) {
        nearest(r);
    }
    for(unsigned int i = 0; i < ndim_; ++i) {
        if(fabs(r[i]) >= 0.5 * node->get_width()) {
            if(r[i] > 0.) {
                r[i] -= 0.5 * node->get_width();
            } else {
                r[i] += 0.5 * node->get_width();
            }
        } else {
            r[i] = 0.;
        }
    }
    return r.norm2() <= _radius;
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

/**
 * @brief Old periodic boundary treatment, method not used
 *
 * @param node TreeNode
 * @param ewald_table EwaldTable
 * @return False, since this method is never used anyway
 */
bool ClosestNgbSearch::periodicsplitnode(TreeNode* node,
                                         EwaldTable& ewald_table) {
    return false;
}

/**
 * @brief Old periodic boundary treatment, method not used
 *
 * @param node TreeNode
 * @param ewald_table EwaldTable
 */
void ClosestNgbSearch::periodicpseudonodeaction(PseudoNode* node,
                                                EwaldTable& ewald_table) {
    // do something
}

/**
 * @brief Old periodic boundary treatment, method not used
 *
 * @param leaf Leaf
 * @param ewald_table EwaldTable
 */
void ClosestNgbSearch::periodicleafaction(Leaf* leaf, EwaldTable& ewald_table) {
    // do something
}

/**
 * @brief Set the closest neighbour of the StarParticle
 */
void ClosestNgbSearch::after_walk() {
    _star->set_closest_gasparticle(_closest, _radius, MPIGlobal::rank);
}

/**
 * @brief Set the closest neighbour of the Import that is then communicated back
 * to the original StarParticle
 *
 * @param import Import to update
 */
void ClosestNgbSearch::after_walk(Import& import) {
    import.set_closest(_closest);
    import.set_radius2(_radius);
}
