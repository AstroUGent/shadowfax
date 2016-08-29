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
#include "NgbSearch.hpp"
#include "Vec.hpp"                      // for Vec, operator-
#include "utilities/GasParticle.hpp"    // for GasParticle
#include "utilities/Particle.hpp"       // for Particle
#include "utilities/ParticleTypes.hpp"  // for ParticleType::PARTTYPE_GAS
#include "utilities/Tree.hpp"           // for Leaf, PseudoNode, TreeNode
#include <cstddef>                      // for NULL
#include <vector>                       // for vector, vector<>::reference
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
