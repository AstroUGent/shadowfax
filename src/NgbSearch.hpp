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
 * @file NgbSearch.hpp
 *
 * @brief TreeWalker implementations for neighbour searching
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef NGBSEARCH_HPP
#define NGBSEARCH_HPP

#include "Vec.hpp"  // for Vec
#include "utilities/TreeWalker.hpp"
#include <vector>  // for vector

class GasParticle;
class Leaf;
class PseudoNode;
class TreeNode;

/**
 * @brief TreeWalker specialization used to search neighbours for a given
 * GasParticle
 */
class NgbSearch : public TreeWalker {
  private:
    /*! @brief Center of the neighbour search */
    Vec _center;

    /*! @brief Radius squared of the neighbour search */
    double _radius2;

    /*! @brief Neighbours found */
    std::vector<GasParticle*> _ngbs;

    /*! @brief Flags for MPI-communication */
    std::vector<bool>& _exportlist;

  public:
    NgbSearch(Vec center, double radius, std::vector<bool>& exportlist,
              unsigned int reserve = 0);

    std::vector<GasParticle*> get_ngbs();

    bool splitnode(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);
};

/**
 * @brief TreeWalker implementation used to find the closest neighbour of a
 * given GasParticle
 */
class ClosestNgbSearch : public TreeWalker {
  private:
    /*! @brief Center of the closest neighbour search */
    Vec _center;

    /*! @brief Radius squared of the closest neighbour search */
    double _radius;

    /*! @brief Current closest particle */
    GasParticle* _closest;

  public:
    ClosestNgbSearch(Vec center, double radius);
    void set_position(Vec position);
    Vec get_position();

    GasParticle* get_closest();

    void increase_radius();

    bool splitnode(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);
};

#endif  // NGBSEARCH_HPP
