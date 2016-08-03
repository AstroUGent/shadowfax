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
 * @file FortuneVorTess.hpp
 *
 * @brief Voronoi tesselation created using Fortune's algorithm: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef FORTUNEVORTESS_HPP
#define FORTUNEVORTESS_HPP

#if ndim_ == 2

#include "DoublyConnectedEdgeList.hpp"  // for DoublyConnectedEdgeList
#include <iosfwd>                       // for ostream
#include <vector>                       // for vector

class BinaryLeaf;
class BinarySearchTree;
class CircleEvent;
class DelCont;
class FortuneSite;
class GasParticle;
class PriorityQueue;
class SiteEvent;

/**
 * @brief Voronoi tesselation constructed using Fortune's algorithm
 *
 * @warning Does not work (yet).
 */
class FortuneVorTess {
  private:
    /*! @brief Generators of the tesselation */
    std::vector<FortuneSite*> _sites;

    /*! @brief DelCont specifying the dimensions of the simulation box */
    DelCont* _container;

    /*! @brief Representation of the Voronoi tesselation */
    DoublyConnectedEdgeList _vorlist;

    void handleSiteEvent(SiteEvent* event, BinarySearchTree& stree,
                         PriorityQueue& queue);
    void handleCircleEvent(CircleEvent* event, BinarySearchTree& stree,
                           PriorityQueue& queue);

    CircleEvent* check_triple(BinaryLeaf* left, BinaryLeaf* middle,
                              BinaryLeaf* right, double ybeach);

  public:
    FortuneVorTess(DelCont* container, unsigned int npart = 0);
    ~FortuneVorTess();

    void add_point(GasParticle* particle);
    void construct();

    void print_tesselation_gnuplot(std::ostream& stream);
    void print_tesselation(std::ostream& stream);

    void test_connections(std::ostream& stream);
};

#endif

#endif  // FORTUNEVORTESS_HPP
