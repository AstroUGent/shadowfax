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
 * @file TreeRoute.hpp
 *
 * @brief Find the way to a node of the tree using a route
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TREEROUTE_HPP
#define TREEROUTE_HPP

#include <ostream>

/**
  * @brief A class representing a route through the tree
  *
  * Every step represents a level in the tree and the corresponding path value
  * gives the index of the next node in the current node. In this way, the
  * TreeRoute effectively tells you how to reach a certain node.
  */
class TreeRoute{
private:
#if ndim_==3
    /*! \brief Indices of the child nodes per level of the tree */
    unsigned int _path[20];
#else
    /*! \brief Indices of the child nodes per level of the tree */
    unsigned int _path[30];
#endif
    /*! \brief Current level in the tree */
    unsigned int _curpos;
    /*! \brief Total length of the route in levels of the tree */
    unsigned int _length;

public:
    /**
     * @brief Empty constructor
     */
    TreeRoute() : _curpos(0), _length(0) {}

    /**
     * @brief Add a step to the route
     *
     * We assume that we always successively call this method from the root of
     * the tree. This method adds the given child on the deepest level and
     * increases the total depth of the route in the tree.
     *
     * @param x Index of the node that is the next step in the route
     */
    inline void add_step(unsigned int x){
        _path[_length] = x;
        _length++;
    }

    /**
     * @brief Get the index of the next child in the route
     *
     * We assume that we successively call this function from the root of the
     * tree. This function then always gives the index of the next child node
     * in the current node.
     *
     * @return Index of the next child node in the route
     */
    inline unsigned int next(){
        return _path[_curpos++];
    }

    /**
     * @brief Check if we are at the deepest level of the route
     *
     * @return True if we are at the end of the route, false otherwise
     */
    inline bool finished(){
        return _curpos == _length;
    }

    /**
     * @brief Print the route to the given stream
     *
     * @param stream std::ostream to write to
     */
    inline void print(std::ostream& stream){
        for(unsigned int i = 0; i < _length; i++){
            stream << _path[i] << " ";
        }
        stream << std::endl;
    }
};

#endif // TREEROUTE_HPP
