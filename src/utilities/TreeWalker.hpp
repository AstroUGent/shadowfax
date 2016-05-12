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
 * @file TreeWalker.hpp
 *
 * @brief General interface for treewalks
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TREEWALKER_HPP
#define TREEWALKER_HPP

#include "EwaldTable.hpp"
#include "Particle.hpp"
#include "Vec.hpp"
#include <vector>

class TreeNode;
class Leaf;
class GasParticle;
class Particle;
class PseudoNode;

/**
 * @brief General interface for treewalkers
 *
 * A TreeWalker is a sort of visitor, which walks the tree recursively based on
 * an opening criterion for its nodes.
 * The opening criterion is provided by the function TreeWalker::splitnode().
 * The action that has to be performed on leaves of the tree is provided in the
 * method TreeWalker::leafaction().
 *
 * The TreeWalker interface itself does not provide any restrictions on the
 * action the implementation performs. It only implements auxiliary functions
 * to work with periodic boxes.
 */
class TreeWalker {
  protected:
    /*! \brief Size of the simulation box */
    double _boxsize;

    /*! \brief Half the size of the simulation box */
    double _boxhalf;

    /**
     * @brief Get the smallest version of the given distance within a periodic
     * box
     *
     * If the distance is larger than half the box size, we subtract the box
     * size. If it is smaller than minus half the box size, we add the box size.
     *
     * @param x Distance to trim
     */
    inline void nearest(double& x) {
        if(x > _boxhalf) {
            x -= _boxsize;
        } else {
            if(x < -_boxhalf) {
                x += _boxsize;
            }
        }
    }

    /**
     * @brief Get the smallest version of the given Vec within a periodic box
     *
     * We call TreeWalker::nearest() on its coordinates.
     *
     * @param v Given Vec
     */
    inline void nearest(Vec& v) {
        nearest(v[0]);
        nearest(v[1]);
#if ndim_ == 3
        nearest(v[2]);
#endif
    }

  public:
    /**
     * @brief Empty constructor
     */
    TreeWalker() {
        _boxsize = 0.;
        _boxhalf = 0.;
    }

    /**
     * @brief Set the size of the simulation box
     *
     * @param boxsize Size of the simulation box
     */
    void set_boxsize(double boxsize) {
        _boxsize = boxsize;
        _boxhalf = 0.5 * boxsize;
    }

    /**
     * @brief Dummy splitnode function for implementations that do not implement
     * this function
     *
     * @param node TreeNode on which to operate
     * @return True: we always open all nodes if no splitnode function is
     * implemented
     */
    virtual bool splitnode(TreeNode* node) {
        return true;
    }

    /**
     * @brief Dummy splitaction
     *
     * Does nothing.
     *
     * @param node TreeNode on which to operate
     */
    virtual void splitaction(TreeNode* node) {}

    /**
     * @brief Dummy nodeaction
     *
     * Does nothing.
     *
     * @param node TreeNode on which to operate
     */
    virtual void nodeaction(TreeNode* node) {}

    /**
     * @brief Dummy leafaction
     *
     * Does nothing.
     *
     * @param leaf Leaf on which to operate
     */
    virtual void leafaction(Leaf* leaf) {}

    /**
     * @brief Dummy pseudonodeaction
     *
     * Does nothing.
     *
     * @deprecated Replaced by export_to_pseudonode()
     *
     * @param pseudonode PseudoNode on which to operate
     */
    virtual void pseudonodeaction(PseudoNode* pseudonode) {}

    /**
     * @brief Dummy treewalk finalize method
     *
     * Does nothing.
     */
    virtual void after_walk() {}
};

/**
 * @brief TreeWalker used for periodic tree walks
 *
 * @deprecated No longer used.
 */
class PeriodicTreeWalker : public TreeWalker {
  public:
    /**
     * @brief Set the position used for the treewalk
     * @param position Position for the treewalk
     */
    virtual void set_position(Vec position) = 0;

    /**
     * @brief Get the position used for the treewalk
     *
     * @return Position used for the treewalk
     */
    virtual Vec get_position() = 0;

    /**
     * @brief Periodic version of TreeWalker::splitnode()
     *
     * @deprecated No longer used
     *
     * @param node TreeNode on which to operate
     * @param ewald_table EwaldTable used for periodic correction terms
     * @return True if the node should be opened, false otherwise
     */
    virtual bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table) = 0;

    /**
     * @brief Periodic version of TreeWalker::pseudonodeaction()
     *
     * @param pseudonode PseudoNode on which to operate
     * @param ewald_table EwaldTable used for periodic correction terms
     */
    virtual void periodicpseudonodeaction(PseudoNode* pseudonode,
                                          EwaldTable& ewald_table) = 0;

    /**
     * @brief Periodic version of TreeWalker::leafaction()
     *
     * @param leaf Leaf on which to operate
     * @param ewald_table EwaldTable used for periodic correction terms
     */
    virtual void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table) = 0;
};

/**
 * @brief TreeWalker specialization used to search neighbours for a given
 * GasParticle
 */
class NgbSearch : public TreeWalker {
  private:
    /*! \brief Center of the neighbour search */
    Vec _center;

    /*! \brief Radius squared of the neighbour search */
    double _radius2;

    /*! \brief Neighbours found */
    std::vector<GasParticle*> _ngbs;

    /*! \brief Flags for MPI-communication */
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
    /*! \brief Center of the closest neighbour search */
    Vec _center;

    /*! \brief Radius squared of the closest neighbour search */
    double _radius;

    /*! \brief Current closest particle */
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

#endif  // TREEWALKER_HPP
