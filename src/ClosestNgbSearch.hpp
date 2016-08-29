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
 * @file ClosestNgbSearch.hpp
 *
 * @brief TreeWalker implementation for closest neighbour searching
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef CLOSESTNGBSEARCH_HPP
#define CLOSESTNGBSEARCH_HPP

#include "MPIMethods.hpp"
#include "Vec.hpp"  // for Vec
#include "utilities/StarParticle.hpp"
#include "utilities/TreeWalker.hpp"
#include <vector>  // for vector

class GasParticle;
class Leaf;
class PseudoNode;
class TreeNode;

/**
 * @brief TreeWalker implementation used to find the closest neighbour of a
 * given GasParticle
 */
class ClosestNgbSearch : public PeriodicTreeWalker {
  private:
    /*! @brief StarParticle for which the closest neighbour search is
     *  performed */
    StarParticle* _star;

    /*! @brief Center of the closest neighbour search */
    Vec _center;

    /*! @brief Radius squared of the closest neighbour search */
    double _radius;

    /*! @brief Current closest particle */
    GasParticle* _closest;

  public:
    /**
     * @brief Auxiliary class used to communicate information between MPI
     * processes during the closest neighbour search tree walk
     */
    class Export {
      private:
        /*! @brief Center of the closest neighbour search */
        Vec _center;

        /*! @brief Radius squared of the closest neighbour search */
        double _radius2;

      public:
        /**
         * @brief Constructor
         *
         * @param center Center of the closest neighbour search
         * @param radius2 Radius squared of the closest neighbour search
         */
        Export(Vec center, double radius2) {
            _center = center;
            _radius2 = radius2;
        }

        /**
         * @brief Add the relevant data to the given communication buffer
         *
         * @param buffer Buffer to write to
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_center[0], ndim_, MPI_DOUBLE, buffer, bufsize,
                       position);
            MyMPI_Pack(&_radius2, 1, MPI_DOUBLE, buffer, bufsize, position);
        }

        /**
         * @brief Read the response of the external communication from the given
         * communication buffer
         *
         * @param buffer Buffer to read from
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void unpack_data(void* buffer, int bufsize, int* position) {
            double radius2;
            int rank;
            GasParticle* closest;
            MyMPI_Unpack(buffer, bufsize, position, &radius2, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &rank, 1, MPI_INT);
            // a pointer should be an unsigned 64-bit integer on 64-bit systems
            MyMPI_Unpack(buffer, bufsize, position, &closest, 1,
                         MPI_UNSIGNED_LONG);
        }
    };

    /**
     * @brief Auxiliary class used to communicate particle information between
     * MPI processes during the gravity treewalk.
     */
    class Import {
      private:
        /*! @brief Position for which the treewalk is performed */
        Vec _pos;

        /*! @brief Old acceleration, used for the relative opening criterion */
        double _olda;

        /*! @brief Local result of the treewalk (is exported) */
        Vec _a;

        /*! @brief Computational cost of the local treewalk (is exported) */
        float _comp_cost;

        /*! @brief Softening length of the Particle for which the treewalk is
         *  performed */
        double _hsoft;

        /*! @brief \f$\eta\f$-parameter for variable softening lengths (is
            exported) */
        double _eta;

      public:
        /**
         * @brief Constructor. Initialize the import based on the given
         * communication stream
         *
         * @param buffer Buffer to read from
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        Import(void* buffer, int bufsize, int* position) {
            MyMPI_Unpack(buffer, bufsize, position, &_pos[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_olda, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_hsoft, 1, MPI_DOUBLE);
        }

        /**
         * @brief Get the position for which the treewalk is performed
         *
         * @return The position of the external Particle
         */
        Vec get_position() {
            return _pos;
        }

        /**
         * @brief Get the old acceleration of the external Particle
         *
         * @return Old acceleration of the external Particle
         */
        double get_olda() {
            return _olda;
        }

        /**
         * @brief Get the softening length of the external Particle
         *
         * @return Softening length of the external Particle
         */
        double get_hsoft() {
            return _hsoft;
        }

        /**
         * @brief set_a Set the acceleration
         *
         * @param a Value of the acceleration
         */
        void set_a(Vec a) {
            _a = a;
        }

        /**
         * @brief Set the computational cost
         * @param comp_cost Value of the computational cost
         */
        void set_comp_cost(float comp_cost) {
            _comp_cost = comp_cost;
        }

        /**
         * @brief Set the \f$\eta\f$-factor for variable softening lengths
         *
         * @param eta Value of the \f$\eta\f$-factor
         */
        void set_eta(double eta) {
            _eta = eta;
        }

        /**
         * @brief Add the relevant data to the given communication buffer to
         * send them back to the original process
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_a[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_comp_cost, 1, MPI_FLOAT, buffer, bufsize, position);
        }
    };

    ClosestNgbSearch(StarParticle* star, Vec center, double radius);
    void set_position(Vec position);
    Vec get_position();

    GasParticle* get_closest();

    void increase_radius();

    bool splitnode(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);

    bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table);
    void periodicpseudonodeaction(PseudoNode* node, EwaldTable& ewald_table);
    void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table);

    void after_walk();

    void after_walk(Import& import);
};

#endif  // CLOSESTNGBSEARCH_HPP
