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
 * @file TimeStepWalker.hpp
 *
 * @brief TreeWalker for hydrodynamical timestep calculation: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TIMESTEPWALKER_HPP
#define TIMESTEPWALKER_HPP

#include "MPIMethods.hpp"             // for MyMPI_Pack, MyMPI_Unpack
#include "utilities/GasParticle.hpp"  // for GasParticle
#include "utilities/TreeWalker.hpp"   // for TreeWalker
#include <algorithm>                  // for min

class DMParticle;
class Leaf;
class PseudoNode;
class StarParticle;
class TreeNode;

/**
 * @brief TreeWalker implemtentation used to calculate the hydrodynamical
 * timestep as discussed in Springel 2010
 */
class TimeStepWalker : public TreeWalker {
  private:
    /*! \brief GasParticle for which the timestep is being calculated */
    GasParticle* _p;

    /*! \brief Position of the particle */
    Vec _position;

    /*! \brief Fluid velocity of the particle */
    Vec _v;

    /*! \brief Magnitude of the fluid velocity of the particle */
    double _vi;

    /*! \brief Current value of the timestep */
    double _t;

    /*! \brief Soundspeed of the particle */
    double _ci;

    /*! \brief Flag indicating if the calculation is done for a particle that
     *  resides on the local MPI process */
    bool _local;

  public:
    /**
     * @brief Auxiliary class used to commmunicate particle information between
     * MPI processes
     */
    class Export {
      private:
        /*! \brief GasParticle for which the timestep is being calculated */
        GasParticle* _p;
        /*! \brief Fluid velocity of the particle (exported) */
        Vec _v;
        /*! \brief Position of the particle (exported) */
        Vec _pos;
        /*! \brief Magnitude of the fluid velocity of the particle (exported) */
        double _vi;
        /*! \brief Soundspeed of the particle (exported) */
        double _ci;
        /*! \brief Current value of the timestep (exported) */
        double _t;

      public:
        /**
         * @brief Constructor
         *
         * @param p GasParticle for which the timestep is being calculated
         * @param v Fluid velocity of the particle
         * @param pos Position of the particle
         * @param vi Magnitude of the velocity of the particle
         * @param ci Soundspeed of the particle
         * @param t Current value of the timestep
         */
        Export(GasParticle* p, Vec v, Vec pos, double vi, double ci, double t)
                : _p(p), _v(v), _pos(pos), _vi(vi), _ci(ci), _t(t) {}

        /**
         * @brief Add data to the given MPI buffer for export
         *
         * @param buffer MPI buffer
         * @param bufsize Buffer size
         * @param position Current position in the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_v[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_pos[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_vi, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_ci, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_t, 1, MPI_DOUBLE, buffer, bufsize, position);
        }

        /**
         * @brief Import data from another MPI process and update the timestep
         * of the particle
         *
         * @param buffer MPI buffer
         * @param bufsize Buffer size
         * @param position Current position in the buffer (is updated)
         */
        void unpack_data(void* buffer, int bufsize, int* position) {
            double t;
            MyMPI_Unpack(buffer, bufsize, position, &t, 1, MPI_DOUBLE);
            // since the particle can be exported to multiple processes,
            // we have no guarantee that this t will be the smallest
            _p->set_real_timestep(std::min(_p->get_real_timestep(), t));
        }
    };

    /**
     * @brief Auxiliary class to communicate particle information between MPI
     * processes
     */
    class Import {
      private:
        /*! \brief Fluid velocity of the particle */
        Vec _v;
        /*! \brief Position of the particle */
        Vec _pos;
        /*! \brief Magnitude of the fluid velocity of the particle */
        double _vi;
        /*! \brief Soundspeed of the particle */
        double _ci;
        /*! \brief Current value of the timestep (exported) */
        double _t;

      public:
        /**
         * @brief Constructor
         *
         * @param buffer MPI buffer
         * @param bufsize Buffer size
         * @param position Current position in the buffer (is updated)
         */
        Import(void* buffer, int bufsize, int* position) {
            MyMPI_Unpack(buffer, bufsize, position, &_v[0], ndim_, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_pos[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_vi, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_ci, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_t, 1, MPI_DOUBLE);
        }

        /**
         * @brief Get the fluid velocity of the particle
         *
         * @return Fluid velocity of the particle
         */
        Vec get_v() {
            return _v;
        }

        /**
         * @brief Get the position of the particle
         *
         * @return Position of the particle
         */
        Vec get_position() {
            return _pos;
        }

        /**
         * @brief Get the magnitude of the fluid velocity of the particle
         *
         * @return Magnitude of the fluid velocity of the particle
         */
        double get_vi() {
            return _vi;
        }

        /**
         * @brief Get the soundspeed of the particle
         *
         * @return Soundspeed of the particle
         */
        double get_ci() {
            return _ci;
        }

        /**
         * @brief Get the current value of the timestep
         *
         * @return Current value of the timestep
         */
        double get_t() {
            return _t;
        }

        /**
         * @brief Set the current value of the timestep
         *
         * @param t Current value of the timestep
         */
        void set_t(double t) {
            _t = t;
        }

        /**
         * @brief Add data to the given MPI buffer for export to the original
         * MPI process
         *
         * @param buffer MPI buffer
         * @param bufsize Buffer size
         * @param position Current position in the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_t, 1, MPI_DOUBLE, buffer, bufsize, position);
        }
    };

    TimeStepWalker(GasParticle* p);

    // needed for compatibility
    /**
     * @brief Dummy constructor
     *
     * Only needed because Tree::walk_tree() needs a DMParticle constructor.
     * This constructor will never be used in practice.
     *
     * @param p DMParticle that compiler-wise could be passed on to the
     * constructor
     */
    TimeStepWalker(DMParticle* p) {}

    /**
     * @brief Dummy constructor
     *
     * Only needed because Tree::walk_tree() needs a StarParticle constructor.
     * This constructor will never be used in practice.
     *
     * @param p StarParticle that compiler-wise could be passed on to the
     * constructor
     */
    TimeStepWalker(StarParticle* p) {}

    TimeStepWalker(Import& import);

    void set_position(Vec position);
    Vec get_position();

    double get_timestep();

    bool splitnode(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);

    bool export_to_pseudonode(PseudoNode* pseudonode);

    void after_walk();

    void after_walk(Import& import);

    Export get_export();
};

#endif  // TIMESTEPWALKER_HPP
