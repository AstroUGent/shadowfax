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
 * @file TimeStepWalker.cpp
 *
 * @brief TreeWalker for hydrodynamical timestep calculation: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "TimeStepWalker.hpp"
#include "StateVector.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/Tree.hpp"
using namespace std;

/**
 * @brief Constructor
 *
 * @param p GasParticle for which the timestep is calculated
 */
TimeStepWalker::TimeStepWalker(GasParticle* p) {
    _p = p;
    _position = p->get_position();
    StateVector W = p->get_Wvec();
#if ndim_ == 3
    _v.set(W[1], W[2], W[3]);
#else
    _v.set(W[1], W[2]);
#endif
    _v -= p->get_velocity();
    _vi = _v.norm();
    _ci = _p->get_soundspeed();
    _t = p->h() / (_vi + p->get_soundspeed());
    _local = true;
}

/**
 * @brief Import constructor
 *
 * @param import Import holding information about the external GasParticle for
 * which the timestep is calculated
 */
TimeStepWalker::TimeStepWalker(Import& import) {
    _p = NULL;
    _position = import.get_position();
    _v = import.get_v();
    _vi = import.get_vi();
    _ci = import.get_ci();
    _t = import.get_t();
    _local = false;
}

/**
 * @brief Set the position for the tree walk
 *
 * @param position Position for the tree walk
 */
void TimeStepWalker::set_position(Vec position) { _position = position; }

/**
 * @brief Get the position for the tree walk
 *
 * @return Position for the tree walk
 */
Vec TimeStepWalker::get_position() { return _position; }

/**
 * @brief Get the resulting timestep
 *
 * @return Resulting particle timestep
 */
double TimeStepWalker::get_timestep() { return _t; }

/**
 * @brief Decide if the given TreeNode should be opened or treated as a whole
 *
 * In the latter case, the required action is also performed by this function.
 *
 * @param node TreeNode on which to act
 * @return True if the node should be opened, false if the tree walk can
 * continue without opening the node
 */
bool TimeStepWalker::splitnode(TreeNode* node) {
    double radius = _t * (_vi + _ci + node->get_cmax() + node->get_vmax());
    radius *= radius;
    double distance = node->get_distance(_position);
    if(_boxsize) {
        nearest(distance);
    }
    return distance <= radius;
}

/**
 * @brief Action to perform on a single Leaf of the tree
 *
 * @param leaf Leaf on which to act
 */
void TimeStepWalker::leafaction(Leaf* leaf) {
    Particle* pj = leaf->get_particle();
    if(pj->type() == PARTTYPE_GAS) {
        Vec rij = _position - pj->get_position();
        if(_boxsize) {
            nearest(rij);
        }
        double r2 = rij.norm2();
        if(r2) {
            double radius =
                    _t * (_vi + _ci + leaf->get_cmax() + leaf->get_vmax());
            radius *= radius;
            if(r2 <= radius) {
                StateVector Wj = ((GasParticle*)pj)->get_Wvec();
#if ndim_ == 3
                Vec vij(Wj[1], Wj[2], Wj[3]);
#else
                Vec vij(Wj[1], Wj[2]);
#endif
                vij = _v - vij;
                double rnrm = sqrt(r2);
                double denom = _ci + ((GasParticle*)pj)->get_soundspeed() -
                               Vec::dot(vij, rij) / rnrm;
                if(denom > 0.) {
                    _t = std::min(_t, rnrm / denom);
                }
            }
        }
    }
}

/**
 * @brief Action to perform when a PseudoNode is encountered
 *
 * @deprecated Replaced by export_to_pseudonode()
 *
 * @param pseudonode PseudoNode on which to act
 */
void TimeStepWalker::pseudonodeaction(PseudoNode* pseudonode) {
    double radius =
            _t * (_vi + _ci + pseudonode->get_cmax() + pseudonode->get_vmax());
    radius *= radius;
    double distance = pseudonode->get_distance(_position);
    if(_boxsize) {
        nearest(distance);
    }
    if(distance <= radius) {
        //        _exportlist[pseudonode->get_source()] = true;
    }
}

/**
 * @brief Decide if the treewalk should be exported to the MPI process
 * corresponding to the given PseudoNode, or if the node can be treated as a
 * whole
 *
 * In the latter case, the required action is performed by this function.
 *
 * @param pseudonode PseudoNode on which to act
 * @return True if the treewalk should be exported, false otherwise
 */
bool TimeStepWalker::export_to_pseudonode(PseudoNode* pseudonode) {
    double radius =
            _t * (_vi + _ci + pseudonode->get_cmax() + pseudonode->get_vmax());
    radius *= radius;
    double distance = pseudonode->get_distance(_position);
    if(_boxsize) {
        nearest(distance);
    }
    return distance <= radius;
}

/**
 * @brief Finalize the treewalk by setting the timestep of the GasParticle
 */
void TimeStepWalker::after_walk() { _p->set_real_timestep(_t); }

/**
 * @brief Finalize the external treewalk by setting the timestep of the Import
 *
 * @param import Import that will be exported back to the original MPI process
 */
void TimeStepWalker::after_walk(Import& import) { import.set_t(_t); }

/**
 * @brief Get an Export to export the treewalk to another MPI process
 *
 * @return Export that can be exported to other MPI processes
 */
TimeStepWalker::Export TimeStepWalker::get_export() {
    return Export(_p, _v, _position, _vi, _ci, _t);
}
