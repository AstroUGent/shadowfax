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
 * @file GravityWalker.cpp
 *
 * @brief Treewalkers used for gravity tree walks: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "GravityWalker.hpp"
#include "utilities/EwaldTable.hpp"
#include "utilities/Particle.hpp"  // for Particle
#include "utilities/Tree.hpp"      // for TreeNode, etc
using namespace std;

// GravityWalker

/**
 * @brief Constructor
 *
 * @param p Particle for which the gravitational acceleration should be
 * calculated
 * @param local Flag indicating whether the Particle resides on the local
 * MPI-process
 */
GravityWalker::GravityWalker(Particle* p, bool local) {
    _p = p;
    _position = p->get_position();
    _olda = p->get_old_acceleration();
    _hsoft = p->get_hsoft();
    _hsoftinv = 1. / _hsoft;
    _hsoftinv3 = _hsoftinv * _hsoftinv * _hsoftinv;
    _local = local;
    _comp_cost = 0;
    if(_p->type() == PARTTYPE_GAS) {
        _kernel = new GasGravityKernel();
        _eta = 0.;
    } else {
        _kernel = new DMGravityKernel();
    }
    _type = _p->type();
}

/**
 * @brief Import constructor.
 *
 * Initialize a GravityWalker for a Particle residing on another MPI-process
 *
 * @param import Import with imported data
 */
GravityWalker::GravityWalker(Import& import) {
    _p = NULL;
    _position = import.get_position();
    _olda = import.get_olda();
    _hsoft = import.get_hsoft();
    _hsoftinv = 1. / _hsoft;
    _hsoftinv3 = _hsoftinv * _hsoftinv * _hsoftinv;
    _local = false;
    _comp_cost = 0;
    if(import.get_type() == PARTTYPE_GAS) {
        _kernel = new GasGravityKernel();
        _eta = 0.;
    } else {
        _kernel = new DMGravityKernel();
    }
    _type = import.get_type();
}

/**
 * @brief Destructor
 *
 * Free GravityKernel memory.
 */
GravityWalker::~GravityWalker() {
    delete _kernel;
}

/**
 * @brief Get the resulting gravitational acceleration
 *
 * @return The result of the treewalk
 */
Vec GravityWalker::get_acceleration() {
    return _a;
}

/**
 * @brief Decide whether the given TreeNode should be opened or can be treated
 * as a whole
 *
 * If it can be treated as a whole, we immediately treat it, as this makes more
 * efficient use of the intermediate quantities.
 *
 * @param node TreeNode on which we operate
 * @return True if the TreeNode should be opened, false otherwise
 */
bool GravityWalker::splitnode(TreeNode* node) {
    Vec r = node->get_center_of_mass_node() - _position;
    if(_boxsize) {
        nearest(r);
    }
    double r2 = r.norm2();
    if(node->get_width() * node->get_width() * node->get_mass_node() <=
       r2 * r2 * _olda) {
        Vec center = node->get_center();
        double nodewidth = 0.6 * node->get_width();
        if(abs(center.x() - _position.x()) < nodewidth) {
            if(abs(center.y() - _position.y()) < nodewidth) {
#if ndim_ == 3
                if(abs(center.z() - _position.z()) < nodewidth) {
                    return true;
                }
#else
                return true;
#endif
            }
        }
        double kernel;
        double rnrm = sqrt(r2);
        if(rnrm >= _hsoft && rnrm >= node->get_hmax()) {
            kernel = node->get_mass_node() / (rnrm * r2);
        } else {
            // we open a node if the force has to be softened
            return true;
        }
        if(_local || node->is_local()) {
            _a += r * kernel;
            _comp_cost++;
        }
        return false;
    } else {
        return true;
    }
}

/**
 * @brief Action to perform on a TreeNode as a whole
 *
 * @deprecated Action is already performed by GravityWalker::splitnode()
 *
 * @param node TreeNode on which we operate
 */
void GravityWalker::nodeaction(TreeNode* node) {
    // obsolete, action is perfomed by splitnode()
}

/**
 * @brief Action to perform on a single Leaf of the Tree
 *
 * @param leaf Leaf on which we operate
 */
void GravityWalker::leafaction(Leaf* leaf) {
    Particle* pj = leaf->get_particle();
    Vec r = pj->get_position() - _position;
    if(_boxsize) {
        nearest(r);
    }
    double m = pj->get_mass();
    double hj = pj->get_hsoft();
    double kernel;
    // we need both r2 and rnrm for all cases
    double r2 = r.norm2();
    double rnrm = sqrt(r2);
    if(rnrm >= _hsoft && rnrm >= hj) {
        kernel = m / (rnrm * r2);
    } else {
        double kernel1, kernel2;
        if(rnrm >= _hsoft) {
            kernel1 = m / (rnrm * r2);
        } else {
            double u = rnrm * _hsoftinv;
            kernel1 = m * (*_kernel)(u, _hsoftinv3);

            if(_type == PARTTYPE_GAS) {
                _eta += m * _kernel->derivative(u, _hsoftinv);
            }
        }
        if(rnrm >= hj) {
            kernel2 = m / (rnrm * r2);
        } else {
            double hjinv = 1. / hj;
            double hjinv2 = hjinv * hjinv;
            double hjinv3 = hjinv2 * hjinv;
            double u = rnrm * hjinv;
            kernel2 = m * (*_kernel)(u, hjinv3);
        }
        kernel = 0.5 * (kernel1 + kernel2);
    }
    _a += r * kernel;
    _comp_cost++;
}

/**
 * @brief Action performed when a PseudoNode is encountered
 *
 * @deprecated Method has been replaced by GravityWalker::export_to_pseudonode
 *
 * @param pseudonode PseudoNode on which we operate
 */
void GravityWalker::pseudonodeaction(PseudoNode* pseudonode) {
    // do nothing
}

/**
 * @brief Decide if the Particle should be exported to the process which holds
 * the actual node represented by the given PseudoNode
 *
 * @param pseudonode PseudoNode on which we operate
 * @return True, since the gravitational force is only approximated for a local
 * TreeNode
 */
bool GravityWalker::export_to_pseudonode(PseudoNode* pseudonode) {
    return true;
}

/**
 * @brief Set the position which is used for the treewalk
 *
 * @deprecated This method should not be used anymore
 *
 * @param position New position for the treewalk
 */
void GravityWalker::set_position(Vec position) {
    _position = position;
}

/**
 * @brief Get the position used for the treewalk
 *
 * @return Position used for the treewalk
 */
Vec GravityWalker::get_position() {
    return _position;
}

/**
 * @brief Periodic version of GravityWalker::splitnode()
 *
 * @deprecated Periodic boundaries are treated in a different way now
 *
 * @param node TreeNode on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 * @return True if the TreeNode should be opened, false otherwise
 */
bool GravityWalker::periodicsplitnode(TreeNode* node, EwaldTable& ewald_table) {
    Vec r = node->get_center_of_mass_node() - _position;
    nearest(r);
    double r2 = r.norm2();
    bool openflag = false;
    if(node->get_width() * node->get_width() * node->get_mass_node() >
       r2 * r2 * _olda) {
        openflag = true;
    } else {
        Vec center = node->get_center();
        double nodewidth = 0.6 * node->get_width();
        if(abs(center.x() - _position.x()) < nodewidth) {
            if(abs(center.y() - _position.y()) < nodewidth) {
#if ndim_ == 3
                if(abs(center.z() - _position.z()) < nodewidth) {
                    openflag = true;
                }
#else
                openflag = true;
#endif
            }
        }
    }
    if(openflag) {
        Vec center = node->get_center();
        for(unsigned int i = 0; i < ndim_; i++) {
            double u = center[i] - _position[i];
            nearest(u);
            if(abs(u) >= 0.5 * (_boxsize - node->get_width())) {
                return true;
            }
        }
        if(node->get_width() > 0.2 * _boxsize) {
            return true;
        }
    }
    _a += node->get_mass_node() * ewald_table.get_correction(r);
    return false;
}

/**
 * @brief Periodic version of GravityWalker::pseudonodeaction()
 *
 * @deprecated Method was never implemented and should not be used, since
 * periodic boundaries are treated differently now
 *
 * @param pseudonode PseudoNode on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void GravityWalker::periodicpseudonodeaction(PseudoNode* pseudonode,
                                             EwaldTable& ewald_table) {}

/**
 * @brief Periodic version of GravityWalker::leafaction()
 *
 * @deprecated Not used anymore
 *
 * @param leaf Leaf on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void GravityWalker::periodicleafaction(Leaf* leaf, EwaldTable& ewald_table) {
    Vec r = leaf->get_particle()->get_position() - _position;
    nearest(r);
    _a += leaf->get_particle()->get_mass() * ewald_table.get_correction(r);
}

/**
 * @brief Finalize the treewalk by setting the gravitational acceleration of the
 * Particle
 */
void GravityWalker::after_walk() {
    _p->set_gravitational_acceleration(_a);
    _p->add_comp_cost(_comp_cost);
    if(_p->type() == PARTTYPE_GAS) {
        GasParticle* gas = (GasParticle*)_p;
        gas->set_eta(_eta);
    }
}

/**
 * @brief Get an Export to export the Particle to another MPI-process
 *
 * @return An Export that can be used to communicate relevant Particle data to
 * another MPI-process
 */
GravityWalker::Export GravityWalker::get_export() {
    return Export(_p, _olda);
}

/**
 * @brief Finalize the local treewalk for another MPI-process by setting the
 * acceleration of the Import
 *
 * @param import Import that will be sent back to the original MPI-process
 */
void GravityWalker::after_walk(Import& import) {
    import.set_a(_a);
    import.set_comp_cost(_comp_cost);
}

// EwaldGravityWalker

/**
 * @brief Constructor
 *
 * @param p Particle for which the Ewald correction force is calculated
 * @param ewald_table Reference to the EwaldTable used to calculate Ewald
 * corrections
 * @param local Flag signaling if the tree walk is done for a particle on the
 * local process
 */
EwaldGravityWalker::EwaldGravityWalker(Particle* p, EwaldTable& ewald_table,
                                       bool local)
        : _ewald_table(ewald_table) {
    _p = p;
    _position = p->get_position();
    _olda = p->get_old_acceleration();
}

/**
 * @brief Import constructor.
 *
 * Initialize an EwaldGravityWalker for a Particle residing on another
 * MPI-process
 *
 * @param import Import with imported data
 * @param ewald_table Reference to the EwaldTable used to calculate Ewald
 * corrections
 */
EwaldGravityWalker::EwaldGravityWalker(Import& import, EwaldTable& ewald_table)
        : _ewald_table(ewald_table) {
    _p = NULL;
    _position = import.get_position();
    _olda = import.get_olda();
}

/**
 * @brief Destructor
 */
EwaldGravityWalker::~EwaldGravityWalker() {}

/**
 * @brief Get the resulting Ewald correction
 *
 * @return The result of the treewalk
 */
Vec EwaldGravityWalker::get_acceleration() {
    return _acorr;
}

/**
 * @brief Decide whether the given TreeNode should be opened or can be treated
 * as a whole
 *
 * If it can be treated as a whole, we immediately treat it, as this makes more
 * efficient use of the intermediate quantities.
 *
 * @param node TreeNode on which we operate
 * @return True if the TreeNode should be opened, false otherwise
 */
bool EwaldGravityWalker::splitnode(TreeNode* node) {
    Vec r = node->get_center_of_mass_node() - _position;
    nearest(r);
    double r2 = r.norm2();
    bool openflag = false;
    if(node->get_width() * node->get_width() * node->get_mass_node() >
       r2 * r2 * _olda) {
        openflag = true;
    } else {
        Vec center = node->get_center();
        double nodewidth = 0.6 * node->get_width();
        if(abs(center.x() - _position.x()) < nodewidth) {
            if(abs(center.y() - _position.y()) < nodewidth) {
#if ndim_ == 3
                if(abs(center.z() - _position.z()) < nodewidth) {
                    openflag = true;
                }
#else
                openflag = true;
#endif
            }
        }
    }
    if(openflag) {
        Vec center = node->get_center();
        for(unsigned int i = 0; i < ndim_; i++) {
            double u = center[i] - _position[i];
            nearest(u);
            if(abs(u) >= 0.5 * (_boxsize - node->get_width())) {
                return true;
            }
        }
        if(node->get_width() > 0.2 * _boxsize) {
            return true;
        }
    }
    _acorr += node->get_mass_node() * _ewald_table.get_correction(r);
    return false;
}

/**
 * @brief Action to perform on a TreeNode as a whole
 *
 * @deprecated Action is already performed by EwaldGravityWalker::splitnode()
 *
 * @param node TreeNode on which we operate
 */
void EwaldGravityWalker::nodeaction(TreeNode* node) {
    // obsolete, action is perfomed by splitnode()
}

/**
 * @brief Action to perform on a single Leaf of the Tree
 *
 * @param leaf Leaf on which we operate
 */
void EwaldGravityWalker::leafaction(Leaf* leaf) {
    Vec r = leaf->get_particle()->get_position() - _position;
    nearest(r);
    _acorr += leaf->get_particle()->get_mass() * _ewald_table.get_correction(r);
}

/**
 * @brief Action performed when a PseudoNode is encountered
 *
 * @deprecated Method has been replaced by
 * EwaldGravityWalker::export_to_pseudonode
 *
 * @param pseudonode PseudoNode on which we operate
 */
void EwaldGravityWalker::pseudonodeaction(PseudoNode* pseudonode) {
    // do something
    //    _exportlist[pseudonode->get_source()] = true;
}

/**
 * @brief Decide if the Particle should be exported to the process which holds
 * the actual node represented by the given PseudoNode
 *
 * @param pseudonode PseudoNode on which we operate
 * @return True, since the gravitational force is only approximated for a local
 * TreeNode
 */
bool EwaldGravityWalker::export_to_pseudonode(PseudoNode* pseudonode) {
    return true;
}

/**
 * @brief Set the position which is used for the treewalk
 *
 * @deprecated This method should not be used anymore
 *
 * @param position New position for the treewalk
 */
void EwaldGravityWalker::set_position(Vec position) {
    _position = position;
}

/**
 * @brief Get the position used for the treewalk
 *
 * @return Position used for the treewalk
 */
Vec EwaldGravityWalker::get_position() {
    return _position;
}

/**
 * @brief Periodic version of EwaldGravityWalker::splitnode()
 *
 * @deprecated Periodic boundaries are treated in a different way now
 *
 * @param node TreeNode on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 * @return True if the TreeNode should be opened, false otherwise
 */
bool EwaldGravityWalker::periodicsplitnode(TreeNode* node,
                                           EwaldTable& ewald_table) {
    return false;
}

/**
 * @brief Periodic version of EwaldGravityWalker::pseudonodeaction()
 *
 * @deprecated Method was never implemented and should not be used, since
 * periodic boundaries are treated differently now
 *
 * @param pseudonode PseudoNode on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void EwaldGravityWalker::periodicpseudonodeaction(PseudoNode* pseudonode,
                                                  EwaldTable& ewald_table) {}

/**
 * @brief Periodic version of EwaldGravityWalker::leafaction()
 *
 * @deprecated Not used anymore
 *
 * @param leaf Leaf on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void EwaldGravityWalker::periodicleafaction(Leaf* leaf,
                                            EwaldTable& ewald_table) {}

/**
 * @brief Finalize the treewalk by setting the gravitational acceleration of the
 * Particle
 */
void EwaldGravityWalker::after_walk() {
    _p->set_gravitational_acceleration(_p->get_gravitational_acceleration() +
                                       _acorr);
}

/**
 * @brief Get an Export to export the Particle to another MPI-process
 *
 * @return An Export that can be used to communicate relevant Particle data to
 * another MPI-process
 */
EwaldGravityWalker::Export EwaldGravityWalker::get_export() {
    return Export(_p, _olda);
}

/**
 * @brief Finalize the local treewalk for another MPI-process by setting the
 * acceleration of the Import
 *
 * @param import Import that will be sent back to the original MPI-process
 */
void EwaldGravityWalker::after_walk(Import& import) {
    import.set_acorr(_acorr);
}

// PotentialWalker

/**
 * @brief Constructor
 *
 * @param p Particle for which to calculate the gravitational potential
 * @param local Flag indicating whether the treewalk is performed for a Particle
 * residing on the local MPI-process
 */
PotentialWalker::PotentialWalker(Particle* p, bool local) {
    _p = p;
    _position = p->get_position();
    _olda = p->get_old_acceleration();
    _hsoft = p->get_hsoft();
    _hsoftinv = 1. / _hsoft;
    _hsoftinv3 = _hsoftinv * _hsoftinv * _hsoftinv;
    _local = local;
    _comp_cost = 0;
    if(p->type() == PARTTYPE_GAS) {
        _kernel = new GasGravityKernel();
    } else {
        _kernel = new DMGravityKernel();
    }

    // we should probably correct for the self-potential of the particle...
    _epot = 0.;
}

/**
 * @brief Import constructor
 *
 * Initialize a PotentialWalker for a Particle residing on another MPI-process.
 *
 * @param import Import containing imported data from another MPI-process
 */
PotentialWalker::PotentialWalker(Import& import) {
    _p = NULL;
    _position = import.get_position();
    _olda = import.get_olda();
    _hsoft = import.get_hsoft();
    _hsoftinv = 1. / _hsoft;
    _hsoftinv3 = _hsoftinv * _hsoftinv * _hsoftinv;
    _local = false;
    _comp_cost = 0;
    if(import.get_type() == PARTTYPE_GAS) {
        _kernel = new GasGravityKernel();
    } else {
        _kernel = new DMGravityKernel();
    }
    _epot = 0.;
}

/**
 * @brief Destructor
 *
 * Free GravityKernel memory.
 */
PotentialWalker::~PotentialWalker() {
    delete _kernel;
}

/**
 * @brief Get the resulting gravitational potential
 *
 * @return The result of the treewalk
 */
double PotentialWalker::get_epot() {
    return _epot;
}

/**
 * @brief Decide if the given TreeNode should be opened
 *
 * If it is not opened, the TreeNode is treated as a whole in this function.
 *
 * @param node TreeNode on which to operate
 * @return True if the TreeNode should be opened, false otherwise
 */
bool PotentialWalker::splitnode(TreeNode* node) {
    Vec r = node->get_center_of_mass_node() - _position;
    if(_boxsize) {
        nearest(r);
    }
    double r2 = r.norm2();
    if(node->get_width() * node->get_width() * node->get_mass_node() <=
       r2 * r2 * _olda) {
        Vec center = node->get_center();
        double nodewidth = 0.6 * node->get_width();
        if(abs(center.x() - _position.x()) < nodewidth) {
            if(abs(center.y() - _position.y()) < nodewidth) {
#if ndim_ == 3
                if(abs(center.z() - _position.z()) < nodewidth) {
                    return true;
                }
#else
                return true;
#endif
            }
        }
        double rnrm = sqrt(r2);
        if(rnrm >= _hsoft && rnrm >= node->get_hmax()) {
            if(_local || node->is_local()) {
                _epot -= node->get_mass_node() / rnrm;
            }
        } else {
            // force opening of node
            return true;
        }
        if(_local || node->is_local()) {
            _comp_cost++;
        }
        return false;
    } else {
        return true;
    }
}

/**
 * @brief Action to perform on a single TreeNode of the Tree
 *
 * This method does nothing, since the actual action is performed by
 * PotentialWalker::splitnode().
 *
 * @param node TreeNode on which we operate
 */
void PotentialWalker::nodeaction(TreeNode* node) {
    // obsolete, action is performed by splitnode()
}

/**
 * @brief Action to perform on a single Leaf of the Tree
 *
 * @param leaf Leaf on which we operate
 */
void PotentialWalker::leafaction(Leaf* leaf) {
    Particle* pj = leaf->get_particle();
    Vec r = pj->get_position() - _position;
    if(_boxsize) {
        nearest(r);
    }
    double m = pj->get_mass();
    double hj = pj->get_hsoft();
    double rnrm = r.norm();
    double kernel;
    if(rnrm >= _hsoft && rnrm >= hj) {
        kernel = -m / rnrm;
    } else {
        double kernel1, kernel2;
        if(rnrm >= _hsoft) {
            kernel1 = -m / rnrm;
        } else {
            double u = rnrm * _hsoftinv;
            kernel1 = m * _kernel->primitive(u, _hsoftinv);
        }
        if(rnrm >= hj) {
            kernel2 = -m / rnrm;
        } else {
            double hjinv = 1. / hj;
            double u = rnrm * hjinv;
            kernel2 = m * _kernel->primitive(u, hjinv);
        }

        kernel = 0.5 * (kernel1 + kernel2);
    }
    _epot += kernel;
    _comp_cost++;
}

/**
 * @brief Action performed on a PseudoNode of the Tree
 *
 * @deprecated This method has been replaced by
 * PotentialWalker::export_to_pseudonode()
 *
 * @param pseudonode PseudoNode on which we operate
 */
void PotentialWalker::pseudonodeaction(PseudoNode* pseudonode) {
    // do nothing
}

/**
 * @brief Decide if the given PseudoNode should be exported to the MPI-process
 * holding the actual node represented by the PseudoNode
 *
 * @param pseudonode PseudoNode on which we operate
 * @return True, since we only approximate the potential for a local TreeNode
 */
bool PotentialWalker::export_to_pseudonode(PseudoNode* pseudonode) {
    return true;
}

/**
 * @brief Set the position for the treewalk
 *
 * @deprecated This method should no longer be used
 *
 * @param position New position for the treewalk
 */
void PotentialWalker::set_position(Vec position) {
    _position = position;
}

/**
 * @brief Get the position used for the treewalk
 *
 * @return Position used for the treewalk
 */
Vec PotentialWalker::get_position() {
    return _position;
}

/**
 * @brief Periodic version of PotentialWalker::splitnode()
 *
 * @deprecated Should no longer be used
 *
 * @param node TreeNode on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 * @return True if the TreeNode should be opened, false otherwise
 */
bool PotentialWalker::periodicsplitnode(TreeNode* node,
                                        EwaldTable& ewald_table) {
    Vec r = node->get_center_of_mass_node() - _position;
    nearest(r);
    double r2 = r.norm2();
    bool openflag = false;
    if(node->get_width() * node->get_width() * node->get_mass_node() >
       r2 * r2 * _olda) {
        openflag = true;
    } else {
        Vec center = node->get_center();
        double nodewidth = 0.6 * node->get_width();
        if(abs(center.x() - _position.x()) < nodewidth) {
            if(abs(center.y() - _position.y()) < nodewidth) {
#if ndim_ == 3
                if(abs(center.z() - _position.z()) < nodewidth) {
                    openflag = true;
                }
#else
                openflag = true;
#endif
            }
        }
    }
    if(openflag) {
        Vec center = node->get_center();
        for(unsigned int i = 0; i < ndim_; i++) {
            double u = center[i] - _position[i];
            nearest(u);
            if(abs(u) >= 0.5 * (_boxsize - node->get_width())) {
                return true;
            }
        }
        if(node->get_width() > 0.2 * _boxsize) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Periodic version of PotentialWalker::pseudonodeaction()
 *
 * @deprecated Method was never implemented and should not be used
 *
 * @param pseudonode PseudoNode on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void PotentialWalker::periodicpseudonodeaction(PseudoNode* pseudonode,
                                               EwaldTable& ewald_table) {}

/**
 * @brief Periodic version of PotentialWalker::leafaction()
 *
 * @deprecated Should not be used anymore
 *
 * @param leaf Leaf on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void PotentialWalker::periodicleafaction(Leaf* leaf, EwaldTable& ewald_table) {
    Vec r = leaf->get_particle()->get_position() - _position;
    nearest(r);
    // not sure if this is right. Pretty sure it is not, since Springel uses
    // a different correction for his potential...
    //    _epot +=
    //    leaf->get_particle()->get_mass()*ewald_table.get_correction(r);
}

/**
 * @brief Finalize the treewalk by setting the gravitational potential of the
 * Particle
 */
void PotentialWalker::after_walk() {
    _p->set_gravitational_potential(_epot);
    _p->add_comp_cost(_comp_cost);
}

/**
 * @brief Get an Export to communicate Particle information to another
 * MPI-process
 *
 * @return An Export that can be sent to another MPI-process
 */
PotentialWalker::Export PotentialWalker::get_export() {
    return Export(_p, _olda);
}

/**
 * @brief Finalize the local treewalk for another MPI-process by setting the
 * data that should be sent back by the Import
 *
 * @param import Import that will be sent back to the original MPI-process
 */
void PotentialWalker::after_walk(Import& import) {
    import.set_epot(_epot);
    import.set_comp_cost(_comp_cost);
}

// Barnes Hut gravity walk

/**
 * @brief Constructor
 *
 * @param p Particle for which to calculate the gravitational acceleration
 * @param local Flag indicating whether the treewalk is done for a Particle
 * residing on the local MPI-process
 */
BHGravityWalker::BHGravityWalker(Particle* p, bool local) {
    _p = p;
    _position = p->get_position();
    _hsoft = p->get_hsoft();
    _hsoftinv = 1. / _hsoft;
    _hsoftinv3 = _hsoftinv * _hsoftinv * _hsoftinv;
    _local = local;
    _comp_cost = 0;
    if(_p->type() == PARTTYPE_GAS) {
        _kernel = new GasGravityKernel();
    } else {
        _kernel = new DMGravityKernel();
    }
}

/**
 * @brief Import constructor
 *
 * Initialize a BHGravityWalker to run a treewalk for a Particle residing on
 * another MPI-process
 *
 * @param import Import containing imported Particle data from another
 * MPI-process
 */
BHGravityWalker::BHGravityWalker(Import& import) {
    _p = NULL;
    _position = import.get_position();
    _hsoft = import.get_hsoft();
    _hsoftinv = 1. / _hsoft;
    _hsoftinv3 = _hsoftinv * _hsoftinv * _hsoftinv;
    _local = false;
    _comp_cost = 0;
    if(import.get_type() == PARTTYPE_GAS) {
        _kernel = new GasGravityKernel();
    } else {
        _kernel = new DMGravityKernel();
    }
}

/**
 * @brief Destructor
 *
 * Free GravityKernel memory.
 */
BHGravityWalker::~BHGravityWalker() {
    delete _kernel;
}

/**
 * @brief Get the resulting gravitational acceleration
 *
 * @return Result of the treewalk
 */
Vec BHGravityWalker::get_acceleration() {
    return _a;
}

/**
 * @brief Decide whether to open the given TreeNode or treat it as a whole
 *
 * When the TreeNode is treated as a whole, the relevant computations are done
 * inside this function.
 *
 * @param node TreeNode on which to operate
 * @return True if the TreeNode should be opened, false otherwise
 */
bool BHGravityWalker::splitnode(TreeNode* node) {
    Vec r = node->get_center_of_mass_node() - _position;
    if(_boxsize) {
        nearest(r);
    }
    double r2 = r.norm2();
    // opening angle = 0.5 (0.5^2 = 0.25)
    if(node->get_width() * node->get_width() <= r2 * 0.25) {
        // if the node is not local, it means that the same treenode is also in
        // the tree on the original process
        // since it is not opened here, it was not opened there and the force
        // was already calculated on the original process
        double kernel;
        double rnrm = sqrt(r2);
        if(rnrm >= _hsoft && rnrm >= node->get_hmax()) {
            kernel = node->get_mass_node() / (rnrm * r2);
        } else {
            // if the node is within the softening length, open it up
            // this way, we never calculate softened forces between particles
            // and nodes
            return true;
        }
        if(_local || node->is_local()) {
            _a += r * kernel;
            _comp_cost++;
        }
        return false;
    } else {
        return true;
    }
}

/**
 * @brief Action performed on a single TreeNode of the Tree
 *
 * This method does nothing, since the relevant action is already performed by
 * BHGravityWalker::splitnode().
 *
 * @param node TreeNode on which to operate
 */
void BHGravityWalker::nodeaction(TreeNode* node) {
    // obsolete, action is perfomed by splitnode()
}

/**
 * @brief Action performed on a single Leaf of the Tree
 *
 * @param leaf Leaf on which to operate
 */
void BHGravityWalker::leafaction(Leaf* leaf) {
    // version with variable softening lengths
    Particle* pj = leaf->get_particle();
    Vec r = pj->get_position() - _position;
    if(_boxsize) {
        nearest(r);
    }
    double m = pj->get_mass();
    // the actual softening is 2.8 times the stored value (see Gadget2)
    // it probably is cleaner to change this, but for now we stick to it
    double hj = pj->get_hsoft();
    double kernel;
    // we need both r2 and rnrm for all cases
    double r2 = r.norm2();
    double rnrm = sqrt(r2);
    if(rnrm >= _hsoft && rnrm >= hj) {
        kernel = m / (rnrm * r2);
    } else {
        // version where we average the contributions of both particles
        double kernel1, kernel2;
        if(rnrm >= _hsoft) {
            kernel1 = m / (rnrm * r2);
        } else {
            double u = rnrm * _hsoftinv;
            kernel1 = m * (*_kernel)(u, _hsoftinv3);
        }
        if(rnrm >= hj) {
            kernel2 = m / (rnrm * r2);
        } else {
            double hjinv = 1. / hj;
            double hjinv2 = hjinv * hjinv;
            double hjinv3 = hjinv2 * hjinv;
            double u = rnrm * hjinv;
            kernel2 = m * (*_kernel)(u, hjinv3);
        }
        kernel = 0.5 * (kernel1 + kernel2);
    }
    _a += r * kernel;
    _comp_cost++;
}

/**
 * @brief Action to perform on a PseudoNode of the Tree
 *
 * @deprecated This method is replaced by BHGravityWalker::export_to_pseudonode
 *
 * @param pseudonode PseudoNode on which to operate
 */
void BHGravityWalker::pseudonodeaction(PseudoNode* pseudonode) {
    // do nothing
}

/**
 * @brief Decide whether the Particle should be exported to the MPI-process
 * containing the actual node represented by this PseudoNode
 *
 * @param pseudonode PseudoNode on which to operate
 * @return True, since the gravitational acceleration is only approximated for
 * a local TreeNode
 */
bool BHGravityWalker::export_to_pseudonode(PseudoNode* pseudonode) {
    return true;
}

/**
 * @brief Set the position used for the treewalk
 *
 * @deprecated This method should no longer be used
 *
 * @param position New value for the position
 */
void BHGravityWalker::set_position(Vec position) {
    _position = position;
}

/**
 * @brief Get the position used for the treewalk
 *
 * @return Position used for the treewalk
 */
Vec BHGravityWalker::get_position() {
    return _position;
}

/**
 * @brief Periodic version of BHGravityWalker::splitnode
 *
 * @deprecated No longer used
 *
 * @param node TreeNode on which to operate
 * @param ewald_table EwaldTable used for periodic correction terms
 * @return
 */
bool BHGravityWalker::periodicsplitnode(TreeNode* node,
                                        EwaldTable& ewald_table) {
    Vec r = node->get_center_of_mass_node() - _position;
    nearest(r);
    double r2 = r.norm2();
    // opening angle = 0.5 (0.5^2 = 0.25)
    if(node->get_width() * node->get_width() > r2 * 0.25) {
        Vec center = node->get_center();
        for(unsigned int i = 0; i < ndim_; i++) {
            double u = center[i] - _position[i];
            nearest(u);
            if(abs(u) >= 0.5 * (_boxsize - node->get_width())) {
                return true;
            }
        }
        if(node->get_width() > 0.2 * _boxsize) {
            return true;
        }
    }
    _a += node->get_mass_node() * ewald_table.get_correction(r);
    return false;
}

/**
 * @brief Periodic version of BHGravityWalker::pseudonodeaction
 *
 * @deprecated Never implemented, should no longer be used
 *
 * @param pseudonode PseudoNode on which to operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void BHGravityWalker::periodicpseudonodeaction(PseudoNode* pseudonode,
                                               EwaldTable& ewald_table) {
    // do nothing
}

/**
 * @brief Periodic version of BHGravityWalker::leafaction
 *
 * @deprecated This method is no longer used
 *
 * @param leaf Leaf on which to operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void BHGravityWalker::periodicleafaction(Leaf* leaf, EwaldTable& ewald_table) {
    Vec r = leaf->get_particle()->get_position() - _position;
    nearest(r);
    _a += leaf->get_particle()->get_mass() * ewald_table.get_correction(r);
}

/**
 * @brief Finalize the treewalk by setting the gravitational acceleration of the
 * Particle
 */
void BHGravityWalker::after_walk() {
    _p->set_gravitational_acceleration(_a);
    _p->add_comp_cost(_comp_cost);
}

/**
 * @brief Get an Export to communicate Particle data to another MPI-process
 *
 * @return An Export that can be sent to another MPI-process
 */
BHGravityWalker::Export BHGravityWalker::get_export() {
    return Export(_p);
}

/**
 * @brief Finalize the local treewalk for another MPI-process by setting the
 * acceleration of the Import
 *
 * @param import Import that will be sent back to the original MPI-process
 */
void BHGravityWalker::after_walk(Import& import) {
    import.set_a(_a);
    import.set_comp_cost(_comp_cost);
}

// BHEwaldGravityWalker

/**
 * @brief Constructor
 *
 * @param p Particle for which the Ewald correction force is calculated
 * @param ewald_table Reference to the EwaldTable used to calculate Ewald
 * corrections
 * @param local Flag signaling if the tree walk is done for a particle on the
 * local process
 */
BHEwaldGravityWalker::BHEwaldGravityWalker(Particle* p, EwaldTable& ewald_table,
                                           bool local)
        : _ewald_table(ewald_table) {
    _p = p;
    _position = p->get_position();
}

/**
 * @brief Import constructor.
 *
 * Initialize an EwaldGravityWalker for a Particle residing on another
 * MPI-process
 *
 * @param import Import with imported data
 * @param ewald_table Reference to the EwaldTable used to calculate Ewald
 * corrections
 */
BHEwaldGravityWalker::BHEwaldGravityWalker(Import& import,
                                           EwaldTable& ewald_table)
        : _ewald_table(ewald_table) {
    _p = NULL;
    _position = import.get_position();
}

/**
 * @brief Destructor
 */
BHEwaldGravityWalker::~BHEwaldGravityWalker() {}

/**
 * @brief Get the resulting Ewald correction
 *
 * @return The result of the treewalk
 */
Vec BHEwaldGravityWalker::get_acceleration() {
    return _acorr;
}

/**
 * @brief Decide whether the given TreeNode should be opened or can be treated
 * as a whole
 *
 * If it can be treated as a whole, we immediately treat it, as this makes more
 * efficient use of the intermediate quantities.
 *
 * @param node TreeNode on which we operate
 * @return True if the TreeNode should be opened, false otherwise
 */
bool BHEwaldGravityWalker::splitnode(TreeNode* node) {
    Vec r = node->get_center_of_mass_node() - _position;
    nearest(r);
    double r2 = r.norm2();
    // opening angle = 0.5 (0.5^2 = 0.25)
    if(node->get_width() * node->get_width() > r2 * 0.25) {
        Vec center = node->get_center();
        for(unsigned int i = 0; i < ndim_; i++) {
            double u = center[i] - _position[i];
            nearest(u);
            if(abs(u) >= 0.5 * (_boxsize - node->get_width())) {
                return true;
            }
        }
        if(node->get_width() > 0.2 * _boxsize) {
            return true;
        }
    }
    _acorr += node->get_mass_node() * _ewald_table.get_correction(r);
    return false;
}

/**
 * @brief Action to perform on a TreeNode as a whole
 *
 * @deprecated Action is already performed by EwaldGravityWalker::splitnode()
 *
 * @param node TreeNode on which we operate
 */
void BHEwaldGravityWalker::nodeaction(TreeNode* node) {
    // obsolete, action is perfomed by splitnode()
}

/**
 * @brief Action to perform on a single Leaf of the Tree
 *
 * @param leaf Leaf on which we operate
 */
void BHEwaldGravityWalker::leafaction(Leaf* leaf) {
    Vec r = leaf->get_particle()->get_position() - _position;
    nearest(r);
    _acorr += leaf->get_particle()->get_mass() * _ewald_table.get_correction(r);
}

/**
 * @brief Action performed when a PseudoNode is encountered
 *
 * @deprecated Method has been replaced by
 * EwaldGravityWalker::export_to_pseudonode
 *
 * @param pseudonode PseudoNode on which we operate
 */
void BHEwaldGravityWalker::pseudonodeaction(PseudoNode* pseudonode) {
    // do something
    //    _exportlist[pseudonode->get_source()] = true;
}

/**
 * @brief Decide if the Particle should be exported to the process which holds
 * the actual node represented by the given PseudoNode
 *
 * @param pseudonode PseudoNode on which we operate
 * @return True, since the gravitational force is only approximated for a local
 * TreeNode
 */
bool BHEwaldGravityWalker::export_to_pseudonode(PseudoNode* pseudonode) {
    return true;
}

/**
 * @brief Set the position which is used for the treewalk
 *
 * @deprecated This method should not be used anymore
 *
 * @param position New position for the treewalk
 */
void BHEwaldGravityWalker::set_position(Vec position) {
    _position = position;
}

/**
 * @brief Get the position used for the treewalk
 *
 * @return Position used for the treewalk
 */
Vec BHEwaldGravityWalker::get_position() {
    return _position;
}

/**
 * @brief Periodic version of EwaldGravityWalker::splitnode()
 *
 * @deprecated Periodic boundaries are treated in a different way now
 *
 * @param node TreeNode on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 * @return True if the TreeNode should be opened, false otherwise
 */
bool BHEwaldGravityWalker::periodicsplitnode(TreeNode* node,
                                             EwaldTable& ewald_table) {
    return false;
}

/**
 * @brief Periodic version of EwaldGravityWalker::pseudonodeaction()
 *
 * @deprecated Method was never implemented and should not be used, since
 * periodic boundaries are treated differently now
 *
 * @param pseudonode PseudoNode on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void BHEwaldGravityWalker::periodicpseudonodeaction(PseudoNode* pseudonode,
                                                    EwaldTable& ewald_table) {}

/**
 * @brief Periodic version of EwaldGravityWalker::leafaction()
 *
 * @deprecated Not used anymore
 *
 * @param leaf Leaf on which we operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void BHEwaldGravityWalker::periodicleafaction(Leaf* leaf,
                                              EwaldTable& ewald_table) {}

/**
 * @brief Finalize the treewalk by setting the gravitational acceleration of the
 * Particle
 */
void BHEwaldGravityWalker::after_walk() {
    _p->set_gravitational_acceleration(_p->get_gravitational_acceleration() +
                                       _acorr);
}

/**
 * @brief Get an Export to export the Particle to another MPI-process
 *
 * @return An Export that can be used to communicate relevant Particle data to
 * another MPI-process
 */
BHEwaldGravityWalker::Export BHEwaldGravityWalker::get_export() {
    return Export(_p, 0.);
}

/**
 * @brief Finalize the local treewalk for another MPI-process by setting the
 * acceleration of the Import
 *
 * @param import Import that will be sent back to the original MPI-process
 */
void BHEwaldGravityWalker::after_walk(Import& import) {
    import.set_acorr(_acorr);
}

// Barnes-Hut potential walk

/**
 * @brief Constructor
 *
 * @param p Particle for which the gravitational potential is calculated
 * @param local Flag indicating if the treewalk is done for a Particle
 * residing on the local MPI-process
 */
BHPotentialWalker::BHPotentialWalker(Particle* p, bool local) {
    _p = p;
    _position = p->get_position();
    _hsoft = p->get_hsoft();
    _hsoftinv = 1. / _hsoft;
    _hsoftinv3 = _hsoftinv * _hsoftinv * _hsoftinv;
    _local = local;
    _comp_cost = 0;
    if(_p->type() == PARTTYPE_GAS) {
        _kernel = new GasGravityKernel();
    } else {
        _kernel = new DMGravityKernel();
    }

    _epot = 0.;
}

/**
 * @brief Import constructor
 *
 * Initialize a BHPotentialWalker for a Particle residing on another
 * MPI-process.
 *
 * @param import Import containing imported Particle data
 */
BHPotentialWalker::BHPotentialWalker(Import& import) {
    _p = NULL;
    _position = import.get_position();
    _hsoft = import.get_hsoft();
    _hsoftinv = 1. / _hsoft;
    _hsoftinv3 = _hsoftinv * _hsoftinv * _hsoftinv;
    _local = false;
    _comp_cost = 0;
    if(import.get_type() == PARTTYPE_GAS) {
        _kernel = new GasGravityKernel();
    } else {
        _kernel = new DMGravityKernel();
    }

    _epot = 0.;
}

/**
 * @brief Destructor
 *
 * Free GravityKernel memory.
 */
BHPotentialWalker::~BHPotentialWalker() {
    delete _kernel;
}

/**
 * @brief Get the resulting gravitational potential
 *
 * @return Result of the treewalk
 */
double BHPotentialWalker::get_epot() {
    return _epot;
}

/**
 * @brief Decide whether to split the given TreeNode or treat it as a whole
 *
 * If it is treated as a whole, this is done inside this function.
 *
 * @param node TreeNode on which to operate
 * @return True if the TreeNode should be opened, false otherwise
 */
bool BHPotentialWalker::splitnode(TreeNode* node) {
    Vec r = node->get_center_of_mass_node() - _position;
    if(_boxsize) {
        nearest(r);
    }
    double r2 = r.norm2();
    // opening angle = 0.5 (0.5^2 = 0.25)
    if(node->get_width() * node->get_width() <= r2 * 0.25) {
        // if the node is not local, it means that the same treenode is also in
        // the tree on the original process
        // since it is not opened here, it was not opened there and the force
        // was already calculated on the original process
        double kernel;
        double rnrm = sqrt(r2);
        if(rnrm >= _hsoft && rnrm >= node->get_hmax()) {
            kernel = -node->get_mass_node() / rnrm;
        } else {
            // if the node is within the softening length, open it up
            // this way, we never calculate softened forces between particles
            // and nodes
            return true;
        }
        if(_local || node->is_local()) {
            _epot += kernel;
            _comp_cost++;
        }
        return false;
    } else {
        return true;
    }
}

/**
 * @brief Action to perform on a single TreeNode of the Tree
 *
 * Does nothing, since the relevant action is done by
 * BHPotentialWalker::splitnode().
 *
 * @param node TreeNode on which to operate
 */
void BHPotentialWalker::nodeaction(TreeNode* node) {
    // obsolete, action is perfomed by splitnode()
}

/**
 * @brief Action to perform on a single Leaf of the Tree
 *
 * @param leaf Leaf on which to operate
 */
void BHPotentialWalker::leafaction(Leaf* leaf) {
    // version with variable softening lengths
    Particle* pj = leaf->get_particle();
    Vec r = pj->get_position() - _position;
    if(_boxsize) {
        nearest(r);
    }
    double m = pj->get_mass();
    // the actual softening is 2.8 times the stored value (see Gadget2)
    // it probably is cleaner to change this, but for now we stick to it (no,
    // we're changing it (finally))
    double hj = pj->get_hsoft();
    double kernel;
    // we need both r2 and rnrm for all cases
    double r2 = r.norm2();
    double rnrm = sqrt(r2);
    if(rnrm >= _hsoft && rnrm >= hj) {
        kernel = -m / rnrm;
    } else {
        // version where we average the contributions of both particles
        double kernel1, kernel2;
        if(rnrm >= _hsoft) {
            kernel1 = -m / rnrm;
        } else {
            double u = rnrm * _hsoftinv;
            kernel1 = m * _kernel->primitive(u, _hsoftinv);
        }
        if(rnrm >= hj) {
            kernel2 = -m / rnrm;
        } else {
            double hjinv = 1. / hj;
            double u = rnrm * hjinv;
            kernel2 = m * _kernel->primitive(u, hjinv);
        }

        kernel = 0.5 * (kernel1 + kernel2);
    }
    _epot += kernel;
    _comp_cost++;
}

/**
 * @brief Action to perform when a PseudoNode is encountered
 *
 * @deprecated Replaced by BHPotentialWalker::export_to_pseudonode()
 *
 * @param pseudonode PseudoNode on which to operate
 */
void BHPotentialWalker::pseudonodeaction(PseudoNode* pseudonode) {
    // nothing
}

/**
 * @brief Decide whether the Particle should be exported to the MPI-process
 * holding the actual node represented by the given PseudoNode
 *
 * @param pseudonode PseudoNode on which to operate
 * @return True, since the potential is only approximated for a local TreeNode
 */
bool BHPotentialWalker::export_to_pseudonode(PseudoNode* pseudonode) {
    return true;
}

/**
 * @brief Set the position used for the treewalk
 *
 * @deprecated Should no longer be used
 *
 * @param position New position used for the treewalk
 */
void BHPotentialWalker::set_position(Vec position) {
    _position = position;
}

/**
 * @brief Get the position used for the treewalk
 *
 * @return Position used for the treewalk
 */
Vec BHPotentialWalker::get_position() {
    return _position;
}

/**
 * @brief Periodic version of BHPotentialWalker::splitnode()
 *
 * @deprecated No longer used
 *
 * @param node TreeNode on which to operate
 * @param ewald_table EwaldTable used for periodic correction terms
 * @return True if the TreeNode should be opened, false otherwise
 */
bool BHPotentialWalker::periodicsplitnode(TreeNode* node,
                                          EwaldTable& ewald_table) {
    Vec r = node->get_center_of_mass_node() - _position;
    nearest(r);
    double r2 = r.norm2();
    // opening angle = 0.5 (0.5^2 = 0.25)
    if(node->get_width() * node->get_width() > r2 * 0.25) {
        Vec center = node->get_center();
        for(unsigned int i = 0; i < ndim_; i++) {
            double u = center[i] - _position[i];
            nearest(u);
            if(abs(u) >= 0.5 * (_boxsize - node->get_width())) {
                return true;
            }
        }
        if(node->get_width() > 0.2 * _boxsize) {
            return true;
        }
    }
    //    _a += node->get_mass_node()*ewald_table.get_correction(r);
    return false;
}

/**
 * @brief Periodic version of BHPotentialWalker::pseudonodeaction()
 *
 * @deprecated Not implemented, should not be used
 *
 * @param pseudonode PseudoNode on which to operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void BHPotentialWalker::periodicpseudonodeaction(PseudoNode* pseudonode,
                                                 EwaldTable& ewald_table) {
    // do nothing
}

/**
 * @brief Periodic version of BHPotentialWalker::leafaction()
 *
 * @deprecated No longer used
 *
 * @param leaf Leaf on which to operate
 * @param ewald_table EwaldTable used for periodic correction terms
 */
void BHPotentialWalker::periodicleafaction(Leaf* leaf,
                                           EwaldTable& ewald_table) {
    Vec r = leaf->get_particle()->get_position() - _position;
    nearest(r);
    //    _a += leaf->get_particle()->get_mass()*ewald_table.get_correction(r);
}

/**
 * @brief Finalize the treewalk by setting the gravitational potential of the
 * Particle
 */
void BHPotentialWalker::after_walk() {
    _p->set_gravitational_potential(_epot);
    _p->add_comp_cost(_comp_cost);
}

/**
 * @brief Get an Export to communicate Particle information to another
 * MPI-process
 *
 * @return An Export that can be exported to another MPI-process
 */
BHPotentialWalker::Export BHPotentialWalker::get_export() {
    return Export(_p);
}

/**
 * @brief Finalize the local treewalk for another MPI-process by setting the
 * gravitational potential of the given Import
 *
 * @param import Import that will be sent back to the original MPI-process
 */
void BHPotentialWalker::after_walk(Import& import) {
    import.set_epot(_epot);
    import.set_comp_cost(_comp_cost);
}
