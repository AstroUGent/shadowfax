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
 * @file GasParticle.cpp
 *
 * @brief Gas particle: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "GasParticle.hpp"
#include "RestartFile.hpp"
#include "VorGen.hpp"
#include <iostream>  // for operator<<, ostream, etc
#include <stddef.h>  // for NULL
using namespace std;

/**
 * @brief Empty constructor
 *
 * Initialize an empty GasParticle with zero properties.
 */
GasParticle::GasParticle() : Particle() {
    _csnd = 0.;
    _vorgen = NULL;
    _vorgenid = 0;
    _h = 0;
    _copies.push_back(0);
    _exports = new bool[MPIGlobal::size];
    for(int i = 0; i < MPIGlobal::size; i++) {
        _exports[i] = false;
    }
    _max_radius = 0.;
    _is_pseudo = false;
    _max_mach = 0.;
    _local_id = 0;
    _total_area = 0.;
    _real_dt = 0.;
}

/**
 * @brief Constructor
 *
 * @param pos Position of the particle
 */
GasParticle::GasParticle(Vec pos) : Particle(pos) {
    _csnd = 0.;
    _vorgen = NULL;
    _vorgenid = 0;
    _h = 0;
    _copies.push_back(0);
    _exports = new bool[MPIGlobal::size];
    for(int i = 0; i < MPIGlobal::size; i++) {
        _exports[i] = false;
    }
    _max_radius = 0.;
    _is_pseudo = false;
    _max_mach = 0.;
    _local_id = 0;
    _total_area = 0.;
    _real_dt = 0.;
}

/**
 * @brief Copy constructor
 *
 * Since we store a raw pointer to a custom allocated buffer (_exports), we
 * cannot rely on the default copy constructor, which just copies everything.
 * If we would use the default constructor, then _exports would point to the
 * same memory for two different particles. This in itself is already dangerous,
 * but even more problematic is the fact that _exports will then effectively be
 * deleted twice by the destructors of the two GasParticle instances, which will
 * cause a segfault.
 *
 * We have to manually assign new memory for the _exports array and then copy
 * all values to the newly created buffer.
 *
 * @param p Reference to the GasParticle that needs to be copied
 */
GasParticle::GasParticle(GasParticle& p) : Particle(p) {
    _a_grav_old = p._a_grav_old;
    _csnd = p._csnd;
    _vorgen = p._vorgen;
    _vorgenid = p._vorgenid;
    _old_Q = p._old_Q;
    _h = p._h;
    _copies = p._copies;
    // the whole reason we declare a copy constructor: we cannot simply copy
    // _exports, we need to assign new memory and copy all elements instead
    _exports = new bool[MPIGlobal::size];
    for(int i = 0; i < MPIGlobal::size; i++) {
        _exports[i] = p._exports[i];
    }
    _W = p._W;
    _Q = p._Q;
    _dQ = p._dQ;
    _delta_E_grav = p._delta_E_grav;
    for(unsigned int i = 0; i < ndim_; i++) {
        _gradientvecs[i] = p._gradientvecs[i];
    }
    _max_radius = p._max_radius;
    _is_pseudo = p._is_pseudo;
    _mesh_v = p._mesh_v;
    _mesh_m = p._mesh_m;
    _max_mach = p._max_mach;
    _centroid = p._centroid;
    _local_id = p._local_id;
    _eta = p._eta;
    _old_v = p._old_v;
    _real_dt = p._real_dt;
    _total_area = p._total_area;
}

/**
 * @brief Destructor
 *
 * Free exports buffer.
 */
GasParticle::~GasParticle() {
    delete[] _exports;
}

/**
 * @brief Set the soundspeed of the particle to the given value
 *
 * @param csnd New soundspeed for the particle
 */
void GasParticle::set_soundspeed(double csnd) {
    _csnd = csnd;
}

/**
 * @brief Get the soundspeed of the particle
 *
 * @return Soundspeed of the particle
 */
double GasParticle::get_soundspeed() {
    return _csnd;
}

/**
 * @brief Set the VorGen corresponding to this GasParticle
 *
 * @param vorgen The VorGen corresponding to this GasParticle
 */
void GasParticle::set_vorgen(VorGen* vorgen) {
    _vorgen = vorgen;
}

/**
 * @brief Get the VorGen corresponding to this GasParticle
 *
 * @return The VorGen corresponding to this GasParticle
 */
VorGen* GasParticle::get_vorgen() {
    return _vorgen;
}

/**
 * @brief Get the VorCell corresponding to this GasParticle
 *
 * @return The VorCell corresponding to this GasParticle
 */
VorCell* GasParticle::get_cell() {
    return _vorgen->get_cell();
}

/**
 * @brief Set the index of the VorGen corresponding to this GasParticle in the
 * DelTess VorGen list
 *
 * @param vorgenid The index of the VorGen corresponding to this GasParticle
 */
void GasParticle::set_vorgenid(unsigned int vorgenid) {
    _vorgenid = vorgenid;
}

/**
 * @brief Get the index of the VorGen corresponding to this GasParticle in the
 * DelTess VorGen list
 *
 * @return The index of the VorGen corresponding to this GasParticle
 */
unsigned int GasParticle::get_vorgenid() {
    return _vorgenid;
}

/**
 * @brief Set the primitive quantities of this particle
 *
 * @param W StateVector containing the primitive quantities for this particle
 */
void GasParticle::set_W(StateVector W) {
    _W = W;
}

/**
 * @brief Get the primitive quantities of this particle
 *
 * @return StateVector containing the primitive quantities for this particle
 */
const StateVector& GasParticle::get_Wvec() {
    return _W;
}

/**
 * @brief Set the conserved quantities of this particle
 *
 * @param Q StateVector containing the conserved quantities for this particle
 */
void GasParticle::set_Q(StateVector Q) {
    _Q = Q;
}

/**
 * @brief Get the conserved quantities of this particle
 *
 * @return StateVector containing the conserved quantities for this particle
 */
StateVector& GasParticle::get_Qvec() {
    return _Q;
}

/**
 * @brief Get the characteristic length of the particle
 *
 * @return The characteristic length of the particle
 */
double GasParticle::h() {
    return _h;
}

/**
 * @brief Set the characteristic length of the particle
 *
 * @param h The new characteristic length for the particle
 */
void GasParticle::set_h(double h) {
    _h = h;
}

/**
 * @brief Get the total surface area of the cell associated with this particle
 *
 * @return Total surface area of the Voronoi cell
 */
double GasParticle::get_total_area() {
    return _total_area;
}

/**
 * @brief Set the total surface area of the cell associated with this particle
 *
 * @param total_area Total surface area of the Voronoi cell
 */
void GasParticle::set_total_area(double total_area) {
    _total_area = total_area;
}

/**
 * @brief Add the gravitational force to the conserved quantities of the
 * particle, using the given timestep
 *
 * @param dt Real timestep of the simulation
 */
void GasParticle::apply_gravity(double dt) {
#define SCHEME3
    // don't add gravity to vacuum cells
    if(!_Q[0]) {
        return;
    }

#ifdef SCHEME1
#if ndim_ == 3
    _Q[4] -= 0.5 * (_Q[1] * _Q[1] + _Q[2] * _Q[2] + _Q[3] * _Q[3]) / _Q[0];
#else
    _Q[3] -= 0.5 * (_Q[1] * _Q[1] + _Q[2] * _Q[2]) / _Q[0];
#endif
    _Q += dt * _Q[0] * _a_grav_new;
#if ndim_ == 3
    _Q[4] += 0.5 * (_Q[1] * _Q[1] + _Q[2] * _Q[2] + _Q[3] * _Q[3]) / _Q[0];
#else
    _Q[3] += 0.5 * (_Q[1] * _Q[1] + _Q[2] * _Q[2]) / _Q[0];
#endif

#endif

#ifdef SCHEME2

#if ndim_ == 3
    _Q[4] += 0.5 * dt * (_Q[1] * _a_grav_new[0] + _Q[2] * _a_grav_new[1] +
                         _Q[3] * _a_grav_new[2]);
#else
    _Q[3] += 0.5 * dt * (_Q[1] * _a_grav_new[0] + _Q[2] * _a_grav_new[1]);
#endif
    _Q += dt * _Q[0] * _a_grav_new;
#if ndim_ == 3
    _Q[4] += 0.5 * dt * (_Q[1] * _a_grav_new[0] + _Q[2] * _a_grav_new[1] +
                         _Q[3] * _a_grav_new[2]);
#else
    _Q[3] += 0.5 * dt * (_Q[1] * _a_grav_new[0] + _Q[2] * _a_grav_new[1]);
#endif

#endif

#ifdef SCHEME3
    _Q[ndim_ + 1] += dt * _Q[0] * Vec::dot(_v, _a_grav_new);

    _Q += dt * _Q[0] * _a_grav_new;

    _Q[ndim_ + 1] -= 0.5 * Vec::dot(_delta_E_grav, _a_grav_new);
#if ndim_ == 3
    _delta_E_grav.set(0., 0., 0.);
#else
    _delta_E_grav.set(0., 0.);
#endif
#endif
}

/**
 * @brief Make an internal backup of the conserved quantities
 */
void GasParticle::save_Q() {
    _old_Q = _Q;
}

/**
 * @brief Set the mass of the particle
 *
 * The mass is the first component of the conserved variables StateVector.
 *
 * @param mass New mass for the particle
 */
void GasParticle::set_mass(double mass) {
    _Q[0] = mass;
}

/**
 * @brief Set the mesh velocity for this particle
 *
 * @deprecated This quantity is not used
 *
 * @param mesh_v New mesh velocity for this particle
 */
void GasParticle::set_mesh_v(Vec& mesh_v) {
    _mesh_v = mesh_v;
}

/**
 * @brief Get the mesh velocity for this particle
 *
 * @deprecated This quantity is not used
 *
 * @return Mesh velocity for this particle
 */
Vec& GasParticle::get_mesh_v() {
    return _mesh_v;
}

/**
 * @brief Set the mesh mass for this particle
 *
 * @deprecated This quantity is not used
 *
 * @param mesh_m New mesh mass for this particle
 */
void GasParticle::set_mesh_m(double mesh_m) {
    _mesh_m = mesh_m;
}

/**
 * @brief Get the mesh mass for this particle
 *
 * @deprecated This quantity is not used
 *
 * @return Mesh mass for this particle
 */
double GasParticle::get_mesh_m() {
    return _mesh_m;
}

/**
 * @brief Flag this GasParticle as a pseudoparticle
 */
void GasParticle::set_pseudo() {
    _is_pseudo = true;
}

/**
 * @brief Check if this GasParticle is a pseudoparticle
 *
 * @return True if this particle is a pseudoparticle, false otherwise
 */
bool GasParticle::is_pseudo() {
    return _is_pseudo;
}

/**
 * @brief Signal that the given copy of the particle was made on the MPI process
 * with the given id
 *
 * @param id unsigned integer flag for the boundary at which the copy was made
 * @param index Rank of the process on which the copy was made (I think)
 */
void GasParticle::make_copy(unsigned int id, unsigned int index) {
    while(index >= _copies.size()) {
        _copies.push_back(0);
    }
    _copies[index] += (1 << id);
}

/**
 * @brief Check if the given copy of the particle was made on the given MPI
 * process
 *
 * @param id unsigned integer flag for the boundary at which the copy was made
 * @param index Rank of the process on which the copy was made (I think)
 * @return True if the given copy exists, false otherwise
 */
bool GasParticle::is_copied(unsigned int id, unsigned int index) {
    bool is_copied = (index < _copies.size());
    if(is_copied) {
        is_copied = (1 & (_copies[index] >> id));
    } else {
        is_copied = false;
    }
    return is_copied;
}

/**
 * @brief Reset copy information
 */
void GasParticle::reset_copies() {
    _copies.clear();
    _copies.push_back(0);
}

/**
 * @brief Signal an export of this particle to the given MPI process
 *
 * @param id Rank of the MPI process to which the particle is exported
 */
void GasParticle::do_export(unsigned int id) {
    _exports[id] = true;
}

/**
 * @brief Check if the particle was exported to the given MPI process
 *
 * @param id Rank of the MPI process
 * @return True if the particle was exported, false otherwise
 */
bool GasParticle::is_exported(unsigned int id) {
    return _exports[id];
}

/**
 * @brief Reset export information
 */
void GasParticle::reset_export() {
    for(int i = 0; i < MPIGlobal::size; i++) {
        _exports[i] = false;
    }
}

/**
 * @brief Set the maximal search radius for the particle
 *
 * @param max_radius New maximal search radius for the particle
 */
void GasParticle::set_max_radius(double max_radius) {
    _max_radius = max_radius;
}

/**
 * @brief Get the maximal search radius for the particle
 *
 * @return Maximal search radius for the particle
 */
double GasParticle::get_max_radius() {
    return _max_radius;
}

/**
 * @brief Update the primitive variables of this particle
 *
 * @deprecated This method has been replaced by a better method
 *
 * @param gamma Adiabatic index of the fluid
 */
void GasParticle::update_W(double gamma) {
    // do nothing
}

/**
 * @brief Get the given component of the gradient in the given coordinate
 * direction
 *
 * @param Windex Index of the primitive quantity component
 * @param dimindex Coordinate direction
 * @return The gradient for the given component in the given direction
 */
double GasParticle::get_gradient(unsigned int Windex, unsigned int dimindex) {
    return _gradientvecs[dimindex][Windex];
}

/**
 * @brief Set the gradients for this particle
 *
 * @param gradients Gradients for this particle
 */
void GasParticle::set_gradients(StateVector* gradients) {
    for(unsigned int i = ndim_; i--;) {
        _gradientvecs[i] = gradients[i];
    }
}

/**
 * @brief Get the gradients of this particle
 *
 * @param gradients Gradients of this particle
 */
void GasParticle::get_gradients(StateVector* gradients) {
    for(unsigned int i = ndim_; i--;) {
        gradients[i] = _gradientvecs[i];
    }
}

/**
 * @brief Increase the conserved quantities of this particle with the given
 * values
 *
 * @param dQ Conserved quantity increase
 */
void GasParticle::increase_dQ(StateVector dQ) {
    _dQ += dQ;
}

/**
 * @brief Get the increase in conserved quantities during the current timestep
 *
 * @return Increase in conserved quantities
 */
StateVector& GasParticle::get_dQvec() {
    return _dQ;
}

/**
 * @brief Reset the increase in conserved quantities during the current timestep
 */
void GasParticle::reset_dQ() {
    _dQ.reset();
}

/**
 * @brief Increase the gravitational energy difference during the current
 * timestep with the given work term
 *
 * @param dE Gravitational work term
 */
void GasParticle::increase_delta_E(Vec dE) {
    _delta_E_grav += dE;
}

/**
 * @brief Add the gravitational work term of the cell
 *
 * @deprecated This method is no longer used
 *
 * @param dt Real timestep of the simulation
 */
void GasParticle::add_dE_grav_cell(double dt) {
    // do nothing
}

/**
 * @brief Reset the gravitational work term during the current timestep
 */
void GasParticle::reset_dE_grav() {
#if ndim_ == 3
    _delta_E_grav.set(0., 0., 0.);
#else
    _delta_E_grav.set(0., 0.);
#endif
}

/**
 * @brief Get the gravitational work term during the current timestep
 *
 * @return Gravitational work term
 */
Vec GasParticle::get_delta_E() {
    return _delta_E_grav;
}

/**
 * @brief Update the conserved quantities of the particle by adding the increase
 * during the current timestep to the actual conserved quantities
 *
 * Errors are thrown if this results in negative mass or energy.
 */
void GasParticle::update_Q() {
    _Q -= _dQ;
    if(_Q.m() < 0. || _Q.e() < 0.) {
        throw UpdateQException(*this);
    }

    _dQ.reset();
}

/**
 * @brief Set the Passively Advected Quantity (PAQ) for this particle
 *
 * @param index Deprecated parameter
 * @param paq New value for the PAQ of this particle
 */
void GasParticle::set_paq(unsigned int index, double paq) {
    _W.set_paq(paq);
}

///**
// * @brief Save the particle properties to the given Block for output
// *
// * @param block Block to write to
// */
// void GasParticle::save_props(Block& block) {
//    vector<double> props;

//    props.push_back(_id);
//    // sometimes, order DOES matter
//    for(unsigned int i = 0; i < ndim_ + 2; i++) {
//        props.push_back(_W[i]);
//    }
//    props.push_back(_Q[0]);

//    block.add_data(props);
//}

///**
// * @brief Save the mesh information of this particle to the given Block for
// * output
// *
// * @param block Block to write to
// */
// void GasParticle::save_grid(Block& block) {
//    vector<double> pos;
//    pos.push_back(_x[0]);
//    pos.push_back(_x[1]);
//#if ndim_ == 3
//    pos.push_back(_x[2]);
//#else
//    pos.push_back(0.);
//#endif
//    block.add_data(pos);
//}

///**
// * @brief Save the particle properties to the given Block for Gadget output
// *
// * @deprecated No longer used
// *
// * @param block Block to write to
// */
// void GasParticle::save_props_gadget(Block& block) {
//    vector<double> props;
//    props.push_back(_x[0]);
//    props.push_back(_x[1]);
//#if ndim_ == 3
//    props.push_back(_x[2]);
//#endif

//    props.push_back(_W.rho());
//    props.push_back(0.);
//    props.push_back(0.);
//    props.push_back(_id);
//    props.push_back(_W.p());
//    props.push_back(0.);
//    props.push_back(0.);
//    props.push_back(0.);
//    props.push_back(0.);
//    props.push_back(_Q.e());

//    props.push_back(_v[0]);
//    props.push_back(_v[1]);
//#if ndim_ == 3
//    props.push_back(_v[2]);
//#endif
//    block.add_row(props);
//}

/**
 * @brief Set the maximal Mach number for this particle
 *
 * @param mach Maximal Mach number of all shocks encountered in faces of the
 * cell associated with this particle
 */
void GasParticle::set_max_mach(double mach) {
    _max_mach = mach;
}

/**
 * @brief Get the maximal Mach number for this particle
 *
 * @return Maximal Mach number of all shocks encountered in faces of the cell
 * associated with this particle
 */
double GasParticle::get_max_mach() {
    return _max_mach;
}

/**
 * @brief MPI constructor
 *
 * @param buffer MPI buffer to read from
 * @param bufsize Buffer size
 * @param position Current position in the buffer (is updated)
 */
GasParticle::GasParticle(void* buffer, int bufsize, int* position)
        : Particle(buffer, bufsize, position) {
    MyMPI_Unpack(buffer, bufsize, position, &_a_grav_old[0], ndim_, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_csnd, 1, MPI_DOUBLE);
    _old_Q = StateVector(buffer, bufsize, position);
    MyMPI_Unpack(buffer, bufsize, position, &_h, 1, MPI_DOUBLE);
    _W = StateVector(buffer, bufsize, position);
    _Q = StateVector(buffer, bufsize, position);
    _dQ = StateVector(buffer, bufsize, position);
    MyMPI_Unpack(buffer, bufsize, position, &_delta_E_grav[0], ndim_,
                 MPI_DOUBLE);
    for(unsigned int i = 0; i < ndim_; i++) {
        _gradientvecs[i] = StateVector(buffer, bufsize, position);
    }
    MyMPI_Unpack(buffer, bufsize, position, &_max_radius, 1, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_centroid[0], ndim_, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_total_area, 1, MPI_DOUBLE);
    MyMPI_Unpack(buffer, bufsize, position, &_real_dt, 1, MPI_DOUBLE);
    _vorgen = NULL;
    _vorgenid = 0;
    _copies.push_back(0);
    _is_pseudo = false;
    _exports = new bool[MPIGlobal::size];
    for(int i = 0; i < MPIGlobal::size; i++) {
        _exports[i] = false;
    }
    _max_mach = 0.;
}

/**
 * @brief Dump particle data to the given MPI buffer for communication
 *
 * @param buffer MPI buffer to write to
 * @param bufsize Buffer size
 * @param position Current position in the buffer (is updated)
 */
void GasParticle::pack_data(void* buffer, int bufsize, int* position) {
    Particle::pack_data(buffer, bufsize, position);
    MyMPI_Pack(&_a_grav_old[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_csnd, 1, MPI_DOUBLE, buffer, bufsize, position);
    _old_Q.pack_data(buffer, bufsize, position);
    MyMPI_Pack(&_h, 1, MPI_DOUBLE, buffer, bufsize, position);
    _W.pack_data(buffer, bufsize, position);
    _Q.pack_data(buffer, bufsize, position);
    _dQ.pack_data(buffer, bufsize, position);
    MyMPI_Pack(&_delta_E_grav[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
    for(unsigned int i = 0; i < ndim_; i++) {
        _gradientvecs[i].pack_data(buffer, bufsize, position);
    }
    MyMPI_Pack(&_max_radius, 1, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_centroid[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_total_area, 1, MPI_DOUBLE, buffer, bufsize, position);
    MyMPI_Pack(&_real_dt, 1, MPI_DOUBLE, buffer, bufsize, position);
}

/**
 * @brief Dump particle data to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void GasParticle::dump(RestartFile& rfile) {
    Particle::dump(rfile);

    rfile.write(_a_grav_old);
    rfile.write(_old_Q);
    rfile.write(_W);
    rfile.write(_Q);
    rfile.write(_dQ);
    rfile.write(_delta_E_grav);
    rfile.write(_mesh_v);
    rfile.write(_centroid);
    rfile.write(_old_v);

    rfile.write(_csnd);
    rfile.write(_vorgenid);
    rfile.write(_h);
    rfile.write(_copies);
    rfile.write(_exports, MPIGlobal::size);
    for(unsigned int i = 0; i < ndim_; i++) {
        rfile.write(_gradientvecs[i]);
    }
    rfile.write(_max_radius);
    rfile.write(_is_pseudo);
    rfile.write(_mesh_m);
    rfile.write(_max_mach);
    rfile.write(_local_id);
    rfile.write(_eta);
    rfile.write(_total_area);

    rfile.write(_real_dt);
}

/**
 * @brief Restart constructor. Initialize the particle based on the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 */
GasParticle::GasParticle(RestartFile& rfile) : Particle(rfile) {
    rfile.read(_a_grav_old);
    rfile.read(_old_Q);
    rfile.read(_W);
    rfile.read(_Q);
    rfile.read(_dQ);
    rfile.read(_delta_E_grav);
    rfile.read(_mesh_v);
    rfile.read(_centroid);
    rfile.read(_old_v);

    rfile.read(_csnd);
    rfile.read(_vorgenid);
    rfile.read(_h);
    rfile.read(_copies);
    _exports = new bool[MPIGlobal::size];
    rfile.read(_exports, MPIGlobal::size);
    for(unsigned int i = 0; i < ndim_; i++) {
        rfile.read(_gradientvecs[i]);
    }
    rfile.read(_max_radius);
    rfile.read(_is_pseudo);
    rfile.read(_mesh_m);
    rfile.read(_max_mach);
    rfile.read(_local_id);
    rfile.read(_eta);
    rfile.read(_total_area);

    rfile.read(_real_dt);
}

/**
 * @brief Dump particle information in ASCII format to the given stream
 *
 * @param stream std::ostream to write to
 */
void GasParticle::dump_ascii(ostream& stream) {
    stream << "key\n";
    stream << get_key() << "\n";

    stream << "x\n";
#if ndim_ == 3
    stream << _x[0] << "\t" << _x[1] << "\t" << _x[2] << "\n";
#else
    stream << _x[0] << "\t" << _x[1] << "\n";
#endif

    stream << "v\n";
#if ndim_ == 3
    stream << _v[0] << "\t" << _v[1] << "\t" << _v[2] << "\n";
#else
    stream << _v[0] << "\t" << _v[1] << "\n";
#endif

    stream << "a_grav_new\n";
#if ndim_ == 3
    stream << _a_grav_new[0] << "\t" << _a_grav_new[1] << "\t" << _a_grav_new[2]
           << "\n";
#else
    stream << _a_grav_new[0] << "\t" << _a_grav_new[1] << "\n";
#endif

    stream << "old_a\n";
    stream << _old_a << "\n";

    stream << "id\n";
    stream << _id << "\n";

    stream << "starttime\n";
    stream << _starttime << "\n";

    stream << "endtime\n";
    stream << _endtime << "\n";

    stream << "comp_cost\n";
    stream << _comp_cost << "\n";

    stream << "epot\n";
    stream << _epot << "\n";

    stream << "hsoft\n";
    stream << _hsoft << "\n";

    stream << "a_grav_old\n";
#if ndim_ == 3
    stream << _a_grav_old[0] << "\t" << _a_grav_old[1] << "\t" << _a_grav_old[2]
           << "\n";
#else
    stream << _a_grav_old[0] << "\t" << _a_grav_old[1] << "\n";
#endif

    stream << "old_Q\n";
#if ndim_ == 3
    stream << _old_Q[0] << "\t" << _old_Q[1] << "\t" << _old_Q[2] << "\t"
           << _old_Q[3] << "\t" << _old_Q[4] << "\n";
#else
    stream << _old_Q[0] << "\t" << _old_Q[1] << "\t" << _old_Q[2] << "\t"
           << _old_Q[3] << "\n";
#endif

    stream << "W\n";
#if ndim_ == 3
    stream << _W[0] << "\t" << _W[1] << "\t" << _W[2] << "\t" << _W[3] << "\t"
           << _W[4] << "\n";
#else
    stream << _W[0] << "\t" << _W[1] << "\t" << _W[2] << "\t" << _W[3] << "\n";
#endif

    stream << "Q\n";
#if ndim_ == 3
    stream << _Q[0] << "\t" << _Q[1] << "\t" << _Q[2] << "\t" << _Q[3] << "\t"
           << _Q[4] << "\n";
#else
    stream << _Q[0] << "\t" << _Q[1] << "\t" << _Q[2] << "\t" << _Q[3] << "\n";
#endif

    stream << "dQ\n";
#if ndim_ == 3
    stream << _dQ[0] << "\t" << _dQ[1] << "\t" << _dQ[2] << "\t" << _dQ[3]
           << "\t" << _dQ[4] << "\n";
#else
    stream << _dQ[0] << "\t" << _dQ[1] << "\t" << _dQ[2] << "\t" << _dQ[3]
           << "\n";
#endif

    stream << "delta_E_grav\n";
#if ndim_ == 3
    stream << _delta_E_grav[0] << "\t" << _delta_E_grav[1] << "\t"
           << _delta_E_grav[2] << "\n";
#else
    stream << _delta_E_grav[0] << "\t" << _delta_E_grav[1] << "\n";
#endif

    stream << "mesh_v\n";
#if ndim_ == 3
    stream << _mesh_v[0] << "\t" << _mesh_v[1] << "\t" << _mesh_v[2] << "\n";
#else
    stream << _mesh_v[0] << "\t" << _mesh_v[1] << "\n";
#endif

    stream << "centroid\n";
#if ndim_ == 3
    stream << _centroid[0] << "\t" << _centroid[1] << "\t" << _centroid[2]
           << "\n";
#else
    stream << _centroid[0] << "\t" << _centroid[1] << "\n";
#endif

    stream << "old_v\n";
#if ndim_ == 3
    stream << _old_v[0] << "\t" << _old_v[1] << "\t" << _old_v[2] << "\n";
#else
    stream << _old_v[0] << "\t" << _old_v[1] << "\n";
#endif

    stream << "csnd\n";
    stream << _csnd << "\n";

    stream << "vorgenid\n";
    stream << _vorgenid << "\n";

    stream << "h\n";
    stream << _h << "\n";

    stream << "copies (" << _copies.size() << ")\n";
    for(unsigned int i = 0; i < _copies.size(); i++) {
        stream << _copies[i] << "\n";
    }

    stream << "exports\n";
    stream << _exports << "\n";

    stream << "gradientvecs\n";
    for(unsigned int i = 0; i < ndim_; i++) {
#if ndim_ == 3
        stream << _gradientvecs[i][0] << "\t" << _gradientvecs[i][1] << "\t"
               << _gradientvecs[i][2] << "\t" << _gradientvecs[i][3] << "\t"
               << _gradientvecs[i][4] << "\n";
#else
        stream << _gradientvecs[i][0] << "\t" << _gradientvecs[i][1] << "\t"
               << _gradientvecs[i][2] << "\t" << _gradientvecs[i][3] << "\n";
#endif
    }

    stream << "max_radius\n";
    stream << _max_radius << "\n";

    stream << "is_pseudo\n";
    stream << _is_pseudo << "\n";

    stream << "mesh_m\n";
    stream << _mesh_m << "\n";

    stream << "max_mach\n";
    stream << _max_mach << "\n";

    stream << "local_id\n";
    stream << _local_id << "\n";

    stream << "eta\n";
    stream << _eta << "\n";

    stream << "real_dt\n";
    stream << _real_dt << "\n";

    stream << "total_area\n";
    stream << _total_area << "\n";
}

/**
 * @brief Set the geometrical centroid of the cell associated with the particle
 *
 * @param centroid Coordinates of the centroid
 */
void GasParticle::set_centroid(Vec centroid) {
    _centroid = centroid;
}

/**
 * @brief Get the geometrical centroid of the last cell associated with this
 * particle
 *
 * @return Reference to the coordinates of the centroid
 */
Vec& GasParticle::get_centroid() {
    return _centroid;
}

/**
 * @brief Set the index of the particle in the local ParticleVector
 *
 * @param local_id Index of the particle in the local ParticleVector
 */
void GasParticle::set_local_id(unsigned int local_id) {
    _local_id = local_id;
}

/**
 * @brief Get the index of the particle in the local ParticleVector
 *
 * @return Index of the particle in the local ParticleVector
 */
unsigned int GasParticle::get_local_id() {
    return _local_id;
}

/**
 * @brief Set the gravitational factor for this particle
 *
 * @param eta Gravitational factor
 */
void GasParticle::set_eta(double eta) {
    _eta = eta;
}

/**
 * @brief Get the gravitational factor for this particle
 *
 * @return Gravitational factor
 */
double GasParticle::get_eta() {
    return _eta;
}

/**
 * @brief Set the real timestep of the particle (as opposed to the integer
 * timestep stored indirectly in Particle)
 *
 * @param real_dt Real timestep of the particle
 */
void GasParticle::set_real_timestep(double real_dt) {
    _real_dt = real_dt;
}

/**
 * @brief Get the real timestep of the particle (as opposed to the integer
 * timestep stored indirectly in Particle)
 *
 * @return Reall timestep of the particle
 */
double GasParticle::get_real_timestep() {
    return _real_dt;
}

/**
 * @brief Constructor
 *
 * @param p GasParticle which throws the exception
 */
UpdateQException::UpdateQException(GasParticle& p) : _p(p) {}

/**
 * @brief Human readable information about the exception
 *
 * @return C-string with more information
 */
const char* UpdateQException::what() const throw() {
    cerr << "Error! Negative mass/energy after update of conserved quantities!"
         << endl;
    cerr << "Particle: " << _p.id() << endl;
    cerr << _p.x() << "\t" << _p.y();
#if ndim_ == 3
    cerr << "\t" << _p.z();
#endif
    cerr << endl;
    StateVector Q = _p.get_Qvec();
    StateVector dQ = _p.get_dQvec();
    cerr << "Q: " << Q.rho() << "\t" << Q.px() << "\t" << Q.py() << "\t";
#if ndim_ == 3
    cerr << Q.pz() << "\t";
#endif
    cerr << Q.e() << endl;
    cerr << "dQ: " << dQ.rho() << "\t" << dQ.px() << "\t" << dQ.py() << "\t";
#if ndim_ == 3
    cerr << dQ.pz() << "\t";
#endif
    cerr << dQ.e() << endl;
    cerr << "Gradients:" << endl;
    for(unsigned int i = 0; i < ndim_; i++) {
        cerr << _p.get_gradient(0, i);
        for(unsigned int j = 1; j < ndim_ + 2; j++) {
            cerr << "\t" << _p.get_gradient(j, i);
        }
        cerr << endl;
    }
    return "UpdateQException";
}
