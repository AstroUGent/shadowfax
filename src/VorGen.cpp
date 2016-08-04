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
 * @file VorGen.cpp
 *
 * @brief Voronoi cell generator: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "VorGen.hpp"
#include "Error.hpp"      // for my_exit
#include "MPIGlobal.hpp"  // for rank
#include <cstdlib>        // for NULL
#include <iostream>       // for cout
#include <iterator>
#include <ostream>  // for operator<<, ostream, etc
using namespace std;

#if ndim_ == 3
/**
 * @brief Constructor
 *
 * Initialize a point with the given coordinates.
 *
 * @param x x-coordinate of the point
 * @param y y-coordinate of the point
 * @param z z-coordinate of the point
 */
VorGen::VorGen(double x, double y, double z) {
    _cell = NULL;
    _mirrored = 0;
    _mirror = NULL;
    _particle = NULL;
    _p.set(x, y, z);
    _original = 0;
    _process = MPIGlobal::rank;
    _flag = false;
}
#else
/**
 * @brief Constructor
 *
 * Initialize a point with the given coordinates.
 *
 * @param x x-coordinate of the point
 * @param y y-coordinate of the point
 */
VorGen::VorGen(double x, double y) {
    _mirror = NULL;
    _particle = NULL;
    _cell = NULL;
    _mirrored = 0;
    _p.set(x, y);
    _original = 0;
    _process = MPIGlobal::rank;
    _flag = false;
}
#endif

/**
 * @brief Destructor
 *
 * Does nothing really...
 */
VorGen::~VorGen() {
    _tetrahedra.clear();
}

// deprecated
/**
 * @brief Move the point to the given coordinates
 *
 * @deprecated No longer used
 *
 * @param x New x-coordinate of the point
 * @param y New y-coordinate of the point
 * @param z New z-coordinate of the point
 */
void VorGen::move(double x, double y, double z) {
    //    _p.set(x,y,z);
}

/**
 * @brief Add the given simplex to the list
 *
 * @param tetrahedron Simplex to add
 */
void VorGen::add_tetrahedron(Simplex* tetrahedron) {
    _tetrahedra.push_back(tetrahedron);
}

/**
 * @brief Add the given simplex index to the list
 *
 * @param index Index of a Simplex in the DelTess simplex list
 */
void VorGen::add_simplex(unsigned int index) {
    _simplices.push_back(index);
}

/**
 * @brief Remove a tetrahedron from the list
 *
 * If the given tetrahedron is not in the list, this function can give rise to
 * undefined behaviour.
 *
 * @param tetrahedron Simplex to remove
 */
void VorGen::remove_tetrahedron(Simplex* tetrahedron) {
    list<Simplex*>::reverse_iterator it = _tetrahedra.rbegin();
    while(*it != tetrahedron) {
        it++;
    }
    _tetrahedra.erase(--(it.base()));
}

/**
 * @brief Remove the given simplex index from the list
 *
 * If the given index is not in the list, the code will abort.
 *
 * @param index Index of a Simplex in the DelTess simplex list
 */
void VorGen::remove_simplex(unsigned int index) {
    list<unsigned int>::iterator it = _simplices.begin();
    while(*it != index && it != _simplices.end()) {
        it++;
    }
    if(it == _simplices.end()) {
        cout << "Trying to remove simplex that does not exist. Error!" << endl;
        my_exit();
    }
    _simplices.erase(it);
}

/**
 * @brief Access the tetrahedra
 *
 * @return Simplex std::list
 */
list<Simplex*> VorGen::get_tetrahedra() {
    return _tetrahedra;
}

/**
 * @brief Get the simplex indices
 *
 * @return std::list with simplex indices
 */
list<unsigned int> VorGen::get_simplices() {
    return _simplices;
}

/**
 * @brief Print the Point (in ascii) to the given stream
 *
 * According to the number of dimensions, the format is
 * p:\\t(c1,c2)\n
 * or
 * p:\\t(c1,c2,c3)\n
 * (with the ci doubles)
 *
 * @param stream std::ostream to write to
 */
void VorGen::print(ostream& stream) {
    stream << _p[0] << "\t" << _p[1];
#if ndim_ == 3
    stream << "\t" << _p[2];
#endif
}

/**
 * @brief Set the VorCell associated with this point
 *
 * @param cell VorCell associated with this generator
 */
void VorGen::set_cell(VorCell* cell) {
    _cell = cell;
}

/**
 * @brief Access the VorCell of this point
 *
 * @return VorCell associated with this generator
 */
VorCell* VorGen::get_cell() {
    return _cell;
}

/**
 * @brief Set the GasParticle associated with this point
 *
 * @param particle GasParticle associated with this generator
 */
void VorGen::set_particle(GasParticle* particle) {
    _particle = particle;
}

/**
 * @brief Get the GasParticle associated with this point
 *
 * @return GasParticle associated with this generator
 */
GasParticle* VorGen::get_particle() {
    return _particle;
}

/**
 * @brief Set the index of the particle associated with this generator in the
 * local ParticleVector
 *
 * @param index Index of a particle in the ParticleVector
 */
void VorGen::set_particle_id(unsigned int index) {
    _particle_index = index;
}

/**
 * @brief Get the index of the particle associated with this generator in the
 * local ParticleVector
 *
 * @return Index of a particle in the ParticleVector
 */
unsigned int VorGen::get_particle_id() {
    return _particle_index;
}

/**
 * @brief Deprecated
 *
 * @param mirror Deprecated
 */
void VorGen::set_mirror(VorGen* mirror) {
    _mirror = mirror;
}

/**
 * @brief Deprecated
 *
 * @return Deprecated
 */
VorGen* VorGen::get_mirror() {
    return _mirror;
}

/**
 * @brief Set the search radius for the DelTess completion
 *
 * @param sr Search radius for the DelTess completion
 */
void VorGen::set_search_radius(double sr) {
    _sr = sr;
}

/**
 * @brief Get the search radius for the DelTess completion
 *
 * @return Search radius for the DelTess completion
 */
double VorGen::get_search_radius() {
    return _sr;
}

/**
 * @brief Signal a mirror copy of this point
 *
 * @param casus Index of the boundary at which this point was mirrored
 */
void VorGen::mirror(unsigned int casus) {
    _mirrored += 1 << casus;
}

/**
 * @brief Check if this point was mirror copied for the given boundary
 *
 * @param casus Index of a boundary
 * @return True if the point has a mirror copy for the given boundary, false
 * otherwise
 */
bool VorGen::mirrored(unsigned int casus) {
    return (_mirrored >> casus) & 1;
}

/**
 * @brief Set the ID of this point
 *
 * We hack the _mirrored variable, because we don't use it when we need ids and
 * it would be overkill to add a special field for this.
 *
 * It could also be that mirror in its original meaning is not used anymore...
 *
 * @param id Index of the point
 */
void VorGen::set_id(unsigned int id) {
    _mirrored = id;
}

/**
 * @brief Get the index of this point
 *
 * @return Index of the point
 */
unsigned int VorGen::get_id() {
    return _mirrored;
}

/**
 * @brief Clear the Simplex list
 */
void VorGen::reset_simplices() {
    _tetrahedra.clear();
}

/**
 * @brief Get the position of this point
 *
 * @return Position of the point
 */
Vec VorGen::get_position() {
    return _p;
}

/**
 * @brief Set the ID of the GasParticle associated to this point
 *
 * @param id ID of the GasParticle associated to this generator
 */
void VorGen::set_original(unsigned long id) {
    _original = id;
}

/**
 * @brief Get the ID of the GasParticle associated to this point
 *
 * @return ID of the GasParticle associated to this generator
 */
unsigned long VorGen::get_original() {
    return _original;
}

/**
 * @brief Set the rank of the process that holds the original copy of this point
 *
 * @param process MPI process that holds the original copy of this generator
 */
void VorGen::set_process(unsigned int process) {
    _process = process;
}

/**
 * @brief Get the rank of the process that holds the original copy of this point
 *
 * @return MPI process that holds the original copy of this generator
 */
unsigned int VorGen::get_process() {
    return _process;
}

/**
 * @brief Set the coordinates of this point, rescaled to the range [1,2] (for
 * exact geometrical tests)
 *
 * @param p12 Coordinates of the point in the range [1,2]
 */
void VorGen::set_p12(Vec p12) {
    _p12 = p12;
}

/**
 * @brief Get the coordinates of this point, rescaled to the range [1,2] (for
 * exact geometrical tests)
 *
 * @return Coordinates of the point in the range [1,2]
 */
Vec& VorGen::get_p12() {
    return _p12;
}

/**
 * @brief Flag this VorGen
 */
void VorGen::flag() {
    _flag = true;
}

/**
 * @brief Unflag this VorGen
 */
void VorGen::unflag() {
    _flag = false;
}

/**
 * @brief Check if this VorGen is flagged
 *
 * @return True if the VorGen is flagged
 */
bool VorGen::flagged() {
    return _flag;
}
