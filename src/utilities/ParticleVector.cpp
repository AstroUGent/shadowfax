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
 * @file ParticleVector.cpp
 *
 * @brief Specialize vector for particles: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ParticleVector.hpp"
#include "ParallelSorter.hpp"
#include "Particle.hpp"
#include "GasParticle.hpp"
#include "DMParticle.hpp"
#include "RestartFile.hpp"
#include <fstream>
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
using namespace std;

/**
 * @brief Constructor
 *
 * @param global_timestep Flag indicating if a global or individual timesteps
 * should be used
 * @param container RectangularBox specifying the dimensions of the box
 * containing all particles
 * @param periodic Flag indicating whether the box is periodic (true) or
 * reflective (false)
 * @param do_ewald Flag indicating if Ewald tables for periodic force
 * corrections should be used
 * @param alpha Ewald \f$\alpha\f$ factor
 * @param size Size of the precomputed Ewald tables
 */
ParticleVector::ParticleVector(bool global_timestep, RectangularBox container,
                               bool periodic, bool do_ewald, double alpha,
                               unsigned int size)
    : _tree(container.get_cuboid(), periodic, do_ewald, alpha, size),
      _container(container){
    double box[2*ndim_] = {0.};
    _container.get_bounding_box(box);
    _header.set_box(box);
    _header.set_global_timestep(global_timestep);
    _header.set_periodic(periodic);
    for(unsigned int i = 0; i < PARTTYPE_COUNTER; i++){
        _offsets[i] = 0;
        _sizes[i] = 0;
    }
    _sizes[PARTTYPE_COUNTER] = 0;
}

/**
 * @brief Destructor
 *
 * Free particle memory and print out the sorting and tree building timers to
 * the stdout.
 */
ParticleVector::~ParticleVector(){
    for(unsigned int i = 0; i < _sizes[PARTTYPE_COUNTER]; i++){
        delete _particles[i];
    }
    cout << "Spent " << _sorttimer.value() << "s sorting particles" << endl;
    cout << "Spent " << _treetimer.value() << "s building tree" << endl;
}

/**
 * @brief Set the dimensions of the box containing all particles
 *
 * @param container RectangularBox specifying the new dimensions of the box that
 * contains all particles
 */
void ParticleVector::set_container(RectangularBox container){
    _container = container;
    double box[ndim_+ndim_] = {0.};
    _container.get_bounding_box(box);
    _header.set_box(box);
    _tree.reset(_container.get_cuboid());
}

/**
 * @brief Indicate if a periodic or reflective box is used
 *
 * @param periodic True if the box is periodic, false if it is reflective
 */
void ParticleVector::set_periodic(bool periodic){
    _header.set_periodic(periodic);
    _tree.set_periodic(periodic);
}

/**
 * @brief Signal that we are finished with adding particles to the vector
 *
 * This will consolidate the total number of particles in the Header.
 */
void ParticleVector::finalize(){
    unsigned int globsizes[PARTTYPE_COUNTER];
    MyMPI_Allreduce(_sizes, globsizes, PARTTYPE_COUNTER, MPI_UNSIGNED, MPI_SUM);
    _header.set_ngaspart(globsizes[PARTTYPE_GAS]);
    _header.set_ndmpart(globsizes[PARTTYPE_DM]);
}

/**
 * @brief Add the given GasParticle to the gas particle list
 *
 * @param particle GasParticle to add
 */
void ParticleVector::add_gas_particle(GasParticle *particle){
    particle->set_key(_container.get_key(particle->get_position()));
    _particles.insert(_particles.begin()+_offsets[PARTTYPE_GAS]+
                      _sizes[PARTTYPE_GAS], particle);
    // shift the offsets of all particles of other types
    for(int i = PARTTYPE_GAS+1; i < PARTTYPE_COUNTER; i++){
        _offsets[i]++;
    }
    _sizes[PARTTYPE_GAS]++;
    _sizes[PARTTYPE_COUNTER]++;
}

/**
 * @brief Add the given DMParticle to the dark matter particle list
 *
 * @param particle DMParticle to add
 */
void ParticleVector::add_DM_particle(DMParticle *particle){
    particle->set_key(_container.get_key(particle->get_position()));
    _particles.insert(_particles.begin()+_offsets[PARTTYPE_DM]+
                      _sizes[PARTTYPE_DM], particle);
    // shift the offsets of all particles of other types
    for(int i = PARTTYPE_DM+1; i < PARTTYPE_COUNTER; i++){
        _offsets[i]++;
    }
    _sizes[PARTTYPE_DM]++;
    _sizes[PARTTYPE_COUNTER]++;
}

/**
 * @brief Construct the octtree containing all particles by adding all particles
 * to it
 */
void ParticleVector::construct_tree(){
    for(unsigned int i = 0; i < _particles.size(); i++){
        _tree.add_particle(_particles[i]);
    }
}

/**
 * @brief Sort the particles in space-filling Hilbert order and construct an
 * octtree containing all particles
 *
 * We first calculate Hilbert keys for all particles. Then we sort the gas
 * particle and dark matter particle lists seperately. The former is just sorted
 * on key, the latter is sorted to provide approximate equal computational cost
 * on all MPI processes.
 *
 * All particles are then added to a local octtree. We then check which parts of
 * the tree are on other MPI processes and complete the local tree with
 * pseudonodes from other processes.
 *
 * When this method exits, we have a complete octtree on all processes and all
 * particles will be redistributed over the processes so that every process
 * holds a compact block of particles with an approximately equal computational
 * cost.
 */
void ParticleVector::sort(){
    _sorttimer.start();
    for(unsigned int i = 0; i < _sizes[PARTTYPE_COUNTER]; i++){
        _particles[i]->set_key(
                    _container.get_key(_particles[i]->get_position())
                    );
    }
    _tree.reset(_container.get_cuboid());
    ParallelSorter sorter;
    // weightsort for DM?
    sorter.sort(_particles);

    // check for duplicate keys
    for(unsigned int i = 1; i < _particles.size(); i++){
        if(_particles[i-1]->get_key() == _particles[i]->get_key()){
            // randomize the key to make sure tree building will succeed
            // we do have to make sure that the particle falls withing the same
            // lowest level node
#if ndim_==3
            unsigned int nodekey = _particles[i]->get_key() & 7;
            if(nodekey == 7){
                _particles[i-1]->set_key(_particles[i-1]->get_key()-1);
            } else {
                _particles[i]->set_key(_particles[i]->get_key()+1);
            }
#else
            unsigned int nodekey = _particles[i]->get_key() & 3;
            if(nodekey == 3){
                _particles[i-1]->set_key(_particles[i-1]->get_key()-1);
            } else {
                _particles[i]->set_key(_particles[i]->get_key()+1);
            }
#endif
        }
    }

    // make sure the particles are ordered per type
    // particles of the same type should keep their relative ordering if we use
    // stable_sort
    std::stable_sort(_particles.begin(), _particles.end(), compare_type);

    // set the offsets and sizes of the local particles
    _sizes[PARTTYPE_COUNTER] = _particles.size();
    for(int type = 0; type < PARTTYPE_COUNTER; type++){
        _sizes[type] = 0;
    }
    for(unsigned int i = 0; i < _sizes[PARTTYPE_COUNTER]; i++){
        _sizes[_particles[i]->type()]++;
    }
    for(int type = 1; type < PARTTYPE_COUNTER; type++){
        _offsets[type] = _offsets[type-1]+_sizes[type-1];
    }

    // set the internal ids for the gas particles
    for(unsigned int i = _offsets[PARTTYPE_GAS]; i < _sizes[PARTTYPE_GAS]; i++){
        ((GasParticle*)_particles[i])->set_local_id(i-_offsets[PARTTYPE_GAS]);
    }

    _sorttimer.stop();
    _treetimer.start();
    construct_tree();
    if(MPIGlobal::size < 2){
        _tree.finalize();
        _treetimer.stop();
        return;
    }

    unsigned long maxkey = 1;
    maxkey <<= 61;
    unsigned long loc_min = _particles[0]->get_key();
    unsigned long loc_max = _particles.back()->get_key();
    vector<unsigned long> glo_min(MPIGlobal::size);
    vector<unsigned long> glo_max(MPIGlobal::size);
    MyMPI_Allgather(&loc_min, 1, MPI_UNSIGNED_LONG, &glo_min[0], 1,
            MPI_UNSIGNED_LONG);
    MyMPI_Allgather(&loc_max, 1, MPI_UNSIGNED_LONG, &glo_max[0], 1,
            MPI_UNSIGNED_LONG);

    glo_min[0] = 0;
    glo_max[MPIGlobal::size-1] = maxkey-1;
    for(unsigned int i = MPIGlobal::size-1; i--;){
        unsigned int count = 0;
        while(((glo_min[i+1]>>1)<<(count+1)) > glo_max[i]){
            glo_min[i+1] >>= 1;
            count++;
        }
        glo_min[i+1] <<= count;
        glo_max[i] = glo_min[i+1]-1;
    }
    // add pseudoparticles to tree
    for(int i = MPIGlobal::size; i--;){
        if(i != MPIGlobal::rank){
            _tree.add_pseudoparticles(glo_min[i], glo_max[i], i);
        }
    }
    // this should not happen before the tree is complete (otherwise it could
    // happen that you still have to split leaves, which can not happen anymore
    // if pseudoparticles have been added)
    _tree.set_keyrange(glo_min[MPIGlobal::rank], glo_max[MPIGlobal::rank]);
    _tree.finalize();
    _tree.exchange_pseudonodes();
    _treetimer.stop();
}

/**
 * @brief Print the local particles to a file with the given name
 *
 * @param filename Name of the file to write
 */
void ParticleVector::print_local_particles(string filename){
    // construct file name
    stringstream name;
    name << filename;
    if(MPIGlobal::size > 1){
        name << "." << MPIGlobal::rank;
    }
    name << ".txt";

    // open file
    ofstream ofile(name.str().c_str());

    // write header
    ofile << "#x\ty";
#if ndim_==3
    ofile << "\tz";
#endif
    ofile << "\n";

    // gas particles
    if(_sizes[PARTTYPE_GAS]){
        ofile << "#gas\n";

        for(unsigned int i = _offsets[PARTTYPE_GAS]; i < _sizes[PARTTYPE_GAS];
            i++){
            Particle* p = _particles[i];
            ofile << p->x() << "\t" << p->y();
#if ndim_==3
            ofile << "\t" << p->z();
#endif
            ofile << "\n";
        }
    }

    // dm particles
    if(_sizes[PARTTYPE_DM]){
        ofile << "#dm\n";

        for(unsigned int i = _offsets[PARTTYPE_DM]; i < _sizes[PARTTYPE_DM];
            i++){
            Particle* p = _particles[i];
            ofile << p->x() << "\t" << p->y();
#if ndim_==3
            ofile << "\t" << p->z();
#endif
            ofile << "\n";
        }
    }
}

/**
 * @brief Get the Header
 *
 * @return Reference to the Header
 */
Header& ParticleVector::get_header(){
    return _header;
}

/**
 * @brief Set the number of active particles at the current system time
 *
 * @param numactive Number of active particles
 */
void ParticleVector::set_numactive(unsigned int numactive){
    _numactive = numactive;
}

/**
 * @brief Get the number of active particles at the current system time
 *
 * @return Number of active particles
 */
unsigned int ParticleVector::get_numactive(){
    return _numactive;
}

/**
 * @brief Dump the particle vector to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void ParticleVector::dump(RestartFile &rfile){
    _header.dump(rfile);
    unsigned int vsize = _sizes[PARTTYPE_GAS];
    rfile.write(vsize);
    for(unsigned int i = _offsets[PARTTYPE_GAS];
        i < _offsets[PARTTYPE_GAS] + _sizes[PARTTYPE_GAS]; i++){
        _particles[i]->dump(rfile);
    }
    vsize = _sizes[PARTTYPE_DM];
    rfile.write(vsize);
    for(unsigned int i = _offsets[PARTTYPE_DM];
        i < _offsets[PARTTYPE_DM] + _sizes[PARTTYPE_DM]; i++){
        _particles[i]->dump(rfile);
    }

}

/**
 * @brief Restart constructor. Initialize the particle vector from the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 * @param box RectangularBox specifying the dimensions of the box containing all
 * particles
 * @param periodic Flag indicating if the box is periodic (true) or reflective
 * (false)
 * @param do_ewald Flag indicating if precomputed Ewald tables for periodic
 * force corrections should be loaded
 * @param alpha Ewald \f$\alpha\f$ factor
 * @param size Size of the Ewald table
 */
ParticleVector::ParticleVector(RestartFile &rfile, RectangularBox& box,
                               bool periodic, bool do_ewald, double alpha,
                               unsigned int size)
    : _tree(box.get_cuboid(), periodic, do_ewald, alpha, size), _container(box),
      _header(rfile){
    unsigned int vsize;
    rfile.read(vsize);
    _particles.resize(vsize, NULL);
    _sizes[PARTTYPE_GAS] = vsize;
    _sizes[PARTTYPE_COUNTER] = vsize;
    _offsets[PARTTYPE_GAS] = 0;
    for(unsigned int i = 0; i < vsize; i++){
        _particles[i] = new GasParticle(rfile);
    }
    rfile.read(vsize);
    _sizes[PARTTYPE_DM] = vsize;
    _sizes[PARTTYPE_COUNTER] += vsize;
    _offsets[PARTTYPE_DM] = _sizes[PARTTYPE_GAS];
    _particles.resize(_sizes[PARTTYPE_COUNTER], NULL);
    for(unsigned int i = _offsets[PARTTYPE_DM];
        i < _offsets[PARTTYPE_DM]+_sizes[PARTTYPE_DM]; i++){
        _particles[i] = new DMParticle(rfile);
    }
}
