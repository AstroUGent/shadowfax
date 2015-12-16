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
 * @file ParticleFactory.hpp
 *
 * @brief Factory to facilitate sending and receiving particles of different
 * types
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PARTICLEFACTORY_HPP
#define PARTICLEFACTORY_HPP

#include "GasParticle.hpp"
#include "DMParticle.hpp"
#include "MPIMethods.hpp"

#include <iostream>

/**
 * @brief Factory to facilitate the sending and receiving of particles of
 * different types
 */
class ParticleFactory{
public:
    /**
     * @brief Load a particle that was packed using ParticleFactory::dump from
     * the given buffer and return a pointer to it
     *
     * @param buffer Buffer to read from
     * @param bufsize Size of the buffer
     * @param position Current position in the buffer (is updated)
     * @return Pointer to the newly created Particle instance
     */
    static Particle* load(void* buffer, int bufsize, int* position){
        int type;
        MyMPI_Unpack(buffer, bufsize, position, &type, 1, MPI_INT);
        if(type == PARTTYPE_GAS){
            return new GasParticle(buffer, bufsize, position);
        }
        if(type == PARTTYPE_DM){
            return new DMParticle(buffer, bufsize, position);
        }
        std::cerr << "Unknown particle type: " << type << std::endl;
        my_exit();
        return NULL;
    }

    /**
     * @brief Write the given Particle to the given buffer. The Particle can
     * unambiguously be recreated on another process by calling
     * ParticleFactory::load
     *
     * @param particle Particle instance to write
     * @param buffer Buffer to write to
     * @param bufsize Size of the buffer
     * @param position Current position in the buffer (is updated)
     */
    static void dump(Particle *particle, void* buffer, int bufsize,
                     int* position){
        int type = particle->type();
        MyMPI_Pack(&type, 1, MPI_INT, buffer, bufsize, position);
        particle->pack_data(buffer, bufsize, position);
    }
};

#endif // PARTICLEFACTORY_HPP
