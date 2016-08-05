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
 * @file DMParticle.hpp
 *
 * @brief Dark matter particle: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef DMPARTICLE_HPP
#define DMPARTICLE_HPP

#include "Particle.hpp"       // for Particle
#include "ParticleTypes.hpp"  // for ParticleType, etc
#include <ostream>            // for ostream

class RestartFile;

/**
  * @brief Representation of a cold dark matter particle
  *
  * A dark matter particle extends Particle with a mass and only interacts
  * through gravity
  */
class DMParticle : public Particle {
  private:
    /*! @brief Mass of the dark matter particle */
    double _mass;

  public:
    DMParticle();
    DMParticle(Vec pos);
    DMParticle(void* buffer, int bufsize, int* position);
    virtual ~DMParticle() {}

    /**
      * @brief Mark this Particle as a dark matter particle
      *
      * @return PARTTYPE_DM
      */
    ParticleType type() {
        return PARTTYPE_DM;
    }

    void set_mass(double mass);

    /**
      * @brief Get the mass of the dark matter particle
      *
      * @return The mass of the dark matter particle
      */
    inline double get_mass() {
        return _mass;
    }

    // NOTE: the "virtual" is really important here
    // it's the thing that makes sure that this particular version of pack_data
    // is called when you call it from a general Particle pointer
    // without the virtual, you call Particle::pack_data instead, which is
    // probably not what you want
    virtual void pack_data(void* buffer, int bufsize, int* position);

    virtual void dump(RestartFile& rfile);
    DMParticle(RestartFile& rfile);

    void dump_ascii(std::ostream& stream);
};

#endif  // DMPARTICLE_HPP
