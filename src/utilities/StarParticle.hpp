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
 * @file StarParticle.hpp
 *
 * @brief Star particle: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef STARPARTICLE_HPP
#define STARPARTICLE_HPP

#include "Particle.hpp"       // for Particle
#include "ParticleTypes.hpp"  // for ParticleType, etc
#include <ostream>            // for ostream

class RestartFile;
class Vec;

/**
  * @brief Representation of a star particle
  *
  * A star particle extends Particle with a mass and only interacts through
  * gravity
  */
class StarParticle : public Particle {
  private:
    /*! @brief Mass of the star particle */
    double _mass;

    /*! @brief Age of the star particle */
    double _age;

  public:
    StarParticle();
    StarParticle(Vec pos);
    StarParticle(void* buffer, int bufsize, int* position);
    virtual ~StarParticle() {}

    /**
      * @brief Mark this Particle as a star particle
      *
      * @return PARTTYPE_STAR
      */
    ParticleType type() {
        return PARTTYPE_STAR;
    }

    void set_mass(double mass);

    /**
      * @brief Get the mass of the star particle
      *
      * @return The mass of the star particle
      */
    inline double get_mass() {
        return _mass;
    }

    /**
     * @brief Set the age of the star particle
     *
     * @param age Age of the star particle
     */
    void set_age(double age) {
        _age = age;
    }

    /**
     * @brief Get the age of the star particle
     *
     * @return Age of the star particle
     */
    double get_age() {
        return _age;
    }

    // NOTE: the "virtual" is really important here
    // it's the thing that makes sure that this particular version of pack_data
    // is called when you call it from a general Particle pointer
    // without the virtual, you call Particle::pack_data instead, which is
    // probably not what you want
    virtual void pack_data(void* buffer, int bufsize, int* position);

    virtual void dump(RestartFile& rfile);
    StarParticle(RestartFile& rfile);

    void dump_ascii(std::ostream& stream);
};

#endif  // STARPARTICLE_HPP
