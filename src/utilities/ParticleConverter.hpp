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
 * @file ParticleConverter.hpp
 *
 * @brief Interface for objects that can convert a particle from one type into a
 * particle of a different type
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PARTICLECONVERTER_HPP
#define PARTICLECONVERTER_HPP

#include "ParticleTypes.hpp"

class Particle;

/**
 * @brief Interface for objects that can convert a particle from one type into a
 * particle of a different type
 */
class ParticleConverter {
  public:
    /**
     * @brief Check if the given ParticleType can be converted
     *
     * @param type ParticleType
     * @return True if the converter applies to this ParticleType
     */
    virtual bool do_conversion(ParticleType type) = 0;

    /**
     * @brief Convert the given Particle to a Particle with a different type
     *
     * The function can also return the original pointer, in which case no
     * conversion takes place. If the Particle is converted, the original
     * Particle is deleted.
     *
     * @param particle Pointer to a Particle that should potentially be
     * converted
     * @return Pointer to a possibly new converted Particle
     */
    virtual Particle* convert(Particle* particle) = 0;
};

#endif  // PARTICLECONVERTER_HPP
