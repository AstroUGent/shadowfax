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
 * @file ParticleTypes.hpp
 *
 * @brief Types of particles in the simulation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PARTICLETYPES_HPP
#define PARTICLETYPES_HPP

/**
 * @brief Particle types
 */
enum ParticleType{
    /*! Gas particle */
    PARTTYPE_GAS = 0,
    /*! Dark matter particle */
    PARTTYPE_DM,
    /*! Counter of the number of types (make sure this stays the last element in
     *  the enum!) */
    PARTTYPE_COUNTER
};

#endif // PARTICLETYPES_HPP
