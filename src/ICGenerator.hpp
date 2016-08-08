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
 * @file ICGenerator.hpp
 *
 * @brief Classes used to generate initial conditions for simulations: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ICGENERATOR_HPP
#define ICGENERATOR_HPP

class ParticleVector;

/**
  * \brief Interface for initial conditions generators
  */
class ICGenerator {
  public:
    virtual ~ICGenerator() {}

    /**
      * \brief Generate initial conditions and return them as a ParticleVector
      *
      * @param conserved_variables Generate conserved variables?
      * @return ParticleVector containing the generated particles
      */
    virtual ParticleVector generate(bool conserved_variables = false) = 0;
};

/**
  * \brief Possible modes for initial condition generation
  */
enum ICMode {
    /*! Set up a regularized uniform random grid */
    IC_RAND,
    /*! Set up a cartesian grid */
    IC_CART
};

#endif  // ICGENERATOR_HPP
