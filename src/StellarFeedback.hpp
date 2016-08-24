/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *                    Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
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
 * @file StellarFeedback.hpp
 *
 * @brief General interface for stellar feedback mechanisms
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#ifndef STELLARFEEDBACK_HPP
#define STELLARFEEDBACK_HPP

class ParticleVector;
class RestartFile;
class StarParticle;

/**
 * @brief General interface for data that needs to be stored in the StarParticle
 * for the stellar feedback
 */
class StellarFeedbackData {
  public:
    // we need to provide an empty destructor to make sure this (empty)
    // interface exists as a well defined entity
    virtual ~StellarFeedbackData() {}
};

/**
 * @brief General interface for stellar feedback mechanisms
 */
class StellarFeedback {
  public:
    virtual ~StellarFeedback() {}

    /**
     * @brief Give stellar feedback from the given StarParticle to the gas in
     * the given ParticleVector
     *
     * @param star StarParticle that does feedback
     * @param particles ParticleVector containing the gas that receives feedback
     * @param dt Time step over which the feedback is done
     */
    virtual void do_feedback(StarParticle* star, ParticleVector& particles,
                             double dt) = 0;

    /**
     * @brief Initialize the StellarFeedbackData for the given StarParticle
     *
     * @param star StarParticle for which the StellarFeedbackData is initialized
     * @return Pointer to a StellarFeedbackData instance for the StarParticle
     */
    virtual StellarFeedbackData* initialize_data(StarParticle* star) = 0;
};

#endif  // STELLARFEEDBACK_HPP
