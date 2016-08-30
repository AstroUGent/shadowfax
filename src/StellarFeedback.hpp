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

class DMParticle;
class GasParticle;
class ParticleVector;
class RestartFile;
class StarParticle;
class StellarFeedbackData;

/**
 * @brief General interface for stellar feedback mechanisms
 */
class StellarFeedback {
  public:
    virtual ~StellarFeedback() {}

    /**
     * @brief Check if the given StarParticle does feedback during the next
     * time step
     *
     * @param star StarParticle
     * @param starttime Physical start time of the next time step
     * @param endtime Physical end time of the next time step
     * @return True if the StarParticle does feedback
     */
    virtual bool does_feedback(StarParticle* star, double starttime,
                               double endtime) = 0;

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

/**
 * @brief TreeFilter used to filter out star particles that give feedback
 */
class StellarFeedbackTreeFilter {
  private:
    /*! @brief StellarFeedback implementation used */
    StellarFeedback* _feedback;

    /*! @brief Physical start time of the current time step */
    double _starttime;

    /*! @brief Physical end time of the current time step */
    double _endtime;

  public:
    /**
     * @brief Constructor
     *
     * @param feedback StellarFeedback implementation used
     * @param starttime Physical start time of the current time step
     * @param endtime Physical end time of the current time step
     */
    inline StellarFeedbackTreeFilter(StellarFeedback* feedback,
                                     double starttime, double endtime) {
        _feedback = feedback;
        _starttime = starttime;
        _endtime = endtime;
    }

    /**
     * @brief Do feedback for the given GasParticle
     *
     * @param gas GasParticle
     * @return False, since a GasParticle does not give stellar feedback
     */
    inline bool do_gas(GasParticle* gas) {
        return false;
    }

    /**
     * @brief Do feedback for the given DMParticle
     *
     * @param dm DMParticle
     * @return False, since a DMParticle does not give stellar feedback
     */
    inline bool do_dm(DMParticle* dm) {
        return false;
    }

    /**
     * @brief Do feedback for the given StarParticle
     *
     * @param star StarParticle
     * @return True if the StellarFeedback implementation signals that this
     * StarParticle does feedback
     */
    inline bool do_star(StarParticle* star) {
        return _feedback->does_feedback(star, _starttime, _endtime);
    }
};

#endif  // STELLARFEEDBACK_HPP
