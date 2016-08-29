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
 * @file StellarFeedbackData.hpp
 *
 * @brief General interface for extra data stored in a StarParticle for stellar
 * feedback purposes
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef STELLARFEEDBACKDATA_HPP
#define STELLARFEEDBACKDATA_HPP

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

#endif  // STELLARFEEDBACKDATA_HPP
