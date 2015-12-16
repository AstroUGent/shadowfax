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
 * @file ShadowfaxSnapshotReader.hpp
 *
 * @brief Shadowfax snapshot format reader: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SHADOWFAXSNAPSHOTREADER_HPP
#define SHADOWFAXSNAPSHOTREADER_HPP

#include "SnapshotHandler.hpp"

/**
 * @brief Specialization of the SnapshotReader interface for Shadowfax snapshots
 *
 * The Shadowfax snapshot format is the native and default snapshot type.
 */
class ShadowfaxSnapshotReader : public SnapshotReader{
public:
    ShadowfaxSnapshotReader(std::string name, UnitSet& units);
    ShadowfaxSnapshotReader(std::string basename, UnitSet& units,
                            unsigned int nr);
    ~ShadowfaxSnapshotReader(){}

    Header read_snapshot(ParticleVector& particles);
};

#endif // SHADOWFAXSNAPSHOTREADER_HPP
