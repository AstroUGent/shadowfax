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
 * @file GadgetSnapshotReader.hpp
 *
 * @brief A SnapshotReader for Gadget/SWIFT/GIZMO HDF5 snapshots: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef GADGETSNAPSHOTREADER_HPP
#define GADGETSNAPSHOTREADER_HPP

#include "SnapshotHandler.hpp"  // for SnapshotReader
#include "Vec.hpp"              // for Vec
#include "io/Header.hpp"        // for Header
#include <string>               // for string

class ParticleVector;
class UnitSet;

/**
 * @brief SnapshotReader specialization to read Gadget snapshots
 */
class GadgetSnapshotReader : public SnapshotReader {
  private:
    /**
     * @brief Make sure the given coordinates are inside the given box by moving
     * them by a box length
     *
     * @param position Coordinates that should be inside the box
     * @param box Origin and sides of a box in 3 dimensions
     */
    void keep_inside(Vec& position, double* box) {
        for(unsigned int i = 0; i < ndim_; i++) {
            if(position[i] < box[i]) {
                position[i] += box[i + ndim_];
            }
            if(position[i] > box[i] + box[i + ndim_]) {
                position[i] -= box[i + ndim_];
            }
        }
    }

  public:
    GadgetSnapshotReader(std::string name, UnitSet& units);

    virtual Header read_snapshot(ParticleVector& particles);
};

#endif  // GADGETSNAPSHOTREADER_HPP
