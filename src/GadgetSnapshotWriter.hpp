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
 * @file GadgetSnapshotWriter.hpp
 *
 * @brief A SnapshotWriter for Gadget/SWIFT/GIZMO snapshots: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef GADGETSNAPSHOTWRITER_HPP
#define GADGETSNAPSHOTWRITER_HPP

#include "SnapshotHandler.hpp"  // for SnapshotWriter
#include <string>               // for string

class ParticleVector;
class RestartFile;
class UnitSet;

/**
 * @brief SnapshotWriter implementation that writes Gadget snapshots
 *
 * This type of snapshots is the same as Gadget2 snapshot type 3 and the same
 * used for Gadget2, SWIFT, GIZMO...
 */
class GadgetSnapshotWriter : public SnapshotWriter {
  public:
    GadgetSnapshotWriter(std::string basename, UnitSet& units,
                         UnitSet& output_units, int lastsnap = 0,
                         bool per_node_output = false);

    virtual ~GadgetSnapshotWriter() {}

    /**
     * @brief Get the std::string tag that distinguishes this SnapshotWriter
     * from other snapshotwriters
     *
     * @return "GADG", because this is GadgetSnapshotWriter
     */
    virtual std::string get_tag() {
        return "GADG";
    }

    virtual void write_snapshot(double t, ParticleVector& particles,
                                bool write_mass = false);

    virtual void dump(RestartFile& rfile);
    GadgetSnapshotWriter(RestartFile& rfile, UnitSet& units,
                         UnitSet& output_units);
};

#endif  // GADGETSNAPSHOTWRITER_HPP
