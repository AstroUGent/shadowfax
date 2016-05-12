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
 * @file ShadowfaxSnapshotWriter.hpp
 *
 * @brief Shadowfax snapshot format writer: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SHADOWFAXSNAPSHOTWRITER_HPP
#define SHADOWFAXSNAPSHOTWRITER_HPP

#include "SnapshotHandler.hpp"

/**
 * @brief SnapshotWriter implementation to write default Shadowfax snapshots
 */
class ShadowfaxSnapshotWriter : public SnapshotWriter {
  private:
    /*! \brief Timer used to quantify time spent writing snapshots */
    Timer _timer;

  public:
    ShadowfaxSnapshotWriter(std::string basename, UnitSet& units,
                            UnitSet& output_units, int lastsnap = 0);
    virtual ~ShadowfaxSnapshotWriter();

    /**
     * @brief Get the tag identifying this SnapshotWriter as a
     * ShadowfaxSnapshotWriter
     *
     * @return SHAD, because this is ShadowfaxSnapshotWriter
     */
    virtual std::string get_tag() { return "SHAD"; }

    void write_snapshot(double currentTime, ParticleVector& particles);

    virtual void dump(RestartFile& rfile);
    ShadowfaxSnapshotWriter(RestartFile& rfile, UnitSet& units,
                            UnitSet& output_units);
};

#endif  // SHADOWFAXSNAPSHOTWRITER_HPP
