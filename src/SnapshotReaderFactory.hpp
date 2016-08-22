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
 * @file SnapshotReaderFactory.hpp
 *
 * @brief Factory to generate SnapshotReader implementations from string type
 * names
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SNAPSHOTREADERFACTORY_HPP
#define SNAPSHOTREADERFACTORY_HPP

#include "Error.hpp"                    // for my_exit
#include "GadgetSnapshotReader.hpp"     // for GadgetSnapshotReader
#include "ShadowfaxSnapshotReader.hpp"  // for ShadowfaxSnapshotReader
#include <cstddef>                      // for NULL
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <string>                       // for operator==, basic_string, etc

class SnapshotReader;
class UnitSet;

/**
 * @brief Factory class for SnapshotReader
 *
 * Can be used to generate a specific implementation of SnapshotReader based on
 * a given std::string.
 * Supported snapshot types:
 *  - Shadowfax (default)
 *  - Gadget
 */
class SnapshotReaderFactory {
  public:
    /**
     * @brief Generate a SnapshotReader with the given type name
     *
     * @param name Type name of the implementation
     * @param snapname Name of the snapshot file to read
     * @param units Internal simulation UnitSet
     * @return
     */
    static SnapshotReader* generate(std::string name, std::string snapname,
                                    UnitSet& units) {
        if(name == "Shadowfax") {
            return new ShadowfaxSnapshotReader(snapname, units);
        }
        if(name == "Gadget") {
            return new GadgetSnapshotReader(snapname, units);
        }
        std::cerr << "Error! Unknown SnapshotWriter: " << name << "!"
                  << std::endl;
        my_exit();
        return NULL;
    }
};

#endif  // SNAPSHOTREADERFACTORY_HPP
