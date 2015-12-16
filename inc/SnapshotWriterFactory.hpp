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
 * @file SnapshotWriterFactory.hpp
 *
 * @brief Factory to generate SnapshotWriter instances from string type names
 * and to load and dump them from/to restart files
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SNAPSHOTWRITERFACTORY_HPP
#define SNAPSHOTWRITERFACTORY_HPP

#include <istream>
#include "SnapshotHandler.hpp"
#include "GadgetSnapshotWriter.hpp"
#include "ShadowfaxSnapshotWriter.hpp"
#include "RestartFile.hpp"
#include "Error.hpp"

/**
 * @brief Factory to generate a SnapshotWriter from a given std::string and to
 * dump or read a SnapshotWriter to/from a RestartFile
 */
class SnapshotWriterFactory{
public:
    /**
     * @brief Load a SnapshotWriter from the given RestartFile
     *
     * We first read the type name from the RestartFile and then call the
     * corresponding restart constructor.
     *
     * @param rfile RestartFile to read from
     * @param units Internal simulation UnitSet
     * @param output_units Output UnitSet
     * @return Pointer to a SnapshotWriter implementation instance
     */
    static SnapshotWriter* load(RestartFile &rfile, UnitSet &units,
                                UnitSet &output_units){
        std::string tag;
        rfile.read(tag);
        if(tag == "SHAD"){
            return new ShadowfaxSnapshotWriter(rfile, units, output_units);
        }
        if(tag == "GADG"){
            return new GadgetSnapshotWriter(rfile, units, output_units);
        }
        return NULL;
    }

    /**
     * @brief Dump the given SnapshotWriter implementation to the given
     * RestartFile
     *
     * We first write the type name to the RestartFile and then call the dump
     * method of the implementation.
     *
     * @param rfile RestartFile to write to
     * @param writer SnapshotWriter instance to dump
     */
    static void dump(RestartFile &rfile, SnapshotWriter *writer){
        std::string tag = writer->get_tag();
        rfile.write(tag);
        writer->dump(rfile);
    }

    /**
     * @brief Generate a SnapshotWriter instance with the given type name
     *
     * @param name Type name of the implementation
     * @param basename Basic name of the snapshot files
     * @param units Internal simulation UnitSet
     * @param output_units Output UnitSet
     * @param lastsnap Counter of the first snapshot to write
     * @param per_node_output Flag indicating if each node should write a
     * separate snapshot file or all nodes should write to the same file (if
     * possible)
     * @return Pointer to a SnapshotWriter implementation instance
     */
    static SnapshotWriter* generate(std::string name, std::string basename,
                                    UnitSet& units, UnitSet& output_units,
                                    int lastsnap = 0,
                                    bool per_node_output = false){
        if(name == "Shadowfax"){
            return new ShadowfaxSnapshotWriter(basename, units, output_units,
                                               lastsnap);
        }
        if(name == "Gadget"){
            return new GadgetSnapshotWriter(basename, units, output_units,
                                            lastsnap, per_node_output);
        }
        std::cerr << "Error! Unknown SnapshotWriter: " << name << "!"
                  << std::endl;
        my_exit();
        return NULL;
    }
};

#endif // SNAPSHOTWRITERFACTORY_HPP
