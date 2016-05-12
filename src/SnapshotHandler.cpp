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
 * @file SnapshotHandler.cpp
 *
 * @brief General methods for snapshot readers and writers: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "SnapshotHandler.hpp"
#include "RestartFile.hpp"
#include <sstream>
using namespace std;

/**
 * @brief Constructor
 *
 * @param name Name of the snapshots. Either a generic name for snapshot writers
 * or an actual filename for snapshot readers
 * @param units Internal simulation UnitSet
 */
SnapshotHandler::SnapshotHandler(std::string name, UnitSet& units)
        : _units(units) {
    _name = name;
}

/**
 * @brief Get the name of the snapshot with given counter
 *
 * The name is composed of the generic name of the snapshot plus a three digit
 * counter plus the .hdf5 extension.
 *
 * @param nr Counter value
 * @param rank Node rank of the snapshot, added to the snapshot name
 * @param size Size of the node space, used to determine the zero padding needed
 * for the node rank part of the name
 * @return Name of the snapshot file
 */
std::string SnapshotHandler::get_snapshot_name(unsigned int nr, int rank,
                                               int size) {
    stringstream name;
    name << _name;
    name.fill('0');
    name.width(3);
    name << nr;
    if(rank >= 0) {
        name << ".";
        int digits = 0;
        int sizecopy = size;
        int check = 1;
        while(sizecopy) {
            sizecopy /= 10;
            check *= 10;
            digits++;
        }
        check /= 10;
        // if size has one digit more than the largest rank in the range, remove
        // the extra digit
        if(check == size) {
            digits--;
        }
        name.fill('0');
        name.width(digits);
        name << rank;
    }
    name << ".hdf5";
    return name.str();
}

/**
 * @brief Dump the snapshot handler to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void SnapshotHandler::dump(RestartFile& rfile) {
    rfile.write(_name);
}

/**
 * @brief Restart constructor. Initialize the snapshot handler using the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 * @param units Internal simulation UnitSet
 */
SnapshotHandler::SnapshotHandler(RestartFile& rfile, UnitSet& units)
        : _units(units) {
    rfile.read(_name);
}
