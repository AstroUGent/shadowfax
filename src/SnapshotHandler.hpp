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
 * @file SnapshotHandler.hpp
 *
 * @brief General interfaces for snapshot readers and writers: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SNAPSHOTHANDLER_HPP
#define SNAPSHOTHANDLER_HPP

#include "MPIGlobal.hpp"    // for rank, local_rank, etc
#include "RestartFile.hpp"  // for RestartFile
#include <ostream>          // for ifstream, ofstream
#include <stdio.h>          // for remove
#include <string>           // for allocator, string, etc

class Header;
class ParticleVector;
class UnitSet;

/**
 * @brief General interface for classes that read or write snapshots
 *
 * Stores the name of the snapshot (or the basic name for snapshot writers) and
 * a reference to the internal UnitSet.
 */
class SnapshotHandler {
  protected:
    /*! \brief Name of the snapshot. Can be either a generic name for snapshot
     *  writers or an actual filename for snapshot readers */
    std::string _name;

    /*! \brief Internal simulation UnitSet */
    UnitSet& _units;

    std::string get_snapshot_name(unsigned int nr, int rank = -1, int size = 1);

  public:
    SnapshotHandler(std::string name, UnitSet& units);
    virtual ~SnapshotHandler() {}

    void dump(RestartFile& rfile);
    SnapshotHandler(RestartFile& rfile, UnitSet& units);
};

/**
 * @brief Abstract interface for snapshot writers
 *
 * The interface implements a counter keeping the index of the last snapshot
 * that was written and keeps track of the internal units, output units and
 * basic name of the snapshots.
 * Actually writing the snapshot should be done in child classes that implement
 * this interface.
 */
class SnapshotWriter : public SnapshotHandler {
  protected:
    /*! \brief UnitSet to be used in the output file */
    UnitSet& _output_units;

    /*! \brief Counter in the name of the next snapshot that will be written */
    int _lastsnap;

    /*! \brief Flag indicating if a separate snapshot file should be written for
     *  different nodes (necessary when e.g. running OpenMPI over SSH) */
    bool _per_node_output;

  public:
    /**
     * @brief Constructor
     *
     * We check if we have to write per node snapshots or not. This needs to be
     * done when the output directory is on a node specific filesystem, in which
     * case processes on different nodes cannot access the files created by each
     * other.
     *
     * @param basename Generic name for the snapshots that will be written
     * @param units Internal simulation UnitSet
     * @param output_units UnitSet used in the output file
     * @param lastsnap Counter in the name of the first snapshot that will be
     * written
     * @param per_node_output Flag indicating if each node should write a
     * separate snapshot file or all nodes should write to the same file (if
     * possible)
     */
    SnapshotWriter(std::string basename, UnitSet& units, UnitSet& output_units,
                   int lastsnap = 0, bool per_node_output = false)
            : SnapshotHandler(basename, units), _output_units(output_units),
              _lastsnap(lastsnap), _per_node_output(per_node_output) {
        // check if we have different nodes
        if(!_per_node_output) {
            if(MPIGlobal::local_size < MPIGlobal::size) {
                // check if all nodes write to the same filesystem
                string testname = basename + ".tmp";
                if(!MPIGlobal::rank) {
                    ofstream testfile(testname);
                    testfile << "Written by process with rank 0\n";
                }
                MyMPI_Barrier();
                if(MPIGlobal::rank && !MPIGlobal::local_rank) {
                    ifstream testfile(testname);
                    if(!testfile) {
                        _per_node_output = true;
                    }
                }
                if(!MPIGlobal::rank) {
                    remove(testname.c_str());
                }

                // communicating is not so easy: we cannot broadcast from a
                // single process, since not all processes know which ranks
                // correspond to local rank 0's
                // we therefore perform an allreduce over all processes
                // since msg_send will be 0 for all processes other than local
                // rank 0's, their contributions will tell us if theirs is 0 or
                // 1
                int msg_send = _per_node_output;
                int msg_recv;
                MyMPI_Allreduce(&msg_send, &msg_recv, 1, MPI_INT, MPI_SUM);
                if(msg_recv) {
                    _per_node_output = true;
                }
            }
        }
    }

    virtual ~SnapshotWriter() {}

    /**
     * @brief Get a tag discriminating different implementations
     *
     * @return A std::string tag that is unique for every implementation
     */
    virtual std::string get_tag() = 0;

    /**
     * @brief Write a snapshot file with the given ParticleVector a the given
     * time
     *
     * @param t Current time of the simulation
     * @param particles ParticleVector to write out
     * @param write_mass Should the mass be written to the snapshot?
     */
    virtual void write_snapshot(double t, ParticleVector& particles,
                                bool write_mass = false) = 0;

    /**
     * @brief Get the current value of the snapshot counter
     *
     * @return The snapshot counter
     */
    unsigned int get_lastsnap() {
        return _lastsnap;
    }

    /**
     * @brief Dump the snapshot writer to the given RestartFile
     *
     * @param rfile RestartFile to write to
     */
    virtual void dump(RestartFile& rfile) {
        SnapshotHandler::dump(rfile);
        rfile.write(_lastsnap);
        rfile.write(_per_node_output);
    }

    /**
     * @brief Restart constructor. Initialize the snapshot writer from the given
     * RestartFile
     *
     * @param rfile RestartFile to read from
     * @param units Internal simulation UnitSet
     * @param output_units UnitSet used in the output file
     */
    SnapshotWriter(RestartFile& rfile, UnitSet& units, UnitSet& output_units)
            : SnapshotHandler(rfile, units), _output_units(output_units) {
        rfile.read(_lastsnap);
        rfile.read(_per_node_output);
    }
};

/**
 * @brief Abstract interface for snapshot readers
 *
 * This interface does nothing in itself but defines the read_snapshot function
 * that should be implemented by child classes.
 */
class SnapshotReader : public SnapshotHandler {
  public:
    /**
     * @brief Constructor
     *
     * @param name Filename to read
     * @param units Internal simulation UnitSet
     */
    SnapshotReader(std::string name, UnitSet& units)
            : SnapshotHandler(name, units) {}

    /**
     * @brief Read the snapshot and store its contents in the given
     * ParticleVector
     *
     * @param particles ParticleVector to fill
     * @param read_mass Should the mass be read from the snapshot?
     * @return Header containing general information about the snapshot
     */
    virtual Header read_snapshot(ParticleVector& particles,
                                 bool read_mass = false) = 0;
};

#endif  // SNAPSHOTHANDLER_HPP
