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
 * @file ShadowfaxSnapshotWriter.cpp
 *
 * @brief Shadowfax snapshot format writer: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ShadowfaxSnapshotWriter.hpp"
#include "Error.hpp"      // for my_exit
#include "MPIGlobal.hpp"  // for rank
#include "RestartFile.hpp"
#include "StateVector.hpp"                // for StateVector
#include "Vec.hpp"                        // for Vec
#include "io/Block.hpp"                   // for Block
#include "io/FileOutput.hpp"              // for FileOutput
#include "io/Header.hpp"                  // for Header
#include "io/Unit.hpp"                    // for Unit, operator/
#include "io/UnitConverter.hpp"           // for UnitConverter
#include "io/UnitSet.hpp"                 // for UnitSet
#include "io/UnitSetGenerator.hpp"        // for UnitSetGenerator
#include "utilities/DMParticle.hpp"       // for DMParticle
#include "utilities/GasParticle.hpp"      // for GasParticle
#include "utilities/HelperFunctions.hpp"  // for make_hdf5_file
#include "utilities/ParticleVector.hpp"   // for ParticleVector
#include <iostream>                       // for operator<<, basic_ostream, etc
#include <vector>                         // for vector
using namespace std;

/**
 * @brief Constructor
 *
 * @param basename Basic name of the snapshots
 * @param units Internal simulation UnitSet
 * @param output_units Output UnitSet
 * @param lastsnap Counter of the first snapshot to write
 */
ShadowfaxSnapshotWriter::ShadowfaxSnapshotWriter(std::string basename,
                                                 UnitSet& units,
                                                 UnitSet& output_units,
                                                 int lastsnap)
        : SnapshotWriter(basename, units, output_units, lastsnap) {
    if(_per_node_output) {
        cerr << "System setup requires per node output, but this is not "
                "supported by ShadowfaxSnapshotWriter. Consider using another "
                "snapshot format!"
             << endl;
        my_exit();
    }
}

/**
 * @brief Destructor
 *
 * Print out the timer to the stdout.
 */
ShadowfaxSnapshotWriter::~ShadowfaxSnapshotWriter() {
    cout << "Time spent writing snapshots: " << _timer.value() << "s" << endl;
}

/**
 * @brief Write a snapshot with the data contained in the given particle vector
 * and with the given simulation time
 *
 * @param currentTime Current time of the simulation, written to the snapshot
 * header
 * @param particles ParticleVector to write out
 * @param write_mass Should the mass be written to the snapshot?
 */
void ShadowfaxSnapshotWriter::write_snapshot(double currentTime,
                                             ParticleVector& particles,
                                             bool write_mass) {
    _timer.start();
    string snapname;
    // if _lastsnap is negative, we do not add an index to the snapshot name
    // this is e.g. done for IC-files
    if(_lastsnap >= 0) {
        snapname = get_snapshot_name(_lastsnap);
    } else {
        snapname = HelperFunctions::make_hdf5_file(_name);
    }
    cout << "Saving snapshot " << snapname << "\n" << endl;

    // quantities in the Header always have SI units!
    UnitSet* si_units = UnitSetGenerator::generate("SI");
    double t = currentTime;
    FileOutput file(snapname);
    UnitConverter time_converter(_units.get_time_unit(),
                                 si_units->get_time_unit());
    t = time_converter.convert(t);
    particles.get_header().set_time(t);
    // convert box to correct units
    double box[6] = {0.};
    particles.get_header().box(box);
    UnitConverter length_converter(_units.get_length_unit(),
                                   si_units->get_length_unit());
    for(unsigned int i = 0; i < 6; i++) {
        box[i] = length_converter.convert(box[i]);
    }
    particles.get_header().set_box(box);
    // convert softening length to correct units
    if(particles.dmsize()) {
        particles.get_header().set_hsoft(
                length_converter.convert(particles.get_header().hsoft()));
    }
    file.write_header(particles.get_header());

    delete si_units;

    if(particles.gassize()) {
        vector<string> headers;
        vector<Unit> units;
        headers.push_back("id");
        units.push_back(Unit());
        headers.push_back("density");
        units.push_back(_output_units.get_density_unit());
        headers.push_back("velocity_x");
        units.push_back(_output_units.get_velocity_unit());
        headers.push_back("velocity_y");
        units.push_back(_output_units.get_velocity_unit());
#if ndim_ == 3
        headers.push_back("velocity_z");
        units.push_back(_output_units.get_velocity_unit());
#endif
        headers.push_back("pressure");
        units.push_back(_output_units.get_pressure_unit());
        headers.push_back("timestep");
        units.push_back(_output_units.get_time_unit());
        headers.push_back("acceleration_x");
        units.push_back(_output_units.get_velocity_unit() /
                        _output_units.get_time_unit());
        headers.push_back("acceleration_y");
        units.push_back(_output_units.get_velocity_unit() /
                        _output_units.get_time_unit());
#if ndim_ == 3
        headers.push_back("acceleration_z");
        units.push_back(_output_units.get_velocity_unit() /
                        _output_units.get_time_unit());
#endif
        headers.push_back("mass");
        units.push_back(_output_units.get_mass_unit());
        vector<unsigned int> dimensions(headers.size(), 1);
        Block block("cells", headers, dimensions, units, particles.gassize());
        vector<double> data(headers.size());
        UnitConverter density_converter(_units.get_density_unit(),
                                        _output_units.get_density_unit());
        UnitConverter velocity_converter(_units.get_velocity_unit(),
                                         _output_units.get_velocity_unit());
        UnitConverter pressure_converter(_units.get_pressure_unit(),
                                         _output_units.get_pressure_unit());
        UnitConverter time_converter(_units.get_time_unit(),
                                     _output_units.get_time_unit());
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            data[0] = particles.gas(i)->id();
            StateVector W = particles.gas(i)->get_Wvec();
            data[1] = density_converter.convert(W.rho());
            data[2] = velocity_converter.convert(W.vx());
            data[3] = velocity_converter.convert(W.vy());
#if ndim_ == 3
            data[4] = velocity_converter.convert(W.vz());
#endif
            data[ndim_ + 2] = pressure_converter.convert(W.p());
            data[ndim_ + 3] =
                    time_converter.convert(particles.gas(i)->get_timestep());
            data[ndim_ + 4] = density_converter.convert(
                    particles.gas(i)->get_gravitational_acceleration().x());
            data[ndim_ + 5] = density_converter.convert(
                    particles.gas(i)->get_gravitational_acceleration().y());
#if ndim_ == 3
            data[ndim_ + 6] = density_converter.convert(
                    particles.gas(i)->get_gravitational_acceleration().z());
#endif
            data[ndim_ + ndim_ + 4] =
                    density_converter.convert(particles.gas(i)->get_mass());
            block.add_data(data);
        }
        file.write(block);
        headers.clear();
        units.clear();
        headers.push_back("x");
        units.push_back(_output_units.get_length_unit());
        headers.push_back("y");
        units.push_back(_output_units.get_length_unit());
        headers.push_back("z");
        units.push_back(_output_units.get_length_unit());
        dimensions.clear();
        dimensions.resize(headers.size(), 1);
        Block block2("grid", headers, dimensions, units, particles.gassize());
        data.resize(headers.size());
        UnitConverter length_converter(_units.get_length_unit(),
                                       _output_units.get_length_unit());
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            Vec& position = particles.gas(i)->get_position();
            data[0] = length_converter.convert(position.x());
            data[1] = length_converter.convert(position.y());
#if ndim_ == 3
            data[2] = length_converter.convert(position.z());
#else
            data[2] = 0.;
#endif
            block2.add_data(data);
        }
        file.write(block2);
    }
    if(particles.dmsize()) {
        vector<string> headers;
        vector<Unit> units;
        headers.push_back("id");
        units.push_back(Unit());
        headers.push_back("mass");
        units.push_back(_output_units.get_mass_unit());
        headers.push_back("x");
        units.push_back(_output_units.get_length_unit());
        headers.push_back("y");
        units.push_back(_output_units.get_length_unit());
        headers.push_back("z");
        units.push_back(_output_units.get_length_unit());
        headers.push_back("vx");
        units.push_back(_output_units.get_velocity_unit());
        headers.push_back("vy");
        units.push_back(_output_units.get_velocity_unit());
        headers.push_back("vz");
        units.push_back(_output_units.get_velocity_unit());
        headers.push_back("process");
        units.push_back(Unit());
        vector<unsigned int> dimensions(headers.size(), 1);
        Block block("DM", headers, dimensions, units, particles.dmsize());
        vector<double> data(headers.size());
        UnitConverter mass_converter(_units.get_mass_unit(),
                                     _output_units.get_mass_unit());
        UnitConverter length_converter(_units.get_length_unit(),
                                       _output_units.get_length_unit());
        UnitConverter velocity_converter(_units.get_velocity_unit(),
                                         _output_units.get_velocity_unit());
        for(unsigned int i = 0; i < particles.dmsize(); i++) {
            double mass = particles.dm(i)->get_mass();
            Vec position = particles.dm(i)->get_position();
            Vec velocity = particles.dm(i)->get_velocity();
            data[0] = particles.dm(i)->id();
            data[1] = mass_converter.convert(mass);
            data[2] = length_converter.convert(position.x());
            data[3] = length_converter.convert(position.y());
#if ndim_ == 3
            data[4] = length_converter.convert(position.z());
#endif
            data[5] = velocity_converter.convert(velocity.x());
            data[6] = velocity_converter.convert(velocity.y());
#if ndim_ == 3
            data[7] = velocity_converter.convert(velocity.z());
#endif
            data[8] = MPIGlobal::rank;
            block.add_data(data);
        }
        file.write(block);
    }
    _lastsnap++;
    _timer.stop();
}

/**
 * @brief Dump the snapshot writer to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void ShadowfaxSnapshotWriter::dump(RestartFile& rfile) {
    SnapshotWriter::dump(rfile);
    _timer.dump(rfile);
}

/**
 * @brief Restart constructor. Initialize the snapshot writer from the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 * @param units Internal simulation UnitSet
 * @param output_units Output UnitSet
 */
ShadowfaxSnapshotWriter::ShadowfaxSnapshotWriter(RestartFile& rfile,
                                                 UnitSet& units,
                                                 UnitSet& output_units)
        : SnapshotWriter(rfile, units, output_units), _timer(rfile) {}
