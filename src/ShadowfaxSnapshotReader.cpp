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
 * @file ShadowfaxSnapshotReader.cpp
 *
 * @brief Shadowfax snapshot format reader: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ShadowfaxSnapshotReader.hpp"
#include "Error.hpp"
#include "io/Block.hpp"
#include "io/Input.hpp"
#include "io/UnitConverter.hpp"
#include "io/UnitSet.hpp"
#include "io/UnitSetGenerator.hpp"
#include "utilities/DMParticle.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/ParticleVector.hpp"
#include <fstream>
using namespace std;

/**
 * @brief Constructor
 *
 * @param name Name of the file to read
 * @param units Internal simulation UnitSet
 */
ShadowfaxSnapshotReader::ShadowfaxSnapshotReader(std::string name,
                                                 UnitSet& units)
    : SnapshotReader(name, units){
    // check if the given snapshot exists
    ifstream file(name.c_str());
    if(!file){
        cerr << "Cannot open " << name << "!" << endl;
        my_exit();
    }
}

/**
 * @brief Constructor
 *
 * Construct a snapshot reader from the given basename and given snapshot
 * counter.
 *
 * @param basename Basic name of the snapshot
 * @param units Internal simulation UnitSet
 * @param nr Counter of the snapshot to read
 */
ShadowfaxSnapshotReader::ShadowfaxSnapshotReader(std::string basename,
                                                 UnitSet& units,
                                                 unsigned int nr)
    : SnapshotReader(basename, units){
    _name = get_snapshot_name(nr);
}

/**
 * @brief Read the snapshot and store its contents in the given particle vector
 *
 * @param particles ParticleVector to fill
 * @return Header containing general information about the snapshot
 */
Header ShadowfaxSnapshotReader::read_snapshot(ParticleVector& particles){
    Header header;
    FileInput file(_name);
    file.read_header(header);

    // convert units in header
    // we assume the header always stores SI-units
    UnitSet* si_units = UnitSetGenerator::generate("SI");
    UnitConverter length_converter(si_units->get_length_unit(),
                                   _units.get_length_unit());
    double box[ndim_+ndim_];
    header.box(box);
    for(unsigned int i = 0; i < ndim_+ndim_; i++){
        box[i] = length_converter.convert(box[i]);
    }
    header.set_box(box);
    UnitConverter time_converter(si_units->get_time_unit(),
                                 _units.get_time_unit());
    double time = header.time();
    time = time_converter.convert(time);
    header.set_time(time);
    double hsoft = header.hsoft();
    hsoft = length_converter.convert(hsoft);
    header.set_hsoft(hsoft);
    delete si_units;

    if(header.ngaspart()){
        vector<string> headers;
        vector<Unit> units;
        headers.push_back("x");
        units.push_back(_units.get_length_unit());
        headers.push_back("y");
        units.push_back(_units.get_length_unit());
        headers.push_back("z");
        units.push_back(_units.get_length_unit());
        vector<unsigned int> dimensions(headers.size(), 1);
        Block block2("grid", headers, dimensions, units);
        file.read(block2, header.ngaspart());
        particles.resizegas(block2.number_of_lines());
        for(unsigned int i = 0; i < block2.number_of_lines(); i++){
            vector<double> coords = block2.get_line(i);
#if ndim_==3
            Vec pcoords(coords[0], coords[1], coords[2]);
#else
            Vec pcoords(coords[0], coords[1]);
#endif
            particles.gas(i) = new GasParticle(pcoords);
        }
        headers.clear();
        units.clear();
        headers.push_back("density");
        units.push_back(_units.get_density_unit());
        headers.push_back("velocity_x");
        units.push_back(_units.get_velocity_unit());
        headers.push_back("velocity_y");
        units.push_back(_units.get_velocity_unit());
#if ndim_==3
        headers.push_back("velocity_z");
        units.push_back(_units.get_velocity_unit());
#endif
        headers.push_back("pressure");
        units.push_back(_units.get_pressure_unit());
        headers.push_back("id");
        units.push_back(Unit());
        dimensions.clear();
        dimensions.resize(headers.size(), 1);
        Block block("cells", headers, dimensions, units);
        file.read(block, header.ngaspart());
        for(unsigned int i = 0; i < particles.gassize(); i++){
            vector<double> Wvec = block.get_line(i);
#if ndim_==3
            StateVector W(Wvec[0], Wvec[1], Wvec[2], Wvec[3], Wvec[4]);
#else
            StateVector W(Wvec[0], Wvec[1], Wvec[2], Wvec[3]);
#endif
            particles.gas(i)->set_W(W);
            particles.gas(i)->set_v(0.,0.,0.);
            particles.gas(i)->set_id(Wvec[ndim_+2]);
        }
    }

    if(header.ndmpart()){
        vector<string> headers;
        vector<Unit> units;
        headers.push_back("id");
        units.push_back(Unit());
        headers.push_back("mass");
        units.push_back(_units.get_mass_unit());
        headers.push_back("x");
        units.push_back(_units.get_length_unit());
        headers.push_back("y");
        units.push_back(_units.get_length_unit());
        headers.push_back("z");
        units.push_back(_units.get_length_unit());
        headers.push_back("vx");
        units.push_back(_units.get_velocity_unit());
        headers.push_back("vy");
        units.push_back(_units.get_velocity_unit());
        headers.push_back("vz");
        units.push_back(_units.get_velocity_unit());
        vector<unsigned int> dimensions(headers.size(), 1);
        Block block("DM", headers, dimensions, units, header.ndmpart());
        file.read(block, header.ndmpart());
        particles.resizedm(block.number_of_lines());
        for(unsigned int i = 0; i < particles.dmsize(); i++){
            vector<double> row = block.get_line(i);
            unsigned long id = row[0];
            double mass = row[1];
#if ndim_==3
            Vec position(row[2], row[3], row[4]);
#else
            Vec position(row[2], row[3]);
#endif
            particles.dm(i) = new DMParticle(position);
            particles.dm(i)->set_v(row[5], row[6], row[7]);
            particles.dm(i)->set_id(id);
            particles.dm(i)->set_mass(mass);
        }
    }

    particles.get_tree().set_periodic(header.periodic());
    return header;
}
