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
 * @file GadgetSnapshotReader.cpp
 *
 * @brief A SnapshotReader for Gadget/SWIFT/GIZMO HDF5 snapshots: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "GadgetSnapshotReader.hpp"
#include "Error.hpp"
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include "io/Header.hpp"
#include "io/HDF5tools.hpp"
#include "io/UnitConverter.hpp"
#include "io/UnitSet.hpp"
#include "io/UnitSetGenerator.hpp"
#include "utilities/ParticleVector.hpp"
#include <fstream>
using namespace std;

/**
 * @brief Constructor
 *
 * Tries to open the file with the given name and spawns an error if this fails.
 *
 * @param name Filename of the snapshot to read
 * @param units Internal UnitSet of the simulation
 */
GadgetSnapshotReader::GadgetSnapshotReader(string name, UnitSet &units)
    : SnapshotReader(name, units) {
    // check if the given snapshot exists
    ifstream file(name.c_str());
    if(!file){
        cerr << "Cannot open " << name << "!" << endl;
        my_exit();
    }
}

/**
 * @brief Read the snapshot and store its contents in the given ParticleVector
 *
 * General information on the snapshot is stored in the ParticleVector Header
 * and all quantities are converted from snapshot units to the internal
 * simulation units.
 *
 * @param particles Reference to a ParticleVector to fill
 * @return A copy of the Header that was read in
 */
Header GadgetSnapshotReader::read_snapshot(ParticleVector &particles){
    Header header;

    for(int irank = 0; irank < MPIGlobal::size; irank++){
        if(irank == MPIGlobal::rank){
            hid_t file = H5Fopen(_name.c_str(), HDF5constants::READONLY,
                                 HDF5constants::DEFAULT_PROPERTY_LIST);
            herr_t status;

            // Units
            hid_t group = H5Gopen(file, "/Units");
            UnitSet *input_units;
            if(group < 0){
                // no units found, assume SI-units
                input_units = UnitSetGenerator::generate("SI");
            } else {
                double unit_mass_in_cgs =
                        HDF5tools::read_attribute_scalar<double>(
                            group, "Unit mass in cgs (U_M)", HDF5types::DOUBLE
                            );
                double unit_length_in_cgs =
                        HDF5tools::read_attribute_scalar<double>(
                            group, "Unit length in cgs (U_L)",
                            HDF5types::DOUBLE
                            );
                double unit_time_in_cgs =
                        HDF5tools::read_attribute_scalar<double>(
                            group, "Unit time in cgs (U_t)", HDF5types::DOUBLE
                            );
                Unit unit_length("length", "U_L", 0.01*unit_length_in_cgs);
                Unit unit_mass("mass", "U_M", 0.001*unit_mass_in_cgs);
                Unit unit_time("time", "U_t", unit_time_in_cgs);
                input_units = new UnitSet(unit_length, unit_mass, unit_time);
                status = H5Gclose(group);
            }

            UnitConverter length_converter(input_units->get_length_unit(),
                                           _units.get_length_unit());
            UnitConverter velocity_converter(input_units->get_velocity_unit(),
                                             _units.get_velocity_unit());
            UnitConverter mass_converter(input_units->get_mass_unit(),
                                         _units.get_mass_unit());
            UnitConverter density_converter(input_units->get_density_unit(),
                                            _units.get_density_unit());
            UnitConverter pressure_converter(input_units->get_pressure_unit(),
                                             _units.get_pressure_unit());
            UnitConverter time_converter(input_units->get_time_unit(),
                                         _units.get_time_unit());

            // Header
            group = H5Gopen(file, "/Header");

            // box
            vector<double> box =
                    HDF5tools::read_attribute_vector<double>(
                        group, "BoxSize", HDF5types::DOUBLE
                        );
            // we allow both 1 and 3 values...
            while(box.size() < 3){
                box.push_back(box[0]);
            }
#if ndim_==3
            double simbox[6];
            simbox[0] = 0.; simbox[1] = 0.; simbox[2] = 0.;
            simbox[3] = box[0]; simbox[4] = box[1]; simbox[5] = box[2];
#else
            double simbox[4];
            simbox[0] = 0.; simbox[1] = 0.;
            simbox[2] = box[0]; simbox[3] = box[1];
#endif
            // convert units
            for(unsigned int i = 0; i < ndim_+ndim_; i++){
                simbox[i] = length_converter.convert(simbox[i]);
            }
            header.set_box(simbox);

            // number of particles
            vector<unsigned int> numpart =
                    HDF5tools::read_attribute_vector<unsigned int>(
                        group, "NumPart_ThisFile", HDF5types::UINT
                        );
            header.set_ngaspart(numpart[0]);
            header.set_ndmpart(numpart[1]);

            // time - can be both float or double
            if(HDF5tools::attribute_has_type(group, "Time", HDF5types::FLOAT)){
                float time =
                        HDF5tools::read_attribute_scalar<float>(
                            group, "Time", HDF5types::FLOAT
                            );
                header.set_time(time_converter.convert(time));
            } else {
                double time =
                        HDF5tools::read_attribute_scalar<double>(
                            group, "Time", HDF5types::DOUBLE
                            );
                header.set_time(time_converter.convert(time));
            }

            status = H5Gclose(group);

            if(status < 0){
                cerr << "ERROR!" << endl;
            }

            // RuntimePars

            if(HDF5tools::exists(file, "RuntimePars")){
                group = H5Gopen(file, "/RuntimePars");

                HDF5tools::read_attribute(group, "PeriodicBoundariesOn",
                                          HDF5types::BOOL,
                                          header.get_periodic());

                if(HDF5tools::exists(group, "GravityOn")){
                    HDF5tools::read_attribute(group, "GravityOn",
                                              HDF5types::BOOL,
                                              header.get_gravity());
                } else {
                    header.set_gravity(false);
                }

                if(HDF5tools::exists(group, "SofteningLength")){
                    HDF5tools::read_attribute(group, "SofteningLength",
                                              HDF5types::DOUBLE,
                                              header.get_hsoft());
                } else {
                    // default value!
                    header.set_hsoft(0.03);
                }

                status = H5Gclose(group);

                if(status < 0){
                    cerr << "ERROR!" << endl;
                }
            } else {
                header.set_periodic(false);
                // irrelevant since this is a parameter in the parameterfile
                header.set_global_timestep(true);
                header.set_gravity(false);
                // irrelevant since gravity is false
                header.set_hsoft(0.03);
            }

            // Particle data

            // gas
            if(header.ngaspart()){
                unsigned int ngaspart = header.ngaspart();
                // npart_other is equal to the number of particles on ranks
                // lower than the current rank
                // if ngaspart is not a multiple of MPIGlobal::size, npart_local
                // will be too large for the process with the highest rank and
                // we correct for this
                unsigned int npart_local = ngaspart/MPIGlobal::size +
                        ((ngaspart%MPIGlobal::size) > 0);
                unsigned int npart_other = MPIGlobal::rank*npart_local;
                while(npart_other + npart_local > ngaspart){
                    npart_local--;
                }
                ngaspart = npart_local;
                particles.resizegas(ngaspart);
                // need: coordinates, velocities, density, pressure, id

                vector<double> coords(ngaspart*3, 0.);
                vector<double> density(ngaspart, 0.);
                vector<double> velocity(ngaspart*3, 0.);
                vector<double> pressure(ngaspart, 0.);
                vector<unsigned int> ids(ngaspart, 0);

                group = H5Gopen(file, "/PartType0");

                // coordinates
                if(HDF5tools::dataset_has_type(group, "Coordinates",
                                               HDF5types::DOUBLE)){
                    coords = HDF5tools::read_dataset_vector_chunk<double>(
                                group, npart_other, ngaspart, "Coordinates",
                                HDF5types::DOUBLE
                                );
                } else {
                    vector<float> fcoords =
                            HDF5tools::read_dataset_vector_chunk<float>(
                                group, npart_other, ngaspart, "Coordinates",
                                HDF5types::FLOAT
                                );
                    for(unsigned int i = 0; i < fcoords.size(); i++){
                        coords[i] = fcoords[i];
                    }
                }
                // convert units
                for(unsigned int i = 0; i < coords.size(); i++){
                    coords[i] = length_converter.convert(coords[i]);
                }

                // velocities
                vector<float> fvelocity =
                        HDF5tools::read_dataset_vector_chunk<float>(
                            group, npart_other, ngaspart, "Velocities",
                            HDF5types::FLOAT
                            );
                for(unsigned int i = 0; i < fvelocity.size(); i++){
                    velocity[i] = fvelocity[i];
                }
                // convert units
                for(unsigned int i = 0; i < velocity.size(); i++){
                    velocity[i] = velocity_converter.convert(velocity[i]);
                }

                // densities and pressures
                if(HDF5tools::dataset_has_type(group, "Density",
                                               HDF5types::DOUBLE)){
                    // if density is double, internal energy will be double too
                    density = HDF5tools::read_dataset_scalar_chunk<double>(
                                group, npart_other, ngaspart, "Density",
                                HDF5types::DOUBLE
                                );
                    pressure = HDF5tools::read_dataset_scalar_chunk<double>(
                                group, npart_other, ngaspart, "InternalEnergy",
                                HDF5types::DOUBLE
                                );
                } else {
                    vector<float> fdensity =
                            HDF5tools::read_dataset_scalar_chunk<float>(
                                group, npart_other, ngaspart, "Density",
                                HDF5types::FLOAT
                                );
                    vector<float> fpressure =
                            HDF5tools::read_dataset_scalar_chunk<float>(
                                group, npart_other, ngaspart, "InternalEnergy",
                                HDF5types::FLOAT
                                );
                    for(unsigned int i = 0; i < fdensity.size(); i++){
                        // density and pressure have the same size
                        density[i] = fdensity[i];
                        pressure[i] = fpressure[i];
                    }
                }
                // convert internal energy to pressure (gamma = 5./3.)
                for(unsigned int i = 0; i < pressure.size(); i++){
                    pressure[i] *= 2.*density[i]/3.;
                }
                // convert units
                for(unsigned int i = 0; i < density.size(); i++){
                    density[i] = density_converter.convert(density[i]);
                    pressure[i] = pressure_converter.convert(pressure[i]);
                }

                // particle IDs
                if(HDF5tools::dataset_has_type(group, "ParticleIDs",
                                               HDF5types::ULONG)){
                    vector<unsigned long> lids =
                            HDF5tools::read_dataset_scalar_chunk<unsigned long>(
                                group, npart_other, ngaspart, "ParticleIDs",
                                HDF5types::ULONG
                                );
                    for(unsigned int i = 0; i < lids.size(); i++){
                        // loss of precision!
                        ids[i] = lids[i];
                    }
                } else {
                    ids = HDF5tools::read_dataset_scalar_chunk<unsigned int>(
                                group, npart_other, ngaspart, "ParticleIDs",
                                HDF5types::UINT
                                );
                }

                for(unsigned int i = 0; i < ngaspart; i++){
#if ndim_==3
                    Vec position(coords[3*i], coords[3*i+1], coords[3*i+2]);
                    StateVector W(density[i], velocity[3*i], velocity[3*i+1],
                                  velocity[3*i+2], pressure[i]);
#else
                    Vec position(coords[3*i], coords[3*i+1]);
                    StateVector W(density[i], velocity[3*i], velocity[3*i+1],
                                  pressure[i]);
#endif
                    keep_inside(position, simbox);
                    particles.gas(i) = new GasParticle(position);
                    particles.gas(i)->set_W(W);
                    particles.gas(i)->set_id(ids[i]);
                    particles.gas(i)->set_v(0., 0., 0.);
                }

                status = H5Gclose(group);

                if(status < 0){
                    cerr << "ERROR!" << endl;
                }
            }

            // dm
            if(header.ndmpart()){
                unsigned int ndmpart = header.ndmpart();
                // npart_other is equal to the number of particles on ranks
                // lower than the current rank
                // if ngaspart is not a multiple of MPIGlobal::size, npart_local
                // will be too large for the process with the highest rank and
                // we correct for this
                unsigned int npart_local = ndmpart/MPIGlobal::size +
                        ((ndmpart%MPIGlobal::size) > 0);
                unsigned int npart_other = MPIGlobal::rank*npart_local;
                while(npart_other + npart_local > ndmpart){
                    npart_local--;
                }
                ndmpart = npart_local;
                particles.resizedm(ndmpart);
                // need: coordinates, velocities, density, pressure, id

                vector<double> coords(ndmpart*3, 0.);
                vector<double> masses(ndmpart, 0.);
                vector<double> velocity(ndmpart*3, 0.);
                vector<unsigned int> ids(ndmpart, 0);

                group = H5Gopen(file, "/PartType1");

                // coordinates
                if(HDF5tools::dataset_has_type(group, "Coordinates",
                                               HDF5types::DOUBLE)){
                    coords = HDF5tools::read_dataset_vector_chunk<double>(
                                group, npart_other, ndmpart, "Coordinates",
                                HDF5types::DOUBLE
                                );
                } else {
                    vector<float> fcoords =
                            HDF5tools::read_dataset_vector_chunk<float>(
                                group, npart_other, ndmpart, "Coordinates",
                                HDF5types::FLOAT
                                );
                    for(unsigned int i = 0; i < fcoords.size(); i++){
                        coords[i] = fcoords[i];
                    }
                }
                // convert units
                for(unsigned int i = 0; i < coords.size(); i++){
                    coords[i] = length_converter.convert(coords[i]);
                }

                // velocities
                vector<float> fvelocity =
                        HDF5tools::read_dataset_vector_chunk<float>(
                            group, npart_other, ndmpart, "Velocities",
                            HDF5types::FLOAT
                            );
                for(unsigned int i = 0; i < fvelocity.size(); i++){
                    velocity[i] = fvelocity[i];
                }
                // convert units
                for(unsigned int i = 0; i < velocity.size(); i++){
                    velocity[i] = velocity_converter.convert(velocity[i]);
                }

                // masses
                if(HDF5tools::dataset_has_type(group, "Masses",
                                               HDF5types::DOUBLE)){
                    masses = HDF5tools::read_dataset_scalar_chunk<double>(
                                group, npart_other, ndmpart, "Masses",
                                HDF5types::DOUBLE
                                );
                } else {
                    vector<float> fmasses =
                            HDF5tools::read_dataset_scalar_chunk<float>(
                                group, npart_other, ndmpart, "Masses",
                                HDF5types::FLOAT
                                );
                    for(unsigned int i = 0; i < fmasses.size(); i++){
                        masses[i] = fmasses[i];
                    }
                }
                // convert units
                for(unsigned int i = 0; i < masses.size(); i++){
                    masses[i] = mass_converter.convert(masses[i]);
                }

                // particle IDs
                if(HDF5tools::dataset_has_type(group, "ParticleIDs",
                                               HDF5types::ULONG)){
                    vector<unsigned long> lids =
                            HDF5tools::read_dataset_scalar_chunk<unsigned long>(
                                group, npart_other, ndmpart, "ParticleIDs",
                                HDF5types::ULONG
                                );
                    for(unsigned int i = 0; i < lids.size(); i++){
                        // loss of precision!
                        ids[i] = lids[i];
                    }
                } else {
                    ids = HDF5tools::read_dataset_scalar_chunk<unsigned int>(
                                group, npart_other, ndmpart, "ParticleIDs",
                                HDF5types::UINT
                                );
                }

                for(unsigned int i = 0; i < ndmpart; i++){
#if ndim_==3
                    Vec position(coords[3*i], coords[3*i+1], coords[3*i+2]);
#else
                    Vec position(coords[3*i], coords[3*i+1]);
#endif
                    keep_inside(position, simbox);
                    particles.dm(i) = new DMParticle(position);
                    particles.dm(i)->set_v(velocity[3*i], velocity[3*i+1],
                                           velocity[3*i+2]);
                    particles.dm(i)->set_id(ids[i]);
                    particles.dm(i)->set_mass(masses[i]);
                }

                status = H5Gclose(group);

                if(status < 0){
                    cerr << "ERROR!" << endl;
                }
            }

            status = H5Fclose(file);

            delete input_units;

            if(status < 0){
                cerr << "ERROR!" << endl;
            }
        }
        MyMPI_Barrier();
    }

    particles.get_tree().set_periodic(header.periodic());

    return header;
}
