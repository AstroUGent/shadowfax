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
 * @file GadgetSnapshotWriter.cpp
 *
 * @brief A SnapshotWriter for Gadget/SWIFT/GIZMO snapshots: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "GadgetSnapshotWriter.hpp"
#include "MPIGlobal.hpp"  // for nodesize, local_rank, etc
#include "RestartFile.hpp"
#include "StateVector.hpp"       // for StateVector
#include "io/HDF5tools.hpp"      // for DOUBLE, etc
#include "io/Unit.hpp"           // for Unit
#include "io/UnitConverter.hpp"  // for UnitConverter
#include "io/UnitSet.hpp"        // for UnitSet
#include "utilities/DMParticle.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/HelperFunctions.hpp"  // for make_hdf5_file
#include "utilities/ParticleVector.hpp"   // for ParticleVector
#include "utilities/StarParticle.hpp"     // for StarParticle
#include <H5Dpublic.h>                    // for H5Dclose
#include <H5Fpublic.h>                    // for H5Fclose, H5Fcreate, etc
#include <H5Gpublic.h>                    // for H5Gclose
#include <H5Ipublic.h>                    // for hid_t
#include <H5public.h>                     // for herr_t
#include <H5version.h>                    // for H5Dopen, H5Gcreate, H5Gopen
#include <iostream>                       // for operator<<, basic_ostream, etc
#include <vector>                         // for vector
using namespace std;

/**
 * @brief Constructor
 *
 * @param basename Basic name of the snapshot, the actual name has a 3-digit
 * counter attached to it
 * @param units Internal UnitSet of the simulation
 * @param output_units Desired UnitSet for the quantities in the snapshots
 * @param lastsnap Index of the first snapshot to write (default: 0)
 * @param per_node_output Flag indicating if each node should write a separate
 * snapshot file or all nodes should write to the same file (if possible)
 */
GadgetSnapshotWriter::GadgetSnapshotWriter(std::string basename, UnitSet& units,
                                           UnitSet& output_units, int lastsnap,
                                           bool per_node_output)
        : SnapshotWriter(basename, units, output_units, lastsnap,
                         per_node_output) {}

/**
 * @brief Write a snapshot containing data from the given ParticleVector at the
 * given simulation time
 *
 * The snapshot has a filename consisting of the basename, then a 3 digit
 * counter, optionally a node rank digit, and a .hdf5-extension.
 *
 * At the end of this method, the internal counter is increased, so that two
 * consective calls of the method will produce two separate snapshot files.
 *
 * @param t Real simulation time (as opposed to the integer internal timeline)
 * @param particles Reference to the ParticleVector that should be dumped in the
 * snapshot
 * @param write_mass Should the mass be written to the snapshot?
 */
void GadgetSnapshotWriter::write_snapshot(double t, ParticleVector& particles,
                                          bool write_mass) {

    const static char* paq_names[NUM_PAQ] = PAQ_NAMES;

    string snapname;
    // if _lastsnap is negative, we do not add an index to the snapshot name
    // this is e.g. done for IC-files
    if(_lastsnap >= 0) {
        if(_per_node_output) {
            snapname = get_snapshot_name(_lastsnap, MPIGlobal::noderank,
                                         MPIGlobal::nodesize);
        } else {
            snapname = get_snapshot_name(_lastsnap);
        }
    } else {
        snapname = HelperFunctions::make_hdf5_file(_name);
    }
    cout << "Saving snapshot " << snapname << "\n" << endl;

    int rank;
    int size;
    MPI_Comm comm;
    if(_per_node_output) {
        rank = MPIGlobal::local_rank;
        size = MPIGlobal::local_size;
        comm = MPIGlobal::nodecomm;
    } else {
        rank = MPIGlobal::rank;
        size = MPIGlobal::size;
        comm = MPI_COMM_WORLD;
    }

    // units
    double uI = 1.;
    double uL = 100. * _output_units.get_length_unit().get_SI_value();
    double uM = 1000. * _output_units.get_mass_unit().get_SI_value();
    double uT = 1.;
    double utime = _output_units.get_time_unit().get_SI_value();

    UnitConverter length_converter(_units.get_length_unit(),
                                   _output_units.get_length_unit());
    UnitConverter velocity_converter(_units.get_velocity_unit(),
                                     _output_units.get_velocity_unit());
    UnitConverter mass_converter(_units.get_mass_unit(),
                                 _output_units.get_mass_unit());
    UnitConverter time_converter(_units.get_time_unit(),
                                 _output_units.get_time_unit());
    UnitConverter density_converter(_units.get_density_unit(),
                                    _output_units.get_density_unit());
    UnitConverter pressure_converter(_units.get_pressure_unit(),
                                     _output_units.get_pressure_unit());

    // collect data over all processes
    unsigned int ngaspart_loc = particles.gassize();
    unsigned int ndmpart_loc = particles.dmsize();
    unsigned int nstarpart_loc = particles.starsize();
    vector<unsigned int> gassizes(size);
    vector<unsigned int> dmsizes(size);
    vector<unsigned int> starsizes(size);
    MyMPI_Allgather(&ngaspart_loc, 1, MPI_UNSIGNED, &gassizes[0], 1,
                    MPI_UNSIGNED, comm);
    MyMPI_Allgather(&ndmpart_loc, 1, MPI_UNSIGNED, &dmsizes[0], 1, MPI_UNSIGNED,
                    comm);
    MyMPI_Allgather(&nstarpart_loc, 1, MPI_UNSIGNED, &starsizes[0], 1,
                    MPI_UNSIGNED, comm);

    // make size vectors cumulative
    for(unsigned int i = 1; i < gassizes.size(); i++) {
        gassizes[i] += gassizes[i - 1];
        dmsizes[i] += dmsizes[i - 1];
        starsizes[i] += starsizes[i - 1];
    }

    unsigned int ngas_glob, ndm_glob, nstar_glob;
    if(_per_node_output) {
        MyMPI_Allreduce(&ngaspart_loc, &ngas_glob, 1, MPI_UNSIGNED, MPI_SUM);
        MyMPI_Allreduce(&ndmpart_loc, &ndm_glob, 1, MPI_UNSIGNED, MPI_SUM);
        MyMPI_Allreduce(&nstarpart_loc, &nstar_glob, 1, MPI_UNSIGNED, MPI_SUM);
    } else {
        ngas_glob = gassizes.back();
        ndm_glob = dmsizes.back();
        nstar_glob = starsizes.back();
    }

    for(int irank = 0; irank < size; irank++) {
        if(irank == rank) {
            hid_t file;
            herr_t status;

            if(!rank) {
                // create file
                file = H5Fcreate(snapname.c_str(),
                                 HDF5constants::OVERWRITE_EXISTING_FILES,
                                 HDF5constants::DEFAULT_PROPERTY_LIST,
                                 HDF5constants::DEFAULT_PROPERTY_LIST);
            } else {
                // open file
                file = H5Fopen(snapname.c_str(), HDF5constants::READWRITE,
                               HDF5constants::DEFAULT_PROPERTY_LIST);
            }

            double box[6];
            // write header, runtimepars and units: only process 0
            if(!rank) {
                particles.get_local_header().box(box);
                vector<double> simbox(3, 0.);
#if ndim_ == 3
                simbox[0] = box[3];
                simbox[1] = box[4];
                simbox[2] = box[5];
#else
                simbox[0] = box[2];
                simbox[1] = box[3];
                simbox[2] = 0.;
#endif
                // convert units
                for(unsigned int i = 0; i < 3; i++) {
                    simbox[i] = length_converter.convert(simbox[i]);
                }

                vector<unsigned int> flag_entropy_ics(6, 0);
                vector<double> mass_table(6, 0.);
                int numfilespersnapshot;
                if(_per_node_output) {
                    numfilespersnapshot = MPIGlobal::nodesize;
                } else {
                    numfilespersnapshot = 1;
                }
                vector<unsigned int> numpart(6, 0);
                numpart[0] = gassizes.back();
                numpart[1] = dmsizes.back();
                numpart[4] = starsizes.back();
                vector<unsigned int> numpart_tot(6, 0);
                numpart_tot[0] = ngas_glob;
                numpart_tot[1] = ndm_glob;
                numpart_tot[4] = nstar_glob;
                vector<unsigned int> numpart_highword(6, 0);
                unsigned int periodic = particles.get_local_header().periodic();
                unsigned int gravity = particles.get_local_header().gravity();
                double hsoft = length_converter.convert(
                        particles.get_local_header().hsoft());

                // write header
                hid_t group = H5Gcreate(file, "Header", -1);

                HDF5tools::write_attribute_array(group, "BoxSize",
                                                 HDF5types::DOUBLE, simbox);
                HDF5tools::write_attribute_array(group, "Flag_Entropy_ICs",
                                                 HDF5types::BOOL,
                                                 flag_entropy_ics);
                HDF5tools::write_attribute_array(group, "MassTable",
                                                 HDF5types::DOUBLE, mass_table);
                HDF5tools::write_attribute_array(group, "NumFilesPerSnapshot",
                                                 HDF5types::INT,
                                                 numfilespersnapshot);
                HDF5tools::write_attribute_array(group, "NumPart_ThisFile",
                                                 HDF5types::UINT, numpart);
                HDF5tools::write_attribute_array(group, "NumPart_Total",
                                                 HDF5types::UINT, numpart_tot);
                HDF5tools::write_attribute_array(
                        group, "NumPart_Total_HighWord", HDF5types::UINT,
                        numpart_highword);
                float time = time_converter.convert(t);
                HDF5tools::write_attribute_array(group, "Time",
                                                 HDF5types::FLOAT, time);

                status = H5Gclose(group);

                // write runtime pars
                group = H5Gcreate(file, "RuntimePars", -1);

                HDF5tools::write_attribute_array(group, "PeriodicBoundariesOn",
                                                 HDF5types::BOOL, periodic);
                HDF5tools::write_attribute_array(group, "GravityOn",
                                                 HDF5types::BOOL, gravity);
                if(gravity) {
                    HDF5tools::write_attribute_array(group, "SofteningLength",
                                                     HDF5types::DOUBLE, hsoft);
                }

                status = H5Gclose(group);

                // write units
                group = H5Gcreate(file, "Units", -1);

                HDF5tools::write_attribute_array(group,
                                                 "Unit current in cgs (U_I)",
                                                 HDF5types::DOUBLE, uI);
                HDF5tools::write_attribute_array(group,
                                                 "Unit length in cgs (U_L)",
                                                 HDF5types::DOUBLE, uL);
                HDF5tools::write_attribute_array(
                        group, "Unit mass in cgs (U_M)", HDF5types::DOUBLE, uM);
                HDF5tools::write_attribute_array(
                        group, "Unit temperature in cgs (U_T)",
                        HDF5types::DOUBLE, uT);
                HDF5tools::write_attribute_array(group,
                                                 "Unit time in cgs (U_t)",
                                                 HDF5types::DOUBLE, utime);

                status = H5Gclose(group);
            } else {
                particles.get_local_header().box(box);
            }

            // gas particles
            if(gassizes.back()) {
                vector<double> coords(particles.gassize() * 3, 0.);
                vector<double> density(particles.gassize(), 0.);
                vector<float> velocity(particles.gassize() * 3, 0.);
                vector<double> internalenergy(particles.gassize(), 0.);
                vector<unsigned long> ids(particles.gassize(), 0);
                vector<vector<double> > paqs(NUM_PAQ);
                for(unsigned int i = 0; i < NUM_PAQ; i++) {
                    paqs[i].resize(particles.gassize(), 0.);
                }
                vector<double> masses;
                if(write_mass) {
                    masses.resize(particles.gassize(), 0);
                }
                for(unsigned int i = 0; i < particles.gassize(); i++) {
                    coords[3 * i] = length_converter.convert(
                            particles.gas(i)->x() - box[0]);
                    coords[3 * i + 1] = length_converter.convert(
                            particles.gas(i)->y() - box[1]);
                    coords[3 * i + 2] = length_converter.convert(
                            particles.gas(i)->z() - box[2]);
                    StateVector W = particles.gas(i)->get_Wvec();
                    density[i] = density_converter.convert(W.rho());
                    velocity[3 * i] = velocity_converter.convert(W.vx());
                    velocity[3 * i + 1] = velocity_converter.convert(W.vy());
#if ndim_ == 3
                    velocity[3 * i + 2] = velocity_converter.convert(W.vz());
#else
                    velocity[3 * i + 2] = 0.;
#endif
                    for(unsigned int j = 0; j < NUM_PAQ; j++) {
                        paqs[j][i] = W.paq(j);
                    }

                    // we assume gamma=5/3
                    if(W.rho()) {
                        internalenergy[i] = 1.5 *
                                            pressure_converter.convert(W.p()) /
                                            density_converter.convert(W.rho());
                    } else {
                        internalenergy[i] = 0.;
                    }
                    ids[i] = particles.gas(i)->id();
                    if(write_mass) {
                        masses[i] = particles.gas(i)->get_Qvec().m();
                    }
                }
                // write particle data
                if(!rank) {
                    // process 0 creates the group and datasets
                    // the datasets that are created have the size of the total
                    // data over all processes
                    hid_t group = H5Gcreate(file, "PartType0", -1);

                    hid_t dataset = HDF5tools::create_dataset_vector(
                            group, "Coordinates", HDF5types::DOUBLE,
                            gassizes.back());
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::DOUBLE, 0, coords);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_scalar(
                            group, "Density", HDF5types::DOUBLE,
                            gassizes.back());
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, 0, density);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_vector(
                            group, "Velocities", HDF5types::FLOAT,
                            gassizes.back());
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::FLOAT, 0, velocity);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_scalar(
                            group, "InternalEnergy", HDF5types::DOUBLE,
                            gassizes.back());
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, 0, internalenergy);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_scalar(
                            group, "ParticleIDs", HDF5types::ULONG,
                            gassizes.back());
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::ULONG, 0, ids);
                    status = H5Dclose(dataset);

                    for(unsigned int i = 0; i < NUM_PAQ; i++) {
                        dataset = HDF5tools::create_dataset_scalar(
                                group, paq_names[i], HDF5types::DOUBLE,
                                gassizes.back());
                        HDF5tools::write_dataset_scalar_chunk(
                                dataset, HDF5types::DOUBLE, 0, paqs[i]);
                        status = H5Dclose(dataset);
                    }

                    if(write_mass) {
                        dataset = HDF5tools::create_dataset_scalar(
                                group, "Masses", HDF5types::DOUBLE,
                                gassizes.back());
                        HDF5tools::write_dataset_scalar_chunk(
                                dataset, HDF5types::DOUBLE, 0, masses);
                        status = H5Dclose(dataset);
                    }

                    status = H5Gclose(group);
                } else {
                    // processes other than 0 open the group and datasets
                    // datasets are appended to in rank order
                    hid_t group = H5Gopen(file, "PartType0");

                    hid_t dataset = H5Dopen(group, "Coordinates");
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::DOUBLE, gassizes[rank - 1],
                            coords);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "Density");
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, gassizes[rank - 1],
                            density);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "Velocities");
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::FLOAT, gassizes[rank - 1],
                            velocity);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "InternalEnergy");
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, gassizes[rank - 1],
                            internalenergy);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "ParticleIDs");
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::ULONG, gassizes[rank - 1], ids);
                    status = H5Dclose(dataset);

                    for(unsigned int i = 0; i < NUM_PAQ; i++) {
                        dataset = H5Dopen(group, paq_names[i]);
                        HDF5tools::write_dataset_scalar_chunk(
                                dataset, HDF5types::DOUBLE, gassizes[rank - 1],
                                paqs[i]);
                        status = H5Dclose(dataset);
                    }

                    if(write_mass) {
                        dataset = H5Dopen(group, "Masses");
                        HDF5tools::write_dataset_scalar_chunk(
                                dataset, HDF5types::DOUBLE, gassizes[rank - 1],
                                masses);
                        status = H5Dclose(dataset);
                    }

                    status = H5Gclose(group);
                }
            }

            // dm particles
            if(dmsizes.back()) {
                vector<double> coords(particles.dmsize() * 3, 0.);
                vector<double> masses(particles.dmsize(), 0.);
                vector<float> velocity(particles.dmsize() * 3, 0.);
                vector<unsigned long> ids(particles.dmsize(), 0);
                for(unsigned int i = 0; i < particles.dmsize(); i++) {
                    coords[3 * i] = length_converter.convert(
                            particles.dm(i)->x() - box[0]);
                    coords[3 * i + 1] = length_converter.convert(
                            particles.dm(i)->y() - box[1]);
                    coords[3 * i + 2] = length_converter.convert(
                            particles.dm(i)->z() - box[2]);
                    masses[i] =
                            mass_converter.convert(particles.dm(i)->get_mass());
                    velocity[3 * i] =
                            velocity_converter.convert(particles.dm(i)->vx());
                    velocity[3 * i + 1] =
                            velocity_converter.convert(particles.dm(i)->vy());
                    velocity[3 * i + 2] =
                            velocity_converter.convert(particles.dm(i)->vz());
                    ids[i] = particles.dm(i)->id();
                }
                // write particle data
                if(!rank) {
                    // process 0 creates the group and datasets
                    // the datasets that are created have the size of the total
                    // data over all processes
                    hid_t group = H5Gcreate(file, "PartType1", -1);

                    hid_t dataset = HDF5tools::create_dataset_vector(
                            group, "Coordinates", HDF5types::DOUBLE,
                            dmsizes.back());
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::DOUBLE, 0, coords);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_scalar(
                            group, "Masses", HDF5types::DOUBLE, dmsizes.back());
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, 0, masses);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_vector(
                            group, "Velocities", HDF5types::FLOAT,
                            dmsizes.back());
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::FLOAT, 0, velocity);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_scalar(
                            group, "ParticleIDs", HDF5types::ULONG,
                            dmsizes.back());
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::ULONG, 0, ids);
                    status = H5Dclose(dataset);

                    status = H5Gclose(group);
                } else {
                    // processes other than 0 open the group and datasets
                    // datasets are appended to in rank order
                    hid_t group = H5Gopen(file, "PartType1");

                    hid_t dataset = H5Dopen(group, "Coordinates");
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::DOUBLE, dmsizes[rank - 1],
                            coords);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "Masses");
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, dmsizes[rank - 1],
                            masses);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "Velocities");
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::FLOAT, dmsizes[rank - 1],
                            velocity);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "ParticleIDs");
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::ULONG, dmsizes[rank - 1], ids);
                    status = H5Dclose(dataset);

                    status = H5Gclose(group);
                }
            }

            // star particles
            if(starsizes.back()) {
                vector<double> coords(particles.starsize() * 3, 0.);
                vector<double> masses(particles.starsize(), 0.);
                vector<float> velocity(particles.starsize() * 3, 0.);
                vector<unsigned long> ids(particles.starsize(), 0);
                vector<double> ages(particles.starsize(), 0.);
                for(unsigned int i = 0; i < particles.starsize(); i++) {
                    coords[3 * i] = length_converter.convert(
                            particles.star(i)->x() - box[0]);
                    coords[3 * i + 1] = length_converter.convert(
                            particles.star(i)->y() - box[1]);
                    coords[3 * i + 2] = length_converter.convert(
                            particles.star(i)->z() - box[2]);
                    masses[i] = mass_converter.convert(
                            particles.star(i)->get_mass());
                    velocity[3 * i] =
                            velocity_converter.convert(particles.star(i)->vx());
                    velocity[3 * i + 1] =
                            velocity_converter.convert(particles.star(i)->vy());
                    velocity[3 * i + 2] =
                            velocity_converter.convert(particles.star(i)->vz());
                    ids[i] = particles.star(i)->id();
                    ages[i] = time_converter.convert(
                            particles.star(i)->get_age());
                }
                // write particle data
                if(!rank) {
                    // process 0 creates the group and datasets
                    // the datasets that are created have the size of the total
                    // data over all processes
                    hid_t group = H5Gcreate(file, "PartType4", -1);

                    hid_t dataset = HDF5tools::create_dataset_vector(
                            group, "Coordinates", HDF5types::DOUBLE,
                            starsizes.back());
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::DOUBLE, 0, coords);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_scalar(
                            group, "Masses", HDF5types::DOUBLE,
                            starsizes.back());
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, 0, masses);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_vector(
                            group, "Velocities", HDF5types::FLOAT,
                            starsizes.back());
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::FLOAT, 0, velocity);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_scalar(
                            group, "ParticleIDs", HDF5types::ULONG,
                            starsizes.back());
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::ULONG, 0, ids);
                    status = H5Dclose(dataset);

                    dataset = HDF5tools::create_dataset_scalar(
                            group, "Ages", HDF5types::DOUBLE, starsizes.back());
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, 0, ages);
                    status = H5Dclose(dataset);

                    status = H5Gclose(group);
                } else {
                    // processes other than 0 open the group and datasets
                    // datasets are appended to in rank order
                    hid_t group = H5Gopen(file, "PartType4");

                    hid_t dataset = H5Dopen(group, "Coordinates");
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::DOUBLE, starsizes[rank - 1],
                            coords);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "Masses");
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, starsizes[rank - 1],
                            masses);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "Velocities");
                    HDF5tools::write_dataset_vector_chunk(
                            dataset, HDF5types::FLOAT, starsizes[rank - 1],
                            velocity);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "ParticleIDs");
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::ULONG, starsizes[rank - 1],
                            ids);
                    status = H5Dclose(dataset);

                    dataset = H5Dopen(group, "Ages");
                    HDF5tools::write_dataset_scalar_chunk(
                            dataset, HDF5types::DOUBLE, starsizes[rank - 1],
                            ages);
                    status = H5Dclose(dataset);

                    status = H5Gclose(group);
                }
            }

            status = H5Fclose(file);

            if(status < 0) {
                std::cerr << "ERROR!" << std::endl;
            }
        }
        // one process writes at a time
        MyMPI_Barrier(comm);
    }

    // put a global barrier as well
    if(_per_node_output) {
        MyMPI_Barrier();
    }

    _lastsnap++;
}

/**
 * @brief Dump this snapshotwriter to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void GadgetSnapshotWriter::dump(RestartFile& rfile) {
    SnapshotWriter::dump(rfile);
}

/**
 * @brief Restart constructor. Initialize the snapshotwriter based on the given
 * RestartFile
 *
 * As for the normal constructor, this constructor does nothing, since all
 * actions are done be the restart constructor of SnapshotWriter.
 *
 * @param rfile RestartFile to read from
 * @param units Internal UnitSet of the simulation
 * @param output_units
 */
GadgetSnapshotWriter::GadgetSnapshotWriter(RestartFile& rfile, UnitSet& units,
                                           UnitSet& output_units)
        : SnapshotWriter(rfile, units, output_units) {}
