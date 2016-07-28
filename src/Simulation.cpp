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
 * @file Simulation.cpp
 *
 * @brief Main simulation program: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Simulation.hpp"
#include "DelCont.hpp"
#include "ExArith.hpp"
#include "GravityWalker.hpp"
#include "LogFiles.hpp"
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include "ParameterFile.hpp"
#include "Physics.hpp"
#include "ProgramLog.hpp"
#include "RestartFile.hpp"
#include "SnapshotReaderFactory.hpp"
#include "TimeLine.hpp"
#include "Vec.hpp"
#include "VorCell.hpp"
#include "VorTessManager.hpp"
#include "io/Block.hpp"
#include "io/HDF5tools.hpp"
#include "io/Header.hpp"
#include "io/Input.hpp"
#include "io/Output.hpp"
#include "io/UnitSet.hpp"
#include "io/UnitSetGenerator.hpp"
#include "riemann/ApproximateSolver.hpp"
#include "riemann/RiemannSolver.hpp"
#include "riemann/RiemannSolverFactory.hpp"
#include "utilities/DMParticle.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/HelperFunctions.hpp"
#include "utilities/ParticleVector.hpp"
#include "utilities/Timer.hpp"
#include "utilities/Tree.hpp"
#include "utilities/TreeWalker.hpp"
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <vector>

//#define VARIABLE_SOFTENING
//#define ENTROPY

#define SIMULATION_DEFAULT_OUTPUTDIR "."
#define SIMULATION_DEFAULT_RESTARTTIME 3600.
#define SIMULATION_DEFAULT_UNITS "SI"
#define SIMULATION_DEFAULT_ICTYPE "Gadget"
#define SIMULATION_DEFAULT_ICNAME "icfile.hdf5"
#define SIMULATION_DEFAULT_MAXMEMORY "1 GB"
#define SIMULATION_DEFAULT_SOFTENING 0.03
#define SIMULATION_DEFAULT_VORONOITOLERANCE 1.e-9
#define SIMULATION_DEFAULT_GRAVALPHA 0.005

// for debug purposes
#include <VorGen.hpp>
using namespace std;

/**
 * @brief The main simulation routine.
 *
 * The routine consists of several consecutive steps:
 *  - read the parameterfile
 *  - set up the RiemannSolver, ParticleVector, TimeLine and initial VorTess
 *  - calculate initial hydrodynamical variables and gradients
 *  - calculate the initial timesteps
 *  - the main simulation loop
 *
 * The main simulation loop itself consists of
 *  - calculating the fluxes around the active cells
 *  - drifting all particles (including the passive ones)
 *  - setting up a new VorTess
 *  - updating hydrodynamical variables and gradients
 *  - updating the timesteps
 *
 * The only parameters to this routine are the command line parameters. If the
 * second command line parameter contains a string, this string is expected to
 * be the path to a valid parameterfile. If no second command line parameter is
 * provided, the default parameterfile "default.ini" is used.
 *
 * @param argc The number of command line arguments
 * @param argv The list of command line arguments
 */
Simulation::Simulation(int argc, char** argv) {
    _starttimer = new Timer();
    // initialize variables
    _periodic = false;

    // suppress output for processes other than the process with rank 0
    if(MPIGlobal::rank) {
        cout.rdbuf(NULL);
    }

    // turn off default HDF5 error handling
    HDF5tools::turn_off_error_handling();

    // reroute cerr output to a different file for every process
    //#warning DISABLE FOR PRODUCTION RUNS!!
    //    stringstream errname;
    //    errname << "error_" << MPIGlobal::rank << ".txt";
    //    ofstream errfile(errname.str());
    //    cerr.rdbuf(errfile.rdbuf());

    // print some output to the terminal
    print_welcome_message();

    if(MPIGlobal::size < 2) {
        // To be fair: a very similar message is written by Gadget-2 :p
        cout << "Running on 1 process\n";
        cout << "Since this code is meant to be run on multiple processes, "
                "there might be some unnecessary overhead\n"
             << endl;
    } else {
        cout << "Running on " << MPIGlobal::size << " processes\n" << endl;
    }

    static struct option long_options[] = {
            {"restart", required_argument, NULL, 'r'},
            {"params", required_argument, NULL, 'p'},
            {0, 0, 0, 0}};

    std::string filename;
    bool restartflag = false;
    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":r:p:", long_options, NULL)) != -1) {
        switch(c) {
            case 'r':
                restartflag = true;
                filename = optarg;
                break;
            case 'p':
                filename = optarg;
                break;
            case ':':
                cerr << "Error! No ";
                if(optopt == 'p') {
                    cerr << "parameterfile ";
                } else {
                    cerr << "restartfile ";
                }
                cerr << "provided. Crashing..." << endl;
                return;
        }
    }

    // initialize the data (based on a restart file or a parameterfile)
    if(!filename.size()) {
        cout << "No parameterfile provided, using default parameterfile "
                "default.ini"
             << endl;
        initialize("default.ini");
    } else {
        if(restartflag) {
            restart(filename);
        } else {
            initialize(filename);
        }
    }
    _starttimer->stop();

    // start the main simulation loop
    main_loop();
}

/**
 * @brief Destructor
 *
 * Clean up all simulation components and print the timers to the stdout.
 */
Simulation::~Simulation() {
    delete _physics;
    delete _parameterfile;
    delete _logfiles;
    delete _solver;
    delete _simulation_units;
    delete _output_units;
    delete _box;
    if(_particles->gassize()) {
        delete _voronoi;
    }
    delete _particles;
    delete _timeline;
    cout << "Initialization time: " << _starttimer->value() << "s" << endl;
    cout << "Timestep calculation: " << _timesteptimer->value() << "s" << endl;
    cout << "MPI communication: " << _mpitimer->value() << "s" << endl;
    cout << "Additional hydro: " << _hydrotimer->value() << "s" << endl;
    cout << "Gravity: " << _gravitytimer->value() << "s" << endl;
    delete _starttimer;
    delete _timesteptimer;
    delete _mpitimer;
    delete _hydrotimer;
    delete _gravitytimer;
    cout << "Total program time: " << _totaltimer.stop() << " seconds" << endl;
}

/**
  * @brief Load a snapshot file and initialize a ParticleVector using its
  * contents
  *
  * @param cells An empty ParticleVector that will be filled
  * @param ictype The type of the initial conditions
  * @param filename The name of the snapshot file containing the particles that
  * should be put in the ParticleVector
  * @param simulation_units Internal Units used
  * @return The time stamp of the snapshot file
  */
double Simulation::load(ParticleVector& cells, string ictype, string filename,
                        UnitSet& simulation_units) {
    LOGS("Start loading particles...");

    cout << "Read initial condition file " << filename << endl;
    SnapshotReader* reader =
            SnapshotReaderFactory::generate(ictype, filename, simulation_units);
    Header header = reader->read_snapshot(cells);
    cout << "Found " << header.npart()
         << " particles (gas: " << header.ngaspart()
         << ", DM: " << header.ndmpart() << ")" << endl;

    _periodic = header.periodic();
    cells.set_periodic(_periodic);

    double box[ndim_ + ndim_] = {0.};
    header.box(box);
    Vec center;
    Vec sides;
    for(unsigned int i = ndim_; i--;) {
        center[i] = box[i] + 0.5 * box[ndim_ + i];
        sides[i] = box[ndim_ + i];
    }
    delete _box;
    _box = new RectangularBox(center, sides);
    cells.set_container(*_box);
    cells.get_header().set_ngaspart(header.ngaspart());
    cells.get_header().set_ndmpart(header.ndmpart());

    cout << "Set up ";
    if(_periodic) {
        cout << "periodic";
    } else {
        cout << "reflective";
    }
    cout << " simulation box with center (";
#if ndim_ == 3
    cout << center[0] << "," << center[1] << "," << center[2];
#else
    cout << center[0] << "," << center[1];
#endif
    cout << ") and sides (" << sides[0] << "," << sides[1];
#if ndim_ == 3
    cout << "," << sides[2];
#endif
    cout << ")" << endl;
    cout << "Initial time is " << header.time() << "\n" << endl;

    delete reader;

    LOGS("Particles loaded");

    return header.time();
}

/**
 * @brief Main simulation loop
 *
 * This method does the actual integration.  It should probably be documented
 * better than this, but then again, I should probably also write a paper that
 * does the same thing and that one has priority.
 */
void Simulation::main_loop() {
    _mpitimer = new Timer();
    _hydrotimer = new Timer();
    _gravitytimer = new Timer();

    // we have to calculate the potential before we calculate the time steps,
    // since active particles are defined by the endtime being equal to 0,
    // which is no longer true after the time steps are calculated
    _gravitytimer->start();
    if(_gravity) {
        LOGS("Start calculating gravitational potential");
        _particles->get_tree().walk_tree<PotentialWalker>(*_particles, true,
                                                          true, 0);
        // multiply with G
        for(unsigned int i = 0; i < _particles->gassize(); i++) {
            GasParticle* p = _particles->gas(i);
            p->set_gravitational_potential(
                    _physics->get_gravitational_constant() *
                    p->get_gravitational_potential());
        }
        for(unsigned int i = 0; i < _particles->dmsize(); i++) {
            DMParticle* p = _particles->dm(i);
            p->set_gravitational_potential(
                    _physics->get_gravitational_constant() *
                    p->get_gravitational_potential());
        }
        LOGS("Finished calculating gravitational potential");
    }
    _gravitytimer->stop();

    // we now have all information needed to calculate the timesteps
    _timesteptimer = new Timer();
    unsigned long dt = _timeline->calculate_timestep();
    _timesteptimer->stop();

    // the actual program loop
    cout << "Starting integration loop\n" << endl;
    // timeline.step_forward() returns true as long as the simulation has not
    // reached its end time
    // snapshots are written by the timeline at appropriate times
    while(_timeline->step_forward()) {
        cout << "t = " << _timeline->get_time()
             << ", dt = " << _timeline->get_realtime_interval(dt) << endl;

        unsigned long currentTime = _timeline->get_integertime();
        if(_particles->gassize()) {
            for(unsigned int i = 0; i < _particles->gassize(); i++) {
                if(_particles->gas(i)->get_starttime() == currentTime) {
                    StateVector Wgas = _particles->gas(i)->get_Wvec();

                    Wgas += 0.5 *
                            _timeline->get_realtime_interval(
                                    _particles->gas(i)->get_endtime() -
                                    _particles->gas(i)->get_starttime()) *
                            _particles->gas(i)
                                    ->get_gravitational_acceleration();
                    _particles->gas(i)->set_W(Wgas);

                    Vec v = _voronoi->get_velocity(i);

                    if(!_periodic) {
                        // make sure the particle stays inside the box
                        double dt = _timeline->get_realtime_interval(
                                _particles->gas(i)->get_endtime() -
                                _particles->gas(i)->get_starttime());
                        Vec newx = _particles->gas(i)->get_position() + dt * v;
                        if(!_particles->get_container().inside(newx)) {
                            // if the particle risks to move outside the box,
                            // we do not move it along with the flow
                            v[0] -= Wgas[1];
                            v[1] -= Wgas[2];
#if ndim_ == 3
                            v[2] -= Wgas[3];
#endif
                        }
                    }

                    _particles->gas(i)->set_velocity(v);
                    if(_gravity) {
                        _particles->gas(i)->apply_gravity(
                                0.5 *
                                _timeline->get_realtime_interval(
                                        _particles->gas(i)->get_endtime() -
                                        _particles->gas(i)->get_starttime()));
                    }
                }
            }
            LOGS("Did first kick");
            // local Ws changed: communicate to other processes
            _mpitimer->start();
            _voronoi->update_Ws();
            _mpitimer->stop();

            for(unsigned int i = 0; i < _particles->gassize(); i++) {
                if(_particles->gas(i)->get_starttime() == currentTime) {
                    StateVector gradients[3];
                    _voronoi->estimate_gradients(i, gradients);
                    _particles->gas(i)->set_gradients(gradients);
                }
            }
            LOGS("Estimated gradients");
            // exchange gradients between processes. This also exchanges
            // centroids
            _mpitimer->start();
            _voronoi->update_gradients();
            _mpitimer->stop();

            // exchange timestep information between processes
            _mpitimer->start();
            _voronoi->update_dts(currentTime);
            _mpitimer->stop();

            // do the actual integration by calculating fluxes between cells
            LOGS("Starting flux exchange");
            _voronoi->hydro(*_timeline, *_solver, *_particles);
            LOGS("Finished flux exchange");
            // exchange fluxes between processes (and gravititional fluxes)
            _mpitimer->start();
            _voronoi->update_dQs();
            _mpitimer->stop();
        }
        _hydrotimer->start();
        for(unsigned int i = _particles->gassize(); i--;) {
            _particles->gas(i)->set_vorgen(NULL);
        }
        for(unsigned int i = 0; i < _particles->gassize(); i++) {
#ifndef STATIC
            // move the mesh generators
            _particles->gas(i)->drift(_timeline->get_realtime_interval(dt));
#endif
            // add the calculated flux for this cell to its conserved quantities
            _particles->gas(i)->update_Q();
        }
        LOGS("Did drift");
        _hydrotimer->stop();
        // update the positions used for the mesh
        // for the AdaptiveVorTess, this effectively moves the cell centers
        // for the normal VorTess, this only makes sure that all particles
        // stay inside the periodic box
        if(_particles->gassize()) {
            _voronoi->update_positions(_periodic);
        }
        _gravitytimer->start();
        if(_gravity) {
            // half a kick and drift for the N-body particles
            for(unsigned int i = _particles->dmsize(); i--;) {
                if(_particles->dm(i)->get_starttime() == currentTime) {
                    _particles->dm(i)->accelerate(
                            _particles->dm(i)
                                    ->get_gravitational_acceleration() *
                            0.5 *
                            _timeline->get_realtime_interval(
                                    _particles->dm(i)->get_endtime() -
                                    _particles->dm(i)->get_starttime()));
                }
                _particles->dm(i)->move(_timeline->get_realtime(dt));
            }
        }
        _gravitytimer->stop();
        LOGS("Did second kick");
        currentTime += dt;
        _timeline->set_timestep(dt);

        // this does a new domain decomposition and recalculates the tree
        _particles->sort();
        unsigned int count = 1;
        if(_particles->gassize()) {
            LOGS("Recalculating grid...");
            // we need to do the reset after the sort, since we need the new
            // number of particles on this process...
            // recalculate the voronoi mesh
            cout << "setting up new mesh" << endl;
            _voronoi->reset(&_particles->get_container(), _periodic);
            count = _voronoi->update(currentTime);
            cout << count << " active particles\n" << endl;
            count = _particles->gassize() - count;
            unsigned int count_glob;
            MyMPI_Allreduce(&count, &count_glob, 1, MPI_INT, MPI_SUM);
            count = count_glob;
            LOGS("Finished recalculating grid.");
        }

#ifdef VARIABLE_SOFTENING
        _gravitytimer->start();
        if(_gravity) {
            // we set the softening lenghts of the gas before we update the tree
            // values
            for(unsigned int i = _particles->gassize(); i--;) {
                if(_particles->gas(i)->get_endtime() == currentTime) {
                    double hsoft = 0.75 * _voronoi->get_volume(i) / M_PI;
                    hsoft = 1.5 * cbrt(hsoft);
                    _particles->gas(i)->set_hsoft(hsoft);
                }
            }
        }
        _gravitytimer->stop();
#endif

        _particles->get_tree().set_velocities();
        if(MPIGlobal::size > 1) {
            _particles->get_tree().exchange_pseudonodes();
        }

        _gravitytimer->start();
        if(_gravity) {
            LOGS("Start calculating gravitational potential");
            _particles->get_tree().walk_tree<PotentialWalker>(
                    *_particles, true, true, currentTime);
            // multiply with G
            for(unsigned int i = 0; i < _particles->gassize(); i++) {
                GasParticle* p = _particles->gas(i);
                p->set_gravitational_potential(
                        _physics->get_gravitational_constant() *
                        p->get_gravitational_potential());
            }
            for(unsigned int i = 0; i < _particles->dmsize(); i++) {
                DMParticle* p = _particles->dm(i);
                p->set_gravitational_potential(
                        _physics->get_gravitational_constant() *
                        p->get_gravitational_potential());
            }
            LOGS("Finished calculating gravitational potential");
        }

        if(_gravity) {
            // calculate the gravitational forces
            LOGS("Start relative treewalk");
            _particles->get_tree().walk_tree<GravityWalker>(*_particles, true,
                                                            true, currentTime);
            // set old acceleration before multiplying with G
            for(unsigned int i = _particles->gassize(); i--;) {
                if(_particles->gas(i)->get_endtime() == currentTime) {
                    _particles->gas(i)->set_old_acceleration(
                            _grav_alpha *
                            _particles->gas(i)
                                    ->get_gravitational_acceleration()
                                    .norm());
                }
            }
            for(unsigned int i = _particles->dmsize(); i--;) {
                if(_particles->dm(i)->get_endtime() == currentTime) {
                    _particles->dm(i)->set_old_acceleration(
                            _grav_alpha *
                            _particles->dm(i)
                                    ->get_gravitational_acceleration()
                                    .norm());
                }
            }
            // multiply with G
            for(unsigned int i = 0; i < _particles->gassize(); i++) {
                _particles->gas(i)->set_gravitational_acceleration(
                        _physics->get_gravitational_constant() *
                        _particles->gas(i)->get_gravitational_acceleration());
                _particles->gas(i)->set_eta(
                        _physics->get_gravitational_constant() *
                        _particles->gas(i)->get_eta());
            }
            for(unsigned int i = 0; i < _particles->dmsize(); i++) {
                _particles->dm(i)->set_gravitational_acceleration(
                        _physics->get_gravitational_constant() *
                        _particles->dm(i)->get_gravitational_acceleration());
            }

            if(_particles->gassize()) {
                _voronoi->update_gravitational_corrections();
            }
            for(unsigned int i = _particles->gassize(); i--;) {
                if(_particles->gas(i)->get_endtime() == currentTime) {
#ifdef VARIABLE_SOFTENING
                    Vec acorr = _voronoi->get_gravitational_correction(i);
                    _particles->gas(i)->set_gravitational_acceleration(
                            _particles->gas(i)
                                    ->get_gravitational_acceleration() +
                            acorr);
#endif
                }
            }
            LOGS("Finished relative treewalk");
        }

        if(_gravity) {
            // update conserved variables with gravitational terms
            for(unsigned int i = _particles->gassize(); i--;) {
                if(_particles->gas(i)->get_endtime() == currentTime) {
                    _particles->gas(i)->apply_gravity(
                            0.5 *
                            _timeline->get_realtime_interval(
                                    _particles->gas(i)->get_endtime() -
                                    _particles->gas(i)->get_starttime()));
                }
            }

            // kick active collisionless particles
            for(unsigned int i = _particles->dmsize(); i--;) {
                if(_particles->dm(i)->get_endtime() == currentTime) {
                    _particles->dm(i)->accelerate(
                            _particles->dm(i)
                                    ->get_gravitational_acceleration() *
                            0.5 *
                            _timeline->get_realtime_interval(
                                    _particles->dm(i)->get_endtime() -
                                    _particles->dm(i)->get_starttime()));
                }
            }
            LOGS("Did gravitational kick");
        }
        _gravitytimer->stop();

        // The volume of the cells has changed, hence we have to recalculate the
        // Ws
        _hydrotimer->start();
        for(unsigned int i = 0; i < _particles->gassize(); i++) {
            // only active particles have received fluxes from all faces and
            // have a reliable value for their conserved quantities
            if(_particles->gas(i)->get_endtime() == currentTime) {
                double V = _voronoi->get_volume(i);
                double A = _voronoi->get_total_area(i);
#ifdef ENTROPY
                StateVector W = _solver->get_W(
                        V, _particles->gas(i)->get_Qvec(),
                        _particles->gas(i)->get_max_mach() == 0. ||
                                _particles->gas(i)->get_max_mach() > 1.1);
#else
                StateVector W =
                        _solver->get_W(V, _particles->gas(i)->get_Qvec());
#endif

                if(W.p() < 0. || W.p() != W.p()) {
                    cerr << "Error! Negative pressure after update of W!"
                         << endl;
                    cerr << W.p() << endl;
                    _particles->gas(i)->dump_ascii(cerr);
                    Vec v = _particles->gas(i)->get_velocity();
                    cerr << v[0] << "\t" << v[1];
#if ndim_ == 3
                    cerr << "\t" << v[2] << endl;
#else
                    cerr << endl;
#endif
                    Vec x = _particles->gas(i)->get_position();
                    cerr << x[0] << "\t" << x[1];
#if ndim_ == 3
                    cerr << "\t" << x[2] << endl;
#else
                    cerr << endl;
#endif
                    my_exit();
                }
                _particles->gas(i)->set_max_mach(0.);
                _particles->gas(i)->set_W(W);
                _particles->gas(i)->set_soundspeed(_solver->get_soundspeed(
                        _particles->gas(i)->get_Wvec()));
                _particles->gas(i)->set_centroid(_voronoi->get_centroid(i));
                _particles->gas(i)->set_total_area(A);
            }
        }
        _hydrotimer->stop();
        LOGS("Primitive variables calculated");

        if(_particles->gassize()) {
            _mpitimer->start();
            _voronoi->update_gradients();
            _mpitimer->stop();
            // no idea why there's a barrier here
            MyMPI_Barrier();
            // calculate new timestep and update the global time
            _voronoi->set_hs(currentTime);
        }

        _logfiles->write_to_logs(_timeline, _particles, _solver);

        check_for_restart();

        _timesteptimer->start();
        // calculate timestep also updates the tree and exchanges pseudonodes
        // we do not have to bother about informing the tree and other processes
        // about the changes in primitive variables
        dt = _timeline->calculate_timestep();
        _timesteptimer->stop();
    }

    // display goodbye message (in ASCII hardcoded green - because I can)
    cout << "\033[1;32mSimulation succesfully finished! Thank you for choosing "
            "Shadowfax!\033[0m"
         << endl;
}

/**
 * @brief Assign memory for the MPI communication buffer
 *
 * @param partsize Total number of particles in the simulation
 * @param maxsize Maximum allowed memory size of the buffer (in bytes)
 */
void Simulation::initialize_MPI_buffer(unsigned int partsize,
                                       unsigned int maxsize) {
    // make sure the buffer is large enough to fit at least one GasParticle and
    // an integer flag for every MPI process (except this one)
    if(maxsize < (sizeof(GasParticle) + sizeof(int)) * (MPIGlobal::size - 1)) {
        cerr << "Communication buffer too small for a single GasParticle!\n"
                "Try increasing the maximum size of the buffer!"
             << endl;
        my_exit();
    }

    unsigned int size = partsize * sizeof(GasParticle) * 100;
    // since size is a 32-bit integer, we have to make sure it is small enough
    // before entering the loop below. If not, we risk overflowing, in which
    // case we enter an endless loop...
    if(size > maxsize) {
        size = maxsize;
    }
    MPIGlobal::sendsize = 1;
    while(MPIGlobal::sendsize < size) {
        MPIGlobal::sendsize <<= 1;
    }
    MPIGlobal::sendsize = std::min(MPIGlobal::sendsize, maxsize);
    MPIGlobal::recvsize = MPIGlobal::sendsize;
    delete[] MPIGlobal::sendbuffer;
    delete[] MPIGlobal::recvbuffer;
    MPIGlobal::sendbuffer = new char[MPIGlobal::sendsize];
    MPIGlobal::recvbuffer = new char[MPIGlobal::recvsize];
    cout << "Assigned "
         << HelperFunctions::human_readable_bytes(MPIGlobal::sendsize)
         << " for MPI communication buffer" << endl;

    LOGS("Initialized MPI buffer");
}

/**
 * @brief Dump the simulation to the given RestartFile
 *
 * @param restartfile RestartFile to write to
 */
void Simulation::dump(RestartFile& restartfile) {
    LOGS("Writing RestartFile");

    _totaltimer.stop();
    restartfile.write(_periodic);
    restartfile.write(_gravity);
    if(_gravity) {
        restartfile.write(_grav_alpha);
    }
    restartfile.write(_outputdir);
    restartfile.write(_restarttime);

    _parameterfile->dump(restartfile);
    _logfiles->dump(restartfile);
    RiemannSolverFactory::dump(restartfile, _solver);
    _simulation_units->dump(restartfile);
    _output_units->dump(restartfile);

    _physics->dump(restartfile);

    _box->dump(restartfile);
    _particles->dump(restartfile);
    _timeline->dump(restartfile);
    // we only have a VorTessManager if there is gas
    if(_particles->gassize()) {
        _voronoi->dump(restartfile);
    }
    _totaltimer.dump(restartfile);
    _totaltimer.start();

    LOGS("RestartFile written");
}

/**
 * @brief Restart the simulation from the file with the given name
 *
 * @param filename Name of the restart file
 */
void Simulation::restart(string filename) {
    RestartFile restartfile(filename);

    restartfile.read(_periodic);
    if(_periodic) {
        cout << "Using periodic box\n" << endl;
    } else {
        cout << "Using reflective box\n" << endl;
    }
    restartfile.read(_gravity);
    if(_gravity) {
        restartfile.read(_grav_alpha);
    }
    restartfile.read(_outputdir);
    restartfile.read(_restarttime);

    _parameterfile = new ParameterFile(restartfile);

    // initialize logging (if enabled)
    LOGINIT(_outputdir);
    LOGS("Parameterfile reconstructed, logfile opened");

    _logfiles = new LogFiles(restartfile);
    _solver = RiemannSolverFactory::load(restartfile);
    _simulation_units = new UnitSet(restartfile);
    _physics = new Physics(*_simulation_units, _parameterfile);
    LOGS("Simulation UnitSet restarted");

    _output_units = new UnitSet(restartfile);
    LOGS("Output UnitSet restarted");

    _physics = new Physics(restartfile);

    _box = new RectangularBox(restartfile);
    _particles = new ParticleVector(restartfile, _parameterfile, *_box,
                                    _periodic, _periodic && _gravity);
    LOGS("ParticleVector restarted");
    _timeline = new TimeLine(restartfile, *_particles, *_simulation_units,
                             *_output_units);

    unsigned int max_memory = HelperFunctions::machine_readable_bytes(
            _parameterfile->get_parameter<string>(
                    "Memory.MaximumSize", SIMULATION_DEFAULT_MAXMEMORY));
    initialize_MPI_buffer(_particles->gassize() + _particles->dmsize(),
                          max_memory);
    // we only build a voronoi mesh if there is gas in the simulation
    // (obviously...)
    // we need to calculate a tree first, since we cannot read it in
    _particles->sort();
    LOGS("Particles resorted");

    if(_particles->gassize()) {
        _voronoi = new VorTessManager(restartfile, *_particles);
    }
    _totaltimer = Timer(restartfile);
    _totaltimer.start();

    LOGS("Simulation restarted");
}

/**
 * @brief Check if it is time to write a restart file and if so, write it
 *
 * @param force_restart If this parameter is true, we write a restart file
 * anyway, irrespective of the outcome of the test
 */
void Simulation::check_for_restart(bool force_restart) {
    LOGS("Checking for restart file...");
    // check if it is time to write a restart file
    // only one process should decide whether or not we write a restart file
    // otherwise, we risk asynchronous restart file writes and a huge mess
    int check = _restarttimer.interval() > _restarttime;
    MyMPI_Bcast(&check, 1, MPI_INT, 0);
    if(check || force_restart) {
        cout << "Writing restart file" << endl;
        RestartFile rfile(_outputdir, _timeline->get_realtime(
                                              _timeline->get_integertime()));
        dump(rfile);

        // restart the timer
        _restarttimer.restart();
    }
}

/**
 * @brief Initialize the simulation from the parameterfile with the given name
 *
 * @param filename Name of the parameter file
 */
void Simulation::initialize(string filename) {
    _parameterfile = new ParameterFile(filename);

    string outputdir = _parameterfile->get_parameter<string>(
            "Snapshots.OutputDir", SIMULATION_DEFAULT_OUTPUTDIR);
    _outputdir = HelperFunctions::get_absolute_path(outputdir);
    _restarttime = _parameterfile->get_parameter<double>(
            "Code.RestartTime", SIMULATION_DEFAULT_RESTARTTIME);

    // initialize global program log (if logging is enabled)
    LOGINIT(_outputdir);
    LOGS("Parameterfile created, programlog opened");

    _logfiles = new LogFiles(_outputdir);
    _solver = RiemannSolverFactory::generate(_parameterfile);
    string internal_units = _parameterfile->get_parameter<string>(
            "Units.InternalUnits", SIMULATION_DEFAULT_UNITS);
    _simulation_units = UnitSetGenerator::generate(internal_units);
    LOGS("Simulation UnitSet created");
    _physics = new Physics(*_simulation_units, _parameterfile);
    string output_units = _parameterfile->get_parameter<string>(
            "Units.OutputUnits", SIMULATION_DEFAULT_UNITS);
    _output_units = UnitSetGenerator::generate(output_units);
    LOGS("Output UnitSet created");
    _box = new RectangularBox();
    _particles = new ParticleVector(_parameterfile, *_box, _periodic, true);
    LOGS("ParticleVector created");
    _timeline = new TimeLine(_parameterfile, _outputdir, *_particles,
                             *_simulation_units, *_output_units, _periodic);
    string icfile_type = _parameterfile->get_parameter<string>(
            "IC.Type", SIMULATION_DEFAULT_ICTYPE);
    string icfile_name = _outputdir + string("/") +
                         _parameterfile->get_parameter<string>(
                                 "IC.FileName", SIMULATION_DEFAULT_ICNAME);
    double time =
            load(*_particles, icfile_type, icfile_name, *_simulation_units);
    _gravity = _parameterfile->check_parameter("Gravity.Gravity");
    // make sure the particles know if gravity is on or off
    _particles->get_header().set_gravity(_gravity);

    unsigned int max_memory = HelperFunctions::machine_readable_bytes(
            _parameterfile->get_parameter<string>(
                    "Memory.MaximumSize", SIMULATION_DEFAULT_MAXMEMORY));
    initialize_MPI_buffer(_particles->gassize() + _particles->dmsize(),
                          max_memory);

    _timeline->set_time(time);

    LOGS("Starting simulation");

    cout << "Starting simulation\n" << endl;

    _particles->sort();
    LOGS("Particles sorted");

    for(unsigned int i = _particles->gassize(); i--;) {
        _particles->gas(i)->set_starttime(_timeline->get_integertime());
    }
    for(unsigned int i = _particles->dmsize(); i--;) {
        _particles->dm(i)->set_starttime(_timeline->get_integertime());
    }

    if(_gravity) {
        double hsoft = _parameterfile->get_parameter<double>(
                "Gravity.Softening", SIMULATION_DEFAULT_SOFTENING);
        for(unsigned int i = _particles->gassize(); i--;) {
            _particles->gas(i)->set_hsoft(hsoft);
        }
        for(unsigned int i = _particles->dmsize(); i--;) {
            _particles->dm(i)->set_hsoft(2.8 * hsoft);
        }
        _particles->get_header().set_gravity(true);
        _particles->get_header().set_hsoft(hsoft);

        _grav_alpha = _parameterfile->get_parameter<double>(
                "Gravity.Alpha", SIMULATION_DEFAULT_GRAVALPHA);
    }
    LOGS("Particles initialized");

    // we only build a voronoi mesh if there is gas in the simulation
    // (obviously...)
    if(_particles->gassize()) {
        _voronoi = new VorTessManager(
                *_particles, _periodic,
                _parameterfile->get_parameter<double>(
                        "Voronoi.Tolerance",
                        SIMULATION_DEFAULT_VORONOITOLERANCE));
    }

    // convert W ic's to Q
    for(unsigned int i = _particles->gassize(); i--;) {
        double V = _voronoi->get_volume(i);
        double A = _voronoi->get_total_area(i);
        StateVector Q = _solver->get_Q(V, _particles->gas(i)->get_Wvec());
        _particles->gas(i)->set_Q(Q);
        _particles->gas(i)->set_soundspeed(
                _solver->get_soundspeed(_particles->gas(i)->get_Wvec()));
        _particles->gas(i)->set_centroid(_voronoi->get_centroid(i));
        _particles->gas(i)->set_total_area(A);
    }
    LOGS("Conserved variables initialized");
    if(_particles->gassize()) {
        _voronoi->update_Ws();
        _voronoi->set_hs(_timeline->get_integertime());
    }

#ifdef VARIABLE_SOFTENING
    if(_gravity) {
        // do this before recalculating the tree (since hmax has to be updated)
        for(unsigned int i = _particles->gassize(); i--;) {
            double hsoft = 0.75 * _voronoi->get_volume(i) / M_PI;
            hsoft = 1.5 * cbrt(hsoft);
            _particles->gas(i)->set_hsoft(hsoft);
        }
    }
#endif

    // this method will recalculate the masses and velocities for all nodes
    _particles->get_tree().set_velocities();

    if(_gravity) {
        // gravity:
        // we first calculate accelerations using a Barnes-Hut treewalk
        // we then use this estimate for the acceleration to walk the
        // tree again with a relative opening criterion (as in Gadget2)

        // Barnes-Hut
        LOGS("Starting Barnes-Hut gravity walk");
        _particles->get_tree().walk_tree<BHGravityWalker>(*_particles, true,
                                                          true, 0);

        // no multiplication by G needed, since it is already included in
        // Springel's relative criterion

        LOGS("Finished Barnes-Hut gravity walk");
        if(_particles->gassize()) {
            _voronoi->update_gravitational_corrections();
        }
        for(unsigned int i = _particles->gassize(); i--;) {
#ifdef VARIABLE_SOFTENING
            Vec acorr = _voronoi->get_gravitational_correction(i);
            _particles->gas(i)->set_gravitational_acceleration(
                    _particles->gas(i)->get_gravitational_acceleration() +
                    acorr);
#endif
            _particles->gas(i)->set_old_acceleration(
                    _grav_alpha *
                    _particles->gas(i)
                            ->get_gravitational_acceleration()
                            .norm());
        }
        for(unsigned int i = _particles->dmsize(); i--;) {
            _particles->dm(i)->set_old_acceleration(
                    _grav_alpha *
                    _particles->dm(i)->get_gravitational_acceleration().norm());
        }

        // relative criterion
        LOGS("Starting relative gravity walk");
        _particles->get_tree().walk_tree<GravityWalker>(*_particles, true, true,
                                                        0);

        // set old acceleration before multiplying by G
        for(unsigned int i = _particles->gassize(); i--;) {
            _particles->gas(i)->set_old_acceleration(
                    _grav_alpha *
                    _particles->gas(i)
                            ->get_gravitational_acceleration()
                            .norm());
        }
        for(unsigned int i = _particles->dmsize(); i--;) {
            _particles->dm(i)->set_old_acceleration(
                    _grav_alpha *
                    _particles->dm(i)->get_gravitational_acceleration().norm());
        }

        // multiply with G
        for(unsigned int i = 0; i < _particles->gassize(); i++) {
            _particles->gas(i)->set_gravitational_acceleration(
                    _physics->get_gravitational_constant() *
                    _particles->gas(i)->get_gravitational_acceleration());
            _particles->gas(i)->set_eta(_physics->get_gravitational_constant() *
                                        _particles->gas(i)->get_eta());
        }
        for(unsigned int i = 0; i < _particles->dmsize(); i++) {
            _particles->dm(i)->set_gravitational_acceleration(
                    _physics->get_gravitational_constant() *
                    _particles->dm(i)->get_gravitational_acceleration());
        }
        LOGS("Finished relative gravity walk");
        if(_particles->gassize()) {
            _voronoi->update_gravitational_corrections();
        }
        for(unsigned int i = _particles->gassize(); i--;) {
#ifdef VARIABLE_SOFTENING
            Vec acorr = _voronoi->get_gravitational_correction(i);
            _particles->gas(i)->set_gravitational_acceleration(
                    _particles->gas(i)->get_gravitational_acceleration() +
                    acorr);
#endif
        }
    }

    for(unsigned int i = _particles->gassize(); i--;) {
#ifndef STATIC
        Vec v = _voronoi->get_velocity(i);
        _particles->gas(i)->set_velocity(v);
#endif
    }
    if(_particles->gassize()) {
        // this communicates gradient information between processes
        _voronoi->update_gradients();
    }
    _particles->get_tree().set_velocities();
    if(MPIGlobal::size > 1) {
        _particles->get_tree().exchange_pseudonodes();
    }
    cout << "Initialization ready" << endl;

    LOGS("Initialization ready");

    check_for_restart(true);
}

/**
 * @brief Print a welcome message to the stdout, containing a fancy wordart
 * "Shadowfax"
 */
void Simulation::print_welcome_message() {
    // show welcome message (made using
    // http://patorjk.com/software/taag/#p=display&h=3&f=Epic&t=Shadowfax,
    // converted to C++ using
    // http://tomeko.net/online_tools/cpp_text_escape.php?lang=en)
    // oh, yeah, it's also scrambled a bit to make it fit in a 70 column window
    cout << " _______         _______ ______  _______         _______ _______  "
            "       \n(  ____ |\\     /(  ___  (  __  \\(  ___  |\\     /(  ___"
            "_ (  ___  |\\     /|\n| (    \\| )   ( | (   ) | (  \\  | (   ) | "
            ")   ( | (    \\| (   ) ( \\   / )\n| (_____| (___) | (___) | |   )"
            " | |   | | | _ | | (__   | (___) |\\ (_) / \n(_____  |  ___  |  __"
            "_  | |   | | |   | | |( )| |  __)  |  ___  | ) _ (  \n      ) | ( "
            "  ) | (   ) | |   ) | |   | | || || | (     | (   ) |/ ( ) \\ \n/"
            "\\____) | )   ( | )   ( | (__/  | (___) | () () | )     | )   ( ( "
            "/   \\ )\n\\_______|/     \\|/     \\(______/(_______(_______|/   "
            "   |/     \\|/     \\|\n"
         << endl;
    cout << "This is Shadowfax, a parallel moving mesh hydrodynamical code\n"
            "written by Bert Vandenbroucke (bert.vandenbroucke@ugent.be)\n"
         << endl;
}
