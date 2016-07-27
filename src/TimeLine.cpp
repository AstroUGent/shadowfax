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
 * @file TimeLine.cpp
 *
 * @brief Simulation timeline: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "TimeLine.hpp"
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include "ParameterFile.hpp"
#include "ProgramLog.hpp"
#include "RestartFile.hpp"
#include "SnapshotHandler.hpp"
#include "SnapshotWriterFactory.hpp"
#include "TimeStepWalker.hpp"
#include "utilities/DMParticle.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/ParticleVector.hpp"
using namespace std;

/**
 * @brief Construct a timeline with a given endtime for the particles in the
 * given ParticleVector.
 *
 * Every snaptime a snapshot with snapname is written out, starting from
 * lastsnap (with default value 0).
 *
 * @param maxtime The endtime of the simulation
 * @param snaptime The time between two subsequent snapshots
 * @param cfl Courant-Friedrich-Levy parameter used for hydrodynamical timestep
 * stability
 * @param grav_eta Parameter setting the size of the gravitational timestep
 * @param particlevector Reference to the simulation ParticleVector
 * @param snaptype Type of the snapshot files
 * @param snapname The name of the snapshots (a three digit index is attached
 * to the end of this name)
 * @param units Internal simulation UnitSet
 * @param output_units Output UnitSet
 * @param gravity Bool specifying if a gravitational timestep criterion should
 * be used
 * @param periodic Bool specifying whether the box containing the
 * ParticleVector is periodic (true) or reflective (false)
 * @param treetime Bool specifying whether the timesteps should be calculated
 * using the expensive treewalk
 * @param max_timestep Maximum particle timestep allowed
 * @param min_timestep Minimum particle timestep allowed
 * @param lastsnap The starting index of the snapshots (default value: 0)
 * @param per_node_output Flag indicating if each node should write a separate
 * snapshot file or all nodes should write to the same file (if possible)
 */
TimeLine::TimeLine(double maxtime, double snaptime, double cfl, double grav_eta,
                   ParticleVector& particlevector, std::string snaptype,
                   std::string snapname, UnitSet& units, UnitSet& output_units,
                   bool gravity, bool periodic, bool treetime,
                   double max_timestep, double min_timestep,
                   unsigned int lastsnap, bool per_node_output)
        : _particles(particlevector) {
    _snapshotwriter = SnapshotWriterFactory::generate(
            snaptype, snapname, units, output_units, lastsnap, per_node_output);
    _maxtime = maxtime;
    _integer_maxtime = 1;
    _integer_maxtime <<= 60;
    _current_time = 0;
    _timestep = 0;
    _snaptime = snaptime;
    _cfl = cfl;
    _gravity = gravity;
    _periodic = periodic;
    _treetime = treetime;
    _grav_eta = grav_eta;
    _max_timestep = _integer_maxtime;
    if(max_timestep) {
        max_timestep /= maxtime;
        while(_max_timestep > max_timestep * _integer_maxtime) {
            _max_timestep >>= 1;
        }
    }
    if(min_timestep) {
        min_timestep /= maxtime;
        _min_timestep = _integer_maxtime;
        while(_min_timestep > max_timestep * _integer_maxtime) {
            _min_timestep >>= 1;
        }
    } else {
        // a timestep of zero seems like a bad idea (although a timestep of 1 is
        // also ridiculously small)
        _min_timestep = 1;
    }

    _iotimer = new Timer();

    LOGS("TimeLine created");
}

/**
 * @brief ParameterFile constructor
 *
 * @param parameters ParameterFile containing the desired parameter values
 * @param outputdir Output directory (already converted to absolute path)
 * @param particlevector Reference to the ParticleVector holding all particles
 * @param units Internal units
 * @param output_units Output units
 * @param periodic Flag indicating if the simulation box is periodic
 * @param cosmology Cosmology used for comoving integration
 */
TimeLine::TimeLine(ParameterFile* parameters, std::string outputdir,
                   ParticleVector& particlevector, UnitSet& units,
                   UnitSet& output_units, bool periodic)
        : TimeLine(parameters->get_parameter<double>("Time.MaxTime", -1000.),
                   parameters->get_parameter<double>(
                           "Snapshots.SnapTime",
                           0.1 *
                                   parameters->get_parameter<double>(
                                           "Time.MaxTime", -1000.)),
                   parameters->get_parameter<double>("RiemannSolver.CFL",
                                                     TIMELINE_DEFAULT_CFL),
                   parameters->get_parameter<double>("Gravity.Eta",
                                                     TIMELINE_DEFAULT_GRAVETA),
                   particlevector,
                   parameters->get_parameter<string>("Snapshots.Type",
                                                     TIMELINE_DEFAULT_SNAPTYPE),
                   outputdir + string("/") +
                           parameters->get_parameter<string>(
                                   "Snapshots.BaseName",
                                   TIMELINE_DEFAULT_BASENAME),
                   units, output_units,
                   parameters->check_parameter("Gravity.Gravity"), periodic,
                   parameters->check_parameter("Time.TreeTime"),
                   parameters->get_parameter<double>(
                           "Time.MaxTimeStep", TIMELINE_DEFAULT_MAXTIMESTEP),
                   parameters->get_parameter<double>(
                           "Time.MinTimeStep", TIMELINE_DEFAULT_MINTIMESTEP),
                   parameters->get_parameter<unsigned int>(
                           "Snapshots.FirstSnap", TIMELINE_DEFAULT_LASTSNAP),
                   parameters->check_parameter("Snapshots.PerNodeOutput")) {}

/**
  * @brief Get the physical current time on the simulation timeline
  *
  * @returns The physical (floating point) time of the simulation
  */
double TimeLine::get_time() {
    return _current_time * _maxtime / (double)_integer_maxtime;
}

/**
  * @brief Get the internal (integer) time of the TimeLine
  *
  * @returns A 64-bit representing the current internal time
  */
unsigned long TimeLine::get_integertime() {
    return _current_time;
}

/**
  * @brief Increase the internal time by the internal timestep and check if a
  * snapshot has to be written out
  *
  * Check if we have reached the end of the simulation.
  *
  * @returns False if we have reached the end of the simulation, true otherwise
  */
bool TimeLine::step_forward() {
    _current_time += _timestep;
    if(get_time() >= _snaptime * _snapshotwriter->get_lastsnap()) {
        _iotimer->start();
        _snapshotwriter->write_snapshot(get_time(), _particles);
        _iotimer->stop();
    }
    LOGS("TimeLine step_forward");

    return _current_time <= _integer_maxtime;
}

/**
  * @brief Set the current internal time of the TimeLine
  *
  * @param time A new 64-bit integer current time for the TimeLine
  */
void TimeLine::set_time(unsigned long time) {
    _current_time = time;
}

/**
 * @brief Set the current internal time of the timeline
 *
 * The floating point time is converted to an integer time on the internal
 * integer timeline.
 *
 * @param time Floating point time
 */
void TimeLine::set_time(double time) {
    _current_time = (unsigned long)(time * _integer_maxtime / _maxtime);
}

/**
 * @brief Calculate new timesteps for all active particles and set the system
 * step to the timestep of the smallest particle
 *
 * @return Integer system timestep for the next system step
 */
unsigned long TimeLine::calculate_timestep() {
    LOGS("Starting timestep calculation");

    unsigned int numactive = 0;

    if(_treetime) {
        _particles.get_tree().set_velocities();
        if(MPIGlobal::size > 1) {
            _particles.get_tree().exchange_pseudonodes();
        }
        // determine hydro timesteps
        _particles.get_tree().walk_tree<TimeStepWalker>(
                _particles, true, false, _current_time + _timestep);
    } else {
        for(unsigned int i = 0; i < _particles.gassize(); i++) {
            if(_particles.gas(i)->get_endtime() == _current_time + _timestep) {
                GasParticle* p = _particles.gas(i);
                StateVector W = p->get_Wvec();
                Vec v;
#if ndim_ == 3
                v.set(W[1], W[2], W[3]);
#else
                v.set(W[1], W[2]);
#endif
                v -= p->get_velocity();
                double vi = v.norm();
                double ci = p->get_soundspeed();
                double dt = p->h() / (vi + ci);
                _particles.gas(i)->set_real_timestep(dt);
            }
        }
    }

    unsigned long Smax = _integer_maxtime;
    unsigned long Smin = Smax;
    double eta = _grav_eta;
    // put them in a power of two hierarchy
    for(unsigned int i = 0; i < _particles.gassize(); i++) {
        if(_particles.gas(i)->get_endtime() == _current_time + _timestep) {
            numactive++;
            double S = _cfl * _particles.gas(i)->get_real_timestep();
            if(_gravity) {
                double a = _particles.gas(i)
                                   ->get_gravitational_acceleration()
                                   .norm();
                // apply gravitational timestep criterion
                double eps = _particles.gas(i)->get_hsoft();
                if(a) {
                    S = std::min(S, sqrt(2. * eps * eta / a));
                }
            }
            unsigned long S2;
            if(get_realtime(Smax) < S) {
                S2 = Smax;
            } else {
                S2 = _integer_maxtime;
                while(get_realtime(S2) > S) {
                    S2 >>= 1;
                }

                if(!S2) {
                    cerr << "Timestep 0!" << endl;
                    my_exit();
                }

                while((_integer_maxtime - _particles.gas(i)->get_endtime()) %
                      S2) {
                    S2 >>= 1;  // cells[i]->get_timestep();
                }
                S2 = std::max(S2, _min_timestep);
                S2 = std::min(S2, _max_timestep);
                Smin = std::min(S2, Smin);
                if(!_particles.global_timestep()) {
                    _particles.gas(i)->set_timestep(S2);
                }
            }
        }
    }

    for(unsigned int i = _particles.dmsize(); i--;) {
        if(_particles.dm(i)->get_endtime() != _current_time + _timestep) {
            continue;
        }
        numactive++;
        double a = _particles.dm(i)->get_gravitational_acceleration().norm();
        // if the acceleration is 0, timestep is formally infinite, so we set it
        // to the maximal timestep
        unsigned long tidt = Smax;
        if(a) {
            double eps = _particles.dm(i)->get_hsoft();
            double dt = sqrt(2. * eps * eta / a);
            tidt = std::min(
                    tidt, (unsigned long)((dt / _maxtime) * _integer_maxtime));
        }
        // make the timestep a power of 2 subdivision of the total simulation
        // time
        unsigned long pidt = _integer_maxtime;
        while(tidt < pidt) {
            pidt >>= 1;
        }

        if(!pidt) {
            cerr << "Timestep 0!" << endl;
            my_exit();
        }

        // if the timestep wants to increase, make sure it synchronizes
        while((_integer_maxtime - _particles.dm(i)->get_endtime()) % pidt) {
            pidt >>= 1;
        }
        pidt = std::max(pidt, _min_timestep);
        pidt = std::min(pidt, _max_timestep);
        if(!_particles.global_timestep()) {
            // increase particle timestep
            _particles.dm(i)->set_starttime(_particles.dm(i)->get_endtime());
            _particles.dm(i)->set_endtime(_particles.dm(i)->get_endtime() +
                                          pidt);
        }
        // set global timestep to smallest timestep of all particles
        Smin = std::min(pidt, Smin);
    }

    _particles.set_numactive(numactive);

    unsigned long Smin_glob;
    MyMPI_Allreduce(&Smin, &Smin_glob, 1, MPI_LONG, MPI_MIN);
    if(_particles.global_timestep()) {
        for(unsigned int i = _particles.gassize(); i--;) {
            _particles.gas(i)->set_timestep(Smin_glob);
        }
        for(unsigned int i = _particles.dmsize(); i--;) {
            _particles.dm(i)->set_timestep(Smin_glob);
        }
    }

    LOGS("Timestep calculation finished");
    return Smin_glob;
}

/**
 * @brief Calculate a new gravitational timestep for all active particles and
 * set the system step to the smallest timestep of all particles
 *
 * @return Integer timestep to be used for the next system step
 */
unsigned long TimeLine::calculate_gravitational_timestep() {
    double eta = _grav_eta;
    unsigned long Smax = _integer_maxtime;
    unsigned long Smin = Smax;

    for(unsigned int i = _particles.gassize(); i--;) {
        if(_particles.gas(i)->get_endtime() != _current_time + _timestep) {
            continue;
        }
        double a = _particles.gas(i)->get_gravitational_acceleration().norm();
        // if the acceleration is 0, timestep is formally infinite, so we set it
        // to the maximal timestep
        unsigned long tidt = Smax;
        if(a) {
            double eps = _particles.gas(i)->get_hsoft();
            double dt = sqrt(2. * eps * eta / a);
            tidt = std::min(
                    tidt, (unsigned long)((dt / _maxtime) * _integer_maxtime));
        }
        // make the timestep a power of 2 subdivision of the total simulation
        // time
        unsigned long pidt = _integer_maxtime;
        while(tidt < pidt) {
            pidt >>= 1;
        }
        // if the timestep wants to increase, make sure it synchronizes
        while((_integer_maxtime - _particles.gas(i)->get_endtime()) % pidt) {
            pidt >>= 1;
        }
        if(!_particles.global_timestep()) {
            // increase particle timestep
            _particles.gas(i)->set_starttime(_particles.gas(i)->get_endtime());
            _particles.gas(i)->set_endtime(_particles.gas(i)->get_endtime() +
                                           pidt);
        }
        // set global timestep to smallest timestep of all particles
        Smin = std::min(pidt, Smin);
    }

    unsigned long Smin_glob;
    MyMPI_Allreduce(&Smin, &Smin_glob, 1, MPI_LONG, MPI_MIN);
    if(_particles.global_timestep()) {
        for(unsigned int i = _particles.gassize(); i--;) {
            _particles.gas(i)->set_timestep(Smin_glob);
        }
    }
    return Smin_glob;
}

/**
 * @brief Convert the given integer time to physical time
 *
 * @param integer_time Integer time on the integer timeline
 * @return Floating point physical time
 */
double TimeLine::get_realtime(unsigned long integer_time) {
    double t = (double)integer_time / (double)_integer_maxtime;
    return _maxtime * t;
}

/**
 * @brief Set the system timestep
 *
 * @param timestep New value for the system timestep
 */
void TimeLine::set_timestep(unsigned long timestep) {
    _timestep = timestep;
}

/**
 * @brief Get the system timestep
 *
 * @return The integer system timestep
 */
unsigned long TimeLine::get_timestep() {
    return _timestep;
}

/**
 * @brief Check if gravity is switched on
 *
 * @return True if gravity is on, false otherwise
 */
bool TimeLine::has_gravity() {
    return _gravity;
}

/**
 * @brief Dump the timeline to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void TimeLine::dump(RestartFile& rfile) {
    SnapshotWriterFactory::dump(rfile, _snapshotwriter);
    rfile.write(_maxtime);
    rfile.write(_integer_maxtime);
    rfile.write(_current_time);
    rfile.write(_timestep);
    rfile.write(_snaptime);
    rfile.write(_cfl);
    rfile.write(_gravity);
    rfile.write(_periodic);
    rfile.write(_grav_eta);
    rfile.write(_max_timestep);
    rfile.write(_min_timestep);
    rfile.write(_treetime);

    _iotimer->dump(rfile);

    LOGS("TimeLine dumped");
}

/**
 * @brief Restart constructor. Initialize the timeline from the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 * @param particlevector Reference to the ParticleVector of the simulation
 * @param units Internal simulation UnitSet
 * @param output_units Output UnitSet
 */
TimeLine::TimeLine(RestartFile& rfile, ParticleVector& particlevector,
                   UnitSet& units, UnitSet& output_units)
        : _particles(particlevector) {
    _snapshotwriter = SnapshotWriterFactory::load(rfile, units, output_units);
    rfile.read(_maxtime);
    rfile.read(_integer_maxtime);
    rfile.read(_current_time);
    rfile.read(_timestep);
    rfile.read(_snaptime);
    rfile.read(_cfl);
    rfile.read(_gravity);
    rfile.read(_periodic);
    rfile.read(_grav_eta);
    rfile.read(_max_timestep);
    rfile.read(_min_timestep);
    rfile.read(_treetime);

    _iotimer = new Timer(rfile);

    LOGS("TimeLine restarted");
}
