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
 * @file TimeLine.hpp
 *
 * @brief Simulation timeline: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TIMELINE_HPP
#define TIMELINE_HPP

#include "GadgetSnapshotWriter.hpp"
#include "ShadowfaxSnapshotWriter.hpp"
#include "utilities/Timer.hpp"
#include <iostream>
#include <string>

class ParticleVector;
class ParameterFile;
class UnitSet;
class RestartFile;

#define TIMELINE_DEFAULT_CFL 0.4
#define TIMELINE_DEFAULT_GRAVETA 0.017857143
#define TIMELINE_DEFAULT_SNAPTYPE "Gadget"
#define TIMELINE_DEFAULT_BASENAME "snapshot"
#define TIMELINE_DEFAULT_TREETIMEFLAG false
#define TIMELINE_DEFAULT_MAXTIMESTEP 0.
#define TIMELINE_DEFAULT_MINTIMESTEP 0.
#define TIMELINE_DEFAULT_LASTSNAP 0
#define TIMELINE_DEFAULT_PERNODEOUTPUTFLAG false

/**
  * \brief Class governing the simulation timeline.
  *
  * The TimeLine keeps track of the current time of the simulation.
  * The simulation is advanced in time by subsequently calling
  * TimeLine::step_forward(), until false is returned.
  *
  * The TimeLine is also responsible for detecting when a snapshot has to be
  * written and calls a SnapshotWriter to do so.
  */
class TimeLine {
  private:
    /*! \brief End time of the simulation */
    double _maxtime;
    /*! \brief End time on the integer timeline */
    unsigned long _integer_maxtime;
    /*! \brief Current time on the integer timeline */
    unsigned long _current_time;
    /*! \brief Current system timestep on the integer timeline */
    unsigned long _timestep;
    /*! \brief Time interval in between successive snapshots */
    double _snaptime;

    /*! \brief Reference to the simulation ParticleVector */
    ParticleVector& _particles;
    /*! \brief SnapshotWriter used to write snapshots */
    SnapshotWriter* _snapshotwriter;
    /*! \brief Courant-Friedrich-Levy constant for hydrodynamical timestep
     *  calculation */
    double _cfl;

    /*! \brief Flag indicating if gravity is switched on */
    bool _gravity;
    /*! \brief Flag indicating if the simulation box is periodic (true) or
     *  reflective (false) */
    bool _periodic;

    /*! \brief Gravitational timestep constant */
    double _grav_eta;

    /*! \brief Maximal allowed timestep on the integer timeline */
    unsigned long _max_timestep;
    /*! \brief Minimal allowed timestep on the integer timeline */
    unsigned long _min_timestep;

    /*! \brief Flag indicating whether the timestep should be calculated using
     *  the expensive treewalk criterion */
    bool _treetime;

    /*! \brief Timer to quantify time spent in writing snapshot files */
    Timer* _iotimer;

  public:
    TimeLine(double maxtime, double snaptime, double cfl, double grav_eta,
             ParticleVector& particlevector, std::string snaptype,
             std::string snapname, UnitSet& units, UnitSet& output_units,
             bool gravity, bool periodic,
             bool treetime = TIMELINE_DEFAULT_TREETIMEFLAG,
             double max_timestep = TIMELINE_DEFAULT_MAXTIMESTEP,
             double min_timestep = TIMELINE_DEFAULT_MINTIMESTEP,
             unsigned int lastsnap = TIMELINE_DEFAULT_LASTSNAP,
             bool per_node_output = TIMELINE_DEFAULT_PERNODEOUTPUTFLAG);
    TimeLine(ParameterFile* parameters, std::string outputdir,
             ParticleVector& particlevector, UnitSet& units,
             UnitSet& output_units, bool periodic);

    /**
     * @brief Destructor
     *
     * Clean up the snapshot writer.
     */
    ~TimeLine() {
        std::cout << "Spent " << _iotimer->value() << "s writing snapshots"
                  << std::endl;
        delete _iotimer;
        delete _snapshotwriter;
    }

    double get_time();
    unsigned long get_integertime();
    bool step_forward();
    void set_time(unsigned long time);
    void set_time(double time);

    double get_realtime(unsigned long integer_time);

    void set_timestep(unsigned long timestep);
    unsigned long get_timestep();

    unsigned long calculate_timestep();
    unsigned long calculate_gravitational_timestep();

    // temporary fix to access gravity in VorFace
    bool has_gravity();

    void dump(RestartFile& rfile);
    TimeLine(RestartFile& rfile, ParticleVector& particlevector, UnitSet& units,
             UnitSet& output_units);
};

#endif  // TIMELINE_HPP
