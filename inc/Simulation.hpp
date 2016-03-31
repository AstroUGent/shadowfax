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
 * @file Simulation.hpp
 *
 * @brief Main simulation program: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <vector>
#include <string>
#include <utilities/Timer.hpp>

class ParameterFile;
class RiemannSolver;
class UnitSet;
class RectangularBox;
class ParticleVector;
class TimeLine;
class VorTessManager;
class RestartFile;
class LogFiles;
class Physics;

/**
  * \brief The main class for actual simulations
  *
  * A simulation object performs an actual simulation. A parameterfile is
  * extracted from the command line arguments or set to "default.ini" if not
  * provided. The constructor then sets up the required objects to run the
  * simulation and starts the time integration.
  */
class Simulation{
private:
    /*! \brief Flag indicating if the simulation box is periodic (true) or
     *  reflective (false) */
    bool _periodic;

    /*! \brief ParameterFile containing simulation parameters */
    ParameterFile* _parameterfile;

    /*! \brief User-friendly log files */
    LogFiles *_logfiles;

    /*! \brief RiemannSolver used to solve the Riemann problem for the
     *  hydrodynamics */
    RiemannSolver* _solver;

    /*! \brief Internal simulation UnitSet */
    UnitSet* _simulation_units;
    /*! \brief Output UnitSet */
    UnitSet* _output_units;

    /*! \brief RectangularBox specifying the dimensions of the simulation box */
    RectangularBox* _box;

    /*! \brief ParticleVector containing the particles of the simulation */
    ParticleVector* _particles;

    /*! \brief TimeLine governing the timeline of the simulation */
    TimeLine* _timeline;

    /*! \brief Physical constants used during the simulation */
    Physics *_physics;

    /*! \brief VorTessManager used to calculate the Voronoi grid used for the
     *  hydrodynamical integration */
    VorTessManager* _voronoi;

    /*! \brief Timer used to time the entire simulation */
    Timer _totaltimer;
    /*! \brief Timer used to time the initialization of the simulation */
    Timer* _starttimer;
    /*! \brief Timer used to quantify the time spent in calculating particle
     *  timesteps */
    Timer* _timesteptimer;
    /*! \brief Timer used to quantify the time spent in MPI operations */
    Timer* _mpitimer;
    /*! \brief Timer used to quantify the time spent in parts of the program
     *  related to the hydrodynamics */
    Timer* _hydrotimer;

    /*! \brief Timer used to quantify the time spent in gravity related
     *  operations */
    Timer* _gravitytimer;

    /*! \brief Timer used to determine when it is time to do restart dumps */
    Timer _restarttimer;

    void main_loop();

    void print_welcome_message();
    void initialize(std::string filename);

    double load(ParticleVector& cells, std::string ictype, std::string filename,
                UnitSet& simulation_units);

    void initialize_MPI_buffer(unsigned int partsize, unsigned int maxsize);
    void dump(RestartFile &restartfile);
    void restart(std::string filename);

    void check_for_restart(bool force_restart = false);

public:
    Simulation(int argc, char** argv);

    ~Simulation();
};

#endif // SIMULATION_HPP
