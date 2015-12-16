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
 * @file ParameterFile.hpp
 *
 * @brief Parameter file: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PARAMETERFILE_HPP
#define PARAMETERFILE_HPP

#include <string>

class RestartFile;

/**
 * @brief Abstraction of the parameter file that contains vital run information
 *
 * The parameter file is a .ini-file which contains parameters that set the code
 * behaviour at runtime. Every parameter has a default value which is used if
 * the parameter is not provided.
 *
 * A complete parameter file would look like this
\verbatim
[Time]
; Total simulation time (no default value, has to be set)
MaxTime = 1.0
; Use a global timestep? (default: false)
GlobalTimestep = false
; Maximum size of the particle timestep (default: MaxTime)
MaxTimeStep = 1.0
; Minimum size of the particle timestep (default: MaxTime/(2**60))
MinTimeStep = 0.867361738e-18
; Use the expensive treewalk timestep criterion? (default: false)
TreeTime = false

[Snapshots]
; Prefix for snapshot name. The snapshot is named $(BaseName)xxx.hdf5
; (default: snapshot)
BaseName = snapshot
; Time between subsequent snapshots (default: 10% of MaxTime)
SnapTime = 0.1
; Index number of first snapshot that is written out (default: 0) Adjust this
; when restarting from a snapshot
FirstSnap = 0
; Directory in which snapshots, log-files and restart files will be stored,
; should not contain trailing '/' (default: .)
OutputDir = .
; Type of the snapshot files (Shadowfax/Gadget, default: Gadget)
Type = Gadget
; Indicate whether every node writes a separate snapshot file, or all nodes
; write (or try to write) to the same file (default: false)
PerNodeOutput = false

[IC]
; Name of the IC-file (default: $(BaseName)$(FirstSnap).hdf5)
FileName = snapshot000.hdf5
; Type of the IC-file (Shadowfax/Gadget, default: Gadget)
Type = Gadget

[RiemannSolver]
; Type of Riemann solver used to solve the Riemann problem (Exact/TRRS, default:
; Exact)
Type = Exact
; Tolerance used for the Newton-Raphson iteration (default: 1.0e-8)
Tolerance = 1.0e-8
; CutOff used to distinguish between Newton-Rapshon regime and Brent regime in
; solution (default: -5.0)
CutOff = -5.0
; Courant-Friedrichs-Lewy parameter for timestep limitation (default: 0.4)
CFL = 0.4

[Hydro]
; Adiabatic index for the polytropic equation of state of the gas
; (default: 1.66667 =~= 5/3)
Gamma = 1.66667

[Gravity]
; Use gravity? (default: true)
Gravity = true
; Softening length used for particles with a fixed softening length, in internal
; code units (default: 0.03)
Softening = 0.03
; Eta factor used in gravitational time step criterion (default: 0.05/2.8)
Eta = 0.017857143
; Alpha factor used in relative tree opening criterion (default: 0.005)
Alpha = 0.005

[Voronoi]
; Tolerance used in the computation of the Delaunay tesselation to distinguish
; between exact and inexact geometric tests (default: 1.e-9)
Tolerance = 1.0e-9

[Tree]
; Side of the grid on which the Ewald correction to the gravitational force or
; mesh movement is precomputed (default: 64) Warning: this scales as 3*Size^3 in
; memory!
EwaldSize = 64
; Parameter alpha determining the cutoff between short range and long range in
; Ewald's method (default: 2.0)
EwaldAlpha = 2.0

[Memory]
; Maximum size of the MPI-buffer in memory (default: 1 GB)
MaximumSize = 1 GB

[Code]
; Computational time interval after which to write a restart-file, in seconds
; (default: 3600.0)
RestartTime = 3600.0

[Units]
; Units used during the run (SI/CGS/galactic, default: SI)
InternalUnits = SI
; Units used for output to snapshots (SI/CGS/galactic, default: SI)
OutputUnits = SI

[Physics]
; Use physical values for physical constants? (default: true)
RealPhysics = true
\endverbatim
 */
class ParameterFile{
private:
    // time options
    /*! \brief End time of the simulation in internal units */
    double _maxtime;
    /*! \brief Flag indicating if we use a global or individual timestep */
    bool _global_timestep;
    /*! \brief Maximal value of the Particle timestep in internal units */
    double _max_timestep;
    /*! \brief Minimal value of the Particle timestep in internal units */
    double _min_timestep;
    /*! \brief Flag indicating if we use the expensive treewalk timestep
     *  criterion */
    bool _treetime;

    // snapshot options
    /*! \brief Basic name of the snapshot dumps */
    std::string _basename;
    /*! \brief Time interval in between successive snapshots in internal
     *  units */
    double _snaptime;
    /*! \brief Index of the first snapshot file that should be written out */
    unsigned int _firstsnap;
    /*! \brief Output directory where the snapshots and other diagnostic files
     *  will be stored */
    std::string _outputdir;
    /*! \brief Type of the snapshot files */
    std::string _snapshot_type;
    /*! \brief Flag indicating if each node writes a separate snapshot file or
     *  all files (try to) write to the same file */
    bool _per_node_output;

    // IC options
    /*! \brief Name of the initial condition file */
    std::string _icfile;
    /*! \brief Type of the initial condition file */
    std::string _icfile_type;

    // Riemann solver options
    /*! \brief Type of the Riemann solver used to solve the Riemann problem */
    std::string _solver_type;
    /*! \brief Tolerance value used in the Riemann solver iteration */
    double _tolerance;
    /*! \brief Cut off used in the exact Riemann solver */
    double _cutoff;
    /*! \brief Courant-Friedrich-Levy constant for time hydrodynamical time
     *  stepping */
    double _CFL;

    // Hydro options
    /*! \brief Adiabatic index of the ideal gas */
    double _gamma;

    // Gravity options
    /*! \brief Flag indicating if we use gravity */
    bool _gravity;
    /*! \brief Softening length used for gravitational force softening when no
     *  individual softening lengths are provided */
    double _hsoft;
    /*! \brief Gravitational \f$\eta\f$ factor */
    double _grav_eta;
    /*! \brief Gravitational \f$\alpha\f$ factor */
    double _grav_alpha;

    // Voronoi options
    /*! \brief Tolerance used to distinguish between exact and approximate
     *  geometric tests during Voronoi grid construction */
    double _voronoi_tolerance;

    // Tree options
    /*! \brief Ewald \f$\alpha\f$ factor */
    double _ewald_alpha;
    /*! \brief Size of the Ewald table */
    unsigned int _ewald_size;

    // Memory options
    /*! \brief Maximum size of the MPI buffer in memory (in bytes) */
    unsigned int _max_memory;

    // Code options
    /*! \brief Time interval (in CPU time) between successive restart file
     *  dumps */
    double _restarttime;

    // Units
    /*! \brief Internal units */
    std::string _units_internal;
    /*! \brief Units used for output files */
    std::string _units_output;

    // Physics
    /*! \brief Physical or idealized values for physical constants? */
    bool _real_physics;

    void print_contents();

public:
    ParameterFile(std::string name);
    ~ParameterFile(){}

    // time options
    double get_maxtime();
    bool get_global_timestep();
    double get_max_timestep();
    double get_min_timestep();
    bool has_treetime();

    // snapshot options
    std::string get_basename();
    double get_snaptime();
    unsigned int get_firstsnap();
    std::string get_outputdir();
    std::string get_snapshot_type();
    bool get_per_node_output();

    // IC options
    std::string get_icfile();
    std::string get_icfile_type();

    // Riemann solver options
    std::string get_solver_type();
    double get_tolerance();
    double get_cutoff();
    double get_CFL();

    // Hydro options
    double get_gamma();

    // Gravity options
    bool has_gravity();
    double get_hsoft();
    double get_grav_eta();
    double get_grav_alpha();

    // Voronoi options
    double get_voronoi_tolerance();

    // Tree options
    double get_ewald_alpha();
    unsigned int get_ewald_size();

    // Memory options
    unsigned int get_max_memory();

    // Code options
    double get_restarttime();

    // Units
    std::string get_units_internal();
    std::string get_units_output();

    // Physics options
    bool has_real_physics();

    void dump(RestartFile &rfile);
    ParameterFile(RestartFile &rfile);
};

#endif // PARAMETERFILE_HPP
