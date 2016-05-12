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
 * @file ParameterFile.cpp
 *
 * @brief Parameter file: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ParameterFile.hpp"
#include "Error.hpp"
#include "RestartFile.hpp"
#include "utilities/HelperFunctions.hpp"
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

// if you want to add a parameter:
//  * make the variable in ParameterFile.hpp
//  * make and implement a getter for the variable
//  * read it in in the constructor
//  * make an entry in print_contents
//  * read and write it in dump and the restart constructor

/**
 * @brief Print the contents of the parameter file in human readable format to
 * the stdout
 */
void ParameterFile::print_contents() {
    cout << "### time options ###" << endl;
    cout << "maxtime: " << _maxtime << endl;
    if(_global_timestep) {
        cout << "using a global timestep" << endl;
    } else {
        cout << "using individual timesteps" << endl;
    }
    cout << "maximal timestep: " << _max_timestep << endl;
    cout << "minimal timestep: " << _min_timestep << endl;
    if(_treetime) {
        cout << "using the treewalk timestep criterion" << endl;
    } else {
        cout << "using the per particle timestep criterion" << endl;
    }
    cout << endl;

    cout << "### snapshot options ###" << endl;
    cout << "basic snapshot name: " << _basename << endl;
    cout << "time between snapshots: " << _snaptime << endl;
    cout << "first snapshot: " << _firstsnap << endl;
    cout << "output directory: " << _outputdir << endl;
    cout << "snapshot type: " << _snapshot_type << endl;
    if(_per_node_output) {
        cout << "using separate output for every node" << endl;
    } else {
        cout << "using one output file for all nodes (if possible)" << endl;
    }
    cout << endl;

    cout << "### ic options ###" << endl;
    cout << "initial condition file: " << _icfile << endl;
    cout << "initial condition file type: " << _icfile_type << endl;
    cout << endl;

    cout << "### riemann solver options ###" << endl;
    cout << "type: " << _solver_type << endl;
    cout << "tolerance: " << _tolerance << endl;
    cout << "cutoff: " << _cutoff << endl;
    cout << "CFL: " << _CFL << endl;
    cout << endl;

    cout << "### hydro options ###" << endl;
    cout << "adiabatic index: " << _gamma << endl;
    cout << endl;

    cout << "### gravity options ###" << endl;
    if(_gravity) {
        cout << "Gravity enabled" << endl;
    } else {
        cout << "No gravity" << endl;
    }
    cout << "softening: " << _hsoft << endl;
    cout << "eta: " << _grav_eta << endl;
    cout << "alpha: " << _grav_alpha << endl;
    cout << endl;

    cout << "### voronoi options ###" << endl;
    cout << "tolerance: " << _voronoi_tolerance << endl;
    cout << endl;

    cout << "### tree options ###" << endl;
    cout << "ewald alpha: " << _ewald_alpha << endl;
    cout << "ewald size: " << _ewald_size << endl;
    cout << endl;

    cout << "### memory options ###" << endl;
    cout << "maximum size: "
         << HelperFunctions::human_readable_bytes(_max_memory) << endl;
    cout << endl;

    cout << "### code options ###" << endl;
    cout << "restart time: " << _restarttime << endl;
    cout << endl;

    cout << "### units ###" << endl;
    cout << "internal units: " << _units_internal << endl;
    cout << "output units: " << _units_output << endl;
    cout << endl;

    cout << "### physics ###" << endl;
    if(_real_physics) {
        cout << "Using real physical values for physical constants" << endl;
    } else {
        cout << "Using idealized values for physical constants" << endl;
    }
    cout << endl;
}

/**
 * @brief Constructor
 *
 * We try to open the file with the given name and initialize a
 * boost::property_tree based on its contents, assuming it is written in
 * .ini-format.
 *
 * We then parse all parameters and initialize them internally. Except for the
 * total simulation time, which needs to be specified, all parameters have
 * reasonable default values and can be omitted from the parameter file.
 *
 * @param name Filename of the parameter file
 */
ParameterFile::ParameterFile(std::string name) {
    ifstream file(name.c_str());
    if(!file) {
        cerr << "Cannot read parameterfile \"" << name << "\"!" << endl;
        my_exit();
    }

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(name, pt);

    // print out some info
    cout << "Read parameters from " << name << endl;

    // time options
    _maxtime = pt.get<double>("Time.MaxTime");
    _global_timestep = pt.get<bool>("Time.GlobalTimestep", false);
    _max_timestep = pt.get<double>("Time.MaxTimeStep", 0.);
    _min_timestep = pt.get<double>("Time.MinTimeStep", 0.);
    _treetime = pt.get<bool>("Time.TreeTime", false);

    // snapshot options
    _basename = pt.get<string>("Snapshots.BaseName", "snapshot");
    _snaptime = pt.get<double>("Snapshots.SnapTime", _maxtime * 0.1);
    _firstsnap = pt.get<unsigned int>("Snapshots.FirstSnap", 0);
    _outputdir = pt.get<string>("Snapshots.OutputDir", ".");
    // translate relative paths to absolute paths
    char* result = realpath(_outputdir.c_str(), NULL);
    _outputdir = string(result);
    // memory allocation is done using alloc, so we have to call free and not
    // delete
    free(result);
    _snapshot_type = pt.get<string>("Snapshots.Type", "Gadget");
    _per_node_output = pt.get<bool>("Snapshots.PerNodeOutput", false);

    // IC options
    stringstream default_ic;
    default_ic << _basename;
    default_ic.fill('0');
    default_ic.width(3);
    default_ic << _firstsnap;
    default_ic << ".hdf5";
    _icfile = pt.get<string>("IC.FileName", default_ic.str());
    _icfile_type = pt.get<string>("IC.Type", "Gadget");

    // Riemann solver options
    _solver_type = pt.get<string>("RiemannSolver.Type", "Exact");
    _tolerance = pt.get<double>("RiemannSolver.Tolerance", 1.e-8);
    _cutoff = pt.get<double>("RiemannSolver.CutOff", -5.);
    _CFL = pt.get<double>("RiemannSolver.CFL", 0.4);

    // Hydro options
    _gamma = pt.get<double>("Hydro.Gamma", 1.66667);

    // Gravity options
    _gravity = pt.get<bool>("Gravity.Gravity", true);
    _hsoft = pt.get<double>("Gravity.Softening", 0.03);
    _grav_eta = pt.get<double>("Gravity.Eta", 0.05 / 2.8);
    _grav_alpha = pt.get<double>("Gravity.Alpha", 0.005);

    // Voronoi options
    _voronoi_tolerance = pt.get<double>("Voronoi.Tolerance", 1.e-9);

    // Tree options
    _ewald_alpha = pt.get<double>("Tree.EwaldAlpha", 2.);
    _ewald_size = pt.get<unsigned int>("Tree.EwaldSize", 64);

    // Memory options
    string max_memory = pt.get<string>("Memory.MaximumSize", "1 GB");
    _max_memory = HelperFunctions::machine_readable_bytes(max_memory);

    // Code options
    _restarttime = pt.get<double>("Code.RestartTime", 3600.);

    // Units
    _units_internal = pt.get<string>("Units.InternalUnits", "SI");
    _units_output = pt.get<string>("Units.OutputUnits", "SI");

    // Physics options
    _real_physics = pt.get<bool>("Physics.RealPhysics", true);

    print_contents();
}

/**
 * @brief Get the total simulation time in internal units
 *
 * @return Total simulation time
 */
double ParameterFile::get_maxtime() {
    return _maxtime;
}

/**
 * @brief Check if we use a global timestep or individual timesteps
 *
 * @return True if we use a global timestep, false otherwise
 */
bool ParameterFile::get_global_timestep() {
    return _global_timestep;
}

/**
 * @brief Get the maximal Particle timestep in internal units
 *
 * @return Maximal Particle timestep
 */
double ParameterFile::get_max_timestep() {
    return _max_timestep;
}

/**
 * @brief Get the minimal Particle timestep in internal units
 *
 * @return Minimal Particle timestep
 */
double ParameterFile::get_min_timestep() {
    return _min_timestep;
}

/**
 * @brief Check if we use the expensive treewalk timestep criterion
 *
 * @return True if we use the treewalk criterion, false otherwise
 */
bool ParameterFile::has_treetime() {
    return _treetime;
}

/**
 * @brief Get the basic name of the snapshots
 *
 * The total name of the snapshot is basic_name + XXX + .hdf5.
 *
 * @return The basic name of the snapshots
 */
string ParameterFile::get_basename() {
    return _basename;
}

/**
 * @brief Get the time interval in between successive snapshots, in internal
 * units
 *
 * @return Time interval in between successive snapshots
 */
double ParameterFile::get_snaptime() {
    return _snaptime;
}

/**
 * @brief Get the index of the first snapshot to be written out
 *
 * The snapshot counter always starts from 0 at the starting time of the
 * simulation. Every time the snapshot time interval is reached, the counter is
 * increased. Only if the counter is equal or larger than the first index, an
 * actual snapshot is written out.
 *
 * @return Index of the first snapshot
 */
unsigned int ParameterFile::get_firstsnap() {
    return _firstsnap;
}

/**
 * @brief Get the output directory where snapshot files and other diagnostic
 * files are stored
 *
 * Also restart files are written to this directory.
 *
 * @return Output directory
 */
string ParameterFile::get_outputdir() {
    return _outputdir;
}

/**
 * @brief Get the snapshot type name
 *
 * @return Snapshot type name
 */
string ParameterFile::get_snapshot_type() {
    return _snapshot_type;
}

/**
 * @brief Check if every node should write a separate snapshot file, or all
 * nodes (try to) write to the same file
 *
 * @return True if nodes write separate snapshots, false otherwise
 */
bool ParameterFile::get_per_node_output() {
    return _per_node_output;
}

/**
 * @brief Get the filename of the initial condition file
 *
 * @return Initial condition filename
 */
string ParameterFile::get_icfile() {
    return _icfile;
}

/**
 * @brief Get the type name of the initial condition file
 *
 * @return Initial condition file type
 */
string ParameterFile::get_icfile_type() {
    return _icfile_type;
}

/**
 * @brief Get the type of the Riemann solver used to solve the Riemann problem
 *
 * @return Riemann solver type
 */
string ParameterFile::get_solver_type() {
    return _solver_type;
}

/**
 * @brief Get the tolerance used for the Riemann solver iteration
 *
 * @return Tolerance for the Riemann solver iteration
 */
double ParameterFile::get_tolerance() {
    return _tolerance;
}

/**
 * @brief Get the cut off value for the exact Riemann solver
 *
 * @return Cut off value for the exact Riemann solver
 */
double ParameterFile::get_cutoff() {
    return _cutoff;
}

/**
 * @brief Get the Courant-Friedrich-Levy constant for hydrodynamical time
 * stepping
 *
 * @return CFL constant for hydrodynamical time stepping
 */
double ParameterFile::get_CFL() {
    return _CFL;
}

/**
 * @brief Get the adiabatic index of the ideal gas
 *
 * @return Adiabatic index of the gas
 */
double ParameterFile::get_gamma() {
    return _gamma;
}

/**
 * @brief Check if we include gravity in the simulation
 *
 * @return True if gravity is enabled, false otherwise
 */
bool ParameterFile::has_gravity() {
    return _gravity;
}

/**
 * @brief Get the softening length used for gravitational force softening if no
 * individual softening lengths are specified
 *
 * @return Gravitational softening length
 */
double ParameterFile::get_hsoft() {
    return _hsoft;
}

/**
 * @brief Get the gravitational \f$\eta\f$ factor
 *
 * The \f$\eta\f$ factor is used in the gravitational time step calculation.
 *
 * @return Gravitational \f$\eta\f$ factor
 */
double ParameterFile::get_grav_eta() {
    return _grav_eta;
}

/**
 * @brief Get the gravitational \f$\alpha\f$ factor
 *
 * The \f$\alpha\f$ factor is used for the relative tree opening criterion
 * during the gravitational force calculation tree walk.
 *
 * @return Gravitational \f$\alpha\f$ factor
 */
double ParameterFile::get_grav_alpha() {
    return _grav_alpha;
}

/**
 * @brief Get the tolerance used to discriminate between exact and approximate
 * geometric tests during Voronoi grid construction
 *
 * @return Voronoi grid construction tolerance
 */
double ParameterFile::get_voronoi_tolerance() {
    return _voronoi_tolerance;
}

/**
 * @brief Get the Ewald \f$\alpha\f$ factor
 *
 * The \f$\alpha\f$ factor sets the cut off between long range and short range
 * forces in the Ewald force approximation.
 *
 * @return Ewald \f$\alpha\f$ factor
 */
double ParameterFile::get_ewald_alpha() {
    return _ewald_alpha;
}

/**
 * @brief Get the size of the Ewald table
 *
 * The size is the number of elements in the precalculated table in one
 * dimensions.
 *
 * @warning For 3D, this scales as N\f$^3\f$!
 *
 * @return Size of the Ewald table
 */
unsigned int ParameterFile::get_ewald_size() {
    return _ewald_size;
}

/**
 * @brief Get the maximal size of the MPI communication buffer in memory (in
 * bytes)
 *
 * @return The maximal size of the MPI communication buffer
 */
unsigned int ParameterFile::get_max_memory() {
    return _max_memory;
}

/**
 * @brief Get the time interval (in CPU time) in between two successive restart
 * file dumps
 *
 * @return Time interval in between restart file dumps
 */
double ParameterFile::get_restarttime() {
    return _restarttime;
}

/**
 * @brief Get the name of the internal UnitSet
 *
 * @return Internal units
 */
std::string ParameterFile::get_units_internal() {
    return _units_internal;
}

/**
 * @brief Get the name of the ouput UnitSet
 *
 * @return Output units
 */
std::string ParameterFile::get_units_output() {
    return _units_output;
}

/**
 * @brief Check if we use physical values for physical constants
 *
 * @return True if we use physical values, false for idealized values
 */
bool ParameterFile::has_real_physics() {
    return _real_physics;
}

/**
 * @brief Dump the parameterfile to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void ParameterFile::dump(RestartFile& rfile) {
    // time options
    rfile.write(_maxtime);
    rfile.write(_global_timestep);
    rfile.write(_max_timestep);
    rfile.write(_min_timestep);
    rfile.write(_treetime);

    // snapshot options
    rfile.write(_basename);
    rfile.write(_snaptime);
    rfile.write(_firstsnap);
    rfile.write(_outputdir);
    rfile.write(_snapshot_type);
    rfile.write(_per_node_output);

    // IC options
    rfile.write(_icfile);
    rfile.write(_icfile_type);

    // Riemann solver options
    rfile.write(_solver_type);
    rfile.write(_tolerance);
    rfile.write(_cutoff);
    rfile.write(_CFL);

    // Hydro options
    rfile.write(_gamma);

    // Gravity options
    rfile.write(_gravity);
    rfile.write(_hsoft);
    rfile.write(_grav_eta);
    rfile.write(_grav_alpha);

    // Voronoi options
    rfile.write(_voronoi_tolerance);

    // Tree options
    rfile.write(_ewald_alpha);
    rfile.write(_ewald_size);

    // Memory options
    rfile.write(_max_memory);

    // Code options
    rfile.write(_restarttime);

    // Units
    rfile.write(_units_internal);
    rfile.write(_units_output);

    // Physics options
    rfile.write(_real_physics);
}

/**
 * @brief Restart constructor. Initialize the restartfile from the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 */
ParameterFile::ParameterFile(RestartFile& rfile) {
    // time options
    rfile.read(_maxtime);
    rfile.read(_global_timestep);
    rfile.read(_max_timestep);
    rfile.read(_min_timestep);
    rfile.read(_treetime);

    // snapshot options
    rfile.read(_basename);
    rfile.read(_snaptime);
    rfile.read(_firstsnap);
    rfile.read(_outputdir);
    rfile.read(_snapshot_type);
    rfile.read(_per_node_output);

    // IC options
    rfile.read(_icfile);
    rfile.read(_icfile_type);

    // Riemann solver options
    rfile.read(_solver_type);
    rfile.read(_tolerance);
    rfile.read(_cutoff);
    rfile.read(_CFL);

    // Hydro options
    rfile.read(_gamma);

    // Gravity options
    rfile.read(_gravity);
    rfile.read(_hsoft);
    rfile.read(_grav_eta);
    rfile.read(_grav_alpha);

    // Voronoi options
    rfile.read(_voronoi_tolerance);

    // Tree options
    rfile.read(_ewald_alpha);
    rfile.read(_ewald_size);

    // Memory options
    rfile.read(_max_memory);

    // Code options
    rfile.read(_restarttime);

    // Units
    rfile.read(_units_internal);
    rfile.read(_units_output);

    // Physics options
    rfile.read(_real_physics);

    // print out info
    print_contents();
}
