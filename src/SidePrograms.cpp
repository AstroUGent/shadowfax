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
 * @file SidePrograms.cpp
 *
 * @brief A collection of programs that need the Voronoi tesselation but are not
 * really related to the main simulation program: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "SidePrograms.hpp"
#include "GravityWalker.hpp"
#include "ShadowfaxSnapshotReader.hpp"  // for ShadowfaxSnapshotReader
#include "ShadowfaxSnapshotWriter.hpp"
#include "SnapshotHandler.hpp"           // for SnapshotReader
#include "SnapshotReaderFactory.hpp"     // for SnapshotReaderFactory
#include "Vec.hpp"                       // for Vec
#include "VorCell.hpp"                   // for VorCell
#include "VorTess.hpp"                   // for VorTess
#include "io/AsciiInput.hpp"             // for AsciiInput
#include "io/Block.hpp"                  // for Block
#include "io/Header.hpp"                 // for Header
#include "io/Unit.hpp"                   // for Unit
#include "io/UnitSet.hpp"                // for UnitSet
#include "utilities/GasParticle.hpp"     // for GasParticle
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include "utilities/Tree.hpp"            // for Tree
#include <algorithm>                     // for max, min
#include <getopt.h>                      // for optarg, required_argument, etc
#include <iostream>                      // for operator<<, basic_ostream, etc
#include <stdlib.h>                      // for NULL, atof
#include <vector>                        // for vector
using namespace std;

/**
 * @brief Read in particle coordinates from the ASCII-file with the given name
 *
 * This method reads coordinates from an ASCII-file and stores them as particles
 * in the given std::vector. The method also calculates the extents of the box
 * containing all the coordinates.
 *
 * @param cells std::vector to store the particles in
 * @param cube_boundaries Array to store the box dimensions in
 * @param name Name of the ASCII-file
 */
void SidePrograms::load_ascii(ParticleVector& cells, double* cube_boundaries,
                              string name) {
    AsciiInput input(name);
    vector<string> headers;
    headers.push_back("x");
    headers.push_back("y");
#if ndim_ == 3
    headers.push_back("z");
#endif
    vector<unsigned int> dims(headers.size(), 1);
    vector<Unit> units(headers.size());
    Block block("position", headers, dims, units);
    input.read(block, 0);
    for(unsigned int i = block.get_size(); i--;) {
        vector<double> row = block.get_line(i);
#if ndim_ == 3
        Vec coords(row[0], row[1], row[2]);
#else
        Vec coords(row[0], row[1]);
#endif
        for(unsigned int j = ndim_; j--;) {
            cube_boundaries[2 * j] = std::min(cube_boundaries[2 * j], row[j]);
            cube_boundaries[2 * j + 1] =
                    std::max(cube_boundaries[2 * j + 1], row[j]);
        }
        cells.add_gas_particle(new GasParticle(coords));
        cells.gasback()->set_id(i);
    }
}

/**
 * @brief Calculate areas for all coordinates in the input file and write them
 * to the file areas.dat
 *
 * @warning Only works in 2 dimensions!
 *
 * This program accepts a single command line argument: the name of a file
 * containing particle coordinates.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
AreaCalculator::AreaCalculator(int argc, char** argv) : SidePrograms() {
#if ndim_ == 3
    cerr << "This option does not work in 3 dimensions. Sell a dimension and "
            "try again!"
         << endl;
#else
    string filename;
    Vec custom_origin;
    double custom_side = 0.;

    static struct option long_options[] = {
            {"filename", required_argument, NULL, 'f'},
            {"origin_x", required_argument, NULL, 'x'},
            {"origin_y", required_argument, NULL, 'y'},
            {"side", required_argument, NULL, 's'},
            {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":f:x:y:s:", long_options, NULL)) !=
          -1) {
        switch(c) {
            case 'f':
                filename = optarg;
                break;
            case 'x':
                custom_origin[0] = atof(optarg);
                break;
            case 'y':
                custom_origin[1] = atof(optarg);
                break;
            case 's':
                custom_side = atof(optarg);
                break;
            case ':':
                cerr << "Error! Missing required argument for option " << optopt
                     << "!" << endl;
                return;
        }
    }

    // command line errors
    if(!filename.size()) {
        cerr << "Error! No filename specified!" << endl;
        return;
    }

    ParticleVector particle_list(false, RectangularBox());
    double cube_boundaries[4] = {0.};
    cout << "Loading ascii file..." << endl;
    load_ascii(particle_list, cube_boundaries, filename);
    Vec side;
    for(unsigned int i = 2; i--;) {
        side[i] = cube_boundaries[2 * i + 1] - cube_boundaries[2 * i];
    }
    for(unsigned int i = 2; i--;) {
        cube_boundaries[2 * i] -= 0.01 * side[i];
    }
    side *= 1.02;
    Vec origin(cube_boundaries[0] + 0.5 * side[0],
               cube_boundaries[2] + 0.5 * side[1]);
    if(custom_side) {
        side.set(custom_side, custom_side);
        origin[0] = custom_origin[0] + 0.5 * side[0];
        origin[1] = custom_origin[1] + 0.5 * side[1];
    }
    cout << "Embedding particles in box with center (" << origin[0] << ","
         << origin[1] << ")";
    cout << " and sides (" << side.x() << "," << side.y() << ")" << endl;
    RectangularBox cube(origin, side);
    particle_list.set_container(cube);
    VorTess voronoi_tesselation(&cube, particle_list.gassize());
    cout << "Start constructing Voronoi mesh" << endl;
    cout << "Added 0 of " << particle_list.gassize() << " particles" << flush;

    // sort particles
    unsigned long width = 1;
    width <<= 10;
    for(unsigned int i = particle_list.gassize(); i--;) {
        particle_list.gas(i)->set_key(
                cube.get_key(particle_list.gas(i)->get_position()));
    }
    particle_list.sort();

    for(unsigned int i = particle_list.gassize(); i--;) {
        voronoi_tesselation.add_point(particle_list.gas(i), i);
        if(i % 10000 == 0) {
            cout << "\rAdded " << particle_list.gassize() - i << " of "
                 << particle_list.gassize() << " particles" << flush;
        }
    }
    cout << "\nFinalizing Voronoi mesh" << endl;
    voronoi_tesselation.complete(particle_list.get_tree());
    voronoi_tesselation.construct();
    cout << "Calculating areas" << endl;
    ofstream areafile("areas.dat");
    for(unsigned int i = particle_list.gassize(); i--;) {
        particle_list.gas(i)->set_key(particle_list.gas(i)->id());
    }
    // have to fix the piece of code below!!
    for(unsigned int i = 0; i < particle_list.gassize(); i++) {
        double dens = particle_list.gas(i)->get_cell()->get_volume();
        areafile << dens << endl;
    }
#endif
}

/**
 * @brief Calculate densities for all particles in the file and write them to
 * a snapshot with filename densities000.hdf5.
 *
 * @warning Only works for 3 dimensions!
 *
 * This program accepts a single command line argument: the name of the file to
 * read coordinates from.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
DensityCalculator::DensityCalculator(int argc, char** argv) : SidePrograms() {
#if ndim_ == 3
    string filename;

    static struct option long_options[] = {
            {"filename", required_argument, NULL, 'f'}, {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":f:", long_options, NULL)) != -1) {
        switch(c) {
            case 'f':
                filename = optarg;
                break;
            case ':':
                cerr << "Error! Missing required argument for option " << optopt
                     << "!" << endl;
                return;
        }
    }

    // command line errors
    if(!filename.size()) {
        cerr << "Error! No filename specified!" << endl;
        return;
    }

    ParticleVector particle_list(false, RectangularBox());
    double cube_boundaries[6] = {0.};
    cout << "Loading ascii file..." << endl;
    load_ascii(particle_list, cube_boundaries, filename);
    Vec side;
    for(unsigned int i = 3; i--;) {
        side[i] = cube_boundaries[2 * i + 1] - cube_boundaries[2 * i];
    }
    for(unsigned int i = 3; i--;) {
        cube_boundaries[2 * i] -= 0.01 * side[i];
    }
    side *= 1.02;
    cout << "Embedding particles in box with origin (" << cube_boundaries[0]
         << "," << cube_boundaries[2] << "," << cube_boundaries[4] << ")";
    cout << " and sides (" << side[0] << "," << side[1] << "," << side[2] << ")"
         << endl;
    Vec origin(cube_boundaries[0] + 0.5 * side[0],
               cube_boundaries[2] + 0.5 * side[1],
               cube_boundaries[4] + 0.5 * side[2]);
    RectangularBox cube(origin, side);
    particle_list.set_container(cube);
    VorTess voronoi_tesselation(&cube, particle_list.gassize(), false);
    cout << "Start constructing Voronoi mesh" << endl;
    cout << "Added 0 of " << particle_list.gassize() << " particles" << flush;

    // sort particles
    unsigned long width = 1;
    width <<= 10;
    for(unsigned int i = particle_list.gassize(); i--;) {
        unsigned long bits[ndim_];
        for(unsigned int j = ndim_; j--;) {
            bits[j] = ((particle_list.gas(i)->pos(j) - cube_boundaries[2 * j]) /
                       side[j]) *
                      width;
        }
        particle_list.gas(i)->set_key(HB::get_key(bits, 30));
    }
    particle_list.sort();

    for(unsigned int i = particle_list.gassize(); i--;) {
        voronoi_tesselation.add_point(particle_list.gas(i), i);
        if(i % 10000 == 0) {
            cout << "\rAdded " << particle_list.gassize() - i << " of "
                 << particle_list.gassize() << " particles" << flush;
        }
    }
    cout << "\nFinalizing Voronoi mesh" << endl;
    voronoi_tesselation.construct();
    unsigned int highdens_index = 0;
    double highdens = 0.;
    cout << "Calculating densities" << endl;
    for(unsigned int i = particle_list.gassize(); i--;) {
        double dens = 1. / particle_list.gas(i)->get_cell()->get_volume();
        if(dens > highdens) {
            highdens = dens;
            highdens_index = i;
        }
        StateVector W(dens, 0., 0., 0., 0.);
        particle_list.gas(i)->set_W(W);
    }
    // save the densities to hdf5
    double box[ndim_ + 1] = {0.};
    cube.get_bounding_box(box);
    Header header;
    header.set_ngaspart(particle_list.gassize());
    header.set_box(box);
    header.set_time(0);
    header.set_global_timestep(false);
    UnitSet units_in;
    UnitSet units_out;
    ShadowfaxSnapshotWriter snapshothandler("densities", units_in, units_out,
                                            0);
    snapshothandler.write_snapshot(0, particle_list);
    cout << "The particle with coordinates ("
         << particle_list.gas(highdens_index)->x() << ","
         << particle_list.gas(highdens_index)->y() << ","
         << particle_list.gas(highdens_index)->z()
         << ") has the highest density" << endl;
    cout << "This is the " << particle_list.gas(highdens_index)->id()
         << "th particle in the file" << endl;
#else
    cerr << "This option does not work in 2 dimensions. Buy another dimension "
            "and try again!"
         << endl;
#endif
}

/**
 * @brief Sort all particles in the given file on Hilbert key and write them in
 * sorted order to a new file with the name FILENAME_sorted.txt
 *
 * @warning Only works in 3 dimensions!
 *
 * This program accepts a single command line argument: the name of a file to
 * read coordinates from.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
HilbertSorter::HilbertSorter(int argc, char** argv) {
#if ndim_ == 3
    string filename;
    string onlyname;

    static struct option long_options[] = {
            {"filename", required_argument, NULL, 'f'}, {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":f:", long_options, NULL)) != -1) {
        switch(c) {
            case 'f':
                filename = optarg;
                // drop the last 5 characters from the name by setting the end
                // of
                // string
                optarg[filename.size() - 5] = '\0';
                onlyname = optarg;
                break;
            case ':':
                cerr << "Error! Missing required argument for option " << optopt
                     << "!" << endl;
                return;
        }
    }

    // command line errors
    if(!filename.size()) {
        cerr << "Error! No filename specified!" << endl;
        return;
    }

    cout << "Sorting " << onlyname << "..." << endl;

    ParticleVector particles(false, RectangularBox());
    double box[6];
    load_ascii(particles, box, filename);
    Vec origin(0., 0., 0.);
    Vec side(1., 1., 1.);
    RectangularBox cube(origin, side);
    particles.set_container(cube);

    particles.sort();

    ofstream ofile((onlyname + string("_sorted.txt")).c_str());
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        GasParticle* p = particles.gas(i);
        ofile << p->x() << "\t" << p->y() << "\t" << p->z() << "\n";
    }
#else
    cerr << "This option does only work in 3 dimensions!" << endl;
#endif
}

/**
 * @brief Calculate masses for all particles in the file and write them to a
 * file with the name FILENAME.mass
 *
 * This program accepts a single command line parameter: the name of a snapshot
 * file.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
MassCalculator::MassCalculator(int argc, char** argv) {
    string filename;
    string onlyname;
    string type = "Gadget";

    static struct option long_options[] = {
            {"filename", required_argument, NULL, 'f'},
            {"type", required_argument, NULL, 't'},
            {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":f:t:", long_options, NULL)) != -1) {
        switch(c) {
            case 'f':
                filename = optarg;
                // drop the last 5 characters from the name by setting the end
                // of
                // string
                optarg[filename.size() - 5] = '\0';
                onlyname = optarg;
                break;
            case 't':
                type = optarg;
                break;
            case ':':
                cerr << "Error! Missing required argument for option " << optopt
                     << "!" << endl;
                return;
        }
    }

    // command line errors
    if(!filename.size()) {
        cerr << "Error! No filename specified!" << endl;
        return;
    }

    cout << "Calculating masses for " << onlyname << "..." << endl;

    RectangularBox container;
    ParticleVector particles(false, container);
    UnitSet simulation_units;
    SnapshotReader* reader =
            SnapshotReaderFactory::generate(type, filename, simulation_units);
    Header header = reader->read_snapshot(particles);
    bool periodic = header.periodic();
    particles.set_periodic(periodic);

    double box[ndim_ + ndim_] = {0.};
    header.box(box);
    Vec center;
    Vec sides;
    for(unsigned int i = ndim_; i--;) {
        center[i] = box[i] + 0.5 * box[ndim_ + i];
        sides[i] = box[ndim_ + i];
    }
    container = RectangularBox(center, sides);
    particles.set_container(container);

    particles.sort();

    VorTess tesselation(&particles.get_container(), particles.gassize(),
                        periodic);
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        tesselation.add_point(particles.gas(i), particles.gas(i)->id());
        particles.gas(i)->reset_copies();
        particles.gas(i)->reset_export();
    }
    tesselation.complete(particles.get_tree());
    tesselation.construct();

    ofstream ofile(onlyname + string(".mass"));
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        double mass = particles.gas(i)->get_cell()->get_volume() *
                      particles.gas(i)->get_Wvec().rho();
        ofile << particles.gas(i)->id() << "\t" << mass << "\n";
    }
}

/**
 * @brief Calculate gravitational potential energies for all particles in the
 * file and write them to a file with the name FILENAME.epot
 *
 * This program accepts a single command line parameter: the name of a snapshot
 * file.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
PotentialCalculator::PotentialCalculator(int argc, char** argv) {
    string filename;
    string onlyname;

    static struct option long_options[] = {
            {"filename", required_argument, NULL, 'f'}, {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":f:", long_options, NULL)) != -1) {
        switch(c) {
            case 'f':
                filename = optarg;
                // drop the last 5 characters from the name by setting the end
                // of
                // string
                optarg[filename.size() - 5] = '\0';
                onlyname = optarg;
                break;
            case ':':
                cerr << "Error! Missing required argument for option " << optopt
                     << "!" << endl;
                return;
        }
    }

    // command line errors
    if(!filename.size()) {
        cerr << "Error! No filename specified!" << endl;
        return;
    }

    cout << "Calculating potential energies for " << onlyname << "..." << endl;

    RectangularBox container;
    ParticleVector particles(false, container);
    UnitSet simulation_units;
    ShadowfaxSnapshotReader reader(filename, simulation_units);
    Header header = reader.read_snapshot(particles);
    bool periodic = header.periodic();
    particles.set_periodic(periodic);

    double box[ndim_ + ndim_] = {0.};
    header.box(box);
    Vec center;
    Vec sides;
    for(unsigned int i = ndim_; i--;) {
        center[i] = box[i] + 0.5 * box[ndim_ + i];
        sides[i] = box[ndim_ + i];
    }
    container = RectangularBox(center, sides);
    particles.set_container(container);

    particles.sort();

    if(particles.gassize()) {
        cout << "Calculating Voronoi grid to obtain masses..." << endl;
        VorTess tesselation(&particles.get_container(), particles.gassize(),
                            periodic);
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            tesselation.add_point(particles.gas(i), particles.gas(i)->id());
            particles.gas(i)->reset_copies();
            particles.gas(i)->reset_export();
        }
        tesselation.complete(particles.get_tree());
        tesselation.construct();
        cout << "Done." << endl;

        // set particle masses and softening lengths
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            double mass = particles.gas(i)->get_cell()->get_volume() *
                          particles.gas(i)->get_Wvec().rho();
            particles.gas(i)->set_mass(mass);
            particles.gas(i)->set_hsoft(0.01);
        }
    }
    if(particles.dmsize()) {
        for(unsigned int i = 0; i < particles.dmsize(); i++) {
            // by using a slightly smaller softening length, we reach higher
            // accuracies
            // this compensates for the fact that we use a different opening
            // criterion during the simulations
            particles.dm(i)->set_hsoft(0.5 * header.hsoft());
        }
    }

    // update the masses in the tree
    particles.get_tree().set_velocities();

    cout << "Performing treewalk to obtain gravitational potentials..." << endl;

    // walk the tree to calculate the potentials
    particles.get_tree().walk_tree<BHPotentialWalker>(particles, true, true,
                                                      true, 0);
    cout << "Done." << endl;

    cout << "Writing file " << onlyname << ".epot..." << endl;
    ofstream ofile(onlyname + string(".epot"));
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        double epot = 0.5 * particles.gas(i)->get_mass() *
                      particles.gas(i)->get_gravitational_potential();
        // need to add G for real simulations!!
        ofile << particles.gas(i)->id() << "\t" << epot << "\n";
    }
    for(unsigned int i = 0; i < particles.dmsize(); i++) {
        double epot = 0.5 * particles.dm(i)->get_mass() *
                      particles.dm(i)->get_gravitational_potential();
        // need to add G for real simulations!!
        ofile << particles.dm(i)->id() << "\t" << epot << "\n";
    }
    cout << "Done." << endl;
}
