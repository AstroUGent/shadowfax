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
 * @file Lloyd.cpp
 *
 * @brief Auxiliary program to compute the centroids of the Voronoi mesh:
 * implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Lloyd.hpp"
#include "MPIMethods.hpp"                // for MyMPI_Finalize, MyMPI_Init
#include "SnapshotHandler.hpp"           // for SnapshotReader
#include "SnapshotReaderFactory.hpp"     // for SnapshotReaderFactory
#include "Vec.hpp"                       // for Vec
#include "VorCell.hpp"                   // for VorCell
#include "VorTess.hpp"                   // for VorTess
#include "io/Header.hpp"                 // for Header
#include "io/UnitSet.hpp"                // for UnitSet
#include "utilities/GasParticle.hpp"     // for GasParticle
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include <cstddef>                       // for NULL
#include <getopt.h>                      // for optarg, getopt_long, etc
#include <iostream>                      // for operator<<, basic_ostream, etc
#include <sstream>
#include <string>  // for string, char_traits, etc
using namespace std;

/**
 * @brief Calculate the centroids of the Voronoi mesh for the given generator
 * coordinates in the given box
 *
 * @param coords Coordinates of the mesh generators, a vector of ndim*N
 * coordinates, with ndim 2 or 3 and N the number of generators. Coordinates
 * should be given in x, y (, z) order
 * @param box_origin Coordinates of the bottom left corner of the box
 * @param box_sides Side lengths of the box in all dimensions
 * @param periodic Flag indicating whether the box is periodic or not
 * @return The ndim*N coordinates of the centroids of the Voronoi mesh
 */
vector<double> Lloyd::calculate_centroids(vector<double>& coords,
                                          vector<double>& box_origin,
                                          vector<double>& box_sides,
                                          bool periodic) {
    RectangularBox container;
    ParticleVector particles(false, container);
    particles.set_periodic(periodic);

    for(unsigned int i = 0; i < coords.size() / ndim_; i++) {
#if ndim_ == 3
        Vec pos(coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
#else
        Vec pos(coords[2 * i], coords[2 * i + 1]);
#endif
        GasParticle* part = new GasParticle(pos);
        part->set_id(i);
        particles.add_gas_particle(part);
    }

#if ndim_ == 3
    double box[6] = {box_origin[0], box_origin[1], box_origin[2],
                     box_sides[0],  box_sides[1],  box_sides[2]};
#else
    double box[4] = {box_origin[0], box_origin[1], box_sides[0], box_sides[1]};
#endif
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

    vector<double> centroids(ndim_ * particles.gassize(), 0.);
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        Vec centroid = tesselation.get_cell(i)->get_centroid();
        centroids[ndim_ * i] = centroid.x();
        centroids[ndim_ * i + 1] = centroid.y();
#if ndim_ == 3
        centroids[3 * i + 2] = centroid.z();
#endif
    }

    return centroids;
}

/**
 * @brief Constructor and main program
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
Lloyd::Lloyd(int argc, char** argv) {
    string readertype = "Gadget";
    string filename;
    string onlyname;

    static struct option long_options[] = {
            {"type", required_argument, NULL, 't'},
            {"filename", required_argument, NULL, 'o'},
            {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":t:o:", long_options, NULL)) != -1) {
        switch(c) {
            case 't':
                readertype = optarg;
                break;
            case 'o':
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

    RectangularBox container;
    ParticleVector particles(false, container);
    UnitSet simulation_units;
    SnapshotReader* reader = SnapshotReaderFactory::generate(
            readertype, filename, simulation_units);
    Header header = reader->read_snapshot(particles);
    delete reader;
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

    stringstream oname;
    oname << onlyname << ".centroid";
    ofstream ofile(oname.str().c_str());
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        Vec centroid = tesselation.get_cell(i)->get_centroid();
        ofile << centroid.x() << "\t" << centroid.y();
#if ndim_ == 3
        ofile << "\t" << centroid.z();
#endif
        ofile << "\n";
    }
}

#ifndef PYTHON_MODULE
/**
 * @brief Entrance point for the lloyd program
 *
 * Calls the Lloyd constructor.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    MyMPI_Init(&argc, &argv);

    Lloyd(argc, argv);

    return MyMPI_Finalize();
}
#endif
