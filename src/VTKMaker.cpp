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
 * @file VTKMaker.cpp
 *
 * @brief Sideprogram to convert snapshots to VTK-files: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "VTKMaker.hpp"
#include "MPIMethods.hpp"             // for MyMPI_Finalize, MyMPI_Init
#include "SnapshotHandler.hpp"        // for SnapshotReader
#include "SnapshotReaderFactory.hpp"  // for SnapshotReaderFactory
#include "StateVector.hpp"            // for StateVector
#include "Vec.hpp"                    // for Vec
#include "VorTess.hpp"                // for VorTess
#include "io/Header.hpp"              // for Header
#include "io/UnitSet.hpp"             // for UnitSet
#include "utilities/GasParticle.hpp"
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include <getopt.h>                      // for optarg, getopt_long, etc
#include <iostream>                      // for cerr
#include <stddef.h>                      // for NULL
#include <string>                        // for string
#include <vector>                        // for vector
using namespace std;

/**
 * @brief Check if the system uses big endian or little endian byte ordering
 *
 * Code based on similar code in the VisIt VTK writer.
 */
void VTKMaker::check_big_endian() {
    // assign an integer value
    int a = 1;
    // cast it to a char to access the first of the 4 bytes
    char* b = reinterpret_cast<char*>(&a);
    // if the first byte contains the value of 1, integers are stored in little
    // endian order
    if(b) {
        _big_endian = false;
    } else {
        _big_endian = true;
    }
}

/**
 * @brief Convert the given byte array from little endian to big endian by
 * reversing it
 *
 * We assume the buffer is subdivided in a given number of groups of 4 bytes.
 *
 * @param buffer Buffer to invert
 * @param length Number of 4-byte groups in the given buffer
 */
void VTKMaker::swap_endian(char* buffer, unsigned int length) {
    for(unsigned i = 0; i < length; i++) {
        char tmp = buffer[4 * i];
        buffer[4 * i] = buffer[4 * i + 3];
        buffer[4 * i + 3] = tmp;
        tmp = buffer[4 * i + 1];
        buffer[4 * i + 1] = buffer[4 * i + 2];
        buffer[4 * i + 2] = tmp;
    }
}

/**
 * @brief Write the given byte array to the given binary stream
 *
 * @param stream Binary std::ostream to write to
 * @param buffer Byte array
 * @param length Number of 4-byte groups in the given array
 */
void VTKMaker::write_big_endian_buffer(std::ostream& stream, char* buffer,
                                       unsigned int length) {
    if(!_big_endian) {
        swap_endian(buffer, length);
    }
    stream.write(buffer, length * 4);
}

/**
 * @brief Write the given float array to the given binary stream
 *
 * @param stream Binary std::ostream to write to
 * @param buffer Float array
 * @param length Length of the array
 */
void VTKMaker::write_big_endian_buffer(std::ostream& stream, float* buffer,
                                       unsigned int length) {
    write_big_endian_buffer(stream, reinterpret_cast<char*>(buffer), length);
}

/**
 * @brief Write the given integer array to the given binary stream
 *
 * @param stream Binary std::ostream to write to
 * @param buffer Integer array
 * @param length Length of the array
 */
void VTKMaker::write_big_endian_buffer(std::ostream& stream, int* buffer,
                                       unsigned int length) {
    write_big_endian_buffer(stream, reinterpret_cast<char*>(buffer), length);
}

/**
 * @brief Constructor
 *
 * Read in the given snapshot, construct a Voronoi tesselation for it and dump
 * the Voronoi tesselation to a .vtk file with the same name as the snapshot.
 *
 * This program accepts two command line arguments:
 *  -# The name of the snapshot file to read
 *  -# The type of the snapshot file (optional, default value: Gadget)
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
VTKMaker::VTKMaker(int argc, char** argv) {
    string readertype = "Gadget";
    string filename;
    string onlyname;
    bool voronoi = true;

    static struct option long_options[] = {
            {"type", required_argument, NULL, 't'},
            {"filename", required_argument, NULL, 'o'},
            {"delaunay", no_argument, NULL, 'd'},
            {"voronoi", no_argument, NULL, 'v'},
            {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":t:o:dv", long_options, NULL)) != -1) {
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
            case 'd':
                voronoi = false;
                break;
            case 'v':
                // not really necessary, but it is nicer to have this option as
                // well
                // for clarity
                voronoi = true;
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
        tesselation.add_point(particles.gas(i),
                              particles.gas(i)->get_local_id());
        particles.gas(i)->reset_copies();
        particles.gas(i)->reset_export();
    }

    onlyname += string(".vtk");
    ofstream stream(onlyname);
    if(voronoi) {
        // calculate the Voronoi mesh and output it
        tesselation.complete(particles.get_tree());
        tesselation.construct();

        vector<float> positions;
        vector<int> connectivity;
        vector<StateVector> data;
        tesselation.get_triangles(positions, connectivity, data);
        stream << "# vtk DataFile Version 2.0\n";
        stream << "some label\n";
        stream << "BINARY\n";
        stream << "DATASET UNSTRUCTURED_GRID\n";
        stream << "POINTS " << (positions.size() / 3) << " float\n";
        write_big_endian_buffer(stream, &positions[0], positions.size());
        stream << "CELLS " << data.size() << " " << connectivity.size() << "\n";
        write_big_endian_buffer(stream, &connectivity[0], connectivity.size());
        stream << "CELL_TYPES " << data.size() << "\n";
        vector<int> celltype(data.size(), 7);
        write_big_endian_buffer(stream, &celltype[0], celltype.size());
        vector<float> densitydata(data.size(), 0.);
        vector<float> velocitydata(data.size() * 3, 0.);
        vector<float> pressuredata(data.size(), 0.);
        for(unsigned int i = 0; i < data.size(); i++) {
            densitydata[i] = data[i].rho();
            velocitydata[3 * i] = data[i].vx();
            velocitydata[3 * i + 1] = data[i].vy();
#if ndim_ == 3
            velocitydata[3 * i + 2] = data[i].vz();
#else
            velocitydata[3 * i + 2] = 0.;
#endif
            pressuredata[i] = data[i].p();
        }
        stream << "CELL_DATA " << data.size() << "\n";
        stream << "SCALARS density float\n";
        stream << "LOOKUP_TABLE default\n";
        write_big_endian_buffer(stream, &densitydata[0], densitydata.size());
        stream << "SCALARS pressure float\n";
        stream << "LOOKUP_TABLE default\n";
        write_big_endian_buffer(stream, &pressuredata[0], pressuredata.size());
        stream << "VECTORS velocity float\n";
        write_big_endian_buffer(stream, &velocitydata[0], velocitydata.size());
        stream << "POINT_DATA " << data.size() << "\n";
    } else {
        // stick with the Delaunay tesselation and output it
        vector<float> positions;
        vector<int> connectivity;
        tesselation.get_delaunay_triangles(positions, connectivity);
        stream << "# vtk DataFile Version 2.0\n";
        stream << "some label\n";
        stream << "BINARY\n";
        stream << "DATASET UNSTRUCTURED_GRID\n";
        stream << "POINTS " << (positions.size() / 3) << " float\n";
        write_big_endian_buffer(stream, &positions[0], positions.size());
        stream << "CELLS " << connectivity.size() / 4 << " "
               << connectivity.size() << "\n";
        write_big_endian_buffer(stream, &connectivity[0], connectivity.size());
        stream << "CELL_TYPES " << connectivity.size() / 4 << "\n";
        vector<int> celltype(connectivity.size() / 4, 7);
        write_big_endian_buffer(stream, &celltype[0], celltype.size());
    }
}

/**
 * @brief Entrance point for the vtkmaker program
 *
 * Calls the VTKMaker constructor.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    MyMPI_Init(&argc, &argv);

    VTKMaker(argc, argv);

    return MyMPI_Finalize();
}
