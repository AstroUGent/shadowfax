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
 * @file Plotter.cpp
 *
 * @brief SideProgram to plot Shadowfax snapshots: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Plotter.hpp"
#include "DelCont.hpp"               // for RectangularBox
#include "Error.hpp"                 // for my_exit
#include "GadgetSnapshotReader.hpp"  // for GadgetSnapshotReader
#include "Image.hpp"                 // for Image, ImageType::BIN_PGM, etc
#include "ImageDescription.hpp"
#include "MPIMethods.hpp"                // for MyMPI_Barrier
#include "Vec.hpp"                       // for Vec
#include "io/Header.hpp"                 // for Header
#include "io/UnitSet.hpp"                // for UnitSet
#include "utilities/GasParticle.hpp"     // for GasParticle
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include <cstdlib>                       // for atof
#include <cstring>                       // for strlen
#include <getopt.h>
#include <iostream>  // for operator<<, ostream, cout, etc
#include <vector>    // for vector
class Particle;
using namespace std;

/**
 * @brief Constructor
 *
 * Load the snapshot specified in the command line arguments and make an Image
 * showing one of the hydrodynamical quantities, which is saved as a PPM-file.
 *
 * This program takes up to 12 command line arguments:
 *  -# The name of the snapshot to be read
 *  -# The name of the variable to be plotted
 *  -# The maximal value of the plot interval (if not set, the maximum value of
 *     the data is used)
 *  -# The minimal value of the plot interval (if not set, the minimal value of
 *     the data is used)
 *  -# Flag indicating if a logarithmic scale should be used (True/true/1). If
 *     no (valid) value is provided, a linear scale is used
 *  -# Flag indicating if the Voronoi grid should be plotted on top of the
 *     image (grid/Grid). If no (valid) value is specified, the grid is not
 *     overplotted
 *  -# Flag indicating if a grayscale should be used or not (gray/Gray). If no
 *     (valid) value is provided, a color scale is used
 *  -# x-coordinate of the offset of the image in the simulation box
 *  -# y-coordinate of the offset of the image in the simulation box
 *  -# Height of the image in physical size
 *  -# Width of the image in physical size
 *  -# Physical side length of a single pixel of the image
 *
 * Only the first argument is required, for the others reasonable default values
 * are provided. If one of the last five is provided, all five should be
 * provided.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
Plotter::Plotter(int argc, char** argv) {
    std::string filename;
    std::string parname;

    static struct option long_options[] = {
            {"filename", required_argument, NULL, 'f'},
            {"parname", required_argument, NULL, 'p'},
            {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":f:p:", long_options, NULL)) != -1) {
        switch(c) {
            case 'f':
                filename = optarg;
                break;
            case 'p':
                parname = optarg;
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

// set up random grid (test)
#if ndim_ == 3
    Vec pos(0.5, 0.5, 0.5);
    Vec sides(1., 1., 1.);
#else
    Vec pos(0.5, 0.5);
    Vec sides(1., 1.);
#endif
    RectangularBox cont(pos, sides);
    vector<Particle*> plist;

    ImageDescription description(parname);

    UnitSet units;
    GadgetSnapshotReader reader(filename, units);
    ParticleVector tempvector(true, cont);
    Header header = reader.read_snapshot(tempvector);
    for(unsigned int i = 0; i < tempvector.gassize(); i++) {
        plist.push_back(tempvector.gas(i));
    }

    double box[ndim_ + ndim_] = {0.};
    header.box(box);
    Vec center;
    Vec side;
    for(unsigned int i = ndim_; i--;) {
        center[i] = box[i] + 0.5 * box[ndim_ + i];
        side[i] = box[ndim_ + i];
    }
    cont = RectangularBox(center, side);

    double snaptime = header.time();
    cout << "snapshot time " << snaptime << endl;
    bool grid = description.plot_grid();
    // full argument example that zooms in on a smaller region
    //    Image img(plist, cont, HYDRO_DENSITY, grid, 0., 0.3, 0.25, 0.25,
    // 0.00025);
    PlotVariable var = HYDRO_DENSITY;
    string varname = description.get_variable_name();
    if(varname == "velocity_x") {
        var = HYDRO_VELOCITY_X;
        cout << "plotting velocity_x" << endl;
    } else if(varname == "velocity_y") {
        var = HYDRO_VELOCITY_Y;
        cout << "plotting velocity_y" << endl;
    } else if(varname == "pressure") {
        var = HYDRO_PRESSURE;
        cout << "plotting pressure" << endl;
    } else {
        cout << "plotting density" << endl;
    }

    double offset[2];
    offset[0] = description.get_offset_x();
    offset[1] = description.get_offset_y();
    double size[2];
    size[0] = description.get_size_x();
    size[1] = description.get_size_y();
    double pixelsize = description.get_pixelsize();

    if(offset[0] == -1.) {
        offset[0] = box[0];
    }
    if(offset[1] == -1.) {
        offset[1] = box[1];
    }
    if(size[0] == -1.) {
        size[0] = side[0];
    }
    if(size[1] == -1.) {
        size[1] = side[1];
    }
    if(pixelsize == -1.) {
        pixelsize = std::max(size[0], size[1]) * 0.001;
    }

    Image img(plist, cont, header.periodic(), var, grid, offset[0], offset[1],
              size[0], size[1], pixelsize);
    string densgridfilename = get_filename(filename);
    double max = description.get_range_max();
    double min = description.get_range_min();
    bool log = description.is_log();
    int imgtype = BIN_PPM | CM_JET;
    if(!description.is_color()) {
        imgtype = BIN_PGM;
    }
    cout << "plotting with " << imgtype << ", " << max << ", " << min << ", "
         << log << endl;
    img.save(densgridfilename, imgtype, max, min, log);
    // we block here to prevent the program to finish before process 0 is
    // finished
    MyMPI_Barrier();
}

/**
 * @brief Subtract the .hdf5 file extension from the given filename
 *
 * @param filename Filename ending with .hdf5
 * @return Filename without the extension
 */
string Plotter::get_filename(string filename) {
    unsigned int i = strlen(filename.c_str());
    return filename.substr(0, i - 5);
}
