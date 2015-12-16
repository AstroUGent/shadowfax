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
#include "Vec.hpp"
#include "DelCont.hpp"
#include "VorTess.hpp"
#include "Image.hpp"
#include "ShadowfaxSnapshotReader.hpp"
#include "Error.hpp"
#include "MPIMethods.hpp"
#include "io/Header.hpp"
#include "io/Input.hpp"
#include "io/UnitSet.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/ParticleVector.hpp"
#include <vector>
#include <sstream>
#include <cstring>
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
Plotter::Plotter(int argc, char **argv){
    // set up random grid (test)
#if ndim_==3
    Vec pos(0.5, 0.5, 0.5);
    Vec sides(1.,1.,1.);
#else
    Vec pos(0.5, 0.5);
    Vec sides(1.,1.);
#endif
    RectangularBox cont(pos, sides);
    vector<Particle*> plist;

    if(argc < 3){
        cout << "not enough arguments, stopping" << endl;
        return;
    }
    if(string(argv[2]) == "help" || string(argv[2]) == "Help"){
        cout << "\nPlot a snapshot.\n\nSYNTAX:\n\n./shadowfax plot NUMBER MAX "
                "MIN LOG GRID GRAYSCALE\n\n";
        cout << "\tNUMBER:\nnumber of snapshot file (assumed to have name "
                "snapshot_small_NUMBER.hdf5, with number having 3 digits)\n\n";
        cout << "\tMAX:\nmaximal value of density (MAX=-1. uses the largest "
                "value in the snapshot as maximal value)\n\n";
        cout << "\tMIN:\nminimal value of density (MIN=-1. uses the smallest "
                "value in the snapshot as minimal value)\n\n";
        cout << "\tLOG:\nboolean specifying if the plot has to be logarithmic "
                "in density. True, true or 1 means logarithmic plot, all other "
                "values give a linear plot\n\n";
        cout << "\tGRID:\nspecify if the grid has to be overplotted on the "
                "image. Accepts Grid or grid, all other values mean no grid is "
                "plotted.\n\n";
        cout << "\tGRAYSCALE:\nspecify if the image has to be plotted in "
                "grayscales. Accepts Gray or gray, all other values give a "
                "colored plot.\n";
        return;
    }
    cout << "plotting snapshot " << argv[2] << endl;

    string filename = string(argv[2]);
    UnitSet units;
    ShadowfaxSnapshotReader reader(filename, units);
    ParticleVector tempvector(true, cont);
    Header header = reader.read_snapshot(tempvector);
    for(unsigned int i = 0; i < tempvector.gassize(); i++){
        plist.push_back(tempvector.gas(i));
    }

    double box[ndim_+ndim_] = {0.};
    header.box(box);
    Vec center;
    Vec side;
    for(unsigned int i = ndim_; i--;){
        center[i] = box[i] + 0.5*box[ndim_+i];
        side[i] = box[ndim_+i];
    }
    cont = RectangularBox(center, side);

    double snaptime = header.time();
    cout << "snapshot time " << snaptime << endl;
    bool grid = false;
    if(argc > 7){
        if(string(argv[7]) == "grid" || string(argv[7]) == "Grid"){
            grid = true;
        }
    }
    // full argument example that zooms in on a smaller region
    PlotVariable var = HYDRO_DENSITY;
    if(argc > 3){
        if(string(argv[3]) == "velocity_x"){
            var = HYDRO_VELOCITY_X;
            cout << "plotting velocity_x" << endl;
        } else
        if(string(argv[3]) == "velocity_y"){
            var = HYDRO_VELOCITY_Y;
            cout << "plotting velocity_y" << endl;
        } else
        if(string(argv[3]) == "pressure"){
            var = HYDRO_PRESSURE;
            cout << "plotting pressure" << endl;
        } else {
            cout << "plotting density" << endl;
        }
    } else {
        cout << "plotting density" << endl;
    }
    double imagepars[5] = {0., 0., 1., 1., 0.001};
    if(argc > 9){
        if(argc < 14){
            cerr << "Error! If you provide one image parameter, you should "
                    "provide all of them!" << endl;
            my_exit();
        }
        imagepars[0] = atof(argv[9]);
        imagepars[1] = atof(argv[10]);
        imagepars[2] = atof(argv[11]);
        imagepars[3] = atof(argv[12]);
        imagepars[4] = atof(argv[13]);
    }
    Image img(plist, cont, header.periodic(), var, grid, imagepars[0],
              imagepars[1], imagepars[2], imagepars[3], imagepars[4]);
    string densgridfilename = get_filename(filename);
    double max = -1.;
    double min = -1.;
    bool log = false;
    int imgtype = BIN_PPM | CM_JET;
    if(argc > 4){
        max = atof(argv[4]);
    }
    if(argc > 5){
        min = atof(argv[5]);
    }
    if(argc > 6){
        log = string(argv[6]) == "True" || string(argv[6]) == "true" ||
                string(argv[6]) == "1";
    }
    if(argc > 8){
        if(string(argv[8]) == "gray" || string(argv[8]) == "Gray"){
            imgtype = BIN_PGM;
        }
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
string Plotter::get_filename(string filename){
    unsigned int i = strlen(filename.c_str());
    return filename.substr(0, i-5);
}
