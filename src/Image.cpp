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
 * @file Image.cpp
 *
 * @brief PPM image: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Image.hpp"
#include "ColorMap.hpp"               // for ColorMap
#include "DelCont.hpp"                // for DelCont
#include "Error.hpp"                  // for my_exit
#include "MPIGlobal.hpp"              // for rank
#include "StateVector.hpp"            // for StateVector
#include "Vec.hpp"                    // for Vec
#include "utilities/GasParticle.hpp"  // for GasParticle
#include "utilities/Particle.hpp"     // for Particle
#include "utilities/Tree.hpp"         // for Tree
#include <algorithm>                  // for min, max
#include <cmath>                      // for log10, sqrt, M_PI
#include <iostream>                   // for operator<<, ofstream, etc
#include <sstream>
#if ndim_ == 2
#include "VorCell.hpp"  // for VorCell
#include "VorFace.hpp"  // for VorFace
#include "VorGen.hpp"   // for VorGen
#include "VorTess.hpp"  // for VorTess
#endif
using namespace std;

// three possible colormaps: CM_JET, CM_AFMHOT, CM_OCEAN, CM_RDBU
// all are shamelessly stolen from matplotlib
// (as is the way in which they are calculated)
/*! @brief Matplotlib jet colormap: red: x-values */
double rxjet[5] = {0., 0.35, 0.66, 0.89, 1.};
/*! @brief Matplotlib jet colormap: red: y-values */
double ryjet[5] = {0., 0., 1., 1., 0.5};
/*! @brief Matplotlib jet colormap: green: x-values */
double gxjet[6] = {0., 0.125, 0.375, 0.64, 0.91, 1.};
/*! @brief Matplotlib jet colormap: green: y-values */
double gyjet[6] = {0., 0., 1., 1., 0., 0.};
/*! @brief Matplotlib jet colormap: blue: x-values */
double bxjet[5] = {0., 0.11, 0.34, 0.65, 1.};
/*! @brief Matplotlib jet colormap: blue: y-values */
double byjet[5] = {0.5, 1., 1., 0., 0.};
/*! @brief Matplotlib jet colormap definition */
static ColorMap _jet(rxjet, ryjet, 5, gxjet, gyjet, 6, bxjet, byjet, 5);

/*! @brief Matplotlib afmhot colormap: red: x-values */
double rxafmhot[3] = {0., 0.5, 1.};
/*! @brief Matplotlib afmhot colormap: red: y-values */
double ryafmhot[3] = {0., 1., 1.};
/*! @brief Matplotlib afmhot colormap: green: x-values */
double gxafmhot[4] = {0., 0.25, 0.75, 1.};
/*! @brief Matplotlib afmhot colormap: green: y-values */
double gyafmhot[4] = {0., 0., 1., 1.};
/*! @brief Matplotlib afmhot colormap: blue: x-values */
double bxafmhot[3] = {0., 0.5, 1.};
/*! @brief Matplotlib afmhot colormap: blue: y-values */
double byafmhot[3] = {0., 0., 1.};
/*! @brief Matplotlib afmhot colormap definition */
static ColorMap _afmhot(rxafmhot, ryafmhot, 3, gxafmhot, gyafmhot, 4, bxafmhot,
                        byafmhot, 3);

/*! @brief Matplotlib ocean colormap: red: x-values */
double rxocean[3] = {0., 0.66, 1.};
/*! @brief Matplotlib ocean colormap: red: y-values */
double ryocean[3] = {0., 0., 1.};
/*! @brief Matplotlib ocean colormap: green: x-values */
double gxocean[3] = {0., 0.33, 1.};
/*! @brief Matplotlib ocean colormap: green: y-values */
double gyocean[3] = {0.5, 0., 1.};
/*! @brief Matplotlib ocean colormap: blue: x-values */
double bxocean[2] = {0., 1.};
/*! @brief Matplotlib ocean colormap: blue: y-values */
double byocean[2] = {0., 1.};
/*! @brief Matplotlib ocean colormap definition */
static ColorMap _ocean(rxocean, ryocean, 3, gxocean, gyocean, 3, bxocean,
                       byocean, 2);

// Red to blue colormap:
/*! @brief Matplotlib rdbu colormap: red: x-values */
double rxrdbu[11] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
/*! @brief Matplotlib rdbu colormap: red: y-values */
double ryrdbu[11] = {0.4039, 0.6980, 0.8392, 0.9569, 0.9922, 0.9686,
                     0.8196, 0.5725, 0.2627, 0.1294, 0.0196};
/*! @brief Matplotlib rdbu colormap: green: x-values */
double gxrdbu[11] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
/*! @brief Matplotlib rdbu colormap: green: y-values */
double gyrdbu[11] = {0.,     0.09412, 0.3765, 0.6471, 0.8588, 0.9686,
                     0.8980, 0.7725,  0.5765, 0.4,    0.1882};
/*! @brief Matplotlib rdbu colormap: blue: x-values */
double bxrdbu[11] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
/*! @brief Matplotlib rdbu colormap: blue: y-values */
double byrdbu[11] = {0.1216, 0.1686, 0.3020, 0.5098, 0.7804, 0.9686,
                     0.9412, 0.8706, 0.7647, 0.6745, 0.3804};
/*! @brief Matplotlib rdbu colormap definition */
static ColorMap _rdbu(rxrdbu, ryrdbu, 11, gxrdbu, gyrdbu, 11, bxrdbu, byrdbu,
                      11);

// we want to be able to use the enum to determine the colormap, so we need an
// array that links integers to colormaps
/*! @brief List of default color maps */
static ColorMap* _default_colormaps[4] = {&_jet, &_afmhot, &_ocean, &_rdbu};

/**
  * @brief Constuct an image for a given Particle vector, plotting the given
  * variable
  *
  * The specifics of the grid of pixels used for the actual image can be
  * specified too.
  *
  * @param plist std::vector with Particle objects to plot
  * @param delcont Reference to a DelCont containing the particles, used to
  * determine the boundaries of the simulation volume
  * @param periodic Bool specifying whether the box containing the particles is
  * periodic (true) or reflective (false)
  * @param variable PlotVariable specifying which quantity to plot
  * @param grid Boolean specifying if the VorTess underlying the particles
  * should be overlayed on the plot
  * @param offset_x,offset_y Coordinates of the origin of the region to plot
  * @param height,width Size of the region to plot
  * @param pixsize Size of one pixel of the image in units of the simulation
  * length
  */
Image::Image(vector<Particle*>& plist, DelCont& delcont, bool periodic,
             PlotVariable variable, bool grid, double offset_x, double offset_y,
             double height, double width, double pixsize)
        : _delcont(delcont) {
    if(MPIGlobal::rank) {
        return;
    }
    unsigned long tottime = 1;
    tottime <<= 60;
    if(variable > 3) {
        //        _maxval = plist[0]->get_paq(variable-4);
        //        _minval = plist[0]->get_paq(variable-4);
        _maxval = 1. * plist[0]->get_timestep() / tottime;
        _minval = 1. * plist[0]->get_timestep() / tottime;
    } else {
        _maxval = ((GasParticle*)plist[0])->get_Wvec()[variable];
        _minval = ((GasParticle*)plist[0])->get_Wvec()[variable];
    }

    _periodic = periodic;
    calculate(plist, variable, grid, offset_x, offset_y, height, width,
              pixsize);
}

/**
  * @brief Calculate pixel values for all pixels in the actual image using a
  * Tree to find the closest cell
  *
  * The hydrodynamical quantities of the closest cell are then used to set the
  * pixel value.
  *
  * @param plist std::vector of Particle objects to plot
  * @param variable PlotVariable specifying which quantity to plot
  * @param grid Boolean specifying if the VorTess underlying the particles
  * should be overlayed on the plot
  * @param offset_x,offset_y Coordinates of the origin of the region to plot
  * @param height,width Size of the region to plot
  * @param pixsize Size of one pixel of the image in units of the simulation
  * length
  */
void Image::calculate(vector<Particle*>& plist, PlotVariable variable,
                      bool grid, double offset_x, double offset_y,
                      double height, double width, double pixsize) {
    unsigned long tottime = 1;
    tottime <<= 60;
    //#if ndim_==3
    //    Vec origin(0.5,0.5,0.5);
    //    Vec sides(1.,1.);
    //#else
    //    Vec origin(0.5,0.5);
    //    Vec sides(1.,1.);
    //#endif
    //    Cuboid box(origin, sides);
    Tree tree(_delcont.get_cuboid(), _periodic);
    for(unsigned int i = plist.size(); i--;) {
        plist[i]->set_key(_delcont.get_key(plist[i]->get_position()));
        tree.add_particle(plist[i]);
    }
    tree.finalize();
    //    ofstream treefile("tree");
    //    tree.print(treefile);
    double pixheight = height / pixsize;
    double pixwidth = width / pixsize;
    double cellheight = pixsize;
    double cellwidth = pixsize;
#if ndim_ == 3
    double sr = cbrt(0.75 / plist.size() / M_PI);
#else
    // the mean radius of a cell, good starting radius for search
    double sr = sqrt(1. / plist.size() / M_PI);
#endif
    StateVector maxW = ((GasParticle*)plist[0])->get_Wvec();
    StateVector minW = ((GasParticle*)plist[0])->get_Wvec();
    for(unsigned int i = 0; i < pixheight; i++) {
        vector<double> rowdens;
        for(unsigned int j = 0; j < pixwidth; j++) {
#if ndim_ == 3
            Vec center(offset_x + (j + 0.5) * cellwidth,
                       offset_y + (i + 0.5) * cellheight, 0.5);
#else
            Vec center(offset_x + (j + 0.5) * cellwidth,
                       offset_y + (i + 0.5) * cellheight);
#endif
            //            ClosestNgbSearch treewalk(center, sr);
            Particle* closest;
            if(_periodic) {
                closest = tree.get_periodic_closest(center, sr);
                //                while(!treewalk.get_closest()){
                //                    tree.walk_tree_periodic(treewalk);
                //                    treewalk.increase_radius();
                //                }
            } else {
                closest = tree.get_closest(center, sr);
                //                while(!treewalk.get_closest()){
                //                    tree.walk_tree(treewalk);
                //                    treewalk.increase_radius();
                //                }
            }
            //            Particle* closest = treewalk.get_closest();
            if(variable > 3) {
                rowdens.push_back(1. * closest->get_timestep() / tottime);
            } else {
                rowdens.push_back(
                        ((GasParticle*)closest)->get_Wvec()[variable]);
            }
            _maxval = max(_maxval, rowdens.back());
            _minval = min(_minval, rowdens.back());
            for(unsigned int maxindex = 0; maxindex < ndim_ + 2; maxindex++) {
                maxW[maxindex] =
                        std::max(maxW[maxindex],
                                 ((GasParticle*)closest)->get_Wvec()[maxindex]);
                minW[maxindex] =
                        std::min(minW[maxindex],
                                 ((GasParticle*)closest)->get_Wvec()[maxindex]);
            }
        }
        _image.push_back(rowdens);
    }
    if(grid) {
#if ndim_ == 3
        cout << "Plotting grid does not work (yet) for 3D" << endl;
#else
        VorTess voronoi(&_delcont, plist.size(), _periodic);
        for(unsigned int i = plist.size(); i--;) {
            voronoi.add_point((GasParticle*)plist[i], i);
        }
        voronoi.complete(tree);
        voronoi.construct();
        //        int type = BIN_PPM | CM_JET;
        for(unsigned int i = plist.size(); i--;) {
            VorCell* cell = ((GasParticle*)plist[i])->get_cell();
            vector<VorFace*> faces = cell->get_faces();
            for(unsigned int j = faces.size(); j--;) {
                if(faces[j]->get_area()) {
                    vector<VorGen*> vertices = faces[j]->get_vertices();
                    double d[4] = {vertices[0]->x(), vertices[0]->y(),
                                   vertices[1]->x(), vertices[1]->y()};
                    d[0] = (d[0] - offset_x) / width;
                    d[2] = (d[2] - offset_x) / width;
                    d[1] = (d[1] - offset_y) / height;
                    d[3] = (d[3] - offset_y) / height;
                    draw_line(d[0] * pixwidth, d[1] * pixheight,
                              d[2] * pixwidth, d[3] * pixheight);
                }
            }
        }
#endif
    }
    cout << "[max,min]: " << _maxval << ", " << _minval << endl;
}

/**
  * @brief Draw a line between the given two points in the pixel grid
  *
  * The pixel values for pixels on the line are set to negative values. The line
  * starts at the coordinates of the first point and stops at the coordinates of
  * the second point.
  *
  * @param x0,y0 Coordinates of the first point
  * @param x1,y1 Coordinates of the second point
  */
void Image::draw_line(int x0, int y0, int x1, int y1) {
    for(double par = 0.; par <= 1.; par += 0.001) {
        unsigned int pix[2] = {
                y0 + static_cast<unsigned int>(par * (y1 - y0) + 0.5),
                x0 + static_cast<unsigned int>(par * (x1 - x0) + 0.5)};
        if(pix[0] < _image[0].size() && pix[1] < _image.size()) {
            _image[pix[0]][pix[1]] = -1.;
        }
    }
}

/**
  * @brief Save the plot to a file with given name
  *
  * Images can be saved in binary or in ascii format, color or grayscale.
  * Minimal and maximal cutoff values can be specified and logarithmic scaling
  * can be applied. The images are saved in the Netpbm format
  * (http://en.wikipedia.org/wiki/Netpbm_format).
  *
  * @param name Name of the image file
  * @param type ImageType specifying which file format and color map to use
  * @param maxval,minval Maximal and minimal cutoff values for the plot
  * @param logarithmic If true, logarithmic scaling is used instead of standard
  * linear scaling
  */
void Image::save(std::string name, int type, double maxval, double minval,
                 bool logarithmic) {
    if(MPIGlobal::rank) {
        return;
    }
    if(maxval == -1.) {
        maxval = _maxval;
    }
    if(minval == -1.) {
        minval = _minval;
    }
    if(maxval == minval) {
        minval -= 0.01;
        maxval += 0.01;
    }
    if(!_image.size()) {
        cerr << "Image does not have a size!" << endl;
        my_exit();
    }
    ofstream file;
    stringstream filename;
    filename << name;
    if(type < 2) {
        filename << ".pgm";
        if(type & 1) {
            // ASCII image
            file.open(filename.str().c_str());
            file << "P2\n";
        } else {
            // binary image
            file.open(filename.str().c_str(), ios::out | ios::binary);
            file << "P5\n";
        }
        file << _image[0].size() << "\t" << _image.size() << "\n255\n";
        // the loop determines the direction of the image, so don't simply
        // change it to boost performance!
        for(unsigned int i = _image.size(); i--;) {
            for(unsigned int j = 0; j < _image[0].size(); j++) {
                double value;
                if(logarithmic) {
                    value = log10(_image[i][j] / minval) /
                            log10(maxval / minval);
                } else {
                    value = (_image[i][j] - minval) / (maxval - minval);
                }
                int grayscale = int(min(1., value) * 255);
                if(type & 1) {
                    if(!j) {
                        file << grayscale << "\n";
                    } else {
                        file << grayscale << "\t";
                    }
                } else {
                    file.write((char*)&grayscale, 1);
                }
            }
        }
        file.close();
    } else {
        if(!type & 2) {
            cout << "Error: can't use a colormap for a grayscale image!"
                 << endl;
            my_exit();
        }
        filename << ".ppm";
        if(type & 1) {
            // ASCII image
            file.open(filename.str().c_str());
            file << "P3\n";
        } else {
            // binary image
            file.open(filename.str().c_str(), ios::out | ios::binary);
            file << "P6\n";
        }
        file << _image[0].size() << "\t" << _image.size() << "\n255\n";
        // the loop determines the direction of the image, so don't simply
        // change it to boost performance!
        for(unsigned int i = _image.size(); i--;) {
            for(unsigned int j = 0; j < _image[0].size(); j++) {
                int colors[3] = {0};
                if(_image[i][j] < _minval) {
                    colors[0] = 0;
                    colors[1] = 0;
                    colors[2] = 0;
                } else {
                    double value;
                    if(logarithmic) {
                        value = log10(_image[i][j] / minval) /
                                log10(maxval / minval);
                    } else {
                        value = (_image[i][j] - minval) / (maxval - minval);
                    }
                    // make sure all values are in the right interval
                    value = min(value, 1.);
                    value = max(value, 0.);
                    _default_colormaps[(type >> 2) - 1]->get_color(
                            min(value, 1.), colors);
                }
                if(type & 1) {
                    if(!j) {
                        file << colors[0] << " " << colors[1] << " "
                             << colors[2] << "\n";
                    } else {
                        file << colors[0] << " " << colors[1] << " "
                             << colors[2] << "\t";
                    }
                } else {
                    file.write((char*)&colors[0], 1);
                    file.write((char*)&colors[1], 1);
                    file.write((char*)&colors[2], 1);
                }
            }
        }
        file.close();
    }
}
