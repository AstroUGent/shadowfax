/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ColorMap.cpp
 *
 * @brief Color Map: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ColorMap.hpp"
#include <fstream>
#include <sstream>  // for operator<<, etc
using namespace std;

/**
  * @brief Construct a ColorMap with the given color functions
  *
  * The color functions are specified by linear chunks in the 3 RGB colors.
  *
  * @param rx,ry Red color chunks
  * @param nr Number of red color chunks
  * @param gx,gy,ng Green color chunks
  * @param bx,by,nb Blue color chunks
  */
ColorMap::ColorMap(double* rx, double* ry, unsigned int nr, double* gx,
                   double* gy, unsigned int ng, double* bx, double* by,
                   unsigned int nb) {
    for(unsigned int i = 0; i < nr; i++) {
        _rx.push_back(rx[i]);
        _ry.push_back(ry[i]);
    }
    for(unsigned int i = 0; i < ng; i++) {
        _gx.push_back(gx[i]);
        _gy.push_back(gy[i]);
    }
    for(unsigned int i = 0; i < nb; i++) {
        _bx.push_back(bx[i]);
        _by.push_back(by[i]);
    }
}

/**
  * @brief Convert a double value in the range [0,1] to RGB colors in the range
  * [0,255]
  *
  * @param value Numerical value to convert from
  * @param colors 3-element array to store the resulting colors in
  */
void ColorMap::get_color(double value, int* colors) {
    colors[0] = int(interpolate(value, _rx, _ry) * 255 + 0.5);
    colors[1] = int(interpolate(value, _gx, _gy) * 255 + 0.5);
    colors[2] = int(interpolate(value, _bx, _by) * 255 + 0.5);
}

/**
  * @brief Linearly interpolate the given value on the given color function
  *
  * @param value A seed value
  * @param cx,cy Color function
  * @return A value in the range [0,1]
  */
double ColorMap::interpolate(double value, vector<double>& cx,
                             vector<double>& cy) {
    unsigned int i = 1;
    while(value > cx[i]) {
        i++;
    }
    if(value == cx[i]) {
        return cy[i];
    } else {
        return cy[i - 1] +
               (cy[i] - cy[i - 1]) / (cx[i] - cx[i - 1]) * (value - cx[i - 1]);
    }
}

/**
  * @brief Save a sample of the ColorMap to a file with given name and type
  *
  * @param filename Name of the file
  * @param type ImageType specifying color map and file type of the image
  */
void ColorMap::save_map(string name, int type) {
    ofstream file;
    stringstream filename;
    filename << name;
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
    file << 1000 << "\t" << 400 << "\n255\n";
    // the loop determines the direction of the image, so don't simply change it
    // to boost performance!
    for(unsigned int i = 0; i < 400; i++) {
        for(unsigned int j = 1000; j--;) {
            int colors[3] = {0};
            get_color(j * 0.001, colors);
            if(type & 1) {
                if(!j) {
                    file << colors[0] << " " << colors[1] << " " << colors[2]
                         << "\n";
                } else {
                    file << colors[0] << " " << colors[1] << " " << colors[2]
                         << "\t";
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
