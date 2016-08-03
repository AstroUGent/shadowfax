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
 * @file ColorMap.hpp
 *
 * @brief Color Map: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef COLORMAP_HPP
#define COLORMAP_HPP

#include "ImageOptions.hpp"
#include <string>
#include <vector>

/**
  * @brief Color map used to convert numerical values to pixel colors
  *
  * The ColorMap consists of a set of 3 color functions that map numerical
  * values in the range [0,1] to the RGB pixel values. The color functions
  * themselves are sampled using a set of points and are evaluated using linear
  * interpolation between this points. The method used here (as well as the
  * pre-defined color maps) are copied from the corresponding Matplotlib code
  * (https://github.com/matplotlib/matplotlib).
  */
class ColorMap {
  private:
    double interpolate(double value, std::vector<double>& cx,
                       std::vector<double>& cy);

    /*! @brief x-components of the red color function */
    std::vector<double> _rx;
    /*! @brief y-components of the red color function */
    std::vector<double> _ry;
    /*! @brief x-components of the green color function */
    std::vector<double> _gx;
    /*! @brief y-components of the green color function */
    std::vector<double> _gy;
    /*! @brief x-components of the blue color function */
    std::vector<double> _bx;
    /*! @brief y-components of the blue color function */
    std::vector<double> _by;

  public:
    ColorMap(double* rx, double* ry, unsigned int nr, double* gx, double* gy,
             unsigned int ng, double* bx, double* by, unsigned int nb);
    ~ColorMap() {}

    void save_map(std::string filename, int type = BIN_PPM);
    void get_color(double value, int* colors);
};

#endif  // COLORMAP_HPP
