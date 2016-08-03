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
 * @file ImageDescription.hpp
 *
 * @brief Parameter file for image generation: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef IMAGEDESCRIPTION_HPP
#define IMAGEDESCRIPTION_HPP

#include <string>

/**
 * @brief Representation of a parameter file for Image creation
 *
 * The parameter file is a .ini file. Every parameter has a default value, which
 * is used if no corresponding parameter is specified in the file.
 *
 * A complete parameter file would look like this:
\verbatim
[Variable]
; Name of the variable to be plotted (density/velocity_x/velocity_y/pressure,
; default: density)
Name = density

[Range]
; Maximal value that can be represented on the color scale or -1 to use the
; maximal value of all cells (default: -1)
Max = -1
; Minimal value that can be represented on the color scale or -1 to use the
; minimal value of all cells (default: -1)
Min = -1
; Use a logarithmic scale? (default: false)
Log = false

[Scale]
; Use a color scale? (default: true)
Color = true

[Grid]
; Plot the Voronoi grid? (default: false)
Active = false

[Offset]
; x-coordinate of the offset of the image within the simulation box or -1 to use
; the origin of the box (default: -1)
X = -1
; y-coordinate of the offset of the image within the simulation box or -1 to use
; the origin of the box (default: -1)
Y = -1

[Dimensions]
; Size of the image in the x-direction, or -1 to plot everything from the offset
; to the boundary of the box (default: -1)
X = -1
; Size of the image in the y-direction, or -1 to plot everything from the offset
; to the boundary of the box (default: -1)
Y = -1

[Pixel]
; Size of a single pixel of the image or -1 to get an image with 1000 pixels in
; the direction with the most pixels (default: -1)
Size = -1
\endverbatim
 */
class ImageDescription {
  private:
    /*! @brief Name of variable to be plotted */
    std::string _variable_name;

    /*! @brief Maximal value of the plot interval */
    double _range_max;

    /*! @brief Minimal value of the plot interval */
    double _range_min;

    /*! @brief Use a logarithmic scale? */
    bool _log;

    /*! @brief Plot the Voronoi grid? */
    bool _grid;

    /*! @brief Color or grayscale? */
    bool _color;

    /*! @brief Offset of the image in the simulation box */
    double _offset[2];

    /*! @brief Dimensions of the image in physical length */
    double _size[2];

    /*! @brief Physical side length of a single pixel */
    double _pixelsize;

    void print_contents();

  public:
    ImageDescription(std::string filename = "");

    std::string get_variable_name();

    double get_range_max();
    double get_range_min();

    bool is_log();
    bool plot_grid();
    bool is_color();

    double get_offset_x();
    double get_offset_y();

    double get_size_x();
    double get_size_y();

    double get_pixelsize();
};

#endif  // IMAGEDESCRIPTION_HPP
