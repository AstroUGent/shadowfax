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
 * @file ImageDescription.cpp
 *
 * @brief Parameter file for image generation: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ImageDescription.hpp"
#include "Error.hpp"                           // for my_exit
#include <boost/property_tree/ini_parser.hpp>  // for read_ini
#include <iostream>                            // for cout, cerr

/**
 * @brief Print the contents of the parameter file to the stdout
 */
void ImageDescription::print_contents() {
    std::cout << "[Variable]" << std::endl;
    std::cout << "Name: " << _variable_name << std::endl;
    std::cout << std::endl;

    std::cout << "[Range]" << std::endl;
    std::cout << "Max: " << _range_max << std::endl;
    std::cout << "Min: " << _range_min << std::endl;
    std::cout << "Log: ";
    if(_log) {
        std::cout << "true" << std::endl;
    } else {
        std::cout << "false" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "[Scale]" << std::endl;
    std::cout << "Color: ";
    if(_color) {
        std::cout << "true" << std::endl;
    } else {
        std::cout << "false" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "[Grid]" << std::endl;
    std::cout << "Active: ";
    if(_grid) {
        std::cout << "true" << std::endl;
    } else {
        std::cout << "false" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "[Offset]" << std::endl;
    std::cout << "X: " << _offset[0] << std::endl;
    std::cout << "Y: " << _offset[1] << std::endl;
    std::cout << std::endl;

    std::cout << "[Dimensions]" << std::endl;
    std::cout << "X: " << _size[0] << std::endl;
    std::cout << "Y: " << _size[1] << std::endl;
    std::cout << std::endl;

    std::cout << "[Pixel]" << std::endl;
    std::cout << "Size: " << _pixelsize << std::endl;
    std::cout << std::endl;
}

/**
 * @brief Constructor
 *
 * Read the parameter file (a .ini file) and initialize the internal variables
 * using its contents.
 *
 * @param filename Name of the .ini file to read
 */
ImageDescription::ImageDescription(std::string filename) {
    if(filename.empty()) {
        // use default values
        _variable_name = "density";

        _range_max = -1.;
        _range_min = -1.;
        _log = false;

        _color = true;

        _grid = false;

        _offset[0] = -1.;
        _offset[1] = -1.;

        _size[0] = -1.;
        _size[1] = -1.;

        _pixelsize = -1.;
    } else {
        std::ifstream file(filename.c_str());
        if(!file) {
            std::cerr << "Cannot read parameterfile \"" << filename << "\"!"
                      << std::endl;
            my_exit();
        }

        boost::property_tree::ptree pt;
        boost::property_tree::ini_parser::read_ini(filename, pt);

        _variable_name = pt.get<std::string>("Variable.Name", "density");

        _range_max = pt.get<double>("Range.Max", -1.);
        _range_min = pt.get<double>("Range.Min", -1.);
        _log = pt.get<bool>("Range.Log", false);

        _color = pt.get<bool>("Scale.Color", true);

        _grid = pt.get<bool>("Grid.Active", false);

        _offset[0] = pt.get<double>("Offset.X", -1.);
        _offset[1] = pt.get<double>("Offset.Y", -1.);

        _size[0] = pt.get<double>("Dimensions.X", -1.);
        _size[1] = pt.get<double>("Dimensions.Y", -1.);

        _pixelsize = pt.get<double>("Pixel.Size", -1.);
    }

    print_contents();
}

/**
 * @brief Get the name of the variable to be plotted in the Image
 *
 * @return Variable name
 */
std::string ImageDescription::get_variable_name() {
    return _variable_name;
}

/**
 * @brief Get the maximal value of the variable that can be represented on the
 * color scale
 *
 * @return Maximal value of the variable to be plotted
 */
double ImageDescription::get_range_max() {
    return _range_max;
}

/**
 * @brief Get the minimal value of the variable that can be represented on the
 * color scale
 *
 * @return Minimal value of the variable to be plotted
 */
double ImageDescription::get_range_min() {
    return _range_min;
}

/**
 * @brief Use a logarithmic scale?
 *
 * @return True if a logscale is used, false for a linear scale
 */
bool ImageDescription::is_log() {
    return _log;
}

/**
 * @brief Plot the Voronoi grid?
 *
 * @return True if the grid should be plotted, false if only the contents of the
 * cell should be plotted and no cell boundaries
 */
bool ImageDescription::plot_grid() {
    return _grid;
}

/**
 * @brief Use color or grayscale?
 *
 * @return True for a color image, false for a grayscale image
 */
bool ImageDescription::is_color() {
    return _color;
}

/**
 * @brief Get the x-coordinate of the offset of the image within the simulation
 * box
 *
 * @return x-coordinate of image offset
 */
double ImageDescription::get_offset_x() {
    return _offset[0];
}

/**
 * @brief Get the y-coordinate of the offset of the image within the simulation
 * box
 *
 * @return y-coordinate of image offset
 */
double ImageDescription::get_offset_y() {
    return _offset[1];
}

/**
 * @brief Get the size of the image (in physical length) in the x-direction
 *
 * @return Size of image in x-direction
 */
double ImageDescription::get_size_x() {
    return _size[0];
}

/**
 * @brief Get the size of the image (in physical length) in the y-direction
 *
 * @return Size of image in y-direction
 */
double ImageDescription::get_size_y() {
    return _size[1];
}

/**
 * @brief Get the side length of a single (square) pixel in physical length
 *
 * @return Size of a pixel
 */
double ImageDescription::get_pixelsize() {
    return _pixelsize;
}
