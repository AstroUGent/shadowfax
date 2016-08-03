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
 * @file Image.hpp
 *
 * @brief PPM image: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef IMAGE_HPP
#define IMAGE_HPP

#include "ImageOptions.hpp"
#include <string>
#include <vector>

class Particle;
class DelCont;

/**
  * @brief Variables that can be plotted
  */
enum PlotVariable {
    /*! @brief The hydrodynamical density */
    HYDRO_DENSITY,
    /*! @brief The x-component of the hydrodynamical velocity */
    HYDRO_VELOCITY_X,
    /*! @brief The y-component of the hydrodynamical velocity */
    HYDRO_VELOCITY_Y,
    /*! @brief The hydrodynamical pressure */
    HYDRO_PRESSURE,
    /*! @brief The Passively Advected Quantity */
    PAQ
};

/**
  * @brief Plot of a hydrodynamical grid variable using a cartesian grid of
  * pixels
  *
  * This class is responsible for the calculation of the grid of pixels and for
  * the saving of the resulting image to a file using the Netpbm format
  * (http://en.wikipedia.org/wiki/Netpbm_format)
  */
class Image {
  private:
    void calculate(std::vector<Particle*>& plist, PlotVariable variable,
                   bool grid, double offset_x, double offset_y, double height,
                   double width, double pixsize);

    /*! @brief The grid of pixels that constitutes the actual image */
    std::vector<std::vector<double> > _image;
    /*! @brief The maximal value of the plotting variable over all pixels */
    double _maxval;
    /*! @brief The minimal value of the plotting variable over all pixels */
    double _minval;
    /*! @brief Reference to the DelCont containing the simulation volume */
    DelCont& _delcont;

    /*! @brief Flag indicating if the simulation box is periodic or
     *  reflective */
    bool _periodic;

  public:
    Image(std::vector<Particle*>& plist, DelCont& delcont,
          bool periodic = false, PlotVariable variable = HYDRO_DENSITY,
          bool grid = false, double offset_x = 0., double offset_y = 0.,
          double height = 1., double width = 1., double pixsize = 0.001);
    ~Image() {}

    void draw_line(int x0, int y0, int x1, int y1);
    void save(std::string name, int type = BIN_PPM | CM_RDBU,
              double maxval = -1., double minval = -1.,
              bool logarithmic = false);
};

#endif  // IMAGE_HPP
