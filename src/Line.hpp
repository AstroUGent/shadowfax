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
 * @file Line.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief Representation of a line (segment): header
 *
 * contains class Line
 *
 * Class to represent a line. A line can be constructed by either indicating two
 * endpoints or by specifying a position and a direction.
 * This class holds some geometrical operations involving lines.
 */
#ifndef HEAD_LINE
#define HEAD_LINE

#include <ostream>  // for ostream
#include <vector>   // for vector

class VorGen;
#if ndim_ == 3
class Plane;
#endif

/**
 * @brief Representation of a geometrical line in 2D or 3D
 *
 * Currently used to convert two locations to the mathematical representation
 * of the line connecting them.
 */
class Line {
  private:
    /*! \brief x-coordinate of the first endpoint */
    double _x1;
    /*! \brief x-coordinate of the second endpoint */
    double _x2;

    /*! \brief y-coordinate of the first endpoint */
    double _y1;
    /*! \brief y-coordinate of the second endpoint */
    double _y2;

    /*! \brief x-component of the direction vector of the line */
    double _xdir;
    /*! \brief y-component of the direction vector of the line */
    double _ydir;
#if ndim_ == 3
    /*! \brief z-coordinate of the first endpoint */
    double _z1;
    /*! \brief z-coordinate of the second endpoint */
    double _z2;

    /*! \brief z-component of the direction vector of the line */
    double _zdir;
#endif

    std::vector<double> get_midpoint();

  public:
    Line(VorGen* point1, VorGen* point2);
#if ndim_ == 3
    Line(double x, double y, double z, double xdir, double ydir, double zdir);
#else
    Line(double x, double y, double xdir, double ydir);
#endif
    ~Line() {}
    double x();
    double xdir();
    double y();
    double ydir();
    std::vector<double> get_coordinates(char coord, double value);
#if ndim_ == 3
    double z();
    double zdir();
#endif
#if ndim_ == 3
    Plane get_bisector();
#endif
    void print(std::ostream& stream);
};

#endif
