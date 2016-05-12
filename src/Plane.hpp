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
 * @file Plane.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief 3D plane: header
 *
 * contains class Plane
 *
 * A Plane represents a plane. It can be constructed by providing 3 points on
 * the plane or by specifying a point and a normal to the plane. This class
 * holds some geometrical operations involving planes.
 */
#ifndef HEAD_PLANE
#define HEAD_PLANE

#if ndim_ == 3

class VorGen;
class Line;

/**
 * @brief Representation of a geometrical plane in 3D
 */
class Plane {
  private:
    /*! \brief x-component of the normal vector of the plane */
    double _nx;
    /*! \brief y-component of the normal vector of the plane */
    double _ny;
    /*! \brief z-component of the normal vector of the plane */
    double _nz;

    /*! \brief x-coordinate of a reference point on the plane */
    double _x;
    /*! \brief y-coordinate of a reference point on the plane */
    double _y;
    /*! \brief z-coordinate of a reference point on the plane */
    double _z;

  public:
    Plane(double nx, double ny, double nz, double x, double y, double z);
    Plane(VorGen* point1, VorGen* point2, VorGen* point3);
    ~Plane(){};
    VorGen intersect(Line* line);
    bool intersect(Plane* plane);
    VorGen intersect_planes(Plane* plane1, Plane* plane2);
    double distance(VorGen* point);
    double distance(Plane* plane);
};

#endif

#endif
