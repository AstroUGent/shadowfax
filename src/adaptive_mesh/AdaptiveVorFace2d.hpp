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
 * @file AdaptiveVorFace2d.hpp
 *
 * @brief 2D Voronoi face, used for geometrical face operations: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_ADAPTIVEVORFACE2D
#define HEAD_ADAPTIVEVORFACE2D

#include <ostream>

class StateVector;
class RestartFile;

/**
 * @brief 3D adaptive Voronoi face
 *
 * Used for the 3D mesh evolution algorithm.
 */
class AdaptiveVorFace2d {
  private:
    /*! @brief Velocity of the face */
    double _v[2];

    /*! @brief Positions of the vertices of the face */
    double _pos[4];

    /*! @brief Position of the midpoint of the face */
    double _mid[2];

    /*! @brief Flag to indicate whether this face is complete or not (probably
     *  no longer used) */
    bool _full;

  public:
    AdaptiveVorFace2d(double* pos);

    void deactivate();
    bool active();
    void set_v(double* v);
    bool add_Lvertex(double* a, double* b, double* c, double* vert);
    void set_Rvertex(double* vert);
    void print(std::ostream& stream);
    void get_points(double* L, double* R);
    void set_points(double* L, double* R);

    //    void save_restart(std::ostream& stream);
    //    AdaptiveVorFace2d(std::istream& stream);

    double get_area();
    void get_midpoint(double* midpoint);
    void transform(StateVector& W, double* left, double* right);
    void invtransform(StateVector& W, double* left, double* right);

    void dump(RestartFile& rfile);
    AdaptiveVorFace2d(RestartFile& rfile);
};

#endif  // HEAD_ADAPTIVEVORFACE2D
