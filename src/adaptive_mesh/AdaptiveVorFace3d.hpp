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
 * @file AdaptiveVorFace3d.hpp
 *
 * @brief 3D Voronoi face, used for geometrical face operations: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_ADAPTIVEVORFACE3D
#define HEAD_ADAPTIVEVORFACE3D

#include <ostream>  // for ostream
#include <vector>   // for vector

class Cuboid;
class RestartFile;
class StateVector;

/**
 * @brief 3D adaptive Voronoi face
 *
 * Used to calculate areas and centroids for the 3D mesh evolution algorithm.
 */
class AdaptiveVorFace3d {
  private:
    /*! @brief Positions of the vertices */
    std::vector<double> _pos;

  public:
    AdaptiveVorFace3d(RestartFile& rfile);
    void dump(RestartFile& rfile);

    AdaptiveVorFace3d();

    void add_vertexpoint(double* vert);
    void print(std::ostream& stream);
    void print_copy(std::ostream& stream, int wall, Cuboid& box);

    std::vector<double> get_vertices();

    double calculate_quantities(double* midpoint);
    void get_midpoint(double* midpoint);
    double get_area();

    void transform(StateVector& W, double* left, double* right);
    void invtransform(StateVector& W, double* left, double* right);
};

#endif
