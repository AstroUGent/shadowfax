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
 * @file Lloyd.hpp
 *
 * @brief Auxiliary program to compute the centroids of the Voronoi mesh: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef LLOYD_HPP
#define LLOYD_HPP

#include <vector>

/**
 * @brief Auxiliary program to calculate the centroids of the Voronoi mesh
 *
 * By moving the mesh generators to the respective centroids, we perform one
 * iteration of Lloyd's algorithm.
 */
class Lloyd {
  public:
    Lloyd(int argc, char** argv);

    static std::vector<double> calculate_centroids(
            std::vector<double>& coords, std::vector<double>& box_origin,
            std::vector<double>& box_sides, bool periodic);
};

#endif  // LLOYD_HPP
