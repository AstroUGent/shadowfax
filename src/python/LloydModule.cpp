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
 * @file LloydModule.cpp
 *
 * @brief Auxiliary program to compute the centroids of the Voronoi mesh:
 * Python module
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Lloyd.hpp"
#include <cstddef>   // for NULL
#include <iostream>  // for operator<<, basic_ostream, etc
#include <vector>
using namespace std;

#include <boost/python/def.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/python/module.hpp>
using namespace boost::python;

/**
 * @brief List wrapper around Lloyd::calculate_centroids
 *
 * @param coords Coordinates of the points
 * @param origin Origin of the box
 * @param sides Side lengths of the box
 * @param periodic Flag indicating if the box is periodic or not
 * @return vector result from Lloyd::calculate_centroids
 */
list python_centroids(list coords, list origin, list sides, bool periodic) {
    // suppress output to cout
    streambuf* backout = cout.rdbuf();
    cout.rdbuf(NULL);

    vector<double> veccoords;
    vector<double> vecorigin;
    vector<double> vecsides;
    for(int i = 0; i < len(coords); i++) {
        veccoords.push_back(extract<double>(coords[i]));
    }
    for(int i = 0; i < len(origin); i++) {
        vecorigin.push_back(extract<double>(origin[i]));
    }
    for(int i = 0; i < len(sides); i++) {
        vecsides.push_back(extract<double>(sides[i]));
    }
    vector<double> centroids = Lloyd::calculate_centroids(veccoords, vecorigin,
                                                          vecsides, periodic);
    list l;
    for(unsigned int i = 0; i < centroids.size(); i++) {
        l.append(centroids[i]);
    }

    // reinstate cout output
    cout.rdbuf(backout);

    return l;
}

#if ndim_ == 3
/*! @brief Register the liblloyd Python module */
BOOST_PYTHON_MODULE(libpython_lloyd3d)
#else
BOOST_PYTHON_MODULE(libpython_lloyd2d)
#endif
{
    def("calculate_centroids", &python_centroids);
}
