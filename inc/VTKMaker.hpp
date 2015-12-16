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
 * @file VTKMaker.hpp
 *
 * @brief Sideprogram to convert snapshots to VTK-files: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef VTKMAKER_HPP
#define VTKMAKER_HPP

#include <ostream>

/**
 * @brief The vtkmaker program
 *
 * This program is used to convert a given snapshot to a .vtk-file that can
 * be used to visualize its Voronoi tesselation with e.g. VisIt.
 *
 * The .vtk-file has the same name as the input snapshot.
 */
class VTKMaker{
private:
    /**
     * \brief Flag indicating if the system stores floating point values in big
     *  endian or little endian order
     *
     * Big endian systems store the most significant byte first (lowest memory
     * address), while little endian systems do it the other way around. This
     * means that a 4-byte integer will be stored in reversed order on little
     * endian systems. VTK binary files use big endian ordering, so we have to
     * invert all 4-byte values on little endian systems.
     */
    bool _big_endian;

    void write_big_endian_buffer(std::ostream& stream, float* buffer,
                                 unsigned int length);
    void write_big_endian_buffer(std::ostream& stream, int* buffer,
                                 unsigned int length);
    void write_big_endian_buffer(std::ostream& stream, char* buffer,
                                 unsigned int length);
    void check_big_endian();
    void swap_endian(char* buffer, unsigned int length);

public:
    VTKMaker(int argc, char** argv);
};

#endif
