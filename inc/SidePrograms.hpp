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
 * @file SidePrograms.hpp
 *
 * @brief A collection of programs that need the Voronoi tesselation but are not
 * really related to the main simulation program: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SIDEPROGRAMS_HPP
#define SIDEPROGRAMS_HPP

#include <string>
#include <vector>

class ParticleVector;

/**
 * @brief General interface for SidePrograms
 *
 * Contains a helper function to load ASCII-files
 */
class SidePrograms {
  public:
    SidePrograms() {}
    ~SidePrograms() {}
    void load_ascii(ParticleVector& cells, double* cube_boundaries,
                    std::string name);
};

/**
 * @brief SideProgram to calculate areas for a file containing coordinates in 2D
 *
 * The program reads coordinates from a given ASCII-file and constructs the
 * Voronoi mesh. It then outputs the areas of the Voronoi cells to the standard
 * output, in the same order as the input coordinates.
 * Only works in 2D.
 */
class AreaCalculator : public SidePrograms {
  public:
    AreaCalculator(int argc, char** argv);
    ~AreaCalculator() {}
};

/**
 * @brief SideProgram to calculate densities for a file containing coordinates
 * in 3D
 *
 * The program reads coordinates from a given ASCII-file and constructs the
 * Voronoi mesh. It then calculates a density for every coordinate triple by
 * taking the inverse of the associated Voronoi cell. These densities, together
 * with the original coordinates, are stored to a Shadowfax snapshot with the
 * name densities.hdf5.
 * The program also outputs the coordinates and line number corresponding to
 * the densest coordinates to the standard output.
 * Only works for 3D.
 */
class DensityCalculator : public SidePrograms {
  public:
    DensityCalculator(int argc, char** argv);
    ~DensityCalculator() {}
};

/**
 * @brief SideProgram to sort a file with particle coordinates in Hilbert order
 *
 * The program takes the name of an ASCII-file as an argument and writes a new
 * ASCII-file containing the coordinates from the input file, but then sorted
 * using the Hilbert key of the input coordinates. Only works for 3D.
 * The output file has the same name as the input file, but with _sorted added
 * to it.
 */
class HilbertSorter : public SidePrograms {
  public:
    HilbertSorter(int argc, char** argv);
};

/**
 * @brief SideProgram to calculate masses for a given Shadowfax snapshot
 *
 * The program reads the coordinates and densities from the snapshot, creates
 * a Voronoi mesh and uses its cell volumes to convert the densities to masses,
 * which are then written to an ASCII-file. The output file has the same name as
 * the input file, but with the extension .mass.
 */
class MassCalculator : public SidePrograms {
  public:
    MassCalculator(int argc, char** argv);
};

/**
 * @brief SideProgram that calculates the gravitational potential energy for a
 * given snapshot
 *
 * The Voronoi tesselation of the gasparticles is calculated and used to convert
 * densities to masses. These masses are then used to calculate the
 * gravitational potential using a Barnes-Hut treewalk.
 * The gravitational potential energies are written to an ASCII-file with the
 * same name as the given snapshot, but with the extension .epot.
 */
class PotentialCalculator : public SidePrograms {
  public:
    PotentialCalculator(int argc, char** argv);
};

#endif  // SIDEPROGRAMS_HPP
