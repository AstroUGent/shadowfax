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
 * @file DelCont.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief Containers for the Delaunay tesselation: general header
 *
 * contains interface Container, class CubicBox, class CircularBox
 *
 * The Container interface defines an object that holds information on the
 * boundaries of an area/volume. The most obvious functions are determining
 * whether a given point/position is inside the area/volume. Furthermore, the
 * Container gives information about the dimensions of the area/volume and
 * gives a prescription to determine mirror points for its boundaries.
 * A Container can not be constructed directly, but has to be constructed by
 * using a child class.
 *
 * The CubicBox represents a cubic (or square) box with given origin and side.
 *
 * The CircularBox represents a round or spheric volume with given central point
 * and radius.
 */
#ifndef HEAD_DELCONT
#define HEAD_DELCONT

#include "utilities/Cuboid.hpp"
#include <vector>

class Particle;
class Vec;
class VorGen;

/**
 * @brief General interface to represent volumes
 *
 * A Container practically defines a volume (box or sphere), by making the
 * geometrical properties of the volume accessible to the rest of the program.
 */
class DelCont {
  public:
    virtual ~DelCont() {}

    /**
      * Assess whether the given VorGen and the region within a sphere with the
      * given radius around this VorGen are inside the DelCont
      *
      * @param point A VorGen specifying the center of the search region
      * @param radius The radius specifying the size of the search region
      * @return true if the region is completely inside the DelCont (borders
      * exclusive), false otherwise
      */
    virtual bool inside(VorGen* point, double radius) = 0;

    /**
      * Assess whether the given Vec is inside the DelCont
      *
      * @param pntpos A Vec specifying the coordinates of a VorGen or Particle
      * @return true if the coordinates are inside the DelCont (borders
      * inclusive), false otherwises
      */
    virtual bool inside(Vec pntpos) = 0;

    /**
      * Get the 4 (or 3) coordinate triplets (/pairs) that make up the vertices
      * of a simplex that contains the entire DelCont volume and ghost regions,
      * to assure that all points that will be added during the construction of
      * the DelTess for the points in the DelCont are inside this simplex.
      *
      * @return A vector containing the coordinates of a Simplex
      */
    virtual std::vector<double> get_bounding_tetrahedron() = 0;

    /**
      * Return the width of the DelCont. For a non-cubic DelCont, this is the
      * side of a cube that contains the entire volume.
      *
      * @return The width of the DelCont
      */
    virtual double get_box_width() = 0;

    /**
      * Get the origin and side of a cubic box containing the entire DelCont
      *
      * @param box A 3- or 4-element array to fill.
      */
    virtual void get_bounding_box(double* box) = 0;

    /**
      * Calculate a hilbert key for the given Vec and taking into account the
      * size of the DelCont
      *
      * @param coords A Vec specifying the coordinates of a Particle or VorGen
      * @return The hilbert key for the given coordinates
      */
    virtual unsigned long get_key(Vec& coords) = 0;

    /**
      * Get the hilbert keys of all blocks that overlap with the region defined
      * by the Vec and radius around it. This is done by systematically
      * subdividing the DelCont volume into blocks. Every block gets a key
      * (similar to the key of a TreeNode). If a block is entirely inside the
      * desired region, we stop subdividing and add the key of the block to the
      * list. If a block is entirely outside, we discard it.
      *
      * @param coords A Vec specifying the coordinates of a Particle or VorGen
      * @param radius A double specifying the size of the region around the
      * given coordinates
      * @param keys An array to fill with keys (length should be at least 9
      * unsigned longs in 2D and 27 in 3D)
      */
    virtual void get_ngb_keys(Vec& coords, double radius,
                              unsigned long* keys) = 0;

    /**
      * Make sure that the given Particle has coordinates inside the DelCont,
      * making use of periodicity (if applicable)
      *
      * @param p A Particle
      */
    virtual void keep_inside(Particle* p) = 0;

    /**
      * Convert a distance Vec between two particles to the smallest possible
      * distance Vec between these two particles, using periodic copies (if
      * applicable). This means that if one of the distance coordinates is
      * larger than half the size of the DelCont, this component is shifted by
      * the size of the DelCont and similar for a coordinate which is too small.
      *
      * @param v A Vec
      */
    virtual void closest_copy(Vec& v) = 0;

    /**
     * @brief Get a Cuboid specifying the exterior extents of the simulation box
     *
     * @return The dimensions of the smallest cuboid containing the entire
     * simulation box
     */
    virtual Cuboid get_cuboid() = 0;
};

#endif
