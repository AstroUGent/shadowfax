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
 * @file CubicBox.cpp
 *
 * @brief Cubic box implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 * @author Sven De Rijcke (sven.derijcke@ugent.be)
 */
#include "DelCont.hpp"
#include "VorGen.hpp"
#include "utilities/Hilbert.hpp"
#include "utilities/Particle.hpp"
#include <cmath>
using namespace std;

/**
 * Initialize a CubicBox with given origin and side
 *
 * @param origin A Vec specifying the coordinates of the center of the box
 * @param side The side of the box
 */
CubicBox::CubicBox( Vec origin, double side){
    // _origin is assumed to be the CENTER of the cube...
    _origin = origin;
    _side = side;
    _bitwidth = 1;
    _bitwidth <<= (60/ndim_);
}

/**
 * @brief Constructor
 *
 * Empty constructor initializing a CubicBox with center (0.5,0.5,0.5) (or
 * (0.5,0.5) in 2D) and side 10
 */
CubicBox::CubicBox(){
#if ndim_==3
    _origin = Vec(0.5, 0.5, 0.5);
#else
    _origin = Vec(0.5, 0.5);
#endif
    _side = 10.;
    _bitwidth = 1;
    _bitwidth <<= (60/ndim_);
}

/**
 * @brief Determine whether the given coordinates define a Point inside the
 * CubicBox
 *
 * Points on the boundary are defined to be inside.
 *
 * @param pntpos A Vec specifying the position of a VorGen or Particle
 * @return true if the coordinates are inside or on the boundary of the
 * CubicBox, false otherwise
 */
bool CubicBox::inside( Vec pntpos ){
    double tst = _side/2;
    for( int i=ndim_; i--; ){
        if( abs( _origin[i]-pntpos[i] ) > tst ){
            return false;
        }
    }
    return true;
}

/**
 * @brief Determine whether the volume defined by Point and radius around Point
 * is completely inside the CubicBox
 *
 * Points on the boundary are defined to be outside.
 *
 * @param point A VorGen that specifies the center of the search region
 * @param radius The radius determining the size of the search region
 * @return true if the specified region is completely inside the CubicBox
 * (borders exclusive), false otherwise
 */
bool CubicBox::inside( VorGen* point, double radius ){
    Vec pntpos = point->get_position( );
    double tst = _side/2;
    for( int i=ndim_; i--; ){
        if( abs( _origin[i]-pntpos[i] ) > tst-radius ){
            return false;
        }
    }
    return true;
}

/**
 * @brief Return the coordinates of the Tetrahedron encompassing the entire cube
 *
 * Because all mirror points have to be inside this Tetrahedron as well, the
 * margins are chosen quite large.
 *
 * @return A vector containing the 4 (or 3) triplets (/doublets) defining the
 * vertices of a Simplex
 */
vector<double> CubicBox::get_bounding_tetrahedron( ){
    vector<double> coords( (ndim_+1)*ndim_ );
#if ndim_==3
    // r = 4 * radius of sphere that contains cube
    double r = sqrt( 12. )*_side;
    coords[0] = _origin[0];
    coords[1] = _origin[1] + 3*r;
    coords[2] = _origin[2];

    coords[3] = _origin[0] + sqrt(8.0)*r;
    coords[4] = _origin[1] - r;
    coords[5] = _origin[2];

    coords[6] = _origin[0] - sqrt(2.0)*r;
    coords[7] = _origin[1] - r;
    coords[8] = _origin[2] + sqrt(6.0)*r;

    coords[9] = _origin[0] - sqrt(2.0)*r;
    coords[10] = _origin[1] - r;
    coords[11] = _origin[2] - sqrt(6.0)*r;
#else
    // r = 4 * radius of circle that contains square
    double r = sqrt( 8. )*_side;
    double l = sqrt( 3. )*r;
    coords[0] = _origin[0] - l;
    coords[1] = _origin[1] - r;

    coords[2] = _origin[0] + l;
    coords[3] = _origin[1] - r;

    coords[4] = _origin[0];
    coords[5] = _origin[1]+ 2*r;
#endif
    return coords;
}

/**
  * @brief Retrieve the origin and side of the box (not the center!)
  *
  * @param box A 3- or 4- dimensional array to fill
  */
void CubicBox::get_bounding_box(double* box){
    for(unsigned i = ndim_; i--;){
        box[i] = _origin[i] - 0.5*_side;
        box[ndim_+i] = _side;
    }
}

/**
 * @brief Return the width of the box
 *
 * Needed for tree construction.
 *
 * @return The side of the CubicBox
 */
double CubicBox::get_box_width( ){
    return _side;
}

/**
  * @brief Calculate a hilbert key for the given coordinates of a Particle or
  * VorGen that is inside the CubicBox, taking into account the size of the box
  *
  * @param coords A Vec specifying the coordinates of a Particle or VorGen
  * inside the box
  * @return The hilbert key for the given coordinates
  */
unsigned long CubicBox::get_key(Vec &coords){
    unsigned long bits[ndim_] = {0};
    for(unsigned int i = ndim_; i--;){
        bits[i] = ((coords[i]-_origin[i]+0.5*_side)/_side)*_bitwidth;
    }
    return HB::get_key(bits, 60);
}

/**
 * @brief Get the Hilbert keys of the neighbouring nodes in a tree for a node
 * with the given size and position
 *
 * @param coords Position inside the node for which we want the neighbours
 * @param radius Size of the node (measure for the level of the tree we're at)
 * @param keys Array to store the Hilbert keys of the neighbouring nodes in
 */
void CubicBox::get_ngb_keys(Vec &coords, double radius, unsigned long *keys){
    unsigned long bits[ndim_] = {0};
    unsigned int level = int(log(_side/radius))-1;
    for(unsigned int i = ndim_; i--;){
        bits[i] = ((coords[i]-_origin[i]+0.5*_side)/_side)*_bitwidth;
    }
#if ndim_==3
    HB::get_ngb_keys(bits, 20, level, keys);
#else
    HB::get_ngb_keys(bits, 30, level, keys);
#endif
}

/**
 * @brief Keep the given Particle inside the box by periodically trimming its
 * coordinates if necessary
 *
 * @param p Particle that can either be inside or outside the CubicBox
 */
void CubicBox::keep_inside(Particle *p){
    for(unsigned int i = ndim_; i--;){
        if(p->pos(i) < _origin[i]-0.5*_side){
            p->get_position()[i] += _side;
        }
        // if both boundaries are inclusive, this gives rise to problems
        // (because the boundaries are their own periodic copies)
        // hence the >=
        if(p->pos(i) >= _origin[i]+0.5*_side){
            p->get_position()[i] -= _side;
        }
    }
}

/**
 * @brief Get the closest copy of the given Vec to the point (0,0,0)
 *
 * @param v Vec specifying coordinates in Cartesian space
 */
void CubicBox::closest_copy(Vec &v){
    for(unsigned int i = ndim_; i--;){
        if(v[i] > 0.5*_side){
            v[i] -= _side;
        }
        if(v[i] < -0.5*_side){
            v[i] += _side;
        }
    }
}

/**
 * @brief Get a cuboid specifying the dimensions of the box
 *
 * @return Cuboid specifying the dimensions of the box
 */
Cuboid CubicBox::get_cuboid(){
#if ndim_==3
    Vec sides(_side, _side, _side);
#else
    Vec sides(_side, _side);
#endif
    Cuboid cuboid(_origin, sides);
    return cuboid;
}
