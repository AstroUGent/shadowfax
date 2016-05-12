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
 * @file Hilbert.cpp
 *
 * @brief Space filling Hilbert curve: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Hilbert.hpp"
#include "MPIMethods.hpp"
#include "RestartFile.hpp"
#include <mpi.h>
using namespace std;

#if ndim_ == 3
static const unsigned int t[12][8][2] = {
        {{5, 0}, {1, 7}, {4, 1}, {2, 6}, {3, 3}, {3, 4}, {4, 2}, {2, 5}},
        {{6, 4}, {2, 7}, {6, 3}, {8, 0}, {0, 5}, {0, 6}, {7, 2}, {7, 1}},
        {{1, 6}, {0, 7}, {1, 5}, {9, 4}, {10, 1}, {11, 0}, {10, 2}, {9, 3}},
        {{9, 2}, {8, 5}, {0, 3}, {0, 4}, {9, 1}, {8, 6}, {6, 0}, {10, 7}},
        {{0, 0}, {5, 1}, {8, 3}, {5, 2}, {11, 7}, {6, 6}, {8, 4}, {6, 5}},
        {{4, 0}, {10, 3}, {9, 7}, {10, 4}, {0, 1}, {0, 2}, {7, 6}, {7, 5}},
        {{11, 6}, {11, 5}, {3, 1}, {3, 2}, {4, 7}, {1, 4}, {9, 0}, {1, 3}},
        {{9, 6}, {8, 1}, {5, 7}, {1, 0}, {9, 5}, {8, 2}, {11, 4}, {11, 3}},
        {{1, 2}, {4, 3}, {1, 1}, {7, 0}, {10, 5}, {4, 4}, {10, 6}, {3, 7}},
        {{2, 4}, {5, 5}, {7, 7}, {5, 6}, {2, 3}, {6, 2}, {3, 0}, {6, 1}},
        {{11, 2}, {11, 1}, {3, 5}, {3, 6}, {5, 3}, {2, 0}, {5, 4}, {8, 7}},
        {{7, 4}, {7, 3}, {4, 5}, {2, 2}, {6, 7}, {10, 0}, {4, 6}, {2, 1}}};

static const unsigned int ti[12][8][4] = {{{5, 0, 0, 0},
                                           {4, 0, 0, 1},
                                           {4, 0, 1, 1},
                                           {3, 0, 1, 0},
                                           {3, 1, 1, 0},
                                           {2, 1, 1, 1},
                                           {2, 1, 0, 1},
                                           {1, 1, 0, 0}},
                                          {{8, 1, 0, 1},
                                           {7, 1, 1, 1},
                                           {7, 0, 1, 1},
                                           {6, 0, 0, 1},
                                           {6, 0, 0, 0},
                                           {0, 0, 1, 0},
                                           {0, 1, 1, 0},
                                           {2, 1, 0, 0}},
                                          {{11, 1, 1, 0},
                                           {10, 0, 1, 0},
                                           {10, 0, 1, 1},
                                           {9, 1, 1, 1},
                                           {9, 1, 0, 1},
                                           {1, 0, 0, 1},
                                           {1, 0, 0, 0},
                                           {0, 1, 0, 0}},
                                          {{6, 0, 1, 1},
                                           {9, 0, 1, 0},
                                           {9, 0, 0, 0},
                                           {0, 0, 0, 1},
                                           {0, 1, 0, 1},
                                           {8, 1, 0, 0},
                                           {8, 1, 1, 0},
                                           {10, 1, 1, 1}},
                                          {{0, 0, 0, 0},
                                           {5, 1, 0, 0},
                                           {5, 1, 0, 1},
                                           {8, 0, 0, 1},
                                           {8, 0, 1, 1},
                                           {6, 1, 1, 1},
                                           {6, 1, 1, 0},
                                           {11, 0, 1, 0}},
                                          {{4, 0, 0, 0},
                                           {0, 0, 1, 0},
                                           {0, 1, 1, 0},
                                           {10, 1, 0, 0},
                                           {10, 1, 0, 1},
                                           {7, 1, 1, 1},
                                           {7, 0, 1, 1},
                                           {9, 0, 0, 1}},
                                          {{9, 0, 1, 1},
                                           {3, 0, 0, 1},
                                           {3, 1, 0, 1},
                                           {1, 1, 1, 1},
                                           {1, 1, 1, 0},
                                           {11, 1, 0, 0},
                                           {11, 0, 0, 0},
                                           {4, 0, 1, 0}},
                                          {{1, 1, 0, 1},
                                           {8, 1, 0, 0},
                                           {8, 1, 1, 0},
                                           {11, 1, 1, 1},
                                           {11, 0, 1, 1},
                                           {9, 0, 1, 0},
                                           {9, 0, 0, 0},
                                           {5, 0, 0, 1}},
                                          {{7, 1, 0, 1},
                                           {1, 0, 0, 1},
                                           {1, 0, 0, 0},
                                           {4, 1, 0, 0},
                                           {4, 1, 1, 0},
                                           {10, 0, 1, 0},
                                           {10, 0, 1, 1},
                                           {3, 1, 1, 1}},
                                          {{3, 0, 1, 1},
                                           {6, 1, 1, 1},
                                           {6, 1, 1, 0},
                                           {2, 0, 1, 0},
                                           {2, 0, 0, 0},
                                           {5, 1, 0, 0},
                                           {5, 1, 0, 1},
                                           {7, 0, 0, 1}},
                                          {{2, 1, 1, 0},
                                           {11, 1, 0, 0},
                                           {11, 0, 0, 0},
                                           {5, 0, 1, 0},
                                           {5, 0, 1, 1},
                                           {3, 0, 0, 1},
                                           {3, 1, 0, 1},
                                           {8, 1, 1, 1}},
                                          {{10, 1, 1, 0},
                                           {2, 1, 1, 1},
                                           {2, 1, 0, 1},
                                           {7, 1, 0, 0},
                                           {7, 0, 0, 0},
                                           {4, 0, 0, 1},
                                           {4, 0, 1, 1},
                                           {6, 0, 1, 0}}};

// xyz
static const unsigned int filltable[2][2][2][27][2] = {
        {{{{0, 7},  {1, 3},  {1, 7},  {4, 5},  {4, 1},  {5, 5},  {5, 7},
           {4, 3},  {4, 7},  {13, 6}, {13, 2}, {12, 6}, {12, 4}, {13, 0},
           {13, 4}, {16, 6}, {16, 2}, {17, 6}, {17, 7}, {16, 3}, {16, 7},
           {13, 5}, {13, 1}, {12, 5}, {12, 7}, {13, 3}, {13, 7}},  // 000
          {{1, 3},  {1, 7},  {2, 3},  {3, 1},  {4, 5},  {4, 1},  {4, 3},
           {4, 7},  {3, 3},  {14, 2}, {13, 6}, {13, 2}, {13, 0}, {13, 4},
           {14, 0}, {15, 2}, {16, 6}, {16, 2}, {16, 3}, {16, 7}, {15, 3},
           {14, 1}, {13, 5}, {13, 1}, {13, 3}, {13, 7}, {14, 3}}},  // 100
         {{{5, 5},  {4, 1},  {4, 5},  {4, 7},  {4, 3},  {5, 7},  {6, 5},
           {7, 1},  {7, 5},  {10, 4}, {10, 0}, {11, 4}, {12, 6}, {13, 2},
           {13, 6}, {13, 4}, {13, 0}, {12, 4}, {12, 5}, {13, 1}, {13, 5},
           {13, 7}, {13, 3}, {12, 7}, {11, 5}, {10, 1}, {10, 5}},  // 010
          {{4, 1},  {4, 5},  {3, 1},  {3, 3},  {4, 7},  {4, 3},  {7, 1},
           {7, 5},  {8, 1},  {9, 0},  {10, 4}, {10, 0}, {13, 2}, {13, 6},
           {14, 2}, {14, 0}, {13, 4}, {13, 0}, {13, 1}, {13, 5}, {14, 1},
           {14, 3}, {13, 7}, {13, 3}, {10, 1}, {10, 5}, {9, 1}}}},  // 110
        {{{{17, 6}, {16, 2}, {16, 6}, {13, 4}, {13, 0}, {12, 4}, {12, 6},
           {13, 2}, {13, 6}, {13, 7}, {13, 3}, {12, 7}, {12, 5}, {13, 1},
           {13, 5}, {16, 7}, {16, 3}, {17, 7}, {18, 6}, {19, 2}, {19, 6},
           {22, 4}, {22, 0}, {23, 4}, {23, 6}, {22, 2}, {22, 6}},  // 001
          {{16, 2}, {16, 6}, {15, 2}, {14, 0}, {13, 4}, {13, 0}, {13, 2},
           {13, 6}, {14, 2}, {14, 3}, {13, 7}, {13, 3}, {13, 1}, {13, 5},
           {14, 1}, {15, 3}, {16, 7}, {16, 3}, {19, 2}, {19, 6}, {20, 2},
           {21, 0}, {22, 4}, {22, 0}, {22, 2}, {22, 6}, {21, 2}}},  // 101
         {{{12, 4}, {13, 0}, {13, 4}, {13, 6}, {13, 2}, {12, 6}, {11, 4},
           {10, 0}, {10, 4}, {10, 5}, {10, 1}, {11, 5}, {12, 7}, {13, 3},
           {13, 7}, {13, 5}, {13, 1}, {12, 5}, {23, 4}, {22, 0}, {22, 4},
           {22, 6}, {22, 2}, {23, 6}, {24, 4}, {25, 0}, {25, 4}},  // 011
          {{13, 0}, {13, 4}, {14, 0}, {14, 2}, {13, 6}, {13, 2}, {10, 0},
           {10, 4}, {9, 0},  {9, 1},  {10, 5}, {10, 1}, {13, 3}, {13, 7},
           {14, 3}, {14, 1}, {13, 5}, {13, 1}, {22, 0}, {22, 4}, {21, 0},
           {21, 2}, {22, 6}, {22, 2}, {25, 0}, {25, 4}, {26, 0}}}}  // 111
};

#define numfill_ 27
#else
/**
 * @brief Basic Hilbert table
 *
 * This table relates the bits on a level to parts of the key and the
 * orientation of the next level.
 */
static const unsigned int t[8][4][2] = {
        {{7, 0}, {0, 1}, {6, 3}, {0, 2}}, {{1, 2}, {7, 3}, {1, 1}, {6, 0}},
        {{2, 1}, {2, 2}, {4, 0}, {5, 3}}, {{4, 3}, {5, 0}, {3, 2}, {3, 1}},
        {{3, 3}, {4, 2}, {2, 0}, {4, 1}}, {{5, 1}, {3, 0}, {5, 2}, {2, 3}},
        {{6, 2}, {6, 1}, {0, 3}, {1, 0}}, {{0, 0}, {1, 3}, {7, 1}, {7, 2}}};

/**
 * @brief Inverse Hilbert table
 *
 * This table relates a part of a key on a level to the bits on that level an
 * the orientation of the next level.
 */
static const unsigned int ti[8][4][3] = {
        {{7, 0, 0}, {0, 0, 1}, {0, 1, 1}, {6, 1, 0}},
        {{6, 1, 1}, {1, 1, 0}, {1, 0, 0}, {7, 0, 1}},
        {{4, 1, 0}, {2, 0, 0}, {2, 0, 1}, {5, 1, 1}},
        {{5, 0, 1}, {3, 1, 1}, {3, 1, 0}, {4, 0, 0}},
        {{2, 1, 0}, {4, 1, 1}, {4, 0, 1}, {3, 0, 0}},
        {{3, 0, 1}, {5, 0, 0}, {5, 1, 0}, {2, 1, 1}},
        {{1, 1, 1}, {6, 0, 1}, {6, 0, 0}, {0, 1, 0}},
        {{0, 0, 0}, {7, 1, 0}, {7, 1, 1}, {1, 0, 1}}};

/**
 * @brief Hilbert neighbour table
 *
 * This table relates the bits on a level to the keys on that level of the
 * neighbouring blocks. Due to the different orientations of blocks on a certain
 * level, it is a priori impossible to tell what the key of a block on the next
 * level will be, if you only know the orientation of one of its neighbours on
 * that level. This is exactly what this table tabulates: it tells you, given
 * the orientation of a block, what the orientation (and hence key part) of its
 * neighbours will be, even if these neighbours are in another block on a
 * coarser level.
 */
static const unsigned int filltable[2][2][9][2] =
        // xy
        // 00
        {{{{0, 3},
           {1, 2},
           {1, 3},
           {8, 1},
           {8, 3},
           {8, 2},
           {7, 3},
           {7, 1},
           {8, 0}},
          // 01
          {{1, 2},
           {1, 3},
           {2, 2},
           {3, 0},
           {3, 2},
           {8, 3},
           {8, 2},
           {8, 0},
           {8, 1}}},
         // 10
         {{{7, 1},
           {8, 0},
           {8, 1},
           {8, 3},
           {5, 1},
           {5, 0},
           {6, 1},
           {7, 3},
           {8, 2}},
          // 11
          {{8, 0},
           {8, 1},
           {3, 0},
           {3, 2},
           {4, 0},
           {5, 1},
           {5, 0},
           {8, 2},
           {8, 3}}}};

#define numfill_ 9
#endif

/*! @brief Placeholder bit added to all keys to mark them as keys on the lowest
 *  level */
static const unsigned long placeholder = ((unsigned long)1) << 60;

// CODE KEPT FOR CONVENIENCE, BECAUSE IT HOLDS VITAL INFORMATION ABOUT THE
// CALCULATION OF HILBERT KEYS
///**
// * Convert a set of integer coordinates to a key with given length.
// * @param bits A vector of integer coordinates (in 3 dimensions)
// * @param nbits The length of the key. Can not be more than 64 (on a 64-bit
// * machine)
// * @return An integer key
// */
//#if ndim_==3
// unsigned long HB::get_key(unsigned long* bits, unsigned int nbits){
//    unsigned long key = 0;
//    unsigned long mask = 1;
//    mask <<= nbits-1;
//    bool x[ndim_];
//    unsigned int ci;
//    unsigned int si = 4;
//    for(unsigned int i = nbits; i--;){
//        key <<= ndim_;
//        for(unsigned int j = ndim_; j--;){
//            x[j] = (bits[j] & mask);
//        }
//        ci = (x[0]<<2) | (x[1]<<1) | x[2];
//        key += t[si][ci][1];
//        si = t[si][ci][0];
//        mask >>= 1;
//    }
//    return key;
//}
//#else
// unsigned long HB::get_key(unsigned long* bits, unsigned int nbits){
//  unsigned long key = 0;
//  unsigned long mask = 1;
//  mask <<= nbits-1;
//  unsigned int si = 0;
//  for(unsigned int i = nbits; i--;){
//    key <<= 2;
//    bool x = (bits[0] & mask);
//    bool y = (bits[1] & mask);
//    unsigned int ci = (x<<1) | y;
//    key += t2[si][ci][1];
//    si = t2[si][ci][0];
//    mask >>= 1;
//  }
//  return key;
//}
//#endif

/**
  * \brief Convert integer coordinates to a key of given length
  *
  * This function takes an array of 2 or 3 integer type coordinates
  * and converts them to a hilbert key with the desired length.
  *
  * @param bits Integer coordinates (64-bit)
  * @param nbits Desired length of the key in number of bits
  * @return An integer type hilbert key
  */
unsigned long HB::get_key(unsigned long* bits, unsigned int nbits) {
    unsigned long key = 0;
    unsigned long mask = 1;
    mask <<= nbits - 1;
    bool x[ndim_];
    unsigned int ci;
#if ndim_ == 3
    unsigned int si = 4;
#else
    unsigned int si = 7;
#endif
    for(unsigned int i = nbits; i--;) {
        key <<= ndim_;
        for(unsigned int j = ndim_; j--;) {
            x[j] = (bits[j] & mask);
        }
#if ndim_ == 3
        ci = (x[0] << 2) | (x[1] << 1) | x[2];
#else
        ci = (x[0] << 1) | x[1];
#endif
        key += t[si][ci][1];
        si = t[si][ci][0];
        mask >>= 1;
    }
    // add placeholder
    key += placeholder;
    return key;
}

/**
  * \brief Convert a key of given length to integer coordinates
  *
  * This function takes a key of given length and calculates integer
  * coordinates from this key, so that evaluation of this function on
  * the result of HB::get_key yields the original integer coordinates.
  *
  * @param key An integer hilbert key
  * @param nbits The length of the given key in number of bits
  * @param coords A 2- or 3-element array of bits, initialized to 0, in which
  * the resulting coordinates will be stored
  */
void HB::get_coords(unsigned long key, unsigned int nbits,
                    unsigned long* coords) {
#if ndim_ == 3
    unsigned int mask = 7;
    unsigned int si = 5;
#else
    unsigned int mask = 3;
    unsigned int si = 7;
#endif
    unsigned int ci;
    for(unsigned int i = nbits + ndim_; i -= ndim_;) {
        ci = (key >> (i - ndim_)) & mask;
        for(unsigned int j = ndim_; j--;) {
            coords[j] <<= 1;
            coords[j] += ti[si][ci][j + 1];
        }
        si = ti[si][ci][0];
    }
}

/**
  * \brief Calculate the neighbouring keys for the given integer coordinates on
  * a given level and with a given length
  *
  * This function is equivalent to HB::get_key, but then not only for the given
  * integer coordinates, but also for its neighbouring blocks on the given
  * level. 9 or 27 keys are calculated at the same time, depending upon the
  * number of dimensions of the code.
  *
  * @param bits Integer coordinates (64-bit)
  * @param nbits The desired length of the resulting keys in number of bits
  * @param level The desired level at which the considered blocks reside (every
  * level divides the 2D or 3D space into 2^(ndim_*level) blocks)
  * @param keys A 9- or 27-element integer array to store the keys in
  * (initialization not required)
  */
void HB::get_ngb_keys(unsigned long* bits, unsigned int nbits,
                      unsigned int level, unsigned long* keys) {
    unsigned int si[numfill_] = {0};
    unsigned int newsi[numfill_] = {0};
    unsigned long ngbs[numfill_] = {0};
    bool x[ndim_];
#if ndim_ == 3
    si[13] = 4;
    keys[13] = 1;
#else
    si[8] = 7;
    keys[8] = 1;
#endif
    unsigned long mask = 1;
    mask <<= nbits - 1;
    for(unsigned int i = level; i--;) {
        for(unsigned int j = ndim_; j--;) {
            x[j] = (bits[j] & mask);
        }
        for(unsigned int j = numfill_; j--;) {
#if ndim_ == 3
            const unsigned int* index = filltable[x[2]][x[1]][x[0]][j];
#else
            const unsigned int* index = filltable[x[0]][x[1]][j];
#endif
            if(keys[index[0]]) {
                ngbs[j] = (keys[index[0]] << ndim_) +
                          t[si[index[0]]][index[1]][1];
                newsi[j] = t[si[index[0]]][index[1]][0];
            }
        }
        for(unsigned int j = numfill_; j--;) {
            keys[j] = ngbs[j];
            si[j] = newsi[j];
        }
        mask >>= 1;
    }
}

/**
  * \brief Convert a hilbert key to a TreeRoute
  *
  * This function takes a hilbert key and converts it to a roadmap for Tree
  * traversal that leads to the block pointed to by the key. This roadmap is
  * stored in a TreeRoute object.
  *
  * @param key An integer 64-bit hilbert key with 60 significant bits
  * @return A TreeRoute that points to the block with the given key
  */
TreeRoute HB::get_route_to_node(unsigned long key) {
    TreeRoute route;
#if ndim_ == 3
    unsigned int mask = 7;
#else
    unsigned int mask = 3;
#endif
    unsigned int nbits = 60;
    while(!(key >> nbits)) {
        nbits -= ndim_;
    }
    for(unsigned int i = nbits + ndim_; i -= ndim_;) {
        route.add_step(((key >> (i - ndim_)) & mask));
    }
    return route;
}

/**
  * \brief Sort hilbert objects along a space filling curve
  *
  * Sort function that can be used in combination with std::sort to sort
  * Hilbert_Object objects. If you want to be able to sort objects in this way,
  * it suffices to let them inherit from the Hilbert_Object class. All required
  * functionalities are provided by this class.
  *
  * @param a,b Two Hilbert_Object objects with a valid hilbert key
  * @return true if the key of a is smaller than that of b, false in all other
  * cases
  */
bool HB::sortfunc(Hilbert_Object* a, Hilbert_Object* b) {
    return a->get_key() < b->get_key();
}

/**
  * \brief Empty Hilbert_Object constructor
  *
  * Initialize a Hilbert_Object with key 0.
  */
Hilbert_Object::Hilbert_Object() {
    _key = 0;
}

/**
  * \brief Set the hilbert key of a Hilbert_Object
  *
  * @param key A hilbert key
  */
void Hilbert_Object::set_key(unsigned long key) {
    _key = key;
}

/**
  * \brief Access the hilbert key of the Hilbert_Object
  *
  * @return The hilbert key of the Hilbert_Object
  */
unsigned long Hilbert_Object::get_key() {
    return _key;
}

/**
 * @brief MPI constructor. Initialize the Hilbert object using the given
 * communication buffer
 *
 * @param buffer Buffer to read from
 * @param bufsize Buffer size
 * @param position Position in the buffer (is updated)
 */
Hilbert_Object::Hilbert_Object(void* buffer, int bufsize, int* position) {
    MyMPI_Unpack(buffer, bufsize, position, &_key, 1, MPI_UNSIGNED_LONG);
}

/**
 * @brief Write object data to the given buffer for MPI communication
 *
 * @param buffer Buffer to write to
 * @param bufsize Buffer size
 * @param position Position in the buffer (is updated)
 */
void Hilbert_Object::pack_data(void* buffer, int bufsize, int* position) {
    MyMPI_Pack(&_key, 1, MPI_UNSIGNED_LONG, buffer, bufsize, position);
}

/**
 * @brief Dump the object to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void Hilbert_Object::dump(RestartFile& rfile) {
    rfile.write(_key);
}

/**
 * @brief Restart constructor. Initialize the object using the given RestartFile
 *
 * @param rfile RestartFile to read from
 */
Hilbert_Object::Hilbert_Object(RestartFile& rfile) {
    rfile.read(_key);
}
