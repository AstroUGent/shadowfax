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
 * \file Hilbert.hpp
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *
 * @brief Space filling Hilbert curve: header
 *
 * contains namespace HB, class Hilbert
 *
 * The namespace HB defines a function that can be used to calculate %Hilbert
 * keys for integer coordinates, as an alternative for the Morton keys
 *
 * The Hilbert class is used to sample Particles along a %Hilbert space-filling
 * curve
 */
#ifndef HEAD_HILBERT
#define HEAD_HILBERT

#include <iostream>
#include "TreeRoute.hpp"

class RestartFile;

/**
  * \brief Abstract class for objects that can be ordered along a hilbert curve
  *
  * This object attaches a key to its children that can be used to sort the
  * children along a space-filling hilbert curve.
  *
  * Every class that inherits from this class can be sorted in this way.
  */
class Hilbert_Object{
private:
    /*! \brief The hilbert key */
    unsigned long _key;

public:
    Hilbert_Object();

    /**
      * Construct a Hilbert_Object from the given key
      *
      * @param key A hilbert key
      */
    Hilbert_Object(unsigned long key) : _key(key) {}
    Hilbert_Object(void* buffer, int bufsize, int* position);
    ~Hilbert_Object(){}

    void set_key(unsigned long key);
    unsigned long get_key();

    void pack_data(void* buffer, int bufsize, int* position);

    void dump(RestartFile &rfile);
    Hilbert_Object(RestartFile &rfile);
};

/**
 * \brief Methods to calculate Hilbert-keys
 *
 * This namespace holds a method to calculate a Hilbert-key for a given set of
 * integer coordinates, using a given key length.
 *
 * It also provides the basic sort function to sort a list of Hilbert_Object
 * objects
 */
namespace HB{
    unsigned long get_key(unsigned long* bits, unsigned int nbits);
    void get_coords(unsigned long key, unsigned int nbits,
                    unsigned long* coords);
    bool sortfunc(Hilbert_Object* a, Hilbert_Object* b);
    TreeRoute get_route_to_node(unsigned long key);
    void get_ngb_keys(unsigned long* bits, unsigned int nbits,
                      unsigned int level, unsigned long* keys);
}

#endif
