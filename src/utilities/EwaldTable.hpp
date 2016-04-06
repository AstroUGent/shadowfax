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
 * @file EwaldTable.hpp
 *
 * @brief Ewald tables used for periodic force corrections: header
 *
 * This code is largely based on similar code in Gadget2 (Springel 2005).
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef EWALDTABLE_HPP
#define EWALDTABLE_HPP

#include "Vec.hpp"

/**
 * \brief Ewald force corrections for N-body gravity with periodic boundary
 * conditions
 *
 * Following Springel (2005) and Hernquist, Bouchet and Suto (1991), we tabulate
 * the Ewald correction that has to be added to the gravitational force (or mesh
 * regularization displacement) to account for periodic boundaries.
 *
 * The precise details are still a bit vague, but we basically copied the
 * expressions from the Gadget2 source file forcetree.c. These expressions
 * correspond to eq. 2.14b in Hernquist et al. (1991), with an extra factor
 * \f$r\f$ in the second term between brackets (which is probably wrong in the
 * article; dimensionally we need a factor \f$r\f$).
 *
 * Since constructing large multidimensional arrays as class members causes
 * problems, we dynamically allocate and deallocate memory for the internal
 * table.
 */
class EwaldTable {
  private:
    void construct();
    bool read_table();

#if ndim_ == 3
    /*! \brief A multidimensional array of Vec force corrections */
    Vec*** _f;
#else
    /*! \brief A multidimensional array of Vec force corrections */
    Vec** _f;
#endif

    /*! \brief Factor that maps a double coordinate to a floating point "index"
     *  in the internal table */
    double _ewald_intp_fac;

    /*! \brief Ewald alpha factor that determines the cutoff between long and
     *  short range forces */
    double _alpha;

    /*! \brief Size of the internal table (in one dimension). The actual size of
     *  the table is a power of this */
    unsigned int _size;

  public:
    EwaldTable(double alpha = 2., unsigned int size = 64);
    ~EwaldTable();

    Vec get_correction(Vec& position);
    void set_boxsize(double L);
};

#endif  // EWALDTABLE_HPP
