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
 * @file TenetCell.hpp
 *
 * @brief Cell used for Tenet implementation: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TENETCELL_HPP
#define TENETCELL_HPP

#include "StateVector.hpp"
#include "Vec.hpp"
#include <vector>

class TenetBasisFunctions;
class TenetCellWeights;

/**
 * @brief Cell of a Cartesian grid, contains the cell center and weights
 */
class TenetCell {
  private:
    /*! @brief Center of the cell */
    Vec _center;

    /*! @brief Weights of the hydrodynamical variables inside the cell */
    std::vector<StateVector> _weights;

    /*! @brief Time derivatives of the weights */
    std::vector<StateVector> _dt_weights;

  public:
    TenetCell(Vec center);

    Vec get_center();

    StateVector get_variables(Vec ksi, TenetCellWeights& weights,
                              TenetBasisFunctions& basis);
};

#endif  // TENETCELL_HPP
