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
 * @file TenetCell.cpp
 *
 * @brief Cell used for Tenet implementation: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "TenetCell.hpp"
#include "TenetBasisFunctions.hpp"
#include "TenetCellWeights.hpp"
using namespace std;

/**
 * @brief Constructor
 *
 * @param center Cell center
 */
TenetCell::TenetCell(Vec center) {
    _center = center;
}

/**
 * @brief Access the center of the cell
 *
 * @return Center of the cell
 */
Vec TenetCell::get_center() {
    return _center;
}

/**
 * @brief Get the hydrodynamical variables at the given (rescaled) position in
 * the cell
 *
 * @param ksi Rescaled position inside the cell
 * @param weights Weights of the hydrodynamical variables inside the cell
 * @param basis Reference to the TenetBasisFunctions
 * @return Reconstructed hydrodynamical variables
 */
StateVector TenetCell::get_variables(Vec ksi, TenetCellWeights& weights,
                                     TenetBasisFunctions& basis) {
    StateVector U;

    for(auto basisfunction = basis.begin(); basisfunction != basis.end();
        ++basisfunction) {
        U += weights[basisfunction.index()] * basisfunction(ksi);
    }

    return U;
}
