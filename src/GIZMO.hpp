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
 * @file GIZMO.hpp
 *
 * @brief Implementation of Phil Hopkins' GIZMO: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef GIZMO_HPP
#define GIZMO_HPP

class ParticleVector;

/**
 * @brief Sideprogram with an implementation of Phil Hopkins' mesh free methods
 */
class GIZMO {
  private:
    double kernel(double r, double h);
    double kernel_d(double r, double h);

    void invert_matrix(double* M, double* IM);

    void test_tree(ParticleVector& particles);

  public:
    GIZMO(int argc, char** argv);
};

#endif  // GIZMO_HPP
