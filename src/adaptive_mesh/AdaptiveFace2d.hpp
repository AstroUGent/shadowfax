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
 * @file AdaptiveFace2d.hpp
 *
 * @brief 2D face for mesh evolution flux calculation: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_ADAPTIVEFACE2D
#define HEAD_ADAPTIVEFACE2D

#include "Vec.hpp"  // for Vec

class GasParticle;
class RiemannSolver;
class StateVector;
class TimeLine;

/**
 * @brief Face of the 2D adaptive Voronoi tesselation
 *
 * Together with AdaptiveVorFace, this class has the same functionality as
 * VorFace for the old algorithm.
 */
class AdaptiveFace2d {
  private:
    /*! @brief GasParticle at the left of the interface */
    GasParticle* _left;

    /*! @brief GasParticle at the right of the interface */
    GasParticle* _right;

    /*! @brief Position of the actual representative of the GasParticle at the
     *  right of the interface */
    Vec _rightpos;

    /*! @brief Velocity of the interface */
    double _v[2];

    /*! @brief Midpoint of the interface */
    double _midpoint[2];

    /*! @brief Area of the interface */
    double _area;

    /*! @brief Sine of the angle between the interface normal and the positive
     *  x-axis */
    double _sint;

    /*! @brief Cosine of the angle between the interface normal and the positive
     *  x-axis */
    double _cost;

    void transform(StateVector& W);
    void invtransform(StateVector& W);

  public:
    AdaptiveFace2d(GasParticle* left, GasParticle* right, Vec& rightpos,
                   double* midpoint, double area);

    void set_v();
    void calculate_flux(TimeLine& timeline, RiemannSolver& solver);
};

#endif  // HEAD_ADAPTIVEFACE2D
