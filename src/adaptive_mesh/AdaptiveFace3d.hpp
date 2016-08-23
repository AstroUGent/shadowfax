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
 * @file AdaptiveFace3d.hpp
 *
 * @brief 3D face for mesh evolution flux calculation: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_ADAPTIVEFACE3D
#define HEAD_ADAPTIVEFACE3D

#include "Vec.hpp"  // for Vec

class GasParticle;
class RiemannSolver;
class StateVector;
class TimeLine;

/**
 * @brief Face of the 3D adaptive Voronoi tesselation
 *
 * Used for the flux calculation in the mesh evolution algorithm. Together with
 * AdaptiveVorFace3d, this class has the same functionality as the VorFace class
 * in the old algorithm.
 */
class AdaptiveFace3d {
  private:
    /*! @brief GasParticle at the left of the interface */
    GasParticle* _left;

    /*! @brief GasParticle at the right of the interface */
    GasParticle* _right;

    /*! @brief Position of the actual representative of the GasParticle at the
     *  right of the interface */
    Vec _rightpos;

    /*! @brief Velocity of the interface */
    double _v[3];

    /*! @brief Midpoint of the interface */
    double _midpoint[3];

    /*! @brief Area of the interface */
    double _area;

    /*! @brief 3D rotation angles: \f$\cos(\phi)\f$ */
    double _cosp;

    /*! @brief 3D rotation angles: \f$\sin(\phi)\f$ */
    double _sinp;

    /*! @brief 3D rotation angles: \f$\cos(\theta)\f$ */
    double _cost;

    /*! @brief 3D rotation angles: \f$\sin(\theta)\f$ */
    double _sint;

    void transform(StateVector& W);
    void invtransform(StateVector& W);
    void get_normal(double* angles);

  public:
    AdaptiveFace3d(GasParticle* left, GasParticle* right, Vec& rightpos,
                   double* midpoint, double area);

    void set_v();
    void calculate_flux(TimeLine& timeline, RiemannSolver& solver);
};

#endif
