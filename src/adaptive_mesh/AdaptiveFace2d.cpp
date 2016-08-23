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
 * @file AdaptiveFace2d.cpp
 *
 * @brief 2D face for mesh evolution flux calculation: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdaptiveFace2d.hpp"
#include "StateVector.hpp"            // for StateVector, operator*
#include "TimeLine.hpp"               // for TimeLine
#include "VorFace.hpp"                // for FluxCalculationException
#include "riemann/RiemannSolver.hpp"  // for RiemannSolver
#include "utilities/GasParticle.hpp"  // for GasParticle
#include <algorithm>                  // for max, min
#include <cmath>                      // for sqrt
using namespace std;

/**
 * @brief Constructor.
 *
 * Unlike for the old algorithm, computation of the geometrical quantities
 * associated to the face is done in a separate class. Here, we just pass them
 * on as parameters to the constructor.
 *
 * The only computation done inside the constructor is that of the angle between
 * the interface normal and the positive x-axis, through its sine and cosine.
 *
 * @param left GasParticle at the left of the interface
 * @param right GasParticle at the right of the interface
 * @param rightpos Position of the actual representative of the GasParticle at
 * the right of the interface
 * @param midpoint Midpoint of the interface
 * @param area Area of the interface
 */
AdaptiveFace2d::AdaptiveFace2d(GasParticle* left, GasParticle* right,
                               Vec& rightpos, double* midpoint, double area) {
    _left = left;
    _right = right;
    _rightpos = rightpos;
    _v[0] = 0.;
    _v[1] = 0.;
    _midpoint[0] = midpoint[0];
    _midpoint[1] = midpoint[1];
    _area = area;

    double d[3];
    d[0] = _left->x() - _rightpos.x();
    d[1] = _left->y() - _rightpos.y();
    d[2] = sqrt(d[0] * d[0] + d[1] * d[1]);
    // t is the angle between the face normal and the x-axis
    _cost = -d[0] / d[2];
    _sint = -d[1] / d[2];
}

/**
 * @brief Rotate the velocity components of the given StateVector to a frame
 * where the x-axis is the interface normal.
 *
 * @param W StateVector that has to be transformed
 */
void AdaptiveFace2d::transform(StateVector& W) {
    double W1 = W[1] * _cost + W[2] * _sint;
    W[2] = -W[1] * _sint + W[2] * _cost;
    W[1] = W1;
}

/**
 * @brief Rotate the velocity components of the given StateVector back from a
 * frame where the x-axis is the interface normal to a static reference frame
 *
 * @param W StateVector that has to be transformed
 */
void AdaptiveFace2d::invtransform(StateVector& W) {
    double W1 = W[1] * _cost - W[2] * _sint;
    W[2] = W[1] * _sint + W[2] * _cost;
    W[1] = W1;
}

/**
 * @brief Set the velocity of the interface based on the velocities of its left
 * and right state
 */
void AdaptiveFace2d::set_v() {
    if(_right) {
        double rRL[3];
        rRL[0] = _rightpos.x() - _left->x();
        rRL[1] = _rightpos.y() - _left->y();
        rRL[2] = rRL[0] * rRL[0] + rRL[1] * rRL[1];
        double fac = (_left->vx() - _right->vx()) *
                     (_midpoint[0] - 0.5 * (_left->x() + _rightpos.x()));
        fac += (_left->vy() - _right->vy()) *
               (_midpoint[1] - 0.5 * (_left->y() + _rightpos.y()));
        fac /= rRL[2];
        _v[0] = 0.5 * (_left->vx() + _right->vx()) + fac * rRL[0];
        _v[1] = 0.5 * (_left->vy() + _right->vy()) + fac * rRL[1];
    }
    // if _right is a ghost, then the face should not move
}

/**
 * @brief Calculate the flux passing through this face
 *
 * Per construction, the _left particle will always be a real particle and
 * correspond to a complete cell around the position of this particle.
 * For _right, there are four different cases:
 *  -# _right is also a real particle, corresponding to a cell (normal case)
 *  -# _right is NULL and corresponds to a reflective ghost particle of _left
 *  -# _right is NULL and corresponds to a periodic copy of an existing cell
 *  -# _right is NULL and corresponds to a periodic boundary ghost
 *
 * Both _left and _right can be active or inactive. However, when _right is
 * a periodic boundary ghost, left cannot be active. In this case, we can skip
 * the face entirely. If both _left and _right are inactive, we can also skip
 * the face (but these faces should not be present in the calculation anyway).
 *
 * Also note that case 2. and cases 3. and 4. are mutually exclusive, since only
 * one type of boundary conditions can be imposed.
 *
 * @param timeline TimeLine of the simulation, needed to convert integer time
 * to real time
 * @param solver Solver used to solve the Riemann problem at the interface
 */
void AdaptiveFace2d::calculate_flux(TimeLine& timeline, RiemannSolver& solver) {
#if ndim_ == 2
    if(!_area) {
        return;
    }
    Vec v(_v[0], _v[1]);

    double dt = timeline.get_realtime_interval(_left->get_timestep());
    if(_right) {
        dt = std::min(dt,
                      timeline.get_realtime_interval(_right->get_timestep()));
    }

    StateVector WL = _left->get_Wvec();
    StateVector WR;
    if(_right) {
        WR = _right->get_Wvec();
    }
    WL -= v;
    WR -= v;
#ifndef NOMUSCL
    for(unsigned int i = 2; i--;) {
        StateVector Wtemp;
        GasParticle* cell;
        Vec pos;
        if(i) {
            Wtemp = WL;
            cell = _left;
            pos = _left->get_position();
        } else {
            if(!_right) {
                continue;
            }
            Wtemp = WR;
            cell = _right;
            pos = _rightpos;
        }
        double d[ndim_];
        Vec centroid = cell->get_centroid();
        // no correction needed, since the centroid for passive cells is
        // translated together with the cell generator
        for(unsigned int n = ndim_; n--;) {
            double correction = cell->pos(n) - pos[n];
            d[n] = _midpoint[n] - centroid[n] + correction;
        }
        double dtp = dt;
        for(unsigned int m = ndim_ + 2; m--;) {
            Wtemp[m] += cell->get_gradient(m, 0) * d[0] +
                        cell->get_gradient(m, 1) * d[1];
        }
        double Ax[4][4] = {{Wtemp[1], Wtemp[0], 0., 0.},
                           {0., Wtemp[1], 0., 1. / Wtemp[0]},
                           {0., 0., Wtemp[1], 0.},
                           {0., solver.get_gamma() * Wtemp[3], 0., Wtemp[1]}};
        double Ay[4][4] = {{Wtemp[2], Wtemp[0], 0., 0.},
                           {0., Wtemp[2], 0., 0.},
                           {0., 0., Wtemp[2], 1. / Wtemp[0]},
                           {0., 0., solver.get_gamma() * Wtemp[3], Wtemp[2]}};
        StateVector deltaW;
        for(unsigned int m = ndim_ + 2; m--;) {
            for(unsigned int n = ndim_ + 2; n--;) {
                deltaW[m] -= 0.5 * dtp * Ax[m][n] * cell->get_gradient(n, 0);
                deltaW[m] -= 0.5 * dtp * Ay[m][n] * cell->get_gradient(n, 1);
            }
        }
        for(unsigned int m = ndim_ + 2; m--;) {
            Wtemp[m] += deltaW[m];
        }
        Wtemp[0] = std::max(0., Wtemp[0]);
        Wtemp[ndim_ + 1] = std::max(0., Wtemp[ndim_ + 1]);

        if(i) {
            WL = Wtemp;
        } else {
            WR = Wtemp;
        }
    }
#endif
    transform(WL);
    transform(WR);
    // no particle means we have a mirror: we mirror the left state
    if(!_right) {
        WR = WL;
        WR[1] = -WR[1];
    }
    StateVector Whalf;
    double maxmach = 0.;
    Vec n;
    n[0] = 1.;
    Whalf = solver.solve(WL, WR, n, maxmach);
    invtransform(Whalf);
    Whalf += v;

    _left->set_max_mach(std::max(_left->get_max_mach(), maxmach));
    if(_right) {
        _right->set_max_mach(std::max(_right->get_max_mach(), maxmach));
    }

    StateVector fluxvec[ndim_];
    for(unsigned int i = ndim_; i--;) {
        fluxvec[i] = solver.get_flux(v, i, Whalf);
    }
    double angles[ndim_];
    angles[0] = _cost;
    angles[1] = _sint;
    StateVector dQL;
    StateVector dQR;
    for(unsigned int i = ndim_; i--;) {
        dQL += dt * _area * angles[i] * fluxvec[i];
        dQR -= dt * _area * angles[i] * fluxvec[i];
    }

    _left->increase_dQ(dQL);
    if(_right) {
        _right->increase_dQ(dQR);
    }

    if(fluxvec[0][0] != fluxvec[0][0]) {
        throw FluxCalculationException(WL, WR, 0, solver);
    }
#endif
}
