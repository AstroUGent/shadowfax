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
 * @file AdaptiveFace3d.cpp
 *
 * @brief 3D face for mesh evolution flux calculation: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdaptiveFace3d.hpp"
#include "StateVector.hpp"            // for StateVector, operator*
#include "TimeLine.hpp"               // for TimeLine
#include "VorFace.hpp"                // for FluxCalculationException
#include "riemann/RiemannSolver.hpp"  // for RiemannSolver
#include "utilities/GasParticle.hpp"  // for GasParticle
#include <algorithm>                  // for max, min
#include <cmath>                      // for sqrt
#include <cstddef>                    // for NULL
using namespace std;

/**
 * @brief Rotate the velocity components of the given StateVector to a frame
 * where the x-axis coincides with the normal of the interface
 *
 * The y- and z-axis are chosen as two mutually orthogonal directions in the
 * plane of the interface.
 *
 * @param W StateVector to transform
 */
void AdaptiveFace3d::transform(StateVector& W) {
    double W1acc = W[1] * _cosp + W[2] * _sinp;
    double W1 = W1acc * _sint + W[3] * _cost;
    double W2 = -W[1] * _sinp + W[2] * _cosp;
    W[3] = -W1acc * _cost + W[3] * _sint;
    W[1] = W1;
    W[2] = W2;
}

/**
 * @brief The inverse of AdaptiveFace3d::transform()
 *
 * If W' = transform(W), then invtransform(W') = W
 *
 * @param W StateVector to transform
 */
void AdaptiveFace3d::invtransform(StateVector& W) {
    double W1acc = W[1] * _sint - W[3] * _cost;
    double W1 = W1acc * _cosp - W[2] * _sinp;
    double W2 = W1acc * _sinp + W[2] * _cosp;
    W[3] = W[1] * _cost + W[3] * _sint;
    W[1] = W1;
    W[2] = W2;
}

/**
 * @brief Get the components of the normal to the face
 *
 * If \f$\theta\f$ is the angle between the normal and the z-axis and \f$\phi\f$
 * the angle between the projection of the normal in the xy-plane and the
 * x-axis, then the components are
 * \f[
 * (\sin(\theta)\cos(\phi), \sin(\theta)\sin(\phi), \cos(\theta))
 * \f]
 *
 * @param angles Array to store the components of the normal vector in
 */
void AdaptiveFace3d::get_normal(double* angles) {
    angles[0] = _sint * _cosp;
    angles[1] = _sint * _sinp;
    angles[2] = _cost;
}

/**
 * @brief Constructor.
 *
 * Unlike for the old algorithm, computation of the geometrical quantities
 * associated to the face is done in a separate class. Here, we just pass them
 * on as parameters to the constructor.
 *
 * The only computation done inside the constructor is that of the angles needed
 * to transform to a reference frame where the x-axis coincides with the
 * interface normal. The y- and z-axis are chosen as arbitrary directions that
 * are mutually orthogonal and lie in the plane of the interface.
 *
 * @param left GasParticle at the left of the interface
 * @param right GasParticle at the right of the interface
 * @param rightpos Position of the actual representative of the GasParticle at
 * the right of the interface
 * @param midpoint Midpoint of the interface
 * @param area Area of the interface
 */
AdaptiveFace3d::AdaptiveFace3d(GasParticle* left, GasParticle* right,
                               Vec& rightpos, double* midpoint, double area)
        : _rightpos(rightpos) {

    _left = left;
    _right = right;
    _midpoint[0] = midpoint[0];
    _midpoint[1] = midpoint[1];
    _midpoint[2] = midpoint[2];
    _v[0] = 0.;
    _v[1] = 0.;
    _v[2] = 0.;
    _area = area;

    double d[4];
    d[0] = _left->x() - _rightpos.x();
    d[1] = _left->y() - _rightpos.y();
    d[2] = _left->z() - _rightpos.z();
    d[3] = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    _cost = -d[2] / d[3];
    _sint = sqrt(1. - _cost * _cost);
    if(_sint) {
        _sinp = -d[1] / d[3] / _sint;
        _cosp = -d[0] / d[3] / _sint;
    } else {
        _sinp = 0.;
        _cosp = 1.;
    }
}

/**
 * @brief Set the velocities of the interface based on the velocities of its
 * left and right state
 */
void AdaptiveFace3d::set_v() {
    double rRL[4];
    rRL[0] = _rightpos.x() - _left->x();
    rRL[1] = _rightpos.y() - _left->y();
    rRL[2] = _rightpos.z() - _left->z();
    rRL[3] = rRL[0] * rRL[0] + rRL[1] * rRL[1] + rRL[2] * rRL[2];
    double fac = (_left->vx() - _right->vx()) *
                 (_midpoint[0] - 0.5 * (_left->x() + _rightpos.x()));
    fac += (_left->vy() - _right->vy()) *
           (_midpoint[1] - 0.5 * (_left->y() + _rightpos.y()));
    fac += (_left->vz() - _right->vz()) *
           (_midpoint[2] - 0.5 * (_left->z() + _rightpos.z()));
    fac /= rRL[3];
    _v[0] = 0.5 * (_left->vx() + _right->vx()) + fac * rRL[0];
    _v[1] = 0.5 * (_left->vy() + _right->vy()) + fac * rRL[1];
    _v[2] = 0.5 * (_left->vz() + _right->vz()) + fac * rRL[2];
}

/**
 * @brief Calculate the flux through the interface, based on the solution of the
 * Riemann problem at the interface
 *
 * We first interpolate and predict the left and right states (\f$W_L\f$ and
 * \f$W_R\f$) from the centroids of the left and right cells to the midpoint of
 * the interface, to obtain the states \f$W'_L\f$ and \f$W'_R\f$. We then rotate
 * and boost these states to the frame connected to the interface to obtain the
 * states \f$W''_L\f$ and \f$W''_R\f$. These states then serve as an input for a
 * given Riemann solver, which returns the state at the interface, \f$W''*\f$.
 * After this StateVector is deboosted and transformed back to the original
 * reference frame, it is ultimately used to retrieve fluxes for the conserved
 * quantities from the Riemann solver, which are in turn used to update the
 * left and right state GasParticle.
 *
 * @param timeline TimeLine of the simulation, used to convert from integer time
 * to real time
 * @param solver Solver used to solve the Riemann problem at the interface and
 * to convert a StateVector to fluxes
 */
void AdaptiveFace3d::calculate_flux(TimeLine& timeline, RiemannSolver& solver) {
    if(!_area) {
        return;
    }
    if(_right == NULL && _left->get_starttime() != timeline.get_integertime()) {
        return;
    }
    if(_right != NULL &&
       _right->get_starttime() != timeline.get_integertime() &&
       _left->get_starttime() != timeline.get_integertime()) {
        return;
    }
    Vec v(_v[0], _v[1], _v[2]);

    double dt = timeline.get_realtime_interval(_left->get_timestep());
    if(_right != NULL) {
        dt = std::min(dt,
                      timeline.get_realtime_interval(_right->get_timestep()));
    }

    StateVector WL = _left->get_Wvec();
    StateVector WR;
    if(_right != NULL) {
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
            if(_right == NULL) {
                continue;
            }
            Wtemp = WR;
            cell = _right;
            pos = _rightpos;
        }
        double d[ndim_];
        Vec centroid = cell->get_centroid();
        for(unsigned int n = ndim_; n--;) {
            double correction = cell->pos(n) - pos[n];
            d[n] = _midpoint[n] - centroid[n] + correction;
        }
        double dtp = dt;
        for(unsigned int m = ndim_ + 2; m--;) {
            Wtemp[m] += cell->get_gradient(m, 0) * d[0] +
                        cell->get_gradient(m, 1) * d[1] +
                        cell->get_gradient(m, 2) * d[2];
        }
        double Ax[5][5] = {
                {Wtemp[1], Wtemp[0], 0., 0., 0.},
                {0., Wtemp[1], 0., 0., 1. / Wtemp[0]},
                {0., 0., Wtemp[1], 0., 0.},
                {0., 0., 0., Wtemp[1], 0.},
                {0., solver.get_gamma() * Wtemp[4], 0., 0., Wtemp[1]}};
        double Ay[5][5] = {
                {Wtemp[2], Wtemp[0], 0., 0., 0.},
                {0., Wtemp[2], 0., 0., 0.},
                {0., 0., Wtemp[2], 0., 1. / Wtemp[0]},
                {0., 0., 0., Wtemp[2], 0.},
                {0., 0., solver.get_gamma() * Wtemp[4], 0., Wtemp[2]}};
        double Az[5][5] = {
                {Wtemp[3], Wtemp[0], 0., 0., 0.},
                {0., Wtemp[3], 0., 0., 0.},
                {0., 0., Wtemp[3], 0., 0.},
                {0., 0., 0., Wtemp[3], 1. / Wtemp[0]},
                {0., 0., 0., solver.get_gamma() * Wtemp[4], Wtemp[3]}};
        StateVector deltaW;
        for(unsigned int m = ndim_ + 2; m--;) {
            for(unsigned int n = ndim_ + 2; n--;) {
                deltaW[m] -= 0.5 * dtp * Ax[m][n] * cell->get_gradient(n, 0);
                deltaW[m] -= 0.5 * dtp * Ay[m][n] * cell->get_gradient(n, 1);
                deltaW[m] -= 0.5 * dtp * Az[m][n] * cell->get_gradient(n, 2);
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
    if(_right == NULL) {
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
    get_normal(angles);
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
}
