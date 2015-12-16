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
 * @file VorFace.hpp
 *
 * @brief Voronoi face: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_VORFACE
#define HEAD_VORFACE

#include "Vec.hpp"
#include "StateVector.hpp"
#include <vector>
#include <ostream>
#include <exception>

class VorGen;
class VorCell;
class StateVector;
class TimeLine;
class RiemannSolver;
class ParticleVector;

/**
 * @brief Voronoi face in 2D or 3D
 *
 * Used for the old algorithm. Represents a Voronoi face during grid
 * construction and performs the flux exchange during the hydro step.
 * This is different from the mesh evolution algorithm, where these steps
 * are represented by separate classes.
 */
class VorFace{
private:
    /*! \brief Vertices that make up the face */
    std::vector<VorGen*> _vertices;

    /*! \brief Right generator of the face */
    VorGen* _right;

    /*! \brief Left generator of the face */
    VorGen* _left;

    /*! \brief Index of the right generator in the DelTess VorGen list */
    unsigned int _vright;

    /*! \brief Index of the left generator in the DelTess VorGen list */
    unsigned int _vleft;

    /*! \brief Face neighbours, necessary for evolution algorithm
     *  initialization */
    std::vector<VorGen*> _ngbs;

    /*! \brief Geometrical area of the face */
    double _area;

    /*! \brief Geometrical centroid of the face */
    Vec _midpoint;

    /*! \brief Velocity of the face */
    Vec _v;

    /**
     * @brief Per face slope limiter for a single quantity
     *
     * Based on the slope limiter described in one of the appendices of Hopkins'
     * GIZMO paper.
     *
     * @param phimid0 Reconstructed value of the quantity at the interface
     * @param phiL Value at the left of the interface
     * @param phiR Value at the right of the interface
     * @param d Coordinates of the vector pointing from the centroid of the left
     * cell to the midpoint of the interface
     * @param r Distance between the centroids of the left and right cell
     * @return Limited value of the quantity at the interface
     */
    static double limit(const double phimid0, const double phiL,
                        const double phiR, const double* d, const double r){
        double psi1 = 0.5;
        double psi2 = 0.25;
#if ndim_==3
        double dnrm = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] );
#else
        double dnrm = sqrt( d[0] * d[0] + d[1] * d[1] );
#endif

        double delta1 = psi1 * fabs( phiL - phiR );
        double delta2 = psi2 * fabs( phiL - phiR );

        double phimin = fmin( phiL, phiR );
        double phimax = fmax( phiL, phiR );

        double phibar = phiL + dnrm / r * ( phiR - phiL );

        /* if sign(phimax+delta1) == sign(phimax) */
        double phiplus;
        if( ( phimax + delta1 ) * phimax > 0. ){
            phiplus = phimax + delta1;
        } else {
            phiplus = phimax / ( 1. + delta1 / fabs(phimax) );
        }

        /* if sign(phimin-delta1) == sign(phimin) */
        double phiminus;
        if( ( phimin - delta1 ) * phimin > 0. ){
            phiminus = phimin - delta1;
        } else {
            phiminus = phimin / ( 1. + delta1 / fabs(phimin) );
        }

        double phimid;
        if( phiL == phiR ){
            phimid = phiL;
        } else {
            if( phiL < phiR ){
                phimid = fmax( phiminus, fmin( phibar + delta2, phimid0 ) );
            } else {
                phimid = fmin( phiplus, fmax( phibar - delta2, phimid0 ) );
            }
        }
        return phimid;
    }

    /**
     * @brief Per face limiter for a StateVector
     *
     * @param W Reconstructed StateVector at the interface
     * @param WL StateVector of the left state
     * @param WR StateVector of the right state
     * @param d Coordinates of the vector pointing from the centroid of the left
     * cell to the midpoint of the interface
     * @param r Distance between the centroids of the left and right cell
     * @return Limited StateVector at the interface
     */
    static StateVector limit(const StateVector& W, const StateVector& WL,
                             const StateVector& WR, const double* d,
                             const double r){
        StateVector Wtemp;
        Wtemp[0] = limit(W[0], WL[0], WR[0], d, r);
        Wtemp[1] = limit(W[1], WL[1], WR[1], d, r);
        Wtemp[2] = limit(W[2], WL[2], WR[2], d, r);
        Wtemp[3] = limit(W[3], WL[3], WR[3], d, r);
#if ndim_==3
        Wtemp[4] = limit(W[4], WL[4], WR[4], d, r);
#endif
        return Wtemp;
    }

    /**
     * @brief Reconstruct the given StateVector using the given gradients and
     * distance between centroid and interface
     *
     * @param W StateVector to reconstruct
     * @param gradW Gradients of the hydrodynamical quantities in the cell
     * @param d Coordinates of the vector pointing from the centroid of the left
     * cell to the midpoint of the interface
     * @return Reconstructed StateVector at the interface
     */
    static StateVector reconstruct(const StateVector& W,
                                   const StateVector* gradW, const double* d){
        StateVector Wtemp(W);
#if ndim_==3
        Wtemp[0] += gradW[0][0]*d[0] + gradW[1][0]*d[1] + gradW[2][0]*d[2];
        Wtemp[1] += gradW[0][1]*d[0] + gradW[1][1]*d[1] + gradW[2][1]*d[2];
        Wtemp[2] += gradW[0][2]*d[0] + gradW[1][2]*d[1] + gradW[2][2]*d[2];
        Wtemp[3] += gradW[0][3]*d[0] + gradW[1][3]*d[1] + gradW[2][3]*d[2];
        Wtemp[4] += gradW[0][4]*d[0] + gradW[1][4]*d[1] + gradW[2][4]*d[2];

#else
        Wtemp[0] += gradW[0][0]*d[0] + gradW[1][0]*d[1];
        Wtemp[1] += gradW[0][1]*d[0] + gradW[1][1]*d[1];
        Wtemp[2] += gradW[0][2]*d[0] + gradW[1][2]*d[1];
        Wtemp[3] += gradW[0][3]*d[0] + gradW[1][3]*d[1];
#endif
        return Wtemp;
    }

public:
    VorFace(unsigned int index, std::vector<VorGen*>& points);
    ~VorFace(){}

    void print(std::ostream &stream);
    void print_pov_frame(std::ostream &stream);
    void print_pov(std::ostream &stream);
    void print_leaflet(std::ostream &vstream, int ox, int oy, VorGen* center);
    bool overlap(double* box);

    void add_vertex(VorGen* point);
    std::vector<VorGen*> get_vertices();

    void add_facengb(VorGen* ngb);
    std::vector<VorGen*> get_facengbs();

    void add_ngb(VorGen* ngb);
    VorGen* get_ngb();
    VorGen* get_left();

    void add_ngb_id(unsigned int ngb);
    unsigned int get_ngb_id();

    Vec& get_midpoint();
    void set_midpoint(double* midpoint);
    void calculate_midpoint();

    double get_area();
    void set_area(double area);
    void calculate_area();

    void transform(StateVector& W, const Vec &left, const Vec &right);
    void invtransform(StateVector& W, const Vec& left, const Vec& right);
    void get_normal(double* angles);

    Vec& get_v();
    void set_v(ParticleVector& particles);

#ifndef ICMAKER
    void calculate_flux(TimeLine& timeline, RiemannSolver& solver);
#endif
    void calculate_advection_flux(double dt);

    unsigned int get_triangles(std::vector<float>& positions,
                               std::vector<int>& connectivity);
};

/**
 * @brief Exception thrown when something goes wrong during the flux calculation
 */
class FluxCalculationException : public std::exception{
private:
    /*! \brief Left StateVector */
    StateVector _WL;

    /*! \brief Right StateVector */
    StateVector _WR;

    /*! \brief Rank of the MPI process */
    unsigned int _rank;

    /*! \brief Reference to the Solver used to solve the Riemann problem */
    RiemannSolver& _solver;

public:
    FluxCalculationException(StateVector WL, StateVector WR, unsigned int rank,
                             RiemannSolver& solver);

    virtual const char* what() const throw();
};

#endif
