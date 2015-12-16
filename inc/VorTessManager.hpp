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
 * @file VorTessManager.hpp
 *
 * @brief General interface for Voronoi tesselations: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_VORTESSMANAGER
#define HEAD_VORTESSMANAGER

class VorTess;
class DelCont;
class GasParticle;
class ParticleVector;
class AdaptiveVorTess;
class AdaptiveVorTess3d;
class FixedGrid;
class RestartFile;

#include "utilities/Tree.hpp"
#include "utilities/Timer.hpp"
#include "TimeLine.hpp"
#include "RiemannSolver.hpp"

/**
 * @brief General Voronoi tesselation that is accessed by Simulation
 *
 * Can either be a classic VorTess that is reconstructed every timestep or an
 * AdaptiveVorTess/AdaptiveVorTess3d that is evolved. Which one is the
 * underlying Voronoi tesselation is set by an internal flag. Currently, only
 * the old VorTess is fully functional.
 *
 * There is also a FixedGrid option, which can be used to test hydro specific
 * stuff without having to bother with the slow Voronoi grid, but this only
 * works in ideal cases.
 */
class VorTessManager{
private:
    /*! \brief Pointer to a VorTess, only used if this tesselation uses the old
     *  algorithm */
    VorTess* _tesselation;
#ifdef FIXED_GRID
    /*! \brief Pointer to a FixedGrid, used for the fixed Cartesian grid */
    FixedGrid* _fixedgrid;
#endif

    /*! \brief ParticleVector for which the tesselation is constructed */
    ParticleVector& _particles;

    /*! \brief Timer to quantify time spent in grid specific operations */
    Timer _gridtimer;
    /*! \brief Timer to quantify time spent in hydro specific operations */
    Timer _hydrotimer;

    /*! \brief Flag indicating if the simulation box is periodic (true) or
     *  reflective (false) */
    bool _periodic;
    /*! \brief Tolerance value used to distinguish between standard floating
     *  point arithmetics and arbitrary precision arithmetics in geometry
     *  tests */
    double _tolerance;

public:
    VorTessManager(ParticleVector& particles, bool periodic = false,
                   double tolerance = 1.e-9, bool evolve = true);
    ~VorTessManager();

    void reset(DelCont* cont, bool periodic = false, double tolerance = 1.e-9);

    unsigned int update(unsigned long currentTime);
    void update_Ws();
    void set_hs(unsigned long currentTime);
    void update_gradients();
    void update_dts(unsigned long currentTime);
    void update_gravitational_corrections();
    void hydro(TimeLine& timeline, RiemannSolver& solver,
               ParticleVector &particles);
    void update_dQs();

    void print_tesselation_gnuplot(std::ostream &stream);

    void dump(RestartFile &rfile);
    VorTessManager(RestartFile &rfile, ParticleVector &particles);

    Vec get_velocity(unsigned int index);
    void estimate_gradients(unsigned int index, StateVector* delta);
    double get_volume(unsigned int index);
    double get_total_area(unsigned int index);
    Vec get_centroid(unsigned int index);

    void update_positions(bool periodic);

    Vec get_gravitational_correction(unsigned int index);

    void dump_connectivity(std::ostream& stream, unsigned long currentTime = 0);
};

#endif
