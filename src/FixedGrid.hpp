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
 * @file FixedGrid.hpp
 *
 * @brief A fixed cartesian grid: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef FIXEDGRID_HPP
#define FIXEDGRID_HPP

#include <vector>

class ParticleVector;
class TimeLine;
class RiemannSolver;
class VorFace;
class VorGen;
class VorCell;
class StateVector;

/**
 * @brief 2D or 3D fixed grid
 *
 * Used to speed up the calculations for a static Cartesian mesh. Does only
 * work for specific initial particle configurations and should only be used
 * to test changes to the hydro solver in controlled test simulations.
 */
class FixedGrid {
  private:
    /*! @brief Number of grid cells in every dimension */
    unsigned int _size[ndim_];

    /*! @brief Volume of a single cell of the grid */
    double _cellvolume;

    /*! @brief Total surface area of a single cell of the grid */
    double _cellarea;

    /*! @brief Characteristic length of a single cell of the grid */
    double _cellh;

    /*! @brief Faces of the grid */
    std::vector<VorFace*> _faces;

    /*! @brief Generators of the grid cells */
    std::vector<VorGen*> _vorgens;

    /*! @brief Ghosts at the boundaries of the grid */
    std::vector<VorGen*> _ghosts;

    /*! @brief Cells of the grid */
    std::vector<VorCell*> _cells;

    /*! @brief Indices of the associated gasparticles in the internal grid */
    std::vector<unsigned int> _idx;

  public:
    FixedGrid(ParticleVector& particles, bool periodic = false);
    ~FixedGrid();

    double get_h();
    double get_volume();
    double get_area();

    void hydro(TimeLine& timeline, RiemannSolver& solver);
    void get_gradients(unsigned int index, StateVector* delta);
};

#endif  // FIXEDGRID_HPP
