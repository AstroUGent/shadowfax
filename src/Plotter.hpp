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
 * @file Plotter.hpp
 *
 * @brief SideProgram to plot Shadowfax snapshots: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PLOTTER_HPP
#define PLOTTER_HPP

#include <string>
#include <vector>

/**
 * @brief SideProgram used to plot snapshots in PPM format
 *
 * All capabilities are now also supported by using VisIt and the vtkmaker.
 */
class Plotter {
  private:
    std::string get_filename(std::string filename);

  public:
    Plotter(int argc, char** argv);
    ~Plotter() {}
};

#endif  // PLOTTER_HPP
