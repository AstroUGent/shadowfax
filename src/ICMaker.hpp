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
 * @file ICMaker.hpp
 *
 * @brief Sideprogram to generate initial conditions: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ICMAKER_H
#define ICMAKER_H

/**
  * @brief Driving class for initial condition generating part of the code
  *
  * This class can be thought of as a stand-alone program which generates
  * initial conditions. The initial conditions are generated using a
  * BlockICGenerator which takes it input from an xml file that is provided as a
  * parameter to the program.
  *
  * Alternatively, a SpecificICGenerator can be used to generate hard-coded
  * initial conditions.
  */
class ICMaker {
  public:
    ICMaker(int argc, char** argv);
    ~ICMaker() {}
};

#endif  // ICMAKER_H
