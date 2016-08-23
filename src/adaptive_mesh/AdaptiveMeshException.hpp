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
 * @file AdaptiveMeshException.hpp
 *
 * @brief Exceptions thrown by the mesh evolution algorithm
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ADAPTIVEMESHEXCEPTION_HPP
#define ADAPTIVEMESHEXCEPTION_HPP

#include <exception>

/**
 * @brief Exception thrown when something goes wrong during the evolution of the
 * mesh
 */
class AdaptiveMeshException : public std::exception {
    /**
     * @brief Return a human-readable explanation of what is going wrong
     *
     * @return The name of the exception, because that's all we know
     */
    virtual const char* what() const throw() {
        return "AdaptiveMeshException";
    }
};

#endif
