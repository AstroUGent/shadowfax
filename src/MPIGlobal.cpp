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
 * @file MPIGlobal.cpp
 *
 * @brief Declaration of global MPI variables
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "MPIGlobal.hpp"

int MPIGlobal::rank;
int MPIGlobal::size;

int MPIGlobal::local_rank;
int MPIGlobal::local_size;
int MPIGlobal::nodesize;
int MPIGlobal::noderank;
MPI_Comm MPIGlobal::nodecomm;

char* MPIGlobal::sendbuffer;
unsigned int MPIGlobal::sendsize;
char* MPIGlobal::recvbuffer;
unsigned int MPIGlobal::recvsize;

Timer MPIGlobal::idletimer;
Timer MPIGlobal::commtimer;
