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
 * @file MPIGlobal.hpp
 *
 * @brief Global MPI variables: size, local rank and communication buffer
 *
 * This header only defines the variables, they have to be initialized
 * elsewhere.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef MPIGLOBAL_HPP
#define MPIGLOBAL_HPP

#include "utilities/Timer.hpp"
#include <mpi.h>

/**
 * \brief Global MPI variables used throughout the program
 */
namespace MPIGlobal {
/*! \brief Number of MPI processes in the system */
extern int size;
/*! \brief Rank of the local MPI process in the system */
extern int rank;

/*! \brief Number of MPI processes on the local node */
extern int local_size;
/*! \brief Rank of the local MPI process on the local node */
extern int local_rank;
/*! \brief Size of the space of all nodes */
extern int nodesize;
/*! \brief Rank of the node in the space of all nodes */
extern int noderank;
/*! \brief Local node communicator */
extern MPI_Comm nodecomm;

/*! \brief Buffer used to send data to other MPI processes */
extern char* sendbuffer;
/*! \brief Size of the send buffer in bytes */
extern unsigned int sendsize;
/*! \brief Buffer used to receive data from other MPI processes */
extern char* recvbuffer;
/*! \brief Size of the receive buffer in bytes */
extern unsigned int recvsize;

/*! \brief Timer to quantify idle time spent waiting on other MPI
 *  processes */
extern Timer idletimer;
/*! \brief Timer to quantify time spent during MPI communication
 *  operations */
extern Timer commtimer;
}

#endif  // MPIGLOBAL_HPP
