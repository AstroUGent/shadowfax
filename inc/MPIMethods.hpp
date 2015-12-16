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
 * @file MPIMethods.hpp
 *
 * @brief Error-checking wrappers around MPI methods
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef MPIMETHODS_HPP
#define MPIMETHODS_HPP

#include "Error.hpp"

/**
 * @brief Wrapper around MPI_Comm_rank
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param rank Pointer to the integer to store the result in
 */
inline void MyMPI_Comm_rank(int *rank){
    int status = MPI_Comm_rank(MPI_COMM_WORLD, rank);
    if(status != MPI_SUCCESS){
        std::cerr << "Failed to obtain process rank!" << std::endl;
        my_exit();
    }
}

/**
 * @brief Wrapper around MPI_Comm_size
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param size Pointer to the integer to store the result in
 */
inline void MyMPI_Comm_size(int *size){
    int status = MPI_Comm_size(MPI_COMM_WORLD, size);
    if(status != MPI_SUCCESS){
        std::cerr << "Failed to obtain MPI world size!" << std::endl;
        my_exit();
    }
}

/**
 * @brief Convert a string to a relatively unique integer value
 *
 * There is no guarantee that the resulting integer will be different for all
 * possible input strings, but it should do a fairly good job.
 *
 * We use this function to provide every MPI node with a different integer value
 * that can then be used to split the MPI communicator into different
 * communicators for every MPI node.
 *
 * We currently use the Adler-32 checksum algorithm, found on Wikipedia.
 *
 * @param name Name of a MPI node
 * @return Semi-unique integer hash of the name
 */
inline int MyMPI_hash(std::string name){
    int a = 1, b = 0;

    /* Process each byte of the data in order */
    for(unsigned int index = 0; index < name.size(); ++index){
        a = (a + name[index]) % 65521;
        b = (b + a) % 65521;
    }

    return (b << 16) | a;
}

/**
 * @brief Wrapper around MPI_Get_processor_name
 *
 * Can be used to optionally emulate multinode behaviour on a single node, in
 * which case two separate nodes are emulated: an even node and an odd node,
 * containing respectively the processes with even and odd ranks.
 *
 * @param name Buffer to store the node name in
 * @param resultlen Length of the node name
 * @return MPI error code
 */
inline int MyMPI_Get_processor_name(char *name, int *resultlen){

//#define EMULATE_MULTINODES

#ifdef EMULATE_MULTINODES
    if(MPIGlobal::rank%2){
        sprintf(name, "%s", "even");
        *resultlen = 4;
    } else {
        sprintf(name, "%s", "odd");
        *resultlen = 3;
    }

    return MPI_SUCCESS;
#else
    return MPI_Get_processor_name(name, resultlen);
#endif
}

/**
 * @brief Method to determine the rank and size on the local node
 *
 * If the program is running on multiple nodes, it can be useful to know the
 * rank and size of the part of the program running on the local node, e.g. to
 * write output files per node.
 *
 * MPI does not provide standard functions to do this, so we have to do it
 * ourselves. The method below is based on a suggestion by Markus Wittmann on
 * his blog
 * (https://blogs.fau.de/wittmann/2013/02/mpi-node-local-rank-determination/).
 *
 * We use MPI_Comm_split to create a new communicator based on a color which is
 * derived from the processor name, obtained from MPI_Get_processor_name. The
 * local rank and size are then obtained by calling MPI_Comm_rank and
 * MPI_Comm_size on this new communicator.
 *
 * Because it is impossible to define a completely unique hash function to
 * convert node names to integer colors, we do a collective gather on the
 * newly created communicator for the node names. We then check if the names
 * are really different. If not, we do a second communicator split, where this
 * time we use the index of the own node name in the node name list as a color.
 * Since this second collective operation is done on a much smaller
 * communicator (which should be the local node for most cases), it is much
 * faster than using the same method at the start.
 *
 * @param rank Variable to store the local rank of this process in
 * @param size Variable to store the local size of the node in
 * @param noderank Rank of the node in the space of all nodes
 * @param nodesize Size of the space of all nodes
 * @param nodecomm Variable to store the local node communicator in
 */
inline void MyMPI_Comm_local_vars(int *rank, int *size, int *noderank,
                                  int *nodesize, MPI_Comm *nodecomm){
    MPIGlobal::commtimer.start();

    char *name = new char[MPI_MAX_PROCESSOR_NAME];
    int namesize;
    int status = MyMPI_Get_processor_name(name, &namesize);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while obtaining MPI processor name!" << std::endl;
        my_exit();
    }

    name[namesize] = '\0';
    std::string procname(name);

    int color = MyMPI_hash(procname);
    MPI_Comm local_world = MPI_COMM_NULL;
    status = MPI_Comm_split(MPI_COMM_WORLD, color, MPIGlobal::rank,
                            &local_world);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while splitting MPI communicator in local MPI "
                     "communicators!" << std::endl;
        my_exit();
    }

    status = MPI_Comm_rank(local_world, rank);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while obtaining local MPI rank!" << std::endl;
        my_exit();
    }
    status = MPI_Comm_size(local_world, size);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while obtaining local MPI size!" << std::endl;
        my_exit();
    }

    // check if the processes in local_world are really on a single node
    char *allnames = new char[(*size)*MPI_MAX_PROCESSOR_NAME];
    status = MPI_Allgather(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, allnames,
                           MPI_MAX_PROCESSOR_NAME, MPI_CHAR, local_world);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while gathering processor names on the presumably "
                     "local communicator!" << std::endl;
        my_exit();
    }

    int owni = 0;
    bool is_single = true;
    for(int i = 0; i < (*size); i++){
        std::string othername(&allnames[i*MPI_MAX_PROCESSOR_NAME]);
        if(othername != procname){
            is_single = false;
        } else {
            owni = i;
        }
    }
    if(!is_single){
        MPI_Comm real_local_world = MPI_COMM_NULL;
        status = MPI_Comm_split(local_world, owni, (*rank), &real_local_world);
        if(status != MPI_SUCCESS){
            std::cerr << "Error while splitting the local communicator in the "
                         "real node communicators!" << std::endl;
            my_exit();
        }

        status = MPI_Comm_rank(real_local_world, rank);
        if(status != MPI_SUCCESS){
            std::cerr << "Error while obtaining real local MPI rank!"
                      << std::endl;
            my_exit();
        }

        status = MPI_Comm_size(real_local_world, size);
        if(status != MPI_SUCCESS){
            std::cerr << "Error while obtaining real local MPI size!"
                      << std::endl;
            my_exit();
        }

        status = MPI_Comm_free(&real_local_world);
        if(status != MPI_SUCCESS){
            std::cerr << "Error while freeing real local MPI communicator!"
                      << std::endl;
            my_exit();
        }

        status = MPI_Comm_dup(real_local_world, nodecomm);
        if(status != MPI_SUCCESS){
            std::cerr << "Error while duplicating local MPI communicator!"
                      << std::endl;
            my_exit();
        }
    } else {
        status = MPI_Comm_dup(local_world, nodecomm);
        if(status != MPI_SUCCESS){
            std::cerr << "Error while duplicating local MPI communicator!"
                      << std::endl;
            my_exit();
        }
    }

    status = MPI_Comm_free(&local_world);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while freeing local MPI communicator!" << std::endl;
        my_exit();
    }

    delete [] name;

    // obtain noderank and nodesize
    MPI_Comm nodeworld;
    int nodemaster = (*rank == 0);
    status = MPI_Comm_split(MPI_COMM_WORLD, nodemaster, MPIGlobal::rank,
                            &nodeworld);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while obtaining node world!" << std::endl;
        my_exit();
    }
    if(*rank == 0){
        status = MPI_Comm_rank(nodeworld, noderank);
        if(status != MPI_SUCCESS){
            std::cerr << "Error while obtaining noderank!" << std::endl;
            my_exit();
        }

        status = MPI_Comm_size(nodeworld, nodesize);
        if(status != MPI_SUCCESS){
            std::cerr << "Error while obtaining nodesize!" << std::endl;
            my_exit();
        }
    }

    status = MPI_Comm_free(&nodeworld);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while freeing nodeworld!" << std::endl;
        my_exit();
    }

    status = MPI_Bcast(noderank, 1, MPI_INT, 0, *nodecomm);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while broadcasting noderank!" << std::endl;
        my_exit();
    }

    status = MPI_Bcast(nodesize, 1, MPI_INT, 0, *nodecomm);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while broadcasting nodesize!" << std::endl;
        my_exit();
    }

    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Allreduce
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param sendbuf Buffer that is being sent
 * @param recvbuf Buffer to receive in
 * @param count Number of elements to be sent
 * @param datatype MPI datatype of the elements
 * @param op Global reduce operation
 */
inline void MyMPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                           MPI_Datatype datatype, MPI_Op op){
    MPIGlobal::commtimer.start();
    int status = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op,
                               MPI_COMM_WORLD);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Allreduce!" << std::endl;
        my_exit();
    }
    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Bcast
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param buffer Buffer to broadcast from root to the other MPI processes
 * @param count Number of elements in the buffer
 * @param datatype MPI datatype of the elements
 * @param root Root process
 * @param comm Communicator identifying the participating processes (default:
 * MPI_COMM_WORLD)
 */
inline void MyMPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
                        int root, MPI_Comm comm = MPI_COMM_WORLD){
    MPIGlobal::commtimer.start();
    int status = MPI_Bcast(buffer, count, datatype, root, comm);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Bcast!" << std::endl;
        my_exit();
    }
    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Barrier
 *
 * We check the error code. However, the only error code that can occur is
 * MPI_ERR_COMM for an invalid communicator and since we use the default
 * communicator MPI_WORLD, this should never happen.
 *
 * @param comm Communicator identifying the participating processes (default:
 * MPI_COMM_WORLD)
 */
inline void MyMPI_Barrier(MPI_Comm comm = MPI_COMM_WORLD){
    MPIGlobal::idletimer.start();
    int status = MPI_Barrier(comm);
    if(status != MPI_SUCCESS){
        std::cerr << "Error in MPI_Barrier!" << std::endl;
        my_exit();
    }
    MPIGlobal::idletimer.stop();
}

/**
 * @brief Wrapper around MPI_Scan
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param sendbuf Buffer that is being sent
 * @param recvbuf Buffer to receive in
 * @param count Number of elements to be sent
 * @param datatype MPI datatype of the elements
 * @param op Global reduce operation
 */
inline void MyMPI_Scan(void *sendbuf, void *recvbuf, int count,
                       MPI_Datatype datatype, MPI_Op op){
    MPIGlobal::commtimer.start();
    int status = MPI_Scan(sendbuf, recvbuf, count, datatype, op,
                          MPI_COMM_WORLD);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Scan!" << std::endl;
        my_exit();
    }
    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Isend
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param buf Buffer that is being sent
 * @param count Number of elements in the buffer
 * @param datatype MPI datatype of the elements
 * @param dest Rank of the destination process
 * @param tag Tag identifying this message
 * @param request Request to store status information in
 */
inline void MyMPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest,
                      int tag, MPI_Request *request){
    MPIGlobal::commtimer.start();
    int status = MPI_Isend(buf, count, datatype, dest, tag, MPI_COMM_WORLD,
                           request);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Isend!" << std::endl;
        my_exit();
    }
    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Recv
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param buf Buffer to receive in
 * @param count Number of elements to receive
 * @param datatype MPI datatype of the elements
 * @param source Rank of the source process
 * @param tag Tag identifying this message
 * @param status Status information about the reception
 */
inline void MyMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source,
                     int tag, MPI_Status *status){
    MPIGlobal::commtimer.start();
    int result = MPI_Recv(buf, count, datatype, source, tag, MPI_COMM_WORLD,
                          status);
    if(result != MPI_SUCCESS){
        std::cerr << "Error during MPI_Recv!" << std::endl;
        my_exit();
    }
    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Get_count
 *
 * We check the error code to detect MPI errors.
 *
 * @param status MPI status object to query
 * @param datatype MPI datatype of the elements that were received
 * @param count Pointer to an integer to store the result in
 */
inline void MyMPI_Get_count(MPI_Status *status, MPI_Datatype datatype,
                            int *count){
    int result = MPI_Get_count(status, datatype, count);
    if(result != MPI_SUCCESS){
        std::cerr << "Error in MPI_Get_count!" << std::endl;
        my_exit();
    }
}

/**
 * @brief Wrapper around MPI_Probe
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param source Rank of the source process
 * @param tag Tag identifying the message
 * @param status MPI status object to store information in
 */
inline void MyMPI_Probe(int source, int tag, MPI_Status *status){
    MPIGlobal::idletimer.start();
    int result = MPI_Probe(source, tag, MPI_COMM_WORLD, status);
    if(result != MPI_SUCCESS){
        std::cerr << "Error in MPI_Probe!" << std::endl;
        my_exit();
    }
    MPIGlobal::idletimer.stop();
}

/**
 * @brief Wrapper around MPI_Pack
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD. We also keep track of the calling function to print out relevant
 * error messages before aborting. This is done by wrapping the function call
 * using the MyMPI_Pack macro.
 *
 * @param inbuf Buffer with elements to pack
 * @param incount Number of elements to be packed
 * @param datatype MPI datatype of the elements
 * @param outbuf Buffer to pack in
 * @param outsize Total size of the pack buffer
 * @param position Current position in the pack buffer (is updated)
 * @param file Filename of the source code file that calls the function
 * @param function Function name of the calling function
 * @param line Line number in the source code file where the call happens
 */
inline void MyMPI_Pack_real(void *inbuf, int incount, MPI_Datatype datatype,
                            void *outbuf, int outsize, int *position,
                            std::string file, std::string function, int line){
    int status = MPI_Pack(inbuf, incount, datatype, outbuf, outsize, position,
                          MPI_COMM_WORLD);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Pack!" << std::endl;
        std::cerr << MPIGlobal::rank << ": " << file << "::" << function << ":"
                  << line << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

#define MyMPI_Pack(buf, count, type, out, outs, pos) {\
    MyMPI_Pack_real(buf, count, type, out, outs, pos, __FILE__, __FUNCTION__,\
    __LINE__);\
    }

/**
 * @brief Wrapper around MPI_Unpack
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param inbuf Buffer to unpack from
 * @param insize Total size of the unpack buffer
 * @param position Current position in the unpack buffer (is updated)
 * @param outbuf Buffer to unpack elements in
 * @param outcount Number of elements to unpack
 * @param datatype MPI datatype of the elements
 * @param file Filename of the source code file that calls the function
 * @param function Function name of the calling function
 * @param line Line number in the source code file where the call happens
 */
inline void MyMPI_Unpack_real(void *inbuf, int insize, int *position,
                              void *outbuf, int outcount, MPI_Datatype datatype,
                              std::string file, std::string function, int line){
    int status = MPI_Unpack(inbuf, insize, position, outbuf, outcount, datatype,
                            MPI_COMM_WORLD);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Unpack!" << std::endl;
        std::cerr << insize << "\t" << *position << "\t" << std::endl;
        std::cerr << MPIGlobal::rank << ": " << file << "::" << function << ":"
                  << line << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

#define MyMPI_Unpack(ib, is, po, ob, oc, dt){\
    MyMPI_Unpack_real(ib, is, po, ob, oc, dt, __FILE__, __FUNCTION__,\
    __LINE__);\
    }

/**
 * @brief Wrapper around MPI_Allgather
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param sendbuf Buffer to send
 * @param sendcount Number of elements in the send buffer
 * @param sendtype MPI type of elements being sent
 * @param recvbuf Buffer to receive in
 * @param recvcount Number of elements to receive
 * @param recvtype MPI type of elements being received
 * @param comm Communicator identifying the participating processes (default:
 * MPI_COMM_WORLD)
 */
inline void MyMPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype,
                            MPI_Comm comm = MPI_COMM_WORLD){
    MPIGlobal::commtimer.start();
    int status = MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                               recvtype, comm);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Allgather!" << std::endl;
        my_exit();
    }
    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Reduce
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param sendbuf Buffer to send
 * @param recvbuf Buffer to receive in
 * @param count Number of elements to communicate
 * @param datatype MPI datatype of the elements
 * @param op Global reduce operation
 * @param root Rank of the process where the result will be stored
 */
inline void MyMPI_Reduce(void *sendbuf, void *recvbuf, int count,
                         MPI_Datatype datatype, MPI_Op op, int root){
    MPIGlobal::commtimer.start();
    int status = MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root,
                            MPI_COMM_WORLD);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Reduce!" << std::endl;
        my_exit();
    }
    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Irecv
 *
 * We check the error code to detect MPI errors and use the default communicator
 * MPI_WORLD.
 *
 * @param buf Buffer to receive in
 * @param count Number of elements to receive
 * @param datatype MPI datatype of the elements
 * @param source Rank of the source process
 * @param tag Tag identifying this message
 * @param request MPI request to store information in
 */
inline void MyMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
                        int tag, MPI_Request *request){
    MPIGlobal::commtimer.start();
    int status = MPI_Irecv(buf, count, datatype, source, tag, MPI_COMM_WORLD,
                           request);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Irecv!" << std::endl;
        my_exit();
    }
    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Waitall
 *
 * We check the error code to detect MPI errors.
 *
 * @param count Number of requests to wait for
 * @param array_of_requests Array of requests
 * @param array_of_statuses Array to store status information about the requests
 */
inline void MyMPI_Waitall(int count, MPI_Request *array_of_requests,
                        MPI_Status *array_of_statuses){
    MPIGlobal::idletimer.start();
    int status = MPI_Waitall(count, array_of_requests, array_of_statuses);
    if(status != MPI_SUCCESS){
        std::cerr << "Error during MPI_Waitall!" << std::endl;
        my_exit();
    }
    MPIGlobal::idletimer.stop();
}

/**
 * @brief Wrapper around MPI_Test
 *
 * We check the error code to detect MPI errors.
 *
 * @param request Request to test
 * @param flag Flag to store the result of the test in
 * @param status MPI_Status object to store status information about the request
 */
inline void MyMPI_Test(MPI_Request *request, int *flag, MPI_Status *status){
    MPIGlobal::commtimer.start();
    int result = MPI_Test(request, flag, status);
    if(result != MPI_SUCCESS){
        std::cerr << "Error during MPI_Test!" << std::endl;
        my_exit();
    }
    MPIGlobal::commtimer.stop();
}

/**
 * @brief Wrapper around MPI_Init
 *
 * We check the error code to detect MPI errors.
 *
 * This method also initializes the rank and size global variables.
 *
 * @param argc Pointer to the number of command line arguments
 * @param argv Pointer to the command line arguments
 */
inline void MyMPI_Init(int *argc, char ***argv){
    int status = MPI_Init(argc, argv);
    if(status != MPI_SUCCESS){
        std::cerr << "Failed to initialize MPI environment!" << std::endl;
        my_exit();
    }

    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    MyMPI_Comm_rank(&MPIGlobal::rank);
    MyMPI_Comm_size(&MPIGlobal::size);
    MyMPI_Comm_local_vars(&MPIGlobal::local_rank, &MPIGlobal::local_size,
                          &MPIGlobal::noderank, &MPIGlobal::nodesize,
                          &MPIGlobal::nodecomm);

    // allocate the communication buffers
    // we initialize them with 0 size, applications needing it should increase
    // the size
    // by initializing it here, we are sure that the buffer will exist, so that
    // we can safely deallocate it again in MyMPI_Finalize
    MPIGlobal::sendsize = 0;
    MPIGlobal::sendbuffer = new char[MPIGlobal::sendsize];
    MPIGlobal::recvsize = 0;
    MPIGlobal::recvbuffer = new char[MPIGlobal::recvsize];
}

/**
 * @brief Wrapper around MPI_Finalize
 *
 * This method prints the total time spent in MPI processes and during idle time
 * to the stdout. It also frees the node communicator.
 *
 * @return The error code of MPI_Finalize
 */
inline int MyMPI_Finalize(){
    int status = MPI_Comm_free(&MPIGlobal::nodecomm);
    if(status != MPI_SUCCESS){
        std::cerr << "Error while freeing node communicator!" << std::endl;
        my_exit();
    }

    // deallocate communication buffers
    delete [] MPIGlobal::sendbuffer;
    delete [] MPIGlobal::recvbuffer;

    double locidletime, loccommtime;
    double totidletime, totcommtime;
    locidletime = MPIGlobal::idletimer.value();
    loccommtime = MPIGlobal::commtimer.value();
    // funny detail: calling MyMPI_Reduce below will increase the commtimer
    // we however do not take this (small) extra time in account, because we
    // already stored the old value and that one is communicated
    MyMPI_Reduce(&locidletime, &totidletime, 1, MPI_DOUBLE, MPI_SUM, 0);
    MyMPI_Reduce(&loccommtime, &totcommtime, 1, MPI_DOUBLE, MPI_SUM, 0);
    // we take the average over all MPI processes
    totidletime /= MPIGlobal::size;
    totcommtime /= MPIGlobal::size;
    std::cout << "Spent " << totidletime << "s waiting for other MPI processes"
              << std::endl;
    std::cout << "Spent " << totcommtime << "s in MPI communications"
              << std::endl;
    return MPI_Finalize();
}

#endif // MPIMETHODS_HPP
