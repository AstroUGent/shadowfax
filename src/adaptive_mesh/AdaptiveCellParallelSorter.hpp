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
 * @file AdaptiveCellParallelSorter.hpp
 *
 * @brief A parallel sorter for Hilbert objects: AdaptiveVorCell2d
 * specialization
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ADAPTIVECELLPARALLELSORTER_HPP
#define ADAPTIVECELLPARALLELSORTER_HPP

#include "utilities/ParallelSorter.hpp"

/**
 * @brief Utility class to sort a list of Hilbert_Objects that is spread out
 * over multiple MPI processes, in parallel
 *
 * Also works in the serial case, in which case it only wraps around std::sort.
 */
class AdaptiveCellParallelSorter : public ParallelSorter {
  public:
    AdaptiveCellParallelSorter() {}

    ~AdaptiveCellParallelSorter() {}

    /**
          * @brief Sort the given list in parallel over all MPI processes
          *
          * We use the algorithm of Siebert & Wolf (2010) to do a parallel
     * vector
          * sort in 4 steps:
          *  -# do a local sort
          *  -# find splitters to subdivide all local vectors into pieces that
     * have
          *     to go to different processes
          *  -# redistribute all pieces to their respective processes
          *  -# merge the individually sorted pieces together in a local p-way
          *     merging
          *
          * @param A The list that needs to be sorted
          */
    void sort(std::vector<AdaptiveVorCell2d*>& A) {
        // (1) local sort
        std::sort(A.begin(), A.end(), HB::sortfunc);

        // (2) splitter finding
        unsigned int p = MPIGlobal::size;
        // if the code is serial, we are done now
        if(p < 2) {
            return;
        }
        vector<Splitter> splitter(p + 1);
        unsigned int loc_n = A.size();

        unsigned long loc_min = A[0]->get_key();
        unsigned long glo_min;
        MyMPI_Allreduce(&loc_min, &glo_min, 1, MPI_UNSIGNED_LONG, MPI_MIN);
        splitter[0].set(glo_min, 0);

        unsigned long loc_max = A[loc_n - 1]->get_key();
        unsigned long glo_max;
        MyMPI_Allreduce(&loc_max, &glo_max, 1, MPI_UNSIGNED_LONG, MPI_MAX);
        unsigned int loc_cnt = loc_n - num_less_than(A, glo_max);
        unsigned int glo_cnt;
        MyMPI_Allreduce(&loc_cnt, &glo_cnt, 1, MPI_UNSIGNED, MPI_SUM);
        splitter[p].set(glo_max, glo_cnt);

        unsigned int target_cnt = loc_n;
        MyMPI_Bcast(&target_cnt, 1, MPI_UNSIGNED, 0);
        for(unsigned int i = 1; i < p; i++) {
            splitter[i] =
                    find_splitter(target_cnt, A, splitter[0], splitter[p]);
            unsigned int next_target_cnt = loc_n;
            MyMPI_Bcast(&next_target_cnt, 1, MPI_UNSIGNED, i);
            target_cnt += next_target_cnt;
        }

        // (3) data redistribution
        unsigned int r = MPIGlobal::rank;
        // find number of occurences of the splitkeys for every process and do
        // an exclusive prefix sum
        vector<unsigned int> numindex(p + 1, 0);
        for(unsigned int i = p - 1; i--;) {
            unsigned int noo = num_of_occurences(A, splitter[i + 1].key());
            unsigned int send_noo = noo;
            // we do not use MPI_Exscan since then the value on rank 0 is
            // undefined (and we want it to be 0)
            MyMPI_Scan(&send_noo, &numindex[i + 1], 1, MPI_UNSIGNED, MPI_SUM);
            numindex[i + 1] -= noo;
        }

        vector<AdaptiveVorCell2d*> newlist;
        newlist.reserve(A.size());
        for(unsigned int i = p; i--;) {
            vector<AdaptiveVorCell2d*> sendlist;
            unsigned int j = 0;
            while(j < A.size() &&
                  (A[j] == NULL ||
                   A[j]->get_key() < splitter[(r + i) % p].key())) {
                j++;
            }
            unsigned int k = numindex[(r + i) % p];
            while(j < A.size() &&
                  (A[j] == NULL ||
                   (A[j]->get_key() == splitter[(r + i) % p].key() &&
                    k < splitter[(r + i) % p].index()))) {
                j++;
                k++;
            }
            while(j < A.size() && A[j] != NULL &&
                  A[j]->get_key() < splitter[(r + i) % p + 1].key()) {
                sendlist.push_back(A[j]);
                A[j] = NULL;
                j++;
            }
            k = numindex[(r + i) % p + 1];
            j -= k;
            while(j + k < A.size() && A[j + k] != NULL &&
                  A[j + k]->get_key() == splitter[(r + i) % p + 1].key() &&
                  k < splitter[(r + i) % p + 1].index()) {
                sendlist.push_back(A[j]);
                A[j + k] = NULL;
                k++;
            }
            if(i) {
                // we first inform the other process on the number of particles
                // we will be sending
                unsigned int sendsize = sendlist.size();
                unsigned int send_sendsize = sendsize;
                MPI_Request reqs[2];
                unsigned int recvsize;
                MyMPI_Isend(&send_sendsize, 1, MPI_INT, (r + i) % p, 0,
                            &reqs[0]);
                MyMPI_Irecv(&recvsize, 1, MPI_INT, (r + p - i) % p, 0,
                            &reqs[1]);
                MPI_Status status[2];
                MyMPI_Waitall(2, reqs, status);

                unsigned int numbuf =
                        MPIGlobal::sendsize /
                        (sizeof(ParticleType) + sizeof(GasParticle));
                unsigned int numsend =
                        sendsize / numbuf + (sendsize % numbuf > 0);
                unsigned int numrecv =
                        recvsize / numbuf + (recvsize % numbuf > 0);
                for(unsigned int bi = 0; bi < std::max(numsend, numrecv);
                    bi++) {
                    int sendbufpos[2] = {0, 0};
                    for(unsigned int si = bi * numbuf;
                        si < std::min((bi + 1) * numbuf, sendsize); si++) {
                        sendlist[si]->pack_data(MPIGlobal::sendbuffer,
                                                MPIGlobal::sendsize,
                                                &sendbufpos[0]);
                        sendbufpos[1]++;
                        delete sendlist[si];
                    }
                    // we have to initialize this in case nothing is received
                    // (otherwise the final loop will fail or --even worse --
                    // add particles for a second time)
                    int recv_size[2] = {0, 0};
                    int recv_pos = 0;
                    if(bi < numsend) {
                        // if nothing is sent, then reqs[0] is still ready from
                        // the previous use. So the Waitall call later on will
                        // succeed.
                        MyMPI_Isend(&sendbufpos, 2, MPI_INT, (r + i) % p, 0,
                                    &reqs[0]);
                    }
                    // the +p is necessary to prevent the term between () to
                    // become negative (negative unsigned ints become giant
                    // unsigned ints)
                    if(bi < numrecv) {
                        MyMPI_Irecv(&recv_size, 2, MPI_INT, (r + p - i) % p, 0,
                                    &reqs[1]);
                    }
                    MyMPI_Waitall(2, reqs, status);
                    if(bi < numsend) {
                        MyMPI_Isend(MPIGlobal::sendbuffer, sendbufpos[0],
                                    MPI_PACKED, (r + i) % p, 1, &reqs[0]);
                    }
                    if(bi < numrecv) {
                        MyMPI_Irecv(MPIGlobal::recvbuffer, recv_size[0],
                                    MPI_PACKED, (r + p - i) % p, 1, &reqs[1]);
                    }
                    MyMPI_Waitall(2, reqs, status);
                    for(int j = 0; j < recv_size[1]; j++) {
                        newlist.push_back(
                                new AdaptiveVorCell2d(MPIGlobal::recvbuffer,
                                                      recv_size[0], &recv_pos));
                    }
                }
            } else {
                newlist.insert(newlist.end(), sendlist.begin(), sendlist.end());
            }
        }
        A = newlist;

        // (4) local p-way merge
        // for the moment, this is a common sort, because c++ doesn't have a
        // standard p-merge algorithm
        // (inplace_merge only merges two vectors, which in this case will be
        // too slow)
        std::sort(A.begin(), A.end(), HB::sortfunc);

        return;
    }
};

#endif  // ADAPTIVECELLPARALLELSORTER_HPP
