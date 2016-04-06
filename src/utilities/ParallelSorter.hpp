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
 * @file ParallelSorter.hpp
 *
 * @brief A parallel sorter for Hilbert objects
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 * @author Sven De Rijcke (sven.derijcke@ugent.be)
 */
#ifndef PARALLELSORTER_HPP
#define PARALLELSORTER_HPP

#include "Hilbert.hpp"
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include "ParticleFactory.hpp"
#include <algorithm>
#include <iostream>
#include <mpi.h>
#include <vector>

/**
 * @brief Index-key pair used for the parallel sorting algorithm
 */
class Splitter {
  private:
    /*! \brief Hilbert key at which the split is performed */
    unsigned long _key;

    /*! \brief Index of the key if multiple objects with the same key exist */
    unsigned int _index;

  public:
    Splitter() : _key(0), _index(0) {}

    /**
     * @brief Constructor
     *
     * @param key Hilbert key at which the split is performed
     * @param index Index of the key if multiple objects have the same key
     */
    Splitter(unsigned long key, unsigned int index)
            : _key(key), _index(index) {}

    /**
     * @brief Set the values of the splitter
     *
     * @param key Hilbert key at which the split is performed
     * @param index Index of the key if multiple objects have the same key
     */
    inline void set(unsigned long key, unsigned int index) {
        _key = key;
        _index = index;
    }

    /**
     * @brief Get the key
     *
     * @return Hilbert key at which the split is performed
     */
    inline unsigned long key() const { return _key; }

    /**
     * @brief Get the index
     *
     * @return Index of the key if multiple objects have the same key
     */
    inline unsigned int index() const { return _index; }
};

/**
 * @brief Utility class to sort a list of Hilbert_Objects that is spread out
 * over multiple MPI processes, in parallel
 *
 * Also works in the serial case, in which case it only wraps around std::sort.
 */
class ParallelSorter {
  private:
    /**
     * @brief Count the number of objects in the local list with Hilbert key
     * lower than the given value
     *
     * @param A The list that is being sorted
     * @param max Hilbert key used as reference for the search
     * @return Number of objects in the list with a key lower than the given key
     */
    template <typename T>
    unsigned int num_less_than(std::vector<T*>& A, unsigned long max) {
        Hilbert_Object maxobj(max);
        return distance(A.begin(),
                        lower_bound(A.begin(), A.end(), &maxobj, HB::sortfunc));
    }

    /**
     * @brief Count the number of objects in the local list with Hilbert key
     * higher than the given value
     *
     * @param A The list that is being sorted
     * @param min Hilbert key used as reference for the search
     * @return Number of objects in the list with a key higher than the given
     * key
     */
    template <typename T>
    unsigned int num_greater_than(std::vector<T*>& A, unsigned long min) {
        Hilbert_Object minobj(min);
        return distance(upper_bound(A.begin(), A.end(), &minobj, HB::sortfunc),
                        A.end());
    }

    /**
     * @brief Count the number of occurences of the given Hilbert key in the
     * list
     *
     * @param A The list that is being sorted
     * @param key Hilbert key we look for
     * @return Number of occurences of the given key in the list
     */
    template <typename T>
    unsigned int num_of_occurences(std::vector<T*>& A, unsigned long key) {
        return A.size() - num_less_than(A, key) - num_greater_than(A, key);
    }

    /**
     * @brief Get the Hilbert key of the median element in the local list
     *
     * @param A The list that is being sorted
     * @param minindex Lower index of the part of the list that is considered
     * @param maxindex Upper index of the part of the list that is considered
     * @return Key of the median element
     */
    template <typename T>
    unsigned long median(std::vector<T*>& A, unsigned int minindex,
                         unsigned int maxindex) {
        return A[(maxindex - minindex - 1) / 2 + minindex]->get_key();
    }

    /**
     * @brief Get the median Hilbert key over all processes
     *
     * @param loc_median Local median key
     * @param have_value Flag indicating if we have a local median
     * @return Global median
     */
    unsigned long sparse_median(unsigned long loc_median, bool have_value) {
        unsigned int r = MPIGlobal::rank;
        unsigned int p = MPIGlobal::size;
        vector<unsigned long> medians;
        for(unsigned int i = p; i--;) {
            if(!i) {
                if(have_value) {
                    medians.push_back(loc_median);
                }
            } else {
                MPI_Request reqs[4];
                char send_have_value = have_value;
                MyMPI_Isend(&send_have_value, 1, MPI_CHAR, (r + i) % p, 0,
                            &reqs[0]);
                unsigned long send_loc_median = loc_median;
                MyMPI_Isend(&send_loc_median, 1, MPI_UNSIGNED_LONG, (r + i) % p,
                            1, &reqs[1]);
                char msg_bool;
                MyMPI_Irecv(&msg_bool, 1, MPI_CHAR, (r + p - i) % p, 0,
                            &reqs[2]);
                unsigned long msg_int;
                MyMPI_Irecv(&msg_int, 1, MPI_UNSIGNED_LONG, (r + p - i) % p, 1,
                            &reqs[3]);
                MPI_Status statuses[4];
                MyMPI_Waitall(4, reqs, statuses);
                if(msg_bool) {
                    medians.push_back(msg_int);
                }
            }
        }
        if(!medians.size()) {
            return 0;
        }
        // for large world sizes, it is probably faster to use the order N
        // algorithm std::nth_element
        std::sort(medians.begin(), medians.end());
        return medians[(medians.size() - 1) / 2];
    }

    /**
     * @brief Find the splitters on all processes so that the requested number
     * of objects on all processes lies below this splitter
     *
     * @param target_cnt Total number of objects in the list that should be
     * below the requested splitter on all processes
     * @param A The list that is being sorted
     * @param splitter0 Lower bound for the interval in the list we consider
     * @param splitterp Upper bound for the interval in the list we consider
     * @return A Splitter for every process so that exactly target_cnt elements
     * in A lie below the splitters on all processes combined
     */
    template <typename T>
    Splitter find_splitter(unsigned int target_cnt, std::vector<T*>& A,
                           Splitter& splitter0, Splitter& splitterp) {
        unsigned long min = splitter0.key();
        unsigned long max = splitterp.key();
        unsigned int n = A.size();
        while(true) {
            unsigned int x = n - num_greater_than(A, min);
            unsigned int y = num_less_than(A, max);
            bool have_value = x < y;
            unsigned long loc_median = 0;
            if(have_value) {
                loc_median = median(A, x, y);
            }
            unsigned long glo_median = sparse_median(loc_median, have_value);

            unsigned int loc_cnt = num_less_than(A, glo_median);
            unsigned int a;
            MyMPI_Allreduce(&loc_cnt, &a, 1, MPI_UNSIGNED, MPI_SUM);
            if(target_cnt < a) {
                max = glo_median;
            } else {
                loc_cnt = n - num_greater_than(A, glo_median);
                unsigned int b;
                MyMPI_Allreduce(&loc_cnt, &b, 1, MPI_UNSIGNED, MPI_SUM);
                if(target_cnt > b) {
                    min = glo_median;
                } else {
                    return Splitter(glo_median, target_cnt - a);
                }
            }
        }
    }

    /**
     * @brief Find the splitters on all processes so that the requested weight
     * of objects on all processes lies below this splitter
     *
     * @param target_weight Total weight of objects in the list that should be
     * below the requested splitter on all processes
     * @param A The list that is being sorted
     * @param cumulative_weights Cumulative weights of the elements in the list
     * @param splitter0 Lower bound for the interval in the list we consider
     * @param splitterp Upper bound for the interval in the list we consider
     * @return A Splitter for every process so that the combined weight of the
     * objects on all processes is approximately target_weight
     */
    template <typename T>
    Splitter find_splitter(unsigned int target_weight, std::vector<T*>& A,
                           std::vector<unsigned int>& cumulative_weights,
                           Splitter& splitter0, Splitter& splitterp) {
        unsigned long min = splitter0.key();
        unsigned long max = splitterp.key();
        unsigned int minweight = 0;
        unsigned int maxweight = 0;
        unsigned int n = A.size();
        while(true) {
            unsigned int x = n - num_greater_than(A, min);
            unsigned int y = num_less_than(A, max);
            bool have_value = x < y;
            unsigned long loc_median = 0;
            if(have_value) {
                loc_median = median(A, x, y);
            }
            unsigned long glo_median = sparse_median(loc_median, have_value);
            if(!glo_median) {
                // normal exit point: the interval can no longer be split. We
                // have to choose a split value
                if(target_weight - minweight < maxweight - target_weight) {
                    return Splitter(min, 0);
                } else {
                    return Splitter(max, 0);
                }
            }

            unsigned int loc_cnt = num_less_than(A, glo_median);
            unsigned int loc_weight;
            if(loc_cnt) {
                loc_weight = cumulative_weights[loc_cnt - 1];
            } else {
                loc_weight = 0;
            }
            unsigned int a;
            MyMPI_Allreduce(&loc_weight, &a, 1, MPI_UNSIGNED, MPI_SUM);
            if(target_weight < a) {
                max = glo_median;
                maxweight = a;
            } else {
                if(target_weight > a) {
                    min = glo_median;
                    minweight = a;
                } else {
                    // we do not handle the case where two keys are equal!!!!
                    // in case a == target_weight (very unlikely)
                    return Splitter(glo_median, 0);
                }
            }
        }
    }

  public:
    ParallelSorter() {}

    ~ParallelSorter() {}

    /**
      * @brief Sort the given list in parallel over all MPI processes
      *
      * We use the algorithm of Siebert & Wolf (2010) to do a parallel vector
      * sort in 4 steps:
      *  -# do a local sort
      *  -# find splitters to subdivide all local vectors into pieces that have
      *     to go to different processes
      *  -# redistribute all pieces to their respective processes
      *  -# merge the individually sorted pieces together in a local p-way
      *     merging
      *
      * @param A The list that needs to be sorted
      */
    template <typename T> void sort(std::vector<T*>& A) {
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

        vector<T*> newlist;
        newlist.reserve(A.size());
        for(unsigned int i = p; i--;) {
            vector<T*> sendlist;
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
                        ParticleFactory::dump(
                                sendlist[si], MPIGlobal::sendbuffer,
                                MPIGlobal::sendsize, &sendbufpos[0]);
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
                        //                        newlist.push_back(new
                        //                        T(MPIGlobal::recvbuffer,
                        //                                                recv_size[0],
                        //                                                &recv_pos));
                        newlist.push_back(
                                ParticleFactory::load(MPIGlobal::recvbuffer,
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

    /**
      * @brief Sort the given list in parallel over all MPI processes, using
      * the computational cost as a weight
      *
      * We use the algorithm of Siebert & Wolf (2010) to do a parallel vector
      * sort in 4 steps:
      *  -# do a local sort
      *  -# find splitters to subdivide all local vectors into pieces that have
      *     to go to different processes
      *  -# redistribute all pieces to their respective processes
      *  -# merge the individually sorted pieces together in a local p-way
      *     merging
      *
      * @param A The list that needs to be sorted
      */
    template <typename T> void weightsort(std::vector<T*>& A) {
        // (1) local sort
        std::sort(A.begin(), A.end(), HB::sortfunc);

        // (2) splitter finding
        unsigned int p = MPIGlobal::size;
        // if the code is serial, we are done now
        if(p < 2) {
            return;
        }

        std::vector<unsigned int> cumulative_weights(A.size());
        unsigned int totweight = 0;
        cumulative_weights[0] = A[0]->get_comp_cost();
        totweight += A[0]->get_comp_cost();
        for(unsigned int i = 1; i < A.size(); i++) {
            cumulative_weights[i] =
                    cumulative_weights[i - 1] + A[i]->get_comp_cost();
            totweight += A[i]->get_comp_cost();
        }

        vector<Splitter> splitter(p + 1);
        unsigned int loc_n = A.size();
        unsigned int globtotweight;
        MyMPI_Allreduce(&totweight, &globtotweight, 1, MPI_UNSIGNED, MPI_SUM);
        if(!globtotweight) {
            sort(A);
            return;
        }
        unsigned int average_weight = globtotweight / p;

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

        unsigned int target_weight = average_weight;
        MyMPI_Bcast(&target_weight, 1, MPI_UNSIGNED, 0);
        for(unsigned int i = 1; i < p; i++) {
            splitter[i] = find_splitter(target_weight, A, cumulative_weights,
                                        splitter[0], splitter[p]);
            unsigned int next_target_weight = average_weight;
            MyMPI_Bcast(&next_target_weight, 1, MPI_UNSIGNED, i);
            target_weight += next_target_weight;
        }

        // (3) data redistribution
        unsigned int r = MPIGlobal::rank;
        // find number of occurences of the splitkeys for every process and do
        // an exclusive prefix sum
        vector<unsigned int> numindex(p + 1, 0);
        for(unsigned int i = p - 1; i--;) {
            unsigned int noo = num_of_occurences(A, splitter[i + 1].key());
            MyMPI_Scan(&noo, &numindex[i + 1], 1, MPI_UNSIGNED, MPI_SUM);
            numindex[i + 1] -= noo;
        }

        vector<T*> newlist;
        newlist.reserve(A.size());
        for(unsigned int i = p; i--;) {
            vector<T*> sendlist;
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
                MPI_Request reqs[2];
                unsigned int recvsize;
                MyMPI_Isend(&sendsize, 1, MPI_INT, (r + i) % p, 0, &reqs[0]);
                MyMPI_Irecv(&recvsize, 1, MPI_INT, (r + p - i) % p, 0,
                            &reqs[1]);
                MPI_Status status[2];
                MyMPI_Waitall(2, reqs, status);

                unsigned int numbuf = MPIGlobal::sendsize / sizeof(T);
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
                        // if nothing is send, then reqs[0] is still ready from
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
                        newlist.push_back(new T(MPIGlobal::recvbuffer,
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

#endif  // PARALLELSORTER_HPP
