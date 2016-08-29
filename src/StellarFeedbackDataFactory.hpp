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
 * @file StellarFeedbackDataFactory.hpp
 *
 * @brief Factory class to pack, unpack, dump and reload StellarFeedbackData
 * instances
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef STELLARFEEDBACKDATAFACTORY_HPP
#define STELLARFEEDBACKDATAFACTORY_HPP

#include "DiscreteStellarFeedbackData.hpp"
#include "MPIMethods.hpp"
#include "RestartFile.hpp"
#include "StellarFeedbackData.hpp"

/**
 * @brief Factory class to pack, unpack, dump and reload StellarFeedbackData
 * instances
 */
class StellarFeedbackDataFactory {
  public:
    /**
     * @brief Dump the given StellarFeedbackData instance to the given MPI
     * buffer for communication
     *
     * We first check if the StellarFeedbackData object exists, since some types
     * of feedback do not require additional variables.
     *
     * @param data StellarFeedbackData object to pack
     * @param buffer Buffer to write to
     * @param bufsize Total size of the buffer
     * @param position Current position in the buffer (is updated)
     */
    static void pack(StellarFeedbackData* data, void* buffer, int bufsize,
                     int* position) {
        int has_data;
        if(data) {
            has_data = true;
            MyMPI_Pack(&has_data, 1, MPI_INT, buffer, bufsize, position);
            DiscreteStellarFeedbackData* discrete_data =
                    (DiscreteStellarFeedbackData*)data;
            discrete_data->pack(buffer, bufsize, position);
        } else {
            has_data = false;
            MyMPI_Pack(&has_data, 1, MPI_INT, buffer, bufsize, position);
        }
    }

    /**
     * @brief Create a StellarFeedbackData instance from the given MPI buffer
     *
     * We first check if the buffer contains data, since not all types of
     * feedback require additional variables.
     *
     * @param buffer Buffer to read from
     * @param bufsize Total size of the buffer
     * @param position Current position in the buffer (is updated)
     * @return A new StellarFeedbackData instance or NULL if no data was found
     */
    static StellarFeedbackData* unpack(void* buffer, int bufsize,
                                       int* position) {
        int has_data;
        MyMPI_Unpack(buffer, bufsize, position, &has_data, 1, MPI_INT);
        if(has_data) {
            return new DiscreteStellarFeedbackData(buffer, bufsize, position);
        } else {
            return NULL;
        }
    }

    /**
     * @brief Dump the given StellarFeedbackData instance to the given
     * RestartFile
     *
     * We first check if the StellarFeedbackData object exists, since not all
     * types of feedback require additional variables.
     *
     * @param data StellarFeedbackData instance to dump
     * @param rfile RestartFile to write to
     */
    static void dump(StellarFeedbackData* data, RestartFile& rfile) {
        bool has_data;
        if(data) {
            has_data = true;
            rfile.write(has_data);
            DiscreteStellarFeedbackData* discrete_data =
                    (DiscreteStellarFeedbackData*)data;
            discrete_data->dump(rfile);
        } else {
            has_data = false;
            rfile.write(has_data);
        }
    }

    /**
     * @brief Create a new StellarFeedbackData object from the given RestartFile
     *
     * We first check if the RestartFile contains data, since not all types of
     * stellar feedback require additional variables.
     *
     * @param rfile RestartFile to read from
     * @return A new StellarFeedbackData instance or NULL if no data was found
     */
    static StellarFeedbackData* restart(RestartFile& rfile) {
        bool has_data;
        rfile.read(has_data);
        if(has_data) {
            return new DiscreteStellarFeedbackData(rfile);
        } else {
            return NULL;
        }
    }
};

#endif  // STELLARFEEDBACKDATAFACTORY_HPP
