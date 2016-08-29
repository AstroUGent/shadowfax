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
 * @file DiscreteStellarFeedbackData.hpp
 *
 * @brief StellarFeedbackData implementation for discrete stellar feedback
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef DISCRETESTELLARFEEDBACKDATA_HPP
#define DISCRETESTELLARFEEDBACKDATA_HPP

#include "MPIMethods.hpp"
#include "RestartFile.hpp"
#include "StellarFeedbackData.hpp"

/**
 * @brief StellarFeedbackData implementation for discrete stellar feedback
 */
class DiscreteStellarFeedbackData : public StellarFeedbackData {
  private:
    /*! @brief Number of PopII SNII explosions that needs to go off */
    unsigned int _PopII_SNII_number;
    /*! @brief Number of PopII SNIa explosions that needs to go off */
    unsigned int _PopII_SNIa_number;
    /*! @brief Number of PopIII SN explosions that needs to go off */
    unsigned int _PopIII_SN_number;

    /*! @brief Total number of PopII SW type stars for the StarParticle */
    double _PopII_SW_fac;
    /*! @brief Total number of PopIII SW type stars for the StarParticle */
    double _PopIII_SW_fac;

    /*! @brief Number of PopII SNII explosions that has gone off */
    unsigned int _PopII_SNII_count;
    /*! @brief Number of PopII SNIa explosions that has gone off */
    unsigned int _PopII_SNIa_count;
    /*! @brief Number of PopIII SN explosions that has gone off */
    unsigned int _PopIII_SN_count;

    /*! @brief Next time a PopII SNII needs to go off */
    double _PopII_SNII_next_time;
    /*! @brief Upper limit for the next mass interval in which a PopII SNII
     *  should go off */
    double _PopII_SNII_interval;
    /*! @brief PopII SNII factor */
    double _PopII_SNII_fac;

    /*! @brief Next time a PopII SNIa needs to go off */
    double _PopII_SNIa_next_time;
    /*! @brief Lower limit for the next time interval in which a PopII SNIa
     *  should go off */
    double _PopII_SNIa_interval;
    /*! @brief PopII SNIa factor */
    double _PopII_SNIa_fac;

    /*! @brief Next time a PopIII SN needs to go off */
    double _PopIII_SN_next_time;
    /*! @brief Upper limit for the next mass interval in which a PopIII SN
     *  should go off */
    double _PopIII_SN_interval;
    /*! @brief PopIII SN factor */
    double _PopIII_SN_fac;

  public:
    /**
     * @brief Empty constructor
     */
    inline DiscreteStellarFeedbackData() {
        _PopII_SNII_number = 0;
        _PopII_SNIa_number = 0;
        _PopIII_SN_number = 0;

        _PopII_SW_fac = 0.;
        _PopIII_SW_fac = 0.;

        _PopII_SNII_count = 0;
        _PopII_SNIa_count = 0;
        _PopIII_SN_count = 0;

        _PopII_SNII_next_time = 0.;
        _PopII_SNII_interval = 0.;
        _PopII_SNII_fac = 0.;

        _PopII_SNIa_next_time = 0.;
        _PopII_SNIa_interval = 0.;
        _PopII_SNIa_fac = 0.;

        _PopIII_SN_next_time = 0.;
        _PopIII_SN_interval = 0.;
        _PopIII_SN_fac = 0.;
    }

    /**
     * @brief Set the total number of PopII SNII that needs to go off
     *
     * @param PopII_SNII_number Total number of PopII SNII
     */
    inline void set_PopII_SNII_number(unsigned int PopII_SNII_number) {
        _PopII_SNII_number = PopII_SNII_number;
    }

    /**
     * @brief Get the total number of PopII SNII that needs to go off
     *
     * @return Total number of PopII SNII
     */
    inline unsigned int get_PopII_SNII_number() {
        return _PopII_SNII_number;
    }

    /**
     * @brief Set the total number of PopII SNIa that needs to go off
     *
     * @param PopII_SNIa_number Total number of PopII SNIa
     */
    inline void set_PopII_SNIa_number(unsigned int PopII_SNIa_number) {
        _PopII_SNIa_number = PopII_SNIa_number;
    }

    /**
     * @brief Get the total number of PopII SNIa that needs to go off
     *
     * @return Total number of PopII SNIa
     */
    inline unsigned int get_PopII_SNIa_number() {
        return _PopII_SNIa_number;
    }

    /**
     * @brief Set the total number of PopIII SN that needs to go off
     *
     * @param PopIII_SN_number Total number of PopIII SN
     */
    inline void set_PopIII_SN_number(unsigned int PopIII_SN_number) {
        _PopIII_SN_number = PopIII_SN_number;
    }

    /**
     * @brief Get the total number of PopIII SN that needs to go off
     *
     * @return Total number of PopIII SN
     */
    inline unsigned int get_PopIII_SN_number() {
        return _PopIII_SN_number;
    }

    /**
     * @brief Set the PopII SW factor
     *
     * @param PopII_SW_fac PopII SW factor
     */
    inline void set_PopII_SW_fac(double PopII_SW_fac) {
        _PopII_SW_fac = PopII_SW_fac;
    }

    /**
     * @brief Get the PopII SW factor
     *
     * @return PopII SW factor
     */
    inline double get_PopII_SW_fac() {
        return _PopII_SW_fac;
    }

    /**
     * @brief Set the PopIII SW factor
     *
     * @param PopIII_SW_fac PopIII SW factor
     */
    inline void set_PopIII_SW_fac(double PopIII_SW_fac) {
        _PopIII_SW_fac = PopIII_SW_fac;
    }

    /**
     * @brief Get the PopIII SW factor
     *
     * @return PopIII SW factor
     */
    inline double get_PopIII_SW_fac() {
        return _PopIII_SW_fac;
    }

    /**
     * @brief Set the number of PopII SNII explosions that has gone off
     *
     * @param PopII_SNII_count Number of PopII SNII explosions that has gone off
     */
    inline void set_PopII_SNII_count(unsigned int PopII_SNII_count) {
        _PopII_SNII_count = PopII_SNII_count;
    }

    /**
     * @brief Get the number of PopII SNII explosions that has gone off
     *
     * @return Number of PopII SNII explosions that has gone off
     */
    inline unsigned int get_PopII_SNII_count() {
        return _PopII_SNII_count;
    }

    /**
     * @brief Increase the number of PopII SNII explosions that has gone off
     */
    inline void increase_PopII_SNII_count() {
        _PopII_SNII_count++;
    }

    /**
     * @brief Set the number of PopII SNIa explosions that has gone off
     *
     * @param PopII_SNIa_count Number of PopII SNIa explosions that has gone off
     */
    inline void set_PopII_SNIa_count(unsigned int PopII_SNIa_count) {
        _PopII_SNIa_count = PopII_SNIa_count;
    }

    /**
     * @brief Get the number of PopII SNIa explosions that has gone off
     *
     * @return Number of PopII SNIa explosions that has gone off
     */
    inline unsigned int get_PopII_SNIa_count() {
        return _PopII_SNIa_count;
    }

    /**
     * @brief Increase the number of PopII SNIa explosions that has gone off
     */
    inline void increase_PopII_SNIa_count() {
        _PopII_SNIa_count++;
    }

    /**
     * @brief Set the number of PopIII SN explosions that has gone off
     *
     * @param PopIII_SN_count Number of PopIII SN explosions that has gone off
     */
    inline void set_PopIII_SN_count(unsigned int PopIII_SN_count) {
        _PopIII_SN_count = PopIII_SN_count;
    }

    /**
     * @brief Get the number of PopIII SN explosions that has gone off
     *
     * @return Number of PopIII SN explosions that has gone off
     */
    inline unsigned int get_PopIII_SN_count() {
        return _PopIII_SN_count;
    }

    /**
     * @brief Increase the number of PopIII SN explosions that has gone off
     */
    inline void increase_PopIII_SN_count() {
        _PopIII_SN_count++;
    }

    /**
     * @brief Set the next time a PopII SNII should go off
     *
     * @param PopII_SNII_next_time Next time a PopII SNII should go off
     */
    inline void set_PopII_SNII_next_time(double PopII_SNII_next_time) {
        _PopII_SNII_next_time = PopII_SNII_next_time;
    }

    /**
     * @brief Get the next time a PopII SNII should go off
     *
     * @return Next time a PopII SNII should go off
     */
    inline double get_PopII_SNII_next_time() {
        return _PopII_SNII_next_time;
    }

    /**
     * @brief Set the upper limit for the next mass interval in which a PopII
     * SNII should go off
     *
     * @param PopII_SNII_interval Upper limit for the next mass interval in
     * which a PopII SNII should go off
     */
    inline void set_PopII_SNII_interval(double PopII_SNII_interval) {
        _PopII_SNII_interval = PopII_SNII_interval;
    }

    /**
     * @brief Get the upper limit for the next mass interval in which a PopII
     * SNII should go off
     *
     * @return Upper limit for the next mass interval in which a PopII SNII
     * should go off
     */
    inline double get_PopII_SNII_interval() {
        return _PopII_SNII_interval;
    }

    /**
     * @brief Set the PopII SNII factor
     *
     * @param PopII_SNII_fac PopII SNII factor
     */
    inline void set_PopII_SNII_fac(double PopII_SNII_fac) {
        _PopII_SNII_fac = PopII_SNII_fac;
    }

    /**
     * @brief Get the PopII SNII factor
     *
     * @return PopII SNII factor
     */
    inline double get_PopII_SNII_fac() {
        return _PopII_SNII_fac;
    }

    /**
     * @brief Set the next time a PopII SNIa should go off
     *
     * @param PopII_SNIa_next_time Next time a PopII SNIa should go off
     */
    inline void set_PopII_SNIa_next_time(double PopII_SNIa_next_time) {
        _PopII_SNIa_next_time = PopII_SNIa_next_time;
    }

    /**
     * @brief Get the next time a PopII SNIa should go off
     *
     * @return Next time a PopII SNIa should go off
     */
    inline double get_PopII_SNIa_next_time() {
        return _PopII_SNIa_next_time;
    }

    /**
     * @brief Set the lower limit for the next time interval in which a PopII
     * SNIa should go off
     *
     * @param PopII_SNIa_interval Upper limit for the next mass interval in
     * which a PopII SNIa should go off
     */
    inline void set_PopII_SNIa_interval(double PopII_SNIa_interval) {
        _PopII_SNIa_interval = PopII_SNIa_interval;
    }

    /**
     * @brief Get the lower limit for the next time interval in which a PopII
     * SNIa should go off
     *
     * @return Upper limit for the next mass interval in which a PopII SNIa
     * should go off
     */
    inline double get_PopII_SNIa_interval() {
        return _PopII_SNIa_interval;
    }

    /**
     * @brief Set the PopII SNIa factor
     *
     * @param PopII_SNIa_fac PopII SNIa factor
     */
    inline void set_PopII_SNIa_fac(double PopII_SNIa_fac) {
        _PopII_SNIa_fac = PopII_SNIa_fac;
    }

    /**
     * @brief Get the PopII SNIa factor
     *
     * @return PopII SNIa factor
     */
    inline double get_PopII_SNIa_fac() {
        return _PopII_SNIa_fac;
    }

    /**
     * @brief Set the next time a PopIII SN should go off
     *
     * @param PopIII_SN_next_time Next time a PopIII SN should go off
     */
    inline void set_PopIII_SN_next_time(double PopIII_SN_next_time) {
        _PopIII_SN_next_time = PopIII_SN_next_time;
    }

    /**
     * @brief Get the next time a PopIII SN should go off
     *
     * @return Next time a PopIII SN should go off
     */
    inline double get_PopIII_SN_next_time() {
        return _PopIII_SN_next_time;
    }

    /**
     * @brief Set the upper limit for the next mass interval in which a PopIII
     * SN should go off
     *
     * @param PopIII_SN_interval Upper limit for the next mass interval in
     * which a PopIII SN should go off
     */
    inline void set_PopIII_SN_interval(double PopIII_SN_interval) {
        _PopIII_SN_interval = PopIII_SN_interval;
    }

    /**
     * @brief Get the upper limit for the next mass interval in which a PopIII
     * SN should go off
     *
     * @return Upper limit for the next mass interval in which a PopIII SN
     * should go off
     */
    inline double get_PopIII_SN_interval() {
        return _PopIII_SN_interval;
    }

    /**
     * @brief Set the PopIII SN factor
     *
     * @param PopIII_SN_fac PopIII SN factor
     */
    inline void set_PopIII_SN_fac(double PopIII_SN_fac) {
        _PopIII_SN_fac = PopIII_SN_fac;
    }

    /**
     * @brief Get the PopIII SN factor
     *
     * @return PopIII SN factor
     */
    inline double get_PopIII_SN_fac() {
        return _PopIII_SN_fac;
    }

    /**
     * @brief Pack the variables to the given MPI buffer for communication
     *
     * @param buffer Buffer to write to
     * @param bufsize Total size of the buffer
     * @param position Current position inside the buffer (is updated)
     */
    inline void pack(void* buffer, int bufsize, int* position) {
        MyMPI_Pack(&_PopII_SNII_number, 1, MPI_UNSIGNED, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopII_SNIa_number, 1, MPI_UNSIGNED, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopIII_SN_number, 1, MPI_UNSIGNED, buffer, bufsize,
                   position);

        MyMPI_Pack(&_PopII_SW_fac, 1, MPI_DOUBLE, buffer, bufsize, position);
        MyMPI_Pack(&_PopIII_SW_fac, 1, MPI_DOUBLE, buffer, bufsize, position);

        MyMPI_Pack(&_PopII_SNII_count, 1, MPI_UNSIGNED, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopII_SNIa_count, 1, MPI_UNSIGNED, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopIII_SN_count, 1, MPI_UNSIGNED, buffer, bufsize,
                   position);

        MyMPI_Pack(&_PopII_SNII_next_time, 1, MPI_DOUBLE, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopII_SNII_interval, 1, MPI_DOUBLE, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopII_SNII_fac, 1, MPI_DOUBLE, buffer, bufsize, position);

        MyMPI_Pack(&_PopII_SNIa_next_time, 1, MPI_DOUBLE, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopII_SNIa_interval, 1, MPI_DOUBLE, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopII_SNIa_fac, 1, MPI_DOUBLE, buffer, bufsize, position);

        MyMPI_Pack(&_PopIII_SN_next_time, 1, MPI_DOUBLE, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopIII_SN_interval, 1, MPI_DOUBLE, buffer, bufsize,
                   position);
        MyMPI_Pack(&_PopIII_SN_fac, 1, MPI_DOUBLE, buffer, bufsize, position);
    }

    /**
     * @brief MPI constructor
     *
     * @param buffer Buffer to read from
     * @param bufsize Total size of the buffer
     * @param position Current position in the buffer (is updated)
     */
    inline DiscreteStellarFeedbackData(void* buffer, int bufsize,
                                       int* position) {
        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNII_number, 1,
                     MPI_UNSIGNED);
        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNIa_number, 1,
                     MPI_UNSIGNED);
        MyMPI_Unpack(buffer, bufsize, position, &_PopIII_SN_number, 1,
                     MPI_UNSIGNED);

        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SW_fac, 1, MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_PopIII_SW_fac, 1, MPI_DOUBLE);

        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNII_count, 1,
                     MPI_UNSIGNED);
        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNIa_count, 1,
                     MPI_UNSIGNED);
        MyMPI_Unpack(buffer, bufsize, position, &_PopIII_SN_count, 1,
                     MPI_UNSIGNED);

        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNII_next_time, 1,
                     MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNII_interval, 1,
                     MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNII_fac, 1,
                     MPI_DOUBLE);

        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNIa_next_time, 1,
                     MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNIa_interval, 1,
                     MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_PopII_SNIa_fac, 1,
                     MPI_DOUBLE);

        MyMPI_Unpack(buffer, bufsize, position, &_PopIII_SN_next_time, 1,
                     MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_PopIII_SN_interval, 1,
                     MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_PopIII_SN_fac, 1, MPI_DOUBLE);
    }

    /**
     * @brief Dump the DiscreteStellarFeedbackData instance to the given
     * RestartFile
     *
     * @param rfile RestartFile to write to
     */
    inline void dump(RestartFile& rfile) {
        rfile.write(_PopII_SNII_number);
        rfile.write(_PopII_SNIa_number);
        rfile.write(_PopIII_SN_number);

        rfile.write(_PopII_SW_fac);
        rfile.write(_PopIII_SW_fac);

        rfile.write(_PopII_SNII_count);
        rfile.write(_PopII_SNIa_count);
        rfile.write(_PopIII_SN_count);

        rfile.write(_PopII_SNII_next_time);
        rfile.write(_PopII_SNII_interval);
        rfile.write(_PopII_SNII_fac);

        rfile.write(_PopII_SNIa_next_time);
        rfile.write(_PopII_SNIa_interval);
        rfile.write(_PopII_SNIa_fac);

        rfile.write(_PopIII_SN_next_time);
        rfile.write(_PopIII_SN_interval);
        rfile.write(_PopIII_SN_fac);
    }

    /**
     * @brief Restart constructor
     *
     * @param rfile RestartFile to read from
     */
    inline DiscreteStellarFeedbackData(RestartFile& rfile) {
        rfile.read(_PopII_SNII_number);
        rfile.read(_PopII_SNIa_number);
        rfile.read(_PopIII_SN_number);

        rfile.read(_PopII_SW_fac);
        rfile.read(_PopIII_SW_fac);

        rfile.read(_PopII_SNII_count);
        rfile.read(_PopII_SNIa_count);
        rfile.read(_PopIII_SN_count);

        rfile.read(_PopII_SNII_next_time);
        rfile.read(_PopII_SNII_interval);
        rfile.read(_PopII_SNII_fac);

        rfile.read(_PopII_SNIa_next_time);
        rfile.read(_PopII_SNIa_interval);
        rfile.read(_PopII_SNIa_fac);

        rfile.read(_PopIII_SN_next_time);
        rfile.read(_PopIII_SN_interval);
        rfile.read(_PopIII_SN_fac);
    }
};

#endif  // DISCRETESTELLARFEEDBACKDATA_HPP
