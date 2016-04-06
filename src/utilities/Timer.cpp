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
 * @file Timer.cpp
 *
 * @brief A simplified interface to the Unix system timer: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Timer.hpp"
#include "RestartFile.hpp"
#include <cstdlib>
using namespace std;

/**
 * @brief Constructor
 *
 * Intialize the internal timeval difference and register the current system
 * time.
 */
Timer::Timer() {
    reset();
    gettimeofday(&_start, NULL);
}

/**
 * @brief Clear the internal timeval difference
 */
void Timer::reset() { timerclear(&_diff); }

/**
 * @brief Record the current system time as starting time
 */
void Timer::start() { gettimeofday(&_start, NULL); }

/**
 * @brief Record the current system time as stopping time and add the difference
 * between start and stop to the internal timeval difference
 *
 * @return The current contents of the internal timeval difference in seconds
 * (with microsecond precision)
 */
double Timer::stop() {
    gettimeofday(&_stop, NULL);
    timeval interval_diff;
    timersub(&_stop, &_start, &interval_diff);
    timeradd(&_diff, &interval_diff, &_diff);
    return _diff.tv_sec + 1.e-6 * _diff.tv_usec;
}

/**
 * @brief Get the current internal timeval difference
 *
 * @return The current contents of the internal timeval difference in seconds
 * (with microsecond precision)
 */
double Timer::value() { return _diff.tv_sec + 1.e-6 * _diff.tv_usec; }

/**
 * @brief Get the current value of the timer without affecting it
 *
 * @return The time in seconds since the timer was last started
 */
double Timer::interval() {
    timeval tempstop;
    gettimeofday(&tempstop, NULL);
    timeval interval;
    timersub(&tempstop, &_start, &interval);
    return interval.tv_sec + 1.e-6 * interval.tv_usec;
}

/**
 * @brief Restart the timer by overwriting the start time.
 */
void Timer::restart() { gettimeofday(&_start, NULL); }

/**
 * @brief Dump the timer to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void Timer::dump(RestartFile& rfile) {
    rfile.write(_start);
    rfile.write(_stop);
    rfile.write(_diff);
}

/**
 * @brief Restart constructor. Initialize the timer based on the given
 * RestartFile
 *
 * @param rfile RestartFile to read from
 */
Timer::Timer(RestartFile& rfile) {
    rfile.read(_start);
    rfile.read(_stop);
    rfile.read(_diff);
}
