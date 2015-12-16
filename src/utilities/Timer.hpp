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
 * @file Timer.hpp
 *
 * @brief A simplified interface to the Unix system timer: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TIMER_HPP
#define TIMER_HPP

#include <sys/time.h>
#include <ostream>
#include <istream>

class RestartFile;

/**
  * \brief A simplified interface to the Unix system timer
  *
  * The Timer automatically registers the current system time when constructed
  * and returns the elapsed time in seconds when it is stopped.
  *
  * The Timer can also be used to time multiple intervals, by using the
  * functions Timer::start and Timer::stop. The function Timer::stop always
  * returns the total registered time, which is the sum of all individual
  * intervals measured.
  */
class Timer{
private:
    /*! \brief Starting time of the timer */
    timeval _start;
    /*! \brief Stop time of the timer */
    timeval _stop;
    /*! \brief Total time interval registered so far */
    timeval _diff;

public:
    Timer();
    ~Timer(){}

    void reset();
    void start();
    double stop();
    double value();

    double interval();
    void restart();

    void dump(RestartFile &rfile);
    Timer(RestartFile &rfile);
};

#endif // TIMER_HPP
