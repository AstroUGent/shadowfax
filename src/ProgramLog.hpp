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
 * @file ProgramLog.hpp
 *
 * @brief Global program logging mechanism: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PROGRAMLOG_HPP
#define PROGRAMLOG_HPP

// log level for the program
//#define LOGLEVEL_ALL

#ifdef LOGLEVEL_ALL

// all logging is activated
#define LOGINIT(foldername) global_logger.open(foldername)
#define LOGS(message) global_logger.log<LOGLEVEL_S>(message)
#define LOGW(message) global_logger.log<LOGLEVEL_W>(message)
#define LOGE(message) global_logger.log<LOGLEVEL_E>(message)
#define DOLOG

#else  // LOGLEVEL_ALL

#ifdef LOGLEVEL_WARNINGS

// only errors and warnings are logged
#define LOGINIT(foldername) global_logger.open(foldername)
#define LOGS(message)
#define LOGW(message) global_logger.log<LOGLEVEL_W>(message)
#define LOGE(message) global_logger.log<LOGLEVEL_E>(message)
#define DOLOG

#else  // LOGLEVEL_WARNINGS

#ifdef LOGLEVEL_ERRORS

// only errors are logged
#define LOGINIT(foldername) global_logger.open(foldername)
#define LOGS(message)
#define LOGW(message)
#define LOGE(message) global_logger.log<LOGLEVEL_E>(message)
#define DOLOG

#else  // LOGLEVEL_ERRORS

// no logging
#define LOGINIT(foldername)
#define LOGS(message)
#define LOGW(message)
#define LOGE(message)

#endif  // LOGLEVEL_ERRORS

#endif  // LOGLEVEL_WARNINGS

#endif  // LOGLEVEL_ALL

#ifdef DOLOG

#include "MPIGlobal.hpp"
#include <fstream>
#include <string>

/**
 * @brief Three possible levels of logging
 */
enum LogLevel {
    /*! @brief Lowest level of logging, provides most information */
    LOGLEVEL_S = 0,
    /*! @brief Warning level. Invoked to log non-fatal warnings */
    LOGLEVEL_W,
    /*! @brief Error level. Invoked to log fatal program errors */
    LOGLEVEL_E
};

/**
 * @brief Class used for global program logging
 *
 * All logs are written to a file programlog.txt in the program output
 * directory and contain a time label, a log level label and a custom message.
 */
class ProgramLog {
  private:
    /*! @brief Output stream to write logs to */
    std::ofstream _ofile;

  public:
    /**
     * @brief Open the log file in the given folder
     *
     * @param foldername Name of the folder where the logfile should be created
     */
    void open(std::string foldername) {
        if(MPIGlobal::rank) {
            return;
        }
        std::string filename = foldername + std::string("/programlog.txt");
        _ofile.open(filename.c_str());
    }

    /**
     * @brief Log a message with the given template log level label
     *
     * @param message Custom message to append to the log line
     */
    template <LogLevel level> void log(std::string message) {
        if(MPIGlobal::rank) {
            return;
        }
        timeval curtime;
        struct tm* date;
        gettimeofday(&curtime, NULL);
        date = localtime(&curtime.tv_sec);
        // convert the date to a string format
        // 11: 2 (hour) + 2 (minutes) + 2 (seconds) + 1 space + 3 colons + '\0'
        char buffer[11];
        strftime(buffer, 11, "%H:%M:%S: ", date);
        _ofile << buffer;

        switch(level) {
            case LOGLEVEL_S:
                _ofile << "S: ";
                break;
            case LOGLEVEL_W:
                _ofile << "W: ";
                break;
            case LOGLEVEL_E:
                _ofile << "E: ";
                break;
        }

        _ofile << message << std::endl;
    }
};

/*! @brief Global ProgramLog instance */
extern ProgramLog global_logger;

#endif

#endif  // PROGRAMLOG_HPP
