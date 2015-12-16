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
 * @file LogFiles.hpp
 *
 * @brief User-friendly diagnostic log files: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef LOGFILES_HPP
#define LOGFILES_HPP

#include <fstream>
#include <string>

#include "utilities/Timer.hpp"

class TimeLine;
class RestartFile;
class ParticleVector;
class RiemannSolver;

/**
 * @brief User-friendly diagnostic log files
 *
 * Information written every timestep to inform the user about what is going on.
 */
class LogFiles{
private:
    // time log
    /*! \brief Path of the time log file */
    std::string _timename;
    /*! \brief Stream to write time log information to */
    std::ofstream _timefile;

    // energy log
    /*! \brief Path of the energy log file */
    std::string _energyname;
    /*! \brief Stream to write energy log information to */
    std::ofstream _energyfile;

    // statistics
    /*! \brief Path of the statistics log file */
    std::string _statname;
    /*! \brief Stream to write statistics log information to */
    std::ofstream _statfile;
    /*! \brief Step counter */
    unsigned int _stat_stepcount;
    /*! \brief Previous number of Riemann solver evaluations */
    unsigned long _stat_solvercount;

    /*! \brief Timer to quantify time spent writing log files */
    Timer _timer;

    void copy_file_single(std::ifstream &ifile, std::ofstream &ofile, long pos);
    void copy_file_buffer(std::ifstream &ifile, std::ofstream &ofile, long pos);
    void copy_file(const std::string &iname, const std::string &oname,
                   long pos);
    void reset_file(std::string name, long pos);

public:
    LogFiles(std::string outputdir);
    ~LogFiles();

    void write_to_logs(TimeLine* &timeline, ParticleVector* &particles,
                       RiemannSolver* &solver);

    LogFiles(RestartFile &rfile);
    void dump(RestartFile &rfile);
};

#endif // LOGFILES_HPP
