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
 * @file IOModule.cpp
 *
 * @brief Expose some input/output methods to Python
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "../io/AsciiInput.hpp"
#include "../io/AsciiOutput.hpp"
#include "../io/Block.hpp"
#include "../io/FileInput.hpp"
#include "../io/FileOutput.hpp"
#include "../io/Header.hpp"
#include "../io/Unit.hpp"
#include "../io/UnitConverter.hpp"
#include "../io/UnitSet.hpp"
#include <string>
#include <vector>
using namespace std;

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/module.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

/**
 * @brief Custom constructor for Python accessible Block instances
 *
 * @param name Name of the Block
 * @param headers Headers of the Block columns
 * @param dimensions Dimensions of the Block
 * @param units Units of the Block columns
 * @param size Number of elements in the Block
 * @return A Python type pointer to a Block instance
 */
boost::shared_ptr<Block> MakeBlock(std::string name, list& headers,
                                   list& dimensions, list& units,
                                   unsigned int size) {
    vector<string> vecheaders;
    for(int i = 0; i < len(headers); i++) {
        vecheaders.push_back(extract<string>(headers[i]));
    }
    vector<unsigned int> vecdimensions;
    for(int i = 0; i < len(dimensions); i++) {
        vecdimensions.push_back(extract<unsigned int>(dimensions[i]));
    }
    vector<Unit> vecunits;
    for(int i = 0; i < len(units); i++) {
        vecunits.push_back(extract<Unit>(units[i]));
    }
    return boost::shared_ptr<Block>(
            new Block(name, vecheaders, vecdimensions, vecunits, size));
}

#if ndim_ == 3
/*! @brief Register the libpython_io Python module */
BOOST_PYTHON_MODULE(libpython_io3d)
#else
BOOST_PYTHON_MODULE(libpython_io2d)
#endif
{
    class_<vector<double> >("vector<double>")
            .def(vector_indexing_suite<vector<double> >());

    class_<FileInput>("FileInput", init<string>())
            .def("read_header", &FileInput::read_header)
            .def("read", &FileInput::read);

    class_<AsciiInput>("AsciiInput", init<string>())
            .def("read_header", &AsciiInput::read_header)
            .def("read", &AsciiInput::read);

    class_<Header>("Header")
            .def("ngaspart", &Header::ngaspart)
            .def("npart", &Header::npart)
            .def("ndmpart", &Header::ndmpart)
            .def("time", &Header::time)
            .def("global_timestep", &Header::global_timestep);

    class_<Block, boost::noncopyable, boost::shared_ptr<Block> >("Block",
                                                                 no_init)
            .def("__init__", make_constructor(&MakeBlock))
            .def("number_of_lines", &Block::number_of_lines)
            .def("get_line", &Block::get_line);

    class_<Unit>("Unit", init<string, string, double>());

    class_<UnitSet>("UnitSet", init<Unit, Unit, Unit>())
            .def("get_length_unit", &UnitSet::get_length_unit)
            .def("get_density_unit", &UnitSet::get_density_unit)
            .def("get_velocity_unit", &UnitSet::get_velocity_unit)
            .def("get_pressure_unit", &UnitSet::get_pressure_unit)
            .def("get_time_unit", &UnitSet::get_time_unit)
            .def("get_mass_unit", &UnitSet::get_mass_unit)
            .def("get_unit", &UnitSet::get_unit);

    class_<AsciiOutput>("AsciiOutput", init<string>())
            .def("write", &AsciiOutput::write);

    class_<FileOutput>("FileOutput", init<string>())
            .def("write", &FileOutput::write)
            .def("write_header", &FileOutput::write_header);

    class_<UnitConverter>("UnitConverter", init<Unit, Unit>())
            .def("convert", &UnitConverter::convert);
}
