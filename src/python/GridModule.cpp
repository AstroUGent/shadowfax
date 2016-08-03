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
 * @file GridModule.cpp
 *
 * @brief Expose some grid building functions to Python
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "../utilities/Cuboid.hpp"
#include "../utilities/GasParticle.hpp"
#include "../utilities/Tree.hpp"
#include "DelCont.hpp"
#include "Vec.hpp"
#include "VorCell.hpp"
#include "VorFace.hpp"
#include "VorGen.hpp"
#include "VorTess.hpp"
#include <string>
#include <vector>
using namespace std;

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/module.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#if ndim_ == 3
/*! @brief Register the libgrid Python module */
BOOST_PYTHON_MODULE(libpython_grid3d)
#else
BOOST_PYTHON_MODULE(libpython_grid2d)
#endif
{
    class_<vector<VorFace*> >("vector<VorFace*>")
            .def(vector_indexing_suite<vector<VorFace*> >());

    class_<vector<VorGen*> >("vector<VorGen*>")
            .def(vector_indexing_suite<vector<VorGen*> >());

    class_<Particle, boost::noncopyable>("Particle", no_init);

#if ndim_ == 3
    class_<GasParticle, boost::noncopyable, bases<Particle> >("GasParticle",
                                                              init<Vec>())
            .def("x", &GasParticle::x)
            .def("y", &GasParticle::y)
            .def("z", &GasParticle::z)
            .def("set_x", &GasParticle::set_x)
            .def("set_y", &GasParticle::set_y)
            .def("set_z", &GasParticle::set_z)
            .def("set_key", &GasParticle::set_key)
            .def("get_id", &GasParticle::id)
            .def("set_id", &GasParticle::set_id);
#else
    class_<GasParticle, boost::noncopyable, bases<Particle> >("GasParticle",
                                                              init<Vec>())
            .def("x", &GasParticle::x)
            .def("y", &GasParticle::y)
            .def("set_x", &GasParticle::set_x)
            .def("set_y", &GasParticle::set_y)
            .def("set_key", &GasParticle::set_key)
            .def("get_id", &GasParticle::id)
            .def("set_id", &GasParticle::set_id);
#endif

#if ndim_ == 3
    class_<Vec>("Vec", init<double, double, double>())
            .def("x", &Vec::x)
            .def("y", &Vec::y)
            .def("z", &Vec::z);
#else
    class_<Vec>("Vec", init<double, double>())
            .def("x", &Vec::x)
            .def("y", &Vec::y);
#endif

    class_<CubicBox>("CubicBox", init<Vec, double>())
            .def("get_key", &CubicBox::get_key);

    class_<VorTess>("VorTess", init<CubicBox*, unsigned int>())
            .def("add_point", &VorTess::add_point)
            .def("get_cell", &VorTess::get_cell,
                 return_value_policy<reference_existing_object>())
            .def("complete", &VorTess::complete)
            .def("construct", &VorTess::construct);

#if ndim_ == 3
    class_<VorGen, VorGen*>("VorGen", init<double, double, double>())
            .def("x", &VorGen::x)
            .def("y", &VorGen::y)
            .def("z", &VorGen::z)
            .def("get_particle", &VorGen::get_particle,
                 return_value_policy<reference_existing_object>())
            .def("get_cell", &VorGen::get_cell,
                 return_value_policy<reference_existing_object>());
#else
    class_<VorGen, VorGen*>("VorGen", init<double, double>())
            .def("x", &VorGen::x)
            .def("y", &VorGen::y)
            .def("get_particle", &VorGen::get_particle,
                 return_value_policy<reference_existing_object>())
            .def("get_cell", &VorGen::get_cell,
                 return_value_policy<reference_existing_object>());
#endif

    class_<VorCell>("VorCell", init<VorGen*>())
            .def("get_volume", &VorCell::get_volume)
            .def("get_centroid", &VorCell::get_centroid,
                 return_value_policy<reference_existing_object>())
            .def("get_faces", &VorCell::get_faces)
            .def("get_ngbs", &VorCell::get_ngbs);

    class_<Cuboid>("Cuboid", init<Vec, Vec>());

    class_<Tree>("Tree", init<Cuboid>())
            .def("add_particle", &Tree::add_particle)
            .def("get_closest", &Tree::get_closest,
                 return_value_policy<reference_existing_object>())
            .def("finalize", &Tree::finalize);

    class_<VorFace, VorFace*>("VorFace", init<unsigned, vector<VorGen*>&>())
            .def("get_vertices", &VorFace::get_vertices)
            .def("get_midpoint", &VorFace::get_midpoint,
                 return_value_policy<reference_existing_object>())
            .def("get_facengbs", &VorFace::get_facengbs)
            .def("get_area", &VorFace::get_area);
}
