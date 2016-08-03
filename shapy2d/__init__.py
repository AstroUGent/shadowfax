################################################################################
# This file is part of Shadowfax
# Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# Shadowfax is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Shadowfax is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Shadowfax. If not, see <http://www.gnu.org/licenses/>.
################################################################################

import sys
if "shapy3d" in sys.modules:
  raise ImportError("3D version of shapy already loaded!\n"
                    "Loading both the 2D and 3D version leads to conflicts!")

from libpython_riemann2d import *
from libpython_io2d import *
from libpython_grid2d import *
from libpython_lloyd2d import *

# by specifying the units in a separate file, we isolate them in their own
# namespace
import si_units

# workaround to make a Header object printable
# this should be possible in C++, but it gave very consistent compile errors
def header_as_string(header):
  string = "Header:\n"
  string += "npart: [" + str(header.ngaspart()) + "," + str(header.ndmpart()) \
            + "]\n"
  string += "time: " + str(header.time())
  return string

# override the default string method of Header
Header.__str__ = header_as_string

def vec_as_string(vec):
  string = "(" + str(vec.x()) + "," + str(vec.y()) + ")"
  return string

Vec.__str__ = vec_as_string
