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

import libpython_io3d as io

unit_length = io.Unit("length", "m", 1.)
unit_mass = io.Unit("mass", "kg", 1.)
unit_time = io.Unit("time", "s", 1.)

unitset = io.UnitSet(unit_length, unit_mass, unit_time)

unit_density = unitset.get_density_unit()
unit_velocity = unitset.get_velocity_unit()
unit_pressure = unitset.get_pressure_unit()
