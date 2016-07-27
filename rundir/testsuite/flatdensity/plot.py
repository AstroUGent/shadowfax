#! /usr/bin/python

################################################################################
# This file is part of Shadowfax
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

import numpy as np
import h5py
import pylab as pl
import glob

t2d = []
E2d = []
for f in sorted(glob.glob("2d/flatdensity*.hdf5")):
  file = h5py.File(f, "r")
  t2d.append(file["/Header"].attrs["Time"])
  u = np.array(file["/PartType0/InternalEnergy"])
  E2d.append(np.mean(u))

t3d = []
E3d = []
for f in sorted(glob.glob("3d/flatdensity*.hdf5")):
  file = h5py.File(f, "r")
  t3d.append(file["/Header"].attrs["Time"])
  u = np.array(file["/PartType0/InternalEnergy"])
  E3d.append(np.mean(u))

pl.plot(t2d, E2d, "rx", markersize = 10, label = "2D")
pl.plot(t3d, E3d, "bo", label = "3D")
pl.legend(loc = "best")
pl.tight_layout()
pl.show()
