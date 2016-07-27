#! /bin/bash

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

testdir="rundir/testsuite/flatdensity"

# compile
cd ../../..
if [ ! -f CMakeCache.txt ];
then
  echo "No build files found. Running cmake without PYTHON_MODULES"
  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=@CMAKE_PREFIX_PATH@ ..
fi
make
cd $testdir
mkdir 2d
mkdir 3d

## 2D simulation

bash run_2d.sh
status=$?

if [ $status -ne 0 ]
then
exit $status
fi

## 3D simulation

bash run_3d.sh
status=$?

exit $status
