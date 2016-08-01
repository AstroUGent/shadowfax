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

n=10000
i=1

# $# is the number of command line arguments passed to the script
if [ $# -gt 0 ]
then
n=$1
fi

if [ $# -gt 1 ]
then
i=$2
fi

cd 2d
# generate initial condition
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $i @MPIEXEC_PREFLAGS@ ../../../icmaker2d \
    @MPIEXEC_POSTFLAGS@ --ncell $n --setup ../vacuum2d.xml \
    --filename ic_vacuum2d.hdf5
status=$?

if [ $status -ne 0 ]
then
exit $status
fi

# run simulation
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $i @MPIEXEC_PREFLAGS@ ../../../shadowfax2d \
    @MPIEXEC_POSTFLAGS@ --params ../vacuum2d.ini 2>&1 | tee vacuum2d.log
status=${PIPESTATUS[0]}
cd ..

exit $status
