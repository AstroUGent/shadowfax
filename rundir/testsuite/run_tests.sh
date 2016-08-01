#! /bin/bash

################################################################################
# This file is part of Shadowfax
# Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2016 Bert Vandenbroucke
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

# registered tests
# every test should be in a subfolder of this folder with the same name as the
# test
# this folder should contain a ./run.sh script with the executable bit set that
# contains all necessary commands to run the test and returns a useful value on
# exit (0 begin success, all other values being a failed test)
testnames=( overdensity nbody evrard gresho_vortex sedov_taylor restarttest \
            flatdensity vacuum )
tests=()
timings=()

for test in "${testnames[@]}"
do
cd $test

echo "Running \"$test\" test"

starttime=$(date +%s.%N)
./run.sh
status=$?
endtime=$(date +%s.%N)
timediff=$(echo "$endtime - $starttime" | bc)
timings+=($timediff)
tests+=($status)

if [ $status -ne 0 ]
then
echo -e "Test \"$test\" \e[31mfailed\e[0m!"
else
echo -e "Test \"$test\" \e[32mpassed\e[0m."
fi

cd ..
done

# summary
echo "### TEST SUMMARY ###"
testcount=$((${#testnames[@]} - 1))
for i in $(eval echo {0..$testcount})
do
flag=${tests[$i]}
testname=${testnames[$i]}
timing=${timings[$i]}
if [ $flag -ne 0 ]
then
echo -e "Test \"$testname\" \e[31mfailed\e[0m!"
else
echo -e "Test \"$testname\" \e[32mpassed\e[0m."
fi
echo $timing
done
