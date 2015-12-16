== Restart test ==

This test is the same as the spherical overdensity test (in 3D even with lower
resolution), and is only used to test the codes capabilities to restart a run
from restart files.

The test is always run with 3 MPI processes, to test the parallel restart
capabilities as well. We first run the entire simulation up to the end, with
a small time interval in between restart file dumps, so that we have at least
one intermediate restart file. After the simulation has finished, we copy the
final snapshot for future reference. We then restart the simulation from the
last restart dumps and run it again until the end. Since restarting a run should
not have any influence on the result, we can then compare the final snapshot
with the reference snapshot to see if restarting indeed works.

Due to stochastic effects, there will be small differences in parallel runs.
These are filtered out by setting the appropriate tolerance parameter for
h5diff.

-- Written by Bert Vandenbroucke, Thursday 3 December 2015 --
