= The N-body test =

The N-body test consists of a cold dark matter cloud with a Plummer density
profile that is evolved under the force of gravity. The cloud has a total mass
of 1000 and a scale parameter of 1, so that its density is given by
rho = 3000/4/pi/(1+r^2)^2.5
We set the gravitational constant G to 1 for this test.

The Plummer model is supposed to be stable, so that the density at a later time
should equal the original density profile. We run the test to time 1, which
corresponds to an order of 10 dynamical times of the cloud.

To run the test, run the command:
bash run.sh

Results can be plotted using the provided python script:
python plot.py

-- Written by Bert Vandenbroucke, Friday 11 December 2015 --
