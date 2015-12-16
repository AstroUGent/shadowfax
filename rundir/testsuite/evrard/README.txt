= The Evrard test =

The Evrard test consists of a self-gravitating cloud with a 1/r density profile
and negligble pressure.
The sphere is centred on the origin (0., 0., 0.) and has radius 1. Inside, the
hydrodynamical variables are:
rho = 1/(2*pi*(r+0.001)), u = 0, v = 0, w = 0, p = 0.05/(3*pi*(r+0.001))
outside the sphere, all variables are set to 0. We use a gravitational constant
G = 1 for this test.

The cloud will start to collapse under its own gravity, causing the pressure in
the center to rise. When the pressure in the center exceeds a critical value, a
strong outward shock is generated that virializes the entire cloud.

This problem is a good test for the coupling between hydrodynamics and gravity,
since gravitational potential energy is converted into thermal energy and then
into kinetic energy. Keeping the total energy constant is challenging.

To run the test, run the command:
bash run.sh

Results can be plotted using the provided python script:
python plot.py

The initial condition was based on Springel (2010).

-- Written by Bert Vandenbroucke, Friday 11 December 2015 --
