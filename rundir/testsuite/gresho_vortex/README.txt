= The Gresho vortex =

This problem consists of a cubic box filled with a homogeneous fluid of density
one. The fluid inside a sphere of radius r = 0.4 w.r.t. the center (0.5, 0.5)
is rotating with an azimuthal velocity profile
        { 5r   for 0.  <= r < 0.2
v_phi = {  2-5r for 0.2 <= r < 0.4
        { 0    for 0.4 <= r
The centrifugal force caused by this rotation is balanced by a pressure profile
    { 5 + 12.5r^2                    for 0.  <= r < 0.2
p = { 9 + 12.5r^2 - 20r + 4ln(r/0.2) for 0.2 <= r < 0.4
    { 3 + 4ln(2)                     for 0.4 <= r
such that the solution is static, i.e. independent of time.

This test demonstrates the local conservation of quantities. Although the
integration method used guarantees total conservation of mass, momentum and
energy for the entire system, the same does not hold for local cell quantities,
even for a static problem. However, Shadowfax does quite a good job in keeping
the velocity profile constant in time, which means that local quantities are
indeed conserved to a good accuracy.

The initial conditions for this test are based on the Gresho vortex problem in
Springel (2010).

To run this problem, run
bash run.sh

The solution, together with the static velocity profile, can be plotted using
python plot.py

== References ==

Springel, V. 2010, MNRAS, 401, 791

-- Written by Bert Vandenbroucke, Tuesday 12 August 2014 --
