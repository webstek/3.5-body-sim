# 3.5-body-sim
The end goal is to create a simulation that makes use of JPL's Horizons Ephemeris to calculate launch windows for Earth to Mars transfers. The inputs would be an arrival date and the maximum delta V. The program is to determine if there is a trajectory that satisfies the input parameters and what the launch window is for it.

Currently, I have only completed the physics for the 3.5 body simulation - although I don't think the current iterative method conserves energy and momentum. Which can be fixed by using the Verlet velocity method.
