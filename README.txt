Larisa Dorman-Gajic & Margot van Laar, March 2019
Computer Modelling 
Project B - Phase Diagram of Argon

This is a program to simulate the solid, liquid, and gas phase of Argon atoms using Lennard Jones pair potential and a velocity Verlet time step algorithm. 

The main body of the program is in the module N-body.py. This inputs the class Particle3D, creating particles for the format to allow for them to be simulated in the simulator VMD.

The reduced particle density, temperature, number of particles, the force calculation cut of distance, the time step, and number of time steps are inputted into the main via a text file called argon_solid.txt, argon_liquid.txt, argon_gas.txt, depending which phase one wishes to simulate.

A module MDUtilities.py sets initial positions, velocities and sets the box length of the particles created by Particle3D. 
PBC.py is a module containing methods for periodic boundary conditions and minimum image convection. 
Observables.py has methods calculating the mean square displacement and radial distribution of the particles. 
forces_energies.py calculates the forces and energies of the particles. 

The program outputs 3 graphs, one containing potential, kinetic and total energy against time, another the mean squared displacement against time and the third the radial distribution of particles plotted as a histogram. 
Respectively, 3 text files are outputted; energy.txt, msd.txt and ref.txt.

A file called trajectory.xyz is outputted with the particles positions over time in a format compatible with VMD to allow for simulation.

The units for all variables in the program are reduced units.

