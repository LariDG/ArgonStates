"""
 CMod Project B: Phase Diagram of Argon 
 Main method
 Authors: L. Dorman-Gajic and M. van Laar
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from particle3D import Particle3D
import MDUtilities as mdu
import PBC as pbc
import forces_energies as fe
import observables as obs

""" Writing positions of all particles at a given time to an output file """
def VMD_output(file_handle, particles, t):
	file_handle.write(str(len(particles)) + "\n")
	file_handle.write("Point = " + str(t) + "\n")

	for i in range(len(particles)):
		file_handle.write(str(particles[i]) + "\n")


def main():

# Defining system parameters
	file_handle = sys.argv[1]
	infile = open(file_handle, "r")
	constant = infile.readline()
	constant_list = constant.split()
	rho = float(constant_list[0])
	temp = float(constant_list[1])
	N = int(constant_list[2])
	r_c = float(constant_list[3])
	numstep = int(constant_list[4])
	dt = float(constant_list[5])
	mass = 1

	# Create open list of particles of length N
	particles = []

	# Initialise particle instance, which will be updated by MDUtilities
	for i in range(N): 
		p_i = Particle3D.create_particle(i)
		particles.append(p_i)
	
	# Open file to write particle trajectories to
	file_handle = open("trajectory.xyz", "w")

	# Set initial particle positions and velocities
	box_size = mdu.set_initial_positions(rho, particles)
	mdu.set_initial_velocities(temp, particles)

	# Saving initial positions
	for i in range(N):
		particles[i].inipos = particles[i].pos

	# Write particle data at t=0 to VMD file
	VMD_output(file_handle, particles, 0)

	# Length of the box of particles
	L = box_size[0]

	# Open list for total energy of system
	energy_list = []

	# Open list for kinetic and potential energy of particles
	energy_K = []
	energy_P = []

	# Open list to for mean squared displacement and radial distribution data
	msd = []
	rdf = []

	# Total force acting on each particle
	force_list = fe.force_init(particles, L, r_c)

	# Velocity Verlet time integration
	time = 0.0
	time_list = []

	for t in range(1, numstep+1): 

		# Update positions, adhering to periodic boundary conditions
		for i in range(N):
			particles[i].secondOstep_pos(dt, force_list[i])
			#particles[i].pos = pbc.pbc(particles[i], L)
		
		# Update force
		f_new = fe.force_init(particles, L, r_c)

		# Update velocities
		for i in range(0,N):
			particles[i].step_velocity(dt, 0.5*(force_list[i]+f_new[i])) 

		# Update time and force
		time += dt
		force_list = f_new

		if t%1==0:
			# Total kinetic and potential and sum of energies of system
			obs.sys_energies(particles, L, mass, energy_list, energy_K, energy_P)

			# MSD data at time step
			msd.append(obs.msd(particles, L))

			# Radial distribution at time step
			obs.radial_distribution(numstep, particles, L, rdf)
			time_list.append(time)
			
		# Writing data to file for VMD
		VMD_output(file_handle, particles, t)

	# Writing data for msd, rdf and energy to output file
	obs.rdf_output(rdf)
	obs.msd_output(msd, time_list)
	obs.energy_output(energy_list, time_list)
	
	#Plotting MSD as function of time
	obs.plot_msd(time_list, msd)

	# Plotting energy as a function of time
	obs.plot_energy(time_list, energy_list, energy_P, energy_K)

	#Plotting the radial distribution function

	obs.plot_rdf(rdf)

main()

