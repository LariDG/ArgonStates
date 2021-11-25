"""
 CMod Project B: Phase Diagram of Argon 
 Module including radial distribution function and mean squared displacement
 Authors: L. Dorman-Gajic and M. van Laar
"""

import sys
import math
import numpy as np
import MDUtilities as mdu
import matplotlib.pyplot as plt
from particle3D import Particle3D
import PBC as pbc
import forces_energies as fe

""" Mean squared displacement function """
def msd(particles, L):
	msd = 0.0
	for i in range(len(particles)):
		separation = np.linalg.norm(particles[i].pos - particles[i].inipos)
		separation = pbc.min_image(separation, L)
		msd += (1/len(particles))*separation**2

	return msd

""" Plot mean squared displacement as a function of time """
def plot_msd(time_list, msd):
	plt.title("MSD averaged over all particles against time")
	plt.xlabel("Time")
	plt.ylabel("Mean squared displacement")
	plt.plot(time_list, msd)
	plt.show()


""" Write MSD values to output file """
def msd_output(msd, time_list):
	file_handle = open("msd.txt", "w")
	zip(msd, time_list)
	file_handle.write(str(msd))

""" Radial distribution function """
def radial_distribution(numstep, particles, L, rdf):
	number_loop = 0

	for i in range(0, len(particles)):
		for j in range(i+1, len(particles)):
			R = np.linalg.norm(Particle3D.vector_separation(particles[i], particles[j], L))
			rdf.append(R)
			number_loop += 1

""" Write RDF values to output file """
def rdf_output(rdf):
	file_handle = open("rdf.txt", "w")
	file_handle.write(str(rdf))

""" Plot the radial distribution function """
def plot_rdf(rdf):
	n_bins = 100
	plt.hist(rdf, n_bins, histtype = "step", density = "True", facecolor = "red")
	plt.title("RDF against radial distance")
	plt.xlabel("Radial distance")
	plt.ylabel("RDF")
	plt.show()

""" Kinetic, potential and total energy of system """
def sys_energies(particles, L, mass, energy_list, energy_K, energy_P):
	energy_kinetic = 0
	energy_potential = 0
	energy = 0
			
	energy_potential = sum(fe.pot_E(particles, L))
	energy_kinetic = sum(fe.kin_E(particles, mass))

	energy = energy_kinetic + energy_potential

	# Append energies to lists defined created main
	energy_list.append(energy)
	energy_K.append(energy_kinetic)
	energy_P.append(energy_potential)

""" Plot energies as a function of time """
def plot_energy(time_list, energy_list, energy_P, energy_K):
	plt.title("Kinetic, potential and total energy of system against time")
	plt.xlabel("Time") #Correct units?
	plt.ylabel("Energy") #Correct units?
	plt.plot(time_list, energy_list, label='Total energy')
	plt.plot(time_list, energy_P, label='Potential energy')
	plt.plot(time_list, energy_K, label='Kinetic energy')
	plt.legend()
	plt.show()

""" Write energy values to output file """
def energy_output(energy_list, time_list):
	file_handle = open("energy.txt", "w")
	file_handle.write(str(energy_list))
