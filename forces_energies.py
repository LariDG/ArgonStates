
import sys
import math
import numpy as np
from particle3D import Particle3D

""" Calculating Lennard-Jones force between two particles """
def force_LJ(p1,p2, L, r_c):
	r = Particle3D.vector_separation(p1, p2, L)

	# Compute force if particle is closer than cut-off radius
	if np.linalg.norm(r) <= r_c:
		if np.linalg.norm(r) == 0.0:
			force = 0.
		else: 
			force = 48*((1/(np.linalg.norm(r))**14)-(1/(2*(np.linalg.norm(r))**8)))*r
	else: 
		force = np.array([0.0,0.0,0.0])
	
	return force

""" Calculating the sum of all forces acting on all particles """
def force_init(particles, L, r_c):
	force_list = []
	for i in range(0, len(particles)):
		force_list.append(np.zeros(3))

	for i in range(0, len(particles)):
		for j in range(i+1, len(particles)):
			F = force_LJ(particles[i], particles[j], L, r_c)
			force_list[i] += F
			force_list[j] -= F

	return force_list
	

""" Calculating potential energy of two interacting particles """
def potential_energy(p1, p2, L):
	r = Particle3D.vector_separation(p1, p2, L)
	if np.linalg.norm(r) == 0:
		pot = 0
	else:
		pot = 4*((1/(np.linalg.norm(r))**12)-(1/(np.linalg.norm(r))**6))

	return pot

""" Total potential energy of a particle due to its interaction with all others """
def pot_E(particles, L):
	potentials = []
	for i in range(0, len(particles)):
		potentials.append(0)
	
	for i in range(0, len(particles)):
		for j in range(i+1, len(particles)):
			P = potential_energy(particles[i], particles[j], L)
			potentials[i] += P

	return potentials

""" Kinetic energy of a particle """
def kinetic_energy(p1, m):
	"""
	Return kinetic energy as
	1/2*mass*vel^2
	"""
	v = np.linalg.norm(p1.vel)
	
	return 0.5*m*v**2.0

""" Creating list of kinetic energies of all particles """
def kin_E(particles, m):
	kinetic = []
	for i in range(0, len(particles)):
		kinetic.append(0)
	for i in range(0,len(particles)):
		energy_kinetic = kinetic_energy(particles[i], m)
		kinetic.append(energy_kinetic)

	return kinetic