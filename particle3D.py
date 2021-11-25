
"""
 CMod Project B: Phase Diagram of Argon 
 Class to describe point-like particles moving in 3D space

 Authors: L. Dorman-Gajic and M. van Laar
"""

import sys
import math
import numpy as np
import MDUtilities as mdu
import PBC as pbc


class Particle3D(object):
	"""
	Class to describe point-like particles moving in 3D space.

	Properties:
	position(array) - positition vector of particle, r
	velocity(array) - velocity vector of particle, v
	mass(float) - mass of particle

	Methods:
	* initialisation of particle object
	* first order velocity update for given timestep and force vector
	* second order update position update for given timestep and force vector
	* create particle from file entry (static)
	* relative vector separation of two particles (static)
	"""

	""" Initialise a particle instance """
	def __init__(self, label, mass, pos, vel):
		"""
		Initialise a Particle3D instance

		:param pos: position as numpy array
		:param vel: velocity as numpy array
		:param mass: mass as float
		"""

		self.pos = pos
		self.inipos = pos
		self.vel = vel
		self.mass = mass
		self.label = label

	""" Formatted output of properties """
	def __str__(self):
		"""
		Format correct output for use in VMD
		<label> <x_pos> <y_pos>, <z_pos> 
		"""
		
		return self.label + " " + str(self.pos[0]) + " " + str(self.pos[1]) + " " + str(self.pos[2]) 

	""" Update the velocity to first order """
	def step_velocity(self, dt, force):
		"""
		first order velocity update given by
		v(t+dt) = v(t) + dt*f(t)/mass.       

		:param dt: time step as float
		:param force: force as numpy array
		"""

		self.vel = self.vel + dt*force/self.mass

	""" Update position to second order """
	def secondOstep_pos(self, dt, force):
		"""
		second order position update given by
		r(t+dt) = r(t) + dt*v(t) + dt**2*force/(2*mass)
		:param dt: time step as float
		:param force: force as numpy array
		"""

		self.pos = self.pos + dt*self.vel + (dt**2)*force/(2.0*self.mass)

	""" Calculate the vector separation between paricles, adhering to minimum image convention """
	@staticmethod
	def vector_separation(p1, p2, L): 
		separation = p1.pos - p2.pos
		new_separation = pbc.min_image(separation,L)
		
		return new_separation

	""" Create a particle instance """
	@staticmethod
	def create_particle(i):
		pos = np.array([0,0,0])
		vel = np.array([0,0,0]) 
		mass = 1
		label = str(i)

		return Particle3D(label, mass, pos, vel)



