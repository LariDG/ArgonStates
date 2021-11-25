"""
 CMod Project B: Phase Diagram of Argon 
 Module including periodic boundary contitions function and minimum image covention function
 Authors: L. Dorman-Gajic and M. van Laar
"""
import sys
import math
import numpy as np


""" Function to correct position of particle according to periodic
boundary conditions
"""
def pbc(p1, L):
	p1.pos[0] = p1.pos[0]%L
	p1.pos[1] = p1.pos[1]%L
	p1.pos[2] = p1.pos[2]%L

	return p1.pos

""" Function to correct separation of two particles according to minimum
image convention 
"""
def min_image(separation, L):
	separation = np.mod(separation + L/2.0, L) - L/2.0

	return separation