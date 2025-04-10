import numpy as np
import matplotlib.pyplot as plt
import magpylib as magpy
from scipy.spatial.transform import Rotation as R

def dipole(m, r, r0):
    """Calculate a field in point r created by a dipole moment m located in r0.
    Spatial components are the outermost axis of r and returned B.
    """
    # we use np.subtract to allow r and r0 to be a python lists, not only np.array
    R = np.subtract(np.transpose(r), r0).T
    
    # assume that the spatial components of r are the outermost axis
    norm_R = np.sqrt(np.einsum("i...,i...", R, R))
    
    # calculate the dot product only for the outermost axis,
    # that is the spatial components
    m_dot_R = np.tensordot(m, R, axes=1)

    # tensordot with axes=0 does a general outer product - we want no sum
    B = 3 * m_dot_R * R / norm_R**5 - np.tensordot(m, 1 / norm_R**3, axes=0)
    
    # include the physical constant
    B *= 1e-7

    return B


# Define the parameters of the magnets array
r = 0.015  # radius in meters
h = 0.1  # height in meters

Br = 1.2  # remanent magnetization in Tesla
mu = 4 * np.pi * 1e-7  # permeability of free space in T*m/A    

V = np.pi * r**2 * h  # volume of the magnet in m^3
m = np.array([0,Br * V / mu,0])  # magnetic moment in Am^2

n_magnets = 11  # number of magnets in the array

angle = np.pi/n_magnets  # angle between magnets in radians

moments = []
positions = []

angle0 = 0*2*np.pi/3
r_array = 200e-2

positions = np.vstack((np.linspace(-r_array, r_array, n_magnets), np.zeros(n_magnets), np.zeros(n_magnets))).T

for i in range(n_magnets):
    # Create a cylinder magnet with the specified parameters
    rot = R.from_euler('z', 2*i * angle + angle0, degrees=False)
    # positions.append(np.array([r * np.cos(i * angle), r * np.sin(i * angle), 0]))
    moments.append(rot.apply(m))

positions = np.array(positions)
moments = np.array(moments)

r_sensor = 100.0
thetas = np.linspace(0, 2 * np.pi, 64)
r_sensor = np.array([[r_sensor * np.cos(theta), r_sensor * np.sin(theta), 0] for theta in thetas])

def N_dipoles(r, positions, moments):
    # Create a list of dipoles
    B = np.array([0., 0., 0.])
    for i in range(len(positions)):
        B += dipole(moments[i], r, positions[i])
    return B

# calculate the field at the sensor positions
B = np.array([N_dipoles(r, positions+np.array([2*r_array,0.,0.]), moments) for r in r_sensor])
B += np.array([N_dipoles(r, positions+np.array([-2*r_array,0.,0.]), moments) for r in r_sensor])
B = np.array([np.linalg.norm(b) for b in B])
plt.polar(thetas,(B), marker='o')

B = np.array([N_dipoles(r, positions[0:1], moments[0:1]) for r in r_sensor])
B = np.array([np.linalg.norm(b) for b in B])
plt.polar(thetas, (B), marker='o')
plt.show()