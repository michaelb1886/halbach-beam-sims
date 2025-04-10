import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import magpylib as magpy

separation = 100e-2
d,h = 30e-3,100e-3

mu0 = 4 * np.pi * 1e-7
m0 = 1.3/mu0
m1 = magpy.magnet.Cylinder(
    magnetization=(0, m0, 0),
    dimension=(d,h),
    position=(0, -separation/2, 0),)

m2 = magpy.magnet.Cylinder(
    magnetization=(0, m0, 0),
    dimension=(d,h),
    position=(0, +separation/2, 0),)

# m2.rotate_from_angax(180, 'z')

# Create a magnet system
magnet_system = magpy.Collection(m1, m2)

# Create a grid of points in space
def plot_magnetic_field(magnet_system):
    ngrid = 128
    x = np.linspace(-10, 10, ngrid)
    y = np.linspace(-10, 10, ngrid)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    # Create a figure and axis

    b = magnet_system.getB(
        np.array((X.flatten(), Y.flatten(), Z.flatten())).T,
    ).reshape(ngrid,ngrid,3)

    c = plt.pcolormesh(X, Y, np.log(np.linalg.norm(b, axis=2)), shading='auto')
    plt.colorbar(c)
    plt.show()

def rotate_both_magnets(magnet_system, p=[1,0,0], ntheta=128, angle0=0):
    # Rotate both magnets around the z-axis
    b = []
    thetas = np.linspace(0, 2*np.pi, ntheta)
    deltatheta = thetas[1] - thetas[0]
    m1.rotate_from_angax(angle0, 'z', degrees=False)
    # magnet_system.add(magpy.Sensor(position=(p,0,0)))
    for i in range(len(thetas)):
        m1.rotate_from_angax(deltatheta, 'z', degrees=False)
        m2.rotate_from_angax(deltatheta, 'z', degrees=False)
        b.append(magnet_system.getB(p))
    b = np.array(b)
    plt.figure()
    plt.polar(thetas, b[:,0], label='Bx')
    plt.polar(thetas, b[:,1], label='By')
    plt.polar(thetas, b[:,2], label='Bz')
    plt.legend()
    plt.figure()
    plt.plot(thetas, b[:,0], label='Bx')
    plt.plot(thetas, b[:,1], label='By')
    plt.plot(thetas, b[:,2], label='Bz')
    plt.show()

def rotate_one_magnet(magnet_system, p=[1,0,0], ntheta=128):
    # Rotate both magnets around the z-axis
    b = []
    thetas = np.linspace(0, 2*np.pi, ntheta)
    deltatheta = thetas[1] - thetas[0]
    # magnet_system.add(magpy.Sensor(position=(p,0,0)))
    # magnet_system.remove(m2)
    for angle in thetas:
        m1.rotate_from_angax(deltatheta, 'z', degrees=False)
        b.append(magnet_system.getB(p))
    b = np.array(b)
    plt.figure()
    plt.polar(thetas, b[:,0], label='Bx')
    plt.polar(thetas, b[:,1], label='By')
    plt.polar(thetas, b[:,2], label='Bz')
    plt.legend()
    plt.figure()
    plt.plot(thetas, b[:,0], label='Bx')
    plt.plot(thetas, b[:,1], label='By')
    plt.plot(thetas, b[:,2], label='Bz')
    plt.show()

# rotate_one_magnet(magnet_system)
# rotate_both_magnets(magnet_system)
p = [0.3,0.7,0.1]
# print(p)
# counter_rotate_both_magnets(magnet_system,p=p)
rotate_both_magnets(magnet_system,p=p,angle0=np.pi*0)