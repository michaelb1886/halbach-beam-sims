import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import magpylib as magpy
from scipy.signal import periodogram


# Create a grid of points in space
def plot_magnetic_field(magnet_system, L=10):
    ngrid = 128
    x = np.linspace(-L, L, ngrid)
    y = np.linspace(-L, L, ngrid)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    # Create a figure and axis

    b = magnet_system.getB(
        np.array((X.flatten(), Y.flatten(), Z.flatten())).T,
    ).reshape(ngrid,ngrid,3)

    c = plt.pcolormesh(X, Y, np.log(np.linalg.norm(b, axis=2)), shading='auto')
    plt.colorbar(c)
    # plt.show()

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
    make_polar_plots(np.array(b), thetas)

def rotate_one_magnet(magnet_system, p=[1,0,0], ntheta=128):
    # Rotate both magnets around the z-axis
    b = []
    thetas = np.linspace(0, 4*np.pi, ntheta)
    deltatheta = thetas[1] - thetas[0]
    # magnet_system.add(magpy.Sensor(position=(p,0,0)))
    # magnet_system.remove(m2)
    for angle in thetas:
        m1.rotate_from_angax(deltatheta, 'z', degrees=False)
        b.append(magnet_system.getB(p))
    make_polar_plots(np.array(b), thetas)

def rotate_both_magnets_one_double(magnet_system, p=[1,0,0], ntheta=256, angle0=0):
    # Rotate both magnets around the z-axis
    b = []
    thetas = np.linspace(0, 20*np.pi, ntheta)
    deltatheta = thetas[1] - thetas[0]
    m1.rotate_from_angax(angle0, 'z', degrees=False)
    # magnet_system.add(magpy.Sensor(position=(p,0,0)))
    for i in range(len(thetas)):
        m1.rotate_from_angax(deltatheta, 'z', degrees=False)
        m2.rotate_from_angax(-deltatheta, 'z', degrees=False)
        b.append(magnet_system.getB(p))
    make_polar_plots(np.array(b), thetas)

def polar_sweep_plot(magnet_system, R=10, ntheta=128):
    # Rotate both magnets around the z-axis
    b = []
    thetas = np.linspace(0, 2*np.pi, ntheta)
    deltatheta = thetas[1] - thetas[0]
    # magnet_system.add(magpy.Sensor(position=(p,0,0)))
    for angle in thetas:
        p = [R*np.cos(angle), R*np.sin(angle), 0]
        b.append(magnet_system.getB(p))
    b = np.array(b)
    plt.figure()
    for i in range(3):
        plt.polar(thetas, b[:,i], label=f'B{i}')
    plt.plot(thetas, np.linalg.norm(b,axis=1), label='|B|')
    plt.legend()
    # make_polar_plots(np.array(b), thetas)

def make_polar_plots(b, thetas):
    plt.figure()
    for i in range(3):
        plt.polar(thetas, b[:,i], label=f'B{i}')
    plt.legend()
    plt.figure()
    for i in range(3):
        plt.plot(thetas, b[:,i], label=f'B{i}')
    plt.legend()
    plt.figure()
    for i in range(3):
        pgm = periodogram(b[:,i], fs=1/(thetas[1]-thetas[0]))
        plt.plot(pgm[0], pgm[1], label=f"B{i}")
    # plt.show()


separation = 2.5e-2
d,h = 5e-2,5e-2
# mx,my,mz = 5e-2,5e-2,5e-2

mu0 = 4 * np.pi * 1e-7
m0 = 1.3/mu0
m1 = magpy.magnet.Cylinder(
    magnetization=(0, m0, 0),
    # dimension=(mx,my,mz),
    dimension = (d,h),
    position=(0, -separation/2, 0),)

m2 = magpy.magnet.Cylinder(
    magnetization=(0, m0, 0),
    # dimension=(mx,my,mz),
    dimension = (d,h),
    position=(0, +separation/2, 0),)

# m2.rotate_from_angax(180, 'z')

# Create a magnet system
magnet_system = magpy.Collection(m1, m2)

# rotate_one_magnet(magnet_system)
# rotate_both_magnets(magnet_system)
# p = np.array([0.3,0.7,0.1])*5
p = np.random.randn(3)*30
# print(p)
# counter_rotate_both_magnets(magnet_system,p=p)
rotate_both_magnets(magnet_system,p=p,angle0=np.pi*0)
rotate_both_magnets(magnet_system,p=p,angle0=np.pi*1)
# plot_magnetic_field(magnet_system, L=10)
# polar_sweep_plot(magnet_system, R=100, ntheta=128)
# m1.rotate_from_angax(180, 'z', degrees=True)
# polar_sweep_plot(magnet_system, R=100, ntheta=128)
# # plot_magnetic_field(magnet_system, L=10)
plt.show()