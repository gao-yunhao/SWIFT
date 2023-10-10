import numpy as np
import h5py
import matplotlib.pyplot as plt

snaps = 12

e_kin = np.zeros(snaps)
e_therm = np.zeros(snaps)
e_mag = np.zeros(snaps)
e_grav = np.zeros(snaps)

e_tot = np.zeros(snaps)

for snap in range(snaps):

    sim = h5py.File("noh_%04d.hdf5" % snap, "r")

    x = sim["/PartType0/Coordinates"][:, 0]
    y = sim["/PartType0/Coordinates"][:, 1]
    z = sim["/PartType0/Coordinates"][:, 2]
    vx = sim["/PartType0/Velocities"][:, 0]
    vy = sim["/PartType0/Velocities"][:, 1]
    vz = sim["/PartType0/Velocities"][:, 2]
    Bx = sim["/PartType0/MagneticFluxDensities"][:, 0]
    By = sim["/PartType0/MagneticFluxDensities"][:, 1]
    Bz = sim["/PartType0/MagneticFluxDensities"][:, 2]
    u = sim["/PartType0/InternalEnergies"][:]

    print(u)
    quit()

    rho = sim["/PartType0/Densities"][:]
    phi = sim["/PartType0/Potentials"][:]

    e_kin_i = 0.5 * (vx ** 2 + vy ** 2 + vz ** 2)
    e_therm_i = u
    e_mag_i = 0.5 * (Bx ** 2 + By ** 2 + Bz ** 2) / rho
    e_grav_i = phi

    e_kin[snap] = np.sum(e_kin_i)
    e_therm[snap] = np.sum(e_therm_i)
    e_mag[snap] = np.sum(e_mag_i)
    e_grav[snap] = np.sum(e_grav_i)

e_tot = e_kin + e_therm + e_mag + e_grav

plt.semilogy(e_kin / e_kin[0], "-g", label="e_kin")
# plt.semilogy(e_therm/e_therm[0], '-p', label="e_therm")
plt.semilogy(e_mag / e_mag[0], "-r", label="e_mag")
# plt.plot(e_grav/e_grav[0], '-b', label="e_grav")
plt.semilogy(e_tot / e_tot[0], "-k", label="e_tot")
plt.legend()
plt.savefig("test.png")
