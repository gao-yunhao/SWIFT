############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

# Parameters
rho = 1.0
# cs = 0.001 # 55
cs2 = 3025.0  # cs**2 #3025.0
L = float(sys.argv[6])  # 1.0
kv = int(sys.argv[4])
kv0 = 2 * np.pi / L * kv
kb = int(sys.argv[5])
kb0 = 2 * np.pi / L * kb
V0 = float(sys.argv[1])  # 22.9
Vz_factor = float(sys.argv[7])
resistive_eta = float(sys.argv[2])
Beq0 = np.sqrt(rho) * V0
B0 = 1e-3 * Beq0
gamma = 5.0 / 3.0
u0 = cs2 / (gamma * (gamma - 1))

# output file
fileOutputName = "RobertsFlow.hdf5"

###---------------------------###

glass = h5py.File(sys.argv[3], "r")
pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

N = len(h)
vol = L ** 3

###---------------------------###

v = np.zeros((N, 3))
B = np.zeros((N, 3))
A = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho * vol / N
u = np.ones(N) * u0

# rescaling the box to size L
pos *= L

# setting up flow
# v[:, 0] = np.sin(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])
# v[:, 1] = -np.sin(kv0 * pos[:, 0]) * np.cos(kv0 * pos[:, 1])
# v[:, 2] = Vz_factor * np.sqrt(2) * np.cos(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])

v[:, 0] = 2 * (np.random.rand(N) - 0.5)
v[:, 1] = 2 * (np.random.rand(N) - 0.5)
v[:, 2] = 2 * (np.random.rand(N) - 0.5)

# Average rms velocity of such field configuration is 1, everything is ok
# vv = np.sqrt(np.mean(v[:,0]**2+v[:,1]**2+v[:,2]**2)), therefore rms in physical field will be Vrms=V0
# print(vv)

v *= V0 * 0.01

B[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
B[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
B[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
B *= B0

A[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
A[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
A[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
A0 = B0 / kb0
A *= A0

###---------------------------###

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]  #####
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")

fileOutput.close()
