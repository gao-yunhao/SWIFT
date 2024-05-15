import numpy as np
import h5py
import argparse
import unyt
import matplotlib.pyplot as plt

from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
#from swiftsimio.visualisation.slice import project_gas

from unyt import kb, g, mh, cm, s, K, statA

background_rho = 0.1 * mh / (cm ** 3)
background_T = 2550 * (5 / 4) * K  # Equilibrium temperature at solar abundance     
mu = 0.5
background_P = (kb / (mu * mh)) * background_rho  * background_T
background_cs = np.sqrt(5 * background_P / (3 * background_rho))
background_B = 1e-8 * g / (statA * s * s)    

# Parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
#argparser.add_argument("output")
args = argparser.parse_args()

# Load snapshot
filename = args.input
data = load(filename)

# Retrieve some information about the simulation run
artDiffusion = data.metadata.hydro_scheme["Artificial Diffusion Constant"]
dedHyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
dedHypDivv = data.metadata.hydro_scheme["Dedner Hyperbolic div(v) Constant"]
dedPar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
eta = data.metadata.hydro_scheme["Resistive Eta"]
git = data.metadata.code["Git Revision"]
gitBranch = data.metadata.code["Git Branch"]
hydroScheme = data.metadata.hydro_scheme["Scheme"]
kernel = data.metadata.hydro_scheme["Kernel function"]
neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]

# Retrieve particle attributes of interest
v = data.gas.velocities
normv = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)
B = data.gas.magnetic_flux_densities
normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
divB = data.gas.magnetic_divergences

'''
print(normB)
print(max(normB))
quit()
'''

divB[normB.value < 0.01*max(normB.value)] = 0.0 * divB.units

h = data.gas.smoothing_lengths

# Generate mass weighted maps of quantities of interest
data.gas.mass_weighted_densities = data.gas.masses * data.gas.densities # / background_rho

data.gas.mass_weighted_pressures = data.gas.masses * np.log10(data.gas.pressures / background_P)

data.gas.mass_weighted_normv = data.gas.masses * normv / background_cs

data.gas.mass_weighted_normB = data.gas.masses * normB / background_B

data.gas.mass_weighted_error = data.gas.masses * np.log10(np.maximum(h * abs(divB) / normB, 1e-6))

boxsize = data.metadata.boxsize
visualise_region = [
    0, boxsize[0],
    0, boxsize[1],
    0.45 * boxsize[2], 0.55 * boxsize[2]
]

#common_arguments = dict(data=data, resolution=512, parallel=True)
common_arguments = dict(data=data, 
                        #z_slice=0.5 * data.metadata.boxsize[2],
                        region = visualise_region,
                        resolution=512,
                        parallel=True)                          

mass_map = project_gas(**common_arguments, project="masses")

mass_weighted_density_map = project_gas(
    **common_arguments, project="mass_weighted_densities"
)

mass_weighted_pressure_map = project_gas(
    **common_arguments, project="mass_weighted_pressures"
)

mass_weighted_normv_map = project_gas(**common_arguments, project="mass_weighted_normv")

mass_weighted_normB_map = project_gas(**common_arguments, project="mass_weighted_normB")

mass_weighted_error_map = project_gas(**common_arguments, project="mass_weighted_error")

# Take out mass dependence
density_map = mass_weighted_density_map / mass_map
pressure_map = mass_weighted_pressure_map / mass_map
normv_map = mass_weighted_normv_map / mass_map
normB_map = mass_weighted_normB_map / mass_map
error_map = mass_weighted_error_map / mass_map

density_map.convert_to_units(0.1 * mh / (cm ** 3))
#density_map /= background_rho

# Plot maps
plt.rcParams.update({"font.size": 16})
fig, ax = plt.subplots(2, 2, figsize=(12.25, 11))

a00 = ax[0, 0].contourf(
    density_map.value.T, cmap="gist_heat", levels=np.linspace(0.0, 1.5, 100), extend='both'
)
a01 = ax[0, 1].contourf(
    pressure_map.value.T, cmap="gist_heat", levels=np.linspace(0.0, 3.0, 100), extend='both'
)
'''
a10 = ax[1, 0].contourf(
    normv_map.value.T, cmap="gist_heat", levels=np.linspace(0.0, 5.0, 100), extend='both'
)
'''
a10 = ax[1, 0].contourf(
    normB_map.value.T, cmap="gist_heat", levels=np.linspace(0.0, 2.0, 100), extend='both'
)
a11 = ax[1, 1].contourf(error_map.value.T, cmap="jet", levels=np.linspace(-5.0, 0.0, 100), extend='both')

'''
# Add panel with infromation about the run
text_common_args = dict(
    fontsize=10, ha="center", va="center", transform=ax[2, 1].transAxes
)

ax[2, 1].text(
    0.5,
    0.8,
    "Feedback event test at time $t=%.4f$" % data.metadata.time,
    **text_common_args,
)
ax[2, 1].text(0.5, 0.7, "SWIFT %s" % git.decode("utf-8"), **text_common_args)
ax[2, 1].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
ax[2, 1].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
ax[2, 1].text(
    0.5,
    0.4,
    kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
    **text_common_args,
)
ax[2, 1].text(
    0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
)
ax[2, 1].text(
    0.5,
    0.2,
    "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
    **text_common_args,
)
ax[2, 1].text(
    0.5, 0.1, "Physical resistivity $\eta$: $%.2f$ " % (eta), **text_common_args
)
ax[2, 1].axis("off")
'''

for axi in ax:
    for axii in axi:
        axii.set_xticks([])
        axii.set_yticks([])
        axii.set_aspect("equal")

# Set appropriate colourbars
fig.colorbar(
    a00,
    ax=ax[0, 0],
    label=r"$\rho / \rho_0$",
    fraction=0.042,
    pad=0.04,
    location="left",
    ticks=np.linspace(0.0, 1.5, 7),
)

fig.colorbar(
    a01,
    ax=ax[0, 1],
    label=r"$\mathrm{log}_{10} \left( P / P_0 \right)$",
    fraction=0.042,
    pad=0.04,
    ticks=np.linspace(0.0, 3.0, 7),
)

'''
fig.colorbar(
    a10,
    ax=ax[1, 0],
    label=r"$|\mathbf{v}| / c_s$",
    fraction=0.042,
    pad=0.04,
    location="left",
    ticks=np.linspace(0.0, 5.0, 6),
)
'''

fig.colorbar(
    a10,
    ax=ax[1, 0],
    label=r"$|\mathbf{B}| / B_0$",
    fraction=0.042,
    pad=0.04,
    location='left',
    ticks=np.linspace(0.0, 2.0, 6),
)

fig.colorbar(
    a11,
    ax=ax[1, 1],
    extend="both",
    label=r"$\mathrm{log}_{10} \frac{h \: \nabla \cdot B}{|B|}$",
    fraction=0.042,
    pad=0.04,
    ticks=np.linspace(-5.0,0.0,6)
)

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig(args.input.replace('hdf5', 'png'))
