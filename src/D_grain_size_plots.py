# -*- coding: utf-8 -*-
"""
Script to create visualizations of soil samples for a survey to assess how well
people can estimate different characteristics of a soil sample.

Author: Georg H. Erharter (georg.erharter@ngi.no)
"""

# importing libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# importing custom libraries from file "X_library.py"
from X_library import laboratory, plotter


def add_grain(ax, x, y, r):
    circle = plt.Circle((x, y), r, facecolor='black', linewidth=0)
    ax.add_patch(circle)


###############################
# fixed values and constant variables
###############################

DENSITY = 2.65  # grain density [g/cm3]
TOT_MASS = 100  # [kg]
MIN_D, MAX_D = 1, 200  # [mm] min & max particle sizes of simulation
N_MESH_SIZES = 50  # number of mesh sizes between MIN_D, MAX_D
BOX_SIZE = 500  # square box size [mm]
SEED = 200  # random seed for reproducibility
MAX_GRAINS = 30_000  # maximum number of grains to place on visualization

###############################
# data loading, instantiations and some data processing
###############################

print(f'running seed: {SEED}')
np.random.seed(SEED)  # fix seed for reproducibility

lab, pltr = laboratory(), plotter()

# create virtual mesh sizes
sieve_sizes = np.exp(np.linspace(np.log(MIN_D), np.log(MAX_D), N_MESH_SIZES))

# generate soil to draw grains from
grain_diameters, grain_masses, grain_ids = lab.make_grains(
    DENSITY, TOT_MASS=TOT_MASS, min_d=MIN_D, max_d=MAX_D)


# draw grains for the plot -> random placement of grains changes overall soil
# distribution as large grains can only be placed on the plot in the beginning
placed_grain_xs = []
placed_grain_ys = []
placed_grain_rs = []

for i in range(len(grain_diameters[:MAX_GRAINS])):

    overlapping = True  # flag indicating if a grain overlaps with others
    count_fail = 0
    while overlapping is True:
        # discard grain if placement attempt failed 30 times
        if count_fail > 30:
            break
        else:
            # generate new grain coordinates
            x, y = np.random.uniform(0, BOX_SIZE, 2)
            r = grain_diameters[i] / 2  # [mm]
            # check if overlapping with other grains
            dx = x - np.array(placed_grain_xs)
            dy = y - np.array(placed_grain_ys)
            distances = np.sqrt(dx**2 + dy**2) - np.array(placed_grain_rs)
            # only place grain if it doesn't overlap
            if len(distances) == 0 or distances.min() >= r:
                overlapping = False
                placed_grain_xs.append(x)
                placed_grain_ys.append(y)
                placed_grain_rs.append(r)
            else:  # increased failed placement attempt counter
                count_fail += 1

print(f'{len(placed_grain_xs)} of {len(grain_diameters)} included')

# get stats of drawn grains that are actually placed on the visualization
placed_grain_ds = np.array(placed_grain_rs)*2
volumes = (4/3)*np.pi*(placed_grain_ds/2)**3  # [mm3]
volumes = volumes / 1000  # [cm3]
masses = volumes * DENSITY / 1000  # [kg]

d_max, d_min = placed_grain_ds.max(), placed_grain_ds.min()
# make sieve analysis
fractions_true = lab.sieve(placed_grain_ds, masses, sieve_sizes)
# calculate geometrical properites of the soil distribution
ds, Cu, Cc, S0 = lab.calc_grading_characteristics(
    placed_grain_ds, masses)


# fig, ax = plt.subplots()
# ax.plot(sieve_sizes, list(fractions_true.values()), label='sieve curve')
# id_sorted = np.argsort(placed_grain_ds)
# ax.plot(placed_grain_ds[id_sorted],
#         np.cumsum(masses[id_sorted])/sum(masses)*100, label='cumulative density function')
# # ax.set_xscale('log')
# # ax.set_xlim(left=100)
# ax.legend()
# ax.grid(alpha=0.5)
# ax.set_xlabel('grain size [mm]')
# ax.set_ylabel('mass percent passing [%]')
# plt.tight_layout()

# save statistics of sample as ground truth
df = pd.DataFrame(columns=['d_min', 'd_max', 'Cu', 'Cc', 'S0']+list(ds.keys()))
for d in ds.keys():
    df.loc[0, d] = ds[d]
df.loc[0, ['Cu', 'Cc', 'S0', 'd_min', 'd_max']] = [Cu, Cc, S0, d_min, d_max]
df.to_excel(fr'../samples/analyses/{SEED}_stats.xlsx')

###############################
# plotting of samples
###############################

# sieve curve of sample
pltr.sieve_curves_plot(sieve_sizes, list(fractions_true.values()),
                       close=True,
                       savepath=fr'../samples/analyses/{SEED}_sieve_curve.jpg')

# visualization of sample
fig = plt.figure(figsize=(9, 10))

ax1 = fig.add_axes([0.10, 0.20, 0.8, 0.8])  # [left, bottom, width, height]
for i in range(len(placed_grain_xs)):
    add_grain(ax1, placed_grain_xs[i], placed_grain_ys[i], placed_grain_rs[i])

# diameter_points = placed_grain_ds * 2.83465
# sizes = (diameter_points / 2) ** 2# * np.pi  # / 4.1
# ax.scatter(placed_grain_xs, placed_grain_ys,
#            s=sizes, alpha=0.5, linewidth=0)
ax1.set_aspect('equal')
ax1.set_xlim(0, BOX_SIZE)
ax1.set_ylim(top=0, bottom=BOX_SIZE)
ax1.set_xticks([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500])
ax1.set_yticks([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500])
ax1.grid(alpha=0.3)
ax1.set_xlabel('mm')
ax1.set_ylabel('mm')

# add second axis with reference grains
ax2 = fig.add_axes([0.1, 0.01, 0.8, 0.2])  # [left, bottom, width, height]

add_grain(ax2, 50, 50, 50)
ax2.text(90, 1, s='Ø = 100 mm', va='bottom')
add_grain(ax2, 160, 75, 25)
ax2.text(178, 55, s='Ø = 50 mm', va='top', ha='left', rotation=-40)
add_grain(ax2, 235, 80, 20)
ax2.text(248, 65, s='Ø = 40 mm', va='top', ha='left', rotation=-40)
add_grain(ax2, 300, 85, 15)
ax2.text(310, 75, s='Ø = 30 mm', va='top', ha='left', rotation=-40)
add_grain(ax2, 350, 90, 10)
ax2.text(356, 84, s='Ø = 20 mm', va='top', ha='left', rotation=-40)
add_grain(ax2, 400, 95, 5)
ax2.text(402, 93, s='Ø = 10 mm', va='top', ha='left', rotation=-40)
add_grain(ax2, 450, 99, 1)
ax2.text(451, 98, s='Ø = 2 mm', va='top', ha='left', rotation=-40)
ax2.set_aspect('equal')
ax2.set_xlim(0, BOX_SIZE)
ax2.set_ylim(0, 100)
ax2.axis('off')

plt.savefig(fr'../samples/analyses/{SEED}_sample.jpg', dpi=600)
plt.close()
