# -*- coding: utf-8 -*-
"""
Script to make individual virtual grain size analyses with different sample
masses and different underlzing distributions.

Author: Georg H. Erharter (georg.erharter@ngi.no)
"""

# importing libraries
import numpy as np
# importing custom libraries from file "X_library.py"
from X_library import laboratory, plotter, statistics


###############################
# fixed values and constant variables
###############################

# constant values and hyperparameters
DENSITY = 2.65  # grain density [g/cm3]
MIN_D, MAX_D = 1, 200  # [mm] min & max particle sizes of simulation
N_MESH_SIZES = 50  # number of mesh sizes between MIN_D, MAX_D
# fractions of the ISO required sample mass that are analyzed
FRACTIONS = [1, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01]
SEED = 6  # random seed for reproducibility

###############################
# main code execution
###############################
print(f'processing seed {SEED}')
# instantiations
lab, pltr, stat = laboratory(), plotter(), statistics()

np.random.seed(SEED)  # fix seed for reproducibility

# create virtual mesh sizes
sieve_sizes = np.exp(np.linspace(np.log(MIN_D), np.log(MAX_D), N_MESH_SIZES))

# initialize ground truth sample
grain_diameters, grain_masses, grain_ids = lab.make_grains(
    DENSITY, TOT_MASS=410, min_d=MIN_D, max_d=MAX_D)

total_mass = grain_masses.sum()

# calculate required sample size
standard_sample_mass = lab.ISO_required_sample_mass(grain_diameters.max()) / 1000
print(f'total mass: {total_mass}')
print(f'required mass: {round(standard_sample_mass, 3)} kg')
# break code if more material is required than available
if standard_sample_mass > total_mass:
    raise ValueError('required sample mass larger than total mass')

fractions_true = lab.sieve(grain_diameters, grain_masses, sieve_sizes)
# calculate grading characteristics of soil distribution
ds, Cu, Cc, S0 = lab.calc_grading_characteristics(grain_diameters,
                                                  grain_masses)

req_sample_masses = []
ks_distances = []
sieved_samples = []
sample_diameters = []

# do analyses for different fractions of required sample mass
for sample_fraction in FRACTIONS:
    print(f'{sample_fraction*100} % fraction in progress')
    req_sample_mass = standard_sample_mass*sample_fraction
    sample_ids, sample_masses, sample_diameter = lab.get_sample(
        req_sample_mass, total_mass, grain_masses, grain_diameters)
    sieved_sample = lab.sieve(sample_diameter, sample_masses, sieve_sizes)
    # calculate distribution - distribution distance metrics
    ks = stat.ks_statistic(list(fractions_true.values()),
                           list(sieved_sample.values()))
    # collect results
    req_sample_masses.append(req_sample_mass)
    sample_diameters.append(sample_diameter)
    ks_distances.append(ks)
    sieved_samples.append(list(sieved_sample.values()))

###############################
# result plotting
###############################

pltr.distances_plot(req_sample_masses, ks_distances,
                    savepath=fr'../figures/{SEED}_distances.jpg')

pltr.sieve_curves_plot(SIEVE_SIZES=sieve_sizes,
                       fractions_true=list(fractions_true.values()),
                       savepath=fr'../figures/{SEED}_sieve_line.jpg',
                       sieved_samples=sieved_samples,
                       req_sample_masses=req_sample_masses,
                       ks_distances=ks_distances)
