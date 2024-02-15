# -*- coding: utf-8 -*-
"""
Script to make individual virtual grain size analyses with different sample
masses and different underlzing distributions.

Author: Georg H. Erharter (georg.erharter@ngi.no)
"""

# importing libraries
import numpy as np
from scipy.stats import wasserstein_distance, energy_distance
# importing custom libraries from file "X_library.py"
from X_library import laboratory, plotter, statistics


###############################
# fixed values and constant variables
###############################

# constant values and hyperparameters
N = 20_000_000  # number of grains of underlying soil distribution
DENSITY = 2.5  # grain density [g/cm3]
SIEVE_SIZES = [0.04, 0.1, 0.25, 0.425, 0.85, 2, 4.75, 10, 20, 50, 100, 150, 200, 300]  # sieve sizes [mm]
# underlying soil distribution shape
DISTRIBUTION = 'lognormal'  #  normal, exponential, beta, uniform, lognormal, combined
# fractions of the ISO required sample mass that are analyzed
FRACTIONS = [1, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01]
SEED = 0  # random seed for reproducibility

###############################
# main code execution
###############################

# instantiations
lab, pltr, stat = laboratory(), plotter(), statistics()

np.random.seed(SEED)  # fix seed for reproducibility

# initialize ground truth sample
grain_diameters, grain_volumes, grain_weights, grain_ids = lab.make_grains(
    N, DENSITY, DISTRIBUTION)

# calculate required sample size
standard_sample_weight = lab.ISO_required_sample_weight(grain_diameters) / 1000
print(DISTRIBUTION)
print(f'total weight: {sum(grain_weights)}')
print(f'required weight: {round(standard_sample_weight, 3)} kg')
# break code if more material is required than available
if standard_sample_weight > sum(grain_weights):
    raise ValueError('required sample weight larger than total weight')

fractions_true = lab.sieve(grain_ids, grain_diameters, grain_weights,
                           SIEVE_SIZES)
# calculate grading characteristics of soil distribution
d10, d12, d30, d50, d60, Cu, Cc, S0 = lab.calc_grading_characteristics(
    fractions_true, SIEVE_SIZES)

req_sample_weights = []
wasserstein_distances = []
energy_distances = []
ks_distances = []
sieved_samples = []
sample_diameters = []

# do analyses for different fractions of required sample mass
for sample_fraction in FRACTIONS:
    print(f'{sample_fraction*100} % fraction in progress')
    req_sample_weight = standard_sample_weight*sample_fraction
    sample_ids, sample_weights, sample_diameter = lab.get_sample(
        req_sample_weight, grain_weights, grain_diameters)
    # calculate distribution - distribution distance metrics
    wd = wasserstein_distance(grain_diameters, sample_diameter)
    ed = energy_distance(grain_diameters, sample_diameter)
    ks = stat.ks_statistic(grain_diameters, sample_diameter)
    sieved_sample = lab.sieve(sample_ids, sample_diameter, sample_weights,
                              SIEVE_SIZES)
    # collect results
    req_sample_weights.append(req_sample_weight)
    sample_diameters.append(sample_diameter)
    wasserstein_distances.append(wd)
    energy_distances.append(ed)
    ks_distances.append(ks)
    sieved_samples.append(sieved_sample)

###############################
# result plotting
###############################

pltr.distances_plot(
    grain_diameters, sample_diameters, req_sample_weights, ks_distances,
    savepath=fr'../figures/{DISTRIBUTION}_{SEED}_distances.jpg')

pltr.sieve_curves_plot(
    SIEVE_SIZES, fractions_true, sieved_samples, req_sample_weights,
    ks_distances, d60, d30, d10, Cu, Cc, grain_diameters,
    standard_sample_weight, DISTRIBUTION,
    savepath=fr'../figures/{DISTRIBUTION}_{SEED}_sample.jpg')

# plot to visualize theoretically required sample mass
pltr.required_weight_plot(300, '../figures/required_weight.svg')
