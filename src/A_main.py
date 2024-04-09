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
sieve_sizes = np.exp(np.linspace(np.log(MIN_D), np.log(MAX_D), 30))

# initialize ground truth sample
grain_diameters, grain_weights, grain_ids = lab.make_grains(
    DENSITY, TOT_MASS=410, min_d=MIN_D, max_d=MAX_D)

total_weight = grain_weights.sum()

# calculate required sample size
standard_sample_weight = lab.ISO_required_sample_weight(grain_diameters.max()) / 1000
print(f'total weight: {total_weight}')
print(f'required weight: {round(standard_sample_weight, 3)} kg')
# break code if more material is required than available
if standard_sample_weight > total_weight:
    raise ValueError('required sample weight larger than total weight')

fractions_true = lab.sieve(grain_ids, grain_diameters, grain_weights,
                           sieve_sizes)
# calculate grading characteristics of soil distribution
ds, Cu, Cc, S0 = lab.calc_grading_characteristics(
    list(fractions_true.values()), sieve_sizes)

req_sample_weights = []
ks_distances = []
sieved_samples = []
sample_diameters = []

# do analyses for different fractions of required sample mass
for sample_fraction in FRACTIONS:
    print(f'{sample_fraction*100} % fraction in progress')
    req_sample_weight = standard_sample_weight*sample_fraction
    sample_ids, sample_weights, sample_diameter = lab.get_sample(
        req_sample_weight, total_weight, grain_weights, grain_diameters)
    sieved_sample = lab.sieve(sample_ids, sample_diameter, sample_weights,
                              sieve_sizes)
    # calculate distribution - distribution distance metrics
    ks = stat.ks_statistic(list(fractions_true.values()),
                           list(sieved_sample.values()))
    # collect results
    req_sample_weights.append(req_sample_weight)
    sample_diameters.append(sample_diameter)
    ks_distances.append(ks)
    sieved_samples.append(list(sieved_sample.values()))

###############################
# result plotting
###############################

pltr.distances_plot(req_sample_weights, ks_distances,
                    savepath=fr'../figures/{SEED}_distances.jpg')

pltr.sieve_curves_plot(
    sieve_sizes, list(fractions_true.values()),
    fr'../figures/{SEED}_sieve_line.jpg',
    sieved_samples, req_sample_weights, ks_distances)

# # make sample preview plot ... very experimental still
# weight = 1.5
# plot_sample = lab.get_sample(weight, total_weight, grain_weights,
#                              grain_diameters)[2]
# if len(plot_sample) < 30_000:
#     pltr.plot_grains(plot_sample, 0.5, SEED, Cu, S0, weight,
#                      savepath=fr'../figures/{SEED}_sample.jpg')
# else:
#     print('too many grains', len(plot_sample))
