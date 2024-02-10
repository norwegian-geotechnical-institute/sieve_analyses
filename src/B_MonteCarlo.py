# -*- coding: utf-8 -*-
"""
Script to make simulated sieve analyses for Monte Carlo analyses to determine
correlations between distribution geometries and sample masses.

Author: Georg H. Erharter (georg.erharter@ngi.no)
"""

# importing libraries
import numpy as np
import pandas as pd
from tqdm import tqdm
# importing custom libraries from file "X_library.py"
from X_library import laboratory, statistics, plotter


###############################
# fixed values and constant variables
###############################

# number of grains of underlying soil distribution -> impacts computation
N = 10_000_000
DENSITY = 2.5  # grain density [g/cm3]
SIEVE_SIZES = [0.04, 0.1, 0.25, 0.425, 0.85, 2, 4.75, 10, 20, 50, 100, 150, 200, 300]  # sieve sizes [mm]
# underlying soil distribution shape
DISTRIBUTION = 'lognormal'  #  normal, exponential, beta, uniform, lognormal, combined
# weight of soil sample, either set to number for sampling always the same mass
# in kg or set to "ISO" to sample the mass as suggested by ISO 17892-4
SAMPLE_WEIGHT = 'ISO'
N_SIMULATIONS = 2  # number of new simulations to do
SAVEPATH = r'../simulations/data.xlsx'  # store simulation results

###############################
# main code execution
###############################

# instantiations
lab, stat, pltr = laboratory(), statistics(), plotter()

# empty lists to collect results of simulations
d10_s, Cu_s, Cc_s, S0_s = [], [], [], []
ks_s, means, sigmas, req_weights, max_diameters = [], [], [], [], []
sample_weights_, soil_classes = [], []

# main simulation loop
for i in tqdm(range(N_SIMULATIONS)):
    # generate underlying soil distribution -> so far only lognormal
    mean = np.random.uniform(0, 3)
    sigma = np.random.uniform(0.1, 2)
    grain_diameters, grain_volumes, grain_weights, grain_ids = lab.make_grains(
        N, DENSITY, DISTRIBUTION, lognorm_mean=0, lognorm_sigma=sigma,
        verbose=False)

    # define sample weight
    match SAMPLE_WEIGHT:  # noqa
        case 'ISO':
            sample_weight = lab.calc_required_sample_weight(grain_diameters) / 1000
        case _:
            sample_weight = SAMPLE_WEIGHT  # kg
    # only process if the sampled weight is lower than the total weight
    # problem may occur when there is a large amount of fine grained sediments
    if sample_weight <= sum(grain_weights):
        max_diameter = max(grain_diameters)
        # make sieve analysis of underlying soil distribution
        fractions_true = lab.sieve(grain_ids, grain_diameters, grain_weights,
                                   SIEVE_SIZES)
        # calculate geometrical properites of underlying soil distribution
        d10, d12, d30, d50, d60, Cu, Cc, S0 = lab.calc_grading_characteristics(
            fractions_true, SIEVE_SIZES)
        # classify underlying soil distribution acc to USCS
        soil_class = lab.USCS_classification(d12, d50, Cu, Cc)

        # get sample out of underlying distribution acc to sample_weight
        sample_ids, sample_weights, sample_diameter = lab.get_sample(
            sample_weight, grain_weights, grain_diameters)
        # compute kolmogorov smirnov distance between sample and real soil
        ks = stat.ks_statistic(grain_diameters, sample_diameter)

        # collect results
        d10_s.append(d10)
        Cu_s.append(Cu)
        Cc_s.append(Cc)
        S0_s.append(S0)
        ks_s.append(ks)
        means.append(mean)
        sigmas.append(sigma)
        req_weights.append(sample_weight)
        max_diameters.append(max_diameter)
        sample_weights_.append(sample_weight)
        soil_classes.append(soil_class)

# make new pandas dataframe with new simulation results
df_new = pd.DataFrame({'lognorm_mean': means,
                       'lognorm_sigma': sigmas,
                       'max diameter true [mm]': max_diameters,
                       'req. weight [kg]': req_weights,
                       'sample weight [kg]': sample_weights_,
                       'd10 true [mm]': d10_s,
                       'Cu true': Cu_s,
                       'Cc true': Cc_s,
                       'S0 true': S0_s,
                       'kolmogorov smirnov distance': ks_s,
                       'USCS soil classes': soil_classes})

# either make new dataframe for new study or add to existing one & save
try:
    df = pd.read_excel(SAVEPATH)
    df = pd.concat([df, df_new])
except FileNotFoundError:
    df = df_new
df.to_excel(SAVEPATH, index=False)

##########################
# plot results
##########################

pltr.monte_carlo_scatterplot(df, savepath=fr'../simulations/monte_carlo.jpg')
