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
from X_library import laboratory, statistics, plotter, utilities


###############################
# fixed values and constant variables
###############################

# number of grains of underlying soil distribution -> impacts computation
DENSITY = 2.65  # grain density [g/cm3]
SIEVE_SIZES = [0.075, 0.105, 0.15, 0.25, 0.425, 0.85, 2, 4.75, 9.5, 19, 25, 37.5, 50, 75, 100, 150, 200, 300]  # sieve sizes [mm]
# weight of soil sample, either set to number for sampling always the same mass
# in kg or set to "ISO" to sample the mass as suggested by ISO 17892-4
N_SIMULATIONS = 2  # number of new simulations to do
STUDY_NAME = '2024_02_23'  # study to work with or to create
PLOT = True  # flag to indicate if plots shall be created
TOT_MASS = 1256.7  # [kg]

###############################
# main code execution
###############################

# instantiations
lab, stat, pltr, utils = laboratory(), statistics(), plotter(), utilities()

# empty lists to collect results of simulations
d_s, Cu_s, Cc_s, S0_s, soil_classes = [], [], [], [], []
max_diameters, total_weights = [], []
ks_ISO, ks_ASTM, ks_new, ks_const = [], [], [], []
req_weights_ISO, req_weights_ASTM = [], []
req_weights_new, req_weights_const = [], []

# main simulation loop
for i in tqdm(range(N_SIMULATIONS)):
    # generate underlying soil distribution -> so far only lognormal
    grain_diameters, grain_weights, grain_ids = lab.make_grains_new(
        DENSITY, TOT_MASS)
    # get maximum grain diameter of soil [mm] and total weight [kg]
    max_diameter = grain_diameters.max()
    total_weight = grain_weights.sum()
    # make sieve analysis of underlying soil distribution
    fractions_true = lab.sieve(grain_ids, grain_diameters, grain_weights,
                               SIEVE_SIZES)
    # calculate geometrical properites of underlying soil distribution
    ds, Cu, Cc, S0 = lab.calc_grading_characteristics(
        fractions_true, SIEVE_SIZES)
    # classify underlying soil distribution acc to USCS
    soil_class = lab.USCS_classification(ds['d12'], ds['d50'], Cu, Cc)

    # define sample weight [kg]
    sample_weight_ISO = lab.ISO_required_sample_weight(max_diameter) / 1000
    sample_weight_ASTM = lab.ASTM_required_sample_weight(max_diameter) / 1000
    sample_weight_new = lab.new_sample_weight(S0)
    sample_weight_const = 10
    req_sample_weights = [sample_weight_ISO, sample_weight_ASTM,
                          sample_weight_new, sample_weight_const]

    # only process if the reqired sample weights are lower than the total
    # weight. problem may occur when there is a large amount of finer grained
    # sediments
    if utils.CheckForLess(req_sample_weights, total_weight):
        # collect parameters of underlying soil distribution
        d_s.append(ds)
        Cu_s.append(Cu)
        Cc_s.append(Cc)
        S0_s.append(S0)
        max_diameters.append(max_diameter)
        total_weights.append(total_weight)
        soil_classes.append(soil_class)

        for i, req_sample_weight in enumerate(req_sample_weights):
            # get sample out of underlying distribution acc to sample_weight
            sample_ids, sample_weights, sample_diameter = lab.get_sample(
                req_sample_weight, total_weight, grain_weights,
                grain_diameters)
            sieved_sample = lab.sieve(sample_ids, sample_diameter,
                                      sample_weights, SIEVE_SIZES)
            # compute kolmogorov smirnov distance between sample and real soil
            ks = stat.ks_statistic(fractions_true, sieved_sample)
            # collect results of sampled soils
            match i:  # noqa
                case 0:
                    ks_ISO.append(ks)
                    req_weights_ISO.append(req_sample_weight)
                case 1:
                    ks_ASTM.append(ks)
                    req_weights_ASTM.append(req_sample_weight)
                case 2:
                    ks_new.append(ks)
                    req_weights_new.append(req_sample_weight)
                case 3:
                    ks_const.append(ks)
                    req_weights_const.append(req_sample_weight)
    else:
        raise ValueError('required sample weight too large')

# make new pandas dataframe with new simulation results
df_new = pd.DataFrame({'Cu': Cu_s,
                       'Cc': Cc_s,
                       'S0': S0_s,
                       'USCS soil classes': soil_classes,
                       'max diameter [mm]': max_diameters,
                       'total weights [kg]': total_weights,
                       'kolmogorov smirnov distance ISO': ks_ISO,
                       'kolmogorov smirnov distance ASTM': ks_ASTM,
                       'kolmogorov smirnov distance new': ks_new,
                       'kolmogorov smirnov distance const': ks_const,
                       'req. weight ISO [kg]': req_weights_ISO,
                       'req. weight ASTM [kg]': req_weights_ASTM,
                       'req. weight new [kg]': req_weights_new,
                       'req. weight const [kg]': req_weights_const,
                       })
df_ds = pd.DataFrame.from_dict(d_s)
df_new = pd.concat((df_new, df_ds), axis=1)

# either make new dataframe for new study or add to existing one & save
try:
    df = pd.read_excel(fr'../simulations/{STUDY_NAME}.xlsx')
    df = pd.concat([df, df_new])
except FileNotFoundError:
    df = df_new
df.to_excel(fr'../simulations/{STUDY_NAME}.xlsx', index=False)

##########################
# plot results
##########################

if PLOT is True:
    weight_modes = ['ISO', 'ASTM', 'new', 'const']

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    for weight_mode in weight_modes:
        ax.scatter(df['S0'],
                   df[f'kolmogorov smirnov distance {weight_mode}'],
                   df[f'req. weight {weight_mode} [kg]'], label=weight_mode)

    ax.set_xlabel('sorting coefficient\nS0 sqrt(d75/d25)')
    ax.set_ylabel('kolmogorov smirnov distance')
    ax.set_zlabel('req. weight [kg]')
    ax.legend()

    plt.tight_layout()

    for weight_mode in weight_modes:
        pltr.monte_carlo_scatterplot(
            df, weight_mode=weight_mode, color_mode='weight',
            savepath=fr'../simulations/{STUDY_NAME}_{weight_mode}_weight.jpg')
        pltr.monte_carlo_scatterplot(
            df, weight_mode=weight_mode, color_mode='USCS',
            savepath=fr'../simulations/{STUDY_NAME}_{weight_mode}_USCS.jpg')
