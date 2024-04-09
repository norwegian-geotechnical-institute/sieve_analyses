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
MIN_D, MAX_D = 1, 200  # [mm] min & max particle sizes of simulation
N_SIMULATIONS = 23  # number of new simulations to do
STUDY_NAME = '2024_03_19'  # study to work with or to create
PLOT = False  # flag to indicate if plots shall be created
TOT_MASS = 1800  # [kg]

###############################
# main code execution
###############################

print(STUDY_NAME)
# instantiations
lab, stat, pltr, utils = laboratory(), statistics(), plotter(), utilities()

# create virtual mesh sizes
sieve_sizes = np.exp(np.linspace(np.log(MIN_D), np.log(MAX_D), 30))

exponents = np.arange(1, 2.6, step=0.1)
standards = ['ISO', 'ASTM', 'const'] + [f'new {round(e, 1)}' for e in exponents]
weight_labels = [f'{s} req. weight [kg]' for s in standards]
ks_labels = [f'{s} ks [%]' for s in standards]

# empty lists to collect results of simulations
d_s, Cu_s, Cc_s, S0_s, soil_classes = [], [], [], [], []
max_diameters, total_weights = [], []
req_weight_p95s, fractions_trues = [], []
req_weights, ks_s = [], []

# main simulation loop
for i in tqdm(range(N_SIMULATIONS)):
    # generate underlying soil distribution -> so far only lognormal
    grain_diameters, grain_weights, grain_ids = lab.make_grains(
        DENSITY, TOT_MASS, min_d=MIN_D, max_d=MAX_D)
    # get maximum grain diameter of soil [mm] and total weight [kg]
    max_diameter = grain_diameters.max()
    total_weight = grain_weights.sum()
    # make sieve analysis of underlying soil distribution
    fractions_true = lab.sieve(grain_ids, grain_diameters, grain_weights,
                               sieve_sizes)
    # calculate geometrical properites of underlying soil distribution
    ds, Cu, Cc, S0 = lab.calc_grading_characteristics(
        list(fractions_true.values()), sieve_sizes)
    # classify underlying soil distribution acc to USCS
    soil_class = lab.USCS_classification(ds['d12'], ds['d50'], Cu, Cc)

    # define sample weight that needs to be taken [kg]
    sample_weight_ISO = lab.ISO_required_sample_weight(max_diameter) / 1000
    sample_weight_ASTM = lab.ASTM_required_sample_weight(max_diameter) / 1000
    sample_weight_const = 10
    new_sample_weights = [lab.new_sample_weight(max_diameter, ds['d90'], exponent=exp) / 1000 for exp in exponents]

    # determine the requ. sample weight to achieve p95 ks of <= 10
    req_weight_p95 = lab.check_required_weight(
        max_error=10, weight_steps=15, tests_per_step=12,
        total_weight=total_weight, grain_weights=grain_weights,
        grain_diameters=grain_diameters, SIEVE_SIZES=sieve_sizes,
        fractions_true=fractions_true, verbose=False)

    req_sample_weights = [sample_weight_ISO, sample_weight_ASTM,
                          sample_weight_const] + new_sample_weights

    # only process if the reqired sample weights are lower than the total
    # weight. problem may occur when there is a large amount of finer grained
    # sediments
    if utils.CheckForLess(req_sample_weights, total_weight):
        # compute errors for all sample weights
        ks_s_temp = []
        for i, req_sample_weight in enumerate(req_sample_weights):
            # get sample out of underlying distribution acc to sample_weight
            sample_ids, sample_weights, sample_diameter = lab.get_sample(
                req_sample_weight, total_weight, grain_weights,
                grain_diameters)
            sieved_sample = lab.sieve(sample_ids, sample_diameter,
                                      sample_weights, sieve_sizes)
            # compute kolmogorov smirnov distance between sample and real soil
            ks_s_temp.append(stat.ks_statistic(list(fractions_true.values()),
                                               list(sieved_sample.values())))

        # collect data of simulation
        d_s.append(ds)
        Cu_s.append(Cu)
        Cc_s.append(Cc)
        S0_s.append(S0)
        max_diameters.append(max_diameter)
        total_weights.append(total_weight)
        soil_classes.append(soil_class)
        req_weight_p95s.append(req_weight_p95)
        fractions_trues.append(fractions_true)
        req_weights.append(dict(zip(weight_labels, req_sample_weights)))
        ks_s.append(dict(zip(ks_labels, ks_s_temp)))

    else:
        raise ValueError('required sample weight too large')

# make new pandas dataframe with new simulation results
df_new = pd.DataFrame({'Cu': Cu_s,
                       'Cc': Cc_s,
                       'S0': S0_s,
                       'USCS soil classes': soil_classes,
                       'max diameter [mm]': max_diameters,
                       'total weights [kg]': total_weights,
                       'req. weight ks_p95 <= 10 [kg]': req_weight_p95s
                       })
df_fractions = pd.DataFrame.from_dict(fractions_trues)
df_ds = pd.DataFrame.from_dict(d_s)
df_req_weights = pd.DataFrame.from_dict(req_weights)
df_ks = pd.DataFrame.from_dict(ks_s)
df_new = pd.concat((df_new, df_fractions, df_ds, df_req_weights, df_ks),
                   axis=1)

df_new.rename(dict(zip(sieve_sizes,
                       [f'{s} mm sieve [m%]' for s in sieve_sizes])),
              axis=1, inplace=True)

# either make new dataframe for new study or add to existing one & save
try:
    df = pd.read_excel(fr'../simulations/{STUDY_NAME}.xlsx')
    df = pd.concat([df, df_new])
except FileNotFoundError:
    df = df_new
df.to_excel(fr'../simulations/{STUDY_NAME}.xlsx', index=False)
print(f'{len(df)} samples processed')

##########################
# plot results
##########################

if PLOT is True:
    weight_modes = ['ISO', 'ASTM', 'new', 'const']

    for weight_mode in weight_modes:
        pltr.monte_carlo_scatterplot(
            df, weight_mode=weight_mode, color_mode='weight',
            savepath=fr'../simulations/{STUDY_NAME}_{weight_mode}_weight.jpg')
        pltr.monte_carlo_scatterplot(
            df, weight_mode=weight_mode, color_mode='USCS',
            savepath=fr'../simulations/{STUDY_NAME}_{weight_mode}_USCS.jpg')
