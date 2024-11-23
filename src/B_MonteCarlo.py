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

# number of grains of soil distribution -> impacts computation
DENSITY = 2.65  # grain density [g/cm3]
MIN_D, MAX_D = 1, 200  # [mm] min & max particle sizes of simulation
N_MESH_SIZES = 50  # number of mesh sizes between MIN_D, MAX_D
N_SIMULATIONS = 25  # number of new simulations to do
STUDY_NAME = '2024_11_22'  # study to work with or to create
TOT_MASS = 1800  # [kg]

###############################
# main code execution
###############################

print(STUDY_NAME)
# instantiations
lab, stat, pltr, utils = laboratory(), statistics(), plotter(), utilities()

# create virtual mesh sizes
sieve_sizes = np.exp(np.linspace(np.log(MIN_D), np.log(MAX_D), N_MESH_SIZES))

exponents = np.arange(1, 2.6, step=0.1)
standards = ['ISO', 'ASTM', 'const'] + [f'new {round(e, 1)}' for e in exponents]
mass_labels = [f'{s} req. mass [kg]' for s in standards]
ks_labels = [f'{s} ks [%]' for s in standards]

# empty lists to collect results of simulations
d_s, Cu_s, Cc_s, S0_s, soil_classes = [], [], [], [], []
min_diameters, max_diameters, total_masses = [], [], []
req_mass_p95s, fractions_trues = [], []
req_masses, ks_s = [], []

# main simulation loop
for i in tqdm(range(N_SIMULATIONS)):
    # generate underlying soil distribution -> so far only lognormal
    grain_diameters, grain_masses, grain_ids = lab.make_grains(
        DENSITY, TOT_MASS, min_d=MIN_D, max_d=MAX_D)
    # get maximum grain diameter of soil [mm] and total mass [kg]
    max_diameter = grain_diameters.max()
    min_diameter = grain_diameters.min()
    total_mass = grain_masses.sum()
    # make sieve analysis of underlying soil distribution
    fractions_true = lab.sieve(grain_diameters, grain_masses, sieve_sizes)
    # calculate geometrical properites of underlying soil distribution
    ds, Cu, Cc, S0 = lab.calc_grading_characteristics(grain_diameters,
                                                      grain_masses)
    # classify underlying soil distribution acc to USCS
    soil_class = lab.USCS_classification(ds['d12'], ds['d50'], Cu, Cc)

    # define sample mass that needs to be taken [kg]
    sample_mass_ISO = lab.ISO_required_sample_mass(max_diameter) / 1000
    sample_mass_ASTM = lab.ASTM_required_sample_mass(max_diameter) / 1000
    sample_mass_const = 10
    new_sample_masses = [lab.new_sample_mass(max_diameter, ds['d90'], exponent=exp) / 1000 for exp in exponents]

    # determine the requ. sample mass to achieve p95 ks of <= 10
    req_mass_p95 = lab.check_required_mass(
        max_error=10, mass_steps=10, tests_per_step=20,
        total_mass=total_mass, grain_masses=grain_masses,
        grain_diameters=grain_diameters, SIEVE_SIZES=sieve_sizes,
        fractions_true=fractions_true, verbose=False)

    req_sample_masses = [sample_mass_ISO, sample_mass_ASTM,
                         sample_mass_const] + new_sample_masses

    # only process if the reqired sample masses are lower than the total mass.
    # Problems may occur when there is a large amount of finer grained
    # sediments
    if utils.CheckForLess(req_sample_masses, total_mass):
        # compute errors for all sample masses
        ks_s_temp = []
        for i, req_sample_mass in enumerate(req_sample_masses):
            # get sample out of underlying distribution acc to sample_mass
            sample_ids, sample_masses, sample_diameter = lab.get_sample(
                req_sample_mass, total_mass, grain_masses, grain_diameters)
            sieved_sample = lab.sieve(sample_diameter, sample_masses,
                                      sieve_sizes)
            # compute kolmogorov smirnov distance between sample and real soil
            ks_s_temp.append(stat.ks_statistic(list(fractions_true.values()),
                                               list(sieved_sample.values())))

        # collect data of simulation
        d_s.append(ds)
        Cu_s.append(Cu)
        Cc_s.append(Cc)
        S0_s.append(S0)
        min_diameters.append(min_diameter)
        max_diameters.append(max_diameter)
        total_masses.append(total_mass)
        soil_classes.append(soil_class)
        req_mass_p95s.append(req_mass_p95)
        fractions_trues.append(fractions_true)
        req_masses.append(dict(zip(mass_labels, req_sample_masses)))
        ks_s.append(dict(zip(ks_labels, ks_s_temp)))

    else:
        raise ValueError('required sample mass too large')

# make new pandas dataframe with new simulation results
df_new = pd.DataFrame({'Cu': Cu_s,
                       'Cc': Cc_s,
                       'S0': S0_s,
                       'USCS soil classes': soil_classes,
                       'min diameter [mm]': min_diameters,
                       'max diameter [mm]': max_diameters,
                       'total masses [kg]': total_masses,
                       'req. mass ks_p95 <= 10 [kg]': req_mass_p95s
                       })
df_fractions = pd.DataFrame.from_dict(fractions_trues)
df_ds = pd.DataFrame.from_dict(d_s)
df_req_masses = pd.DataFrame.from_dict(req_masses)
df_ks = pd.DataFrame.from_dict(ks_s)
df_new = pd.concat((df_new, df_fractions, df_ds, df_req_masses, df_ks),
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
df['ID'] = np.arange(len(df))
df.to_excel(fr'../simulations/{STUDY_NAME}.xlsx', index=False)
print(f'{len(df)} samples processed')
