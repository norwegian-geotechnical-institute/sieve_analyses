# -*- coding: utf-8 -*-
"""
Script to primarily create visualizations as a means to build equations of
relationships between parameters.

Author: Georg H. Erharter (georg.erharter@ngi.no)
"""

import numpy as np
import pandas as pd

from X_library import utilities, laboratory, plotter


###############################
# fixed values and constant variables
###############################

STUDY_NAME = '2024_07_07'  # study to work with or to create

###############################
# data loading, instantiations and some data processing
###############################

df = pd.read_excel(fr'../simulations/{STUDY_NAME}.xlsx')
print(len(df))
lab, pltr, utils = laboratory(), plotter(), utilities()

# apply new functions to compute sample mass recommendations with ISO and
# ASTM p95 errors: ISO 7.7, ASTM 5.3
df['new. req. mass ISO [kg]'] = (df['d90']/10)**((np.log(7.7)-np.log(123.65))/-1.29)
df['new. req. mass ASTM [kg]'] = (df['d90']/10)**((np.log(5.3)-np.log(123.65))/-1.29)
# ratio
df['comparison ISO'] = df['ISO req. mass [kg]'] / df['new. req. mass ISO [kg]']
df['comparison ASTM'] = df['ASTM req. mass [kg]'] / df['new. req. mass ASTM [kg]']
print(f"ISO mass = {round(df['comparison ISO'].mean(), 2)} larger on avg")
print(f"ISO mass = {round(df['comparison ISO'].max(), 2)} larger at max")
print(f"ASTM mass = {round(df['comparison ASTM'].mean(), 2)} larger on avg")
print(f"ASTM mass = {round(df['comparison ASTM'].max(), 2)} larger at max")
# move ID column to the front
ID_column = df.pop('ID')
df.insert(0, 'ID', ID_column)
df.to_excel(r'../simulations/MonteCarloSimulations.xlsx', index=False)

###############################
# plotting
###############################

# plot to visualize theoretically required sample mass
pltr.required_mass_plot(300, '../figures/required_mass.svg')

# plot showing errors of ASTM and ISO
pltr.error_violin_plot(df, r'../figures/error_violin.jpg')

# scatterplot showing required sample mass against max diameter
pltr.req_sample_mass_vs_dmax_plot(df, annotate_all=False,
                                  annotate_some=[310, 210, 94, 52],
                                  close=True,
                                  savepath=r'../figures/req_mass_dmax.jpg')

# plot individual simplified sieve curves for combined plots
# 310 ... max req. mass, 210 ... max dmax & max req. mass, 94... min dmax
# 52 ... center
pltr.simple_sieve_plot(df, ids=[310, 210, 94, 52], close=True,
                       savepath=r'../figures/req_mass_dmax_samples.jpg')

# plot showing required sample mass against d90
pltr.req_sample_mass_vs_d90_plot(df, annotate_some=[310, 94, 210, 52],
                                 savepath=r'../figures/req_mass_d90.jpg')

# plot showing new m_min functions with different epsilon
# and plot showing relationship between exponent and errors
pltr.exponents_plot(df, savepath=r'../figures/exponents.jpg')

# plot showing a comparison between standard required sample masses and new one
pltr.comparison_plot(df, savepath=r'../figures/comparison.jpg')

# plot showing exemplary sieve curves
sieve_cols = [c for c in df.columns if 'mm sieve [m%]' in c]
sieve_sizes = [eval(s.split(' ')[0]) for s in sieve_cols]
fractions_trues = [list(df[sieve_cols].iloc[i].values) for i in range(100)]

pltr.sieve_curves_plot(SIEVE_SIZES=sieve_sizes, fractions_true=fractions_trues,
                       color=df['S0'],
                       savepath=r'../figures/sieve_samples_new.jpg',
                       close=True)

# plot showing real lab results: sieve curves
pltr.real_sieve_curves_plot(savepath=r'../laboratory/lab_tests_PSDs.jpg')

# plot showing real lab results: scatter plot
pltr.real_sieve_curves_scatter(r'../laboratory/lab_tests_CcCu_scatter.jpg')

###############################
# other analyses of the simulations
###############################

# check correlation between parameters and "bottom up" determined required
# sample mass
for f in ['Cu', 'Cc', 'S0', 'max diameter [mm]', 'd10', 'd12', 'd25',
          'd30', 'd50', 'd60', 'd75', 'd90']:
    cc = np.corrcoef(df[f], df['req. mass ks_p95 <= 10 [kg]'])[0][1]
    print(f'{round(cc, 2)} = corr. coeff. {f} - required mass')

# # new parameter for grading characterization
# df['grading'] = np.mean(np.vstack(((np.log(df['d30']) - np.log(df['d10'])),
#                                   (np.log(df['d50']) - np.log(df['d30'])),
#                                   (np.log(df['d70']) - np.log(df['d50'])),
#                                   (np.log(df['d90']) - np.log(df['d70'])))).T,
#                        axis=1)
# fig, ax = plt.subplots()
# ax.scatter(df['Cc'], df['grading'])
# df['grading'] = np.where(df['grading'] > 2, 1, 0)
