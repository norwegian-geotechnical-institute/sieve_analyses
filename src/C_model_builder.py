# -*- coding: utf-8 -*-
"""
Script to primarily create visualizations as a means to build equations of
relationships between parameters.

Author: Georg H. Erharter (georg.erharter@ngi.no)
"""


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd

from X_library import utilities, laboratory, plotter


###############################
# fixed values and constant variables
###############################

STUDY_NAME = '2024_07_07'  # study to work with or to create

###############################
# data loading and instantiations
###############################

df = pd.read_excel(fr'../simulations/{STUDY_NAME}.xlsx')
print(len(df))
lab, pltr, utils = laboratory(), plotter(), utilities()

###############################
# plotting
###############################

# plot to visualize theoretically required sample mass
pltr.required_weight_plot(300, '../figures/required_weight.svg')

# plot showing errors of ASTM and ISO
pltr.error_violin_plot(df, r'../figures/error_violin.jpg')


# scatterplot showing required sample mass against max diameter
pltr.req_sample_mass_vs_dmax_plot(df, annotate_all=False,
                                  annotate_some=[310, 210, 94, 52],
                                  close=True,
                                  savepath=r'../figures/req_weight_dmax.jpg')

# plot individual simplified sieve curves for combined plots
# 310 ... max req. weight, 210 ... max dmax & max req. weight, 94... min dmax
# 52 ... center
pltr.simple_sieve_plot(df, ids=[310, 210, 94, 52], close=True,
                       savepath=r'../figures/req_weight_dmax_samples.jpg')

# plot showing required sample mass against d90
pltr.req_sample_mass_vs_d90_plot(df, annotate_some=[310, 94, 210, 52],
                                 savepath=r'../figures/req_weight_d90.jpg')

# plot showing new m_min functions with different epsilon
# and plot showing relationship between exponent and errors
pltr.exponents_plot(df, savepath=r'../figures/exponents.jpg')

print(ghjklø)

###############################
# plot showing new m_min functions with different epsilon
###############################

# ISO 7.6, ASTM 5.1

df['new. req. weight ISO [kg]'] = (df['d90']/10)**((np.log(8.2)-np.log(125.31))/-1.29)
df['new. req. weight ASTM [kg]'] = (df['d90']/10)**((np.log(5.3)-np.log(125.31))/-1.29)

df['ISO req. weight [kg]'] / df['new. req. weight ISO [kg]']
df['ASTM req. weight [kg]'] / df['new. req. weight ASTM [kg]']

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))

ax1.scatter(df['ISO req. weight [kg]'], df['new. req. weight ISO [kg]'],
            color='grey', edgecolor='black', alpha=0.4)
ax1.plot([0, 400], [0, 400], color='black', lw=4, ls='--')
ax1.text(.05, .94, 'a', fontsize=20, transform=ax1.transAxes)
ax1.grid(alpha=0.5)
ax1.set_xlabel('ISO req. weight [kg]')
ax1.set_ylabel('req. weight acc. new criterium [kg]')

ax2.scatter(df['ASTM req. weight [kg]'], df['new. req. weight ASTM [kg]'],
            color='grey', edgecolor='black', alpha=0.4)
ax2.plot([0, 1200], [0, 1200], color='black', lw=4, ls='--')
ax2.text(.05, .94, 'b', fontsize=20, transform=ax2.transAxes)
ax2.grid(alpha=0.5)
ax2.set_xlabel('ASTM req. weight [kg]')
ax2.set_ylabel('req. weight acc. new criterium [kg]')

plt.tight_layout()
plt.savefig(r'../figures/comparison.jpg')
plt.close()


sieve_cols = [c for c in df.columns if 'mm sieve [m%]' in c]

sieve_sizes = [eval(s.split(' ')[0]) for s in sieve_cols]

fractions_trues = [list(df[sieve_cols].iloc[i].values) for i in range(100)]

# new parameter for grading characterization
df['grading'] = np.mean(np.vstack(((np.log(df['d30']) - np.log(df['d10'])),
                                  (np.log(df['d50']) - np.log(df['d30'])),
                                  (np.log(df['d70']) - np.log(df['d50'])),
                                  (np.log(df['d90']) - np.log(df['d70'])))).T,
                       axis=1)
fig, ax = plt.subplots()
ax.scatter(df['Cc'], df['grading'])
df['grading'] = np.where(df['grading'] > 2, 1, 0)

pltr.sieve_curves_plot(SIEVE_SIZES=sieve_sizes, fractions_true=fractions_trues,
                       color=df['S0'],
                       savepath=r'../figures/sieve_samples_new.jpg',
                       close=True)

###############################
# plot showing real lab results: sieve curves
###############################

fp = r'../laboratory/LabResults.xlsx'
df = pd.read_excel(fp, header=2, nrows=13, usecols=list(range(1, 9)))
headers = list(df.columns)

fig, ax = pltr.make_sieve_plot()

ax.plot(df['Sieve size Ø [mm]'], df['Soil C (ISO)'],
        lw=3, color='C1', label='Soil C (ISO), 50 kg')
ax.plot(df['Sieve size Ø [mm]'], df['Soil C'],
        lw=1.5, color='C1', alpha=0.5, label='Soil C, 20 kg')

ax.plot(df['Sieve size Ø [mm]'], df['Soil A (ISO)'],
        lw=3, color='C0', label='Soil A (ISO), 200 g')
ax.plot(df['Sieve size Ø [mm]'], df['Soil A (5g)'],
        lw=1.5, color='C0', alpha=0.5, label='Soil A, 5 g')
ax.legend(loc='upper left')

plt.tight_layout()
plt.savefig(r'../laboratory/lab_tests_PSDs.jpg', dpi=600)
plt.close()

###############################
# plot showing real lab results: scatter plot
###############################

fp = r'../laboratory/LabResults.xlsx'
df = pd.read_excel(fp, skiprows=20, nrows=2, usecols=list(range(1, 9)))
headers[0] = 'parameter'
df.columns = headers
df.set_index('parameter', inplace=True)

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))

ms = 100
min_, max_ = df.loc['Cc'].min(), df.loc['Cc'].max()
min_, max_ = min_-min_*0.1, max_+max_*0.1

ax1.scatter(df['Soil A (ISO)'].loc['Cc'], df['Soil A (100g)'].loc['Cc'],
            label='Soil A (100g)', s=ms, marker='o', color='C0')
ax1.scatter(df['Soil A (ISO)'].loc['Cc'], df['Soil A (75g)'].loc['Cc'],
            label='Soil A (75g)', s=ms, marker='o', color='C1')
ax1.scatter(df['Soil A (ISO)'].loc['Cc'], df['Soil A (50g)'].loc['Cc'],
            label='Soil A (50g)', s=ms, marker='o', color='C2')
ax1.scatter(df['Soil A (ISO)'].loc['Cc'], df['Soil A (5g)'].loc['Cc'],
            label='Soil A (5g)', s=ms, marker='o', color='C3')

ax1.scatter(df['Soil C (ISO)'].loc['Cc'], df['Soil C'].loc['Cc'],
            label='Soil C (20kg)', s=ms, marker='P', color='C0')
ax1.plot([min_, max_], [min_, max_], color='black')
middle = min_ + (max_ - min_) / 2
ax1.text(x=middle, y=middle, s='1 : 1 line', rotation=45, va='top',
         ha='center')

ax1.set_xlim(left=min_, right=max_)
ax1.set_ylim(bottom=min_, top=max_)
ax1.set_aspect('equal')
ax1.grid(alpha=0.5)
ax1.set_title('Coefficient of Curvature')
ax1.set_xlabel('tests with ISO sample weight')
ax1.set_ylabel('tests with lower sample weight')
ax1.legend()


min_, max_ = df.loc['Cu'].min(), df.loc['Cu'].max()
min_, max_ = min_-min_*0.1, max_+max_*0.1

ax2.scatter(df['Soil A (ISO)'].loc['Cu'], df['Soil A (100g)'].loc['Cu'],
            label='Soil A (100g)', s=ms, marker='o', color='C0')
ax2.scatter(df['Soil A (ISO)'].loc['Cu'], df['Soil A (75g)'].loc['Cu'],
            label='Soil A (75g)', s=ms, marker='o', color='C1')
ax2.scatter(df['Soil A (ISO)'].loc['Cu'], df['Soil A (50g)'].loc['Cu'],
            label='Soil A (50g)', s=ms, marker='o', color='C2')
ax2.scatter(df['Soil A (ISO)'].loc['Cu'], df['Soil A (5g)'].loc['Cu'],
            label='Soil A (5g)', s=ms, marker='o', color='C3')

ax2.scatter(df['Soil C (ISO)'].loc['Cu'], df['Soil C'].loc['Cu'],
            label='Soil C (20kg)', s=ms, marker='P', color='C0')
ax2.plot([min_, max_], [min_, max_], color='black')
middle = min_ + (max_ - min_) / 2
ax2.text(x=10, y=10, s='1 : 1 line', rotation=45, va='top',
         ha='center')

ax2.set_xlim(left=min_, right=max_)
ax2.set_ylim(bottom=min_, top=max_)
ax2.set_aspect('equal')
ax2.grid(alpha=0.5)
ax2.set_title('Coefficient of Uniformity')
ax2.set_xlabel('tests with ISO sample weight')
ax2.set_ylabel('tests with lower sample weight')
ax2.legend()
ax2.set_xscale('log')
ax2.set_yscale('log')

plt.tight_layout()
plt.savefig(r'../laboratory/lab_tests_CcCu_scatter.jpg', dpi=600)
plt.close()

###############################
# other analyses of the simulations
###############################

for f in ['Cu', 'Cc', 'S0', 'max diameter [mm]', 'd10', 'd12', 'd25',
          'd30', 'd50', 'd60', 'd75', 'd90']:
    cc = np.corrcoef(df[f], df['req. weight ks_p95 <= 10 [kg]'])[0][1]
    print(f'{round(cc, 2)} = corr. coeff. {f} - required weight')
