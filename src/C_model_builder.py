# -*- coding: utf-8 -*-
"""
Script to primarily create visualizations as a means to build equations of
relationships between parameters.

Author: Georg H. Erharter (georg.erharter@ngi.no)
"""

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd

from X_library import laboratory, plotter


def exponential(x, a, b):
    return a * np.exp(b*x)


def fit_exponential(x, y):
    # Perform the curve fit
    p0 = [100, -1]
    params, _ = curve_fit(exponential, x, y, p0=p0)

    # Extract the fitted parameters
    a_fit, b_fit = params
    return a_fit, b_fit


STUDY_NAME = '2024_03_17'  # study to work with or to create

df = pd.read_excel(fr'../simulations/{STUDY_NAME}.xlsx')
print(len(df))
lab, pltr = laboratory(), plotter()

# plot to visualize theoretically required sample mass
pltr.required_weight_plot(300, '../figures/required_weight.svg')

###############################
# plot showing errors of ASTM and ISO
###############################

med_ISO = round(df['ISO ks [%]'].median(), 1)
p95_ISO = round(np.percentile(df['ISO ks [%]'], 95), 1)
med_ASTM = round(df['ASTM ks [%]'].median(), 1)
p95_ASTM = round(np.percentile(df['ASTM ks [%]'], 95), 1)

fig, ax = plt.subplots()

parts = ax.violinplot([df['ISO ks [%]'], df['ASTM ks [%]']], showextrema=False,
                      points=100, widths=0.7)
for pc in parts['bodies']:
    pc.set_facecolor('grey')
    pc.set_edgecolor('black')
    pc.set_alpha(1)
ax.plot([1, 1], np.percentile(df['ISO ks [%]'], [5, 95]), color='black')
ax.plot([2, 2], np.percentile(df['ASTM ks [%]'], [5, 95]), color='black')
ax.scatter([1, 2], [med_ISO, med_ASTM], color='white', edgecolor='black',
           s=100, zorder=10)
ax.scatter([1, 2], [p95_ISO, p95_ASTM], color='black', edgecolor='black',
           s=50, zorder=10)
ax.text(x=1.05, y=med_ISO, s=f'{med_ISO}%')
ax.text(x=1.05, y=p95_ISO, s=f'{p95_ISO}%')
ax.text(x=2.05, y=med_ASTM, s=f'{med_ASTM}%')
ax.text(x=2.05, y=p95_ASTM, s=f'{p95_ASTM}%')

ax.grid(alpha=0.3)
ax.set_ylabel('Kolmogorov-Smirnov statisic [mass %]')
ax.set_ylim(top=-1, bottom=15)
ax.set_xticks([1, 2], labels=['ISO', 'ASTM'])

plt.tight_layout()
plt.savefig(r'../figures/error_violin.jpg')
plt.close()


###############################
# plot showing required sample mass against max diameter
###############################

x = np.arange(200)
req_ISO = [lab.ISO_required_sample_weight(dmax)/1000 for dmax in x]
req_ASTM = [lab.ASTM_required_sample_weight(dmax)/1000 for dmax in x]

fig, ax = plt.subplots()

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)

im = ax.scatter(df['max diameter [mm]'], df['req. weight ks_p95 <= 10 [kg]'],
                c=df['S0'], cmap='Greys', edgecolor='black',
                label='simulated samples')
ax.plot(x, req_ISO, color='black', label='ISO req. sample mass')
ax.plot(x, req_ASTM, color='black', ls='--', label='ASTM req. sample mass')
ax.set_xlabel('max diameter [mm]')
ax.set_ylabel('req. sample weight [kg]')
ax.grid(alpha=0.5)
ax.set_ylim(bottom=0, top=400)
ax.legend()

fig.colorbar(im, cax=cax, orientation='vertical',
             label='sorting coefficient - S0')
plt.tight_layout()
plt.savefig(r'../figures/req_weight_dmax.jpg')
plt.close()


###############################
# plot showing required sample mass against d90
###############################

for f in ['Cu', 'Cc', 'S0', 'max diameter [mm]', 'd10', 'd12', 'd25',
          'd30', 'd50', 'd60', 'd75', 'd90']:
    cc = np.corrcoef(df[f], df['req. weight ks_p95 <= 10 [kg]'])[0][1]
    print(f'{round(cc, 2)} = corr. coeff. {f} - required weight')

fig, ax = plt.subplots()

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)

im = ax.scatter(df['d90'], df['req. weight ks_p95 <= 10 [kg]'],
                c=df['max diameter [mm]'], cmap='Greys', edgecolor='black')
ax.set_xlabel('d90 [mm]')
ax.set_ylabel('req. sample weight [kg]')
ax.grid(alpha=0.5)
ax.set_ylim(bottom=0, top=400)

fig.colorbar(im, cax=cax, orientation='vertical', label='max diameter [mm]')
plt.tight_layout()
plt.savefig(r'../figures/req_weight_d90.jpg')
plt.close()

###############################
# plot showing relationship between exponent and errors
###############################

# compute statistics
results = {}
for exp in np.arange(1, 2.6, step=0.1):
    exp = round(exp, 1)
    results[exp] = [df[f'new {exp} ks [%]'].median(),
                    np.percentile(df[f'new {exp} ks [%]'], 95)]

exponents = list(results.keys())
exponents.reverse()
med_errors = [results[x][0] for x in exponents]
p95_errors = [results[x][1] for x in exponents]

med_a, med_b = fit_exponential(exponents, med_errors)
p95_a, p95_b = fit_exponential(exponents, p95_errors)


fig, ax = plt.subplots()

ax.scatter(exponents, med_errors, label='median error', color='C0')
label = f'error = {round(med_a, 2)}e^({round(med_b, 2)}*exp)'
ax.plot(exponents, exponential(np.array(exponents), med_a, med_b),
        label=label, color='C0')
ax.scatter(exponents, p95_errors, label='p95 error', color='C1')
label = f'error = {round(p95_a, 2)}e^({round(p95_b, 2)}*exp)'
ax.plot(exponents, exponential(np.array(exponents), p95_a, p95_b),
        label=label, color='C1')

ax.grid()
ax.legend()
ax.set_xlabel('exponent')
ax.set_ylabel('error')

plt.tight_layout()
plt.savefig(r'../figures/exponent_vs_error.jpg')
plt.close()

###############################
# plot showing new m_min functions with different epsilon
###############################

x = np.arange(200)

fig, ax = plt.subplots()

for exp in list(results.keys()):
    y = ((x / 10)**exp)
    ax.plot(x, y,
            label=f'ε {exp}; med.e.:{round(results[exp][0], 1)}, p95 e.:{round(results[exp][1], 1)}')

ax.grid(alpha=0.4)
ax.legend()
ax.set_ylabel('req. sample weight (m) [kg]')
ax.set_xlabel('d90 [mm]')
plt.tight_layout()
plt.savefig(r'../figures/exponents.jpg')
plt.close()

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
            color='grey', edgecolor='black', alpha=0.5)
ax1.plot([0, 400], [0, 400], color='black', lw=4, ls='--')
ax1.grid(alpha=0.5)
ax1.set_xlabel('ISO req. weight [kg]')
ax1.set_ylabel('req. weight acc. new criterium [kg]')

ax2.scatter(df['ASTM req. weight [kg]'], df['new. req. weight ASTM [kg]'],
            color='grey', edgecolor='black', alpha=0.5)
ax2.plot([0, 1200], [0, 1200], color='black', lw=4, ls='--')
ax2.grid(alpha=0.5)
ax2.set_xlabel('ASTM req. weight [kg]')

plt.tight_layout()
plt.savefig(r'../figures/comparison.jpg')
plt.close()


sieve_cols = [c for c in df.columns if 'mm sieve [m%]' in c]

sieve_sizes = [eval(s.split(' ')[0]) for s in sieve_cols]

fractions_trues = [list(df[sieve_cols].iloc[i].values) for i in range(500)]

pltr.sieve_curves_plot(sieve_sizes, fractions_trues,
                       savepath=r'../figures/sieve_samples_new.jpg',
                       close=True)
