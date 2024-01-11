# -*- coding: utf-8 -*-
"""
Script to make virtual grain size analyses with different sample sizes

Author: Georg H. Erharter (georg.erharter@ngi.no)
"""

# importing libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

from X_library import laboratory, plotter


###############################
# main code execution
###############################

# constant values and hyperparameters
N = 10_000_000  # number of grains to choose from
DENSITY = 2.5  # grain density [g/cm3]
SIEVE_SIZES = [0.04, 0.1, 0.25, 0.425, 0.85, 2, 4.75, 10, 20, 50, 100, 150, 200, 300]  # sieve sizes [mm]
DISTRIBUTION = 'exponential'  # 'normal', 'exponential', 'beta', 'uniform', 'lognormal'
SEED = 3

# instantiations
lab, pltr = laboratory(), plotter()

# plot for GBV to visualize theoretically required sample mass
# pltr.required_weight_plot('../figures/required_weight.svg')

np.random.seed(SEED)  # fix seed for reproducibility

# initialize ground truth sample
grain_diameters, grain_volumes, grain_weights, grain_ids = lab.make_grains(
    N, DENSITY, DISTRIBUTION)

# calculate required sample size
required_sample_weight = lab.calc_required_sample_weight(grain_diameters) / 1000
print(f'total weight: {sum(grain_weights)}')
print(f'required weight: {round(required_sample_weight, 3)} kg')
# break if more material is required than available
if required_sample_weight > sum(grain_weights):
    raise ValueError('required sample weight larger than total weight')


fractions_true = lab.sieve(grain_ids, grain_diameters, grain_weights, SIEVE_SIZES)

# initialize plot
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(SIEVE_SIZES, fractions_true, label="sieve curve real", color='black',
        lw=5)

d60 = np.interp(60, fractions_true, SIEVE_SIZES)
d30 = np.interp(30, fractions_true, SIEVE_SIZES)
d10 = np.interp(10, fractions_true, SIEVE_SIZES)

Cu = d60/d10
Cc = (d30**2)/(d60*d10)

for i in [1, 2, 10, 20, 100]:
    sample_weight = required_sample_weight/i
    sample_ids, sample_weights, sample_diameters = lab.get_sample(
        sample_weight, grain_weights, grain_diameters)
    fractions_sample = lab.sieve(
        sample_ids, sample_diameters, sample_weights, SIEVE_SIZES)
    r2 = r2_score(np.log(fractions_true), np.log(fractions_sample))
    mae = mean_absolute_error(np.log(fractions_true), np.log(fractions_sample))
    mse = mean_squared_error(np.log(fractions_true), np.log(fractions_sample))
    ax.plot(SIEVE_SIZES, fractions_sample,
            label=f"sample {round(sample_weight, 1)} kg, R2: {round(r2, 3)}, mae: {round(mae, 3)}, mse: {round(mse, 3)}",
            alpha=0.8)

ax.plot([0.002, d60, d60], [60, 60, 0], color='black')
ax.plot([0.002, d30, d30], [30, 30, 0], color='black')
ax.plot([0.002, d10, d10], [10, 10, 0], color='black')

text = f'd60: {round(d60, 1)}\nd30: {round(d30, 1)}\nd10: {round(d10, 1)}\nCu: {round(Cu, 1)}\nCc: {round(Cc, 1)}'
ax.text(x=0.003, y=5, s=text, backgroundcolor='white')

ax.set_xscale('log')
ax.set_xlim(left=0.002, right=630)
ax.set_ylim(bottom=0, top=101)
ax.set_xticks([0.006, 0.02, 0.06, 2, 63])
ax.set_xlabel('grain size [mm]')
ax.set_ylabel('[%]')
ax.set_title(f'max grain size: {round(max(grain_diameters), 1)} mm, ISO required sample weight: {round(required_sample_weight, 1)} kg, distribution: {DISTRIBUTION}')
ax.grid(alpha=0.5)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
ax.legend(loc='upper left')
plt.tight_layout()
plt.savefig(fr'../figures/sample_{DISTRIBUTION}_{SEED}.png')
plt.close()
