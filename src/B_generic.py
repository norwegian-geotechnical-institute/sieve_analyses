# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 12:08:40 2024

@author: GEr
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import wasserstein_distance, energy_distance

from X_library import laboratory


N = 10_000
DISTRIBUTION = 'normal'  # normal, exponential, beta, uniform, lognormal, combined
DENSITY = 2.5  # grain density [g/cm3]

# instantiate random sampler
rng = np.random.default_rng()

lab = laboratory()

grain_diameters = lab.make_grains(
    n=N, density=DENSITY, distribution=DISTRIBUTION)[0]

percentages = np.arange(0.01, 1, step=0.02)
samples = []
distances = []

for percentage in percentages:
    print(percentage)
    sample_size = int(N*percentage)
    sample = rng.choice(grain_diameters, size=sample_size, replace=False)
    samples.append(sample)
    distances.append(energy_distance(grain_diameters, sample))

fig, (ax1, ax2) = plt.subplots(figsize=(10, 6), nrows=1, ncols=2)

ax1.hist(grain_diameters, bins=30, edgecolor='black', label='data',
         density=True)
ax1.hist(samples[0], bins=30, edgecolor='black', label='10% sample',
         alpha=0.6, density=True)
ax1.legend()
ax1.grid(alpha=0.5)

ax2.scatter(percentages*100, distances)
ax2.grid(alpha=0.5)
ax2.set_xlabel('sample sizes [%]')
ax2.set_ylabel('distance between distributions')

plt.tight_layout()
plt.savefig(fr'../figures/generic_{DISTRIBUTION}.svg')
plt.close()
