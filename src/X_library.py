# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:16:46 2024

@author: GEr
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import numpy as np


class laboratory:

    def __init__(self):
        pass

    def calc_required_sample_weight(self, diameters: np.array) -> float:
        '''calculate required sample weight acc. to ISO 17892-4'''
        d_max = max(diameters)  # [mm]
        if d_max <= 2:
            weight = 100
        elif d_max <= 6.3:
            weight = 300
        elif d_max <= 10:
            weight = 500
        elif d_max <= 20:
            weight = 2000
        else:
            weight = ((d_max / 10)**2) * 1000
        return weight  # [g]

    def get_sample(self, required_sample_weight: float, grain_weight: np.array,
                   grain_diameters: np.array) -> list:
        '''get a sample with a specific weight from the total number of grains
        by gathering single grains until th desired weight is reached'''
        sample_ids = []
        sample_weight, sample_id = 0, 0
        while sample_weight < required_sample_weight:
            sample_ids.append(sample_id)
            sample_weight += grain_weight[sample_id]
            sample_id += 1
        return np.array(sample_ids), grain_weight[sample_ids], grain_diameters[sample_ids]

    def make_grains(self, n: int, density: float, distribution: str) -> list:
        '''make the initial distribution of grains following different
        statistical dsitributions'''
        # initializing random values
        if distribution == 'normal':
            diameters = np.random.normal(loc=0, scale=0.2, size=n)  # [mm]
            diameters = np.where(diameters < 0, diameters*-1, diameters)
        elif distribution == 'exponential':
            diameters = np.random.exponential(scale=3, size=n)  # [mm]
        elif distribution == 'beta':
            diameters = np.random.beta(1, 3, size=n) * 100
        elif distribution == 'uniform':
            diameters = np.random.uniform(0, 100, size=n)
        elif distribution == 'lognormal':
            diameters = np.random.lognormal(mean=0, sigma=1, size=n)
        elif distribution == 'combined':
            clay_to_silt_samples = np.random.lognormal(mean=1.5, sigma=0.5,
                                                       size=n)
            sand_to_gravel_samples = np.random.lognormal(mean=1.0, sigma=0.5,
                                                         size=n)
            gravel_samples = np.random.exponential(scale=1 / 0.2, size=n)
            diameters = np.concatenate((clay_to_silt_samples,
                                        sand_to_gravel_samples,
                                        gravel_samples))
            np.random.shuffle(diameters)

        diameters = np.where(diameters < 0.002, 0.002, diameters)
        diameters = np.where(diameters > 300, 300, diameters)
        print(f'grain size min: {diameters.min()}, max: {diameters.max()} mm')
        volumes = (4/3)*np.pi*(diameters/2)**3  # [mm3]
        volumes = volumes / 1000  # [cm3]
        weights = volumes * density / 1000  # [kg]
        ids = np.arange(n)
        return diameters, volumes, weights, ids

    def sieve(self, ids: np.array, diameters: np.array, weights: np.array,
              sieve_sizes: list) -> list:
        '''make a virtual sieve analysis of the sample'''
        tot_weight = sum(weights)
        fractions = []
        for size in sieve_sizes:
            fraction_ids = np.where(diameters < size)[0]
            fraction_weight = sum(weights[fraction_ids])
            try:
                fraction = (fraction_weight / tot_weight) * 100
            except ZeroDivisionError:
                fraction = 0
            fractions.append(fraction)
        return fractions


class statistics:

    def __init__(self):
        pass

    def calc_cdf(self, data):
        '''calculate cumulative density function of data'''
        count, bins_count = np.histogram(data, bins=50)
        # finding the PDF of the histogram using count values
        pdf = count / sum(count)
        # using numpy np.cumsum to calculate the CDF
        cdf = np.cumsum(pdf)
        return cdf

    def ks_statistic(self, data1, data2):
        '''calculate Kolmogorov-Smirnov (KS) statistic'''
        # Calculate the empirical CDFs of the two datasets
        cdf1 = self.calc_cdf(data1)
        cdf2 = self.calc_cdf(data2)
        # Calculate the maximum vertical distance between the CDFs
        ks_statistic = np.max(np.abs(cdf1 - cdf2))
        return ks_statistic


class plotter(laboratory):

    def __init__(self):
        pass

    def required_weight_plot(self, savepath: str) -> None:
        '''plot that shows the theoretically required sample weight acc. to the
        standards'''
        sizes = np.arange(200)  # [mm]
        weights = [self.calc_required_sample_weight([ds])/1000 for ds in sizes]

        fig, ax = plt.subplots(figsize=(3, 3))
        ax.plot(sizes, weights)
        ax.grid(alpha=0.5)
        ax.set_xlabel('max. grain diameter [mm]')
        ax.set_ylabel('required sample weight [kg]\nacc. to ISO 17892-4')
        plt.tight_layout()
        plt.savefig(savepath)
        plt.close()

    def distances_plot(self, grain_diameters: np.array, sample_diameters: list,
                       req_sample_weights: list, wasserstein_distances: list,
                       energy_distances: list, ks_distances: list,
                       savepath: str) -> None:
        fig, (ax1, ax2) = plt.subplots(figsize=(12, 6), nrows=1, ncols=2)

        ax1.hist(grain_diameters, bins=30, edgecolor='black', color='C0',
                 label='true soil', density=True)
        ax1.hist(sample_diameters[0], bins=30, edgecolor='black', color='C1',
                 label=f'sample acc. standards: {round(req_sample_weights[0], 1)}kg',
                 alpha=0.5, density=True)
        ax1.hist(sample_diameters[-3], bins=30, edgecolor='black', color='C2',
                 label=f'{round(req_sample_weights[-3], 1)}kg',
                 alpha=0.3, density=True)
        ax1.set_xlabel('grain diameters [mm]')
        ax1.legend()
        ax1.grid(alpha=0.5)

        ax2.scatter(req_sample_weights, ks_distances, color='C0',
                    edgecolor='black', s=60)
        ax2.grid(alpha=0.5)
        ax2.set_xlabel('sample weight [kg]')
        ax2.set_ylabel('kolmogorov smirnov distance', color='C0')

        # ax2.scatter(req_sample_weights, wasserstein_distances, color='C0',
        #             edgecolor='black', s=60)
        # ax2.grid(alpha=0.5)
        # ax2.set_xlabel('sample weight [kg]')
        # ax2.set_ylabel('wasserstein distance', color='C0')

        # ax2_2 = ax2.twinx()
        # ax2_2.scatter(req_sample_weights, energy_distances, color='C1',
        #               edgecolor='black', s=60, marker='v')
        # ax2_2.set_ylabel('energy distance', color='C1')

        plt.tight_layout()
        plt.savefig(savepath)
        plt.close()

    def sieve_curves_plot(self, SIEVE_SIZES: list, fractions_true: list,
                          sieved_samples: list, req_sample_weights: list,
                          wasserstein_distances: list, energy_distances: list,
                          d60: float, d30: float, d10: float, Cu: float,
                          Cc: float, grain_diameters: np.array,
                          standard_sample_weight: float, DISTRIBUTION: str,
                          savepath: str) -> None:
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(SIEVE_SIZES, fractions_true, label="sieve curve real",
                color='black', lw=5)

        for i in range(len(sieved_samples)):
            ax.plot(SIEVE_SIZES, sieved_samples[i],
                    label=f"sample {round(req_sample_weights[i], 1)}kg,\
                    wd: {round(wasserstein_distances[i], 3)},\
                    ed: {round(energy_distances[i], 3)}",
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
        ax.set_title(f'max grain size: {round(max(grain_diameters), 1)} mm, ISO required sample weight: {round(standard_sample_weight, 1)} kg, distribution: {DISTRIBUTION}')
        ax.grid(alpha=0.5)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax.legend(loc='upper left')
        plt.tight_layout()
        plt.savefig(savepath)
        plt.close()
