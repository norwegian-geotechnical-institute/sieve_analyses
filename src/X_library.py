# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:16:46 2024

@author: GEr
"""

import matplotlib.pyplot as plt
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
            diameters = np.random.normal(loc=0, scale=0.1, size=n)  # [mm]
            diameters = np.where(diameters < 0, diameters*-1, diameters)
        elif distribution == 'exponential':
            diameters = np.random.exponential(scale=3, size=n)  # [mm]
        elif distribution == 'beta':
            diameters = np.random.beta(1, 3, size=n) * 100
        elif distribution == 'uniform':
            diameters = np.random.uniform(0, 100, size=n)
        elif distribution == 'lognormal':
            diameters = np.random.lognormal(mean=-2, sigma=1.5, size=n)
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


class plotter(laboratory):

    def __init__(self):
        pass

    def required_weight_plot(self, savepath):
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
