# -*- coding: utf-8 -*-
"""
Custom libraries for simulated sieve analyses to develop a new way to determine
minimum required soil sample masses.
Libraries:
    laboratory ... contains functions for simulated lab analyses
    statistics ... contains statistical functions for analyses
    plotter ... makes figures

Author: Georg H. Erharter (georg.erharter@ngi.no)
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd


class laboratory:

    def ISO_required_sample_weight(self, d_max: float) -> float:
        '''calculate required sample weight acc. to ISO 17892-4'''
        if d_max <= 2:  # [mm]
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

    def ASTM_required_sample_weight(self, d_max: float) -> float:
        '''calculate required sample weight acc. to ASTM D6913/D6913M − 17,
        method A'''
        density = 3  # [g/mm3] ... ASTM seems to use this high grain density

        if d_max <= 4.75:
            weight = 50
        elif d_max <= 9.5:
            weight = 75
        elif d_max <= 76.2:
            v = (4/3)*np.pi*(d_max/2)**3  # [mm3]
            v = v / 1000  # [cm3]
            m = v * density / 1000  # [kg]
            factor = 1.2
            weight = m * 100 * factor * 1000
        else:
            v = (4/3)*np.pi*(d_max/2)**3  # [mm3]
            v = v / 1000  # [cm3]
            m = v * density / 1000  # [kg]
            weight = m * 100 * 1000
        return weight  # [g]

    def new_sample_weight(self, S0) -> float:
        '''tests for new sample weight calculations'''
        return (np.log(np.log(S0)) + 3.22) * 25

    def get_sample(self, required_sample_weight: float, total_weight: float,
                   grain_weights: np.array, grain_diameters: np.array) -> list:
        '''get a sample with a specific weight from the total number of grains
        by gathering single grains until th desired weight is reached'''
        fraction = required_sample_weight / total_weight
        split = int(len(grain_diameters)*fraction)
        sample_ids = np.arange(split)
        sample_weights = grain_weights[:split]
        sample_diameters = grain_diameters[:split]
        return sample_ids, sample_weights, sample_diameters

    def make_grains(self, DENSITY: float, TOT_MASS: float,
                        verbose: bool = False) -> list:
        '''function generates a new soil sample'''

        def make_beta(lower, upper, n):
            x = np.random.beta(a=np.random.uniform(1, 5),
                               b=np.random.uniform(1, 5), size=n)
            x = x*(upper-lower)
            return x + lower

        COMPONENTS = np.array(['Sa', 'fGr', 'mGr', 'cGr', 'Co'])

        fractions = np.random.uniform(0, 1, np.random.randint(
            2, len(COMPONENTS)+1))
        fractions = fractions / sum(fractions)
        fractions = np.sort(fractions)

        components = np.random.choice(COMPONENTS, len(fractions),
                                      replace=False)
        if verbose is True:
            print(components)
            print(fractions)

        diameters = []
        tot_weight = 0

        for i, component in enumerate(components):
            while tot_weight < sum(fractions[:i+1]) * TOT_MASS:
                match component:  # noqa
                    case 'Sa':
                        diameter = make_beta(0.063, 2, 1000)
                    case 'fGr':
                        diameter = make_beta(2, 6.3, 100)
                    case 'mGr':
                        diameter = make_beta(6.3, 20, 50)
                    case 'cGr':
                        diameter = make_beta(20, 63, 10)
                    case 'Co':
                        diameter = make_beta(63, 200, 1)
                volume = ((4/3)*np.pi*((diameter/2)**3)) / 1000  # [cm3]
                try:
                    tot_weight += volume.sum() * DENSITY / 1000
                except AttributeError:
                    tot_weight += volume * DENSITY / 1000
                diameters.append(diameter)

        diameters = np.hstack(diameters)
        np.random.shuffle(diameters)
        volumes = (4/3)*np.pi*(diameters/2)**3  # [mm3]
        volumes = volumes / 1000  # [cm3]
        weights = volumes * DENSITY / 1000  # [kg]
        ids = np.arange(len(diameters))
        return diameters, weights, ids

    def sieve(self, ids: np.array, diameters: np.array, weights: np.array,
              sieve_sizes: list) -> list:
        '''make a virtual sieve analysis of the sample'''
        tot_weight = weights.sum()
        fractions = []
        for size in sieve_sizes:
            fraction_ids = np.where(diameters < size)[0]
            fraction_weight = weights[fraction_ids].sum()
            try:
                fraction = (fraction_weight / tot_weight) * 100
            except ZeroDivisionError:
                fraction = 0
            fractions.append(fraction)
        return fractions

    def calc_grading_characteristics(self, fractions_true: list,
                                     SIEVE_SIZES: list) -> list:
        '''function computes different metrics about the grading of a soil'''
        d_keys = ['d10', 'd12', 'd25', 'd30', 'd50', 'd60', 'd75', 'd90']
        d_vals = np.interp([10, 12, 25, 30, 50, 60, 75, 90], fractions_true,
                           SIEVE_SIZES)
        ds = dict(zip(d_keys, d_vals))

        Cu = ds['d60']/ds['d10']
        Cc = (ds['d30']**2)/(ds['d60']*ds['d10'])
        S0 = np.sqrt(ds['d75']/ds['d25'])
        return ds, Cu, Cc, S0

    def USCS_classification(self, d12: float, d50: float, Cu: float,
                            Cc: float) -> str:
        '''function determines the soil class according to the unified soil
        classification system'''
        # 50% or more passes the no. 200 Sieve
        if d50 <= 0.075:  # [mm]
            soil_class = 'fines'
        # 50% or more of coarse fraction passes No.4 sieve
        elif d50 <= 4.75:
            if d12 <= 0.075:
                soil_class = 'S_dirty'
            else:
                if Cu >= 6 and Cc <= 3 and Cc >= 1:
                    soil_class = 'SW'
                elif Cu < 6 or Cc < 1 or Cc > 3:
                    soil_class = 'SP'
        # More than 50% of coarse fraction on No. 4 Sieve
        elif d50 > 4.75:
            if d12 <= 0.075:
                soil_class = 'G_dirty'
            else:
                if Cu >= 4 and Cc <= 3 and Cc >= 1:
                    soil_class = 'GW'
                elif Cu < 4 or Cc < 1 or Cc > 3:
                    soil_class = 'GP'
        else:
            soil_class = 'undefined'

        return soil_class


class utilities:

    def CheckForLess(self, list1, val):
        '''function checks if any of values in list1 is smaller than val'''
        for x in list1:
            if val < x:
                return False
        return True


class statistics:

    def ks_statistic(self, sieve_curve1, sieve_curve2):
        '''calculate Kolmogorov-Smirnov (KS) statistic for 2 sieve curves'''
        # Calculate the maximum vertical distance between the sieve curves
        return np.max(np.abs(np.array(sieve_curve1) - np.array(sieve_curve2)))


class plotter(laboratory):

    def required_weight_plot(self, max_grain_size: float, savepath: str,
                             close: bool = True) -> None:
        '''plot that shows the theoretically required sample weight acc. to the
        standards'''
        sizes = np.arange(max_grain_size)  # [mm]
        ms_ISO = [self.ISO_required_sample_weight(ds)/1000 for ds in sizes]
        ms_ASTM = [self.ASTM_required_sample_weight(ds)/1000 for ds in sizes]

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.plot(sizes, ms_ISO, label='ISO 17892-4')
        ax.plot(sizes, ms_ASTM, label='ASTM D6913/D6913M − 17')
        ax.grid(alpha=0.5)
        ax.set_xlabel('max. grain diameter of soil [mm]')
        ax.set_ylabel('required sample weight [kg]')
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
        ax.legend()

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def distances_plot(self, req_sample_weights: list, ks_distances: list,
                       savepath: str, close: bool = True) -> None:
        '''plot an underlying soil distribution vs. distributions of different
        samples with reduced mass and also their kolmogorov smirnov distance'''
        fig, ax = plt.subplots(figsize=(6, 6), nrows=1, ncols=1)

        ax.scatter(req_sample_weights, ks_distances, color='C0',
                   edgecolor='black', s=60)
        ax.grid(alpha=0.5)
        ax.set_xlabel('sample weight [kg]')
        ax.set_ylabel('kolmogorov smirnov distance', color='C0')

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def sieve_curves_plot(self, SIEVE_SIZES: list, fractions_true: list,
                          savepath: str,
                          sieved_samples: list = None,
                          req_sample_weights: list = None,
                          ks_distances: list = None,
                          close: bool = True) -> None:
        '''plot sieve curves of underlying soil distribution and taken
        samples'''
        fig, ax = plt.subplots(figsize=(10, 5))
        if len(np.array(fractions_true).shape) == 1:
            ax.plot(SIEVE_SIZES, fractions_true, label="sieve curve real",
                    color='black', lw=3)
        else:
            for f in fractions_true:
                ax.plot(SIEVE_SIZES, f, color='black', alpha=0.2)

        if sieved_samples is not None:
            for i in range(len(sieved_samples)):
                ax.plot(SIEVE_SIZES, sieved_samples[i],
                        label=f"sample {round(req_sample_weights[i], 1)}kg,\
                        ks distance: {round(ks_distances[i], 3)}",
                        alpha=0.8)

        ax.set_xscale('log')
        ax.set_xlim(left=0.002, right=630)
        ax.set_ylim(bottom=0, top=101)
        ax.set_xticks([0.006, 0.02, 0.06, 2, 63])
        ax.set_xlabel('grain size [mm]')
        ax.set_ylabel('[%]')
        ax.grid(alpha=0.5)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        ax.legend(loc='upper left')

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def monte_carlo_scatterplot(self, df: pd.DataFrame, weight_mode: str,
                                color_mode: str, savepath: str,
                                close: bool = True) -> None:
        '''scatterplot that shows results of the monte carlo simulations in the
        form of different soil distribution parameters Cu, Cc etc. vs. a metric
        of how well the sample fits the underlying distribution'''
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=1, ncols=4,
                                                 figsize=(20, 6))
        match color_mode:  # noqa
            case 'USCS':
                for sc in df.groupby('USCS soil classes'):
                    ax1.scatter(
                        sc[1]['Cu'],
                        sc[1][f'kolmogorov smirnov distance {weight_mode}'],
                        alpha=0.5, label=sc[0])
                    ax2.scatter(
                        sc[1]['Cc'],
                        sc[1][f'kolmogorov smirnov distance {weight_mode}'],
                        alpha=0.5, label=sc[0])
                    ax3.scatter(
                        sc[1]['S0'],
                        sc[1][f'kolmogorov smirnov distance {weight_mode}'],
                        alpha=0.5, label=sc[0])
                    ax4.scatter(
                        sc[1]['max diameter [mm]'],
                        sc[1][f'kolmogorov smirnov distance {weight_mode}'],
                        alpha=0.5, label=sc[0])
                ax1.legend()
                ax2.legend()
                ax3.legend()
                ax4.legend()
            case 'weight':
                ax1.scatter(df['Cu'],
                            df[f'kolmogorov smirnov distance {weight_mode}'],
                            c=df[f'req. weight {weight_mode} [kg]'], alpha=0.5)
                ax2.scatter(df['Cc'],
                            df[f'kolmogorov smirnov distance {weight_mode}'],
                            c=df[f'req. weight {weight_mode} [kg]'], alpha=0.5)
                ax3.scatter(df['S0'],
                            df[f'kolmogorov smirnov distance {weight_mode}'],
                            c=df[f'req. weight {weight_mode} [kg]'], alpha=0.5)
                im = ax4.scatter(
                    df['max diameter [mm]'],
                    df[f'kolmogorov smirnov distance {weight_mode}'],
                    c=df[f'req. weight {weight_mode} [kg]'], alpha=0.5)
                divider = make_axes_locatable(ax4)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                fig.colorbar(im, cax=cax, orientation='vertical',
                             label='required sample weight [kg]')

        ax1.grid(alpha=0.5)
        ax1.set_ylabel(f'kolmogorov smirnov distance {weight_mode}')
        ax1.set_xlabel('coefficient of uniformity\nCu  d60/d10')
        ax1.set_yscale('log')
        ax1.set_ylim(top=100, bottom=0.01)

        ax2.grid(alpha=0.5)
        ax2.set_xlabel('coefficient of curvature\nCc  d30**2/(d60*d10)')
        ax2.set_yscale('log')
        ax2.set_ylim(top=100, bottom=0.01)

        ax3.grid(alpha=0.5)
        ax3.set_xlabel('sorting coefficient\nS0 sqrt(d75/d25)')
        ax3.set_yscale('log')
        ax3.set_ylim(top=100, bottom=0.01)

        ax4.grid(alpha=0.5)
        ax4.set_xlabel('maximum soil grainsize [mm]')
        ax4.set_yscale('log')
        ax4.set_ylim(top=100, bottom=0.01)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def identify_overlapping(self, grain_xs, grain_ys, grain_rs):
        # check for overlapping grains
        dx = grain_xs - grain_xs[:, np.newaxis]
        dy = grain_ys - grain_ys[:, np.newaxis]
        distances = np.sqrt(dx**2 + dy**2)
        radii_sum = grain_rs + grain_rs[:, np.newaxis]
        overlapping = distances < radii_sum
        # remove self overlapping
        np.fill_diagonal(overlapping, False)
        # get indices
        return np.where(np.any(overlapping, axis=1) == True)[0]

    def shake(self, grain_xs, grain_ys, grain_rs, BOX_SIZE):
        delta_x = np.random.normal(loc=0, scale=grain_rs/10, size=len(grain_xs))
        delta_y = np.random.normal(loc=0, scale=grain_rs/10, size=len(grain_ys))
        grain_xs += delta_x
        grain_ys += delta_y
        grain_xs = np.where(grain_xs < 0, grain_xs*-1, grain_xs)
        grain_xs = np.where(grain_xs > BOX_SIZE,
                            BOX_SIZE + (BOX_SIZE - grain_xs), grain_xs)
        grain_ys = np.where(grain_ys < 0, grain_ys*-1, grain_ys)
        grain_ys = np.where(grain_ys > BOX_SIZE,
                            BOX_SIZE + (BOX_SIZE - grain_ys), grain_ys)
        return grain_xs, grain_ys

    def adjust_grain(self, overlap_ids, grain_xs, grain_ys, grain_rs,
                     BOX_SIZE):
        for overlap_id in overlap_ids:
            overlapping = True
            # update position until it fits
            counter = 0
            while overlapping is True:
                counter += 1
                if counter > 40_000:  # safety to avoid infinite loops
                    break
                # choose new position from grain randomly
                delta_x = np.random.normal(loc=0, scale=grain_rs[overlap_id])
                delta_y = np.random.normal(loc=0, scale=grain_rs[overlap_id])
                new_x = grain_xs[overlap_id] + delta_x
                new_y = grain_ys[overlap_id] + delta_y
                if new_x > 0 and new_x < BOX_SIZE:
                    grain_xs[overlap_id] = new_x
                if new_y > 0 and new_y < BOX_SIZE:
                    grain_ys[overlap_id] = new_y
                # grain_xs[overlap_id] = np.random.uniform(0, BOX_SIZE)
                # grain_ys[overlap_id] = np.random.uniform(0, BOX_SIZE)
                # check if new position overlaps or not
                dx_grain = grain_xs[overlap_id] - grain_xs
                dy_grain = grain_ys[overlap_id] - grain_ys
                distances_grain = np.sqrt(dx_grain**2 + dy_grain**2) - grain_rs
                overlaps_grain = distances_grain < grain_rs[overlap_id]
                overlaps_grain[overlap_id] = False
                if True in overlaps_grain:
                    overlapping = True
                else:
                    overlapping = False

        return grain_xs, grain_ys

    def plot_grains(self, grain_diameters, BOX_SIZE, SEED, Cu, S0, weight,
                    savepath: str, close: bool = True) -> None:
        # preprocessing of grains and initialize positions
        grain_diameters = grain_diameters / 1000
        grain_rs = grain_diameters / 2

        if (grain_rs**2 * np.pi).sum() > BOX_SIZE**2 * 0.66:
            raise ValueError('too small plotting area')

        grain_xs = np.random.uniform(0, BOX_SIZE, len(grain_diameters))
        grain_ys = np.random.uniform(0, BOX_SIZE, len(grain_diameters))
        # identify overlapping ones
        overlapping_indexes = self.identify_overlapping(grain_xs, grain_ys,
                                                        grain_rs)
        n_overlapping = len(overlapping_indexes)
        print(f'plotting grains; there are {n_overlapping} overlapping grains')

        n_overlappings = []  # track how many grains still overlap

        # adjust overlapping grains
        while n_overlapping > 0:
            # number of smallest overlapping grains to relocate
            if len(n_overlappings) > 20 and np.mean(np.diff(n_overlappings[-10:])) > -1:
                n_smallest_grains = 1
            else:
                n_smallest_grains = int(n_overlapping * 0.1)
                if n_smallest_grains < 1: n_smallest_grains = 1

            print(n_overlapping, n_smallest_grains,
                  np.mean(np.diff(n_overlappings[-10:])))
            id_smallest = np.argsort(grain_diameters[overlapping_indexes])[:n_smallest_grains]
            # id_smallest = np.argmin(grain_diameters[overlapping_indexes])
            id_smallest = overlapping_indexes[id_smallest]
            grain_xs, grain_ys = self.adjust_grain(id_smallest, grain_xs,
                                                   grain_ys, grain_rs,
                                                   BOX_SIZE)

            overlapping_indexes = self.identify_overlapping(grain_xs,
                                                            grain_ys,
                                                            grain_rs)
            n_overlapping = len(overlapping_indexes)
            n_overlappings.append(n_overlapping)

            if np.mean(np.diff(n_overlappings[-5:])) >= 0:
                grain_xs, grain_ys = self.shake(grain_xs, grain_ys,
                                                grain_rs, BOX_SIZE)
                print('shaked')
        print('plotting')
        # plot grains
        fig, ax = plt.subplots(figsize=(9, 9))
        for i in range(len(grain_xs)):
            if i in overlapping_indexes:
                circle = plt.Circle((grain_xs[i], grain_ys[i]), grain_rs[i],
                                    facecolor='none', edgecolor='red')
            else:
                circle = plt.Circle((grain_xs[i], grain_ys[i]), grain_rs[i],
                                    color='black', linewidth=0)
            ax.add_patch(circle)
        ax.set_xlim(0, BOX_SIZE)
        ax.set_ylim(0, BOX_SIZE)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.set_title(f'seed {SEED}; Cu: {round(Cu, 2)}, S0: {round(S0, 2)}, weight: {weight} kg')
        ax.set_xlabel(f'{BOX_SIZE} meters')

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()


if __name__ == '__main__':
    DENSITY = 2.65
    TOT_MASS = 50
    SIEVE_SIZES = [0.075, 0.105, 0.15, 0.25, 0.425, 0.85, 2, 4.75, 9.5, 19, 25, 37.5, 50, 75, 100, 150, 200, 300]  # sieve sizes [mm]

    lab, pltr = laboratory(), plotter()

    # make plot of random sieve curves to demonstrate sampler
    fractions_trues = []

    for _ in range(200):
        grain_diameters, grain_weights, grain_ids = lab.make_grains_new(
            DENSITY, TOT_MASS=TOT_MASS)
        print(TOT_MASS, grain_weights.sum())
        fractions_true = lab.sieve(grain_ids, grain_diameters, grain_weights,
                                   SIEVE_SIZES)
        fractions_trues.append(fractions_true)

    pltr.sieve_curves_plot(SIEVE_SIZES, fractions_trues,
                           savepath=fr'../figures/sieve_samples.jpg')
