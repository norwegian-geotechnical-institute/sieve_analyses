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

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from tqdm import tqdm


class statistics:

    def ks_statistic(self, sieve_curve1: list, sieve_curve2: list) -> float:
        '''calculate Kolmogorov-Smirnov (KS) statistic for 2 sieve curves'''
        # Calculate the maximum vertical distance between the sieve curves
        return np.max(np.abs(np.array(sieve_curve1) - np.array(sieve_curve2)))


class laboratory(statistics):

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

    def new_sample_weight(self, d_max: float, d90: float,
                          exponent: float = 2) -> float:
        '''new equation to determine sample weight based on d90'''
        if d_max <= 2:  # [mm]
            weight = 100
        elif d_max <= 6.3:
            weight = 300
        elif d_max <= 10:
            weight = 500
        elif d_max <= 20:
            weight = 2000
        else:
            weight = ((d90 / 10)**exponent) * 1000
        return weight  # [g]

    def get_sample(self, required_sample_weight: float, total_weight: float,
                   grain_weights: np.array, grain_diameters: np.array,
                   strategy: str = 'random choice') -> list:
        '''get a sample with a specific weight from the total number of grains
        by gathering single grains until th desired weight is reached'''

        fraction = required_sample_weight / total_weight
        match strategy:  # noqa
            case 'from start':
                split = int(len(grain_diameters)*fraction)
                sample_ids = np.arange(split)
                sample_weights = grain_weights[:split]
                sample_diameters = grain_diameters[:split]
            case 'random choice':
                n_grains = int(len(grain_diameters)*fraction)
                sample_ids = np.random.choice(np.arange(len(grain_diameters)),
                                              size=n_grains, replace=False)
                sample_weights = grain_weights[sample_ids]
                sample_diameters = grain_diameters[sample_ids]
        return sample_ids, sample_weights, sample_diameters

    def make_grains(self, DENSITY: float, TOT_MASS: float,
                    min_d: float = 1, max_d: float = 200,
                    verbose: bool = False) -> list:
        '''function generates a new soil sample'''

        def make_beta(lower, upper, n):
            x = np.random.beta(a=np.random.uniform(1, 5),
                               b=np.random.uniform(1, 5), size=n)
            x = x*(upper-lower)
            return x + lower

        fractions = np.random.uniform(0, 1, np.random.randint(1, 6))
        fractions = fractions / sum(fractions)
        fractions = np.sort(fractions)

        diameters = []
        tot_weight = 0

        for i in range(len(fractions)):
            lower, upper = sorted(np.exp(np.random.uniform(np.log(min_d),
                                                           np.log(max_d), 2)))
            while tot_weight < sum(fractions[:i+1]) * TOT_MASS:
                if upper < 3:
                    n_grains = 2000
                elif upper < 6:
                    n_grains = 200
                elif upper < 18:
                    n_grains = 50
                elif upper < 60:
                    n_grains = 10
                else:
                    n_grains = 1
                diameter = make_beta(lower, upper, n_grains)
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
              sieve_sizes: list) -> dict:
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
        return dict(zip(sieve_sizes, fractions))

    def calc_grading_characteristics(self, fractions_true: list,
                                     SIEVE_SIZES: list) -> list:
        '''function computes different metrics about the grading of a soil'''
        d_s = [10, 12, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90]
        d_keys = [f'd{d}' for d in d_s]
        d_vals = np.interp(d_s, fractions_true, SIEVE_SIZES)
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

    def check_required_weight(self, max_error: float,  # [%]
                              weight_steps: float,  # steps to increase weight
                              tests_per_step: int, total_weight: float,
                              grain_weights: np.array,
                              grain_diameters: np.array, SIEVE_SIZES: list,
                              fractions_true: dict,
                              verbose: bool = False):
        '''function computes the theoretically required mass to achieve a
        defined error by gathering sequentially more soil. Criterium is p95
        percentile of errors below max_error'''
        ks_p95, weight = 2*max_error, 0  # initialize parameters
        while ks_p95 > max_error:
            weight += weight_steps
            weights_temp, ks_s_temp = [], []
            # make multiple tests to reduce randomnes
            for _ in range(tests_per_step):
                sample_ids, sample_weights, sample_diameter = self.get_sample(
                    weight, total_weight, grain_weights, grain_diameters,
                    strategy='random choice')
                fractions = self.sieve(sample_ids, sample_diameter,
                                       sample_weights, SIEVE_SIZES)
                weights_temp.append(sample_weights.sum())
                ks_s_temp.append(self.ks_statistic(
                    list(fractions_true.values()), list(fractions.values())))
            weight_avg = np.mean(weights_temp)
            ks_p95 = np.percentile(ks_s_temp, 95)
            if verbose is True:
                print(f'{round(weight_avg, 1)} kg\t\tp95 {round(ks_p95, 1)}')

        return weight  # required weight [kg]


class utilities:

    def CheckForLess(self, list1, val):
        '''function checks if any of values in list1 is smaller than val'''
        for x in list1:
            if val < x:
                return False
        return True

    def exponential(self, x, a, b):
        return a * np.exp(b*x)

    def fit_exponential(self, x, y):
        # Perform the curve fit
        p0 = [100, -1]
        params, _ = curve_fit(self.exponential, x, y, p0=p0)

        # Extract the fitted parameters
        a_fit, b_fit = params
        return a_fit, b_fit


class plotter(laboratory, utilities):

    def __init__(self):
        self.fsize = 8  # main font size for plots
        self.lss = ['-', '--', ':', '-.'] * 2  # choice of linestyles

    def required_weight_plot(self, max_grain_size: float, savepath: str,
                             close: bool = True) -> None:
        '''plot that shows the theoretically required sample weight acc. to the
        standards'''
        sizes = np.arange(max_grain_size)  # [mm]
        ms_ISO = [self.ISO_required_sample_weight(ds)/1000 for ds in sizes]
        ms_ASTM = [self.ASTM_required_sample_weight(ds)/1000 for ds in sizes]

        fig, ax = plt.subplots(figsize=(3.54331, 3.54331))
        ax.plot(sizes, ms_ISO, color='black', ls='--',
                label='ISO 17892-4', zorder=10)
        ax.plot(sizes, ms_ASTM, color='black', ls='-',
                label='ASTM D6913/D6913M − 17', zorder=10)
        vlines = (20, 63, 200)
        ax.vlines(vlines, ymin=0, ymax=max(ms_ASTM), color='grey',
                  zorder=5)
        for vline in vlines:
            ax.text(x=vline+2, y=max(ms_ASTM), s=f'{vline}\nmm', ha='left',
                    va='top', fontsize=self.fsize)
        ax.grid(alpha=0.5, axis='y')
        ax.set_xlabel('max. grain diameter of soil [mm]', fontsize=self.fsize)
        ax.set_ylabel('required sample weight [kg]', fontsize=self.fsize)
        ax.set_xlim(left=2)
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
        ax.legend(loc='lower right', framealpha=1, fontsize=self.fsize)
        ax.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def error_violin_plot(self, df: pd.DataFrame, savepath: str,
                          close: bool = True) -> None:
        '''violin plots showing the KS-errors of sieve curves of soils sampled
        according to ISO and ASTM standard'''
        med_ISO = round(df['ISO ks [%]'].median(), 1)
        p95_ISO = round(np.percentile(df['ISO ks [%]'], 95), 1)
        med_ASTM = round(df['ASTM ks [%]'].median(), 1)
        p95_ASTM = round(np.percentile(df['ASTM ks [%]'], 95), 1)

        fig, ax = plt.subplots(figsize=(3.54331, 3))

        parts = ax.violinplot([df['ISO ks [%]'], df['ASTM ks [%]']],
                              showextrema=False, points=100, widths=0.7)
        for pc in parts['bodies']:
            pc.set_facecolor('grey')
            pc.set_edgecolor('black')
            pc.set_alpha(1)
        ax.plot([1, 1], np.percentile(df['ISO ks [%]'], [5, 95]),
                color='black')
        ax.plot([2, 2], np.percentile(df['ASTM ks [%]'], [5, 95]),
                color='black')
        ax.scatter([1, 2], [med_ISO, med_ASTM], color='white',
                   edgecolor='black', s=100, zorder=10)
        ax.scatter([1, 2], [p95_ISO, p95_ASTM], color='black',
                   edgecolor='black', s=50, zorder=10)
        ax.text(x=1.07, y=med_ISO, s=r'$KS_{med}$' + f'\n{med_ISO}%',
                fontsize=self.fsize)
        ax.text(x=1.07, y=p95_ISO, s=r'$KS_{p95}$' + f'\n{p95_ISO}%',
                fontsize=self.fsize, va='top')
        ax.text(x=2.07, y=med_ASTM, s=r'$KS_{med}$' + f'\n{med_ASTM}%',
                fontsize=self.fsize)
        ax.text(x=2.07, y=p95_ASTM, s=r'$KS_{p95}$' + f'\n{p95_ASTM}%',
                fontsize=self.fsize, va='top')

        ax.grid(alpha=0.3)
        ax.set_ylabel('Kolmogorov-Smirnov statisic [mass %]',
                      fontsize=self.fsize)
        ax.set_ylim(top=-1, bottom=15)
        ax.set_xticks([1, 2], labels=['ISO 17892-4', 'ASTM D6913/D6913M – 17'],
                      fontsize=self.fsize)
        ax.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def req_sample_mass_vs_dmax_plot(self, df: pd.DataFrame, savepath: str,
                                     annotate_all=False, annotate_some=None,
                                     close: bool = True) -> None:
        dmaxs = np.arange(200)
        req_ISO = [self.ISO_required_sample_weight(dmax)/1000 for dmax in dmaxs]
        req_ASTM = [self.ASTM_required_sample_weight(dmax)/1000 for dmax in dmaxs]

        fig, ax = plt.subplots(figsize=(3.54331, 3))

        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        im = ax.scatter(df['max diameter [mm]'],
                        df['req. weight ks_p95 <= 10 [kg]'],
                        c=df['S0'], cmap='Greys', edgecolor='black', lw=0.2,
                        label='simulated samples', s=8)

        ax.plot(dmaxs, req_ISO, color='dimgrey', ls='--',
                label='ISO req. sample mass')
        ax.plot(dmaxs, req_ASTM, color='dimgrey', ls='-',
                label='ASTM req. sample mass')
        if annotate_all is True:
            for i in range(len(df)):
                ax.text(x=df['max diameter [mm]'].iloc[i],
                        y=df['req. weight ks_p95 <= 10 [kg]'].iloc[i],
                        s=df['ID'].iloc[i])
        if annotate_some is not None:
            for id_ in annotate_some:
                ax.scatter(df['max diameter [mm]'].iloc[id_],
                           df['req. weight ks_p95 <= 10 [kg]'].iloc[id_],
                           edgecolor='black', facecolors='none', s=20, lw=1,
                           zorder=200)
                t = ax.text(x=df['max diameter [mm]'].iloc[id_],
                            y=df['req. weight ks_p95 <= 10 [kg]'].iloc[id_]+10,
                            s=df['ID'].iloc[id_], ha='center',
                            fontsize=self.fsize, weight='bold', zorder=100)
                t.set_bbox(dict(facecolor='white', alpha=0.5, lw=0))

        ax.set_xlabel('max diameter [mm]', fontsize=self.fsize)
        ax.set_ylabel('req. sample weight [kg]', fontsize=self.fsize)
        ax.grid(alpha=0.5)
        ax.set_ylim(bottom=0, top=400)
        ax.legend(fontsize=self.fsize, loc='upper left')

        cbar = fig.colorbar(im, cax=cax, orientation='vertical')

        cbar.set_label(label=r'sorting coefficient - $S_0$',
                       fontsize=self.fsize)
        ticks = cbar.ax.get_yticks()
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticks(ticks)
        cbar.ax.set_yticklabels(ticklabs, fontsize=self.fsize)
        ax.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def simple_sieve_plot(self, df: pd.DataFrame, ids: list, savepath: str,
                          close: bool = True) -> None:
        '''simple sieve curves plot to combine with
        req_sample_mass_vs_dmax_plot'''
        mpercs = [c for c in df.columns if ' mm sieve [m%]' in c]
        ssizes = [eval(s.split(' ')[0]) for s in mpercs]

        fig, ax = plt.subplots(figsize=(3.54331, 3))
        for i, id_ in enumerate(ids):
            ax.plot(ssizes, df[mpercs].iloc[id_].values, color='black',
                    ls=self.lss[i], label=id_)

        ax.set_xscale('log')
        ax.set_xlabel('grain size [mm]', fontsize=self.fsize)
        ax.set_ylabel('mass percentage passing [%]', fontsize=self.fsize)
        ax.set_xticks([2, 63, 200])
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
        ax.grid(alpha=0.5)
        ax.tick_params(axis='both', labelsize=self.fsize)
        ax.legend(fontsize=self.fsize, loc='lower right')

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def req_sample_mass_vs_d90_plot(self, df, savepath: str,
                                    annotate_some=None,
                                    close: bool = True) -> None:
        fig, ax = plt.subplots(figsize=(3.54331, 3))

        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        im = ax.scatter(df['d90'], df['req. weight ks_p95 <= 10 [kg]'],
                        c=df['max diameter [mm]'], cmap='Greys',
                        edgecolor='black', lw=0.2, s=8)
        if annotate_some is not None:
            for id_ in annotate_some:
                ax.scatter(df['d90'].iloc[id_],
                           df['req. weight ks_p95 <= 10 [kg]'].iloc[id_],
                           edgecolor='black', facecolors='none', s=20, lw=1,
                           zorder=200)
                if id_ == 210:
                    t = ax.text(
                        x=df['d90'].iloc[id_],
                        y=df['req. weight ks_p95 <= 10 [kg]'].iloc[id_]+10,
                        s=df['ID'].iloc[id_], ha='left', fontsize=self.fsize,
                        weight='bold', zorder=100)
                else:
                    t = ax.text(
                        x=df['d90'].iloc[id_],
                        y=df['req. weight ks_p95 <= 10 [kg]'].iloc[id_]+10,
                        s=df['ID'].iloc[id_], ha='center', fontsize=self.fsize,
                        weight='bold', zorder=100)
                t.set_bbox(dict(facecolor='white', alpha=0.5, lw=0))

        ax.set_xlabel('d90 [mm]', fontsize=self.fsize)
        ax.set_ylabel('req. sample weight [kg]', fontsize=self.fsize)
        ax.grid(alpha=0.5)
        ax.set_ylim(bottom=0, top=400)

        cbar = fig.colorbar(im, cax=cax, orientation='vertical')

        cbar.set_label(label='max diameter [mm]', fontsize=self.fsize)
        ticks = cbar.ax.get_yticks()
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticks(ticks)
        cbar.ax.set_yticklabels(ticklabs, fontsize=self.fsize)
        ax.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def exponents_plot(self, df: pd.DataFrame, savepath: str,
                       close: bool = True) -> None:
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

        med_a, med_b = self.fit_exponential(exponents, med_errors)
        p95_a, p95_b = self.fit_exponential(exponents, p95_errors)

        # plotting
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1,
                                       figsize=(3.54331, 6.5))

        x = np.arange(200)
        colors = ['black'] * 4 + ['grey'] * 4
        print(colors)

        for i, exp in enumerate(list(results.keys())):
            if i % 2 == 0:  # only show every second function
                y = ((x / 10)**exp)
                ax1.plot(x, y, ls=self.lss[int(i/2)], color=colors[int(i/2)],
                         label=fr'$\epsilon$ {exp}; $KS_{{med}}$ {round(results[exp][0], 1)}, $KS_{{p95}}$ {round(results[exp][1], 1)}')

        ax1.grid(alpha=0.4)
        ax1.legend(fontsize=self.fsize)
        ax1.set_ylabel(r'$m_{min}$ [kg]', fontsize=self.fsize)
        ax1.set_xlabel('d90 [mm]', fontsize=self.fsize)
        ax1.tick_params(axis='both', labelsize=self.fsize)

        ax2.scatter(exponents, med_errors,  # label=r'$KS_{med}$',
                    color='black', alpha=0.7, s=12)
        label = r'$KS_{{med}}={}*e^{{({}*\epsilon)}}$'.format(round(med_a, 2), round(med_b, 2))
        ax2.plot(exponents,
                 self.exponential(np.array(exponents), med_a, med_b),
                 label=label, color='black')
        ax2.scatter(exponents, p95_errors,  # label=r'$KS_{p95}$',
                    color='black', alpha=0.7, s=12)
        label = r'$KS_{{p95}}={}*e^{{({}*\epsilon)}}$'.format(round(p95_a, 2), round(p95_b, 2))
        ax2.plot(exponents,
                 self.exponential(np.array(exponents), p95_a, p95_b),
                 label=label, color='black', ls='--')

        ax2.grid(alpha=0.5)
        ax2.legend(fontsize=self.fsize)
        ax2.set_xlabel(r'exponent ($\epsilon$)', fontsize=self.fsize)
        ax2.set_ylabel('Kolmogorov-Smirnov statisic [mass %]',
                       fontsize=self.fsize)
        ax2.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def comparison_plot(self, df, savepath: str, close: bool = True) -> None:
        '''plot that compares ISO and ASTM sample requirements with new ones'''

        fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2,
                                       figsize=(3.54331, 6.5))

        ax1.scatter(df['ISO req. weight [kg]'],
                    df['new. req. weight ISO [kg]'],
                    color='grey', alpha=0.4, s=10)
        ax1.plot([0, 400], [0, 400], color='black', lw=2, ls='--')
        ax1.text(400, 400, '1:1', weight='bold', fontsize=self.fsize,
                 ha='right')
        ax1.grid(alpha=0.5)
        ax1.set_xlabel('ISO req. weight [kg]', fontsize=self.fsize)
        ax1.set_ylabel('req. weight acc. new criterium [kg]',
                       fontsize=self.fsize)
        ax1.tick_params(axis='both', labelsize=self.fsize)

        ax2.scatter(df['ASTM req. weight [kg]'],
                    df['new. req. weight ASTM [kg]'],
                    color='grey', alpha=0.4, s=10)
        ax2.plot([0, 1200], [0, 1200], color='black', lw=2, ls='--')
        ax2.text(1200, 1200, '1:1', weight='bold', fontsize=self.fsize,
                 ha='right')
        ax2.grid(alpha=0.5)
        ax2.set_xlabel('ASTM req. weight [kg]', fontsize=self.fsize)
        ax2.set_ylabel('req. weight acc. new criterium [kg]',
                       fontsize=self.fsize)
        ax2.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def real_sieve_curves_plot(self, savepath: str,
                               close: bool = True) -> None:
        fp = r'../laboratory/LabResults.xlsx'
        df = pd.read_excel(fp, header=2, nrows=13, usecols=list(range(1, 9)))
        self.headers = list(df.columns)

        fig, ax = self.make_sieve_plot()

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
        plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def real_sieve_curves_scatter(self, savepath: str,
                                  close: bool = True) -> None:
        fp = r'../laboratory/LabResults.xlsx'
        df = pd.read_excel(fp, skiprows=20, nrows=2, usecols=list(range(1, 9)))
        self.headers[0] = 'parameter'
        df.columns = self.headers
        df.set_index('parameter', inplace=True)

        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1,
                                       figsize=(3.54331, 6.5))

        ms = 100
        min_, max_ = df.loc['Cc'].min(), df.loc['Cc'].max()
        min_, max_ = min_-min_*0.1, max_+max_*0.1

        ax1.scatter(df['Soil A (ISO)'].loc['Cc'],
                    df['Soil A (100g)'].loc['Cc'],
                    label='Soil A (100g)', s=ms, marker='o', color='C0')
        ax1.scatter(df['Soil A (ISO)'].loc['Cc'],
                    df['Soil A (75g)'].loc['Cc'],
                    label='Soil A (75g)', s=ms, marker='o', color='C1')
        ax1.scatter(df['Soil A (ISO)'].loc['Cc'],
                    df['Soil A (50g)'].loc['Cc'],
                    label='Soil A (50g)', s=ms, marker='o', color='C2')
        ax1.scatter(df['Soil A (ISO)'].loc['Cc'],
                    df['Soil A (5g)'].loc['Cc'],
                    label='Soil A (5g)', s=ms, marker='o', color='C3')

        ax1.scatter(df['Soil C (ISO)'].loc['Cc'], df['Soil C'].loc['Cc'],
                    label='Soil C (20kg)', s=ms, marker='P', color='C0')
        ax1.plot([min_, max_], [min_, max_], color='black')
        middle = min_ + (max_ - min_) / 2
        ax1.text(x=middle, y=middle, s='1 : 1 line', rotation=45, va='top',
                 ha='center', fontsize=self.fsize)

        ax1.set_xlim(left=min_, right=max_)
        ax1.set_ylim(bottom=min_, top=max_)
        ax1.set_aspect('equal')
        ax1.grid(alpha=0.5)
        ax1.set_title('Coefficient of Curvature', fontsize=self.fsize)
        ax1.set_xlabel('tests with ISO sample weight', fontsize=self.fsize)
        ax1.set_ylabel('tests with lower sample weight', fontsize=self.fsize)
        ax1.legend(fontsize=self.fsize)
        ax1.tick_params(axis='both', labelsize=self.fsize)

        min_, max_ = df.loc['Cu'].min(), df.loc['Cu'].max()
        min_, max_ = min_-min_*0.1, max_+max_*0.1

        ax2.scatter(df['Soil A (ISO)'].loc['Cu'],
                    df['Soil A (100g)'].loc['Cu'],
                    label='Soil A (100g)', s=ms, marker='o', color='C0')
        ax2.scatter(df['Soil A (ISO)'].loc['Cu'],
                    df['Soil A (75g)'].loc['Cu'],
                    label='Soil A (75g)', s=ms, marker='o', color='C1')
        ax2.scatter(df['Soil A (ISO)'].loc['Cu'],
                    df['Soil A (50g)'].loc['Cu'],
                    label='Soil A (50g)', s=ms, marker='o', color='C2')
        ax2.scatter(df['Soil A (ISO)'].loc['Cu'],
                    df['Soil A (5g)'].loc['Cu'],
                    label='Soil A (5g)', s=ms, marker='o', color='C3')

        ax2.scatter(df['Soil C (ISO)'].loc['Cu'], df['Soil C'].loc['Cu'],
                    label='Soil C (20kg)', s=ms, marker='P', color='C0')
        ax2.plot([min_, max_], [min_, max_], color='black')
        middle = min_ + (max_ - min_) / 2
        ax2.text(x=10, y=10, s='1 : 1 line', rotation=45, va='top',
                 ha='center', fontsize=self.fsize)

        ax2.set_xlim(left=min_, right=max_)
        ax2.set_ylim(bottom=min_, top=max_)
        ax2.set_aspect('equal')
        ax2.grid(alpha=0.5)
        ax2.set_title('Coefficient of Uniformity', fontsize=self.fsize)
        ax2.set_xlabel('tests with ISO sample weight', fontsize=self.fsize)
        ax2.set_ylabel('tests with lower sample weight', fontsize=self.fsize)
        ax2.legend(fontsize=self.fsize)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.tick_params(axis='both', labelsize=self.fsize)

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

    def make_sieve_plot(self, main_fontsize=12, iso_fontsize=10):
        '''function creates the background of a sieve plot ... does not plot
        data!'''
        # create sieve curves background
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set_xscale('log')
        ax.set_xlim(left=0.002, right=630)
        ax.set_ylim(bottom=0, top=101)
        ax.set_xticks([0.06, 2, 63, 200])
        ax.vlines([0.002, 0.06, 2, 63, 200], ymin=0, ymax=101,
                  color='lightgrey')
        ax.vlines([0.006, 0.02, 0.2, 0.6, 6.3, 20], ymin=0, ymax=101,
                  color='lightgrey', ls='--')
        ax.set_xlabel('grain size [mm]', fontsize=main_fontsize)
        ax.set_ylabel('mass percentage passing [%]', fontsize=main_fontsize)
        ax.grid(alpha=0.5)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        ax.tick_params(axis='both', labelsize=main_fontsize)

        # grain size labels
        ax.text(0.009, 118, 'SILT', fontsize=main_fontsize, va='top')
        ax.text(0.003, 112, 'fine', fontsize=main_fontsize, va='top')
        ax.text(0.007, 112, 'medium', fontsize=main_fontsize, va='top')
        ax.text(0.025, 112, 'coarse', fontsize=main_fontsize, va='top')
        ax.text(0.25, 118, 'SAND', fontsize=main_fontsize, va='top')
        ax.text(0.1, 112, 'fine', fontsize=main_fontsize, va='top')
        ax.text(0.22, 112, 'medium', fontsize=main_fontsize, va='top')
        ax.text(0.7, 112, 'coarse', fontsize=main_fontsize, va='top')
        ax.text(7, 118, 'GRAVEL', fontsize=main_fontsize, va='top')
        ax.text(3, 112, 'fine', fontsize=main_fontsize, va='top')
        ax.text(7.1, 112, 'medium', fontsize=main_fontsize, va='top')
        ax.text(25, 112, 'coarse', fontsize=main_fontsize, va='top')
        ax.text(70, 115, 'COBBEL', fontsize=main_fontsize, va='top')
        ax.text(220, 115, 'BLOCKS', fontsize=main_fontsize, va='top')
        # ISO lablling
        yiso = 101.5
        ax.text(0.005, yiso, 'ISO Standard sieves:', fontsize=iso_fontsize,
                va='bottom')
        ax.text(0.063, yiso, '.063', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(0.125, yiso, '.125', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(0.25, yiso, '.25', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(0.5, yiso, '.5', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(1, yiso, '1', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(2, yiso, '2', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(4, yiso, '4', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(8, yiso, '8', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(16, yiso, '16', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(31.5, yiso, '31.5', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(45, yiso, '45', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(63, yiso, '63', fontsize=iso_fontsize,
                va='bottom', ha='center')
        ax.text(90, yiso, '90', fontsize=iso_fontsize,
                va='bottom', ha='center')

        # plot grid around labels
        trans = ax.get_xaxis_transform()
        # horizontal lines
        ax.plot([0.002, 630], [1.06, 1.06], color="k", transform=trans,
                clip_on=False)
        ax.plot([0.002, 63], [1.12, 1.12], color="k", transform=trans,
                clip_on=False)
        ax.plot([0.002, 630], [1.18, 1.18], color="k", transform=trans,
                clip_on=False)
        # bold vertical lines
        ax.plot([0.002, 0.002], [1.06, 1.18], color="k", transform=trans,
                clip_on=False)
        ax.plot([0.063, 0.063], [1.06, 1.18], color="k", transform=trans,
                clip_on=False)
        ax.plot([2, 2], [1.06, 1.18], color="k", transform=trans,
                clip_on=False)
        ax.plot([63, 63], [1.06, 1.18], color="k", transform=trans,
                clip_on=False)
        ax.plot([200, 200], [1.06, 1.18], color="k", transform=trans,
                clip_on=False)
        ax.plot([630, 630], [1.06, 1.18], color="k", transform=trans,
                clip_on=False)
        # thin vertical lines
        ax.plot([0.0063, 0.0063], [1.06, 1.11], color="grey", linewidth=0.5,
                transform=trans, clip_on=False)
        ax.plot([0.63, 0.63], [1.06, 1.11], color="grey", linewidth=0.5,
                transform=trans, clip_on=False)
        ax.plot([6.3, 6.3], [1.06, 1.11], color="grey", linewidth=0.5,
                transform=trans, clip_on=False)
        ax.plot([0.02, 0.02], [1.06, 1.11], color="grey", linewidth=0.5,
                transform=trans, clip_on=False)
        ax.plot([0.2, 0.2], [1.06, 1.11], color="grey", linewidth=0.5,
                transform=trans, clip_on=False)
        ax.plot([20, 20], [1.06, 1.11], color="grey", linewidth=0.5,
                transform=trans, clip_on=False)

        return fig, ax

    def sieve_curves_plot(self, SIEVE_SIZES: list, fractions_true: list,
                          color: pd.Series = None,
                          savepath: str = None,
                          sieved_samples: list = None,
                          req_sample_weights: list = None,
                          ks_distances: list = None,
                          close: bool = True) -> None:
        '''plot sieve curves of underlying soil distribution and taken
        samples'''
        cmap = matplotlib.colormaps['inferno']

        fig, ax = self.make_sieve_plot()

        # plot data
        if len(np.array(fractions_true).shape) == 1:
            ax.plot(SIEVE_SIZES, fractions_true, label="underlying soil",
                    color='black', lw=3)
        else:
            for i, f in enumerate(fractions_true):
                if color is None:
                    ax.plot(SIEVE_SIZES, f, color='black', alpha=0.2)
                else:
                    c = cmap(int(color[i] / color[:len(fractions_true)].max() * 255))
                    ax.plot(SIEVE_SIZES, f, color=c)

        if sieved_samples is not None:
            for i in range(len(sieved_samples)):
                ax.plot(
                    SIEVE_SIZES, np.array(sieved_samples[i]),
                    label=f"sample {round(req_sample_weights[i], 1)}kg  ks: {round(ks_distances[i], 1)} %",
                    alpha=0.8)

        ax.legend(loc='upper left', fontsize=12)

        plt.tight_layout()
        if savepath is not None:
            plt.savefig(savepath, dpi=600)
        if close is True:
            plt.close()

    def monte_carlo_scatterplot_legacy(self, df: pd.DataFrame, weight_mode: str,
                                color_mode: str, savepath: str,
                                close: bool = True) -> None:
        # old plot that is not required anymore -> TODO remove
        '''scatterplot that shows results of the monte carlo simulations in the
        form of different soil distribution parameters Cu, Cc etc. vs. a metric
        of how well the sample fits the underlying distribution'''

        ks_median = df[f'kolmogorov smirnov distance {weight_mode}'].median()
        ks_p95 = np.percentile(df[f'kolmogorov smirnov distance {weight_mode}'], 95)

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

                ax3.axhline(y=ks_median, color='black', ls='-',
                            label=f'med. error {round(ks_median, 1)}, p95 {round(ks_p95, 1)}')
                ax3.axhline(y=ks_p95, color='black', ls='--')

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


class sample_preview:
    '''class with functions to make a preview of a soil sample - very
    experimental, use not recommended'''

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
        delta_x = np.random.normal(loc=0, scale=grain_rs/10,
                                   size=len(grain_xs))
        delta_y = np.random.normal(loc=0, scale=grain_rs/10,
                                   size=len(grain_ys))
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
    TOT_MASS = 60
    SIEVE_SIZES = np.exp(np.linspace(np.log(1), np.log(200), 30))

    lab, pltr = laboratory(), plotter()

    # make plot of random sieve curves to demonstrate sampler
    fractions_trues, max_diameters = [], []

    for i in range(400):
        grain_diameters, grain_weights, grain_ids = lab.make_grains(
            DENSITY, TOT_MASS=TOT_MASS)
        Dmax = grain_diameters.max()
        print(i, grain_weights.sum())
        fractions_true = lab.sieve(grain_ids, grain_diameters, grain_weights,
                                   SIEVE_SIZES)
        fractions_trues.append(list(fractions_true.values()))
        max_diameters.append(Dmax)

    pltr.sieve_curves_plot(list(fractions_true.keys()), fractions_trues,
                           savepath=fr'../figures/sieve_samples_new.jpg',
                           close=False)

    fig, ax = plt.subplots()
    ax.hist(max_diameters, color='grey', edgecolor='black', bins=30)
