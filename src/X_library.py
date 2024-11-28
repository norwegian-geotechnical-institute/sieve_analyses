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

# importing libraries
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from tqdm import tqdm


class statistics:
    '''class with use case specfic statistical functions'''

    def ks_statistic(self, sieve_curve1: list, sieve_curve2: list) -> float:
        '''calculate Kolmogorov-Smirnov (KS) statistic for 2 sieve curves'''
        # Calculate the maximum vertical distance between the sieve curves
        return np.max(np.abs(np.array(sieve_curve1) - np.array(sieve_curve2)))


class laboratory(statistics):
    '''class with functions for simulating laboratory sieve tests'''

    def ISO_required_sample_mass(self, d_max: float) -> float:
        '''calculate required sample mass acc. to ISO 17892-4'''
        if d_max <= 2:  # [mm]
            mass = 100
        elif d_max <= 6.3:
            mass = 300
        elif d_max <= 10:
            mass = 500
        elif d_max <= 20:
            mass = 2000
        else:
            mass = ((d_max / 10)**2) * 1000
        return mass  # [g]

    def ASTM_required_sample_mass(self, d_max: float) -> float:
        '''calculate required sample mass acc. to ASTM D6913/D6913M,
        method A'''
        density = 3  # [g/mm3] ... ASTM seems to use this high grain density

        if d_max <= 4.75:
            mass = 50
        elif d_max <= 9.5:
            mass = 75
        elif d_max <= 76.2:
            v = (4/3)*np.pi*(d_max/2)**3  # [mm3]
            v = v / 1000  # [cm3]
            m = v * density / 1000  # [kg]
            factor = 1.2
            mass = m * 100 * factor * 1000
        else:
            v = (4/3)*np.pi*(d_max/2)**3  # [mm3]
            v = v / 1000  # [cm3]
            m = v * density / 1000  # [kg]
            mass = m * 100 * 1000
        return mass  # [g]

    def new_sample_mass(self, d_max: float, d90: float,
                        exponent: float = 2) -> float:
        '''new equation to determine sample mass based on d90'''
        if d_max <= 2:  # [mm]
            mass = 100
        elif d_max <= 6.3:
            mass = 300
        elif d_max <= 10:
            mass = 500
        elif d_max <= 20:
            mass = 2000
        else:
            mass = ((d90 / 10)**exponent) * 1000
        return mass  # [g]

    def get_sample(self, required_sample_mass: float, total_mass: float,
                   grain_masses: np.array, grain_diameters: np.array,
                   strategy: str = 'random choice') -> list:
        '''get a sample with a specific mass from the total number of grains
        by gathering single grains until the desired mass is reached'''

        fraction = required_sample_mass / total_mass
        match strategy:  # noqa
            case 'from start':
                split = int(len(grain_diameters)*fraction)
                sample_ids = np.arange(split)
                sample_masses = grain_masses[:split]
                sample_diameters = grain_diameters[:split]
            case 'random choice':
                n_grains = int(len(grain_diameters)*fraction)
                sample_ids = np.random.choice(np.arange(len(grain_diameters)),
                                              size=n_grains, replace=False)
                sample_masses = grain_masses[sample_ids]
                sample_diameters = grain_diameters[sample_ids]
        return sample_ids, sample_masses, sample_diameters

    def make_grains_from_distribution(self, DENSITY: float,  # [g/cm3]
                                      TOT_MASS: float,  # [kg]
                                      distribution: str) -> list:
        '''Function generates a new soil sample, based on a standard
        statistical distribution.'''

        diameters = []  # [mm]
        tot_mass = 0  # [kg]

        match distribution:
            case 'normal':
                while tot_mass < TOT_MASS:
                    diameter = np.random.normal(loc=32.5, scale=15,
                                                size=2)  # [mm]
                    diameter = np.where(diameter < 0, diameter * -1, diameter)
                    volume = ((4/3)*np.pi*((diameter/2)**3)) / 1000  # [cm3]
                    tot_mass += volume.sum() * DENSITY / 1000  # [kg]
                    diameters.append(diameter)
            case _:
                raise ValueError(f'{distribution} not implemented')

        volumes = (4/3)*np.pi*(np.array(diameters)/2)**3  # [mm3]
        volumes = volumes / 1000  # [cm3]
        masses = volumes * DENSITY / 1000  # [kg]
        ids = np.arange(len(diameters))
        return diameters, masses, ids

    def make_grains(self, DENSITY: float, TOT_MASS: float,
                    min_d: float = 1, max_d: float = 200,
                    verbose: bool = False) -> list:
        '''Function generates a new soil sample, based on generation process
        described in the paper that mimics real soils that are typically
        composed of multiple soil fractions due to their geological history.'''

        def make_beta(lower: float, upper: float, n: int) -> np.array:
            '''function generates an array of grain sizes based on a beta
            distribution and scaled between lower and upper [mm]'''
            x = np.random.beta(a=np.random.uniform(1, 5),
                               b=np.random.uniform(1, 5), size=n)
            x = x*(upper-lower)
            return x + lower

        fractions = np.random.uniform(0, 1, np.random.randint(1, 6))
        fractions = fractions / sum(fractions)
        fractions = np.sort(fractions)

        diameters = []  # [mm]
        tot_mass = 0  # [kg]

        for i in range(len(fractions)):
            # lower, upper = size boundaries of soil fractions in [mm]
            lower, upper = sorted(np.exp(np.random.uniform(np.log(min_d),
                                                           np.log(max_d), 2)))
            while tot_mass < sum(fractions[:i+1]) * TOT_MASS:
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
                    tot_mass += volume.sum() * DENSITY / 1000
                except AttributeError:
                    tot_mass += volume * DENSITY / 1000
                diameters.append(diameter)

        diameters = np.hstack(diameters)
        np.random.shuffle(diameters)
        volumes = (4/3)*np.pi*(diameters/2)**3  # [mm3]
        volumes = volumes / 1000  # [cm3]
        masses = volumes * DENSITY / 1000  # [kg]
        ids = np.arange(len(diameters))
        return diameters, masses, ids

    def sieve(self, diameters: np.array, masses: np.array,
              sieve_sizes: list) -> dict:
        '''make a virtual sieve analysis of the sample'''
        tot_mass = masses.sum()
        fractions = []
        for size in sieve_sizes:
            fraction_ids = np.where(diameters < size)[0]
            fraction_mass = masses[fraction_ids].sum()
            try:
                fraction = (fraction_mass / tot_mass) * 100
            except ZeroDivisionError:
                fraction = 0
            fractions.append(fraction)
        return dict(zip(sieve_sizes, fractions))

    def calc_grading_characteristics_legacy(self, fractions_true: list,
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

    def calc_grading_characteristics(self, grain_diameters: np.array,
                                     grain_masses: np.array) -> list:
        '''function computes different metrics about the grading of a soil
        based on the real cumulative density curve of the soil mass - not based
        on the sieve curve, as this introduces error esp. in large grain due to
        mesh widths sizes'''
        d_s = [10, 12, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90]
        d_keys = [f'd{d}' for d in d_s]

        id_sorted = np.argsort(grain_diameters)
        sorted_grain_diameters = grain_diameters[id_sorted]
        cumulative_relative_grain_masses = np.cumsum(grain_masses[id_sorted])/sum(grain_masses)*100
        d_vals = np.interp(d_s, cumulative_relative_grain_masses,
                           sorted_grain_diameters)
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

    def check_required_mass(self, max_error: float,  # [%]
                            mass_steps: float,  # steps to increase mass
                            tests_per_step: int, total_mass: float,
                            grain_masses: np.array,
                            grain_diameters: np.array, SIEVE_SIZES: list,
                            fractions_true: dict,
                            verbose: bool = False):
        '''function computes the theoretically required mass to achieve a
        defined error by gathering sequentially more soil. Criterium is p95
        percentile of errors below max_error'''
        ks_p95, mass = 2*max_error, 0  # initialize parameters
        while ks_p95 > max_error:
            mass += mass_steps
            masses_temp, ks_s_temp = [], []
            # make multiple tests to reduce randomnes
            for _ in range(tests_per_step):
                sample_ids, sample_masses, sample_diameter = self.get_sample(
                    mass, total_mass, grain_masses, grain_diameters,
                    strategy='random choice')
                fractions = self.sieve(sample_diameter, sample_masses,
                                       SIEVE_SIZES)
                masses_temp.append(sample_masses.sum())
                ks_s_temp.append(self.ks_statistic(
                    list(fractions_true.values()), list(fractions.values())))
            mass_avg = np.mean(masses_temp)
            ks_p95 = np.percentile(ks_s_temp, 95)
            if verbose is True:
                print(f'{round(mass_avg, 1)} kg\t\tp95 {round(ks_p95, 1)}')

        return mass  # required mass [kg]


class utilities:
    '''class with general purpose functions that do not fit anywhere else'''

    def CheckForLess(self, list1, val):
        '''function checks if any of values in list1 is smaller than val'''
        for x in list1:
            if val < x:
                return False
        return True

    def exponential(self, x: np.array, a: float, b: float) -> np.array:
        '''exponential function'''
        return a * np.exp(b*x)

    def fit_exponential(self, x: np.array, y: np.array) -> tuple:
        '''functions fits an exponential function to x and y data and yields
        fitted parameters'''
        p0 = [100, -1]  # initial guess
        params, _ = curve_fit(self.exponential, x, y, p0=p0)  # curve fitting
        a_fit, b_fit = params  # Extract the fitted parameters
        return a_fit, b_fit


class plotter(laboratory, utilities):
    '''class with functions for data and result plotting'''

    def __init__(self):
        self.fsize = 8  # main font size for plots
        self.lss = ['-', '--', ':', '-.'] * 2  # choice of linestyles

    def required_mass_plot(self, max_grain_size: float, savepath: str,
                           save_pdf: bool = False, close: bool = True) -> None:
        '''plot that shows the theoretically required sample mass acc. to the
        standards'''
        sizes = np.arange(max_grain_size)  # [mm]
        ms_ISO = [self.ISO_required_sample_mass(ds)/1000 for ds in sizes]
        ms_ASTM = [self.ASTM_required_sample_mass(ds)/1000 for ds in sizes]

        fig, ax = plt.subplots(figsize=(3.54331, 3.54331))
        ax.plot(sizes, ms_ISO, color='black', ls='--',
                label='ISO 17892-4', zorder=10)
        ax.plot(sizes, ms_ASTM, color='black', ls='-',
                label='ASTM D6913/D6913M', zorder=10)
        vlines = (20, 63, 200)
        ax.vlines(vlines, ymin=0, ymax=max(ms_ASTM), color='grey',
                  zorder=5)
        for vline in vlines:
            ax.text(x=vline+2, y=max(ms_ASTM), s=f'{vline}\nmm', ha='left',
                    va='top', fontsize=self.fsize)
        ax.grid(alpha=0.5, axis='y')
        ax.set_xlabel('max. grain diameter\n$D_{max}$ [mm]', fontsize=self.fsize)
        ax.set_ylabel('required sample mass\n$m_{min}$ [kg]', fontsize=self.fsize)
        ax.set_xlim(left=2)
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
        ax.legend(loc='lower right', framealpha=1, fontsize=self.fsize)
        ax.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if save_pdf is True:
            plt.savefig(savepath[:-3]+'pdf')
        if close is True:
            plt.close()

    def error_violin_plot(self, df: pd.DataFrame, savepath: str,
                          save_pdf: bool = False, close: bool = True) -> None:
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
        ax.set_xticks([1, 2], labels=['ISO 17892-4', 'ASTM D6913/D6913M'],
                      fontsize=self.fsize)
        ax.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if save_pdf is True:
            plt.savefig(savepath[:-3]+'pdf')
        if close is True:
            plt.close()

    def req_sample_mass_vs_dmax_plot(self, df: pd.DataFrame, savepath: str,
                                     save_pdf: bool = False,
                                     annotate_all=False, annotate_some=None,
                                     close: bool = True) -> None:
        '''plot shows results from Monte-Carlo simulation where the required
        sample mass to achieve a certain KS error is scattered against the soil
        d_max'''

        def add_scatter_axis():
            '''add axis to the figure that contains the scatterplot of dmax vs
            required sample mass'''
            ax = plt.gca()

            dmaxs = np.arange(200)
            req_ISO = [self.ISO_required_sample_mass(dmax)/1000 for dmax in dmaxs]
            req_ASTM = [self.ASTM_required_sample_mass(dmax)/1000 for dmax in dmaxs]
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)

            im = ax.scatter(df['max diameter [mm]'],
                            df['req. mass ks_p95 <= 10 [kg]'],
                            c=df['S0'], cmap='Greys', edgecolor='black',
                            lw=0.2, label='simulated samples', s=8, vmax=9)

            ax.plot(dmaxs, req_ISO, color='dimgrey', ls='--',
                    label='ISO req. sample mass')
            ax.plot(dmaxs, req_ASTM, color='dimgrey', ls='-',
                    label='ASTM req. sample mass')
            if annotate_all is True:
                for i in range(len(df)):
                    ax.text(x=df['max diameter [mm]'].iloc[i],
                            y=df['req. mass ks_p95 <= 10 [kg]'].iloc[i],
                            s=df['ID'].iloc[i])

            ax.set_xlabel('max. grain diameter\n$D_{max}$ [mm]',
                          fontsize=self.fsize)
            ax.set_ylabel('required sample mass\n$m_{min}$ [kg]',
                          fontsize=self.fsize)
            ax.grid(alpha=0.5)
            ax.set_ylim(bottom=0, top=400)
            ax.legend(fontsize=self.fsize, loc='upper left')

            cbar = fig.colorbar(im, cax=cax, orientation='vertical')

            cbar.set_label(label=r'sorting coefficient - $S_0$',
                           fontsize=self.fsize)
            # ticks = cbar.ax.get_yticks()
            # ticklabs = cbar.ax.get_yticklabels()
            # cbar.ax.set_yticks(ticks)
            # cbar.ax.set_yticklabels(ticklabs, fontsize=self.fsize)
            ax.tick_params(axis='both', labelsize=self.fsize)

        def add_simple_sieve_plot_axis():
            '''simple sieve curves plot to combine with
            req_sample_mass_vs_dmax_plot'''

            ax = plt.gca()

            mpercs = [c for c in df.columns if ' mm sieve [m%]' in c]
            ssizes = [eval(s.split(' ')[0]) for s in mpercs]

            for i, id_ in enumerate(annotate_some):
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

        fig = plt.figure(figsize=(3.54331, 3))

        if annotate_some is None:
            ax = fig.add_subplot(111)
            add_scatter_axis()
        else:
            fig.set_size_inches(3.54331, 6)

            ax = fig.add_subplot(211)
            add_scatter_axis()
            for id_ in annotate_some:
                ax.scatter(df['max diameter [mm]'].iloc[id_],
                           df['req. mass ks_p95 <= 10 [kg]'].iloc[id_],
                           edgecolor='black', facecolors='none', s=20, lw=1,
                           zorder=200)
                t = ax.text(x=df['max diameter [mm]'].iloc[id_],
                            y=df['req. mass ks_p95 <= 10 [kg]'].iloc[id_]+10,
                            s=df['ID'].iloc[id_], ha='center',
                            fontsize=self.fsize, weight='bold', zorder=100)
                t.set_bbox(dict(facecolor='white', alpha=0.5, lw=0))

            ax = fig.add_subplot(212)
            add_simple_sieve_plot_axis()

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if save_pdf is True:
            plt.savefig(savepath[:-3]+'pdf')
        if close is True:
            plt.close()

    def req_sample_mass_vs_d90_plot(self, df, savepath: str,
                                    annotate_some=None, save_pdf: bool = False,
                                    close: bool = True) -> None:
        '''plot shows results from Monte-Carlo simulation where the required
        sample mass to achieve a certain KS error is scattered against the soil
        d90'''
        fig, ax = plt.subplots(figsize=(3.54331, 3))

        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        im = ax.scatter(df['d90'], df['req. mass ks_p95 <= 10 [kg]'],
                        c=df['max diameter [mm]'], cmap='Greys',
                        edgecolor='black', lw=0.2, s=8)
        if annotate_some is not None:
            for id_ in annotate_some:
                ax.scatter(df['d90'].iloc[id_],
                           df['req. mass ks_p95 <= 10 [kg]'].iloc[id_],
                           edgecolor='black', facecolors='none', s=20, lw=1,
                           zorder=200)
                if id_ == 210:
                    t = ax.text(
                        x=df['d90'].iloc[id_],
                        y=df['req. mass ks_p95 <= 10 [kg]'].iloc[id_]+10,
                        s=df['ID'].iloc[id_], ha='left', fontsize=self.fsize,
                        weight='bold', zorder=100)
                else:
                    t = ax.text(
                        x=df['d90'].iloc[id_],
                        y=df['req. mass ks_p95 <= 10 [kg]'].iloc[id_]+10,
                        s=df['ID'].iloc[id_], ha='center', fontsize=self.fsize,
                        weight='bold', zorder=100)
                t.set_bbox(dict(facecolor='white', alpha=0.5, lw=0))

        ax.set_xlabel('$D_{90}$ [mm]', fontsize=self.fsize)
        ax.set_ylabel('required sample mass\n$m_{min}$ [kg]', fontsize=self.fsize)
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
        if save_pdf is True:
            plt.savefig(savepath[:-3]+'pdf')
        if close is True:
            plt.close()

    def exponents_plot(self, df: pd.DataFrame, savepath: str,
                       save_pdf: bool = False, close: bool = True) -> None:
        '''plot that shows the different error exponents in relation to
        parameters like minimum required smaple mass, d90, and KS error'''
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

        for i, exp in enumerate(list(results.keys())):
            if i % 2 == 0:  # only show every second function
                y = ((x / 10)**exp)
                ax1.plot(x, y, ls=self.lss[int(i/2)], color=colors[int(i/2)],
                         label=fr'$\epsilon$ {exp}; $KS_{{med}}$ {round(results[exp][0], 1)}, $KS_{{p95}}$ {round(results[exp][1], 1)}')

        ax1.grid(alpha=0.4)
        ax1.legend(fontsize=self.fsize)
        ax1.set_ylabel('required sample mass\n$m_{min}$ [kg]', fontsize=self.fsize)
        ax1.set_xlabel('$D_{90}$ [mm]', fontsize=self.fsize)
        ax1.tick_params(axis='both', labelsize=self.fsize)

        ax2.scatter(exponents, med_errors,  # label=r'$KS_{med}$',
                    color='black', alpha=0.7, s=12)
        label = r'$KS_{{med}}={}*e^{{({}*\epsilon)}}$'.format(round(med_a, 2),
                                                              round(med_b, 2))
        ax2.plot(exponents,
                 self.exponential(np.array(exponents), med_a, med_b),
                 label=label, color='black')
        ax2.scatter(exponents, p95_errors,  # label=r'$KS_{p95}$',
                    color='black', alpha=0.7, s=12)
        label = r'$KS_{{p95}}={}*e^{{({}*\epsilon)}}$'.format(round(p95_a, 2),
                                                              round(p95_b, 2))
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
        if save_pdf is True:
            plt.savefig(savepath[:-3]+'pdf')
        if close is True:
            plt.close()

    def comparison_plot(self, df, savepath: str, save_pdf: bool = False,
                        close: bool = True) -> None:
        '''plot that compares ISO and ASTM sample requirements with new ones'''

        fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2,
                                       figsize=(3.54331, 6.5))

        ax1.scatter(df['ISO req. mass [kg]'],
                    df['new. req. mass ISO [kg]'],
                    color='grey', alpha=0.4, s=10)
        ax1.plot([0, 400], [0, 400], color='black', lw=2, ls='--')
        ax1.text(400, 400, '1:1', weight='bold', fontsize=self.fsize,
                 ha='right')
        ax1.grid(alpha=0.5)
        ax1.set_xlabel('$m_{min}$ acc. ISO [kg]', fontsize=self.fsize)
        ax1.set_ylabel('$m_{min}$ acc. new criterium [kg]',
                       fontsize=self.fsize)
        ax1.tick_params(axis='both', labelsize=self.fsize)

        ax2.scatter(df['ASTM req. mass [kg]'],
                    df['new. req. mass ASTM [kg]'],
                    color='grey', alpha=0.4, s=10)
        ax2.plot([0, 1200], [0, 1200], color='black', lw=2, ls='--')
        ax2.text(1200, 1200, '1:1', weight='bold', fontsize=self.fsize,
                 ha='right')
        ax2.grid(alpha=0.5)
        ax2.set_xlabel('$m_{min}$ acc. ASTM [kg]', fontsize=self.fsize)
        ax2.set_ylabel('$m_{min}$ acc. new criterium [kg]',
                       fontsize=self.fsize)
        ax2.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if save_pdf is True:
            plt.savefig(savepath[:-3]+'pdf')
        if close is True:
            plt.close()

    def real_sieve_curves_plot(self, savepath: str, save_pdf: bool = False,
                               close: bool = True) -> None:
        '''function that plots the real laboratory sieve tests based on an
        excel file'''
        fp = r'../laboratory/LabResults.xlsx'
        df = pd.read_excel(fp, header=2, nrows=13, usecols=list(range(1, 12)))
        self.headers = list(df.columns)

        fig, ax = self.make_sieve_plot()

        ax.plot(df['Sieve size Ø [mm]'], df['Soil A (ISO)'],
                lw=3, color='C0', label='Soil A (ISO), 200 g')
        ax.plot(df['Sieve size Ø [mm]'], df['Soil A (5g)'],
                lw=1.5, color='C0', alpha=0.5, label='Soil A, 5 g')

        ax.plot(df['Sieve size Ø [mm]'], df['Soil B (ISO)'],
                lw=3, color='C1', label='Soil B (ISO), 9 kg')
        ax.plot(df['Sieve size Ø [mm]'], df['Soil B (1000g)'],
                lw=1.5, color='C1', alpha=0.5, label='Soil B, 1 kg')
        ax.plot(df['Sieve size Ø [mm]'], df['Soil B (300g)'],
                lw=1.5, color='C1', ls='--', alpha=0.5, label='Soil B, 0.3 kg')

        ax.plot(df['Sieve size Ø [mm]'], df['Soil C (ISO)'],
                lw=3, color='C2', label='Soil C (ISO), 50 kg')
        ax.plot(df['Sieve size Ø [mm]'], df['Soil C'],
                lw=1.5, color='C2', alpha=0.5, label='Soil C, 20 kg')

        ax.legend(loc='upper left')

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if save_pdf is True:
            plt.savefig(savepath[:-3]+'pdf')
        if close is True:
            plt.close()

    def real_sieve_curves_scatter(self, savepath: str, save_pdf: bool = False,
                                  close: bool = True) -> None:
        '''function that creates a scatter plot that scatters Cc and Cu from
        ISO sample masses against smaller sample masses'''
        fp = r'../laboratory/LabResults.xlsx'
        df = pd.read_excel(fp, skiprows=20, nrows=2,
                           usecols=list(range(1, 12)))
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

        ax1.scatter(df['Soil B (ISO)'].loc['Cc'],
                    df['Soil B (1000g)'].loc['Cc'],
                    label='Soil B (1 kg)', s=ms, marker='v', color='C0')
        ax1.scatter(df['Soil B (ISO)'].loc['Cc'],
                    df['Soil B (300g)'].loc['Cc'],
                    label='Soil B (0.3 kg)', s=ms, marker='v', color='C1')

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
        ax1.set_xlabel('tests with ISO sample mass', fontsize=self.fsize)
        ax1.set_ylabel('tests with lower sample mass', fontsize=self.fsize)
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

        ax2.scatter(df['Soil B (ISO)'].loc['Cu'],
                    df['Soil B (1000g)'].loc['Cu'],
                    label='Soil B (1 kg)', s=ms, marker='v', color='C0')
        ax2.scatter(df['Soil B (ISO)'].loc['Cu'],
                    df['Soil B (300g)'].loc['Cu'],
                    label='Soil B (0.3 kg)', s=ms, marker='v', color='C1')

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
        ax2.set_xlabel('tests with ISO sample mass', fontsize=self.fsize)
        ax2.set_ylabel('tests with lower sample mass', fontsize=self.fsize)
        ax2.legend(fontsize=self.fsize)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.tick_params(axis='both', labelsize=self.fsize)

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if save_pdf is True:
            plt.savefig(savepath[:-3]+'pdf')
        if close is True:
            plt.close()

    def distances_plot(self, req_sample_masses: list, ks_distances: list,
                       savepath: str, save_pdf: bool = False,
                       close: bool = True) -> None:
        '''plot a soil distribution vs. distributions of different samples with
        reduced mass and also their kolmogorov smirnov distance'''
        fig, ax = plt.subplots(figsize=(6, 6), nrows=1, ncols=1)

        ax.scatter(req_sample_masses, ks_distances, color='C0',
                   edgecolor='black', s=60)
        ax.grid(alpha=0.5)
        ax.set_xlabel('sample mass [kg]')
        ax.set_ylabel('kolmogorov smirnov distance', color='C0')

        plt.tight_layout()
        plt.savefig(savepath, dpi=600)
        if save_pdf is True:
            plt.savefig(savepath[:-3]+'pdf')
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

        plt.tight_layout()

        return fig, ax

    def sieve_curves_plot(self, SIEVE_SIZES: list, fractions_true: list,
                          color: pd.Series = None,
                          savepath: str = None,
                          sieved_samples: list = None,
                          req_sample_masses: list = None,
                          ks_distances: list = None,
                          show_dxxs: dict = None,
                          save_pdf: bool = False,
                          close: bool = True) -> None:
        '''plot sieve curves of soil distributions and taken samples'''
        cmap = matplotlib.colormaps['inferno']

        fig, ax = self.make_sieve_plot()

        # plot data
        if len(np.array(fractions_true).shape) == 1:
            ax.plot(SIEVE_SIZES, fractions_true, label="soil", color='black',
                    lw=3)
            ax.legend(loc='upper left', fontsize=12)
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
                    label=f"sample {round(req_sample_masses[i], 1)} kg  $KS$: {round(ks_distances[i], 1)} %",
                    alpha=0.8)
            ax.legend(loc='upper left', fontsize=12)

        # if show_dxxs is not None:
        #     for d_x in show_dxxs.keys():
        #         d = 
        #         ax.plot([])

        plt.tight_layout()
        if savepath is not None:
            plt.savefig(savepath, dpi=600)
            if save_pdf is True:
                plt.savefig(savepath[:-3]+'pdf')
        if close is True:
            plt.close()


if __name__ == '__main__':
    DENSITY = 2.65
    TOT_MASS = 60
    SIEVE_SIZES = np.exp(np.linspace(np.log(1), np.log(200), 30))

    lab, pltr = laboratory(), plotter()

    diameters, masses, ids = lab.make_grains_from_distribution(
        DENSITY, TOT_MASS, distribution='normal')
    fractions_true = lab.sieve(ids, diameters, masses, SIEVE_SIZES)

    pltr.sieve_curves_plot(SIEVE_SIZES, list(fractions_true.values()),
                           close=False)
