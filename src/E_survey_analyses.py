# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 11:55:18 2024

@author: GEr
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# load ground truth data
fp_gt = r'C:\Users\GEr\OneDrive - NGI\Research\Internal_Funding\GBV_GrainSizes\sieve_analyses\samples\analyses\Sample_selection_map.xlsx'
df_gt = pd.read_excel(fp_gt)

# load data from survey
fp_s = r'C:\Users\GEr\Downloads\Particle size distribution characterization survey(1-36).xlsx'
df_s = pd.read_excel(fp_s)

for c in df_s.columns:
    if ' ... ' in c:
        if '\n4' in c:
            c_new = f'S3_{c}'
        elif '\n3' in c:
            c_new = f'S2_{c}'
        elif '\n2' in c:
            c_new = f'S1_{c}'
        else:
            c_new = f'S0_{c}'
        c_new = c_new.split(' ... ')[0]
        df_s.rename({c: c_new}, axis=1, inplace=True)

df_s.drop(['ID', 'Start time', 'Completion time', 'Name2',
           'Email address (to send you results after the survey)\n'],
          axis=1, inplace=True)

# for c in df_s.columns:
#     print(c, df_s[c].dtype)

# compute Cu and Cc for suggestions
for sample in [0, 1, 2, 3]:
    df_s[f'S{sample}_Cu'] = df_s[f'S{sample}_d60'] / df_s[f'S{sample}_d10']
    df_s[f'S{sample}_Cc'] = (df_s[f'S{sample}_d30']**2) / (df_s[f'S{sample}_d60'] * df_s[f'S{sample}_d10'])
average_guess_Cu = df_s[['S0_Cu', 'S1_Cu', 'S2_Cu', 'S3_Cu']].mean(axis=0)
average_guess_Cc = df_s[['S0_Cc', 'S1_Cc', 'S2_Cc', 'S3_Cc']].mean(axis=0)

# correct labels in participant data
df_s['What is your main area of expertise?\n'] = df_s['What is your main area of expertise?\n'].replace(
    {'Engineering geology': 'Engineering\ngeology',
     'Quaternary geology': 'Quaternary\ngeology',
     'geophysics': 'Geophysics',
     'Digital engineering applied to geotechnics and strctures': 'Geotechnics',
     'Geotechnical engineering': 'Geotechnics',
     'Earthquake geology & sedimentology': 'Sedimentology',
     'Planetology / Astronautics': 'Other',
     'Geomechanics': 'Other'})

df_s['What is currently your main field of work?\n'] = df_s['What is currently your main field of work?\n'].replace(
    {'Industry (consulting, contractors, technology development,...)': 'Industry'})

df_s['How many years of experience post master do you have?'] = df_s['How many years of experience post master do you have?'].replace(
    {'None (still student or not from this field)': 'None'})

###############################
# plotting of survey results
###############################

# plot to show participant statistics
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(3.54331, 8))
colors = ['white', 'lightgrey', 'silver', 'darkgrey', 'grey'] * 10

professions, p_counts = np.unique(
    df_s['What is your main area of expertise?\n'], return_counts=True)
sort_ids = np.flip(np.argsort(p_counts))
axs[0].pie(p_counts[sort_ids], labels=professions[sort_ids], autopct='%1.0f%%',
           pctdistance=0.8, startangle=90, colors=colors,
           wedgeprops={'edgecolor': 'black'})
axs[0].set_title('Main area of expertise')

fields, f_counts = np.unique(
    df_s['What is currently your main field of work?\n'], return_counts=True)
sort_ids = [1, 0, 2]
axs[1].pie(f_counts[sort_ids], labels=fields[sort_ids], autopct='%1.0f%%',
           pctdistance=0.8, startangle=90, colors=colors,
           wedgeprops={'edgecolor': 'black'})
axs[1].set_title('Main field of work')

expertise, e_counts = np.unique(
    df_s['How many years of experience post master do you have?'],
    return_counts=True)
sort_ids = [5, 0, 3, 1, 2, 4]
axs[2].pie(e_counts[sort_ids], labels=expertise[sort_ids], autopct='%1.0f%%',
           pctdistance=0.8, startangle=90, colors=colors,
           wedgeprops={'edgecolor': 'black'})
axs[2].set_title('Years of experience post master')

plt.tight_layout()
plt.savefig(r'../samples/analyses/result_participants.svg')
plt.savefig(r'../samples/analyses/result_participants.pdf')
plt.close()

# plots to show results for Cu, Cc
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(3.54331, 7))
x = np.arange(1, 5)

axs[0].scatter(x, df_gt['Cu'], zorder=10, marker='P', s=80, color='black',
               edgecolor='white', label='true value')
axs[0].scatter(x, average_guess_Cu, zorder=10, marker='o', color='black',
               s=80, edgecolor='white', label='average estimations')

parts = axs[0].violinplot([df_s['S0_Cu'], df_s['S1_Cu'], df_s['S2_Cu'],
                           df_s['S3_Cu']],
                          showextrema=False, showmeans=False, points=10,
                          widths=0.7)
for pc in parts['bodies']:
    pc.set_facecolor('grey')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

axs[0].legend()
axs[0].grid(alpha=0.5)
axs[0].set_ylabel('$C_u$')
axs[0].set_xticks(x)
axs[0].set_xticklabels(['sample 1', 'sample 2', 'sample 3', 'sample 4'])


axs[1].scatter(x, df_gt['Cc'], zorder=10, marker='P', s=80, color='black',
               edgecolor='white', label='true value')
axs[1].scatter(x, average_guess_Cc, zorder=10, marker='o', color='black',
               s=80, edgecolor='white', label='average estimations')

parts = axs[1].violinplot([df_s['S0_Cc'], df_s['S1_Cc'], df_s['S2_Cc'],
                           df_s['S3_Cc']],
                          showextrema=False, showmeans=False, points=50,
                          widths=0.7)
for pc in parts['bodies']:
    pc.set_facecolor('grey')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

axs[1].legend()
axs[1].grid(alpha=0.5)
axs[1].set_ylabel('$C_c$')
axs[1].set_xticks(x)
axs[1].set_xticklabels(['sample 1', 'sample 2', 'sample 3', 'sample 4'])

plt.tight_layout()
plt.savefig(r'../samples/analyses/result_Cc_Cu.svg')
plt.savefig(r'../samples/analyses/result_Cc_Cu.pdf')
plt.close()

# plot to show results for every sample
fig = plt.figure(figsize=(7.08662, 6))
x = np.arange(1, 8)
for sample in [0, 1, 2, 3]:
    ax = fig.add_subplot(2, 2, sample+1)
    ax.scatter(
        x,
        df_gt[['d_min', 'd10', 'd30', 'd50', 'd60', 'd90', 'd_max']].iloc[sample].values,
        zorder=10, marker='P', s=80, color='black', edgecolor='white',
        label='true value')

    average_guesses = df_s[[f'S{sample}_d_min', f'S{sample}_d10',
                            f'S{sample}_d30', f'S{sample}_d50',
                            f'S{sample}_d60', f'S{sample}_d90',
                            f'S{sample}_d_max']].mean(axis=0)

    ax.scatter(x, average_guesses, zorder=10, marker='o', color='black', s=80,
               edgecolor='white', label='average estimations')

    parts = ax.violinplot([df_s[f'S{sample}_d_min'], df_s[f'S{sample}_d10'],
                           df_s[f'S{sample}_d30'], df_s[f'S{sample}_d50'],
                           df_s[f'S{sample}_d60'], df_s[f'S{sample}_d90'],
                           df_s[f'S{sample}_d_max']],
                          showextrema=False, showmeans=False, points=20,
                          widths=0.7)
    for pc in parts['bodies']:
        pc.set_facecolor('grey')
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    ax.set_xticks(x)
    ax.set_xticklabels(['$D_{min}$', '$D_{10}$', '$D_{30}$', '$D_{50}$',
                        '$D_{60}$', '$D_{90}$', '$D_{max}$'], fontsize=8)
    ax.tick_params(axis='both', labelsize=8)
    ax.set_title(f'Sample {sample+1}', fontsize=8)
    ax.set_ylabel('particle diameter [mm]', fontsize=8)
    ax.grid(alpha=0.5)
    ax.legend(loc='upper left', fontsize=8)

plt.tight_layout()
plt.savefig(r'../samples/analyses/result_samples.svg')
plt.savefig(r'../samples/analyses/result_samples.pdf')
plt.close()
