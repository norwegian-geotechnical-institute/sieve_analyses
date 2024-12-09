# Monte-Carlo-Simulation based revision of sample quantities for grain size analyses
Determining particle size distributions (PSD) of soils is a basic first step in many geotechnical analyses and guidance is given in different national standards. For ambiguous reasons, the recommended required minimum sample mass (m_min) for the PSD-analyses of soils with a main component of gravel or greater is based on equations including the soil's maximum grain diameter (D_max). We claim that the recommended m_min is overestimated in many cases as D_max does not represent the relevant large soil fraction but only the PSD's uppermost outlier. Furthermore, the recommended m_min is not based on a specific sampling confidence (i.e. how closely does the sample’s PSD need to approximate the soil’s PSD?) and thus it is not clear why the m_min should even be necessary. We conducted Monte-Carlo simulation-based sieve analyses of coarse-grained soils and developed a new, practically applicable framework to determine m_min based on D_90 that also includes explicit consideration of sampling confidence. A survey was conducted that shows that there is no difference in how well operators are able to assess parameters like D_90 or D_max. Real sieve tests performed on three different sands and gravels corroborate the theoretical results and show that substantially lower sample masses yield PSDs with only marginal differences to PSDs from samples according to the standards. While the results are promising, they open up for new research questions about which geotechnical application requires which soil sampling confidence. 

## Publication
This repository contains the code for the above described Monte-Carlo simulation and is published alongside a paper that is submitted to the journal **Engineering Geology**. The paper has the title *"Uncertainty aware sample mass determination of coarse-grained soils for particle size analyses"* and was originally submitted in August 2024, then rejected and a major revision was submitted in December 2024.

A preprint of the paper can be found on EarthArXiv: https://eartharxiv.org/repository/view/7592/

## Repository and code structure

The code framework depends on a certain folder structure. The main functionality is in the src folder. Here are mainly two types of files:

- "<A, B, C, D, E>"_description - scripts to be run
- X_library - custom library with functions to run scripts.

```
sieve_analyses
├── figures                       - Figures that are produced from the simulations.
├── laboratory                    - Real laboratory test results and figures from them.
├── samples                       - Folder containing simulated soil samples and results for a survey that investigates how well people can visually assess grain sizes.
├── simulations                   - Collections of simulation data / recordings.
├── src                           - Source code files.
│   ├── A_main.py                 - Script to make individual virtual grain size analyses with different sample masses and different underlzing distributions.
│   ├── B_MonteCarlo.py           - Script to make simulated sieve analyses for Monte Carlo analyses to determine correlations between distribution geometries and sample masses.
│   ├── C_model_builder.py        - Script to primarily create visualizations as a means to build equations of relationships between parameters.
│   ├── D_grain_size_plots.py     - Script to create visualizations of simulated soil samples for a survey.
│   ├── E_survey_analyses.py      - Script to analyse the survey results.
│   ├── X_library.py              - Custom libraries for simulated sieve analyses to develop a new way to determine minimum required soil sample masses.
├── .gitignore                    - Gitignore file that specifies what should not be synchronized.
├── environment.yaml              - Dependency file to use with conda.
├── LICENSE                       - File clarifying the repository's license.
├── README.md                     - Descritpion of the repository.
```

To clone the repository and have everything set up for you, run:

```bash
git clone https://github.com/norwegian-geotechnical-institute/sieve_analyses.git
```

## Project organisation
The investigation is conducted by the **Norwegian Geotechnical Institute (NGI)** together with the **Graz University of Technology - Institute of Soil Mechanics, Foundation Engineering and Computational Geotechnics (IBG)**
- **NGI**: Georg H. Erharter, Diana Cordeiro, Santiago Quinteros
- **IBG**: Franz Tschuchnigg, Matthias Rebhan

contact: georg.erharter@ngi.no

## Requirements

The environment is set up using `conda`.

To do this create an environment called `sieve_analyses` using `environment.yaml` with the help of `conda`. If you get pip errors, install pip libraries manually, e.g. `pip install pandas`
```bash
conda env create --file environment.yaml
```

Activate the new environment with:

```bash
conda activate sieve_analyses
```
