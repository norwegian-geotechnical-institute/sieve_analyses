# Monte-Carlo-Simulation based revision of sample quantities for grain size analyses
Standards and guidelines that regulate grain size analyses of gravelly soils, such as ISO 17892-4 or ASTM D6913, suggest the use of large quantities of soil (by mass). The rationale for the soil mass to be used is based on formulas that link the maximum grain diameter of the whole sample and the required mass. Such formulas usually lead to a rapid increase in required sample quantities as the max. grain diameter increases (e.g. a max. diameter 200 mm requires 400 kg of soil sample mass).
In practice, implementation of the standard(s) leads to unpractical, massive sample quantities which is especially problematic if the sampling is done based on boreholes which only yield a finite amount of soil.
We hypothesize that the suggestions in the standards vastly overestimate the required sample quantities in coarse grained sediments and the required sample quantity should rather be based on a preliminary estimation of the soil grading. For example, consider the extreme case of a soil which contains only one grain size. In this case a single grain would be sufficient to fully represent the soil’s grain size distribution irrespective of the max. diameter.
Consequently, much smaller sample quantities would be enough to sufficiently represent the underlying grain size distribution in most cases.
The goal of this project is to develop a new recommendation for minimum required sample quantities based on laboratory testing, and statistical methods such as monte-carlo-simulations of virtual sieve analyses.

## project organisation
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
conda activate MLpFEM
```
