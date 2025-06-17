# generic_human_analysis
Supplemental code for generic human CRAC MS. This repository provides step-by-step guidance to reproduce data in the MS.
It is divided into two stages (cluster and desktop) according to computational demand. 
All steps require git and conda installed

# Cluster processing
## Getting started

Clone this repository
```
git clone git@github.com:tturowski/generic_human_analysis.git
```

Create and activate conda environment (you can use mamba instead)
```
conda env create -f envs/snakemake.yml
conda activate snakemake
```
Run SnakeMake file to process raw files (-c determines number of CPUs to use)
```
snakemake -c64 --use-conda -s SM_CRACprocessing3end_all.py
```
OR using slurm
```
snakemake  -c64 --use-conda -s SM_CRACprocessing.smk --slurm -j12
```
**NOTE:** Very first run of the SnakeMake file ```SM_CRACprocessing3end_all.py``` will initialize new conda environments. This may take a several minutes.
