# Cyclon trxtools example
Supplemental code for XXXXXX et al. 20XX. This repository provides step-by-step guidance to reproduce data in the MS and serves as a manual for using the trxtools package.
The workflow is divided into two stages (cluster and desktop) according to computational demand, however it can run entirely on as single machine. 
All steps require git and conda installed.

# Cluster processing
## Getting started

### Clone this repository
```
git clone git@github.com:TurowskiLab/example-trxtools-Cyclon.git
```

### Prepare STAR index
```
conda env create -f envs/processing.yml
conda activate processing
STAR --runThreadN 30 --runMode genomeGenerate --genomeDir hg41_STAR_index/ --genomeFastaFiles GRCh38.primary_assembly.genome.cleaned.fa --sjdbGTFfile gencode.v41.annotation.gtf --limitGenomeGenerateRAM 33524399488
```
**NOTE:** Adjust --runThreadN (number of CPU threads to use) and --limitGenomeGenerateRAM (available RAM) to your systems capabilities.

**NOTE:** The STAR index can be saved to any location (--genomeDir).

**IMPORTANT:** You need to specify the path to your STAR index in the ```SM_CRACprocessing3end_all.smk``` file. Set ```STAR_INDEX``` to an absolute path to your index, e.g. ```STAR_INDEX = "/home/user/seq_references/hg41/hg41_STAR_index/" ```

### Create and activate conda environment (you can use mamba instead)
```
conda env create -f envs/snakemake.yml
conda activate snakemake
```
### Run Snakemake file to process raw files 
```
snakemake -c64 --use-conda -s SM_CRACprocessing_read_3end.smk
```
**NOTE:** -c determines number of CPUs to use.

OR if using slurm:
```
snakemake  -c64 --use-conda -s SM_CRACprocessing_read_3end.smk --slurm -j12
```
**NOTE:** The first run of the Snakemake file ```SM_CRACprocessing3end_all.smk``` will initialize new conda environments. This may take several minutes.

After the run finishes continue to the analysis steps below. If you're performing the analysis stage on a different machine (e.g. a desktop after running the pipeline on a cluster), copy the whole repository folder there, including the output files produced.

# Analysis using Jupyter notebooks

## Create and activate jupyter environment
```
conda env create -f envs/jupyter.yml -n jupyter-trxtools-05.yml
conda activate jupyter-trxtools-0.5
```


## Open Jupyter Lab and run notebooks
```
jupyter lab .
```

Afterwards open the subsequent notebooks in Jupyter and run them to perform the analysis.