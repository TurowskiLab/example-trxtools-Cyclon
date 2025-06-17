#!/bin/bash
#SBATCH --job-name=polyAfraction
#SBATCH --output=polyAfraction%j.out
#SBATCH --error=polyAfraction%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=28-00
#SBATCH --partition=long

# Load any required modules
# module load your_module_name

# Activate the Conda environment
source /home/${USER}/.bashrc
source activate processing

# Navigate to the directory containing your SAM files
cd $PWD/04_BigWig

# Loop through each SAM file and run SAM2profilesGenomic.py
for f in *polyA_fwd.bw; do
    polyAfraction.py -f $f &
done

# Wait for all background jobs to finish
wait

# Deactivate the Conda environment
conda deactivate
