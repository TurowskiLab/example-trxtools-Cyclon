#!/bin/bash
#SBATCH --job-name=BigWig_3end_Rpo21
#SBATCH --output=BigWig_3end_Rpo21%j.out
#SBATCH --error=BigWig_3end_Rpo21%j.err
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
for f in *sam; do
    SAM2profilesGenomic.py -f $f -u 3end -n -s polyA &
done

# Wait for all background jobs to finish
wait

# Deactivate the Conda environment
conda deactivate
