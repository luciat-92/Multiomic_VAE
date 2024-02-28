#!/bin/bash
#SBATCH --job-name=prepr_methyl
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucia.trastulla@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/preproc_methyl_%j.out.log
#SBATCH --error=logs/preproc_methyl_%j.err.log
#SBATCH --mem=250G

module load R/4.1.0
Rscript data/preprocessing_methyl.R

# gzip output
gzip data/preprocessed/methylation_var1000_tcgaonly.csv
gzip data/preprocessed/methylation_var5000_tcgaonly.csv
gzip data/preprocessed/methylation.csv