#!/bin/bash
#SBATCH --job-name=train_adv
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucia.trastulla@fht.org
#SBATCH --partition=gpuq
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/train_%j.out.log
#SBATCH --error=logs/train_%j.err.log
#SBATCH --mem=20G
#SBATCH --gres=gpu:1

module load cuda11.7/toolkit/11.7.1
module load cudnn8.5-cuda11.7/8.5.0.96
source ${HOME}/miniforge3/etc/profile.d/conda.sh
conda activate VAE_momics_v2

# test 1
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_var1000_tcgaonly.csv.gz' \
    --samples='tcgaonly' \
    --ngene='var1000' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --no-norm_feat 

# test 2
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_var1000_tcgaonly.csv.gz' \
    --samples='tcgaonly' \
    --ngene='var1000' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --norm_feat 

# test 3
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_var5000_tcgaonly.csv.gz' \
    --samples='tcgaonly' \
    --ngene='var5000' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --no-norm_feat 

# test 4
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_var5000_tcgaonly.csv.gz' \
    --samples='tcgaonly' \
    --ngene='var5000' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --norm_feat 

# test 5
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_all.csv.gz' \
    --samples='tcgaonly' \
    --ngene='all' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --no-norm_feat 

# test 6
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_all.csv.gz' \
    --samples='tcgaonly' \
    --ngene='all' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --norm_feat 

# all samples
# test 1
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_var1000.csv.gz' \
    --ngene='var1000' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --no-norm_feat 

# test 2
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_var1000.csv.gz' \
    --ngene='var1000' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --norm_feat 

# test 3
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_var5000.csv.gz' \
    --ngene='var5000' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --no-norm_feat 

# test 4
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_var5000.csv.gz' \
    --ngene='var5000' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --norm_feat 

# test 5
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_all.csv.gz' \
    --ngene='all' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --no-norm_feat 

# test 6
python code/train_codeae.py \
    --method='code_adv' \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data_preprocessing/gene_expression_all.csv.gz' \
    --ngene='all' \
    --root='/group/iorio/lucia/Multiomic_VAE/' \
    --norm_feat 

