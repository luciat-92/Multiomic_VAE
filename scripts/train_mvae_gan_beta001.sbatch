#!/bin/bash
#SBATCH --job-name=train_adv
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucia.trastulla@fht.org
#SBATCH --partition=gpuq
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/train_%j.out.log
#SBATCH --error=logs/train_%j.err.log
#SBATCH --mem=90G
#SBATCH --gres=gpu:1

module load cuda11.7/toolkit/11.7.1
module load cudnn8.5-cuda11.7/8.5.0.96
source ${HOME}/miniforge3/etc/profile.d/conda.sh
conda activate VAE_momics_v2

GROUP_PATH='/group/iorio/lucia/'

# var1000
# test 1
python train/training.py \
    --folder=experiment_3 \
    --gex_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/gene_expression_var1000_tcgaonly.csv.gz' \
    --meth_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/methylation_var1000_tcgaonly.csv.gz' \
    --nfeat='var1000' \
    --norm_feat \
    --no-only_shared \
    --norm_type='zscore' \
    --method='mvae_gan' \
    --beta=0.001 \
    --samples='tcgaonly'

# test 2
python train/training.py \
    --folder=experiment_3 \
    --gex_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/gene_expression_var1000_tcgaonly.csv.gz' \
    --meth_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/methylation_var1000_tcgaonly.csv.gz' \
    --nfeat='var1000' \
    --norm_feat \
    --no-only_shared \
    --norm_type='minmax' \
    --method='mvae_gan' \
    --beta=0.001 \
    --samples='tcgaonly'

# test 3
python train/training.py \
    --folder=experiment_3 \
    --gex_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/gene_expression_var1000_tcgaonly.csv.gz' \
    --meth_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/methylation_var1000_tcgaonly.csv.gz' \
    --nfeat='var1000' \
    --no-norm_feat \
    --no-only_shared \
    --method='mvae_gan' \
    --beta=0.001 \
    --samples='tcgaonly'

# var5000
# test 4
python train/training.py \
    --folder=experiment_3 \
    --gex_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/gene_expression_var5000_tcgaonly.csv.gz' \
    --meth_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/methylation_var5000_tcgaonly.csv.gz' \
    --nfeat='var5000' \
    --norm_feat \
    --no-only_shared \
    --norm_type='zscore' \
    --method='mvae_gan' \
    --beta=0.001 \
    --samples='tcgaonly'

# test 5
python train/training.py \
    --folder=experiment_3 \
    --gex_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/gene_expression_var5000_tcgaonly.csv.gz' \
    --meth_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/methylation_var5000_tcgaonly.csv.gz' \
    --nfeat='var5000' \
    --norm_feat \
    --no-only_shared \
    --norm_type='minmax' \
    --method='mvae_gan' \
    --beta=0.001 \
    --samples='tcgaonly'

# test 6
python train/training.py \
    --folder=experiment_3 \
    --gex_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/gene_expression_var5000_tcgaonly.csv.gz' \
    --meth_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/methylation_var5000_tcgaonly.csv.gz' \
    --nfeat='var5000' \
    --no-norm_feat \
    --no-only_shared \
    --method='mvae_gan' \
    --beta=0.001 \
    --samples='tcgaonly'

# # all
# # test 7
# python train/training.py \
#     --folder=experiment_3 \
#     --gex_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/gene_expression_all.csv.gz' \
#     --meth_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/methylation.csv.gz' \
#     --nfeat='all' \
#     --norm_feat \
#     --no-only_shared \
#     --norm_type='zscore' \
#     --method='mvae_gan' \
#     --beta=0.001 \
#     --samples='tcgaonly'

# # test 8
# python train/training.py \
#     --folder=experiment_3 \
#     --gex_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/gene_expression_all.csv.gz' \
#     --meth_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/methylation.csv.gz' \
#     --nfeat='all' \
#     --norm_feat \
#     --no-only_shared \
#     --norm_type='minmax' \
#     --method='mvae_gan' \
#     --beta=0.001 \
#     --samples='tcgaonly'

# # test 9
# python train/training.py \
#     --folder=experiment_3 \
#     --gex_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/gene_expression_all.csv.gz' \
#     --meth_feature_file=${GROUP_PATH}'Multiomic_VAE/data/preprocessed/methylation.csv.gz' \
#     --nfeat='all' \
#     --no-norm_feat \
#     --no-only_shared \
#     --method='mvae_gan' \
#     --beta=0.001 \
#     --samples='tcgaonly'
