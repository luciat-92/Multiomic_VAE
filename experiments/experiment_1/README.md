# Autoencoder with Generative Adversarial Network (ae_gan)
1. only gene expression 
2. autoencoder

## Run
```bash
python train/training.py \
    --folder experiment_1 \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data/preprocessed/gene_expression_all.csv.gz' \
    --ngene='all' \
    --no-norm_feat
```

## Run with sbatch
All the different settings are tested
```bash
sbatch scripts/train_ae_gan.sbatch
```