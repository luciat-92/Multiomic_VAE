# Variational Autoencoder with Generative Adversarial Network (vae_gan)
1. only gene expression
2. variational autoencoder

## Run
```bash
python train/training.py \
    --folder experiment_2 \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data/preprocessed/gene_expression_all.csv.gz' \
    --ngene='all' \
    --no-norm_feat \
    --only_shared \
    --method='vae_gan'
```

## Run with sbatch
To run all the configurations (samples / feature normalization / n. of genes) use:
```bash
sbatch scripts/train_vae_gan.sbatch
```
