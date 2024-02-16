# Variational Autoencoder with Generative Adversarial Network (vae_gan)
1. only gene expression
2. variational autoencoder

## Run
```bash
python train/training.py \
    --folder experiment_2 \
    --gex_feature_file='/group/iorio/lucia/Multiomic_VAE/data/preprocessed/gene_expression_var1000.csv.gz' \
    --ngene='var1000' \
    --no-norm_feat \
    --no-only_shared \
    --method='vae_gan' \
    --beta=0.001
```

## Run with sbatch
To run all the configurations (samples / feature normalization / n. of genes) use:
```bash
sbatch scripts/train_vae_gan.sbatch # beta not specified: is based on 0.0005
```
