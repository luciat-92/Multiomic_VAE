### Creation enviroment:
```bash
# inside gpu (HPC)
module load cuda11.7/toolkit/11.7.1
module load cudnn8.5-cuda11.7/8.5.0.96
mamba create --name VAE_momics_v2 python=3.9
mamba activate VAE_momics_v2
mamba install pytorch=2.0.1=py3.9_cuda11.7_cudnn8.5.0_0 -c pytorch -c nvidia
mamba install ipykernel notebook
mamba install pandas scikit-learn matplotlib seaborn
mamba install -c conda-forge umap-learn
mamba install datashader bokeh holoviews colorcet scikit-image # for umap plot
mamba install -c plotly plotly=5.18.0
mamba install -c conda-forge python-kaleido
mamba install openpyxl # to open excel
# export the environment
mamba env export | grep -v "prefix" > environment.yml # to export 
# to set-up the environment after cloning the repo:
git clone https://gitlab.fht.org/lucia.trastulla/Multiomic_VAE.git
mamba env create -f environment.yml
```

## Structure:
- `data/`: contains the data used for training
- `experiments/`: each subdirectory corresponds to a different experiment or model variation 
- `models/`: contains files for the model architecture
- `utils/`: contains utility functions for configuration, logging, and visualization and additional utilities (e.g. seed fixing)
- `train/`: contains training and evaluation scripts
- `tests/`: contains unit tests
- `scripts/`: contains bash scripts for running experiments on HPC

## Usage:
1. *data preparation*:
    - python data/preprocessing_gex.py: preprocess the gene expression and save it in data/processed/ (this is done only once)
    - sbatch scripts/preproc_methyl.sbatch: preprocess the methylation and save it in data/processed/ (this is done only once)
2. *training*:
    Each sbatch include multiple configurations with different hyperparameters, the hyperparameter beta controlling VAE is set externally
    - sbatch scripts/train_ae_gan.sbatch  # train with gene expression data and autoencoder
    - sbatch scripts/train_vae_gan_beta<>.sbatch # train with gene expression data and variational autoencoder
    - sbatch scripts/train_mvae_gan_beta<>.sbatch  # train with gene expression and methylation data and variational autoencoder
3. *evaluation*: 
    Evaluate results based on tissue information, tissue type (solid, blood and non-cancerous) and study type (depmap, tgca, gtex)
    - sbatch scripts/eval_ae_gan.sbatch  # train with gene expression data and autoencoder
    - sbatch scripts/eval_vae_gan_beta<>.sbatch # train with gene expression data and variational autoencoder
    - sbatch scripts/eval_mvae_gan_beta<>.sbatch  # train with gene expression and methylation data and variational autoencoder
    - summary of evaluations at experiments/compare_experiments.R # compare results from ex1 and ex2
    - summary of evaluations at experiments/compare_experiments_ex3.R # compare results from ex3
    - plot of correlation between latent space and original features: experiments/plot_corr_encoded_space.R # 

#### Additional information:
In OpenTarget presentation of 28/02/2024 shown:
- VAE gene expression: experiments/experiments_2/vae_gan/samples__ngene_var5000_norm_feat_flag_False_only_shared_True_beta_0.001/
- VAE gene expression and methylation: experiments/experiment_3/mvae_gan/samples_tcgaonly_nfeat_var5000_norm_feat_flag_False_norm_type_zscore_only_shared_True_beta_0.0005/
