### Creation enviroment:
```
# inside gpu
module load cuda11.7/toolkit/11.7.1
module load cudnn8.5-cuda11.7/8.5.0.96
mamba create --name VAE_momics_v2 python=3.9
mamba activate VAE_momics_v2
mamba install pytorch=2.0.1=py3.9_cuda11.7_cudnn8.5.0_0 -c pytorch -c nvidia
mamba install ipykernel notebook
# to make visualizations
mamba install pandas scikit-learn matplotlib seaborn
mamba install -c plotly plotly=5.18.0
mamba install -c conda-forge python-kaleido
mamba env export | grep -v "prefix" > environment.yml # to export 
```