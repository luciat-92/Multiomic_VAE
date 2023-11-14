### Creation enviroment:
```
conda create --name VAE python=3.9
conda activate VAE
conda install pandas numpy scikit-learn
conda install pytorch torchvision -c pytorch
conda env export | grep -v "prefix" > environment.yml # to export 
```


