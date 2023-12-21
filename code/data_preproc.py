
import random
import torch
import numpy as np
import os
import gzip
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
from torch.utils.data import TensorDataset, DataLoader


# load data sample
ROOT = "/group/iorio/lucia/"
FOLD_PROJECT = os.path.join(ROOT, "Multiomic_VAE/")
xena_folder = os.path.join(ROOT, "datasets/XENA/TCGA_TARGET_GTEx/")
xena_sample_file = os.path.join(xena_folder, 'TcgaTargetGTEX_phenotype.txt.gz')
depmap_folder = os.path.join(ROOT, "datasets/DEPMAP_PORTAL/version_23Q2/")
depmap_sample_file = os.path.join(depmap_folder, "Model.csv")
# load tissue data
with gzip.open(xena_sample_file) as f:
    xena_sample_df = pd.read_csv(f, sep='\t', index_col=0, encoding='ISO-8859-1')
xena_sample_df.index.name = 'sample_id'
# load depmap data
depmap_sample_df = pd.read_csv(depmap_sample_file, sep=',', index_col=0)
depmap_sample_df.index.name = 'sample_id'

# correct site annotation
map_common_annotation = {"Vagina":"Vulva/Vagina", 
                         "Thyroid Gland": "Thyroid", 
                         "SympatheticÊNervous System": "Peripheral Nervous System", 
                         "Stomach" : "Esophagus/Stomach",
                         "Soft tissue,Bone" : "Soft Tissue", 
                         "Head and Neck region" : "Head and Neck", 
                         "Paraganglia" : "Peripheral Nervous System", 
                         "Ovary" : "Ovary/Fallopian Tube", 
                         "Lymphoid" : "White blood cell", 
                         "Fallopian Tube" : "Ovary/Fallopian Tube", 
                         "Esophagus" : "Esophagus/Stomach", 
                         "Cervix Uteri" : "Cervix", 
                         "Brain" : "CNS/Brain", 
                         "Adrenal gland" : "Adrenal Gland",
                         "Myeloid" : "White blood cell", 
                         "Bladder/Urinary Tract": "Bladder", 
                         "Bile duct": "Biliary Tract"}
depmap_sample_df['OncotreePrimaryDisease'] = depmap_sample_df['OncotreePrimaryDisease'].replace(map_common_annotation)
xena_sample_df["_primary_site"] = xena_sample_df["_primary_site"].replace(map_common_annotation)

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.device_count() > 0:
        torch.cuda.manual_seed_all(seed)
        
def get_unlabeled_dataloaders(gene_df, depmap_sample_df, xena_sample_df, seed, batch_size, normalize_features=False):
    
    set_seed(seed)
    # filter samples based on the combined dataset intersection
    xena_sample_df = xena_sample_df.dropna(subset=['primary disease or tissue'])
    depmap_sample_df = depmap_sample_df.dropna(subset = ['OncotreePrimaryDisease'])
    
    depmap_samples = depmap_sample_df.index.intersection(gene_df.index)
    xena_samples = xena_sample_df.index.intersection(gene_df.index)

    # match all the resources
    depmap_sample_df = depmap_sample_df.loc[depmap_samples]
    xena_sample_df = xena_sample_df.loc[xena_samples]

    # filter gene_df for both depmap and xena samples
    gene_df = gene_df.loc[depmap_samples.union(xena_samples)]
    # normalize features
    if normalize_features:
        gene_df = (gene_df - gene_df.mean()) / gene_df.std()
    depmap_df = gene_df.loc[depmap_samples]
    xena_df = gene_df.loc[xena_samples]

    # exclude samples with only one sample per disease 
    excluded_depmap_samples = []
    excluded_depmap_samples.extend(depmap_df.index.difference(depmap_sample_df.index))
    excluded_depmap_diseases = depmap_sample_df.OncotreePrimaryDisease.value_counts()[
        depmap_sample_df.OncotreePrimaryDisease.value_counts() < 2].index
    excluded_depmap_samples.extend(
        depmap_sample_df[depmap_sample_df.OncotreePrimaryDisease.isin(excluded_depmap_diseases)].index)
    to_split_depmap_df = depmap_df[~depmap_df.index.isin(excluded_depmap_samples)]

    # split train and test (?? IS IT NECESSARY ??)
    train_depmap_df, test_depmap_df = train_test_split(
        to_split_depmap_df, 
        test_size=0.1,
        stratify=depmap_sample_df.loc[to_split_depmap_df.index].OncotreePrimaryDisease,
        random_state=seed)
    test_depmap_df = pd.concat([test_depmap_df, depmap_df.loc[excluded_depmap_samples]], ignore_index=True)
   
    train_xena_df, test_xena_df = train_test_split(
        xena_df, 
        test_size=len(test_depmap_df) / len(xena_df),
        stratify=xena_sample_df['primary disease or tissue'],
        random_state=seed)

    # prepare data for VAE, torch.from_numpy converts numpy array to torch tensor without copying the underlying data
    xena_dataset = TensorDataset(
        torch.from_numpy(xena_df.values.astype('float32'))
    )
    depmap_dataset = TensorDataset(
        torch.from_numpy(depmap_df.values.astype('float32'))
    )
    train_xena_dateset = TensorDataset(
        torch.from_numpy(train_xena_df.values.astype('float32')))
    test_xena_dateset = TensorDataset(
        torch.from_numpy(test_xena_df.values.astype('float32')))
    train_depmap_dateset = TensorDataset(
        torch.from_numpy(train_depmap_df.values.astype('float32')))
    test_depmap_dateset = TensorDataset(
        torch.from_numpy(test_depmap_df.values.astype('float32')))

    xena_dataloader = DataLoader(xena_dataset,
                                 batch_size=batch_size,
                                 shuffle=True)
    train_xena_dataloader = DataLoader(train_xena_dateset,
                                       batch_size=batch_size,
                                       shuffle=True)
    test_xena_dataloader = DataLoader(test_xena_dateset,
                                      batch_size=batch_size,
                                      shuffle=True)

    depmap_data_loader = DataLoader(depmap_dataset,
                                  batch_size=batch_size,
                                  shuffle=True,
                                  drop_last=True) # drop the last batch if the sample size is lower than the settled batch size
    train_depmap_dataloader = DataLoader(train_depmap_dateset,
                                       batch_size=batch_size,
                                       shuffle=True, 
                                       drop_last=True)
    test_depmap_dataloader = DataLoader(test_depmap_dateset,
                                      batch_size=batch_size,
                                      shuffle=True)
    return (depmap_data_loader, test_depmap_dataloader, depmap_sample_df), (xena_dataloader, test_xena_dataloader, xena_sample_df)
