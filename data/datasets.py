import sys, os
sys.path.insert(0, os.getcwd())
import torch
import numpy as np
import gzip
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
from torch.utils.data import TensorDataset, DataLoader

from utils.utils import *
import utils.config as config
from data.preprocessing_gex import load_data

# load data sample
output_gex = load_data(load_gex = False)
depmap_sample_df = output_gex[2]
xena_sample_df = output_gex[3]

def get_dataloaders(depmap_sample_df,
                    xena_sample_df, 
                    seed, 
                    batch_size, 
                    gene_df = None,
                    meth_df = None,
                    data_type = 'both',
                    normalize_features=False, 
                    norm_type='zscore'):
    
    print("Loading data")
    # concatenate gene and methylation data (if requested)
    if (gene_df is None) & (meth_df is None):
        raise ValueError('At least one of gene_df and meth_df should be provided')
    elif (gene_df is not None) & (meth_df is not None):
        # get common samples located in the first column
        common_samples = gene_df.index.intersection(meth_df.index)
        gene_df = gene_df.loc[common_samples]
        meth_df = meth_df.loc[common_samples]
        meth_df = meth_df.groupby(meth_df.index).mean() # remove duplicates if present
        if data_type == 'both':
            print('Concatenating gene and methylation data')
            tot_df = pd.concat([gene_df, meth_df], axis=1)
        elif data_type == 'gene':
            print('Consider only gene data')
            tot_df = gene_df
        elif data_type == 'methyl':
            print('Consider only methylation data')
            tot_df = meth_df
        else:
            raise ValueError(f'Data type {data_type} not supported')
    elif gene_df is not None:
        tot_df = gene_df
    else:
        meth_df = meth_df.groupby(meth_df.index).mean()
        tot_df = meth_df

    set_seed(seed)
    # filter samples based on the combined dataset intersection
    xena_sample_df = xena_sample_df.dropna(subset=['primary disease or tissue'])
    depmap_sample_df = depmap_sample_df.dropna(subset = ['OncotreePrimaryDisease'])
    
    depmap_samples = depmap_sample_df.index.intersection(tot_df.index)
    xena_samples = xena_sample_df.index.intersection(tot_df.index)

    # match all the resources
    depmap_sample_df = depmap_sample_df.loc[depmap_samples]
    xena_sample_df = xena_sample_df.loc[xena_samples]

    # filter tot_df for both depmap and xena samples
    tot_df = tot_df.loc[depmap_samples.union(xena_samples)]
    # # normalize features
    if normalize_features:
        if norm_type == 'zscore':
            tot_df = (tot_df - tot_df.mean()) / tot_df.std()
        elif norm_type == 'minmax':
             tot_df = (tot_df - tot_df.min()) / (tot_df.max() - tot_df.min())
        else:
            raise ValueError(f'Normalization type {norm_type} not supported')
        
    depmap_df = tot_df.loc[depmap_samples]
    xena_df = tot_df.loc[xena_samples]

    print(f'Loaded CL data with shape: {depmap_df.shape} and patient data with shape: {xena_df.shape}')

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

    # prepare data for AE, torch.from_numpy converts numpy array to torch tensor without copying the underlying data
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
    # not used now! trained on the entire set
    #train_xena_dataloader = DataLoader(train_xena_dateset,
    #                                   batch_size=batch_size,
    #                                   shuffle=True)
    test_xena_dataloader = DataLoader(test_xena_dateset,
                                      batch_size=batch_size,
                                      shuffle=True)

    depmap_data_loader = DataLoader(depmap_dataset,
                                  batch_size=batch_size,
                                  shuffle=True,
                                  drop_last=True) # drop the last batch if the sample size is lower than the settled batch size
    # not used now! trained on the entire set
    #train_depmap_dataloader = DataLoader(train_depmap_dateset,
    #                                   batch_size=batch_size,
    #                                   shuffle=True, 
    #                                   drop_last=True)
    test_depmap_dataloader = DataLoader(test_depmap_dateset,
                                      batch_size=batch_size,
                                      shuffle=True)
    return (depmap_data_loader, test_depmap_dataloader, depmap_sample_df), (xena_dataloader, test_xena_dataloader, xena_sample_df)
