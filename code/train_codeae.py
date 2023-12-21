import pandas as pd
import torch
import json
import os
import sys
import argparse
import random
import pickle
import itertools
import gzip
import umap
import matplotlib.pyplot as plt
import seaborn as sns

import data_preproc
import train_code_adv
import train_code_base

sys.path.append("code/")

# 12 different models to run: 3 with TCGAonly *2 (adv and base) and 3 with all *2 (adv and base)
# we could also add another flag to initially normalize the data (check publication)
# this would lead in total 24 different models to run

# different from the original code
def generate_encoded_features(encoder, dataloader, normalize_flag=False):
    encoder.eval()
    raw_feature_tensor = dataloader[0].dataset.tensors[0].cpu()
    sample_info = dataloader[2]
    with torch.no_grad():
        encoded_feature_tensor = encoder.cpu()(raw_feature_tensor)
    if normalize_flag:
        encoded_feature_tensor = torch.nn.functional.normalize(encoded_feature_tensor, p=2, dim=1)
    encoded_feature_tensor = encoded_feature_tensor.cpu().detach().numpy()
    enc_out = pd.DataFrame(encoded_feature_tensor, index = sample_info.index)
    return enc_out

def safe_make_dir(new_folder_name):
    if not os.path.exists(new_folder_name):
        os.makedirs(new_folder_name)
    else:
        print(new_folder_name, 'exists!')

def dict_to_str(d):
    return "_".join(["_".join([k, str(v)]) for k, v in d.items()])

def get_umap(encoded_features, depmap_sample_df, xena_sample_df, save_folder):
    depmap_sample_df['type'] = 'depmap'
    depmap_sample_df['study'] = 'CL_depmap'
    xena_sample_df['type'] = 'xena'
    t1 = depmap_sample_df[["type", "Sex", "OncotreeLineage", "study"]].rename(columns={"Sex":"sex", "OncotreeLineage":"site"})
    t2 = xena_sample_df[["type", "_primary_site", "_gender", "_study", "_sample_type"]].rename(columns={"_primary_site": "site", "_gender":"sex", "_study":"study"})
    tot_sample_df = pd.concat([t1, t2], axis=0, sort=False)
    # order tot_sample_df according tmp
    tot_sample_df = tot_sample_df.loc[encoded_features.index]
    embedding = umap.UMAP(random_state=42).fit_transform(encoded_features)
    umap_df = pd.DataFrame(embedding, index = encoded_features.index, columns = ['umap_1', 'umap_2'])
    umap_df = pd.concat([umap_df, tot_sample_df], axis = 1)
    umap_df.to_csv(os.path.join(save_folder, f'umap.csv'))
    return umap_df

def plot_umap(umap_df, save_folder):

    # plot
    plt.figure(figsize=(5, 5))
    sns.scatterplot(data=umap_df ,x='umap_1', y='umap_2', hue='site', style='type', markers = ['o', '.'], alpha = 0.4, linewidth=0.1)
    #Set plot title and labels
    plt.xlabel("UMAP Dimension 1", fontsize=10)
    plt.ylabel("UMAP Dimension 2", fontsize=10)
    plt.legend(fontsize=6, title_fontsize=6, loc="upper left", bbox_to_anchor=(1.05, 1))
    # save plot
    plt.savefig(os.path.join(save_folder, f'umap_lineage.pdf'), format = 'pdf', bbox_inches='tight')

    plt.figure(figsize=(5, 5))
    sns.scatterplot(data=umap_df, x='umap_1', y='umap_2', hue='study', style='type', markers = ['o', '.'], alpha = 0.4, linewidth=0.1)
    # Set plot title and labels
    plt.xlabel("UMAP Dimension 1", fontsize=10)
    plt.ylabel("UMAP Dimension 2", fontsize=10)
    plt.legend(fontsize=6, title_fontsize=6, loc="upper left", bbox_to_anchor=(1.05, 1))
    plt.savefig(os.path.join(save_folder, f'umap_study.pdf'), format = 'pdf', bbox_inches='tight')


def main(args, update_params_dict):
    if args.method == 'code_base':
        train_fn = train_code_base.train_code_base
    else:
        train_fn = train_code_adv.train_code_adv

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(device)
    with gzip.open(args.gex_feature_file) as f:
        gene_df = pd.read_csv(f, sep=',', index_col=0)

    training_params = {"batch_size": 64, 
    "lr": 0.0001, 
    "pretrain_num_epochs": 500, 
    "train_num_epochs": 1000, 
    "alpha": 1.0, 
    "classifier_hidden_dims": [64, 32], 
    "encoder_hidden_dims": [512, 256], 
    "latent_dim": 128, 
    "dop": 0.1, 
    "device": device,
    "input_dim": gene_df.shape[-1],
    "model_save_folder": os.path.join(args.root, "output/"),
    "es_flag": False,
    "norm_flag": args.norm_flag, 
    "retrain_flag": True}

    training_params.update(update_params_dict)
    param_str = dict_to_str(update_params_dict)

    if not args.norm_flag:
        method_save_folder = os.path.join(training_params['model_save_folder'], 'model_save', f'{args.method}')
    else:
        method_save_folder = os.path.join(training_params['model_save_folder'], 'model_save', f'{args.method}_norm')
    training_params.update(
        {
            'model_save_folder': os.path.join(method_save_folder, param_str),
        })

    safe_make_dir(training_params['model_save_folder'])
    print(training_params['model_save_folder'])
    random.seed(2020)

    if training_params['samples'] == 'tcgaonly':
        data_preproc.xena_sample_df = data_preproc.xena_sample_df[data_preproc.xena_sample_df['_study'] == 'TCGA']

    s_dataloaders, t_dataloaders = data_preproc.get_unlabeled_dataloaders(
        gene_df = gene_df,
        depmap_sample_df = data_preproc.depmap_sample_df,
        xena_sample_df = data_preproc.xena_sample_df,
        seed=2020,
        batch_size=training_params['batch_size'], 
        normalize_features = args.norm_feat
    )

    # start unlabeled training
    encoder, historys = train_fn(s_dataloaders=s_dataloaders,
                                 t_dataloaders=t_dataloaders,
                                 **training_params)
    with open(os.path.join(training_params['model_save_folder'], f'unlabel_train_history.pickle'),
              'wb') as f:
        for history in historys:
            pickle.dump(dict(history), f)
    
    enc_depmap = generate_encoded_features(encoder, s_dataloaders)
    enc_xena = generate_encoded_features(encoder, t_dataloaders)
    enc_tot = pd.concat([enc_depmap, enc_xena], axis = 0, sort = False)
    enc_tot.to_csv(os.path.join(training_params['model_save_folder'], f'encoded_features.csv'))

    # plot umap
    umap_df = get_umap(enc_tot, data_preproc.depmap_sample_df, data_preproc.xena_sample_df, training_params['model_save_folder'])
    plot_umap(umap_df, training_params['model_save_folder'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser('training and evaluation')
    parser.add_argument('--method', dest='method', nargs='?', default='code_adv',
                        choices=['code_adv', 'code_base'])
    norm_group = parser.add_mutually_exclusive_group(required=False)
    norm_group.add_argument('--norm', dest='norm_flag', action='store_true')
    norm_group.add_argument('--no-norm', dest='norm_flag', action='store_false')
    parser.set_defaults(norm_flag=True)
    parser.add_argument('--gex_feature_file', dest='gex_feature_file', nargs='?')
    parser.add_argument('--samples', dest='samples', nargs='?', default='',
                       choices=['', 'tcgaonly'])
    parser.add_argument('--ngene', dest='ngene', nargs='?', default='all')
    parser.add_argument('--root', dest='root', nargs='?', default='.')
    
    norm_feat_group = parser.add_mutually_exclusive_group(required=False)
    norm_feat_group.add_argument('--norm_feat', dest='norm_feat', action='store_true')
    norm_feat_group.add_argument('--no-norm_feat', dest='norm_feat', action='store_false')
    parser.set_defaults(norm_feat=True)
    
    args = parser.parse_args()
    print(f'current config is {args}')
    # params_grid = {
    #     "pretrain_num_epochs": [0, 100, 300],
    #     "train_num_epochs": [100, 200, 300, 500, 750, 1000, 1500, 2000, 2500, 3000],
    #     "dop": [0.0, 0.1]
    # }
    params_grid = {
        "pretrain_num_epochs": [500],
        "train_num_epochs": [1000],
        "dop": [0.0],
        "samples": [args.samples],
        "ngene": [args.ngene],
        "norm_feat_flag": [args.norm_feat]
    }
    if args.method not in ['code_adv', 'adsn', 'adae', 'dsnw']:
        params_grid.pop('pretrain_num_epochs')

    keys, values = zip(*params_grid.items())
    update_params_dict_list = [dict(zip(keys, v)) for v in itertools.product(*values)]

    for param_dict in update_params_dict_list:
        main(args=args, update_params_dict=param_dict)
