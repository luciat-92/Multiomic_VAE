import os, sys
sys.path.insert(0, os.getcwd())
import pandas as pd
import torch
import json
import argparse
import pickle
import itertools
import gzip
import logging

from utils.utils import *
import utils.config as config
import data.datasets as data_loader
from train.helper_plots import *
from train.helper_utils import *
import train.helper_train_GAN as train_code
from utils.logger import setup_logging, get_logger

@timeit
def main(args, update_params_dict):

    logger = get_logger()
    logger.info(f'current config is {args}')

    if args.method == 'ae_gan':
        train_fn = train_code.train_code_adv
        var_flag = False
    else:
        if args.method == 'vae_gan':
            train_fn = train_code.train_code_adv_vae
            var_flag = True
        else:
            logger.error(f'Invalid method: {args.method}')
    
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    logger.info(f'Using device: {device}')

    with gzip.open(args.gex_feature_file) as f:
        gene_df = pd.read_csv(f, sep=',', index_col=0)
    logger.info(f'Loaded gene_df with shape {gene_df.shape}')
    
    with open(str2path('experiments/train_params.json'), 'r') as f:
        training_params = json.load(f)

    training_params.update(update_params_dict)
    param_str = dict_to_str(update_params_dict)
    method_save_folder = config.RESULTS_DIR / args.folder / f'{args.method}'

    training_params.update(
        {
            "device": device,
            "input_dim": gene_df.shape[-1],
            'model_save_folder': method_save_folder / param_str
        })
    
    print(training_params['model_save_folder'])
    safe_create_dir(training_params['model_save_folder'])

    set_seed(2024)

    if training_params['samples'] == 'tcgaonly':
        data_loader.xena_sample_df = data_loader.xena_sample_df[data_loader.xena_sample_df['_study'] == 'TCGA']

    logger.info(f'Feature normalization: {args.norm_feat}')
    s_dataloaders, t_dataloaders = data_loader.get_dataloaders(
        gene_df = gene_df,
        depmap_sample_df = data_loader.depmap_sample_df,
        xena_sample_df = data_loader.xena_sample_df,
        seed=2024,
        batch_size=training_params['batch_size'], 
        normalize_features = args.norm_feat
    )
    logger.info(f'Loaded dataloaders with s_dataloaders: {len(s_dataloaders[0])} and t_dataloaders: {len(t_dataloaders[0])}')
    # check the shape of the data (TODO: remove this part later on)
    for data in s_dataloaders[0]:  
        print(data[0].shape)
        break
    for data in t_dataloaders[0]:  
        print(data[0].shape)
        break

    # start unlabeled training
    logger.info(training_params)
    encoder, historys = train_fn(s_dataloaders=s_dataloaders,
                                 t_dataloaders=t_dataloaders,
                                 **training_params)
    logger.info(f'Trained model')
    # save history (hardcoded info of the historys list)
    with open(training_params['model_save_folder'] / 'pretrain_history.pickle', 'wb') as f:
        pickle.dump(dict(historys[0]), f)
    with open(training_params['model_save_folder'] / 'critic_history.pickle', 'wb') as f:
        pickle.dump(dict(historys[2]), f)
    with open(training_params['model_save_folder'] / 'trainGAN_history.pickle', 'wb') as f:
        pickle.dump(dict(historys[3]), f)

    # plot history
    plot_histories_gan(histories = historys, 
                  save_folder = training_params['model_save_folder'] / 'plots/', 
                  pretrain_num_epochs=training_params['pretrain_num_epochs'], 
                  train_num_epochs=training_params['train_num_epochs'], 
                  variational_flag = var_flag)
    
    enc_depmap = generate_encoded_features(encoder, s_dataloaders, variational_flag=var_flag)
    enc_xena = generate_encoded_features(encoder, t_dataloaders, variational_flag=var_flag)
    enc_tot = pd.concat([enc_depmap, enc_xena], axis = 0, sort = False)
    enc_tot.to_csv(training_params['model_save_folder'] / 'encoded_features.csv')
    logger.info(f'Saved new encoded features')

    # plot umap
    umap_df = get_umap(encoded_features = enc_tot, 
                       depmap_sample_df = data_loader.depmap_sample_df, 
                       xena_sample_df = data_loader.xena_sample_df, 
                       save_folder = training_params['model_save_folder'] / 'plots', 
                       make_plot = True)
    logger.info(f'Saved UMAP plot')

if __name__ == '__main__':
    parser = argparse.ArgumentParser('training and evaluation')
    parser.add_argument('--gex_feature_file', dest='gex_feature_file', nargs='?')
    parser.add_argument('--samples', dest='samples', nargs='?', default='',
                       choices=['', 'tcgaonly'])
    parser.add_argument('--ngene', dest='ngene', nargs='?', default='all')
    parser.add_argument('--folder', dest='folder', nargs='?', default='.')
    norm_feat_group = parser.add_mutually_exclusive_group(required=False)
    norm_feat_group.add_argument('--norm_feat', dest='norm_feat', action='store_true')
    norm_feat_group.add_argument('--no-norm_feat', dest='norm_feat', action='store_false')
    parser.set_defaults(norm_feat=True)
    shared_group = parser.add_mutually_exclusive_group(required=False)
    shared_group.add_argument('--only_shared', dest='only_shared', action='store_true')
    shared_group.add_argument('--no-only_shared', dest='only_shared', action='store_false')
    parser.set_defaults(only_shared=False)
    parser.add_argument('--method', dest='method', nargs='?', default='ae_gan',
                        choices=['ae_gan', 'vae_gan'])
    
    args = parser.parse_args()
    # params_grid = {
    #     "pretrain_num_epochs": [0, 100, 300],
    #     "train_num_epochs": [100, 200, 300, 500, 750, 1000, 1500, 2000, 2500, 3000],
    #     "dop": [0.0, 0.1]
    # }
    params_grid = {
        "samples": [args.samples],
        "ngene": [args.ngene],
        "norm_feat_flag": [args.norm_feat], 
        "only_shared": [args.only_shared]
    }

    keys, values = zip(*params_grid.items())
    update_params_dict_list = [dict(zip(keys, v)) for v in itertools.product(*values)]

    for param_dict in update_params_dict_list:
        
        # logging
        logging.getLogger().handlers = [] # to stop logging on the file
        param_dict_log = param_dict.copy() # to avoid changing the original dict
        param_dict_log.update({"method": args.method})
        param_str_log = dict_to_str(param_dict_log)
        print(param_str_log)
        setup_logging(log_file=f"job_{os.environ.get('SLURM_JOBID')}_train_{param_str_log}.log")
        logger = get_logger()

        main(args=args, update_params_dict=param_dict)
