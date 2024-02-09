import umap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import pandas as pd
from utils.utils import *
from utils.logger import setup_logging, get_logger


def plot_training_loss(minibatch_losses, 
                       num_epochs, 
                       save_folder = None,
                       averaging_iterations=100, 
                       custom_label=''):
    
    logger = get_logger()
    iter_per_epoch = len(minibatch_losses) // num_epochs
    if len(minibatch_losses) < 100:
            num_losses = len(minibatch_losses) // 2
    else:
        num_losses = 100

    try:
        plt.figure(figsize=(7, 5))
        ax1 = plt.subplot(1, 1, 1)
        ax1.plot(range(len(minibatch_losses)),
                 (minibatch_losses), label=f'Minibatch Loss')
        ax1.set_xlabel('Iterations')
        ax1.set_ylabel('Loss')

        ax1.set_ylim([
            min(0, np.min(minibatch_losses[num_losses:])*1.01), np.max(minibatch_losses[num_losses:])*1.01
            ])
        ax1.plot(np.convolve(minibatch_losses,
                             np.ones(averaging_iterations,)/averaging_iterations,
                             mode='valid'),
                label='Running Average')
        # add title to ax1
        ax1.set_title(f'{custom_label}')
        ax1.legend()

        ###################
        # Set scond x-axis
        ax2 = ax1.twiny()
        newlabel = list(range(num_epochs+1))
        newpos = [e*iter_per_epoch for e in newlabel]

        ax2.set_xticks(newpos[::100])
        ax2.set_xticklabels(newlabel[::100])

        ax2.xaxis.set_ticks_position('bottom')
        ax2.xaxis.set_label_position('bottom')
        ax2.spines['bottom'].set_position(('outward', 45))
        ax2.set_xlabel('Epochs')
        ax2.set_xlim(ax1.get_xlim())
        ###################
        plt.tight_layout()
        # remove space in custum_label string and substitute with _
        custom_label = custom_label.replace(' ', '_')
        if save_folder is not None:
            plt.savefig(save_folder / f'{custom_label}.pdf', format = 'pdf', bbox_inches='tight')
    except:
        logger.info(f'Error in plotting loss {custom_label}')

def plot_histories_gan(histories, 
                   save_folder, 
                   pretrain_num_epochs, 
                   train_num_epochs, 
                   variational_flag = False):
    
    safe_create_dir(save_folder)
    
    # plot error pretrain adversarial network
    if variational_flag:
        name_title_loss = "Recons Ortho and KL Div"
        plot_training_loss(histories[0]['kl_div'], pretrain_num_epochs, save_folder = save_folder, custom_label="KL Div (Pretrain)")
    else:
        name_title_loss = "Recons and Ortho"

    plot_training_loss(histories[0]['loss'], pretrain_num_epochs, save_folder = save_folder, custom_label=f'{name_title_loss} (Pretrain)')
    plot_training_loss(histories[0]['recons_loss'], pretrain_num_epochs, save_folder = save_folder, custom_label="Recons (Pretrain)")
    plot_training_loss(histories[0]['ortho_loss'], pretrain_num_epochs, save_folder = save_folder, custom_label="Ortho (Pretrain)")

    # plot error train adversarial network
    if variational_flag:
        plot_training_loss(histories[3]['kl_div'], train_num_epochs, save_folder = save_folder, custom_label="KL Div (Train GAN)")
    
    plot_training_loss(histories[2]['critic_loss'], train_num_epochs, save_folder = save_folder, custom_label="Critic (Train GAN)")    
    histories[3]['tot_loss'] = [x + y for x, y in zip(histories[3]['gen_loss'], histories[3]['loss'])]
    plot_training_loss(histories[3]['loss'], train_num_epochs, save_folder = save_folder, custom_label=f'{name_title_loss} (Train GAN)')
    plot_training_loss(histories[3]['recons_loss'], train_num_epochs, save_folder = save_folder, custom_label="Recons (Train GAN)")
    plot_training_loss(histories[3]['ortho_loss'], train_num_epochs, save_folder = save_folder, custom_label="Ortho (Train GAN)")
    plot_training_loss(histories[3]['gen_loss'], train_num_epochs, save_folder = save_folder, custom_label="Generative (Train GAN)")
    plot_training_loss(histories[3]['tot_loss'], train_num_epochs, save_folder = save_folder, custom_label="Total (Train GAN)") 



def plot_umap(umap_df, save_folder):
    """
    The function `plot_umap` takes a DataFrame `umap_df` and a save folder path `save_folder` as input,
    and plots two UMAP scatterplots based on different columns of the DataFrame, saving the plots as PDF
    files in the specified folder.
    
    :param umap_df: The `umap_df` parameter is a DataFrame that contains the data to be plotted on the
    UMAP plot. It should have the following columns: umap_1, umap_2, site, type, and study. 
    :param save_folder: The `save_folder` parameter is the path to the folder where you want to save the
    generated plots and umap in csv format.
    """

    # plot
    plt.figure(figsize=(5, 5))
    sns.scatterplot(data=umap_df ,x='umap_1', y='umap_2', hue='site', style='type', markers = ['o', '.'], alpha = 0.4, linewidth=0.1)
    #Set plot title and labels
    plt.xlabel("UMAP Dimension 1", fontsize=10)
    plt.ylabel("UMAP Dimension 2", fontsize=10)
    plt.legend(fontsize=6, title_fontsize=6, loc="upper left", bbox_to_anchor=(1.05, 1))
    # save plot
    plt.savefig(save_folder / 'umap_lineage.pdf', format = 'pdf', bbox_inches='tight')

    plt.figure(figsize=(5, 5))
    sns.scatterplot(data=umap_df, x='umap_1', y='umap_2', hue='study', style='type', markers = ['o', '.'], alpha = 0.4, linewidth=0.1)
    # Set plot title and labels
    plt.xlabel("UMAP Dimension 1", fontsize=10)
    plt.ylabel("UMAP Dimension 2", fontsize=10)
    plt.legend(fontsize=6, title_fontsize=6, loc="upper left", bbox_to_anchor=(1.05, 1))
    plt.savefig(save_folder / 'umap_study.pdf', format = 'pdf', bbox_inches='tight')

def get_umap(encoded_features, depmap_sample_df, xena_sample_df, save_folder, make_plot = True):

    depmap_sample_df['type'] = 'depmap'
    depmap_sample_df['study'] = 'CL_depmap'
    xena_sample_df['type'] = 'xena'
    t1 = depmap_sample_df[["type", "Sex", "OncotreeLineage", "study"]].rename(columns={"Sex":"sex", "OncotreeLineage":"site"})
    t2 = xena_sample_df[["type", "_primary_site", "_gender", "_study", "_sample_type"]].rename(columns={"_primary_site": "site", "_gender":"sex", "_study":"study"})
    tot_sample_df = pd.concat([t1, t2], axis=0, sort=False)
    # order tot_sample_df according tmp
    tot_sample_df = tot_sample_df.loc[encoded_features.index]
    embedding = umap.UMAP().fit_transform(encoded_features)
    umap_df = pd.DataFrame(embedding, index = encoded_features.index, columns = ['umap_1', 'umap_2'])
    umap_df = pd.concat([umap_df, tot_sample_df], axis = 1)
    # save
    safe_create_dir(save_folder)
    umap_df.to_csv(save_folder / 'umap.csv')
    # plot
    if make_plot:
        plot_umap(umap_df = umap_df, save_folder = save_folder)

    return umap_df



