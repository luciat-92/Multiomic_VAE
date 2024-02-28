import sys, os
sys.path.insert(0, os.getcwd())
import gzip
import pandas as pd
import numpy as np
import re
from utils.utils import *
from utils.config import *
from utils.logger import setup_logging, get_logger
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go

def harmonize_names(df, var, mapping_dict):
    df[var] = df[var].replace(mapping_dict)
    return df

# filter based on most variable genes
def get_most_variable_genes(df, samples = None, n_genes = 1000):
    # filter based on samples
    if samples is not None:
        df = df.loc[samples, :]
    std_top = df.std(axis=0).sort_values(ascending=False)
    # get the top n genes
    std_top = std_top.iloc[0:n_genes]
    # get the gene names
    gene_names = std_top.index.tolist()
    return gene_names

def filter_gene_expression(df, samples = None, n_genes = 1000):
    # samples is a dictorionary with keys: 'depmap' and 'xena'
    if samples is None:
        genes = get_most_variable_genes(df, n_genes = n_genes)
    else:
        genes = []
        for key in samples.keys():
            samples_id = samples[key]
            genes_tmp = get_most_variable_genes(df, samples = samples_id, n_genes = n_genes)
            genes.extend(genes_tmp)
        genes = list(set(genes))
    logger.info(f'Number of common genes: {len(genes)} from top {n_genes} most variable genes')
    # filter based on genes
    df = df.loc[:, genes]
    return df

# create a function to create standard deviation and mean per gene (colum)
def get_mean_std_cv(df):
    std = df.std(axis=0)
    mean = df.mean(axis=0)
    df_summ = pd.DataFrame({'id': df.columns, 'mean': mean, 'std': std})
    df_summ['cv'] =  df_summ['std'] /  df_summ['mean']
    # remove index
    df_summ = df_summ.reset_index(drop=True)
    return df_summ 

def generate_sankey_chart_data(df, columns, sankey_link_weight):
    # list of list: each list are the set of nodes in each tier/column
    column_values = [df[col] for col in columns]
    # this generates the labels for the sankey by taking all the unique values
    labels = sum([list(node_values.unique()) for node_values in column_values],[])
    # initializes a dict of dicts (one dict per tier) 
    link_mappings = {col: {} for col in columns}

    # each dict maps a node to unique number value (same node in different tiers
    # will have different nubmer values
    i = 0
    for col, nodes in zip(columns, column_values):
        for node in nodes.unique():
            link_mappings[col][node] = i
            i = i + 1
    # specifying which coluns are serving as sources and which as sources
    # ie: given 3 df columns (col1 is a source to col2, col2 is target to col1 and 
    # a source to col 3 and col3 is a target to col2
    source_nodes = column_values[: len(columns) - 1]
    target_nodes = column_values[1:]
    source_cols = columns[: len(columns) - 1]
    target_cols = columns[1:]
    links = []
    # loop to create a list of links in the format [((src,tgt),wt),(),()...]
    for source, target, source_col, target_col in zip(source_nodes, target_nodes, source_cols, target_cols):
        for val1, val2, link_weight in zip(source, target, df[sankey_link_weight]):
            links.append(((link_mappings[source_col][val1],link_mappings[target_col][val2]),link_weight,))

    # creating a dataframe with 2 columns: for the links (src, tgt) and weights
    df_links = pd.DataFrame(links, columns=["link", "weight"])

    # aggregating the same links into a single link (by weight)
    df_links = df_links.groupby(by=["link"], as_index=False).agg({"weight": sum})

    # generating three lists needed for the sankey visual
    sources = [val[0] for val in df_links["link"]]
    targets = [val[1] for val in df_links["link"]]
    weights = df_links["weight"]

    return labels, sources, targets, weights

# plot sunkey chart
def plot_sankey_chart(link_sankey_plot, file_name):
    # make plot
    node = dict(label = link_sankey_plot[0], pad = 50, thickness = 20)
    link = dict(source = link_sankey_plot[1], target = link_sankey_plot[2], value = link_sankey_plot[3])
    data_plot = go.Sankey(link=link, node=node)
    fig = go.Figure(data_plot)
    fig.update_layout(
        autosize=False,
        width=1200,
        height=1000, 
        font_size=15
    )
    safe_create_dir(file_name.parent)
    fig.write_image(file_name.with_suffix(".pdf"))

# make plot of distributions
def plot_summary_CL_vs_tissue(depmap_gene_df_match, xena_gene_df_match, save_folder):
    # plot the distribution of the coefficient of variation
    df_summary = pd.merge(
        get_mean_std_cv(depmap_gene_df_match), 
        get_mean_std_cv(xena_gene_df_match), 
        on = 'id', 
        suffixes=('_depmap', '_xena'))
    # change the size of the plot
    fig_size = (12, 4) 
    fig, axes = plt.subplots(1,3, figsize=fig_size)
    # plt.rcParams['figure.figsize'] = [14, 4]
    sns.scatterplot(data = df_summary, x = "cv_depmap", y = "cv_xena", size = 0.5, alpha = 0.5, ax = axes[0])
    sns.scatterplot(data = df_summary, x = "mean_depmap", y = "mean_xena", size = 0.5, alpha = 0.5, ax = axes[1])
    sns.scatterplot(data = df_summary, x = "std_depmap", y = "std_xena", size = 0.5, alpha = 0.5, ax = axes[2])
    # change labels
    axes[0].set_xlabel("Cell lines")
    axes[0].set_ylabel("Solid tissues")
    axes[1].set_xlabel("Cell lines")
    axes[1].set_ylabel("Solid tissues")
    axes[2].set_xlabel("Cell lines")
    axes[2].set_ylabel("Solid tissues")
    # remove legend
    axes[0].get_legend().remove()
    axes[1].get_legend().remove()
    axes[2].get_legend().remove()
    # title
    axes[0].set_title("Coefficient of variation")
    axes[1].set_title("Mean")
    axes[2].set_title("Standard deviation")
    plt.tight_layout()
    safe_create_dir(save_folder)
    plt.savefig(save_folder / 'GEX_compare_distributions.pdf', format = 'pdf')


def load_data(load_gex = True):
    XENA_FOLDER = PERSONAL_DATA_DIR / 'XENA/TCGA_TARGET_GTEx/' 
    XENA_GENE_FILE = XENA_FOLDER / 'TcgaTargetGtex_rsem_gene_tpm'
    # XENA_GENE_FILE = XENA_FOLDER / 'TEST_rsem_gene_tpm'
    XENA_GENE_MAPPING_FILE = XENA_FOLDER / 'probeMap%2Fgencode.v23.annotation.gene.probemap'
    XENA_SAMPLE_FILE = XENA_FOLDER / 'TcgaTargetGTEX_phenotype.txt.gz'
    DEPMAP_FOLDER = PERSONAL_DATA_DIR / 'DEPMAP_PORTAL/version_23Q2/'
    DEPMAP_GENE_FILE = DEPMAP_FOLDER / 'OmicsExpressionProteinCodingGenesTPMLogp1.csv'  
    DEPMAP_SAMPLE_FILE = DEPMAP_FOLDER / 'Model.csv'
    HGNC_GENE_ANN_FILE = PERSONAL_DATA_DIR / 'hgnc_complete_set_202311.txt'
    ## load data ##
    if load_gex:
        # gene annotation
        hgnc_gene_ann_df = pd.read_csv(HGNC_GENE_ANN_FILE, sep="\t", low_memory=False)
        # depmap data
        depmap_gene_df = pd.read_csv(DEPMAP_GENE_FILE, sep=',', index_col=0)
        depmap_gene_df.index.name = 'sample_id'
        # xena data
        xena_gene_mapping_df = pd.read_csv(XENA_GENE_MAPPING_FILE, sep='\t')
        xena_gene_df = pd.read_csv(XENA_GENE_FILE, sep='\t')    
        xena_gene_df.index = xena_gene_df['sample']
        xena_gene_df = xena_gene_df.drop(columns=['sample'])
        xena_gene_df.index.name = 'xena_gene_id'
    else:
        hgnc_gene_ann_df = None
        depmap_gene_df = None
        xena_gene_mapping_df = None
        xena_gene_df = None
    depmap_sample_df = pd.read_csv(DEPMAP_SAMPLE_FILE, sep=',', index_col=0)
    depmap_sample_df.index.name = 'sample_id'
    # load tissue data
    with gzip.open(XENA_SAMPLE_FILE) as f:
        xena_sample_df = pd.read_csv(f, sep='\t', index_col=0, encoding='ISO-8859-1')
    xena_sample_df.index.name = 'sample_id'
    # correct site and type annotation
    type_harmonization_names = {
        "Solid Tissue Normal": "Normal Tissue", 
        "Primary Solid Tumor": "Primary Tumor", 
        "Recurrent Solid Tumor": "Recurrent Tumor",
        "Additional - New Primary": "Primary Tumor", 
        "Additional Metastatic": "Metastatic", 
        "Primary Blood Derived Cancer - Peripheral Blood": "Primary Blood Derived Cancer", 
        "Primary Blood Derived Cancer - Bone Marrow": "Primary Blood Derived Cancer", 
        "Recurrent Blood Derived Cancer - Bone Marrow": "Recurrent Blood Derived Cancer",
        "Recurrent Blood Derived Cancer - Peripheral Blood": "Recurrent Blood Derived Cancer", 
        "Post treatment Blood Cancer - Bone Marrow": "Post treatment Blood Cancer", 
        "Post treatment Blood Cancer - Peripheral Blood": "Post treatment Blood Cancer", 
        "Post treatment Blood Cancer - Blood": "Post treatment Blood Cancer"
    }
    site_harmonization_names = {
        "Vagina":"Vulva/Vagina", 
        "Thyroid Gland": "Thyroid", 
        "SympatheticÊNervous System": "Peripheral Nervous System", 
        "Stomach" : "Esophagus/Stomach",
        "Soft tissue,Bone" : "Soft Tissue", 
        "Head and Neck region" : "Head and Neck", 
        "Paraganglia" : "Peripheral Nervous System", 
        "Ovary" : "Ovary/Fallopian Tube", 
        "Lymphoid" : "Immune cells", 
        "Myeloid" : "Immune cells", 
        "White blood cell" : "Immune cells",
        "Fallopian Tube" : "Ovary/Fallopian Tube", 
        "Esophagus" : "Esophagus/Stomach", 
        "Cervix Uteri" : "Cervix", 
        "Brain" : "CNS/Brain", 
        "Adrenal gland" : "Adrenal Gland",
        "Bladder": "Bladder/Urinary Tract", 
        "Bile duct": "Biliary Tract", 
        "Rectum" : "Colon/Rectum",
        "Colon" : "Colon/Rectum",
        "Bowel" : "Colon/Rectum"
    }
    depmap_sample_df = harmonize_names(depmap_sample_df, 'OncotreeLineage', site_harmonization_names)
    xena_sample_df = harmonize_names(xena_sample_df, '_primary_site', site_harmonization_names)
    xena_sample_df = harmonize_names(xena_sample_df, '_sample_type', type_harmonization_names)
    return hgnc_gene_ann_df, depmap_gene_df, depmap_sample_df, xena_sample_df, xena_gene_mapping_df, xena_gene_df


def match_and_filter_genes(depmap_gene_df, hgnc_gene_ann_df, xena_gene_mapping_df, xena_gene_df):
    # match data frames by gene id (use hgnc annotation)
    # to get the new gene names, assign to a dataframe and match with the two available annotations!
    gene_names = []
    gene_ids = []
    # Split each element in the array and extract gene name and ID
    for column_name in depmap_gene_df.columns:
        parts = re.split(r'\s*\(|\)\s*', column_name)
        gene_name = parts[0].strip()
        gene_id = parts[1].strip()
        gene_names.append(gene_name)
        gene_ids.append(gene_id)
    # Create a DataFrame with the gene names and gene IDs as columns
    df = pd.DataFrame({'gene_name': gene_names, 'entrez_id': gene_ids})
    # match with hgnc
    df['entrez_id'] = df['entrez_id'].astype('float64')
    match_df = pd.merge(df, 
                    hgnc_gene_ann_df[["hgnc_id", "entrez_id", "symbol", "ensembl_gene_id"]], 
                    on = 'entrez_id')
    # match with xena gene mapping
    xena_gene_mapping_df['ensembl_gene_id'] = xena_gene_mapping_df['id'].str.split('.').str[0]
    match_df = pd.merge(match_df, xena_gene_mapping_df, on = 'ensembl_gene_id')
    # remvoe duplicates and modify names
    match_df = match_df.drop_duplicates(subset=["ensembl_gene_id"])
    match_df['depmap_gene_id'] = match_df['gene_name'] + " (" + match_df['entrez_id'].astype('int').astype('str') + ")"
    match_df['xena_gene_id'] = match_df['id']
    match_df['complete_id'] = match_df['symbol'] + " (" + match_df['ensembl_gene_id'] + ")"
    # create a new data frame from depmap_gene_df, keep only those genes that have a value in depmap_id in match_df
    depmap_gene_df_match = depmap_gene_df[match_df['depmap_gene_id']]
    depmap_gene_df_match.columns = match_df['complete_id']
    depmap_gene_df_match.columns.name = None
    # create a new data frame from xena_gene_df, keep only those genes that have a value in xena_id in match_df
    xena_gene_df_match = xena_gene_df.loc[match_df['xena_gene_id'], :]
    # transpose xena_gene_df_match
    xena_gene_df_match = xena_gene_df_match.transpose()
    xena_gene_df_match.index.name = 'sample_id'
    xena_gene_df_match.columns = match_df['complete_id']
    xena_gene_df_match.columns.name = None
    # convert entries in xena_gene_df_match: computed as log2(x+0.001) but should be log2(x+1)
    xena_gene_df_match = np.log2(np.power(2, xena_gene_df_match) - 0.001 + 1)
    xena_gene_df_match[xena_gene_df_match.abs() < 10**-7] = 0
    # create final table (samples x genes)
    tot_gene_df = pd.concat([depmap_gene_df_match, xena_gene_df_match], axis=0, sort=False)
    #### PLOT ####
    plot_summary_CL_vs_tissue(depmap_gene_df_match, xena_gene_df_match, DATA_DIR / 'preprocessed/plots/')
    return tot_gene_df

def main():
    # load
    hgnc_gene_ann_df, depmap_gene_df, depmap_sample_df, xena_sample_df, xena_gene_mapping_df, xena_gene_df = load_data()
    logger.info("All data loaded")
    
    # filter for common genes
    tot_gene_df = match_and_filter_genes(depmap_gene_df, hgnc_gene_ann_df, xena_gene_mapping_df, xena_gene_df)
    logger.info("Genes matched and filtered")

    # save based on variability
    depmap_samples = depmap_sample_df.index.intersection(tot_gene_df.index)
    xena_samples = xena_sample_df.index.intersection(tot_gene_df.index)
    tcga_samples = xena_sample_df[xena_sample_df["_study"] == "TCGA"].index.intersection(tot_gene_df.index)
    depmap_sample_df = depmap_sample_df.loc[depmap_samples]
    xena_sample_df = xena_sample_df.loc[xena_samples]
    tcga_samples_df = xena_sample_df.loc[tcga_samples]

    # most 1000 variable genes
    logger.info("Filtering based on variability: top 1000 genes, all samples")
    var1000_gene_df = filter_gene_expression(
        df = tot_gene_df, samples = {'depmap': depmap_samples, 'xena': xena_samples}, 
        n_genes = 1000)
    var1000_gene_df.to_csv( DATA_DIR / 'preprocessed/gene_expression_var1000.csv.gz', compression='gzip')

    # most 5000 variable genes
    logger.info("Filtering based on variability: top 5000 genes, all samples")
    var5000_gene_df = filter_gene_expression(
        df = tot_gene_df, samples = {'depmap': depmap_samples, 'xena': xena_samples}, 
        n_genes = 5000)
    var5000_gene_df.to_csv(DATA_DIR / 'preprocessed/gene_expression_var5000.csv.gz', compression='gzip')

    # most 1000 variable genes using tcga
    logger.info("Filtering based on variability: top 1000 genes, only tcga")
    var1000_tcgaonly_gene_df = filter_gene_expression(
        df = tot_gene_df, samples = {'depmap': depmap_samples, 'tcga': tcga_samples}, 
        n_genes = 1000)
    var1000_tcgaonly_gene_df.to_csv(DATA_DIR / 'preprocessed/gene_expression_var1000_tcgaonly.csv.gz', compression='gzip')

    # most 5000 variable genes using tcga
    logger.info("Filtering based on variability: top 5000 genes, only tcga")
    var5000_tcgaonly_gene_df = filter_gene_expression(
        df = tot_gene_df, samples = {'depmap': depmap_samples, 'tcga': tcga_samples}, 
        n_genes = 5000)
    var5000_tcgaonly_gene_df.to_csv(DATA_DIR / 'preprocessed/gene_expression_var5000_tcgaonly.csv.gz', compression='gzip')

    # save including rownames and colnames, compress
    logger.info("Save all genes, all samples")
    tot_gene_df.to_csv(DATA_DIR / 'preprocessed/gene_expression_all.csv.gz', compression='gzip')
    logger.info("Combined matrices saved")

    # plot samples info (sankey chart)
    tmp = xena_sample_df
    tmp["weight"] = 1
    link_sankey_plot = generate_sankey_chart_data(tmp, ['_gender', '_study', '_sample_type', '_primary_site'], 'weight')
    plot_sankey_chart(link_sankey_plot, DATA_DIR / 'preprocessed/plots/samples_TCGA_TARGET_GTEX')

    tmp = depmap_sample_df
    tmp["weight"] = 1
    link_sankey_plot = generate_sankey_chart_data(tmp, ['Sex', 'GrowthPattern', 'PrimaryOrMetastasis', 'OncotreeLineage'], 'weight')
    plot_sankey_chart(link_sankey_plot, DATA_DIR / 'preprocessed/plots/samples_DEPMAP')
    logger.info("Plots: samples info")

if __name__ == "__main__":
    logging.getLogger().handlers = [] # to stop logging on the file
    setup_logging(log_file=f"job_{os.environ.get('SLURM_JOBID')}_gexpreproc.log")
    logger = get_logger()
    main()
