# compare results
library(umap)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(ggpubr)

create_list_folds <- function(exp_folder, beta_val = NULL){
  
  if(is.null(beta_val)){
    folder_models <- c(
      # all 
      sprintf('%ssamples__ngene_all_norm_feat_flag_False_only_shared_True/', exp_folder), 
      sprintf('%ssamples__ngene_all_norm_feat_flag_False_only_shared_False/', exp_folder), 
      sprintf('%ssamples__ngene_all_norm_feat_flag_True_only_shared_True/', exp_folder), 
      sprintf('%ssamples__ngene_all_norm_feat_flag_True_only_shared_False/', exp_folder), 
      # var1000
      sprintf('%ssamples__ngene_var1000_norm_feat_flag_False_only_shared_True/', exp_folder), 
      sprintf('%ssamples__ngene_var1000_norm_feat_flag_False_only_shared_False/', exp_folder), 
      sprintf('%ssamples__ngene_var1000_norm_feat_flag_True_only_shared_True/', exp_folder), 
      sprintf('%ssamples__ngene_var1000_norm_feat_flag_True_only_shared_False/', exp_folder), 
      # var5000
      sprintf('%ssamples__ngene_var5000_norm_feat_flag_False_only_shared_True/', exp_folder), 
      sprintf('%ssamples__ngene_var5000_norm_feat_flag_False_only_shared_False/', exp_folder), 
      sprintf('%ssamples__ngene_var5000_norm_feat_flag_True_only_shared_True/', exp_folder), 
      sprintf('%ssamples__ngene_var5000_norm_feat_flag_True_only_shared_False/', exp_folder),
      # all (tcga only)
      sprintf('%ssamples_tcgaonly_ngene_all_norm_feat_flag_False_only_shared_True/', exp_folder), 
      sprintf('%ssamples_tcgaonly_ngene_all_norm_feat_flag_False_only_shared_False/', exp_folder), 
      sprintf('%ssamples_tcgaonly_ngene_all_norm_feat_flag_True_only_shared_True/', exp_folder), 
      sprintf('%ssamples_tcgaonly_ngene_all_norm_feat_flag_True_only_shared_False/', exp_folder), 
      # var1000 (tcga only)
      sprintf('%ssamples_tcgaonly_ngene_var1000_norm_feat_flag_False_only_shared_True/', exp_folder), 
      sprintf('%ssamples_tcgaonly_ngene_var1000_norm_feat_flag_False_only_shared_False/', exp_folder), 
      sprintf('%ssamples_tcgaonly_ngene_var1000_norm_feat_flag_True_only_shared_True/', exp_folder), 
      sprintf('%ssamples_tcgaonly_ngene_var1000_norm_feat_flag_True_only_shared_False/', exp_folder), 
      # var5000 (tcga only)
      sprintf('%ssamples_tcgaonly_ngene_var5000_norm_feat_flag_False_only_shared_True/', exp_folder), 
      sprintf('%ssamples_tcgaonly_ngene_var5000_norm_feat_flag_False_only_shared_False/', exp_folder), 
      sprintf('%ssamples_tcgaonly_ngene_var5000_norm_feat_flag_True_only_shared_True/', exp_folder), 
      sprintf('%ssamples_tcgaonly_ngene_var5000_norm_feat_flag_True_only_shared_False/', exp_folder))
  }else{
    
    folder_models <- c(
      # all 
      sprintf('%ssamples__ngene_all_norm_feat_flag_False_only_shared_True_beta_%s/', exp_folder, beta_val), 
      sprintf('%ssamples__ngene_all_norm_feat_flag_False_only_shared_False_beta_%s/', exp_folder, beta_val), 
      sprintf('%ssamples__ngene_all_norm_feat_flag_True_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples__ngene_all_norm_feat_flag_True_only_shared_False_beta_%s/', exp_folder,beta_val), 
      # var1000
      sprintf('%ssamples__ngene_var1000_norm_feat_flag_False_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples__ngene_var1000_norm_feat_flag_False_only_shared_False_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples__ngene_var1000_norm_feat_flag_True_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples__ngene_var1000_norm_feat_flag_True_only_shared_False_beta_%s/', exp_folder,beta_val), 
      # var5000
      sprintf('%ssamples__ngene_var5000_norm_feat_flag_False_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples__ngene_var5000_norm_feat_flag_False_only_shared_False_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples__ngene_var5000_norm_feat_flag_True_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples__ngene_var5000_norm_feat_flag_True_only_shared_False_beta_%s/', exp_folder,beta_val), 
      # all (tcga only)
      sprintf('%ssamples_tcgaonly_ngene_all_norm_feat_flag_False_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_ngene_all_norm_feat_flag_False_only_shared_False_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_ngene_all_norm_feat_flag_True_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_ngene_all_norm_feat_flag_True_only_shared_False_beta_%s/', exp_folder,beta_val), 
      # var1000 (tcga only)
      sprintf('%ssamples_tcgaonly_ngene_var1000_norm_feat_flag_False_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_ngene_var1000_norm_feat_flag_False_only_shared_False_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_ngene_var1000_norm_feat_flag_True_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_ngene_var1000_norm_feat_flag_True_only_shared_False_beta_%s/', exp_folder,beta_val), 
      # var5000 (tcga only)
      sprintf('%ssamples_tcgaonly_ngene_var5000_norm_feat_flag_False_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_ngene_var5000_norm_feat_flag_False_only_shared_False_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_ngene_var5000_norm_feat_flag_True_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_ngene_var5000_norm_feat_flag_True_only_shared_False_beta_%s/', exp_folder, beta_val))
  }
  return(folder_models)
  
}


get_model_info <- function(folder_model){
  
  tmp <- str_split(folder_model, "/")[[1]]  
  model <- tmp[3]
  tmp1 <- str_split(tmp[4], "_")[[1]]
  samples <- ifelse(tmp1[2] == "", "all", tmp1[2])
  ngene <- tmp1[4]
  norm_feat <- tmp1[8]
  only_shared <- tmp1[11]
  
  return(data.frame(model = model, 
                    samples = samples, 
                    ngene = ngene, 
                    norm_feat = norm_feat, 
                    only_shared = only_shared))
}


# total results: 24 ae_gan + 24*3 vae_gan (beta=0.0005, beta=0.001, beta=0.0001)
### ae_gan ###
setwd("/Volumes/iorio/lucia/Multiomic_VAE/")
exp1 <- 'experiments/experiment_1/ae_gan/'
folder_models_exp1 <- create_list_folds(exp1)
exp2 <- 'experiments/experiment_2/vae_gan/'
folder_models_exp2_b0001 <- create_list_folds(exp2, beta_val = "0.0001")
folder_models_exp2_b0005 <- create_list_folds(exp2, beta_val = "0.0005")
folder_models_exp2_b001 <- create_list_folds(exp2, beta_val = "0.001")

df_auc_site <- list()
df_auc_study <- list()
df_auc_type <- list()
df_auc_site_p <- list()
df_auc_study_p <- list()
df_auc_type_p <- list()
for(i in 1:length(folder_models_exp1)){
  
  folder_model <- folder_models_exp1[i]
  info <- get_model_info(folder_model)
  
  df_auc_site[[i]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = NA) %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s", model, samples, ngene, norm_feat, only_shared))
  
  df_auc_site_p[[i]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = NA) %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s", model, samples, ngene, norm_feat, only_shared))
  
  df_auc_study[[i]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = NA) %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s", model, samples, ngene, norm_feat, only_shared))
  
  df_auc_study_p[[i]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = NA) %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s", model, samples, ngene, norm_feat, only_shared))
  
  df_auc_type[[i]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = NA) %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s", model, samples, ngene, norm_feat, only_shared))
  
  df_auc_type_p[[i]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = NA) %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s", model, samples, ngene, norm_feat, only_shared))
}


for(i in 1:length(folder_models_exp2_b0001)){
  
  folder_model <- folder_models_exp2_b0001[i]
  info <- get_model_info(folder_model)
  
  df_auc_site[[i+length(folder_models_exp1)]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0001"))
  
  
  df_auc_site_p[[i+length(folder_models_exp1)]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0001"))
  
  df_auc_study[[i+length(folder_models_exp1)]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0001"))
  
  df_auc_study_p[[i+length(folder_models_exp1)]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0001"))
  
  df_auc_type[[i+length(folder_models_exp1)]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0001"))
  
  df_auc_type_p[[i+length(folder_models_exp1)]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0001"))
}

n <- length(folder_models_exp1) + length(folder_models_exp2_b0001)
for(i in 1:length(folder_models_exp2_b0005)){
  
  folder_model <- folder_models_exp2_b0005[i]
  info <- get_model_info(folder_model)
  
  df_auc_site[[i+n]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0005") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0005"))
  
  df_auc_site_p[[i+n]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0005") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0005"))
  
  df_auc_study[[i+n]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0005") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0005"))
  
  df_auc_study_p[[i+n]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0005") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0005"))
  
  df_auc_type[[i+n]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0005") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0005"))
  
  df_auc_type_p[[i+n]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.0005") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.0005"))
}

n <- length(folder_models_exp1) + length(folder_models_exp2_b0001) + length(folder_models_exp2_b0005)
for(i in 1:length(folder_models_exp2_b001)){
  
  folder_model <- folder_models_exp2_b001[i]
  info <- get_model_info(folder_model)
  
  df_auc_site[[i+n]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.001"))
  
  df_auc_site_p[[i+n]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.001"))
  
  df_auc_study[[i+n]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.001"))
  
  df_auc_study_p[[i+n]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.001"))
  
  df_auc_type[[i+n]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.001"))
  
  df_auc_type_p[[i+n]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  ngene = info$ngene, 
                  norm_feat = info$norm_feat, 
                  only_shared = info$only_shared, 
                  beta = "0.001") %>%
    dplyr::mutate(tot_id = sprintf("%s_samples%s_ngene%s_normfeat%s_onlyshared%s_beta%s", 
                                   model, samples, ngene, norm_feat, only_shared, "0.001"))
}

df_auc_site <- do.call(rbind, df_auc_site)
df_auc_study <- do.call(rbind, df_auc_study)
df_auc_type <- do.call(rbind, df_auc_type)

df_auc_site_p <- do.call(rbind, df_auc_site_p)
df_auc_study_p <- do.call(rbind, df_auc_study_p)
df_auc_type_p <- do.call(rbind, df_auc_type_p)

ord_tot_id <- df_auc_site %>%
  dplyr::group_by(tot_id) %>%
  dplyr::summarize(median_auc = median(AUC)) %>%
  dplyr::arrange(desc(median_auc)) %>%
  dplyr::pull(tot_id)
pl1 <- ggplot(subset(df_auc_site), 
              aes(x = factor(tot_id, levels = ord_tot_id), y = AUC, fill = ngene)) + 
  geom_boxplot() + 
  facet_wrap(.~samples, nrow = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") + 
  ylab("AUC (Same tissue)")
pl1

ord_tot_id <- df_auc_site_p %>%
  dplyr::group_by(tot_id) %>%
  dplyr::summarize(median_auc = median(AUC)) %>%
  dplyr::arrange(desc(median_auc)) %>%
  dplyr::pull(tot_id)
pl1_p <- ggplot(subset(df_auc_site_p), 
              aes(x = factor(tot_id, levels = ord_tot_id), y = AUC, fill = ngene)) + 
  geom_boxplot() + 
  facet_wrap(.~samples, nrow = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") + 
  ylab("AUC (Same tissue)")
pl1_p

ord_tot_id <- df_auc_study %>%
  dplyr::mutate(auc_min05 = abs(AUC - 0.5)) %>%
  dplyr::group_by(tot_id) %>%
  dplyr::summarize(median_auc = median(auc_min05)) %>%
  dplyr::arrange((median_auc)) %>%
  dplyr::pull(tot_id)
pl2 <- ggplot(subset(df_auc_study), aes(x = factor(tot_id, levels = ord_tot_id), 
                                y = abs(AUC - 0.5), fill = ngene)) +
  geom_boxplot() + 
  facet_wrap(.~samples, nrow = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") + 
  ylab("AUC (Same study: DepMap, GTEx, TCGA, Target)")
pl2  

ord_tot_id <- df_auc_study_p %>%
  dplyr::mutate(auc_min05 = abs(AUC - 0.5)) %>%
  dplyr::group_by(tot_id) %>%
  dplyr::summarize(median_auc = median(auc_min05)) %>%
  dplyr::arrange((median_auc)) %>%
  dplyr::pull(tot_id)
pl2_p <- ggplot(subset(df_auc_study_p), aes(x = factor(tot_id, levels = ord_tot_id), 
                                        y = abs(AUC - 0.5), fill = ngene)) +
  geom_boxplot() + 
  facet_wrap(.~samples, nrow = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") + 
  ylab("AUC (Same study: DepMap, GTEx, TCGA, Target)")
pl2_p  

ord_tot_id <- df_auc_type %>%
  dplyr::group_by(tot_id) %>%
  dplyr::summarize(median_auc = median(AUC)) %>%
  dplyr::arrange(desc(median_auc)) %>%
  dplyr::pull(tot_id)
pl3 <- ggplot(subset(df_auc_type), aes(x = factor(tot_id, levels = ord_tot_id), 
                               y = AUC, fill = ngene)) +
  geom_boxplot() + 
  facet_wrap(.~samples, nrow = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") +
  ylab("AUC (Same class: Solid T, Blood T. or Non-Cancerous)")
pl3

ord_tot_id <- df_auc_type_p %>%
  dplyr::group_by(tot_id) %>%
  dplyr::summarize(median_auc = median(AUC)) %>%
  dplyr::arrange(desc(median_auc)) %>%
  dplyr::pull(tot_id)
pl3_p <- ggplot(subset(df_auc_type_p), aes(x = factor(tot_id, levels = ord_tot_id), 
                                       y = AUC, fill = ngene)) +
  geom_boxplot() + 
  facet_wrap(.~samples, nrow = 2, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") +
  ylab("AUC (Same class: Solid T, Blood T. or Non-Cancerous)")
pl3_p

### get medians 
# create summary medians 
summary_auc <- function(df_auc_site, df_auc_study, df_auc_type, samples_class = "tcgaonly") {
  df_auc_site_medians <- df_auc_site %>%
    dplyr::filter(samples == samples_class) %>%
    dplyr::group_by(tot_id) %>%
    dplyr::summarize(median_auc = median(AUC)) %>%
    dplyr::arrange(desc(median_auc)) %>%
    dplyr::mutate(tot_id = factor(tot_id, levels = .$tot_id)) %>%
    dplyr::select(tot_id, median_auc) %>%
    dplyr::mutate(type = "Same tissue") %>%
    dplyr::rename(median_auc_site = median_auc) %>%
    dplyr::mutate(order = 1:length(tot_id))
  
  df_auc_study_medians <- df_auc_study %>%
    dplyr::filter(samples == samples_class) %>%
    dplyr::group_by(tot_id) %>%
    dplyr::mutate(auc_min05 = abs(AUC - 0.5)) %>%
    dplyr::summarize(median_auc = median(auc_min05)) %>%
    dplyr::arrange(median_auc) %>%
    dplyr::mutate(tot_id = factor(tot_id, levels = .$tot_id)) %>%
    dplyr::select(tot_id, median_auc) %>%
    dplyr::mutate(type = "Same study") %>%
    dplyr::rename(median_auc_study = median_auc) %>%
    dplyr::mutate(order = 1:length(tot_id))
  
  df_auc_type_medians <- df_auc_type %>%
    dplyr::filter(samples == samples_class) %>%
    dplyr::group_by(tot_id) %>%
    dplyr::summarize(median_auc = median(AUC)) %>%
    dplyr::arrange(desc(median_auc)) %>%
    dplyr::mutate(tot_id = factor(tot_id, levels = .$tot_id)) %>%
    dplyr::select(tot_id, median_auc) %>%
    dplyr::mutate(type = "Same class") %>%
    dplyr::rename(median_auc_type = median_auc)  %>%
    dplyr::mutate(order = 1:length(tot_id))
  
  df_auc_medians <- full_join(df_auc_site_medians, df_auc_study_medians, by = "tot_id", 
                              suffix = c("_site", "_study"))
  df_auc_medians <- full_join(df_auc_medians, df_auc_type_medians, by = "tot_id", 
                              suffix = c("", "_type")) %>%
    mutate(order_sum = order_site + order_study + order) %>%
    arrange(order_sum)
  colnames(df_auc_medians) <- c("tot_id", "median_auc_site", "type_site", "order_site", 
                                "median_auc_study", "type_study", "order_study", 
                                "median_auc_type", "type_type", "order_type", "order_sum")
  
  return(df_auc_medians)
}
df_auc_medians_tcgaonly <- summary_auc(df_auc_site, df_auc_study, df_auc_type, samples_class = "tcgaonly")

#View(df_auc_medians)
write.table(df_auc_medians_tcgaonly, 
            file = "experiments/auc_medians_tcgaonly_ex1_ex2.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

#df_auc_medians_tcga = read.table("experiments/auc_medians_tcgaonly_ex1_ex2.txt", header = TRUE, sep = "\t")
df_auc_medians_all <- summary_auc(df_auc_site, df_auc_study, df_auc_type, samples_class = "all")
write.table(df_auc_medians_all, 
            file = "experiments/auc_medians_all_ex1_ex2.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

df_auc_p_medians_tcgaonly <- summary_auc(df_auc_site_p, df_auc_study_p, df_auc_type_p, samples_class = "tcgaonly")
write.table(df_auc_medians_all, 
            file = "experiments/auc_medians_private_tcgaonly_ex1_ex2.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

df_auc_p_medians_all <- summary_auc(df_auc_site_p, df_auc_study_p, df_auc_type_p, samples_class = "all")
write.table(df_auc_medians_all, 
            file = "experiments/auc_medians_private_all_ex1_ex2.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

