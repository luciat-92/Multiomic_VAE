# compare results
library(umap)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(ggpubr)

create_list_folds <- function(exp_folder, beta_val = NULL){
  
    folder_models <- c(
      # var1000 (tcga only)
      sprintf('%ssamples_tcgaonly_nfeat_var1000_norm_feat_flag_False_norm_type_zscore_only_shared_True_beta_%s/', exp_folder,beta_val),
      sprintf('%ssamples_tcgaonly_nfeat_var1000_norm_feat_flag_True_norm_type_minmax_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_nfeat_var1000_norm_feat_flag_True_norm_type_zscore_only_shared_True_beta_%s/', exp_folder,beta_val),
      sprintf('%ssamples_tcgaonly_nfeat_var1000_norm_feat_flag_False_norm_type_zscore_only_shared_False_beta_%s/', exp_folder,beta_val),
      sprintf('%ssamples_tcgaonly_nfeat_var1000_norm_feat_flag_True_norm_type_minmax_only_shared_False_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_nfeat_var1000_norm_feat_flag_True_norm_type_zscore_only_shared_False_beta_%s/', exp_folder,beta_val),
      # var5000 (tcga only)
      sprintf('%ssamples_tcgaonly_nfeat_var5000_norm_feat_flag_False_norm_type_zscore_only_shared_True_beta_%s/', exp_folder,beta_val),
      sprintf('%ssamples_tcgaonly_nfeat_var5000_norm_feat_flag_True_norm_type_minmax_only_shared_True_beta_%s/', exp_folder,beta_val), 
      sprintf('%ssamples_tcgaonly_nfeat_var5000_norm_feat_flag_True_norm_type_zscore_only_shared_True_beta_%s/', exp_folder,beta_val), 
    sprintf('%ssamples_tcgaonly_nfeat_var5000_norm_feat_flag_False_norm_type_zscore_only_shared_False_beta_%s/', exp_folder,beta_val),
    sprintf('%ssamples_tcgaonly_nfeat_var5000_norm_feat_flag_True_norm_type_minmax_only_shared_False_beta_%s/', exp_folder,beta_val), 
    sprintf('%ssamples_tcgaonly_nfeat_var5000_norm_feat_flag_True_norm_type_zscore_only_shared_False_beta_%s/', exp_folder,beta_val))
  return(folder_models)
}


get_model_info <- function(folder_model){
  
  tmp <- str_split(folder_model, "/")[[1]]  
  model <- tmp[3]
  tmp1 <- str_split(tmp[4], "_")[[1]]
  samples <- ifelse(tmp1[2] == "", "all", tmp1[2])
  nfeat <- tmp1[4]
  norm_feat <- tmp1[8]
  norm_type <- ifelse(norm_feat == "False", "none", tmp1[11])
  only_shared <- tmp1[14]
  beta <- tmp1[16]
  
  return(data.frame(model = model, 
                    samples = samples, 
                    nfeat = nfeat, 
                    norm_feat = norm_feat, 
                    norm_type = norm_type,
                    only_shared = only_shared, 
                    beta = beta))
}


# total results: 6*3 mvae_gan (beta=0.0005, beta=0.001, beta=0.0001)
### mvae_gan ###
setwd("/Volumes/iorio/lucia/Multiomic_VAE/")
exp3 <- 'experiments/experiment_3/mvae_gan/'
folder_models_exp3 <- c(create_list_folds(exp3, beta_val = "0.0001"), 
                              create_list_folds(exp3, beta_val = "0.0005"), 
                              create_list_folds(exp3, beta_val = "0.001"))

df_auc_site <- list()
df_auc_study <- list()
df_auc_type <- list()
df_auc_site_p <- list()
df_auc_study_p <- list()
df_auc_type_p <- list()
for(i in 1:length(folder_models_exp3)){
  
  folder_model <- folder_models_exp3[i]
  info <- get_model_info(folder_model)
  
  df_auc_site[[i]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  nfeat = info$nfeat, 
                  norm_feat = info$norm_feat, 
                  norm_type = info$norm_type,
                  only_shared = info$only_shared, 
                  beta = info$beta) %>%
    dplyr::mutate(tot_id = sprintf("%s_nfeat%s_normfeat%s_normtype%s_onlyshared%s_beta%s", model, nfeat, norm_feat, norm_type, only_shared, beta))
  
  df_auc_site_p[[i]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_tissue.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  nfeat = info$nfeat, 
                  norm_feat = info$norm_feat, 
                  norm_type = info$norm_type,
                  only_shared = info$only_shared, 
                  beta = info$beta) %>%
    dplyr::mutate(tot_id = sprintf("%s_nfeat%s_normfeat%s_normtype%s_onlyshared%s_beta%s", model, nfeat, norm_feat, norm_type, only_shared, beta))
  
  df_auc_study[[i]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  nfeat = info$nfeat, 
                  norm_feat = info$norm_feat, 
                  norm_type = info$norm_type,
                  only_shared = info$only_shared, 
                  beta = info$beta) %>%
    dplyr::mutate(tot_id = sprintf("%s_nfeat%s_normfeat%s_normtype%s_onlyshared%s_beta%s", model, nfeat, norm_feat, norm_type, only_shared,  beta))
  
  df_auc_study_p[[i]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_study.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  nfeat = info$nfeat, 
                  norm_feat = info$norm_feat, 
                  norm_type = info$norm_type,
                  only_shared = info$only_shared, 
                  beta = info$beta) %>%
    dplyr::mutate(tot_id = sprintf("%s_nfeat%s_normfeat%s_normtype%s_onlyshared%s_beta%s", model, nfeat, norm_feat, norm_type, only_shared, beta))
  
  df_auc_type[[i]] <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  nfeat = info$nfeat, 
                  norm_feat = info$norm_feat, 
                  norm_type = info$norm_type,
                  only_shared = info$only_shared, 
                  beta = info$beta) %>%
    dplyr::mutate(tot_id = sprintf("%s_nfeat%s_normfeat%s_normtype%s_onlyshared%s_beta%s", model, nfeat, norm_feat, norm_type, only_shared, beta))
  
  df_auc_type_p[[i]] <- read.csv(sprintf("%s/plots/private/AUC_CLs_pertissue_Same_type.csv", folder_model)) %>%
    dplyr::mutate(model = info$model, 
                  samples = info$samples, 
                  nfeat = info$nfeat, 
                  norm_feat = info$norm_feat, 
                  norm_type = info$norm_type,
                  only_shared = info$only_shared, 
                  beta = info$beta) %>%
    dplyr::mutate(tot_id = sprintf("%s_nfeat%s_normfeat%s_normtype%s_onlyshared%s_beta%s", model, nfeat, norm_feat, norm_type, only_shared, beta))
  
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
              aes(x = factor(tot_id, levels = ord_tot_id), y = AUC, fill = nfeat)) + 
  geom_boxplot() + 
  facet_wrap(.~only_shared, ncol = 1, scales = "free_y") +
  theme_bw() +
  # scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") + 
  ylab("AUC (Same tissue)")
pl1

# ord_tot_id <- df_auc_site_p %>%
#   dplyr::group_by(tot_id) %>%
#   dplyr::summarize(median_auc = median(AUC)) %>%
#   dplyr::arrange(desc(median_auc)) %>%
#   dplyr::pull(tot_id)
pl1_p <- ggplot(subset(df_auc_site_p), 
                aes(x = factor(tot_id, levels = ord_tot_id), y = AUC, fill = nfeat)) + 
  geom_boxplot() + 
  facet_wrap(.~only_shared, ncol = 1, scales = "free_y") +
  theme_bw() +
  # scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") + 
  ylab("AUC (Same tissue)")
pl1_p

# ord_tot_id <- df_auc_study %>%
#   dplyr::mutate(auc_min05 = AUC) %>%
#   dplyr::group_by(tot_id) %>%
#   dplyr::summarize(median_auc = median(auc_min05)) %>%
#   dplyr::arrange((median_auc)) %>%
#   dplyr::pull(tot_id)
pl2 <- ggplot(subset(df_auc_study), aes(x = factor(tot_id, levels = ord_tot_id), 
                                y = AUC, fill = nfeat)) +
  geom_boxplot() + 
  facet_wrap(.~only_shared, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") + 
  ylab("AUC (Same study: DepMap, TCGA)")
pl2  

# ord_tot_id <- df_auc_study_p %>%
#   dplyr::mutate(auc_min05 = abs(AUC - 0.5)) %>%
#   dplyr::group_by(tot_id) %>%
#   dplyr::summarize(median_auc = median(auc_min05)) %>%
#   dplyr::arrange((median_auc)) %>%
#   dplyr::pull(tot_id)
pl2_p <- ggplot(subset(df_auc_study_p), aes(x = factor(tot_id, levels = ord_tot_id), 
                                        y = AUC, fill = nfeat)) +
  geom_boxplot() + 
  facet_wrap(.~only_shared, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "bottom") + 
  coord_flip() + 
  xlab("") + 
  ylab("AUC (Same study: DepMap, TCGA)")
pl2_p  


### get medians 
# create summary medians 
summary_auc <- function(df_auc_site, df_auc_study, samples_class = "tcgaonly") {
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
    dplyr::summarize(median_auc_min05 = median(auc_min05), 
                     median_auc = median(AUC)) %>%
    dplyr::arrange(median_auc, decreasing = TRUE) %>%
    dplyr::mutate(tot_id = factor(tot_id, levels = .$tot_id)) %>%
    dplyr::select(tot_id, median_auc_min05, median_auc) %>%
    dplyr::mutate(type = "Same study") %>%
    dplyr::rename(median_auc_min05_study = median_auc_min05, 
                  median_auc_study = median_auc) %>%
    dplyr::mutate(order = 1:length(tot_id))
  
  df_auc_medians <- full_join(df_auc_site_medians, df_auc_study_medians, by = "tot_id", 
                              suffix = c("_site", "_study")) %>%
    mutate(order_sum = order_site + order_study) %>%
    arrange(order_sum)

  return(df_auc_medians)
}
df_auc_medians_tcgaonly <- summary_auc(df_auc_site, df_auc_study, samples_class = "tcgaonly")

#View(df_auc_medians_tcgaonly)
write.table(df_auc_medians_tcgaonly, 
            file = "experiments/auc_medians_tcgaonly_ex3.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

df_auc_p_medians_tcgaonly <- summary_auc(df_auc_site_p, df_auc_study_p, samples_class = "tcgaonly")
write.table(df_auc_medians_all, 
            file = "experiments/auc_medians_private_tcgaonly_ex3.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

