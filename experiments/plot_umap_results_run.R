library(umap)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(argparse)

# load input with args
parser <- ArgumentParser(description = "")
parser$add_argument("--folder_model", type = "character", help = "Input and Output fold")
parser$add_argument("--depmap_meta_file", type = "character", help = "file depmap info")
args <- parser$parse_args()
print(args)

folder_model <- args$folder_model
depmap_meta_file <- args$depmap_meta_file

###########
folder_model <- "experiments/experiment_2/vae_gan/samples_tcgaonly_ngene_var1000_norm_feat_flag_False_only_shared_True/"
depmap_meta_file <- "/Volumes/iorio/lucia/datasets/DEPMAP_PORTAL/version_23Q2/Model.csv"
###########

print("Start")

umap_df <- readr::read_csv(sprintf("%s/plots/umap.csv", folder_model))
depmap_sample_df <- readr::read_csv(depmap_meta_file)
enc_df <- readr::read_csv(sprintf("%s/encoded_features.csv", folder_model))

# refactor type, first xena then depmap
umap_df <- umap_df %>% 
  mutate(type = factor(type, levels = c("xena", "depmap"))) %>%
  arrange(type) %>%
  mutate(study = ifelse(study == "CL_depmap", "Cell Line - DepMap", paste0("Tissue - ", study))) %>%
  filter(!is.na(site)) %>%
  rename(sample_type = !!("_sample_type"))

pl1 <- ggplot(subset(umap_df), 
             aes(x = umap_1, y = umap_2, color = study, size = type)) + 
  geom_point(alpha = 0.5) +
  scale_size_manual(values = c(0.2, 1)) +
  scale_color_manual(values = c("#e8702a", "#6bd2db", "#9ed670", "#0c457d")) +
  theme_classic() + 
  theme(legend.title = element_blank(),
        legend.position = c(0.99, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 10)) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(size = 'none')
ggsave(sprintf("%s/plots/umap_study_ggplot2.png", folder_model), pl1, width = 6, height = 6)

# select random colors
n_colors <- length(unique(umap_df$site))
all_c <- colors()
set.seed(12344)
colors <- sample(all_c[!grepl("grey", all_c) & !grepl("gray", all_c)], n_colors)

pl2 <- ggplot(subset(umap_df), 
             aes(x = umap_1, y = umap_2, size = type)) + 
  geom_point(aes(fill = site, color = type), shape = 21, alpha = 0.8) + 
  scale_size_manual(values = c(0.7, 1)) +
  scale_color_manual(values = c("transparent", "black")) +
  # scale_fill_manual(values = colors) +
  theme_classic() + 
  theme(legend.title = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10), 
        # reduce legend spacing
        legend.key.height = unit(0.1, "cm"),  
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  xlab("UMAP 1") +
  ylab("UMAP 2") + 
  guides(fill = guide_legend(ncol = 1), size = 'none', color = 'none')
ggsave(sprintf("%s/plots/umap_site_ggplot2.png", folder_model), pl2, width = 8, height = 6.2)
# ggarrange(plotlist = list(pl1, pl2), nrow = 1)

# # compute correlation among encoded features
# cor_enc_df <- cor(t(enc_df[, -1]))
# colnames(cor_enc_df) <- enc_df$sample_id
# rownames(cor_enc_df) <- enc_df$sample_id
# 
# sample_df <- umap_df %>%
#   dplyr::select(-umap_1, -umap_2) %>%
#   dplyr::rename(sample_type = !!("_sample_type")) %>%
#   dplyr::mutate(sample_type = case_when(
#     sample_type == "Solid Tissue Normal" ~ "Normal Tissue",
#     sample_type %in% c("Primary Solid Tumor", "Additional - New Primary") ~ "Primary Tumor",
#     sample_type == "Recurrent Solid Tumor" ~ "Recurrent Tumor",
#     sample_type == "Additional Metastatic" ~ "Metastatic",
#     sample_type %in% c("Primary Blood Derived Cancer - Peripheral Blood", "Primary Blood Derived Cancer - Bone Marrow") ~ "Primary Blood Derived Cancer",
#     sample_type %in% c("Recurrent Blood Derived Cancer - Bone Marrow", "Recurrent Blood Derived Cancer - Peripheral Blood") ~ "Recurrent Blood Derived Cancer",
#     sample_type %in% c("Post treatment Blood Cancer - Bone Marrow", "Post treatment Blood Cancer - Peripheral Blood", "Post treatment Blood Cancer - Blood") ~ "Post treatment Blood Cancer",
#     TRUE ~ sample_type
#   ))
# 
# 
# sample_df <- sample_df[match(colnames(cor_enc_df), sample_df$sample_id), ]
# id_depmap <- sample_df$type == "depmap"
# id_all_CL <- which(id_depmap)
# df <- list()
# corr_threshold <- 0.9
# for(i in 1:length(id_all_CL)){
#   
#   idx <- id_all_CL[i]
#   id_CL <- sample_df$sample_id[idx]
#   curr_corr <- cor_enc_df[id_CL, ]
#   
#   tmp_ann <- sample_df[curr_corr >= corr_threshold & !id_depmap, ]
#   max_study_perc <-  best_study <- max_sample_type_perc <- best_sample_type <- max_site_perc <- best_site <- NA
#   if(nrow(tmp_ann) > 0){
#     max_study_perc <- max(table(tmp_ann$study) / sum(table(tmp_ann$study)))
#     best_study <- paste0(names(table(tmp_ann$study)[table(tmp_ann$study) == max(table(tmp_ann$study))]), collapse = ",")
#     max_sample_type_perc <- max(table(tmp_ann$sample_type) / sum(table(tmp_ann$sample_type)))
#     best_sample_type <- paste0(names(table(tmp_ann$sample_type)[table(tmp_ann$sample_type) == max(table(tmp_ann$sample_type))]), collapse = ",")
#     max_site_perc <- max(table(tmp_ann$site) / sum(table(tmp_ann$site)))
#     best_site <- paste0(names(table(tmp_ann$site)[table(tmp_ann$site) == max(table(tmp_ann$site))]), collapse = ",")
#   }
#  
#   df[[i]] <- data.frame(sample_df[idx, ], 
#                         n_threshold = sum(curr_corr >= corr_threshold & !id_depmap),
#                         best_perc_study = max_study_perc, 
#                         best_study = best_study,
#                         best_perc_sample_type = max_sample_type_perc,
#                         best_sample_type = best_sample_type,
#                         best_perc_site = max_site_perc,
#                         best_site = best_site)
#   
# }
# 
# df <- do.call(rbind, df)
# # use color scale of reds
# # 
# png(sprintf("%s/heatmap_count_CL_site0.9corr.png", folder_model), width = 8, height = 10, units = "in", res = 300)
# pheatmap::pheatmap(t(table(df$site, df$best_site)), 
#                    cluster_rows = FALSE, 
#                    cluster_cols = FALSE, 
#                    color = c("grey80", colorRampPalette(brewer.pal(9, "Reds"))(100)))
# dev.off()
# 
# png(sprintf("%s/heatmap_count_CL_study0.9corr.png", folder_model), width = 5, height = 3, units = "in", res = 300)
# pheatmap::pheatmap(t(table(df$site, df$best_study)), 
#                    cluster_rows = FALSE, 
#                    cluster_cols = FALSE, 
#                    color = c("grey80", colorRampPalette(brewer.pal(9, "Reds"))(100)))
# # save
# dev.off()
# 
# png(sprintf("%s/heatmap_count_CL_sample_type0.9corr.png", folder_model), width = 5, height = 3, units = "in", res = 300)
# pheatmap::pheatmap(t(table(df$site, df$best_sample_type)), 
#                    cluster_rows = FALSE, 
#                    cluster_cols = FALSE, 
#                    fontsize = 7, 
#                    cellwidth = 6, 
#                    cellheight = 6, 
#                    color = c("grey80", colorRampPalette(brewer.pal(9, "Reds"))(100)))
# dev.off()

# compute euclidean distance
dist_enc <- as.matrix(dist(enc_df[, -1]))
colnames(dist_enc) <- enc_df$sample_id
rownames(dist_enc) <- enc_df$sample_id

sample_df <- umap_df %>%
  dplyr::select(-umap_1, -umap_2)

sample_df <- sample_df %>%
  dplyr::mutate(sample_type_macro = case_when(
    sample_type %in% c("Normal Tissue", "Cell Line") ~ "Non-cancerous",
    sample_type %in% c("Primary Tumor", "Recurrent Tumor", "Metastatic") ~ "Solid Tumor",
    sample_type %in% c("Primary Blood Derived Cancer", "Recurrent Blood Derived Cancer", "Post treatment Blood Cancer") ~ "Blood Tumor",
    TRUE ~ sample_type
  ))
id_depmap <- sample_df$type == "depmap"
id_all_CL <- which(id_depmap)

sample_df$sample_type_macro[id_depmap] <- "Solid Tumor"
sample_df$sample_type_macro[sample_df$sample_id %in% depmap_sample_df$ModelID[depmap_sample_df$OncotreePrimaryDisease == "Non-Cancerous"]] <- "Non-cancerous"
sample_df$sample_type_macro[id_depmap & sample_df$site == "Immune cells"] <- "Blood Tumor"

common_s <- intersect(sample_df$sample_id, enc_df$sample_id)
sample_df <- sample_df[match(common_s, sample_df$sample_id), ]
dist_enc <- dist_enc[common_s, common_s]

df_auc_site <- list()
df_auc_study <- list()
df_auc_type <- list()

for(i in 1:length(id_all_CL)){
  
  # print(i)
  idx <- id_all_CL[i]
  id_CL <- sample_df$sample_id[idx]
  curr_dist <- dist_enc[id_CL, ]
  
  # categorize based on study
  curr_s <- sample_df$study[idx]
  bin_study <- as.numeric(sample_df$study == curr_s)
  
  # compute AURC with pROC
  roc_obj <- pROC::roc(response = bin_study, predictor = curr_dist,
                       # arguments for ci
                       ci=TRUE, boot.n=1000, ci.alpha=0.95, stratified=TRUE, 
                       direction = ">", quiet = TRUE)
  ci <- roc_obj$ci
  # get auc and CI
  df_auc_study[[i]] <- data.frame(CI_low = as.numeric(ci)[1], 
                       AUC = as.numeric(ci)[2], 
                       CI_up = as.numeric(ci)[3], 
                       comp = "study", 
                       sample_id = id_CL, 
                       sample_info = sample_df$site[idx])
  
  # categorize based on site
  curr_t <- sample_df$site[idx]
  bin_site <- as.numeric(sample_df$site == curr_t)
  
  # compute AURC with pROC
  if(length(unique(bin_site[!id_depmap])) > 1){
    roc_obj <- pROC::roc(response = bin_site[!id_depmap], predictor = curr_dist[!id_depmap],
                         # arguments for ci
                         ci=TRUE, boot.n=1000, ci.alpha=0.95, stratified=TRUE, 
                         direction = ">", quiet = TRUE)
    ci <- roc_obj$ci
    # get auc and CI
    df_auc_site[[i]] <- data.frame(CI_low = as.numeric(ci)[1], 
                                   AUC = as.numeric(ci)[2], 
                                   CI_up = as.numeric(ci)[3], 
                                   comp = "site", 
                                   sample_id = id_CL, 
                                   sample_info = curr_t)
  }else{
    df_auc_site[[i]] <- data.frame(CI_low = NA, 
                                   AUC = NA, 
                                   CI_up = NA, 
                                   comp = "site", 
                                   sample_id = id_CL, 
                                   sample_info = curr_t)
  }
 
  # categorize based on type
  curr_tm <- sample_df$sample_type_macro[idx]
  bin_type <- as.numeric(sample_df$sample_type_macro == curr_tm)
  
  # compute AURC with pROC
  if(length(unique(bin_type[!id_depmap])) > 1){
    roc_obj <- pROC::roc(response = bin_type[!id_depmap], predictor = curr_dist[!id_depmap],
                         # arguments for ci
                         ci=TRUE, boot.n=1000, ci.alpha=0.95, stratified=TRUE, 
                         direction = ">", quiet = TRUE)
    ci <- roc_obj$ci
    # get auc and CI
    df_auc_type[[i]] <- data.frame(CI_low = as.numeric(ci)[1], 
                                   AUC = as.numeric(ci)[2], 
                                   CI_up = as.numeric(ci)[3], 
                                   comp = "type", 
                                   sample_id = id_CL, 
                                   sample_info = curr_t, 
                                   sample_type_macro = curr_tm)
  }else{
    df_auc_type[[i]] <- data.frame(CI_low = NA, 
                                   AUC = NA, 
                                   CI_up = NA, 
                                   comp = "type", 
                                   sample_id = id_CL, 
                                   sample_info = curr_t, 
                                   sample_type_macro = curr_tm)
  }
}
df_auc_site <- do.call(rbind, df_auc_site)
df_auc_study <- do.call(rbind, df_auc_study)
df_auc_type <- do.call(rbind, df_auc_type)

plot_auc_pertissue <- function(df_auc, title_pl, folder_model, save_csv = TRUE){
  
  # drop NA (not the same site name between tissue and cell lines)
  df_auc <- df_auc[!is.na(df_auc$AUC), ]
  # order class based on the median AUC in each sample info
  df_auc$sample_info <- factor(df_auc$sample_info, 
                               levels = names(sort(tapply(df_auc$AUC, df_auc$sample_info, median))))
  pl <- ggplot(df_auc, aes(x = sample_info, y = AUC, fill = sample_info)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.7) + 
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none") + 
    # facet_wrap(.~comp, scales = "free_x") + 
    labs(x = "", y = "AUCROC") + 
    ggtitle(title_pl)
  
  print(pl)
  ggsave(sprintf("%s/boxplotAUC_CLs_pertissue_%s.pdf", folder_model, gsub(" ", "_", title_pl)), 
         pl, width = 6, height = 4.5)
  if(save_csv){
    write.csv(df_auc, sprintf("%s/AUC_CLs_pertissue_%s.csv", folder_model, gsub(" ", "_", title_pl)), 
              row.names = F)
  }
  
  if( "sample_type_macro" %in% colnames(df_auc)){
    
  df_auc$sample_type_macro <- factor(df_auc$sample_type_macro, 
                                 levels = names(sort(tapply(df_auc$AUC, df_auc$sample_type_macro, median))))
   pl2 <- ggplot(df_auc, aes(x = sample_type_macro, y = AUC, fill = sample_type_macro)) + 
     geom_boxplot(outlier.shape = NA) + 
     geom_jitter(width = 0.2, alpha = 0.5, size = 0.7) + 
     geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
     theme_bw() + 
     theme(axis.text.x = element_text(angle = 45, hjust = 1), 
           legend.position = "none") + 
     # facet_wrap(.~comp, scales = "free_x") + 
     labs(x = "", y = "AUCROC") + 
     ggtitle(title_pl)
   
   print(pl2)
   ggsave(sprintf("%s/boxplotAUC_CLs_pertype_%s.pdf", folder_model, gsub(" ", "_", title_pl)), 
          pl2, width =4.56, height = 4.5)
   if(save_csv){
     write.csv(df_auc, sprintf("%s/AUC_CLs_pertype_%s.csv", folder_model, gsub(" ", "_", title_pl)), 
               row.names = F)
   }
  }
  
}

plot_auc_pertissue(df_auc_site, "Same tissue", sprintf("%s/plots/", folder_model))
plot_auc_pertissue(df_auc_study, "Same study", sprintf("%s/plots/", folder_model))
plot_auc_pertissue(df_auc_type, "Same type", sprintf("%s/plots/", folder_model))

## load results
#folder_model_1 <- "experiments/experiment_1/ae_gan/samples__ngene_all_norm_feat_flag_True_only_shared_True/"
#folder_model_2 <- "experiments/experiment_2/vae_gan/samples__ngene_all_norm_feat_flag_True_only_shared_False/"
#
#df_auc_site_1 <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_tissue.csv", folder_model_1)) %>%
#  dplyr::mutate(model = "ae_nonorm_only_shared")
#
#df_auc_site_2 <- read.csv(sprintf("%s/plots/AUC_CLs_pertissue_Same_tissue.csv", folder_model_2)) %>%
#  mutate(model = "vae_norm_noonly_shared")
#
#df_auc_site_tot <- rbind(df_auc_site_1, df_auc_site_2)
#boxplot(df_auc_site_tot$AUC ~ df_auc_site_tot$model)

