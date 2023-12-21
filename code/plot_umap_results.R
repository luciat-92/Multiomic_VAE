library(umap)
library(ggplot2)
library(dplyr)
library(tidyverse)

# load results
setwd("/Volumes/iorio/lucia/Multiomic_VAE/output/")

folder_model <- "model_save/code_adv_norm/pretrain_num_epochs_500_train_num_epochs_1000_dop_0.0_samples__ngene_all_norm_feat_flag_False/"
umap_df <- read_csv(sprintf("%s/umap.csv", folder_model))
depmap_sample_df <- read_csv("/Volumes/iorio/lucia/datasets/DEPMAP_PORTAL/version_23Q2/Model.csv")
enc_df <- read_csv(sprintf("%s/encoded_features.csv", folder_model))

# refactor type, first xena then depmap
umap_df <- umap_df %>% 
  mutate(type = factor(type, levels = c("xena", "depmap"))) %>%
  arrange(type) %>%
  mutate(study = ifelse(study == "CL_depmap", "Cell Line - DepMap", paste0("Tissue - ", study))) %>%
  filter(!is.na(site))

pl <- ggplot(subset(umap_df), aes(x = umap_1, y = umap_2, color = study, size = type)) + 
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
pl
ggsave(sprintf("%s/umap_study_ggplot2.png", folder_model), pl, width = 6, height = 6)

umap_df <- umap_df %>% 
  mutate(
    site = case_when(site %in% c("White blood cell", "Myeloid", "Lymphoid") ~ "Immune cells", 
                     site %in% c("Bladder", "Bladder/Urinary Tract") ~ "Bladder/Urinary Tract",
                     site %in% c("Rectum", "Colon", "Bowel") ~ "Colon/Rectum",
                     TRUE ~ site)
  )

# select random colors
n_colors <- length(unique(umap_df$site))
all_c <- colors()
set.seed(12344)
colors <- sample(all_c[!grepl("grey", all_c) & !grepl("gray", all_c)], n_colors)

pl <- ggplot(subset(umap_df, umap_1 < -3 & umap_1 > -7 & umap_2 < 2 & umap_2 > -2), 
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
pl
ggsave(sprintf("%s/umap_site_ggplot2.png", folder_model), pl, width = 8, height = 6.2)

# compute correlation among encoded features
cor_enc_df <- cor(t(enc_df[, -1]))
colnames(cor_enc_df) <- enc_df$sample_id
rownames(cor_enc_df) <- enc_df$sample_id

sample_df <- umap_df %>%
  dplyr::select(-umap_1, -umap_2) %>%
  dplyr::rename(sample_type = !!("_sample_type")) %>%
  dplyr::mutate(sample_type = case_when(
    sample_type == "Solid Tissue Normal" ~ "Normal Tissue",
    sample_type %in% c("Primary Solid Tumor", "Additional - New Primary") ~ "Primary Tumor",
    sample_type == "Recurrent Solid Tumor" ~ "Recurrent Tumor",
    sample_type == "Additional Metastatic" ~ "Metastatic",
    sample_type %in% c("Primary Blood Derived Cancer - Peripheral Blood", "Primary Blood Derived Cancer - Bone Marrow") ~ "Primary Blood Derived Cancer",
    sample_type %in% c("Recurrent Blood Derived Cancer - Bone Marrow", "Recurrent Blood Derived Cancer - Peripheral Blood") ~ "Recurrent Blood Derived Cancer",
    sample_type %in% c("Post treatment Blood Cancer - Bone Marrow", "Post treatment Blood Cancer - Peripheral Blood", "Post treatment Blood Cancer - Blood") ~ "Post treatment Blood Cancer",
    TRUE ~ sample_type
  ))


sample_df <- sample_df[match(colnames(cor_enc_df), sample_df$sample_id), ]
id_depmap <- sample_df$type == "depmap"
id_all_CL <- which(id_depmap)
df <- list()
corr_threshold <- 0.9
for(i in 1:length(id_all_CL)){
  
  idx <- id_all_CL[i]
  id_CL <- sample_df$sample_id[idx]
  curr_corr <- cor_enc_df[id_CL, ]
  
  tmp_ann <- sample_df[curr_corr >= corr_threshold & !id_depmap, ]
  max_study_perc <-  best_study <- max_sample_type_perc <- best_sample_type <- max_site_perc <- best_site <- NA
  if(nrow(tmp_ann) > 0){
    max_study_perc <- max(table(tmp_ann$study) / sum(table(tmp_ann$study)))
    best_study <- paste0(names(table(tmp_ann$study)[table(tmp_ann$study) == max(table(tmp_ann$study))]), collapse = ",")
    max_sample_type_perc <- max(table(tmp_ann$sample_type) / sum(table(tmp_ann$sample_type)))
    best_sample_type <- paste0(names(table(tmp_ann$sample_type)[table(tmp_ann$sample_type) == max(table(tmp_ann$sample_type))]), collapse = ",")
    max_site_perc <- max(table(tmp_ann$site) / sum(table(tmp_ann$site)))
    best_site <- paste0(names(table(tmp_ann$site)[table(tmp_ann$site) == max(table(tmp_ann$site))]), collapse = ",")
  }
 
  df[[i]] <- data.frame(sample_df[idx, ], 
                        n_threshold = sum(curr_corr >= corr_threshold & !id_depmap),
                        best_perc_study = max_study_perc, 
                        best_study = best_study,
                        best_perc_sample_type = max_sample_type_perc,
                        best_sample_type = best_sample_type,
                        best_perc_site = max_site_perc,
                        best_site = best_site)
  
}

df <- do.call(rbind, df)
# use color scale of reds
# library(RColorBrewer)
png(sprintf("%s/heatmap_count_CL_site0.9corr.png", folder_model), width = 8, height = 10, units = "in", res = 300)
pheatmap::pheatmap(t(table(df$site, df$best_site)), 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   color = c("grey80", colorRampPalette(brewer.pal(9, "Reds"))(100)))
dev.off()

png(sprintf("%s/heatmap_count_CL_study0.9corr.png", folder_model), width = 5, height = 3, units = "in", res = 300)
pheatmap::pheatmap(t(table(df$site, df$best_study)), 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   color = c("grey80", colorRampPalette(brewer.pal(9, "Reds"))(100)))
# save
dev.off()

png(sprintf("%s/heatmap_count_CL_sample_type0.9corr.png", folder_model), width = 5, height = 3, units = "in", res = 300)
pheatmap::pheatmap(t(table(df$site, df$best_sample_type)), 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   fontsize = 7, 
                   cellwidth = 6, 
                   cellheight = 6, 
                   color = c("grey80", colorRampPalette(brewer.pal(9, "Reds"))(100)))
dev.off()



# remove comparisons with other cell lines from depmap

#id_normal <- sample_df$sample_type %in% c("Normal Tissue", "Cell Line")
#id_same_site <- sample_df$site %in% sample_df$site[idx]
#curr_corr_samesite_can <- curr_corr[!id_depmap & id_same_site & !id_normal]
#curr_corr_diffsite_can <- curr_corr[!id_depmap & !id_same_site & !id_normal]
#curr_corr_samesite_normal <- curr_corr[!id_depmap & id_same_site & id_normal]
#curr_corr_diffsite_normal <- curr_corr[!id_depmap & !id_same_site & id_normal]









