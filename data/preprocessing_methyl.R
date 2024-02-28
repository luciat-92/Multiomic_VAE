# preprocess methylation data
library(tidyverse)

#####
# input
file_CL_meth <- "/group/iorio/Alessandro.D/EpiClock/data/CLs_methylation_data.csv"
file_CL_meta <- "/group/iorio/Alessandro.D/EpiClock/data/annotations/IorioCell2016-MethylAccesionCellLines.txt"
file_CL_depmap <- "/group/iorio/lucia/datasets/DEPMAP_PORTAL/version_23Q2/Model.csv"
fold_SAMPLE_meth <- "/group/iorio/Irene/epiclock/data/"
fold_output <- "/group/iorio/lucia/Multiomic_VAE/data/preprocessed/"
#####

# load TCGA
# get all folders in fold_SAMPLE_meth
fold_percancer_meth <- list.dirs(path = fold_SAMPLE_meth, full.names = TRUE)
fold_percancer_meth  <- fold_percancer_meth[grepl("TCGA", fold_percancer_meth)]
percancer_meth <- list()
# sample_meta <- list()
for (i in 1:length(fold_percancer_meth)) {
  print(i)
  try(percancer_meth[[i]] <- get(load(sprintf("%s/data_met_preprocessed.RData", fold_percancer_meth[i]))))
  # sample_meta[[i]] <- get(load(sprintf("%s/data_samples.RData", fold_percancer_meth[i])))
}
# remove null one
id_null <- sapply(percancer_meth, is.null)
percancer_meth <- percancer_meth[!id_null]
print("loaded TCGA")

# load CL
CL_meth <- read_csv(file_CL_meth)
cpgs_names <- CL_meth[,1]
CL_meth <- as.matrix(CL_meth[, -1])
rownames(CL_meth) <- pull(cpgs_names)

CL_meta <- read_tsv(file_CL_meta) %>%
  dplyr::mutate(StrippedCellLineName = gsub("[^[:alnum:]]", "", Title.2)) %>%
  dplyr::distinct(StrippedCellLineName, .keep_all = TRUE)
CL_depmap_ann <- read_csv(file_CL_depmap)
CL_meta_tot <- inner_join(CL_depmap_ann, CL_meta, by = "StrippedCellLineName")
# filter samples 
CL_meth_filt <- CL_meth[, CL_meta_tot$CAccession]
colnames(CL_meth_filt) <- CL_meta_tot$ModelID
print("loaded CLs")

# get common list of variables
common_cpg <- Reduce(intersect, lapply(percancer_meth, rownames))
# merge with cpgs in cell line
common_cpg <- intersect(common_cpg, rownames(CL_meth_filt))
percancer_meth_common <- lapply(percancer_meth, function(x) x[common_cpg,])
percancer_meth_tot <- do.call(cbind, percancer_meth_common)
meta_tcga <- data.frame(complete_id = colnames(percancer_meth_tot))
meta_tcga$sample_id <- unname(sapply(meta_tcga$complete_id, function(x) substr(x, 1, 15)))
colnames(percancer_meth_tot) <- meta_tcga$sample_id

CL_meth_common <- CL_meth_filt[common_cpg,]
rm(percancer_meth)
rm(percancer_meth_common)

all_meth_tot <- t(cbind(CL_meth_common, percancer_meth_tot))
samples_list <- list(depmap = colnames(CL_meth_common), tcga = colnames(percancer_meth_tot))
rm(percancer_meth_tot)
rm(CL_meth_common)
print("got common list of cpgs")

# get the most variable meths
get_most_variable_meth <- function(df, samples = NULL, n_meth = 1000) {
  # Filter based on samples
  if (!is.null(samples)) {
    df <- df[samples, , drop = FALSE]
  }
  # Calculate standard deviation for each gene
  std <- apply(df, 2, sd)
  # Sort standard deviations in descending order
  std_top <- sort(std, decreasing = TRUE)
  # Get the top n meths
  std_top <- std_top[1:n_meth]
  # Get the gene names
  meth_names <- names(std_top)
  return(meth_names)
}

filter_methylation <- function(df, samples = NULL, n_meth = 1000) {
  # If samples is NULL, use all samples
  if (is.null(samples)) {
    meths <- get_most_variable_meth(df, n_meth)
  } else {
    meths <- c()
    for (key in names(samples)) {
      samples_id <- samples[[key]]
      print(length(samples_id))
      meths_tmp <- get_most_variable_meth(df = df[samples_id, ], n_meth = n_meth)
      meths <- c(meths, meths_tmp)
    }
    meths <- unique(meths)
  }
  
  # Log number of common meths
  cat(paste("Number of common methylation sites:", length(meths), "from top", n_meth, "most variable sites\n"))
  # Filter based on meths
  df <- df[, meths, drop = FALSE]
  return(df)
}

gc(); gc(); gc()

all_samples_names <- rownames(all_meth_tot)
var1000_meths_df <- filter_methylation(
  df = all_meth_tot, samples = samples_list, 
  n_meth = 1000)

# save
save_var1000_meths_df <- as.data.frame(var1000_meths_df) %>%
  dplyr::mutate(sample_id = all_samples_names, .before = colnames(var1000_meths_df)[1])
write_csv(save_var1000_meths_df, file = sprintf("%smethylation_var1000_tcgaonly.csv", fold_output))

var5000_meths_df <- filter_methylation(
  df = all_meth_tot, samples = samples_list, 
  n_meth = 5000)
save_var5000_meths_df <- as.data.frame(var5000_meths_df) %>%
  dplyr::mutate(sample_id = all_samples_names, .before = colnames(var5000_meths_df)[1])
# save
write_csv(save_var5000_meths_df, file = sprintf("%smethylation_var5000_tcgaonly.csv", fold_output))

# save
save_all <- as.data.frame(all_meth_tot) %>%
  dplyr::mutate(sample_id = all_samples_names, .before = colnames(all_meth_tot)[1])
write_csv(save_all, file = sprintf("%smethylation.csv", fold_output))


