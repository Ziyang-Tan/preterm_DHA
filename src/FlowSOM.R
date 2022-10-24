# FlowSOM clustering

# 0. Load packages ----------

library(FlowSOM)
library(flowCore)
library(tidyverse)
library(magrittr)
library(paletteer)

# 0. Load data ----------

ff = read.FCS('data/cytof_outlier99_arcsinh5_combat_data.fcs', 
              transformation = FALSE, truncate_max_range = FALSE)
label = readr::read_csv('data/cytof_outlier99_arcsinh5_combat_label.csv.gz', 
                        col_names = F) 
#sample_info = readxl::read_excel('data/210825_database_for_petter_sent_210826.xlsx') 

# 1. FlowSOM 1: fish out the Neutrophils ----------

common_marker = colnames(ff@exprs) %>% unname()

# Run FlowSOM
fSOM_30 = FlowSOM(ff,
                  compensate = F,
                  transform = F,
                  scale = F,
                  colsToUse = common_marker,
                  nClus = 10,
                  xdim = 5, ydim = 6,
                  seed = 824)

# Median marker expression of FlowSOM cell clusters
MFI_30 = GetClusterMFIs(fSOM_30, colsUsed = T)
clusters_30 = GetClusters(fSOM_30)
freqClusters_30 = data.frame(clusters = clusters_30) %>%
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_30 = cbind(MFI_30, freqClusters_30)

# Heatmap of FlowSOM clusters
q99 = quantile(as.matrix(MFI_30), 0.99)
pheatmap::pheatmap(as.matrix(MFI_30),
                   scale = 'none',
                   labels_row = paste(rownames(MFI_30),' (', round(freqClusters_30$percentage,1), '%', ')',sep = ''),
                   display_numbers = TRUE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "30 clusters: Median marker expression per cluster")

# 2. Separate Neutrophil clusters ----------

# Identify Neutrophil clusters according to heatmap
nrow_neutrophils = which(clusters_30 %in% c(22,23,24,26,27,28,29,30))

# Separate Neutrophil data
label_neutrophils = label[nrow_neutrophils,]
clusters_neutrophils = clusters_30[nrow_neutrophils]
label_neutrophils$subtype1 = 'Neutrophils'
label_neutrophils$clusters = clusters_neutrophils
data_neutrophils = ff@exprs[nrow_neutrophils,]
res_neutrophils = res_30 %>% dplyr::filter(clusters %in% c(22,23,24,26,27,28,29,30))
res_neutrophils$subtype1 = 'Neutrophils'

# Remaining cells after removing Neutrophils
data_remaining = ff@exprs[-nrow_neutrophils,]
label_remaining = label[-nrow_neutrophils,]

# 3. FlowSOM 100 clusters on remaining cells ----------

ff_remaining = ff
ff_remaining@exprs = data_remaining

# Run FlowSOM
fSOM_100 = FlowSOM(ff_remaining,
                   compensate = F,
                   transform = F,
                   scale = F,
                   colsToUse = common_marker,
                   nClus = 10,
                   xdim = 10, ydim = 10,
                   seed = 824)

# FlowSOM results summary
# FlowSOMmary(fSOM_100, plotFile = 'flowsom_summary_100clusters.pdf')

# Median marker expression of FlowSOM cell clusters
MFI_100 = GetClusterMFIs(fSOM_100, colsUsed = T)
MFI_100 <- scale(MFI_100) # zscore transfer
clusters_100 = GetClusters(fSOM_100)
freqClusters_100 = data.frame(clusters = clusters_100) %>%
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()

# Heatmap of FlowSOM clusters
q99 = quantile(as.matrix(MFI_100), 0.99)
pheatmap::pheatmap(as.matrix(MFI_100),
                   scale = 'none',
                   labels_row = paste(rownames(MFI_100),' (', round(freqClusters_30$percentage,1), '%', ')',sep = ''),
                   display_numbers = TRUE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "100 clusters: Median marker expression per cluster")

# Manual annotation

# load data (skip if already)

#load(file = 'data/flowsom_results.RData')
#all_bind = mmR::mm.fastread('data/flowSOM_all_data.csv.gz')

annotation = read.csv2('data/flowSOM100_annotation.csv', check.names = F, stringsAsFactors = F)

res_100 <- cbind(MFI_100, freqClusters_100) %>%
  left_join(annotation, by='clusters')

label_remaining$clusters = clusters_100
label_remaining %<>% left_join(annotation, by='clusters')

# 4. Bind remaining cells with neutrophils ----------

data_bind = rbind(data_neutrophils, data_remaining)
label_neutrophils$clusters = label_neutrophils$clusters+100 # avoid confusion with results of flowSOM100
res_neutrophils$clusters = res_neutrophils$clusters+100 # avoid confusion with results of flowSOM100
label_bind = rbind(label_neutrophils, label_remaining)
label_bind %<>% mutate(cell_cluster = paste(subtype1, '_', clusters, sep = '')) %>% 
  dplyr::rename(ID_unique = X5)
all_bind = cbind(data_bind, label_bind) 
res_bind = rbind(res_neutrophils, res_100)
res_bind %<>% mutate(frequency = (n/sum(res_bind$n))*100)
res_bind %<>% mutate(cell_cluster = paste(subtype1, '_', clusters, sep = ''))

# Heatmap of all cell clusters except T cells (will do second flowSOM for T cells)
d <- res_bind[res_bind$subtype1 != 'T cell',]
q99 = quantile(as.matrix(d[,1:length(common_marker)]), 0.99)
pheatmap::pheatmap(as.matrix(d[,1:length(common_marker)]),
                   color = paletteer::paletteer_d("rcartocolor::PurpOr", n = 100, type = 'continuous'), # New color
                   border_color = 'white',
                   scale = 'none',
                   labels_row = paste(d$cell_cluster,' (', round(d$percentage,1), '%', ')',sep = ''),
                   display_numbers = FALSE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "Median marker expression per cluster",
                   filename = 'figures/flowSOM_annotated_remove_Tcells.pdf',
                   height = 15, width = 18)

# Calculate relative frequency of each celltype in each sample
# freq_lineage = label_bind %>%
#   group_by(ID_unique) %>%
#   count(.data$lineage) %>%
#   mutate(percentage = .data$n / sum(.data$n) * 100) %>%
#   dplyr::select(-n) %>%
#   spread(lineage, percentage)
# freq_subtype = label_bind %>%
#   group_by(ID_unique) %>%
#   count(.data$subtype) %>%
#   mutate(percentage = .data$n / sum(.data$n) * 100) %>%
#   dplyr::select(-n) %>%
#   spread(subtype, percentage)
# freq_cellcluster = label_bind %>%
#   group_by(ID_unique) %>%
#   count(.data$cell_cluster) %>%
#   mutate(percentage = .data$n / sum(.data$n) * 100) %>%
#   dplyr::select(-n) %>%
#   spread(cell_cluster, percentage)

# 5. Save results ----------

#save(all_bind, file = 'data/flowsom_results.RData')
all_bind <- all_bind %>% select(common_marker, ID_unique, subtype1)
mmR::mm.fastwrite(all_bind, path = 'data/flowSOM_all_data.csv', compress = T, overwrite = T)
#write.csv(freq_lineage, file = 'data/flowsom_freq_lineage.csv', row.names = F)
#write.csv(freq_subtype, file = 'data/flowsom_freq_subtype.csv', row.names = F)
#write.csv(freq_cellcluster, file = 'data/flowsom_freq_cellcluster.csv', row.names = F)

# Save median marker expression data for Network plot
write.table(res_bind, file = 'data/flowsom_clustered.txt', 
            sep="\t", row.names=F, col.names = T, quote = F)


