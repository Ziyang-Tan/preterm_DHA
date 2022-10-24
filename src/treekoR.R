# treekoR analysis

# 0. Load packages ----------

library(tidyverse)
library(magrittr)
library(data.table)
library(treekoR)
library(ggsignif)

# 0. Load data ----------

load(file = 'data/flowsom_results.RData')
sample_info = sample_info = readxl::read_excel('data/210825_database_for_petter_sent_210826.xlsx') 

# 1. Define cohort and data ----------

# Comparision1: Neuroblastoma: Metastasis vs. Non-metastasis
cohort = sample_info %>% dplyr::filter(tumor_level2 == 'Neuroblastoma')
data = all_bind %>% dplyr::filter(study_id %in% cohort$study_id) 
comparison = 'metastasis'
comparison_level = c('No', 'Yes')
title = 'Neuroblastoma - Metastasis vs. Non-metastasis'

# 2. Input for treekoR ----------

common_marker = colnames(all_bind)[1:32]
exprs = data[,1:length(common_marker)] %>% as.matrix()
clusters = data$cell_cluster %>% as.factor()
samples = data$study_id %>% as.factor()
classes = data$metastasis %>% as.character() %>% factor(levels = c('No', 'Yes')) 

# 3. Manually fixed lineage tree ----------

getFixedClusterTree = function(exprs, clusters, samples, classes){
  # calculate median marker expression
  clust_med_dt = as.data.table(exprs)
  clust_med_dt[, cluster_id := clusters]
  res = clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
  res2 = res[,.SD, .SDcols = !c('cluster_id')]
  rownames(res2) = res[["cluster_id"]]
  res2[, (colnames(res2)) := lapply(.SD, scale), .SDcols=colnames(res2)]
  order = order(rownames(res2))
  res2 = res2[order,]
  rownames(res2) = res[["cluster_id"]][order]
  
  # label the position of cell clusters in the tree
  cell_level3 = unique(rownames(res2))
  cell_level2 = lapply(strsplit(cell_level3,'_'), function(x){x[1]}) %>% unlist()
  cell_order = data.frame(cell_level2 = unique(cell_level2), node = 1:length(unique(cell_level2)))
  cell_nodes = cell_order %>% right_join(data.frame(cell_level2 = cell_level2, cell_level3 = cell_level3))
  
  cutree_1 = rep(1,length(cell_level3))
  names(cutree_1) = 1:length(cell_level3)
  cutree_2 = cell_nodes$node
  names(cutree_2) = 1:length(cell_level3)
  cutree_3 = 1:length(cell_level3)
  names(cutree_3) = 1:length(cell_level3)
  cutree_list = list(cutree_1, cutree_2, cutree_3)
  hp_dend = list()
  hp_dend$cutree_list = cutree_list
  
  # generate phylogram
  hc_phylo = hopachToPhylo(hp_dend)
  hc_phylo$tip.label = rownames(res2)[as.numeric(hc_phylo$tip.label)]
  
  clust_tree=list(
    median_freq = res2,
    clust_tree = hc_phylo
  )
  
  return(clust_tree)
}

clust_tree = getFixedClusterTree(exprs, clusters, samples, classes)

# 4. Significance testing of cell subpopulations ----------

tested_tree = testTree(phylo=clust_tree$clust_tree,
                       clusters=clusters,
                       samples=samples,
                       classes=classes,
                       sig_test = 'wilcox',
                       pos_class_name=NULL)

# 5. Interactive heatmap ----------

plotInteractiveHeatmap(tested_tree,
                       clust_med_df = clust_tree$median_freq,
                       clusters=clusters)

# 6. Feature extraction ----------

# Significance testing results
res_df = getTreeResults(tested_tree)
res_df$clusters = res_df$clusters %>% as.character()
sig_df = res_df %>% dplyr::filter(pval_total < 0.05 | pval_parent < 0.05)

# Percentage of each cell cluster
prop_df = getCellProp(phylo=clust_tree$clust_tree,
                      clusters=clusters,
                      samples=samples,
                      classes=classes)

# Mean marker expression in each cell cluster
means_df = getCellGMeans(clust_tree$clust_tree,
                         exprs=exprs,
                         clusters=clusters,
                         samples=samples,
                         classes=classes)

# Log2fc of mean percentage of each cell cluster between two comparison groups
log2fc_df = prop_df %>% gather(variable, value, -c(class, sample_id)) 
log2fc_df$class = as.integer(log2fc_df$class)
log2fc_df[is.na(log2fc_df)] = 0
log2fc_df$perc = lapply(strsplit(log2fc_df$variable, '_'), function(x){x[2]}) %>% unlist()
log2fc_df$cell_cluster = lapply(strsplit(log2fc_df$variable, paste(log2fc_df$perc,'_',sep='')), 
                                function(x){x[2]}) %>% unlist()
log2fc_df %<>% group_by(class, cell_cluster) %>%
  mutate(mean = mean(value)) %>%
  dplyr::select(class, perc, cell_cluster, mean) %>%
  distinct() %>%
  mutate(cell_cluster = ifelse(substring(cell_cluster, 1, 2) == '1_', yes = substring(cell_cluster, 3), no = cell_cluster)) %>%
  spread(class, mean) %>% 
  mutate(log2fc = log2(`2`/`1`)) %>%
  dplyr::select(perc, cell_cluster, log2fc) %>%
  spread(perc, log2fc)
colnames(log2fc_df)[c(2,3)] = c('log2fc_parent', 'log2fc_total')

# Add log2fc into res_df
colnames(res_df)[which(colnames(res_df) == 'clusters')] = 'cell_cluster'
res_df %<>% left_join(log2fc_df) %>% 
  mutate(log2fc_parent_direction = ifelse(log2fc_parent > 0, yes = 'pos', no = 'neg')) %>%
  mutate(log2fc_total_direction = ifelse(log2fc_total > 0, yes = 'pos', no = 'neg'))

# 7. Boxplot of significant cell clusters ----------

cluster1 = sig_df %>% dplyr::filter(isTip == TRUE)
cluster2 = sig_df %>% dplyr::filter(isTip == FALSE)
cluster_sig = c(cluster1$clusters, as.character(cluster2$node)) 

for (i in 1:length(cluster_sig)) {
  cluster = cluster_sig[i]
  prop_df %>%
    dplyr::select(contains(c(cluster)), class, sample_id) %>%
    gather(variable, value, -c(class, sample_id)) %>%
    dplyr::filter(variable %in% c(paste0("perc_total_", cluster),
                                  paste0("perc_parent_1_", cluster))) %>%
    mutate(variable = factor(variable, 
                             levels = sapply(c("perc_total_", "perc_parent_1_"), 
                                             function(x) paste0(x, c(cluster))) %>% 
                               as.vector)) %>%
    ggplot(aes(x = class, y = value, col = class, fill = class)) +
    geom_boxplot(outlier.size = -1, width = 0.3, alpha = 0.2) +
    geom_jitter(position = position_jitter(seed = 1),alpha = 0.5) +
    geom_text(aes(label=sample_id), position = position_jitter(seed = 1), color='grey')+
    facet_wrap(~variable, scales="free", nrow = 2) +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(color = 'black')) +
    geom_signif(
      color = "grey40",
      comparisons = list(comparison_level),
      margin_top = 0.15,
      map_signif_level = function(p) sprintf("p = %.2g", p),
      test = "wilcox.test"
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background =element_blank(),
          axis.line = element_line(size = 0.35),
          legend.position = "right") +
    labs(title="Comparison between proportion to all and proportion to parent",
         subtitle="Proportions measured relative to parent cells per sample",
         fill= comparison, 
         col = comparison,
         y="Proportion") 
  ggsave(filename = paste0(cluster, '.pdf'), width = 10, height = 15)
}

#8. Save treekoR results ----------

write.csv(res_df, paste('treekoR_res_', title, '.csv', sep = ''), row.names = F)

