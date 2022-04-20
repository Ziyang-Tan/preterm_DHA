library(SingleCellExperiment)
library(dplyr)
library(zellkonverter)
library(FlowSOM)
library(flowCore)
library(readr)
library(sva)
library(ggplot2)
library(ggpubr)
source('src/treekoR_src/analyseTree.R')
source('src/treekoR_src/visualiseTree.R')
'%nin%' <-Negate('%in%') 

sample_info <- read_csv('/Users/tan/preterm_DHA/data/refined_fileInfo.csv')
processed_data_path <- 'data/CyTOF_baby_fatty_acid_arcsinh_combat_eachSample50k.csv'
flowSOM_data_path <- 'data/CyTOF_baby_fatty_acid_flowSOM.csv'

if (file.exists(paste0(flowSOM_data_path, '.gz'))){
  data_flowSOM <- mmR::mm.fastread(paste0(flowSOM_data_path, '.gz'))
} else {
  if (file.exists(paste0(processed_data_path, '.gz'))){
    data_merge_combat <- mmR::mm.fastread(paste0(processed_data_path, '.gz'))
  } else {
    data_paths <- Sys.glob(file.path('/Users/tan/cytof_data/*/renamed/*', paste0(sample_info$Sample_ID, '.fcs')))
    label_paths <- Sys.glob(file.path('/Users/tan/cytof_data/*/classifiedV3/*', paste0(sample_info$Sample_ID, '.csv')))
    
    data_files <- lapply(data_paths, function(x){read.FCS(x)})
    names(data_files) <- sample_info$Sample_ID
    label_files <- lapply(label_paths, function(x){
      return(read.csv(x) %>% 
               tibble::add_column(Sample_ID = tools::file_path_sans_ext(basename(x)),
                                  batch = strsplit(x, '/')[[1]][5]) # EXP-XX-XXXXXX
      )
    })
    names(label_files) <- sample_info$Sample_ID
    
    # Remove empty channels
    channelFilter = function(flowframe){
      exprs = exprs(flowframe)[,-1]
      colnames(exprs) = markernames(flowframe)
      nonEmptyChannels = pData(parameters(flowframe))[c("name", "desc")] %>% dplyr::filter(name != desc)
      exprs = exprs[,colnames(exprs) %in% nonEmptyChannels$desc]
      return(exprs)
    }
    data = lapply(data_files, channelFilter)
    common_markers <- lapply(data, function(x){colnames(x)}) %>% Reduce(f=intersect)
    common_panel_markers <- common_markers[common_markers %nin% c("102Pd", "104Pd", "105Pd", "106Pd", "108Pd", "116Cd", "131Xe", "133Cs", "191Ir", "193Ir")]
    data <- lapply(data, function(x){x = x[, common_panel_markers]})
    
    # Arcsinh transformation
    data <- lapply(data, function(x){asinh(x/5)})
    
    # attach label to data
    data_label <- lapply(c(1:length(data)), function(x){
      cbind(data[[x]], label_files[[x]]) 
    })
    
    # relabel the cells, see details in doc/
    data_label <- lapply(data_label, function(x){
      x %>% dplyr::filter(level0 == 'cells') %>% #filter out non cells
        mutate(subpop = case_when(
          level1 %in% c("Neutrophils", "NK-cells", "Monocytes", "B-cells", "Eosinophils") ~ level1,
          level3 %in% c("CD39 Memory Tregs", "Memory Tregs", "Naive Tregs") ~ 'Tregs',
          level2 %in% c("pDC", "CD4 T", "CD8 T") ~ level2,
          TRUE ~ 'Others'
        ))
    })
    
    # subsample per sample to keep their cell counts the same
    data_label_trim <- lapply(data_label, function(x){
      full_size <- dim(x)[1]
      sample_size <- min(full_size, 50000)
      return(x[sample(full_size, sample_size),])
    })
    data_merge <- do.call(rbind, data_label_trim) %>% 
      left_join(sample_info %>% select(Sample_ID, Timepoint_string, randtrt), by = 'Sample_ID')
    data_matrix <- data_merge[, common_panel_markers] %>% as.matrix()
    data_combat <- sva::ComBat(t(data_matrix), batch = data_merge$batch, mod=model.matrix(~as.factor(Timepoint_string)+ as.factor(randtrt), data=data_merge)) %>% t()
    data_merge_combat = cbind(data_combat, data_merge[,colnames(data_merge) %nin% common_panel_markers])
    mmR::mm.fastwrite(data_merge_combat, path = processed_data_path,
                      compress = TRUE)
  }
  
  # FlowSOM
  # Run FlowSOM within each cell lineage --------------------------------------------------------------
  # read a fcs file as a flowframe container
  ff = read.FCS('/Users/tan/cytof_data/EXP-20-DG3616/renamed/1_con_0_cells/334 d0.fcs')
  # table of marker names and isotopes
  marker_df = ff@parameters@data 
  # parameters for ff
  parameter = data.frame(desc = unique(colnames(data_merge_combat[,1:49])),
                         range = 500, minRange = 0, maxRange = 499) %>% 
    left_join(marker_df %>% select(name, desc)) 
  ff@parameters@data = parameter
  # common markers 
  common_marker = colnames(data_merge_combat)[1:49] 
  # major cell pops by Grid
  grid_pops = unique(data_merge_combat$subpop)
  # run flowSOM for each major cell pop
  data_flowSOM_list <- lapply(grid_pops, function(pop){
    if (pop %in% c('CD4 T')) {nClus=10
    }else if (pop %in% c('pDC', 'Tregs')) {nClus=5
    }else {nClus=10}
    
    tmp = data_merge_combat %>% dplyr::filter(subpop == pop) 
    ff@exprs = tmp[,common_marker] %>% as.matrix()
    fSOM = FlowSOM(ff, 
                   compensate = F,
                   transform = F,
                   scale = F,
                   colsToUse = common_marker,
                   xdim = 4, ydim = 4,
                   nClus = nClus) 
    tmp$fSOM = GetMetaclusters(fSOM)
    filename = paste0('figures/flowSOM/FlowSOMmary_', pop, '.pdf')
    FlowSOMmary(fsom = fSOM, plotFile = filename)
    return(tmp)
  })
  
  data_flowSOM <- do.call(rbind, data_flowSOM_list)
  data_flowSOM <- data_flowSOM %>% mutate(cluster = paste0(subpop, '-', fSOM))
  
  mmR::mm.fastwrite(data_flowSOM, path = flowSOM_data_path,
                    compress = TRUE)
}

for (time_point in c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks')){
  #time_point = 'Day 1'
  #  data_dir = './PAGA_result_data/'
  #  sub_names <- sapply(dir(path=data_dir), 
  #                      function(x){strsplit(x,split = '_')[[1]][2]})
  #  nec_info <- readr::read_delim('/Users/tan/preterm_DHA/data/NEC_info.csv', delim = ';') %>%
  #    mutate(Inf_ID = as.character(Inf_ID))
  # data_merge <- lapply(sub_names, function(sub_name){
  #   file_path = file.path('./PAGA_result_data', 
  #                         paste0('all_', sub_name, '_sample1000_raw.h5ad'))
  #   datah5 = readH5AD(file = file_path)
  #   data = t(assay(datah5)) %>% as_tibble()
  #   meta = colData(datah5) %>% 
  #     as_tibble() %>% 
  #     select(timepoint, Inf_ID, group, leiden) %>%
  #     rename(randtrt = group) %>%
  #     left_join(nec_info, by='Inf_ID') %>%
  #     tidyr::replace_na(list(NEC = FALSE))
  #   
  #   data_sub <- data %>% filter(meta$timepoint == time_point)
  #   meta_sub <- meta %>% filter(timepoint == time_point) %>% mutate(cluster = paste0(sub_name, '-', leiden))
  #   return(list(data=data_sub, meta=meta_sub))
  # })
  # 
  # meta <- lapply(data_merge, function(x){
  #   return(x$meta)
  # }) %>% do.call(what=rbind) %>%
  #   mutate(cluster_inf = paste0(cluster, '-', Inf_ID))
  # data <- lapply(data_merge, function(x){
  #   return(x$data)
  # }) %>% do.call(what=dplyr::bind_rows) %>%
  #   replace(is.na(.), 0)
  
  # filter out clusters with low cell counts
  data <- data_flowSOM
  filter_out <- data %>% group_by(cluster) %>% summarise(count = n()) %>% dplyr::filter(count<200) %>% select(cluster) %>% unlist()
  data <- data %>% dplyr::filter(cluster %nin% filter_out)
  
  data <- data %>% dplyr::filter(Timepoint_string == time_point) %>%
    left_join(sample_info %>% select(-any_of(colnames(data)), Sample_ID), by = 'Sample_ID')
  # filter out clusters with only 1 subject
  filter_out2 <- data %>% group_by(cluster, randtrt) %>% summarise(unique_Subject = n_distinct(Subject_ID)) %>% dplyr::filter(unique_Subject==1) %>% select(cluster) %>% unlist()
  data <- data %>% dplyr::filter(cluster %nin% filter_out2)
  
  common_marker = colnames(data)[1:49] 
  exprs <- data[,common_marker]
  meta <- data %>% select(!all_of(common_marker)) 
  
  #clust_tree0 <- getClusterTree(exprs = exprs,
  #                             clusters = meta$cluster,
  #                             hierarchy_method = "hopach")
  ##### fix the tree-------------------------------------------------------------------------------------------
  # calculate median marker expression
  clust_med_dt = data.table::as.data.table(exprs)
  clust_med_dt[, cluster_id := meta$cluster]
  res = clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
  res2 = res[,.SD, .SDcols = !c('cluster_id')]
  rownames(res2) = res[["cluster_id"]]
  res_unscaled = res2
  res2[, (colnames(res2)) := lapply(.SD, scale), .SDcols=colnames(res2)]
  
  # label the position of cell clusters in the tree
  cell_level3 = unique(rownames(res2))
  cell_level2 = lapply(strsplit(cell_level3,'-'),function(x){x[1]}) %>% unlist()
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
    median_freq = res_unscaled,
    clust_tree = hc_phylo
  )
  
  #####----------------------------------------------------------------------------------------------------------------
  
  tested_tree <- testTree(phylo = clust_tree$clust_tree,
                          clusters = meta$cluster,
                          samples = meta$Subject_ID,
                          classes = meta$randtrt,
                          sig_test = 'ttest',
                          pos_class_name='2')
  plotInteractiveHeatmap(tested_tree,
                         clust_med_df = clust_tree$median_freq,
                         clusters=meta$cluster,
                         hm_offset=5)
  ggsave(filename = paste0('/Users/tan/preterm_DHA/figures/treeKoR_SOM/', time_point, '.pdf'), width = 14, height = 8)
  
  res_df = getTreeResults(tested_tree)
  head(res_df, 10)
  
  sig = res_df %>% mutate(sig = ifelse(pval_total>0.05 & pval_parent>0.05, yes='unsig',no='sig'))
  sig <- sig %>% dplyr::filter(sig == 'sig')
  if (dim(sig)[1]==0) next
  # Feature extraction
  prop_df = getCellProp(phylo=clust_tree$clust_tree,
                        clusters=meta$cluster,
                        samples=meta$Subject_ID,
                        classes=meta$randtrt) %>%
    mutate(class = as.factor(class))
  head(prop_df[,1:8])
  
  #means_df = getCellGMeans(clust_tree$clust_tree,
  #                         exprs=data,
  #                         clusters=meta$cluster,
  #                         samples=meta$Inf_ID,
  #                         classes=meta$randtrt)
  #head(means_df[,1:8])
  
  # Boxplot
  cluster1 = sig %>% dplyr::filter(isTip == TRUE)
  cluster2 = sig %>% dplyr::filter(isTip == FALSE)
  glist <- lapply(cluster1$clusters, function(cluster){
    prop_df %>%
      select(contains(c(cluster)), class, sample_id) %>%
      tidyr::gather(variable, value, -c(class, sample_id)) %>%
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
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background =element_blank(),
            axis.line = element_line(size = 0.35),
            legend.position = "right") +
      labs(title=paste0('total pval = ', signif((sig %>% dplyr::filter(clusters==cluster))$pval_total, digits = 2),
                        ', ',
                        'parent pval = ', signif((sig %>% dplyr::filter(clusters==cluster))$pval_parent, digits = 2)),
           fill="DHA_AA", 
           col = "DHA_AA",
           y="Proportion")
  })
  ggarrange(plotlist = glist, ncol = 3, nrow = 3) %>%
    ggexport(filename = paste0('/Users/tan/preterm_DHA/figures/treeKoR_SOM/', time_point, '_sig.pdf'), height = 20, width = 20)
}


