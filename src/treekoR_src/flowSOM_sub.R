library(FlowSOM)
library(flowCore)
library(tidyverse)
library(magrittr)

ff = read.FCS('data/cytof_outlier99_arcsinh5_combat_data.fcs', 
              transformation = FALSE, truncate_max_range = FALSE)
label = readr::read_csv('data/cytof_outlier99_arcsinh5_combat_label.csv.gz', 
                        col_names = F) 
sample_info = readxl::read_excel('data/210825_database_for_petter_sent_210826.xlsx') 
annotation = read.csv2('data/flowSOM100_annotation.csv', check.names = F, stringsAsFactors = F)
all_bind = mmR::mm.fastread('data/flowSOM_all_data.csv.gz')

data_sub <- all_bind %>% filter(lineage == 'T cell')
ff@exprs <- as.matrix(data_sub[,1:49])

fSOM_sub = FlowSOM(ff,
                   compensate = F,
                   transform = F,
                   scale = F,
                   #colsToUse = common_marker,
                   nClus = 10,
                   xdim = 5, ydim = 6,
                   seed = 824)
FlowSOMmary(fSOM_sub, plotFile = 'flowsom_summary_Tcell.pdf')

MFI_sub = GetClusterMFIs(fSOM_sub, colsUsed = T)
clusters_sub = GetClusters(fSOM_sub)
freqClusters_sub = data.frame(clusters = clusters_sub) %>%
  dplyr::count(.data$clusters) %>%
  dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
  as.data.frame()
res_sub = cbind(MFI_sub, freqClusters_sub)

q99 = quantile(as.matrix(MFI_sub), 0.99)
pheatmap::pheatmap(as.matrix(MFI_sub[,c('CD4', 'CD8a', 'gdTCR', 'CD25', 'CD127', 'pan-KIR')]),
                   scale = 'column',
                   labels_row = paste(rownames(MFI_sub),' (', round(freqClusters_sub$percentage,1), '%', ')',sep = ''),
                   display_numbers = TRUE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "all T cells")

data_sub$T_cluster <- clusters_sub
tmp <- data_sub %>% mutate(
  T_cell_subtype_1 = case_when(
    T_cluster %in% c(19, 15, 20, 7, 28, 10, 9, 24, 17, 23, 12) ~ 'CD4T',
    T_cluster %in% c(4, 29, 25, 30, 21) ~ 'CD8T',
    T_cluster %in% c(14, 26, 18, 27, 13, 8, 22, 5) ~ 'DNT',
    T_cluster %in% c(6, 16) ~ 'DPT',
    T_cluster %in% c(3, 2, 11, 1) ~ 'gdT'
  ),
  T_cell_subtype_2 = case_when(
    T_cluster %in% c(9) ~ 'CD4Treg',
  )
)

tmp2 <- tmp %>% filter(T_cell_subtype_1 == 'CD4T') %>% mutate(ID_unique == as.factor(ID_unique))

Treg_freq = tmp2 %>% group_by(ID_unique) %>% tally(name = 'CD4_count') %>% left_join(
  tmp2 %>% filter(T_cell_subtype_2 == 'CD4Treg') %>% group_by(ID_unique) %>% tally(name = 'Treg_count'),
  by = 'ID_unique'
) %>% replace_na(list(Treg_count = 0)) %>%
  mutate(Treg_freq = Treg_count/CD4_count) 

ggplot(
  data = sample_info %>% 
    left_join(Treg_freq, by='ID_unique') %>%
    select(ID_unique, PN_days_sample, acttrt, Treg_freq, Inf_ID),
  aes(x=log(PN_days_sample+1), y=Treg_freq, group=Inf_ID, color = acttrt)
) +
  geom_point() +
  geom_smooth(method = 'loess', aes(group=acttrt))

all_T_flowsom <- tmp

mmR::mm.fastwrite(all_T_flowsom, 'data/all_T_flowsom.csv.gz')










