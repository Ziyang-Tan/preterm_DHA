library(FlowSOM)
library(flowCore)
library(dplyr)
library(mmR)
library(pheatmap)
library(paletteer)
library(ggpubr)
library(readr)

ff = read.FCS('data/cytof_outlier99_arcsinh5_combat_data.fcs', 
              transformation = FALSE, truncate_max_range = FALSE)

all_bind <- mmR::mm.fastread('data/flowSOM_all_data.csv.gz')
sample_info <- readxl::read_excel('data/210825_database_for_petter_sent_210826.xlsx') %>%
  select(ID_unique, Inf_ID, Timepoint_string, PN_days_sample, acttrt) %>%
  mutate(group = case_when(acttrt==1 ~ 'Control',
                           acttrt==2 ~ 'DHA_AA'),
         Timepoint_string = factor(Timepoint_string, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks')))

data <- all_bind[,colnames(ff)]
labels <- all_bind[,c('ID_unique', 'subtype1')] %>% left_join(sample_info, by='ID_unique')

subtype_index <- labels$subtype1 == 'T cell'

data <- scale(data)

data_sub <- data[subtype_index,]
labels_sub <- labels[subtype_index,]

ff@exprs <- data_sub
fsub.SOM = FlowSOM(ff,
                   compensate = F,
                   transform = F,
                   scale = F,
                   #colsToUse = setdiff(colnames(ff),c('IgD', 'CD1c', 'gdTCR', 'Siglec-8', 'CD20', 'CD14')), # exclude markers not for CD4T
                   colsToUse = colnames(ff),
                   nClus = 10,
                   xdim = 5, ydim = 6,
                   seed = 824)

# Median marker expression of FlowSOM cell clusters
fsub.res <- cbind(
  GetClusterMFIs(fsub.SOM, colsUsed = T), # median marker expression
  data.frame(clusters = GetClusters(fsub.SOM)) %>%
    dplyr::count(.data$clusters) %>%
    dplyr::mutate(percentage = .data$n / sum(.data$n) * 100)  # cluster and cluster frequency
)


# Heatmap of FlowSOM clusters
tmp <- fsub.res %>% select(-clusters,-n,-percentage) %>% as.matrix()
q99 = quantile(tmp, 0.99)

pheatmap::pheatmap(tmp[,c('CD3', 'CD4', 'CD8a', 'CD25', 'CD45RA', 'CD127', 'gdTCR', 'CD39')],
                   color = paletteer::paletteer_d("rcartocolor::PurpOr", n = 100, type = 'continuous'), # New color
                   border_color = 'white',
                   scale = 'none',
                   labels_row = paste(rownames(fsub.res),' (', round(fsub.res$percentage,1), '%', ')',sep = ''),
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   filename = 'figures/flowSOM_Tcell_for_subset_determination.pdf')

##
## manual annotaion here
##

annotation <- read_delim(file = 'data/flowSOM_T cell_annotation.csv', delim=';')
fsub.res_anno <- left_join(fsub.res, annotation, by='clusters')

pheatmap::pheatmap(tmp,
                   color = paletteer::paletteer_d("rcartocolor::PurpOr", n = 100, type = 'continuous'), # New color
                   border_color = 'white',
                   scale = 'none',
                   labels_row = paste(fsub.res_anno$subtype3, '_', fsub.res_anno$clusters,' (', round(fsub.res_anno$percentage,1), '%', ')',sep = ''),
                   display_numbers = FALSE,
                   angle_col = 45,
                   breaks = seq(0, q99, q99/90),
                   main = "Median marker expression per cluster (Z-score of all subtypes)",
                   #filename = 'figures/flowSOM_Tcell.pdf',
                   height = 8, width = 10)

write.table(fsub.res_anno, file = 'data/flowsom_clustered_Tcell.txt', 
            sep="\t", row.names=F, col.names = T, quote = F)

# cluster frequency

labels_sub <- cbind(labels_sub, tibble(clusters2 = GetClusters(fsub.SOM))) %>%
  left_join(annotation %>% rename(clusters2 = clusters), by='clusters2')

# save results for PAGA

mm.fastwrite(cbind(data_sub, labels_sub),  path = 'data/flowSOM_Tcell_data.csv', compress = T, overwrite = T)


#

raw <- mm.fastread('data/flowSOM_Tcell_data.csv.gz')

data_sub <- raw[,1:49]
labels_sub <- raw[,50:59]

data_freq <- labels_sub %>%
  filter(subtype2 != 'Exclude') %>%
  group_by(ID_unique, subtype2, subtype3) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(ID_unique) %>%
  mutate(freq = n/sum(n)) %>%
  ungroup() %>%
  #tidyr::complete(ID_unique, subtype3, fill = list(freq=0, n=0)) %>%
  left_join(labels %>% select(ID_unique, Inf_ID, Timepoint_string, PN_days_sample, group) %>% distinct(), by='ID_unique')


lapply(unique(data_freq$subtype3), function(x){
  d <- data_freq %>% filter(subtype3 == x)
  ggplot(data = d, 
         aes(x = Timepoint_string, y = freq, color = group)) +
    geom_line(aes(group=Inf_ID), alpha=0.3) +
    geom_point() +
    geom_smooth(aes(x = as.numeric(Timepoint_string), y = freq), method='loess') +
    labs(title = x) +
    theme_minimal()
}) %>%
  ggarrange(plotlist = ., ncol = 2, nrow = 2, common.legend = T) %>% 
  ggexport(filename = 'figures/T cell subtype3 freq changes.pdf',
           width = 10,
           height = 6)

data_freq_subtype2 <- data_freq %>% group_by(ID_unique, subtype2) %>% mutate(freq = sum(freq)) %>% select(-subtype3, -n) %>% distinct()

lapply(unique(data_freq_subtype2$subtype2), function(x){
  ggplot(data = data_freq_subtype2 %>% filter(subtype2 == x), 
         aes(x = Timepoint_string, y = freq, color = group)) +
    geom_line(aes(group=Inf_ID), alpha=0.3) +
    geom_point(alpha=0.3) +
    geom_smooth(aes(x = as.numeric(Timepoint_string), y = freq), method='loess') +
    labs(title = paste0(x, ' freq of all T cells')) +
    theme_minimal()
}) %>%
  ggarrange(plotlist = ., ncol = 2, nrow = 2, common.legend = T) %>% 
  ggexport(filename = 'figures/T cell subtype2 freq changes.pdf',
           width = 10,
           height = 6)
