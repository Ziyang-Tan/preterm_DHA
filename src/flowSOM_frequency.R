library(dplyr)
library(ggplot2)
library(robCompositions)
library(ggpubr)
library(rstatix)


info_path <- file.path('data', 'sample_info_refined.csv')
data_path <- file.path('data', 'flowSOM_all_data.csv.gz')
data = mmR::mm.fastread(data_path)
olink <- readr::read_delim('data/Olink_breastmilk_group_72_inds.txt', delim = '\t') %>% rename(Subject_ID = Infant_ID)
group2 <- olink %>% select(Subject_ID, Breastmilk_group) %>% distinct()
fileInfo <- read.csv(info_path)

data_path <- Sys.glob(file.path('/Users/tan/cytof_data', '*', 'renamed', '*', paste0('*', fileInfo$ID_unique, '*.fcs'))) %>% unique()
fileInfo$cytof_batch <- sapply(stringr::str_split(data_path, '/'), function(x) x[5])

data_T <- mmR::mm.fastread('data/flowSOM_Tcell_data.csv.gz')

labels <- data[,50:51]
T_index <- which(labels$subtype1 == 'T cell')
labels$merge_subtype <- labels$subtype1
labels[T_index, 'merge_subtype'] <- data_T$subtype3

# QC
low_cell <- labels %>% group_by(ID_unique) %>% tally() %>% filter(n<3000)

# MDS by subtype annotations
subpop_freq <- labels %>% 
  filter(!ID_unique %in% low_cell$ID_unique) %>%
  group_by(merge_subtype, ID_unique) %>%
  tally() %>% 
  group_by(ID_unique) %>% 
  mutate(freq = n/sum(n))

tmp <- subpop_freq %>% select(-n) %>%
  tidyr::pivot_wider(names_from = merge_subtype, values_from = freq, values_fill = 1e-9) %>%
  ungroup()
m <- tmp %>% select(-ID_unique) %>% as.matrix()
datDist <- aDist(m)
datMDS = cmdscale(datDist)
colnames(datMDS) = c('MDS1', 'MDS2')
datMDS <- datMDS %>% 
  as_tibble() %>% 
  mutate(ID_unique = tmp$ID_unique) %>%
  left_join(fileInfo, by='ID_unique') %>%
  arrange(`Subject_ID`, `PN_days`) # for geom_path to connect in order

nameList <- c('group', 'SEX', 'WAYOFDE', 'Breastmilk_group', 'GA_days_at_birth', 'center', 'Death_PNA_days', 'GA_bin', 'PNA_bin', 'cytof_batch')
lapply(nameList, function(x){
  g <- ggplot(datMDS, aes_string(x='MDS1',y='MDS2',colour = paste0('`', x, '`'))) + 
    geom_point() +
    geom_path(aes_string(group='`Subject_ID`'))
  return(g)
}) %>% ggarrange(plotlist = ., ncol = 2, align = 'hv') %>%
  ggexport(filename = 'figures/MDS_by_subtype_annotations.pdf',
           width = 15, height = 6)

# frequency changes over time (group merged)
subpop_freq <- labels %>% 
  filter(!ID_unique %in% low_cell$ID_unique) %>%
  filter(merge_subtype != 'Exclude') %>%
  group_by(merge_subtype, ID_unique) %>%
  tally() %>% 
  group_by(ID_unique) %>% 
  mutate(freq = n/sum(n))
freq_anno <- left_join(subpop_freq, fileInfo, by='ID_unique')


subset_freq_change <- function(data){
  lapply(unique(data$merge_subtype), function(x){
    ggplot(data %>% 
             filter(merge_subtype == x) %>% 
             group_by(Timepoint_string) %>%
             mutate(median_freq = median(freq)), aes(x=Timepoint_string, y=freq)) +
      lemon::geom_pointline(aes(group=Subject_ID), position = position_jitter(), alpha=0.3) +
      geom_boxplot(aes(fill=median_freq), outlier.shape = NA) +
      viridis::scale_fill_viridis() +
      theme_bw() +
      ggtitle(x) +
      theme(axis.text.x = element_text(angle = -30, vjust = 0.5, hjust=0),
            axis.title.x = element_blank())
  }) %>%
    ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = TRUE, legend = "right", align = 'hv')
}

subset_freq_change(freq_anno) %>%
  ggexport(filename = 'figures/flowSOM_subset_frequency.pdf', width = 10, height = 10)
# separate groups
subset_freq_change(freq_anno %>% filter(group=='control')) %>%
  ggexport(filename = 'figures/flowSOM_subset_frequency_control.pdf', width = 10, height = 10)
subset_freq_change(freq_anno %>% filter(group=='DHA_AA')) %>%
  ggexport(filename = 'figures/flowSOM_subset_frequency_DHA_AA.pdf', width = 10, height = 10)

# compare the groups
lapply(unique(freq_anno$merge_subtype), function(x){
  df <- freq_anno %>% filter(merge_subtype == x)
  tb <- df %>% 
    group_by(Timepoint_string) %>% 
    t_test(freq ~ group) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj") %>% 
    add_xy_position(x = 'Timepoint_string')
  ggboxplot(df, x = "Timepoint_string", y = "freq", color = "group", palette = "jco", outlier.shape = NA) + 
    stat_pvalue_manual(tb, tip.length = 0) +
    geom_jitter(aes(x=Timepoint_string, y=freq, fill=group, color=group), position = position_jitterdodge()) +
    ggtitle(x) +
    theme(axis.text.x = element_text(angle = -30, vjust = 0.5, hjust=0),
          axis.title.x = element_blank())
}) %>%
  ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = TRUE, legend = "right", align = 'hv') %>%
  ggexport(filename = 'figures/flowSOM_subset_frequency_compare_group.pdf', width = 13, height = 13)

# permanova
lapply(unique(freq_anno$Timepoint_string), function(x){
  df <- freq_anno %>% filter(Timepoint_string == x) %>% 
    ungroup() %>%
    select(merge_subtype, ID_unique, freq) %>%
    tidyr::pivot_wider(names_from = merge_subtype, values_from = freq, values_fill = 1e-9) %>%
    left_join(fileInfo %>% select(ID_unique, SEX, group) %>% distinct(), by='ID_unique')
  subfreq <- df %>% select(-c("ID_unique", "SEX", "group"))
  tmp.stat <- vegan::adonis2(formula=subfreq ~ group + SEX, data=df[,c("ID_unique", "SEX", "group")], method = 'aitchison')
  return(tibble(var_name = rownames(tmp.stat), R2 = tmp.stat$R2, timepoint=x))
}) %>% do.call(what = rbind) %>%
  filter(var_name != 'Total') %>%
  ggplot(aes(x=timepoint, y=R2, fill=var_name)) +
  geom_bar(stat = 'identity', position = 'stack') +
  geom_text(aes(label=round(R2*100)), position=position_stack(vjust = 0.5)) +
  theme_classic()
ggsave(filename = 'figures/permanova_group_gender.pdf', width = 6, height = 6)


# postnatal age (environmental) vs gestational age (maturity)

subpop_freq_anno <- labels %>% 
  filter(merge_subtype != 'Exclude') %>%
  group_by(merge_subtype, ID_unique) %>%
  tally() %>% 
  group_by(ID_unique) %>% 
  mutate(freq = n/sum(n)) %>% select(-n) %>%
  ungroup() %>%
  left_join(fileInfo %>% select(ID_unique, Subject_ID, GA_days, PN_days, group, Timepoint_string, GA_bin, PNA_bin), by='ID_unique')


# by PNA

d <- subpop_freq_anno 

lapply(unique(d$merge_subtype), function(x){
  ggplot(d %>% 
           filter(merge_subtype == x) %>% 
           group_by(PNA_bin) %>%
           mutate(median_freq = median(freq)), aes(x=PNA_bin, y=freq)) +
    lemon::geom_pointline(aes(group=Subject_ID), position = position_jitter(), alpha=0.3) +
    geom_boxplot(aes(fill=median_freq), outlier.shape = NA) +
    viridis::scale_fill_viridis() +
    theme_bw() +
    ggtitle(x) +
    theme(axis.text.x = element_text(angle = -30, vjust = 0.5, hjust=0),
          axis.title.x = element_blank())
}) %>%
  ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = TRUE, legend = "none", align = 'hv') %>%
  ggexport(filename = 'figures/PMA/subset_frequency_PN days.pdf', width = 10, height = 10)

lapply(unique(d$merge_subtype), function(x){
  ggplot(d %>% 
           filter(merge_subtype == x) %>% 
           group_by(GA_bin) %>%
           mutate(median_freq = median(freq)), aes(x=GA_bin, y=freq)) +
    lemon::geom_pointline(aes(group=Subject_ID), position = position_jitter(), alpha=0.3) +
    geom_boxplot(aes(fill=median_freq), outlier.shape = NA) +
    viridis::scale_fill_viridis() +
    theme_bw() +
    ggtitle(x) +
    theme(axis.text.x = element_text(angle = -30, vjust = 0.5, hjust=0),
          axis.title.x = element_blank())
}) %>%
  ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = TRUE, legend = "none", align = 'hv') %>%
  ggexport(filename = 'figures/PMA/subset_frequency_GA.pdf', width = 10, height = 10)

# compare the groups based on new binning
lapply(unique(d$merge_subtype), function(x){
  df <- d %>% filter(merge_subtype == x)
  tb <- df %>% 
    group_by(PNA_bin) %>% 
    t_test(freq ~ group) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj") %>% 
    add_xy_position(x = 'PNA_bin')
  ggboxplot(df, x = "PNA_bin", y = "freq", color = "group", palette = "jco", outlier.shape = NA) + 
    stat_pvalue_manual(tb, tip.length = 0) +
    geom_jitter(aes(x=PNA_bin, y=freq, fill=group, color=group), position = position_jitterdodge()) +
    ggtitle(x) +
    theme(axis.text.x = element_text(angle = -30, vjust = 0.5, hjust=0),
          axis.title.x = element_blank())
}) %>%
  ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = TRUE, legend = "right", align = 'hv') %>%
  ggexport(filename = 'figures/PMA/compare_group_PNA.pdf', width = 13, height = 13)

lapply(unique(d$merge_subtype), function(x){
  df <- d %>% filter(merge_subtype == x)
  tb <- df %>% 
    group_by(GA_bin) %>% 
    t_test(freq ~ group) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj") %>% 
    add_xy_position(x = 'GA_bin')
  ggboxplot(df, x = "GA_bin", y = "freq", color = "group", palette = "jco", outlier.shape = NA) + 
    stat_pvalue_manual(tb, tip.length = 0) +
    geom_jitter(aes(x=GA_bin, y=freq, fill=group, color=group), position = position_jitterdodge()) +
    ggtitle(x) +
    theme(axis.text.x = element_text(angle = -30, vjust = 0.5, hjust=0),
          axis.title.x = element_blank())
}) %>%
  ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = TRUE, legend = "right", align = 'hv') %>%
  ggexport(filename = 'figures/PMA/compare_group_GA.pdf', width = 13, height = 13)
