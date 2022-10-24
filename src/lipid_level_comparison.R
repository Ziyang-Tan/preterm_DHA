library(dplyr)
library(ggplot2)
library(ggpubr)

figurePath <- '/Users/tan/preterm_DHA/figures/lipids_related/'

lipid_data <- xlsx::read.xlsx('/Users/tan/preterm_DHA/data/AUC fatty acids and nutrition for Petter 2022-02-22.xlsx',
                              sheetIndex = 1,
                              password = 'Formulaid') %>% 
  rename(Subject_ID = Inf_ID)
# normalize by mean daily lipids
lipid_data_normed <- lipid_data %>% mutate_at(vars(contains('AUC_d1to28')), ~./Lipids_mean)

fileInfo <- readr::read_csv('/Users/tan/preterm_DHA/data/refined_fileInfo.csv') %>%
  mutate(randtrt = case_when(
    randtrt==1 ~ 'control',
    randtrt==2 ~ 'AA_DHA'
  )) %>%
  left_join(lipid_data_normed, by = 'Subject_ID')
#left_join(lipid_data, by = 'Subject_ID')
#data <- readr::read_csv('/Users/tan/preterm_DHA/data/refined_frequency_data.csv') %>%
#  left_join(fileInfo, by = 'Sample_ID')










frequency_cols <- grep('/', colnames(data), value=TRUE)
lipid_cols <- grep('AUC_d1to28', colnames(data), value=TRUE)

# correlation heatmap
cor_data_treated <- cor(data %>% filter(randtrt=='AA_DHA') %>% select(frequency_cols),
                        data %>% filter(randtrt=='AA_DHA') %>% select(lipid_cols),
                        use = 'complete.obs')

cor_data_control <- cor(data %>% filter(randtrt=='control') %>% select(frequency_cols),
                        data %>% filter(randtrt=='control') %>% select(lipid_cols),
                        use = 'complete.obs')

cor_data_all <- cor(data %>% select(frequency_cols),
                    data %>%  select(lipid_cols),
                    use = 'complete.obs')

library(ComplexHeatmap)
pdf(file = file.path(figurePath, 'cor_heatmap_treated_normed.pdf'), width = 16, height = 12)
Heatmap(cor_data_treated)
dev.off()
pdf(file = file.path(figurePath, 'cor_heatmap_control_normed.pdf'), width = 16, height = 12)
Heatmap(cor_data_control)
dev.off()
pdf(file = file.path(figurePath, 'cor_heatmap_all_normed.pdf'), width = 16, height = 12)
Heatmap(cor_data_all)
dev.off()


# fatty acid level
lipid_data <- data %>% select(lipid_cols, randtrt, Subject_ID) %>% distinct()
names(lipid_cols) <- lipid_cols
pvals <- lapply(lipid_cols, function(x){
  tmp <- t.test((lipid_data%>%filter(randtrt=='control'))[x], (data%>%filter(randtrt=='AA_DHA'))[x])
  return(tmp$p.value)
})
p.adj <- p.adjust(unlist(pvals), method='BH')

fatty_acid_g_list <- lapply(lipid_cols, function(lipid_name){
  ggplot(lipid_data, aes_string(x='randtrt', y=lipid_name, group='randtrt', color='randtrt')) +
    geom_violin() +
    geom_jitter() +
    theme(legend.position = "none") +
    labs(title=paste0('p.adj = ', signif(p.adj[lipid_name], digits = 2)))
})
ggarrange(plotlist = fatty_acid_g_list, ncol = 4, nrow = 10) %>%
  ggexport(filename = file.path(figurePath, 'fatty_acid_levels_normed.pdf'), width = 10, height = 20)
