library(FlowSOM)
library(flowCore)
library(dplyr)
library(mmR)
library(pheatmap)
library(paletteer)
library(ggpubr)
library(readr)
library(anndata)
library(umap)

#
ff = read.FCS('data/cytof_outlier99_arcsinh5_combat_data.fcs', 
              transformation = FALSE, truncate_max_range = FALSE)
data = as_tibble(ff@exprs)
labels = mm.fastread('data/cytof_outlier99_arcsinh5_combat_label.csv.gz', header = F)
colnames(labels) <- c('level1', 'level2', 'level3', 'label_path', 'ID_unique')

sample_info <- readxl::read_excel('data/210825_database_for_petter_sent_210826.xlsx') %>%
  select(ID_unique, Inf_ID, Timepoint_string, PN_days_sample, acttrt) %>%
  mutate(group = case_when(acttrt==1 ~ 'Control',
                           acttrt==2 ~ 'DHA_AA'),
         Timepoint_string = factor(Timepoint_string, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks')))

data_anno <- cbind(data,labels)
data_sample <- data_anno %>% group_by(level1) %>% sample_n(5000)
pca.res <- prcomp(scale(data_sample[,1:49]))
data.umap <- umap(pca.res$x[,1:20])
d <- cbind(data.umap$layout, data_sample) %>% rename(umap1 = `...1`, umap2 = `...2`)
ggplot(d %>% filter(level1 != 'Other'), aes(x=umap1, y=umap2, color=level1)) +
  geom_point(size=0.2)

# 
all_bind <- mmR::mm.fastread('data/flowSOM_all_data.csv.gz')
data <- all_bind[,colnames(ff)]
labels <- all_bind[,c('ID_unique', 'subtype1')] %>% left_join(sample_info, by='ID_unique')
grid_label = readr::read_csv('data/cytof_outlier99_arcsinh5_combat_label.csv.gz', 
                             col_names = F)
colnames(grid_label) <- c('level1', 'level2', 'level3', 'label_path', 'ID_unique')
data_anno <- cbind(data,grid_label)
data_sample <- data_anno %>% group_by(level1) %>% sample_n(5000)
pca.res <- prcomp(scale(data_sample[,1:49]))
data.umap <- umap(pca.res$x[,1:20])
d <- cbind(data.umap$layout, data_sample) %>% rename(umap1 = `...1`, umap2 = `...2`)
ggplot(d %>% filter(level1 != 'Other'), aes(x=umap1, y=umap2, color=level1)) +
  geom_point(size=0.2)

## all_bind has been permutated. We fish out Neutrophils and then append them in the end so the order is different

#
all_bind <- mmR::mm.fastread('data/flowSOM_all_data.csv.gz')
data <- all_bind[,1:49]
labels <- all_bind[,c('ID_unique', 'subtype1')] %>% left_join(sample_info, by='ID_unique')
data_anno <- cbind(data, labels)
data_sample <- data_anno %>% group_by(ID_unique) %>% sample_n(1000)
pca.res <- prcomp(scale(data_sample[,1:49]))
data.umap <- umap(pca.res$x[,1:20])
d <- cbind(data.umap$layout, data_sample) %>% 
  rename(umap1 = `...1`, umap2 = `...2`) %>% 
  mutate(subtype1 = as.factor(subtype1))

library(RColorBrewer)
myColors <- brewer.pal(nlevels(d$subtype1),"Set3")
names(myColors) <- levels(d$subtype1)

ggplot(d, aes(x=umap1, y=umap2, color=subtype1)) +
  geom_point(size=0.8) +
  scale_colour_manual(name = "subtype1",values = myColors)
ggsave('figures_project_rescue/flowSOM_umap_all.pdf')

ggplot(d, aes(x=umap1, y=umap2, color=CD3)) +
  geom_point(size=0.8)
ggsave('figures_project_rescue/flowSOM_umap_all_CD3.pdf')

l <- data_anno %>% group_by(ID_unique) %>% tally() %>% arrange(n)
l1 <- l %>% arrange(n) %>% slice(1:20)
l2 <- l %>% arrange(desc(n)) %>% slice(1:20)

lapply(l1$ID_unique, function(x){
  ggplot(d %>% filter(ID_unique == x), aes(x=umap1, y=umap2, color=subtype1)) +
    geom_point(size=0.8) +
    scale_colour_manual(name = "subtype1",values = myColors) +
    ggtitle(x)
}) %>%
  ggarrange(plotlist = ., ncol = 4, nrow = 5, common.legend = T) %>%
  ggexport(filename = 'figures_project_rescue/flowSOM_umap_per_sample_low_count.pdf', width = 14, height = 18)

lapply(l2$ID_unique, function(x){
  ggplot(d %>% filter(ID_unique == x), aes(x=umap1, y=umap2, color=subtype1)) +
    geom_point(size=0.8) +
    scale_colour_manual(name = "subtype1",values = myColors) +
    ggtitle(x)
}) %>%
  ggarrange(plotlist = ., ncol = 4, nrow = 5, common.legend = T) %>%
  ggexport(filename = 'figures_project_rescue/flowSOM_umap_per_sample_high_count.pdf', width = 14, height = 18)


lapply(l$ID_unique, function(x){
  ggplot(d %>% filter(ID_unique == x), aes(x=umap1, y=umap2, color=subtype1)) +
    geom_point(size=0.8) +
    scale_colour_manual(name = "subtype1",values = myColors) +
    ggtitle(x)
}) %>%
  ggarrange(plotlist = ., ncol = 4, nrow = 5, common.legend = T) %>%
  ggexport(filename = 'figures_project_rescue/flowSOM_umap_per_sample_all.pdf', width = 14, height = 18)






#files_fcs = lapply(data_path, function(x){read.FCS(x, truncate_max_range = FALSE)})
#files_csv = lapply(grid_path, function(x){read.csv(x, row.names = 1)})

anndata <- read_h5ad('/Users/tan/preterm_DHA/adata/all_data_preterm_DHA.h5ad')

df <- as_tibble(anndata$X)
label <- as_tibble(anndata$obs)

data_anno <- cbind(df, label %>% select(level1, ID_unique))
data_sample <- data_anno %>% group_by(ID_unique) %>% sample_n(500)

pca.res <- prcomp(data_sample[,1:49])
data.umap <- umap(pca.res$x[,1:20])
d <- cbind(data.umap$layout, data_sample) %>% rename(umap1 = `...1`, umap2 = `...2`)
ggplot(d %>% filter(umap2<10), aes(x=umap1, y=umap2, color=level1)) +
  geom_point(size=0.2)
ggsave('figures_project_rescue/umap_500perSample_from_adata_all.pdf')


# from script 'CyTOF_preprocessing.R'
data_path <- Sys.glob(file.path('/Users/tan/cytof_data', '*', 'renamed', '*', paste0('*', sample_info$ID_unique, '*.fcs'))) %>% unique()
grid_path <- Sys.glob(file.path('/Users/tan/cytof_data', '*', 'ClassifiedV3', '*', paste0('*', sample_info$ID_unique, '*.csv'))) %>% unique()

files_fcs = lapply(data_path, function(x){read.FCS(x, truncate_max_range = FALSE)})
files_csv = lapply(grid_path, function(x){
  tmp <- read.csv(x, row.names = 1)
  tmp$ID_unique <- sub('.csv', '', sub('^.*/', '', x))
  return(tmp)
})

channelFilter = function(flowframe){
  exprs = exprs(flowframe)[,-1]
  colnames(exprs) = markernames(flowframe)
  nonEmptyChannels = pData(parameters(flowframe))[c("name", "desc")] %>% dplyr::filter(name != desc)
  exprs = exprs[,colnames(exprs) %in% nonEmptyChannels$desc]
  return(exprs)
}

data <- lapply(files_fcs, channelFilter) %>% do.call(what = rbind)
data <- as_tibble(data)

label <- files_csv %>% do.call(what = rbind)
label <- as_tibble(label)

data_anno <- cbind(data %>% select(-c("102Pd", "104Pd", "105Pd", "106Pd", "108Pd",
                                      "116Cd", "131Xe", "133Cs", "191Ir", "193Ir")), 
                   label %>% select(level1, ID_unique))

data_sample <- data_anno %>% group_by(ID_unique) %>% sample_n(500)
pca.res <- prcomp(data_sample[,1:49], scale.=T)
data.umap <- umap(pca.res$x[,1:20])
d <- cbind(data.umap$layout, data_sample) %>% rename(umap1 = `...1`, umap2 = `...2`)
ggplot(d %>% filter(umap2>-10), aes(x=umap1, y=umap2, color=level1)) +
  geom_point(size=0.2)

#

df <- data_T %>% 
  filter(subtype3 != 'Exclude',
         subtype3 != 'gdT') %>% 
  group_by(ID_unique, subtype3) %>%
  tally() %>% 
  group_by(ID_unique) %>%
  mutate(freq = n/sum(n)) %>%
  left_join(fileInfo, by='ID_unique')

ggplot(df%>% filter(subtype3=='gdT'), aes(x=log(PN_days+1), y=freq, color=group, group=group)) +
  geom_point() +
  geom_smooth()

