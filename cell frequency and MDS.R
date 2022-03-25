library(stats)
library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(robCompositions)
library(ggpubr)

info_path <- file.path('data', '210825_database_for_petter_sent_210826.xlsx')
nec_info <- readr::read_delim('/Users/tan/preterm_DHA/data/NEC_info.csv', delim = ';')
fileInfo = readxl::read_excel(info_path, sheet=1) %>%
  select(-Sample_ID) %>%
  left_join(nec_info, by = 'Inf_ID') %>%
  tidyr::replace_na(list(NEC = FALSE)) %>%
  rename(Sample_ID = ID_unique,
         timepoint = PN_days_sample, 
         Subject_ID = Inf_ID) %>%
  mutate(randtrt = as.factor(randtrt)) # 1: control 2: DHA
  

projs<- c('EXP-20-DG3616', 'EXP-20-DG3619', 'EXP-20-DG3620', 'EXP-20-DG3621')
data_paths <- file.path('/Users/tan/cytof_data/', projs, 'classifiedV3', 'abundance.csv')
raw <- lapply(data_paths, function(x){
  read_csv(x, show_col_types = FALSE) %>%
    mutate(subpop = paste(...1,...2,sep='/')) %>%
    select(-c('...1', '...2')) %>%
    select(-contains('Nil'))
}) %>% Reduce(f = function(...){left_join(..., by = 'subpop')})


tmp <- raw %>% select(-subpop) %>% t()
colnames(tmp) = raw$subpop
dat <- tmp %>% as_tibble() %>% add_column(Sample_ID = rownames(tmp)) %>%
  replace(is.na(.), 0)
#mutate(`abT-cells/CD4 T/` = `abT-cells/CD4 T/` + `abT-cells/CD4 T`) %>% # fix the abundance matrix
#select(-c('abT-cells/CD4 T'))

# save processed tables
write_csv(dat, '/Users/tan/preterm_DHA/data/refined_frequency_data.csv')
write_csv(fileInfo, '/Users/tan/preterm_DHA/data/refined_fileInfo.csv')

# changes of T reg
subpopDat <- dat %>%
  select(contains('abT')) %>%
  select(contains('CD4')) %>%
  mutate(`T regs` = rowSums(select(.,contains('Tregs')))) %>%
  select(-contains('Tregs')) %>%
  mutate(`T regs of CD4` = `T regs`/rowSums(.)) %>%
  add_column(Sample_ID = dat$Sample_ID) %>%
  tidyr::pivot_longer(!Sample_ID, names_to = 'subpop', values_to = 'frequency') %>%
  left_join(fileInfo, by='Sample_ID')

# plot Treg freq changes of randtrt
g1 <- ggplot(subpopDat %>% filter(subpop == 'T regs of CD4'),
             aes(x = log(timepoint+1), y = frequency, color = randtrt)) +
  geom_point() +
  geom_line(aes(group = Subject_ID)) +
  geom_smooth(aes(x = log(as.numeric(timepoint)+1), y = frequency), method='loess') +
  stat_cor(aes(x = log(as.numeric(timepoint)+1), y = frequency), 
           label.x = 3,
           label.x.npc = "right",
           method = 'spearman') +
  ylab('T reg frequency of CD4 T cells') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggarrange(g1, ncol=1) %>%
  ggexport(filename = '/Users/tan/preterm_DHA/figures/frequency changes/Treg changes_randtrt.pdf',
           width = 10, height = 10)

# plot Treg freq changes of nec
g1 <- ggplot(subpopDat %>% filter(subpop == 'T regs of CD4'),
             aes(x = log(timepoint+1), y = frequency, color = NEC)) +
  geom_point() +
  geom_line(aes(group = Subject_ID)) +
  geom_smooth(aes(x = log(as.numeric(timepoint)+1), y = frequency), method='loess') +
  stat_cor(aes(x = log(as.numeric(timepoint)+1), y = frequency), 
           label.x = 3,
           label.x.npc = "right",
           method = 'spearman') +
  ylab('T reg frequency of CD4 T cells') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggarrange(g1, ncol=1) %>%
  ggexport(filename = '/Users/tan/preterm_DHA/figures/frequency changes/Treg changes_nec.pdf',
           width = 10, height = 10)

# plot Treg freq changes of nec x randtrt
g1 <- ggplot(subpopDat %>% filter(subpop == 'T regs of CD4'),
             aes(x = log(timepoint+1), y = frequency, color = NEC, linetype = randtrt)) +
  geom_point() +
  geom_line(aes(group = Subject_ID)) +
  geom_smooth(aes(x = log(as.numeric(timepoint)+1), y = frequency), method='loess') +
  stat_cor(aes(x = log(as.numeric(timepoint)+1), y = frequency), 
           label.x = 3,
           label.x.npc = "right",
           method = 'spearman') +
  ylab('T reg frequency of CD4 T cells') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggarrange(g1, ncol=1) %>%
  ggexport(filename = '/Users/tan/preterm_DHA/figures/frequency changes/Treg changes_nec_n_randtrt.pdf',
           width = 10, height = 10)

ggplot(subpopDat %>% 
         filter(Timepoint_string == 'Day 7') %>%
         filter(subpop == 'T regs of CD4'),
       aes(x = randtrt, y = frequency, fill = randtrt)) +
  geom_violin() +
  geom_point() +
  labs(title='Day 7, p=0.02519')

ggplot(subpopDat %>% 
         filter(Timepoint_string == 'Day 14') %>%
         filter(subpop == 'T regs of CD4'),
       aes(x = randtrt, y = frequency, fill = randtrt)) +
  geom_violin() +
  geom_point() +
  labs(title='Day 14, p=0.005956')

ggplot(subpopDat %>% 
         filter(Timepoint_string == 'Day 28') %>%
         filter(subpop == 'T regs of CD4'),
       aes(x = randtrt, y = frequency, fill = randtrt)) +
  geom_violin() +
  geom_point() +
  labs(title='Day 28, p=0.04828')

# changes of other subpop

sub_names = c('B-cells', 'CD4 T', 'CD8 T', 'Eosinophils', 'Monocytes', 'Neutrophils', 
              'NK-cells', 'pDC', 'Plasmablasts', 'Basophils')
g_list <- lapply(sub_names, function(sub_name){
  subpopDat <- dat %>%
    select(contains(sub_name)) %>%
    mutate(current_sub = rowSums(select(.,contains(sub_name)))) %>%
    add_column(Sample_ID = dat$Sample_ID) %>%
    tidyr::pivot_longer(!Sample_ID, names_to = 'subpop', values_to = 'frequency') %>%
    left_join(fileInfo, by='Sample_ID')
  
  ggplot(subpopDat %>% filter(subpop == 'current_sub'),
         aes(x = log(timepoint+1), y = frequency, color = randtrt)) +
    geom_point() +
    geom_line(aes(group = Subject_ID)) +
    geom_smooth(aes(x = log(as.numeric(timepoint)+1), y = frequency), method='loess') +
    stat_cor(aes(x = log(as.numeric(timepoint)+1), y = frequency), 
             label.x = 3,
             label.x.npc = "right",
             method = 'spearman') +
    ylab(paste0(sub_name, ' frequency')) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
})
ggarrange(plotlist = g_list, ncol = 2, nrow = 3) %>%
  ggexport(filename = '/Users/tan/preterm_DHA/figures/frequency changes/frequency changes.pdf',
           width = 12, height = 10)

#ggplot(d, aes(fill=subpop, y=frequency, x=timepoint)) + 
#  geom_bar(position="fill", stat="identity")

# MDS
m <- (dat %>% select(-Sample_ID) %>% replace(.==0, 1e-9) %>% as.matrix())
datDist <- aDist(m)
datMDS = cmdscale(datDist)
colnames(datMDS) = c('MDS1', 'MDS2')
datMDS <- datMDS %>% 
  as_tibble() %>% 
  add_column(Sample_ID = dat$Sample_ID) %>%
  left_join(fileInfo, by='Sample_ID') %>%
  arrange(`Inf_ID`, `timepoint`)

nameList <- c('Timepoint_string', 'randtrt', 'timepoint')
gList <- lapply(nameList, function(x){
  g <- ggplot(datMDS, aes_string(x='MDS1',y='MDS2',colour = paste0('`', x, '`'))) + 
    geom_point() +
    geom_path(aes_string(group='`Inf_ID`'))
  return(g)
})
ggarrange(plotlist = gList, ncol = 2, nrow = 2, align = 'hv')
ggexport(filename = '/Users/tan/Ionctura-collab/results/cell_frequency_mds.pdf',
         width = 12, height = 6)


