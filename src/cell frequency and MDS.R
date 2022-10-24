library(stats)
library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(robCompositions)
library(ggpubr)
library(rstatix)

info_path <- file.path('data', '210825_database_for_petter_sent_210826.xlsx')
nec_info <- readr::read_delim('/Users/tan/preterm_DHA/data/NEC_info.csv', delim = ';')
fileInfo = readxl::read_excel(info_path, sheet=1) %>%
  select(-Sample_ID) %>%
  left_join(nec_info, by = 'Inf_ID') %>%
  tidyr::replace_na(list(NEC = FALSE)) %>%
  dplyr::rename(Sample_ID = ID_unique,
         timepoint = PN_days_sample, 
         Subject_ID = Inf_ID) %>%
  mutate(group = as.factor(case_when(
    acttrt == 1 ~ 'control',
    acttrt == 2 ~ 'DHA_AA'
  )),
  Timepoint_string = factor(Timepoint_string, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks'))) # 1: control 2: DHA


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

d <- subpopDat %>% filter(subpop == 'T regs of CD4')

ggplot(d, aes(x = Timepoint_string, y = frequency, color = group)) +
  geom_point(alpha=0.2) +
  geom_line(aes(group = Subject_ID), alpha=0.2) +
  geom_smooth(aes(x = as.numeric(Timepoint_string), y = frequency), method='loess') +
  theme_minimal() + 
  ggtitle('T regs of CD4')
ggsave('figures/manu/raw/Treg of CD4 frequency changes.pdf')

# two way mixed anova (from Day 7 to Day 28)
res.aov <- anova_test(
  data = d %>% filter(!Timepoint_string %in% c('PMA 40 weeks')), 
  dv = frequency, wid = Subject_ID,
  between = group, 
  within = Timepoint_string
)
write.csv(get_anova_table(res.aov), 'figures/manu/raw/Treg of CD4 anova (Day 1 to 28).csv')
# pair wise t test
tb <- d %>% 
  group_by(Timepoint_string) %>% 
  pairwise_t_test(
    frequency ~ group, p.adjust.method = 'BH'
  ) %>% add_xy_position(x = 'Timepoint_string')
ggboxplot(
  d, x = "Timepoint_string", y = "frequency",
  color = "group", palette = "jco"
) + stat_pvalue_manual(tb, tip.length = 0, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE))
ggsave(filename = 'figures/manu/raw/Treg of CD4 boxplot.pdf')
write_csv(tb, 'figures/manu/raw/Treg of CD4 pairwise t.csv')

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
  d <- subpopDat %>% filter(subpop == 'current_sub')
  g <- ggplot(d,
              aes(x = Timepoint_string, y = frequency, color = group)) +
    geom_point(alpha=0.2) +
    geom_line(aes(group = Subject_ID), alpha=0.2) +
    geom_smooth(aes(x = as.numeric(Timepoint_string), y = frequency), method='loess') +
    ggtitle(sub_name)+
    theme_minimal()
  res.aov <- anova_test(
    data = d %>% filter(!Timepoint_string %in% c('PMA 40 weeks')), 
    dv = frequency, wid = Subject_ID,
    between = group, 
    within = Timepoint_string
  )
  write.csv(get_anova_table(res.aov), paste0('figures/manu/raw/', sub_name, ' anova (Day 1 to 28).csv'))
  tb <- d %>% 
    group_by(Timepoint_string) %>% 
    pairwise_t_test(
      frequency ~ group, p.adjust.method = 'BH'
    ) %>% add_xy_position(x = 'Timepoint_string')
  g_tmp <- ggboxplot(
    d, x = "Timepoint_string", y = "frequency",
    color = "group", palette = "jco"
  ) + stat_pvalue_manual(tb, tip.length = 0, hide.ns = TRUE) +
    labs(title = sub_name,
         subtitle = get_test_label(res.aov, detailed = TRUE))
  ggsave(g_tmp, filename = paste0('figures/manu/raw/', sub_name, ' boxplot.pdf'))
  write_csv(tb, paste0('figures/manu/raw/', sub_name, ' pairwise t.csv'))
  return(g)
})
ggarrange(plotlist = g_list, ncol = 2, nrow = 3) %>%
  ggexport(filename = 'figures/manu/raw/other subpop frequency changes.pdf',
           width = 12, height = 10)




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


