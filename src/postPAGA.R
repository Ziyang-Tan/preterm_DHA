library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(tibble)
library(wesanderson)
library(rstatix)

data_dir = './PAGA_result_data_flowsom/'
figure_dir = './figures_flowsom_PAGA/'
sub_names <- sapply(dir(path=data_dir), 
                    function(x){strsplit(x,split = '_')[[1]][2]})
info_path <- file.path('data', '210825_database_for_petter_sent_210826.xlsx')
nec_info <- readr::read_delim('/Users/tan/preterm_DHA/data/NEC_info.csv', delim = ';')
fileInfo = readxl::read_excel(info_path, sheet=1) %>%
  select(-Sample_ID) %>%
  left_join(nec_info, by = 'Inf_ID') %>%
  tidyr::replace_na(list(NEC = FALSE)) %>%
  dplyr::rename(#Sample_ID = ID_unique,
                #timepoint = PN_days_sample, 
                Subject_ID = Inf_ID) %>%
  mutate(group = as.factor(case_when(
    acttrt == 1 ~ 'control',
    acttrt == 2 ~ 'AA:DHA'
  )),
  Timepoint_string = factor(Timepoint_string, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks'))) # 1: control 2: DHA
group_info <- fileInfo %>% select(Subject_ID, group) %>% distinct()


for (sub_name in sub_names){
#for (sub_name in c('NK-cells', 'pDC', 'Tregs')){
  #sub_name <- 'Memory CD8 T'
  file_path = file.path(data_dir, paste0('all_', sub_name, '_sample1000_raw.h5ad'))
  datah5 = readH5AD(file = file_path)
  data = t(assay(datah5)) %>% as_tibble()
  meta = colData(datah5) %>% 
    as_tibble() %>% 
    select(timepoint, leiden, ID_unique) %>%
    left_join(fileInfo, by='ID_unique') %>%
    tidyr::replace_na(list(NEC = FALSE))
  
  # subpop freq change by group
  # data_freq = meta %>%
  #   group_by(timepoint, group, leiden) %>%
  #   summarise(n = n()) %>%
  #   ungroup() %>%
  #   group_by(timepoint, group) %>% # for each timepoint x randtrt, frequency add up to 1
  #   mutate(freq = n / sum(n))
  # subpop freq by subject
  data_freq <- meta %>%
    group_by(timepoint, group, leiden) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(timepoint, group) %>%
    mutate(freq = n/sum(n))
  
  g_list <- lapply(sort(unique(meta$leiden)), function(x){
    #d <- data_freq %>% filter(leiden == x) %>% left_join(group_info, by='Subject_ID')
    d <- data_freq %>% filter(leiden == x) %>% ungroup()
    #------------------ aov not yet finished
    #res.aov <- anova_test(
    #  data = d %>% filter(!timepoint %in% c('PMA 40 weeks')), 
    #  dv = frequency, wid = Subject_ID,
    #  between = group, 
    #  within = timepoint
    #)
    #---------------------
    ggplot(data = d, 
           aes(x = timepoint, y = freq, color = group, group=group)) +
      #geom_line(aes(group=Subject_ID), alpha=0.3) +
      geom_point() +
      geom_line() +
      #geom_smooth(aes(x = as.numeric(timepoint), y = freq), method='loess') +
      labs(title = paste0(sub_name, '_cluster_', x)) +
      theme_minimal()
  })
  ggarrange(plotlist = g_list, ncol = 4, nrow = 3) %>% 
    ggexport(filename = file.path(figure_dir, sub_name, paste0(sub_name, '_cluster_changes.pdf')),
             width = 20,
             height = 10)
  

  ################################
  
  # subpop marker expression
  data_anno <- data %>% add_column(leiden = meta$leiden)
  
  g_list <- lapply(colnames(data), function(x){
    d <- data_anno %>% select(x, leiden) %>% dplyr::rename(marker_expression = x)
    ggplot(d, 
           aes(x = marker_expression, y = leiden, fill = leiden)) +
      geom_density_ridges(scale = 4, rel_min_height=.01) +
      scale_fill_manual(values = wes_palette("Rushmore1", nlevels(meta$leiden), type = "continuous")) +
      theme_ridges() + 
      theme(legend.position = "none") +
      xlim(quantile(d$marker_expression, 0.01),
           quantile(d$marker_expression, 0.99)) +
      labs(title = x)
  })
  ggarrange(plotlist = g_list, ncol = 5, nrow = 9) %>% 
    ggexport(filename = file.path(figure_dir, sub_name, paste0(sub_name, '_cluster_marker_expressions.pdf')),
             width = 20,
             height = 40)
}

