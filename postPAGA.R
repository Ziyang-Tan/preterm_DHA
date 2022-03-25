library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(tibble)
library(wesanderson)

data_dir = './PAGA_result_data/'
sub_names <- sapply(dir(path=data_dir), 
                    function(x){strsplit(x,split = '_')[[1]][2]})
nec_info <- readr::read_delim('/Users/tan/preterm_DHA/data/NEC_info.csv', delim = ';') %>%
  mutate(Inf_ID = as.character(Inf_ID))
for (sub_name in sub_names){
  file_path = file.path('./PAGA_result_data', 
                        paste0('all_', sub_name, '_sample1000_raw.h5ad'))
  datah5 = readH5AD(file = file_path)
  data = t(assay(datah5)) %>% as_tibble()
  meta = colData(datah5) %>% 
    as_tibble() %>% 
    select(timepoint, Inf_ID, group, leiden) %>%
    rename(randtrt = group) %>%
    left_join(nec_info, by='Inf_ID') %>%
    tidyr::replace_na(list(NEC = FALSE))
  
  # subpop freq change by randtrt
  data_freq = meta %>% 
    group_by(timepoint, randtrt, leiden) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(timepoint, randtrt) %>% # for each timepoint x randtrt, frequency add up to 1
    mutate(freq = n / sum(n))
  
  g_list <- lapply(sort(unique(meta$leiden)), function(x){
    ggplot(data = data_freq %>% filter(leiden == x), 
           aes(x = timepoint, y = freq, color = randtrt, group = randtrt)) +
      geom_line() +
      geom_point() +
      labs(title = paste0(sub_name, '_cluster_', x))
  })
  ggarrange(plotlist = g_list, ncol = 4, nrow = 3) %>% 
    ggexport(filename = file.path('./figures', sub_name, paste0(sub_name, '_cluster_changes_timepoint_trt.pdf')),
             width = 20,
             height = 10)
  
  # subpop freq change by nec situation
  data_freq = meta %>% 
    group_by(timepoint, NEC, leiden) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(timepoint, NEC) %>%
    mutate(freq = n / sum(n))
  
  g_list <- lapply(sort(unique(meta$leiden)), function(x){
    ggplot(data = data_freq %>% filter(leiden == x), 
           aes(x = timepoint, y = freq, color = NEC, group = NEC)) +
      geom_line() +
      geom_point() +
      labs(title = paste0(sub_name, '_cluster_', x))
  })
  ggarrange(plotlist = g_list, ncol = 4, nrow = 3) %>% 
    ggexport(filename = file.path('./figures', sub_name, paste0(sub_name, '_cluster_changes_timepoint_nec.pdf')),
             width = 20,
             height = 10)
  
  # subpop freq change by nec x trt
  data_freq = meta %>% 
    group_by(timepoint, NEC, randtrt, leiden) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(timepoint, NEC, randtrt) %>%
    mutate(freq = n / sum(n))
  
  g_list <- lapply(sort(unique(meta$leiden)), function(x){
    ggplot(data = data_freq %>% filter(leiden == x), 
           aes(x = timepoint, y = freq, color = NEC, linetype = randtrt,
               group = interaction(NEC, randtrt))) +
      geom_line() +
      geom_point() +
      labs(title = paste0(sub_name, '_cluster_', x))
  })
  ggarrange(plotlist = g_list, ncol = 4, nrow = 3) %>% 
    ggexport(filename = file.path('./figures', sub_name, paste0(sub_name, '_cluster_changes_timepoint_nec_n_trt.pdf')),
             width = 20,
             height = 10)
  ################################
  
  # subpop marker expression
  data_anno <- data %>% add_column(leiden = meta$leiden)
  
  g_list <- lapply(colnames(data), function(x){
    d <- data_anno %>% select(x, leiden) %>% rename(marker_expression = x)
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
    ggexport(filename = file.path('./figures', sub_name, paste0(sub_name, '_cluster_marker_expressions.pdf')),
             width = 20,
             height = 40)
}

