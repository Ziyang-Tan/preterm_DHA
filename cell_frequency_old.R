library(dplyr)
library(tidyr)
library(mmR)
library(ggplot2)
library(ggpubr)


sample_info_path = file.path('data', '210825_database_for_petter_sent_210826.xlsx')
data_path = Sys.glob('data/*/abundance.csv')


resolve_frequency_table <- function(data){
  tmp1 <- data %>% 
    group_by(level1) %>% 
    summarise(across(!any_of(c('level1', 'level2', 'level3')), sum)) %>%
    rename(cell_type = level1)
  tmp2 <- data %>% 
    group_by(level2) %>%
    summarise(across(!any_of(c('level1', 'level2', 'level3')), sum)) %>%
    filter(level2 != '') %>%
    rename(cell_type = level2)
  tmp3 <- data %>%
    group_by(level3) %>%
    summarise(across(!any_of(c('level1', 'level2', 'level3')), sum)) %>%
    filter(level3 != '') %>%
    rename(cell_type = level3)
  
  return(rbind(tmp1, tmp2, tmp3))
}


sample_info = mm.fastread(sample_info_path) %>%
  rename(sample_name = ID_unique)
data_list = lapply(data_path, function(x){
  mm.fastread(x) %>% 
    rename(level1 = V1, level2 = V2)
})
data = Reduce(
  function(x, y) merge(x, y, all = T),
  data_list
) %>% mutate_all(~replace(., is.na(.), 0)) %>%
  separate(level2, c('level2', 'level3'), sep = '/')

data_mod = resolve_frequency_table(data)

data_long = data_mod %>% pivot_longer(!cell_type, 
                                      names_to = 'sample_name', 
                                      values_to = 'frequency') %>%
  filter(sample_name != 'Nil') %>% 
  left_join(sample_info, by = 'sample_name') %>%
  mutate(Timepoint_string = factor(Timepoint_string, 
                                   levels = c('Day 1', 'Day 3', 'Day 7', 
                                              'Day 14', 'Day 28', 'PMA 40 weeks'))) %>%
  mutate(randtrt = factor(randtrt)) %>%
  mutate(log_PN_days = log10(PN_days_sample))


g_list = lapply(data_mod$cell_type, function(x){
  ggplot(data_long %>% filter(cell_type == x), 
         aes(x = log_PN_days, y = frequency, group = Inf_ID, color = randtrt)) +
    geom_line(alpha=0.5) + 
    geom_point(size=0.5, alpha = 0.3) +
    geom_smooth(aes(group = randtrt, color = randtrt), method = 'loess') +
    labs(title=x)
})
ggarrange(plotlist = g_list, ncol = 4, nrow = 15) %>% 
  ggexport(filename = file.path(getwd(), 'figures', 'cell_frequency_log_PN_days.pdf'),
           width = 18,
           height = 50)

