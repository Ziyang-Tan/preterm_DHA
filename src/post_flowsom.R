library(tidyverse)
library(ggpubr)
library(rstatix)


all_data <- mmR::mm.fastread('data/flowSOM_all_data.csv.gz')
t_data <- mmR::mm.fastread('data/all_T_flowsom.csv.gz')
sample_info = readxl::read_excel('data/210825_database_for_petter_sent_210826.xlsx') %>%
  mutate(group = case_when(
    acttrt == 1 ~ 'control',
    acttrt == 2 ~ 'DHA_AA'
  ),
  Timepoint_string = as.factor(Timepoint_string))

data <- rbind(all_data %>% 
                filter(lineage != 'T cell'),
              t_data %>%
                mutate(T_cell_subtype_1 = case_when(
                  T_cell_subtype_2 == 'CD4Treg' ~ 'CD4Treg',
                  TRUE ~ T_cell_subtype_1
                )) %>% 
                mutate(subtype = T_cell_subtype_1) %>%
                select(-T_cell_subtype_1, -T_cell_subtype_2, -T_cluster)
)
rm(all_data, t_data)
label_info <- data %>% select(ID_unique, lineage, subtype)

# calculate frequency  -----

freq <- label_info %>% 
  group_by(ID_unique, subtype) %>% 
  summarise(count = n()) %>% 
  group_by(ID_unique) %>% 
  mutate(freq = count/sum(count)) %>%
  left_join(sample_info, by='ID_unique')

# mixed ANOVA (two-way) -----

# 
for(sub_name in unique(freq$subtype)){
  print(sub_name)
  d <- freq %>% filter(subtype == !!sub_name) %>% ungroup()
  res.aov <- anova_test(
    data = d, dv = freq, wid = Inf_ID,
    between = group, within = Timepoint_string
  )
  print(get_anova_table(res.aov))
}

d<- Treg_freq %>% left_join(sample_info, by='ID_unique')
res.aov <- anova_test(
  data = d, dv = Treg_freq, wid = Inf_ID,
  between = group, within = Timepoint_string
)
print(get_anova_table(res.aov))
ggboxplot(
  d, x = "Timepoint_string", y = "Treg_freq",
  color = "group", palette = "jco"
)










