library(dplyr)
library(ggplot2)
library(ggpubr)

info_path <- file.path('data', '210825_database_for_petter_sent_210826.xlsx')
fileInfo = readxl::read_excel(info_path, sheet=1) %>%
  select(-Sample_ID) %>%
  tidyr::replace_na(list(NEC = FALSE)) %>%
  dplyr::rename(timepoint = PN_days_sample, 
                Subject_ID = Inf_ID) %>%
  mutate(group = as.factor(case_when(
    acttrt == 1 ~ 'control',
    acttrt == 2 ~ 'DHA_AA'
  )),
  Timepoint_string = factor(Timepoint_string, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks')),
  gender = if_else(SEX==1, 'boy', 'girl')
  )


timepoint_summary <- lapply(unique(fileInfo$Timepoint_string), function(x) {
  fileInfo %>% filter(Timepoint_string == x) %>% group_by(timepoint) %>% tally() %>% mutate(Timepoint_string = x)
}) %>% do.call(what=rbind)

lapply(unique(timepoint_summary$Timepoint_string), function(x){
  ggplot(timepoint_summary %>% filter(Timepoint_string == x), aes(x=timepoint, y = n)) +
    geom_bar(stat='identity') +
    ggtitle(x) +
    theme_classic() +
    xlab('PN days')
}) %>%
  ggarrange(plotlist = ., ncol = 2, nrow = 3, align = 'hv') %>%
  ggexport(filename = 'figures/timepoint_distribution_within_group.pdf')

# gender and group
ggplot(fileInfo %>% select(Subject_ID, gender, group) %>% distinct(),
       aes(x=group, fill=gender)) +
  geom_bar(stat='count', position = 'fill') +
  geom_text(aes(label=after_stat(count)), stat='count', position=position_fill(vjust = 0.5)) +
  theme_classic()
ggsave(filename = 'figures/cohort_gender_group.pdf', width = 4, height = 4)

# sampling info
ggplot(fileInfo %>% mutate(Subject_ID = factor(Subject_ID, levels = unique((fileInfo %>% arrange(group))$Subject_ID))), 
       aes(x=timepoint, y=Subject_ID)) +
  lemon::geom_pointline(aes(group=Subject_ID, color=group)) +
  geom_vline(xintercept = c(1.5, 5.5, 12, 20, 45), color = "blue", size=0.5, alpha=0.5) +
  theme_classic() +
  xlab('PN days') +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(filename = 'figures/sampling_info.pdf', width = 8, height = 8)

ggplot(fileInfo %>% select(Subject_ID, Breastmilk_group, preterm_class) %>% distinct(), 
       aes(x=Breastmilk_group, fill=preterm_class)) + 
  geom_bar(stat='count', position = 'stack') +
  geom_text(aes(label=after_stat(count)), stat='count', position=position_stack(vjust = 0.5)) +
  theme_classic()
ggsave(filename = 'figures_milk_group/milk_preterm_distribution.pdf', width = 8, height = 8)

