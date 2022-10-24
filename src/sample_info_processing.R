library(dplyr)
library(readr)


info_path <- file.path('data', '210825_database_for_petter_sent_210826.xlsx')
olink <- readr::read_delim('data/Olink_breastmilk_group_72_inds.txt', delim = '\t') %>% rename(Subject_ID = Infant_ID)
group2 <- olink %>% select(Subject_ID, Breastmilk_group) %>% distinct()

fileInfo = readxl::read_excel(info_path, sheet=1) %>%
  select(-Sample_ID) %>%
  dplyr::rename(PN_days = PN_days_sample, 
                Subject_ID = Inf_ID) %>%
  left_join(group2, by= 'Subject_ID') %>%
  tidyr::replace_na(list(NEC = FALSE, Breastmilk_group='Unknown')) %>%
  mutate(group = as.factor(case_when(
    acttrt == 1 ~ 'control',
    acttrt == 2 ~ 'DHA_AA'
  )),
  GA_days_at_birth = GAGEW * 7 + GAGED,
  GA_days = GA_days_at_birth + PN_days,
  Timepoint_string = factor(Timepoint_string, levels = c('Day 1', 'Day 3', 'Day 7', 'Day 14', 'Day 28', 'PMA 40 weeks'))) %>% 
  mutate(PNA_bin = case_when(
    PN_days == 0 ~ '0 days',
    PN_days > 0 & PN_days <= 2 ~ '0-2 days',
    PN_days > 2 & PN_days <= 4 ~ '2-4 days',
    PN_days > 4 & PN_days <= 7 ~ '4-7 days',
    PN_days > 7 & PN_days <= 14 ~ '7-14 days',
    PN_days > 14 & PN_days <= 28 ~ '14-28 days',
    PN_days > 28 & PN_days <= 42 ~ '28-42 days',
    PN_days > 75 & PN_days <= 100 ~ '75-100 days',
    PN_days > 100 & PN_days <= 115 ~ '100-115 days',
    PN_days > 115 ~ '>115 days',
  )) %>% 
  mutate(PNA_bin = factor(PNA_bin, levels = c('0 days', '0-2 days', '2-4 days', '4-7 days', '7-14 days',
                                              '14-28 days', '28-42 days', '75-100 days', '100-115 days',
                                              '>115 days'))) %>%
  mutate(GA_bin = case_when(
    GA_days > 160 & GA_days <= 170 ~ '160-170 days',
    GA_days > 170 & GA_days <= 175 ~ '170-175 days',
    GA_days > 175 & GA_days <= 180 ~ '175-180 days',
    GA_days > 180 & GA_days <= 185 ~ '180-185 days',
    GA_days > 185 & GA_days <= 190 ~ '185-190 days',
    GA_days > 190 & GA_days <= 195 ~ '190-195 days',
    GA_days > 195 & GA_days <= 200 ~ '195-200 days',
    GA_days > 200 & GA_days <= 230 ~ '200-230 days',
    GA_days > 260 & GA_days <= 280 ~ '260-280 days',
    GA_days > 280 & GA_days <= 285 ~ '280-285 days',
    GA_days > 285 & GA_days <= 290 ~ '285-290 days',
    GA_days > 290 ~ '>290 days'
  )) %>%
  mutate(GA_bin = factor(GA_bin, levels = c('160-170 days', '170-175 days', '175-180 days', '180-185 days', '185-190 days',
                                            '190-195 days', '195-200 days', '200-230 days', '260-280 days', '280-285 days',
                                            '285-290 days', '>290 days'))) %>%
  mutate(preterm_class = if_else(GAW < 25, 'Periviability (<25)', 'Extremely preterm'))

write_csv(fileInfo, file = 'data/sample_info_refined.csv')



