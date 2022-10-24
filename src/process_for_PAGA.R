library(dplyr)



load(file = 'data/flowsom_results.RData')
sample_info = sample_info = readxl::read_excel('data/210825_database_for_petter_sent_210826.xlsx') 

mmR::mm.fastwrite(all_bind, path = 'data/flowSOM_all_data.csv.gz', compress = TRUE)
