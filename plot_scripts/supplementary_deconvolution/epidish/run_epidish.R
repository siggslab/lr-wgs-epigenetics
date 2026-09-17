library(EpiDISH)
library(readr)
library(tidyverse)

d = read_csv(joined.csv)
d = d %>% column_to_rownames('ilmn')
data(centDHSbloodDMC.m)
outputs = epidish(beta.m=d, ref.m=centDHSbloodDMC.m)

out_ests = as.data.frame(outputs$estF) %>% rownames_to_column('sample')
write_csv(out_ests, 'output.csv')

