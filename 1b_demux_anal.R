#
library(tidyverse)
library(parallel)
##library(data.table)

### 
outFolder="./1_demux_output/"

opfn <- "./1_demux_output/1_demux_New.SNG.rds"
demux <- read_rds(opfn)

# filter 
aa <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10,NUM.READS<20000) %>%
  select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,Sample_ID=SNG.BEST.GUESS)%>%
    separate(Sample_ID,c('Pregnancy_ID','Origin'),'-',remove=FALSE)



libLoc <- read.csv("./LibLocation.csv")
aa <- left_join(aa,libLoc)

cv <- read.csv("./parturition_cv.csv")
aa <- left_join(aa,cv)

aa$Condition <- aa$Labor

cell.counts <- aa %>% group_by(EXP,Location,Sample_ID,Pregnancy_ID,Origin,Condition) %>%
  summarize(n=n()) 

write_tsv(cell.counts,paste0(outFolder,"cell.counts.tsv"))

cc.wider <- cell.counts %>% ungroup() %>% select(EXP,Location,Pregnancy_ID,Origin,Condition,n,Sample_ID) %>%
    filter(!is.na(Origin)) %>% 
  group_by(EXP,Pregnancy_ID,Origin) %>%  
  pivot_wider(names_from=Origin,values_from=c(n,Sample_ID),values_fill=list(n=0,Sample_ID="UNK")) %>%
  mutate(n_T=n_M+n_F+n_NA) %>% 
  filter(n_T>200)  %>%
  arrange(EXP,-n_T)

write_tsv(cc.wider,paste0(outFolder,"cc.wider.tsv"))


