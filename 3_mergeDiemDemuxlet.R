## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

###########################################
sc <- read_rds("./2_kb_diem_Output/kb_diem_Seurat.obj.rds")

sc@meta.data$Library <- (gsub("_[ACGT]{6,}","",colnames(sc)))

table(sc$Library)

expNames <- unique(sc@meta.data$Library)
expNames


outFolder= "./3_MergeDemux_Output/"
system(paste0("mkdir -p ",outFolder))

future::plan(strategy = 'multicore', workers = 12)
options(future.globals.maxSize = 20 * 1024 ^ 3)


## Demuxlet output. 
dd <- read_rds("./1_demux_output/1_demux_New.ALL.rds") %>% filter(DIFF.LLK.BEST.NEXT>1) %>%
    select(BARCODE=NEW_BARCODE,SNG.BEST.GUESS,DROPLET.TYPE,NUM.SNPS,NUM.READS,DIFF.LLK.BEST.NEXT,EXP)

table(dd$DROPLET.TYPE)
table(dd$EXP)

md = sc@meta.data %>% rownames_to_column("BARCODE") %>% left_join(dd) %>%
    column_to_rownames("BARCODE")

stopifnot(identical(rownames(md),rownames(sc@meta.data)))

## Repeat with soup or cell
dd <- read_rds("./1_souporcell_output/1_souporcell.ALL.rds")

md = md %>% rownames_to_column("barcode") %>% left_join(dd) %>%
    column_to_rownames("barcode")

md$bc <- rownames(md)

md$Pregnancy_ID <- gsub("-.*","",md$SNG.BEST.GUESS)

# Todo: merge with other tables with covariates.
cc <- read.csv("./parturition_cv.csv")

md <- md %>% left_join(cc)
head(md)


# begin (azam)
libLoc <- read.csv("./LibLocation.csv")
md <- left_join(md,libLoc)
# 
# cv <- read.csv("./parturition_cv.csv")
# md <- left_join(md,cv)
# end (azam)

## Have a file with all the individuals and the correct libraries they should be in. 
cc <- read_csv("./ParturitionLibraries.csv")
tail(cc)          

okbatches <- paste0("",cc$Pregnancy_ID,"_",cc$EXP)
md$okbatch <- paste0(md$Pregnancy_ID,"_",md$EXP) %in% okbatches 
mean(md$okbatch)

md2 <- md

md <- md2

## Filter cells in the right batch.
md <- filter(md,okbatch)
mean(md$okbatch)



table(md$status,md$DROPLET.TYPE)

table(md$assig2==md$SNG.BEST.GUESS)

table(md$Library,is.na(md$SNG.BEST.GUESS))

table(md$Library,md$SNG.BEST.GUESS)


sc <- sc[,md$bc]

sc@meta.data <- md %>% column_to_rownames("bc")

## Annotate genes with annotables or other and filtering chrM?
cmd <- paste0("zcat ",
              "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz",
              " | awk '$3~/gene/'",
              " | sed 's/ID=.*gene_id=//;s/;gene_type=/\\t/;s/;gene_name=/\\t/;s/;.*//'")
cat(cmd,"\n")

## Check 0 or 1-based coordinates. 

aux <- data.table::fread(cmd=cmd) %>% mutate(TSS=ifelse(V7=="+",V4,V5)) %>%
    select(Chr=V1,Min=V4,Max=V5,kbid=V9,TSS,Strand=V7,Type=V10,gene_name=V11) 
##anno <- tibble(kbid=rownames(sc)) %>% mutate(ensgene=gsub("\\..*","",kbid)) %>% left_join(grch38)

anno <- tibble(kbid=rownames(sc),rs=rowSums(sc@assays$RNA@data)) %>% filter(rs>0) %>% left_join(aux) %>%
    filter(!is.na(Chr))

table(is.na(anno$Chr))

table(anno$Chr)

table(anno$Type)

head(anno)

head(aux)


## Subset to genes in anno 
sc <- sc[anno$kbid,]

## Calculate features 
sc[["percent.mt"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrM",]$kbid)

sc[["percent.mt"]] %>% summary()

## Filter sc for things matching genotype. chrM or number of RNAs.  

sc[["percent.Y"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrY",]$kbid)
sc[["percent.Y"]] %>% summary()

##anno[anno$gene_name=="XIST",]$kbid 
sc[["percent.XIST"]] <- PercentageFeatureSet(sc,features=anno[anno$gene_name=="XIST",]$kbid )
sc[["percent.XIST"]] %>% summary()

table(sc[["percent.XIST"]]>0.01,md$SNG.BEST.GUESS)


## Save merged object without any filter.
fname=paste0(outFolder,"scFullSeurat.Rdata")
cat(fname,"\n")
save(sc,anno,file=fname)

##    ##  ## 
table(sc[["percent.mt"]]<25,sc@meta.data$nFeature_RNA > 100) 

table(sc@meta.data$Library,sc[["percent.mt"]]<25 & sc@meta.data$DIFF.LLK.BEST.NEXT > 3 & sc@meta.data$nFeature_RNA > 100 & sc$DROPLET.TYPE=="SNG")

table(sc@meta.data$SNG.BEST.GUESS,sc[["percent.mt"]]<25 & sc@meta.data$DIFF.LLK.BEST.NEXT > 3 & sc@meta.data$nFeature_RNA > 100 & sc$DROPLET.TYPE=="SNG") 

table(sc$status,sc$nFeature_RNA > 10000)

sc <- subset(sc, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)

table(sc@meta.data$Library, sc@meta.data$SNG.BEST.GUESS) 

## Save merged object without any filter.
fname=paste0(outFolder,"scFilteredSeurat.Rdata")
cat(fname,"\n")
save(sc,anno,file=fname)

