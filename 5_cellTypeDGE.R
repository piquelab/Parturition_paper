#############################################################
### cluster marker identification
### 
#############################################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
#library(clusterProfiler)
##library(SingleR)





future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################

#outFolder="./5_harmonyClustersDGE/"
#outFolder="./5_harmonyClusters_SoupX_DGE/"
outFolder="./5_harmonyClusters_withcovid19control_DGE/"
system(paste0("mkdir -p ", outFolder))


#sc <- read_rds("4_harmony_res0.8/sc.NormByLocationRep.Harmony.rds")
#sc <- read_rds("4_harmony_SoupX_res1.0/sc.NormByLocationRep.Harmony.rds")
sc <- read_rds("4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.rds")
#sc <- subset(sc, Condition=="Control")
# new_names <- read_tsv("ClusterAssignment.tsv")
# clust2Names <- new_names$scLabor_ID
# names(clust2Names) <- new_names$seurat_clusters
# cc <- new_names %>% dplyr::select(scLabor_ID,color) %>% unique
# cluster.Colors <- cc$color
# names(cluster.Colors) <- cc$scLabor_ID
# sc@meta.data$cluster_name <- clust2Names[sc@meta.data$seurat_clusters]
# Idents(sc) <- "cluster_name"



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
#anno <- read_rds("3_MergeDemux_Output/")




dim(sc)

table(sc$Library)

table(sc$Location) 


# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


m2 <- markers %>% left_join(select(anno,gene=kbid,symbol=gene_name,rs)) 

m2 <- m2 %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)

fname=paste0(outFolder,"ClusterDEG.tsv");
write_tsv(m2,fname)
m2<-read_tsv(fname)
write.csv(m2,file=paste0(outFolder,"ClusterDEG.csv"))
write_rds(m2,file=paste0(outFolder,"ClusterDEG.rds"))

#m2$celltype<-clust2Names[as.character(m2$cluster)]
#write.csv(m2,file=paste0(outFolder,"ClusterDEG_withcelltypes.csv"))
write.csv(m2,file=paste0("5_harmonyClusters_withcovid19control_DGE/ClusterDEG_withcelltypes.csv"))





top20 <- m2 %>% top_n(n = 20, wt = avg_log2FC)


# fname=paste0(outFolder,"ClusterDEGtop20.tsv");
# write_tsv(top20,fname)
top20<-read_tsv(fname)
write.csv(top20,file=paste0(outFolder,"ClusterDEGtop20.csv"))
write_rds(top20,file=paste0(outFolder,"ClusterDEGtop20.rds"))



top20$celltype<-clust2Names[as.character(top20$cluster)]
write.csv(top20,file=paste0(outFolder,"ClusterDEGtop20_withcelltypes.csv"))


