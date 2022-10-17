#############################################################
### sub-type marker identification
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

outFolder="./5_harmonyClusters_withcovid19control_subtypes_DGE/"
system(paste0("mkdir -p ", outFolder))


sc <- read_rds("4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.rds")

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



 
clusters1<-as.character(c(0,2,8,20,22,27)) #Stromal
clusters2<-as.character(c(4,6,26,29)) #Macrophage
clusters3<-as.character(c(5,9,15,32))#Tcell
clusters4<-as.character(c(10,12,17)) #Decidual
clusters5<-as.character(c(3,13,16,19,30)) #CTB

subtypes_list<-list(clusters1,clusters2,clusters3,clusters4,clusters5)
names(subtypes_list)<-c("Stromal","Macrophage","Tcell","Decidual","CTB")
subtypes<-names(subtypes_list)


 mlist<-lapply(1:length(subtypes),function(x){
   
   print(subtypes[x])
   sc1 <- subset(sc, subset = seurat_clusters %in% subtypes_list[[x]])
   markers <- FindAllMarkers(sc1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
   markers$Celltype<-as.character(clust2Names[as.character(markers$cluster)])
   m2 <- markers %>% left_join( dplyr::select(anno,gene=kbid,symbol=gene_name,rs)) 
   m2 <- m2 %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)
   
   outFolder<-paste0(outFolder,subtypes[x],"/")
   system(paste0("mkdir -p ", outFolder,subtypes[x],"/"))
   
   fname=paste0(outFolder,"ClusterDEG.tsv");
   write_tsv(m2,fname)
   #write.csv(m2,file=paste0(outFolder,"ClusterDEG.csv"))
   
   print(dim(m2))
   top100 <- m2 %>% top_n(n = 100, wt = avg_log2FC)
   
   fname=paste0(outFolder,"ClusterDEGtop100.tsv");
   write_tsv(top100,fname)
   #write.csv(top100,file=paste0(outFolder,"ClusterDEGtop100.csv"))
   m2
   })
 
 m2final<-do.call(rbind,mlist)
 fname=paste0(outFolder,"ClusterDEG.tsv");
 write_tsv(m2final,fname)
 