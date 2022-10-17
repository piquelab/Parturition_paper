##############################################
### cell type classification using SingleR ###
### reference: parturition elife           ### 
##############################################

## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)

#################
library(SingleR)



future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################
## Testing sc transform           
## 2and3_Diem_Output

# no soupx

#load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
#sc1 <- sc


# with control samples from covid19
load("/wsu/home/groups/prbgenomics/covid19/covid19analysis_public_repo/3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc_covid<-sc
sc_covid <- subset(sc_covid, subset = Condition=="Control")
sc_covid@meta.data$Location <- "CAM"
sc_covid@meta.data$Location [grepl("PVBP",sc_covid@meta.data$EXP)] <- "PVBP" 
load("./3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc@meta.data$Condition<-sc$Labor
sc1 <- merge(sc_covid,sc, project="parturition")



# with soupx
#sc1<-read_rds("1_soupx/seuratObj-newfilter-merge.2022-03-10.rds")

sc1 <- subset(sc1, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)



sc3 <- read_rds("/nfs/rprdata/scilab/novogene/Analyses/Roger_20200218/3_scTransferLabel_scLabor/ST_Integrated.scLabor.obj.rds")

# Subset a Seurat object
sc3 <- subset(sc3, subset = nFeature_RNA > 100)

# to identify the cell types of sc1, another study with known celltypes will be merged to this study and the cell types will be identified

sc <- merge(sc1,list(sc3))

dim(sc)

table(sc$Library)

## table(sc$Location) 

##table(sc$sclabor.tlabel)

##table(sc$Location,sc$sclabor.tlabel)

## Harmony

DefaultAssay(sc) <- "RNA"

sc <- NormalizeData(sc, verbose=TRUE) 

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE) 

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

##sc <- RunHarmony(sc,c("Location","percent.mt","Rep"),reduction="pca")
sc <- RunHarmony(sc,c("Library"),reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)
sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)


#sapply(c(0.5,0.8,1.5),function(x)
#sapply(c(0.6,"1.0"),function(x)
#sapply(c(0.8,"1.0",1.5),function(x)
sapply(c( 0.8, "1.0", 1.5,2),function(x)
{

  #outFolder=paste0("./4_harmony_cellClass_elife",x,"/")
  
  #outFolder=paste0("./4_harmony_cellClass_SoupX_elife",x,"/")
  
  outFolder=paste0("./4_harmony_cellClass_with_covidcontrol_elife",x,"/")
  system(paste0("mkdir -p ", outFolder))


###### Cluster



sc <- FindClusters(sc, verbose = TRUE,resolution=as.numeric(x))


fname=paste0(outFolder,"sc.NormIntegrated.ref.Harmony.rds")
write_rds(sc,fname)

################

he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

#unknown cell types
query.he <- he[,is.na(sc@meta.data$FinalName)]

#known cell types
ref.he <- he[,!is.na(sc@meta.data$FinalName)]

ref.labels <- sc@meta.data$FinalName[!is.na(sc@meta.data$FinalName)]

pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.rds")
write_rds(pred.labels,fname)

md <- pred.labels %>% as.data.frame() %>% 
  rownames_to_column("BARCODES") %>%
  left_join(sc@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.csv")
write_csv(md,fname)
})
## save object.





