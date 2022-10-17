#############################################################
### cell type classification using SingleR 
### reference: PBMC
#############################################################
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)




#outFolder="./4_harmony_cellClass_PBMC/"
#outFolder="./4_harmony_cellClass_SoupX_PBMC/"
outFolder="./4_harmony_cellClass__with_covidcontrol_PBMC/"
system(paste0("mkdir -p ", outFolder))


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)



##############################################
# loading query data
##############################################



# load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
# sc1 <- sc
# DefaultAssay(sc1) <- "RNA"

# including covid control samples
load("/wsu/home/groups/prbgenomics/covid19/covid19analysis_public_repo/3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc_covid<-sc
sc_covid <- subset(sc_covid, subset = Condition=="Control")
sc_covid@meta.data$Location <- "CAM"
sc_covid@meta.data$Location [grepl("PVBP",sc_covid@meta.data$EXP)] <- "PVBP" 
load("./3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc@meta.data$Condition<-sc$Labor
sc1 <- merge(sc_covid,sc, project="parturition")




# including SoupX
# sc1<-read_rds("1_soupx/seuratObj-after-mt-filtering.2022-03-10.rds")
# sc1 <- subset(sc1, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)
# DefaultAssay(sc1) <- "RNA"



##############################################
# loading the PBMC reference 
##############################################

# reference
# more details here: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
# https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub#bib6
# https://www.sciencedirect.com/science/article/pii/S0031320310005819


# reference <- LoadH5Seurat("Seurat_PBMC_multimodal/pbmc_multimodal.h5seurat")
# save(reference,file="Seurat_PBMC_multimodal/pbmc_multimodal.RData")
load("../../labor_myo/myometrium_analysis/Seurat_PBMC_multimodal/pbmc_multimodal.RData")


DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()



# seurat query example
#library(SeuratData)
# InstallData('pbmc3k')
# data("pbmc3k")
#pbmc3k <- SCTransform(pbmc3k, verbose = TRUE)



##############################################
# Converting ensamble to symbol (qurery)
##############################################

library(tidyverse)
# ## Annotate genes with annotables or other and filtering chrM?
cmd <- paste0("zcat ",
              "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz",
              " | awk '$3~/gene/'",
              " | sed 's/ID=.*gene_id=//;s/;gene_type=/\\t/;s/;gene_name=/\\t/;s/;.*//'")
cat(cmd,"\n")

## Check 0 or 1-based coordinates.

aux <- data.table::fread(cmd=cmd) %>% mutate(TSS=ifelse(V7=="+",V4,V5)) %>%
    select(Chr=V1,Min=V4,Max=V5,kbid=V9,TSS,Strand=V7,Type=V10,gene_name=V11)
##anno <- tibble(kbid=rownames(sc)) %>% mutate(ensgene=gsub("\\..*","",kbid)) %>% left_join(grch38)

anno <- tibble(kbid=rownames(sc1),rs=rowSums(sc1@assays$RNA@data)) %>% filter(rs>0) %>% left_join(aux) %>%
    filter(!is.na(Chr))
# 
# 
sc1 <- sc1[anno$kbid,]

cn<-sc1@assays$RNA@counts
md<-sc1@meta.data

cn<-cn[anno$kbid,]
rownames(cn)<-anno$gene_name


#soupx

#cn<-sc1@assays$RNA@counts

selectedgene <- intersect(rownames(cn),rownames(reference))
cn<-cn[selectedgene,]
newquery <- CreateSeuratObject(counts = cn, meta.data = md)



#sc1<-RenameGenesSeurat(sc1,newnames=anno$gene_name)
sc1 <- SCTransform(newquery, verbose = TRUE)

## Subset to genes in anno
# sc1<-sc1[selectedgene,]
# reference<-reference[selectedgene,]

# sc1@assays$RNA@data<-sc1@assays$RNA@data[selectedgene,]
# reference@assays$SCT@data<-reference@assays$SCT@data[selectedgene,]
#  sc1 <- FindVariableFeatures(sc1, selection.method = "vst")
#  sc1 <- ScaleData(sc1, verbose = TRUE) 
#  sc1 <- RunPCA(sc1,pc.genes = sc1@var.genes, npcs = 100, verbose = TRUE)
# # DefaultAssay(sc1) <- "RNA"
#  sc1 <- RunUMAP(object = sc1, reduction = "pca", dims = 1:30)
#  
#  sc1 <- FindNeighbors(object = sc1, reduction = "pca", dims = 1:30)
#  
#  sc1 <- FindClusters(sc1, resolution=0.8)
 


##############################################
#  finding anchors between reference and query
##############################################


anchors <- FindTransferAnchors(
    reference = reference,
    query = sc1, #pbmc3k, 
    normalization.method = "SCT",
    reference.reduction = "spca",
    recompute.residuals=FALSE,
    dims = 1:50)


##############################################
#transfer cell type labels and protein data from the reference to the query
##############################################
# https://satijalab.org/seurat/reference/mapquery
# MapQuery() is a wrapper around three functions: TransferData(), IntegrateEmbeddings(), and ProjectUMAP(). 
# TransferData() is used to transfer cell type labels and impute the ADT values. IntegrateEmbeddings() and ProjectUMAP() are used to project the query data onto the UMAP structure of the reference. The equivalent code for doing this with the intermediate functions is below:
    
    
# Transfer data
#TransferData()
# https://satijalab.org/seurat/reference/transferdata 
# The main difference between label transfer (classification) and feature imputation is what gets multiplied by the weights matrix. For label transfer, we perform the following steps: 
#     Create a binary classification matrix, the rows corresponding to each possible class and the columns corresponding to the anchors. If the reference cell in the anchor pair is a member of a certain class, that matrix entry is filled with a 1, otherwise 0.


# IntegrateEmbeddings()
# https://satijalab.org/seurat/reference/integrateembeddings

# ProjectUMAP()
# https://satijalab.org/seurat/reference/projectumap

# This function will take a query dataset and project it into the coordinates of a provided reference UMAP. This is essentially a wrapper around two steps:
#     
#     FindNeighbors - Find the nearest reference cell neighbors and their distances for each query cell.
# 
# RunUMAP - Perform umap projection by providing the neighbor set calculated above and the umap model previously computed in the reference.

sc1 <- MapQuery(
    anchorset = anchors,
    query = sc1,
    reference = reference,
    refdata = list(
        celltype.l1 = "celltype.l1",
        celltype.l2 = "celltype.l2",
        predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
)

#Explore the mapping results

fname=paste0(outFolder,"UMAP_predicted.celltypes.l1.l2.png");
png(fname,width=1000,height=600)
p1 = DimPlot(sc1, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(sc1, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2
dev.off()


# predicted.celltype.l1<-sc1$predicted.celltype.l1
# predicted.celltype.l2<-sc1$predicted.celltype.l2


predicted.celltype<-sc1@meta.data %>% select(predicted.celltype.l1,predicted.celltype.l2,predicted.celltype.l1.score,predicted.celltype.l2.score)
    
    
fname=paste0(outFolder,"sc.NormByLocation.ref.Anchors.rds")
write_rds(predicted.celltype,fname)

#ref
 md <- predicted.celltype %>% as.data.frame() %>%
    rownames_to_column("BARCODES") %>%
     left_join(sc1@meta.data %>% rownames_to_column("BARCODES"))

 fname=paste0(outFolder,"sc.NormByLocation.ref.Anchors.csv")
 write_csv(md,fname)


fname=paste0(outFolder,"sc.Integrated.Harmony.rds")
write_rds(sc1,fname)



