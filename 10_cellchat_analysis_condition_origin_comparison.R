library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)



#######################################################################
# Location- specific analysis
#######################################################################




cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)

outFolder="./10_CellChat_analysis_withcovidcontrol_origin_or_500/"

system(paste0("mkdir -p ", outFolder))


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

# CellChat requires two user inputs: 
# one is the gene expression data of cells, 
# and the other is either user assigned cell labels (i.e., label-based mode) 
# or a low-dimensional representation of the single-cell data

####################################
# Load data
####################################

# load("./3_MergeDemux_Output/scFilteredSeurat.Rdata")
# sc2<-sc


load("/wsu/home/groups/prbgenomics/covid19/covid19analysis_public_repo/3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc_covid<-sc
sc_covid <- subset(sc_covid, subset = Condition=="Control")
sc_covid@meta.data$Location <- "CAM"
sc_covid@meta.data$Location [grepl("PVBP",sc_covid@meta.data$EXP)] <- "PVBP" 

sc_covid<-subset(sc_covid, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)


md <- sc_covid@meta.data
md<-md %>% separate(SNG.BEST.GUESS,into=c("Pcase", "Origin",sep="-",remove=FALSE))
md<-md %>% filter(!is.na(Condition))
md<-md %>% filter (!Pregnancy_ID %in% c("HPL20888", "HPL20874", "HPL20922"))
md<-md[,which(!colnames(md)%in% c("Condition"))]
parturion_samples<-read_delim("parturition_cv.txt")
parturion_samples_data<- parturion_samples %>% select(Pregnancy_ID,Condition=Labor)
md<-md %>% rownames_to_column("barcode") %>% inner_join(parturion_samples_data)
rownames(md)<-md$barcode
sc_covid@meta.data<-md
sc_covid$Condition<-md$Condition

#table(sc_covid$status)
#table(sc_covid$DROPLET.TYPE)


load("./3_MergeDemux_Output/scFilteredSeurat.Rdata")

sc@meta.data$Condition<-sc$Labor
sc<-subset(sc, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)
table(sc$status)
table(sc$DROPLET.TYPE)
sc$Origin<-unlist(str_split(sc$SNG.BEST.GUESS,"-"))[seq(2,2*length(sc$SNG.BEST.GUESS),by=2)]
sc2 <- merge(sc_covid,sc, project="parturition")


cmd <- paste0("zcat ",
              "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz",
              " | awk '$3~/gene/'",
              " | sed 's/ID=.*gene_id=//;s/;gene_type=/\\t/;s/;gene_name=/\\t/;s/;.*//'")
cat(cmd,"\n")

## Check 0 or 1-based coordinates. 

aux <- data.table::fread(cmd=cmd) %>% mutate(TSS=ifelse(V7=="+",V4,V5)) %>%
  select(Chr=V1,Min=V4,Max=V5,kbid=V9,TSS,Strand=V7,Type=V10,gene_name=V11) 
##anno <- tibble(kbid=rownames(sc)) %>% mutate(ensgene=gsub("\\..*","",kbid)) %>% left_join(grch38)

anno <- tibble(kbid=rownames(sc2),rs=rowSums(sc2@assays$RNA@data)) %>% filter(rs>0) %>% left_join(aux) %>%
  filter(!is.na(Chr))



sc2 <- sc2[anno$kbid,]
sc2 <- NormalizeData(sc2, verbose=TRUE)



sc <- read_rds("4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.rds")
metadata <- sc@meta.data
metadata<-metadata %>% separate(SNG.BEST.GUESS,into=c("Pcase", "Origin2",sep="-",remove=FALSE))
metadata$Origin[is.na(metadata$Origin)]<-metadata$Origin2[is.na(metadata$Origin)]
metadata<-metadata %>% filter (!Pregnancy_ID %in% c("HPL20888", "HPL20874", "HPL20922"))
metadata<-metadata[,which(colnames(metadata)!="Condition")]
parturion_samples_data<- parturion_samples %>% select(Pregnancy_ID,Condition=Labor)
metadata<-metadata %>% rownames_to_column("barcode") %>% inner_join(parturion_samples_data)
rownames(metadata)<-metadata$barcode
metadatacopy<-metadata


locations<- c("CAM"  ,     "PVBP")
conditions<-c("TIL" ,"TNL")
mclapply(locations, function(xlocation)
  {
  
  for (xcondition in conditions)
  {
    
    
    metadata<-metadatacopy
    metadata<-metadata %>% filter(Origin!="NA")
    
    sc2_filter<-sc2
    nextcondition<-conditions[which(conditions!=xcondition)]
    
    sc2_filter2<-subset(sc2,Location==xlocation & Labor==nextcondition) 
    mapping<-clust2Names[metadata$seurat_clusters]
    names(mapping)<-rownames(metadata)
    sc2_filter2$barcode<-colnames(sc2_filter2)
    sc2_filter2$seurat_clusters<-mapping[sc2_filter2$barcode]
    sc2_filter2$celltype<-paste0(sc2_filter2$seurat_clusters,"_",sc2_filter2@meta.data$Origin)
    cluster_filter2<-names(which(table(sc2_filter2$celltype)>500))
    
    #subset(sc2,Location==xlocation & Labor==xcondition) 
    
    sc2_filter<-subset(sc2,Location==xlocation & Labor==xcondition) 
    mapping<-clust2Names[metadata$seurat_clusters]
    names(mapping)<-rownames(metadata)
    sc2_filter$barcode<-colnames(sc2_filter)
    sc2_filter$seurat_clusters<-mapping[sc2_filter$barcode]
    sc2_filter$celltype<-paste0(sc2_filter$seurat_clusters,"_",sc2_filter@meta.data$Origin)
    cluster_filter<-names(which(table(sc2_filter$celltype)>500))
    
    additional<-cluster_filter2 [ which(cluster_filter2 %in% names(table(sc2_filter$celltype)))]
    cluster_filter<-unique(c(cluster_filter,additional))
    metadata<-metadata %>% filter (Location==xlocation & Condition==xcondition) 
    #metadata_cluster_count <- metadata %>% dplyr::count(seurat_clusters)%>% filter(n>100)
    
    metadata$seurat_clusters<-paste0(clust2Names[metadata$seurat_clusters] ,"_",metadata$Origin)
    metadata<-metadata %>% filter( seurat_clusters %in% cluster_filter )
    sc2_filter<-subset(sc2_filter,barcode %in% rownames(metadata) & celltype %in% cluster_filter) 
    
    metadata$labels<-metadata$seurat_clusters
    data<-sc2_filter@assays$RNA@data
    
    convertensymbol<-anno$gene_name
    names(convertensymbol)<-anno$kbid
    rownames(data)<-convertensymbol[rownames(data)]
    data<-data[,colnames(data)%in% rownames(metadata)]
    
    rw<-intersect(colnames(data),rownames(metadata))
    data<-data[,rw]
    metadata<-metadata[rw,]
    
    cellchat <- createCellChat(object = data, meta = metadata, group.by = "labels")
    cellchat <- addMeta(cellchat, meta = metadata)
    cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
    levels(cellchat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
    
    
    ####################################
    # Set the ligand-receptor interaction database
    ####################################
    
    CellChatDB <- CellChatDB.human
    
    # pdf(paste0(outFolder,"showDatabaseCategory.pdf"))
    # showDatabaseCategory(CellChatDB)
    # dev.off()
    
    # Show the structure of the database
    dplyr::glimpse(CellChatDB$interaction)
    
    
    # use a subset of CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
    # use all CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    
    # simply use the default CellChatDB
    
    CellChatDB.use <- CellChatDB 
    
    
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
    
    ########################################################################
    # Preprocessing the expression data for cell-cell communication analysis
    ########################################################################
    
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    # project gene expression data onto PPI network (optional)
    cellchat <- projectData(cellchat, PPI.human)
    
    
    # Part II: Inference of cell-cell communication network
    
    # Compute the communication probability and infer cellular communication network
    cellchat <- computeCommunProb(cellchat)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    
    
    ########################################################################
    # Infer the cell-cell communication at a signaling pathway level
    ########################################################################
    
    # NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
    cellchat <- computeCommunProbPathway(cellchat)
    
    ########################################################################
    # Calculate the aggregated cell-cell communication network
    ########################################################################
    
    cellchat <- aggregateNet(cellchat)
    
    write_rds(cellchat,file=paste0(outFolder,"cellchat_",xlocation,"_",xcondition,"_",Sys.Date(),".rds"))
    
  }
  
},mc.cores=6)









