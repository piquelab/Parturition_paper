

# covid control samples

# these samples were not used because there were no details: "HPL20888" "HPL20874" "HPL20922"
load("/wsu/home/groups/prbgenomics/covid19/covid19analysis_public_repo/3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc_covid<-sc
sc_covid <- subset(sc_covid, subset = Condition=="Control")
sc_covid@meta.data$Location <- "CAM"
sc_covid@meta.data$Location [grepl("PVBP",sc_covid@meta.data$EXP)] <- "PVBP" 

sc_covid<-subset(sc_covid, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)

table(sc_covid$status)
table(sc_covid$DROPLET.TYPE)


load("./3_MergeDemux_Output/scFilteredSeurat.Rdata")

sc@meta.data$Condition<-sc$Labor
sc<-subset(sc, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)
table(sc$status)
table(sc$DROPLET.TYPE)

sc <- merge(sc_covid,sc, project="parturition")




### Merge


#LogNormalize
sc <- NormalizeData(sc, verbose=TRUE) 

#Identifies features that are outliers on a 'mean variability plot'.

#vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

#Scales and centers features in the dataset. If variables are provided in vars.to.regress, they are individually regressed against each feautre, and the resulting residuals are then scaled and centered.
sc <- ScaleData(sc, verbose = TRUE) 

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)


# remove the influence of dataset-of-origin from the embedding. By default, Harmony accepts a normalized gene expression matrix and performs PCA. Since here we already have the PCs, we specify do_pca=FALSE. The matrix harmony_embeddings is the matrix of Harmony corrected PCA embeddings.
sc <- RunHarmony(sc,c("Library"),reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)


###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)


sapply(c(0.8, "1.0", 1.5,2),function(x)
{
  #outFolder=paste0("./4_harmony_res",x,"/")
  outFolder=paste0("./4_harmony_with_covidcontrol_res",x,"/")
  system(paste0("mkdir -p ", outFolder))
  
  sc2 <- FindClusters(sc, verbose = TRUE,resolution=as.numeric(x))
  fname=paste0(outFolder,"sc.NormByLocationRep.Harmony.rds")
  write_rds(sc2,fname)
  
  sc2 <- subset(sc2, subset = Condition!="NA")
  
  #sc2<-read_rds(paste0(outFolder,"sc.NormByLocationRep.Harmony.rds"))
  # Make a simple plot here:
  fname=paste0(outFolder,"UMAP_LocationHarmony.pdf");
  pdf(fname,width=10,height=5)
  aa <- FetchData(sc2,c("UMAP_1","UMAP_2","seurat_clusters","Location","Labor","Condition","SNG.BEST.GUESS","Library"))
  
  p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
    geom_point(size=0.1) +
    facet_grid(Condition ~ Location) +
    theme_bw()
  p1
  ##    theme_black()
  dev.off()
  
  
})




#################

## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.pdf");
pdf(fname,width=12,height=5)
aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","SNG.BEST.GUESS")) 
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
  geom_point(size=0.1) +
  facet_grid(Condition ~ Location) +
  theme_bw()
p1
##    theme_black()
dev.off()


### END- HERE ###
########################################################
md<-sc@meta.data
md<-md[,c("Library","nCount_RNA")]
md  %>% dplyr::count(Library,nCount_RNA )

md %>% 
  group_by(Library) %>% 
  summarise(UMI = sum(nCount_RNA))


sc<-read_rds("./4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.rds")
aa <- FetchData(sc,c("UMAP_1","UMAP_2","Location","Condition","Labor","Origin","status","FetalSex","seurat_clusters","cluster_name","SNG.BEST.GUESS","Status")) 
aa<-aa %>% filter (!is.na(status))
aa<-aa %>% arrange(desc(status))
aa$status <- factor(aa$status,levels=unique(aa$status))

outFolder=paste0("./4_harmony_with_covidcontrol_res1.0/")
system(paste0("mkdir -p ", outFolder))
fname=paste0(outFolder,"UMAP_LocationHarmony.Status.pdf");
pdf(fname,width=7,height=5)
p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=status)) +
  geom_point(size=0.1) +
  ##    scale_color_manual(values=group.colors) +
  guides(colour = guide_legend(override.aes = list(size=10),title="Status")) +
  scale_color_manual("Status",values=c("doublet"="#A61BB5" ,"singlet"="#D1D1D1","unassigned"="blue"))+
  #facet_wrap(~Location) +
  theme_bw()
p2
##    theme_black()
dev.off()