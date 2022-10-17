
#############################################################
### UMAP plots  
 
#############################################################


library(Seurat)
library(Matrix)
library(tidyverse)
library(dplyr)
library(future)
library(harmony)
##library(SingleR)



###########################################
future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)



outFolder="./5_harmony_cellClass_resolutions0.6_plot/"
system(paste0("mkdir -p ", outFolder))


cell.type.annotation<-read.delim("../../labor_myo/myometrium_analysis/data.colors.mouse.human.txt")
cluster.Colors.rest<-cell.type.annotation$cluster.Colors_human2
names(cluster.Colors.rest)<-cell.type.annotation$cluster_human
cluster.Colors.rest<-cluster.Colors.rest[cluster.Colors.rest!=""]




future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")


cluster.Colors.rest<-cluster.Colors.rest[which(cluster.Colors.rest %in% cluster.Colors)]

cluster.Colors<-c(cluster.Colors,cluster.Colors.rest)
names(cluster.Colors)<-as.character(c(0:43))

outFolder=paste0("./5_harmony_cellClass_thresholds0.6/")
system(paste0("mkdir -p ", outFolder))    
sc <- read_rds(paste0("4_harmony_res0.8/sc.NormByLocationRep.Harmony.rds"))


md <- read_rds("./4_harmony_cellClass_thresholds0.6/sc.NormByLibrary.ref.Harmony.singler.thresholds.rds") %>%
    as.data.frame %>%
    rownames_to_column("BARCODES") %>%
    dplyr::select(BARCODES,pruned.labels,scLabor_ID=pruned.labels,scLabor_Score=tuning.scores.second)

# md <- read_rds("./4_harmony_cellClass_PBMC/sc.NormByLocation.ref.Anchors.rds") %>%
#     as.data.frame %>%
#     rownames_to_column("BARCODES") %>%
#     dplyr::select(BARCODES,scLabor_ID=predicted.celltype.l1,scLabor_Score=predicted.celltype.l1.score)

md$BARCODES<-sapply(md$BARCODES, function(x){return(unlist(strsplit(x,"_1"))[1])})

md <- sc@meta.data %>% rownames_to_column("BARCODES") %>%
    left_join(md) 

identical(md$BARCODES,rownames(sc@meta.data))





##setwd(outFolder)

dim(sc)

table(sc$Library)

table(sc$Location) 

table(md$scLabor_ID, md$Location)


cluster.Colors2<-cluster.Colors[unique(sc$seurat_clusters)]

table(md$scLabor_ID, md$Location)

table(md$SNG.BEST.GUESS, md$scLabor_ID)


cc <- md %>% dplyr::select(seurat_clusters,scLabor_ID) %>%
    group_by(seurat_clusters,scLabor_ID) %>%
    summarize(n=n()) %>%
    group_by(seurat_clusters) %>%
    mutate(prop=n/sum(n)) %>%
    arrange(-n)
# 
cc2 <- cc %>% top_n(n=1) %>% mutate(cname=paste0(seurat_clusters,":",scLabor_ID))
clust2name <- cc2$cname
names(clust2name)=cc2$seurat_clusters
clust2name


tt <- table(md$seurat_clusters,md$scLabor_ID)
tt <- tt/rowSums(tt)

colnames(tt)<-sapply(colnames(tt), function(x){return(unlist(strsplit(x,"res0.6_"))[2])})
library(pheatmap)

# wd<-8
# if (length(unique(sc$seurat_clusters))>=12)
#     wd<-12

tt<-tt[,order(as.numeric(colnames(tt)))]
fname=paste0(outFolder,"Singler.HeatMap.pdf");
pdf(fname,width=6,height=7)
#pheatmap(t(tt),cluster_rows=TRUE,cluster_cols=TRUE,scale="none")
pheatmap(t(tt),cluster_rows=FALSE,cluster_cols=FALSE,scale="none")
dev.off()

fname=paste0(outFolder,"UMAP_Harmony.png");
png(fname,width=1000,height=1000)
DimPlot(sc, reduction = "umap", cols = cluster.Colors2, label = TRUE, pt.size = 0.5,label.size = 8) #+ NoLegend()
dev.off()

aa <- FetchData(sc,c("UMAP_1","UMAP_2","Location","Labor","Origin","status","FetalSex","seurat_clusters","cluster_name","SNG.BEST.GUESS")) 
head(aa)

aa<-aa %>% separate(SNG.BEST.GUESS,into=c("Pcase", "Origin",sep="-",remove=FALSE))


fname=paste0(outFolder,"UMAP_LocationHarmony.Origin.pdf");
pdf(fname,width=14,height=5)
p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Origin)) +
    geom_point(size=0.1) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=10),title="Cell origin")) +
    scale_color_manual("Origin",values=c("M"="#D1D1D1","F"="#A61BB5"))+
    facet_wrap(~Location) +
    theme_bw()
p2
##    theme_black()
dev.off()

#aa<-aa %>% filter(Labor!="NA" & status!="NA")


## Make a simple plot here:

fname=paste0(outFolder,"UMAP_LocationHarmony.png");
png(fname,width=1600,height=1200)
#aa$seurat_clusters <- clust2name[aa$seurat_clusters]
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Location)) +
    geom_point(size=0.1) +
    ##    scale_color_manual(values=group.colors) +
    theme_bw()+
    scale_colour_manual(values=c("CAM"="#F4B183","PVBP"="#8AD2CD"))+
    #scale_color_manual(values=c("Control"="#333399","E. coli"="#A50021"))+
    #theme(legend.text=element_text(size=30,face="bold"), axis.text=element_text(size=30,face="bold"), axis.title=element_text(size=20,face="bold"))+
    ##    facet_wrap(~LocTime) +
    guides(colour = guide_legend(override.aes = list(size=10)),title="Location")+
    theme(text = element_text(size=30,face = "bold"),
          plot.title = element_text(size = 25, face = "bold"),
          legend.title=element_text(size=25,face="bold"),
          legend.text=element_text(size=25,face="bold"))
p1
##    theme_black()
dev.off()



## Make a simple plot here:
fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Labor)) +
    geom_point(size=0.1) +
    scale_color_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=10)),title="Condition") +
    ##    facet_wrap(~LocTime) +
    theme_bw()+
    theme(text = element_text(size=30,face = "bold"),
          plot.title = element_text(size = 25, face = "bold"),
          legend.title=element_text(size=25,face="bold"), 
          legend.text=element_text(size=25,face="bold"))

p1
##    theme_black()
dev.off()


aa<-aa %>% filter(Labor!="NA" & status!="NA")

## Make a simple plot here:
fname=paste0(outFolder,"UMAP_StatusSoC_Harmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=status)) +
    geom_point(size=0.1) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=10),title="Status")) +
    theme_bw()+
    theme(text = element_text(size=30,face = "bold"),
          plot.title = element_text(size = 25, face = "bold"),
          legend.title=element_text(size=25,face="bold"), 
          legend.text=element_text(size=25,face="bold"))
##    facet_wrap(~LocTime) +

p1
##    theme_black()
dev.off()




aa<-aa %>% filter(Labor!="NA" & status!="NA")
## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.Cell_annotation.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
    geom_point(size=0.1) +
    ##    scale_color_manual(values=group.colors) +
    scale_color_manual(values=cluster.Colors2) +
    guides(colour = guide_legend(override.aes = list(size=15),title="Cell type")) +
    #facet_wrap(~Location) +
    facet_grid(Labor ~ Location) +
    theme_bw()+
    theme(text = element_text(size=30,face = "bold"),
          plot.title = element_text(size = 25, face = "bold"),
          legend.title=element_text(size=25,face="bold"), 
          legend.text=element_text(size=25,face="bold"))
p1
##    theme_black()
dev.off()




# fname=paste0(outFolder,"UMAP_Location.Barplot.pdf");
# pdf(fname,width=10,height=6)
# p2 <- ggplot(aa,aes(x=seurat_clusters,fill=SNG.BEST.GUESS)) +
#     geom_bar(position="stack") +
#     ##    scale_color_manual(values=group.colors) +
#     guides(colour = guide_legend(override.aes = list(size=10),title="Cell origin")) +
#     facet_grid(.~Location) + coord_flip() +
#     theme_bw()
# p2
# ##    theme_black()
# dev.off()


fname=paste0(outFolder,"barplot_pregnancyID.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=seurat_clusters,fill=Pcase)) +
    geom_bar(position="stack") +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=15),title="Pregnancy ID")) +
    facet_grid(.~Location) + coord_flip() +
    theme_bw()
p2
##    theme_black()
dev.off()



pregnancy.Colors<-cluster.Colors[1:length(unique(aa$Pcase))]
names(pregnancy.Colors)<-unique(aa$Pcase)

fname=paste0(outFolder,"UMAP_PregnancyID.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Pcase)) +
    scale_color_manual(values=pregnancy.Colors) +
    geom_point(size=0.1) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=10),title="Pregnancy ID")) +
    theme_bw()+
    theme(text = element_text(size=30,face = "bold"),
          plot.title = element_text(size = 25, face = "bold"),
          legend.title=element_text(size=25,face="bold"), 
          legend.text=element_text(size=25,face="bold"))
##    facet_wrap(~LocTime) +

p1
##    theme_black()
dev.off()


#aa<-aa %>% filter(Origin!="NA")




# cell_counts<-table(aa$cluster_name)
# write.csv(cell_counts,file=paste0(outFolder,"cluster_cell_counts.csv"))


### END- HERE ###
########################################################
sapply(as.character(c(0.5,0.6,0.7,0.8,"1.0",1.1)),function(x)
{
    
    outFolder=paste0("./5_harmony_cellClass_PBMC_plots_res",x,"/")
    system(paste0("cp -r ", outFolder," /nfs/rprscratch/wwwShare/azam/parturition_project_01.19.2022/"))})


