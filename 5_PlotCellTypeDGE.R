#############################################################
###  Plot top 20 cluster markers (umap plot / dot plot)
### 
#############################################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
#################


# clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# 
# names(clust2Names)<-as.character(c(0:23))



#outFolder='./5_harmonyClustersDGE_plots/'
#outFolder='./5_harmonyClustersDGE_SoupX_plots/'
outFolder="./5_harmonyClusters_withcovid19control_DGE_plots/"
system(paste0("mkdir -p ", outFolder))

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)



#sc <- read_rds("4_harmony_res0.8/sc.NormByLocationRep.Harmony.rds")
#sc <- read_rds("4_harmony_SoupX_res1.0/sc.NormByLocationRep.Harmony.rds")
sc <- read_rds("4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.rds")


# sc <- read_rds("4_harmony_res0.6/sc.NormByLibrary.Harmony.StringentFiltering.rds")
# sc <- subset(sc, Condition=="Control")
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

# anno<- anno %>%select(kbid,gene=gene_name)
# top20<-top20 %>% left_join(anno)
#top20<-read_rds("./5_harmonyClustersDGE/top20.rds")



# m2<-read_tsv("./5_harmonyClusters_SoupX_DGE/ClusterDEG.tsv")
# m2 <- m2 %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)
# top20 <- m2 %>% top_n(n = 20, wt = avg_log2FC)


#outFolder="./5_PlotClustersDGE/"
# outFolder="./5_PlotClustersSubtyDGE/"
# system(paste0("mkdir -p ", outFolder))



dim(sc)

table(sc$Library)

table(sc$Location) 



#m2 = read_tsv("5_harmonyClustersDGE/ClusterDEG.tsv")
#m2<-read_tsv("./5_harmonyClusters_SoupX_DGE/ClusterDEG.tsv")
m2<-read_tsv("./5_harmonyClusters_withcovid19control_DGE/ClusterDEG.tsv")
Hmax=log2(max(m2$cluster)+1)

# m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.2) %>%
#   group_by(gene) %>%
#   mutate(H=log2(length(cluster))) %>%
#   filter(H<=1) %>%
#   ungroup()

m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.25) %>%
  group_by(gene) %>%
  mutate(H=log2(length(cluster))) %>%
  filter(H<=1) %>%
  ungroup()

table(m3$cluster)


top20 <- m3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) %>% ungroup()
table(top20$cluster)


aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Rep","SNG.BEST.GUESS","cluster_name"))

myscale = 1/colSums(sc@assays$RNA@counts)*1000000

for(i in unique(top20$cluster)){
  ##    
  genesel <- filter(top20,cluster==i)
  ##        FetchData(sc,genesel$gene)
  genexpr <- map_dfr(1:nrow(genesel),function(i){
    aux <- FetchData(sc,genesel$gene[i])
    colnames(aux) <- "expr"
    aux <- cbind(aa,aux*myscale)
    #aux$gene=genesel$gene[i]
    aux$symbol=genesel$symbol[i] #genesel$symbol[i]
    aux
  })
  ##    
  cat(i,dim(genexpr),"\n")    
  ##
  fname=paste0(outFolder,i,"_umap.png");
  #fname=paste0(outFolder,"cluster_",i,"_umap.png");
  p1 <- genexpr %>% arrange(symbol,expr) %>%
    ggplot(aes(UMAP_1,UMAP_2,color=log10(0.1+expr))) +
    geom_point(size=0.1) +
    ##        scale_fill_viridis_c(option = "plasma") +
    scale_color_gradient(low = "lightblue", high = "darkred") +
    facet_wrap(~symbol) +
    theme_bw()
  ggsave(fname,p1,width=15,height=10)
  ##dev.off()
  cat(fname,"\n")
}


# top20 <- m3 %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC) %>% ungroup()
# table(top20$cluster)
# 
# 
# for(i in unique(top20$cluster)){
#   ##    
#   genesel <- filter(top20,cluster==i)
#   genexpr <- map_dfr(1:nrow(genesel),function(i){
#     aux <- FetchData(sc,genesel$gene[i])
#     colnames(aux) <- "expr"
#     aux <- cbind(aa,aux*myscale)
#     aux$gene=genesel$gene[i]
#     #aux$symbol=genesel$symbol[i]
#     aux
#   })
#   ##    
#   cat(i,dim(genexpr),"\n")
#   ##
#   sum_rec <- genexpr %>% dplyr::group_by(symbol,seurat_clusters) %>% dplyr::summarize(Prop=mean(expr>0),Expr=mean(expr[expr>0]))
#   ##
#   #fname=paste0(outFolder,"cluster_",i,"_dotplot2.png");
#   fname=paste0(outFolder,i,"_dotplot.png");
#   p0 <- ggplot(sum_rec,aes(x=seurat_clusters,y=symbol,color=Expr,size=Prop)) +
#     geom_point() +
#     scale_size_area(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),max_size=6,na.value=0) +
#     scale_color_distiller(trans='log10', palette = "RdYlBu",direction = -1,na.value=0,limits=c(10,NA)) +
#     ##facet_grid( .~ symbol) +
#     theme_classic() +
#     labs(y="Gene Symbol",x="Cell type",size="Prop.",color="Expr.(TPM)") + 
#     theme(axis.text.x = element_text(angle = 45,hjust=1)) 
#   ggsave(fname,p0,width=8,height=10)
#   ##dev.off()
#   cat(fname,"\n")
# }




# clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#               "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast",
#               "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
# 
# names(clust2Name)<-c(0:23)
# clust2Name<-paste0(names(clust2Name),"_",clust2Name)
# 
# 
# 
# fname=paste0(outFolder,"UMAP_Harmony.pdf");
# pdf(fname,width=5,height=5)
# DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# dev.off()
# 
# 
# 
# 
# 
# scale.data<-sc@assays$RNA@scale.data
# rw<-anno$gene_name[which(anno$kbid %in% rownames(scale.data))]
# names(rw)<-anno$kbid[which(anno$kbid %in% rownames(scale.data))]
# scale.data<-scale.data[names(rw),]
# rownames(scale.data)<-as.character(rw)
# sc@assays$RNA@scale.data<-scale.data
# #colnames(scale.data)<-newsc$seurat_clusters[colnames(scale.data)]
# 
# m2 = read_tsv("5_harmonyClustersDGE/ClusterDEG.tsv")
# m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.5) %>%
#   group_by(gene) %>%
#   mutate(H=log2(length(cluster))) %>%
#   filter(H<=1) %>%
#   ungroup()
# 
# top10 <- m3 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# 
# fname=paste0(outFolder,"Heatmap-top10.pdf");
# pdf(fname,width=35,height=25)
# DoHeatmap(sc, features = top10$symbol)+ #+ NoLegend()
# theme(text = element_text(size = 10),axis.text.y = element_text(size = 10)) + NoLegend()#
# dev.off()
# 
# 
# 


#system(paste0("mkdir -p ", outFolder,"CAM/"))


for(i in unique(top20$cluster)){
  ##    
  genesel <- filter(top20,cluster==i)
  genexpr <- map_dfr(1:nrow(genesel),function(i){
    aux <- FetchData(sc,genesel$gene[i])
    colnames(aux) <- "expr"
    aux <- cbind(aa,aux*myscale)
    #aux$gene=genesel$gene[i]
    aux$symbol=genesel$symbol[i]
    aux
  })
  ##    
  cat(i,dim(genexpr),"\n")
  ##
  sum_rec <- genexpr %>% dplyr::group_by(symbol,seurat_clusters) %>% dplyr::summarize(Prop=mean(expr>0),Expr=mean(expr[expr>0]))
  ##
  #fname=paste0(outFolder,"cluster_",i,"_dotplot2.png");
  fname=paste0(outFolder,i,"_dotplot.png");
  p0 <- ggplot(sum_rec,aes(x=seurat_clusters,y=symbol,color=Expr,size=Prop)) +
    geom_point() +
    scale_size_area(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),max_size=6,na.value=0) +
    scale_color_distiller(trans='log10', palette = "RdYlBu",direction = -1,na.value=0,limits=c(10,NA)) +
    ##facet_grid( .~ symbol) +
    theme_classic() +
    labs(y="Gene Symbol",x="Cell type",size="Prop.",color="Expr.(TPM)") + 
    theme(axis.text.x = element_text(angle = 45,hjust=1)) 
  ggsave(fname,p0,width=10,height=12)
  ##dev.off()
  cat(fname,"\n")
}



