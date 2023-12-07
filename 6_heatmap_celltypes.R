library(gplots)
library(tidyverse)
library(reshape2)
library(ggplot2)

outFolder="6_celltype_cor_heatmap_plot_Roger/"
system(paste0("mkdir -p ", outFolder))

# cell type labels


cell.type.annotation<-read.delim("cell.type.annotation.v2.txt")
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)

#cell.type.annotation$color<-as.character(cluster.Colors[1:33])
#write_tsv(cell.type.annotation,file="cell.type.annotation.v2.tsv")



# union of DE genes across cell types
res_de <-read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/SIG.combined.2022-03-29.tsv")
res_de <- res_de %>% separate(cname,c("Location","Cluster","Origin"),sep="_",remove=FALSE)
de_gene<-unique(res_de$kbid)
res_de$Cell_type<-clust2Names[res_de$Cluster]
res_de$cname<-paste0(res_de$Cell_type,"_",res_de$Location,"_",res_de$Origin)  

# calculating correlation based on union of DE genes 
res <-read_tsv("7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv")
res<-res %>% filter(!is.na(log2FoldChange) & !is.na(padj))
res<-res  %>% filter (kbid %in% de_gene)
res <- res %>% separate(cname,c("Location","Cluster","Origin"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cluster]
res$cname<-paste0(res$Cell_type,"_",res$Location,"_",res$Origin)  
### RPR changed this line. I think this will solve the problem of not getting the same plot. 
##clusters<-unique(res_de$cname)
###clusters<-unique(res$cname)
### RPR 2023/02/21 Additionally changed this to remove small clusters. 
tt <- table(res_de$cname)
clusters <- names(which(tt>10))

# clusters_location<-sapply(clusters,function(x){
#   y<-unlist(strsplit(x,"_"))
#   return(y[length(y)])
#   
# })
# clusters_location<-clusters_location[order(clusters_location,decreasing = TRUE)]
# clusters<-names(clusters_location)


cor_matrix<-matrix(NA,nrow=length(clusters),ncol=length(clusters))
for (i in 1:length(clusters))
{
  for (j in 1: length(clusters))
  {
      resi<-res %>% dplyr::filter(cname == clusters[i]) %>% dplyr::select(kbid,log2FoldChange,lfcSE)
      resj<-res %>% dplyr::filter(cname == clusters[j])%>% dplyr::select(kbid,log2FoldChange,lfcSE)
      colnames(resj)<-c("kbid","log2FoldChange2","lfcSE2")
      res_intersect<-resi%>% inner_join(resj)
      res_intersect<-res_intersect%>% dplyr::select(log2FoldChange,log2FoldChange2)
      if(nrow(res_intersect)>5)
      {
        cr<-cor(res_intersect)[1,2]
        cor_matrix[i,j]<-cr
        cor_matrix[j,i]<-cr
      }
      }
  }
 
rownames(cor_matrix)<-clusters
colnames(cor_matrix)<-clusters
write_rds(cor_matrix,file=paste0(outFolder,"cor_matrix.rds"))

write_rds(cor_matrix,file=paste0(outFolder,"cor_matrix_all.rds"))

# heatmap plot

library(pheatmap)

fname=paste0(outFolder,"heatmap_celltype_cor.pdf");
pdf(fname,width=8,height=8)
paletteLength<-50
my_palette <- colorRampPalette(colors = c("#E8F4F8","white", "#A50021"))(n = paletteLength)
myBreaks <- c(seq(min(cor_matrix), 0, length.out=ceiling(paletteLength/2) ), 
              seq(max(cor_matrix)/paletteLength, max(cor_matrix), length.out=floor(paletteLength/2)))

#myBreaks <- c(seq( min(cor_matrix),0,length.out=20), seq(0,max(cor_matrix), length.out=30) )


pheatmap(cor_matrix,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize=12)
dev.off()


fname=paste0(outFolder,"heatmap_celltype_cor_all.pdf");
pdf(fname,width=15,height=15)
paletteLength<-50
my_palette <- colorRampPalette(colors = c("#E8F4F8","white", "#A50021"))(n = paletteLength)
myBreaks <- c(seq(min(cor_matrix), 0, length.out=ceiling(paletteLength/2) ), 
              seq(max(cor_matrix)/paletteLength, max(cor_matrix), length.out=floor(paletteLength/2)))

#myBreaks <- c(seq( min(cor_matrix),0,length.out=20), seq(0,max(cor_matrix), length.out=30) )

pheatmap(cor_matrix,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize=10)
dev.off()

# fname=paste0(outFolder,"heatmap_celltype_cor2.pdf");
# pdf(fname,width=30,height=30)
# 
# pheatmap(cor_matrix,cluster_rows=TRUE,fontsize=20)
# dev.off()



### new plot #### 

# library(pheatmap)
# 
# fname=paste0(outFolder,"heatmap_celltype_cor3.pdf");
# pdf(fname,width=38,height=28)
# #pheatmap(cor_matrix,cluster_rows=TRUE,cluster_cols=TRUE,scale="none")
# #my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 201)
# paletteLength<-30
# # my_palette <- colorRampPalette(colors = c("#333399", "white", "#A50021"))(n = paletteLength+1)
# #  myBreaks <- c(seq(min(cor_matrix),0, length.out=ceiling(paletteLength/2) + 1), 
# #                seq(max(cor_matrix)/paletteLength, max(cor_matrix), length.out=floor(paletteLength/2)))
# 
# my_palette <- colorRampPalette(colors = c("#D8EAF0","white", "#A50021"))(n = paletteLength+1)
# myBreaks <- c(seq(min(cor_matrix),0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(cor_matrix)/paletteLength, max(cor_matrix), length.out=floor(paletteLength/2)))
# 
# pheatmap(cor_matrix,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize =10)
# dev.off()


############



clusters_location<-sapply(clusters,function(x){
  y<-unlist(strsplit(x,"_"))
  return(y[2])

})


clusters_number<-sapply(clusters,function(x){
  y<-unlist(strsplit(x,"_"))
  y<-unlist(strsplit(y,":"))
  return(as.numeric(y[1]))
  
})

# clusters_number<-clusters_number[order(clusters_number,decreasing = FALSE)]
# 
# cor_matrix<-cor_matrix[names(clusters_number),]
# cor_matrix_cervix<-cor_matrix[,clusters_location[colnames(cor_matrix)]=="Cervix"]



########################################################
# Same cell type across locations
# Decidua-Cervix
# Decidua-Myometrium
# Myometrium-Cervix
########################################################
CAM<-names(clusters_location)[which(clusters_location =="CAM")]
PVBP<-names(clusters_location)[which(clusters_location =="PVBP")]


#D_C<-names(clusters_location)[which(clusters_location !="Myometrium")]
cor_matrix_CAM_PVBP<-cor_matrix[CAM,PVBP]


library(pheatmap)

fname=paste0(outFolder,"heatmap_celltype_cor_CAM_PVBP.pdf");
pdf(fname,width=5,height=6)
paletteLength<-50
my_palette <- colorRampPalette(colors = c("#E8F4F8","white", "#A50021"))(n = paletteLength)
myBreaks <- c(seq(min(cor_matrix_CAM_PVBP), 0, length.out=ceiling(paletteLength/2) ), 
              seq(max(cor_matrix_CAM_PVBP)/paletteLength, max(cor_matrix_CAM_PVBP), length.out=floor(paletteLength/2)))

#myBreaks <- c(seq( min(cor_matrix),0,length.out=20), seq(0,max(cor_matrix), length.out=30) )

pheatmap(cor_matrix_CAM_PVBP,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize=10)
dev.off()



# pdf(paste0(outFolder,"boxplot_celltype_locations.pdf"),width=3,height = 3)
# ggplot(data_locations,aes(x=Location, y=Correlation,fill=Location)) + 
# geom_boxplot()+
# xlab("")+
# scale_fill_manual(values=c("Decidua-Cervix"="#74AF81","Decidua-Myometrium"="#B46FB5","Myometrium-Cervix"="#99BE7F") )+
# theme_bw()+
# theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.2))
# #theme(axis.text=element_text(size=30))
# #scale_color_manual(values=c("red", "green", "blue"))
# #facet_wrap(~Location, scale="free")
# dev.off()


########################################################
# within location across cell types

cor_matrix_CAM<-cor_matrix[CAM,CAM]
cor_matrix_PVBP<-cor_matrix[PVBP,PVBP]



library(pheatmap)

fname=paste0(outFolder,"heatmap_celltype_cor_CAM.pdf");
pdf(fname,width=12,height=13)
paletteLength<-50
my_palette <- colorRampPalette(colors = c("#E8F4F8", "white", "#A50021"))(n = paletteLength)
myBreaks <- c(seq(min(cor_matrix_CAM), 0, length.out=ceiling(paletteLength/2) ), 
              seq(max(cor_matrix_CAM)/paletteLength, max(cor_matrix_CAM), length.out=floor(paletteLength/2)))
pheatmap(cor_matrix_CAM,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize=10)

dev.off()


# fname=paste0(outFolder,"heatmap_celltype_cor_PVBP.pdf");
# pdf(fname,width=5,height=5)
# paletteLength<-10
# my_palette <- colorRampPalette(colors = c("#E8F4F8", "white", "#A50021"))(n = paletteLength)
# myBreaks <- c(seq(min(cor_matrix_PVBP), 0, length.out=ceiling(paletteLength/2) ), 
#               seq(max(cor_matrix_PVBP)/paletteLength, max(cor_matrix_PVBP), length.out=floor(paletteLength/2)))
# pheatmap(cor_matrix_PVBP,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize=10)
# dev.off()




