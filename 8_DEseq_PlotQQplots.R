library(tidyverse)
library(DESeq2)
library(qqman)
library(org.Hs.eg.db)
library(clusterProfiler)


outFolder <- paste0("./8_outputs_DESeq_batch_library_with_covidcontrol_Plots_v2/")
system(paste0("mkdir -p ",outFolder))

cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")

clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)



# cell.type.annotation<-read.delim("cell.type.annotation.txt")
# clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
# clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
# names(clust2Names)<-as.character(cell.type.annotation$Cluster)
 

#write.csv(cell.type.annotation,file="cell.type.annotation.csv")

########################################################
# load single cell data 

########################################################

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library_v2/ALL.combined.2022-05-04.tsv" )

res <- res %>% separate(cname,c("Location","Cluster","Origin"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cluster]
res$cname<-paste0(res$Cell_type,"_",res$Location,"_",res$Origin)
#write_tsv(res, paste0("7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.withcelltypes",Sys.Date(),".tsv"))


## cluster colors 

# clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# 
# names(clust2Names)<-as.character(c(0:23))
# cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
#                   "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
#                   "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF")#,"#6F00FF")
# names(cluster.Colors)<-clust2Names #as.character(c(0:29))#clust2Names #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")


cluster.Colors<-cell.type.annotation$color
names(cluster.Colors)<-clust2Names
# 
# 
# # cluster labels
# res$Cell_type<-clust2Names[res$Cell_type]



#ENTREZID
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))


res2 <- res %>% filter(!is.na(pvalue)) %>%
    arrange(pvalue) %>%
    group_by(Location,Cell_type,Origin) %>%
    mutate(r=rank(pvalue, ties.method = "random"),pexp=r/length(pvalue))


# new_names <- read_tsv("./ClusterAssignment.tsv")
# clust2Names <- new_names$scLabor_ID
# names(clust2Names) <- new_names$seurat_clusters
# cc <- new_names %>% dplyr::select(scLabor_ID,color) %>% unique 
# cluster.Colors <- cc$color
# names(cluster.Colors) <- cc$scLabor_ID


fname=paste0(outFolder,"all.qqplot.png");
png(fname,width=800,height=800)
qq(res$pvalue)
dev.off()


###################################################
# qqplot to show the p-values splited by Origin and Location  
###################################################
fname=paste0(outFolder,"split.qqplot.png");
p1 <- res2 %>%
    ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
    geom_point() +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    geom_abline(slope=1,intercept=0) +
    facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) + 
    theme_bw()

ggsave(fname,p1,width=9,height=6)


#############################################

# res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-16.tsv")
# res %>% filter(padj<=0.05)

