library(tidyverse)


##################################################################
# comparison between between single cell and bulk datasets reanalysis: 
##################################################################

outFolder="12_comparison_with_bulk_Roger3/"
system(paste0("mkdir -p ", outFolder))

# cell type labels
cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
#cell.type.annotation$color[31]<-"#8B0000"
#write_tsv(cell.type.annotation,file="cell.type.annotation.v2.tsv")
#rownames(cell.type.annotation)<-NULL
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)


# single cell DEGs
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv")
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
res$cname<-paste0(res$Cell_type,"_",res$Origin)
#res <-res %>% filter(!is.na(log2FoldChange) & !is.na(padj))
res <-res %>% filter(!is.na(log2FoldChange))
dim(res)

cluster.Colors<-cell.type.annotation$color
names(cluster.Colors)<-unique(res$cname)
cluster.Colors<-cluster.Colors[1:length(unique(res$cname))]
#celltype_DE<-table(res$Cell_type,res$Location)

#ENTREZID id 
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

#colnames(res)[1]<-"gene_name"

res <- res %>% left_join(eg) ## %>% filter(!is.na(ENTREZID))
dim(res)


## Where was this file produced!!!????
##Open reslist 
bulkres <- read_rds(file="2022-08-29-NewBulk/reslist.rds");

scres <- mutate(res,cnameLoc=paste(Location,cname,sep="_"),zscore=log2FoldChange/lfcSE)

scnames <- names(table(scres$cnameLoc))

##scres <- res %>% group_by(cnameLoc) %>% group_split()

scnames2 <- scnames[grepl("CAM",scnames)]

cor.mat = sapply(bulkres,function(bulk){
  ##bulk <- bulkres[[1]]
  sapply(scnames2,function(scn){
    ##scn="CAM_4:Macrophage-2 (Hofbauer)_M"
    scdeg <- filter(scres,cnameLoc==scn)
    aux<-inner_join(bulk,scdeg,by="ENTREZID")
    ct = cor.test(aux$zscore.x,aux$zscore.y,method = "pearson")
    est=ct$estimate
    if(ct$p.value>0.05){est=NA}
    est
  })
})

cor.mat2 <- cor.mat[,grepl("PTLDT",colnames(cor.mat))]

rownames(cor.mat2) <- gsub(".cor","",rownames(cor.mat2))

## probably drop the PVBP
pdf("test2.pdf",width =4,height=8)
pheatmap::pheatmap(cor.mat2,cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()
