rm(list=ls())

library(tidyverse)
library(factoextra)
library(clusterProfiler)
library(annotate)


outFolder<-"./15_ROC_Signatures_DEGs_CAM_Blood_v4/"
system(paste0("mkdir -p ",outFolder))

# cell type labels
cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)

##### DEGs single cell data to build signatures
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv")
# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","cluster","Origin"),sep="_",remove=FALSE) #Cell_type
res<-res %>% filter(!is.na(padj) & padj<0.1)
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
dat <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

# up-regulated DEGs.  This is with absolute value. I may also want to include more than the top 20. 
#dat <- dat %>% group_by(cluster) %>% top_n(n = 20, wt = log2FoldChange) %>% ungroup()
# both up-regulated and down-regulated DEGs. 
dat<- dat %>%filter( Location=="CAM",padj<0.1,abs(log2FoldChange)>0.5)
##dat <- dat %>% group_by(cluster) %>% top_n(n = 20, wt = abs(log2FoldChange)) %>% ungroup()
table(dat$cluster)
dat$cluster<-clust2Names[as.character(dat$cluster)]
dat$cluster<-paste0(dat$cluster,"_",dat$Origin)
# filter small clusters
#dat<-dat %>% filter(cluster %in% cluster_filter )
dat<-dat %>% dplyr::select(cluster,ENTREZID,gene=gene_name,avg_logFC=log2FoldChange,p_val=pvalue,p_val_adj=padj)


tc <- table(dat$cluster) 
tc
tc <- names(which(tc>10))

dat <- dat %>% filter(cluster %in% tc)


### Generage weights signed with the direction of labor, to build a metagene per cell-type. 
generate_gene_cell_matrix_v2<-function(dat) 
{
  dat<-as.data.frame(dat)
  #edge_list<-dat[,c("ENTREZID","cluster","avg_logFC")]
  edge_list<-dat[,c("gene","cluster","avg_logFC")]
  edge_list<-as.data.frame(edge_list)
  colnames(edge_list)[3]<-"weight"
  edge_list<-edge_list[!duplicated(edge_list), ]
  #edge_list$weight<-1
  cl<-unique(as.character(edge_list[,"cluster"]))
  rw<-unique(as.character(edge_list[,"gene"]))
  edge_list$cluster<-as.character(edge_list$cluster)
  gene_cell_matrix<-matrix(0,nrow=length(rw),ncol = length(cl))
  rownames(gene_cell_matrix)<-rw
  colnames(gene_cell_matrix)<-cl
  gene_cell_matrix[as.matrix(edge_list[,1:2])] <- edge_list[,3]
  gene_cell_matrix[which(gene_cell_matrix>0)]<- 1
  gene_cell_matrix[which(gene_cell_matrix<0)]<- -1
  return(gene_cell_matrix)
}


# ##### Need to update to the other dataset. 
# #anoSC2,esetSC2
# load("GSE96083/ano_eset_GSE96083.RData")
# eset<-esetSC2
# anpack="org.Hs.eg.db"
# rownames(eset)<-gsub("_at","",rownames(eset))
# SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
# eset=eset[!is.na(SYMBOLS),]
# SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
# eset=eset[order(apply(eset,1,mean),decreasing=TRUE),]
# SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
# eset=eset[!duplicated(SYMBOLS),]
# SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
# rownames(eset)<-SYMBOLS
# count_bulk<-eset
# meta<-anoSC2 
# 
# rownames(meta)<-meta$SampleID

################################
###############################
load("./preterm_adi/DREAM_PPROM_toroger.RData")
##anoSC2 and esetSC2
GA=factor(ifelse(anoSC2$GA<=25,"T1","T2"))
anoSC2$Time <- GA 
## T1 earlier time point. 
## batch is the two dataswets, HTA20 is the PRB, HG21ST is the other. 
ano <- anoSC2 ##%>% filter(Time=="T2",Platform=="HG21ST")
eset <- esetSC2[,ano$SampleID]

rownames(ano) <- ano$SampleID 

## change rowname to SYMBOL to merge data
anpack="org.Hs.eg.db"
rownames(eset)<-gsub("_at","",rownames(eset))
SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
eset=eset[!is.na(SYMBOLS),]
SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
eset=eset[order(apply(eset,1,mean),decreasing=TRUE),]
SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
eset=eset[!duplicated(SYMBOLS),]
SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
rownames(eset)<-SYMBOLS
bulk<-eset


sample_names<-ano$Group[which(ano$SampleID %in% colnames(bulk))]
names(sample_names)<-ano$SampleID[which(ano$SampleID %in% colnames(bulk))]

# CREATE SINGLE CELL SIGNATURES ASSOCIATED WITH LABOR
# both up and down regulated DEGs 
gene_cell_matrix<-generate_gene_cell_matrix_v2(dat) #gene = "ENSG",clustertype = "cluster",

## Find genes in common bulk and sc signatures matrix. 
rw1<-rownames(bulk)
rw2<-rownames(gene_cell_matrix)
rw<-intersect(rw1,rw2)
gene_cell_matrix<-gene_cell_matrix[rw,]
bulk2 <- bulk[rw,]
gene_cell_matrix<-t(gene_cell_matrix)  # cell types * gene



# Building the actual metagene of the bulk data using the gene signatures. 
Metagene<-gene_cell_matrix  %*% as.matrix(bulk2)
Metagene<-t(Metagene)


#Metgene <- log10(as.matrix(Metgene)+1)
groups<-unique(ano$Group)
groups

##Rcountcopy<-Rcount ##bulk
##Metgenecopy<-Metgene ##Metagene

library(pROC)

stopifnot(rownames(ano) == rownames(Metagene))

system(paste0("mkdir -p ",outFolder,"/","sep/"))

## Subset metagene for timepoint, batch, and condition
tp <- "T1"
batch <- "HG21ST"
condition <- "PPROM"


cc <- expand.grid(Time = c("T1", "T2"), 
                    Platform = c("HG21ST", "HTA20"), 
                    Condition = c("PPROM", "sPTD"),stringsAsFactors = FALSE)



rocall <- map_dfr(seq_len(nrow(cc)), function(i) {
  tp <- cc$Time[i]
  batch <- cc$Platform[i]
  condition <- cc$Condition[i]
    ##
  anocopy <- ano %>% 
    filter(Time == tp, Platform == batch, Group %in% c("Control", condition))
  cat(paste("##",condition, batch, tp),"\n")
  cat(unique(anocopy$Group),"\n")
  ##
  Metagenecopy <- Metagene[anocopy$SampleID,]
  ##
  roc_celltype <- map(colnames(Metagenecopy), function(celltype) {
    cat(paste(condition, batch, tp, celltype), "\n")
    pred <- Metagenecopy[,celltype]
    aucs <- round(ci.auc(roc(response = as.numeric(anocopy$Group == condition), predictor = pred)), 3)
    ##
      RC <- roc(response = as.numeric(anocopy$Group == condition), predictor = pred)
      # pdf(paste0(outFolder,"/sep/","ROC_",paste(condition,batch,tp,celltype,sep="_"),".pdf"))
      # plot(1-RC$specificities, RC$sensitivities, col="black", type="l", lwd=2, xlab="False positive rate", ylab="Sensitivity")
      # abline(0,1)
      # text(0.5, 0.05, paste("AUC: ", aucs[1], aucs[2], aucs[3]))
      # title(paste(condition, batch, tp, celltype))
      # dev.off()
##    tibble(auc.l = aucs[1], auc = aucs[2], auc.u = aucs[3],
##           celltype = celltype, condition = condition, batch = batch, tp = tp)
      tibble(FPR=rev(1-RC$specificities),Sensitivity=rev(RC$sensitivities),
            celltype = celltype, condition = condition, cohort = batch, tp = tp,
            auc.l = aucs[1], auc = aucs[2], auc.u = aucs[3])
  })
  ##
  roc_celltype
})

##write_tsv(rocall,paste0(outFolder,"/sep/","auc_table_all.tsv"))
##mauc %>% filter(auc.l>0.5)
pdf(paste0(outFolder,"ROC_sep_gg.pdf"),height=26,width=7)
rocall %>% ggplot(aes(x=FPR,y=Sensitivity,color=cohort,lty=tp)) +
  scale_color_manual(values = c("#1144DD", "#AA5522")) + ## "#E7B800"
  geom_step() + 
  geom_abline(slope=1,intercept=0) + 
  facet_grid(celltype ~ condition) +
  theme_bw()
dev.off()



### Combined batches
system(paste0("mkdir -p ",outFolder,"/","combined/"))
cc <- expand.grid(Time = c("T1", "T2"), 
                  Condition = c("PPROM", "sPTD"),stringsAsFactors = FALSE)



rocall <- map_dfr(seq_len(nrow(cc)), function(i) {
  tp <- cc$Time[i]
  condition <- cc$Condition[i]
  ##
  anocopy <- ano %>% 
    filter(Time == tp, Group %in% c("Control", condition))
  cat(paste("##",condition, tp),"\n")
  cat(unique(anocopy$Group),"\n")
  ##
  Metagenecopy <- Metagene[anocopy$SampleID,]
  ##
  roc_celltype <- map(colnames(Metagenecopy), function(celltype) {
    cat(paste(condition, tp, celltype), "\n")
    pred <- Metagenecopy[,celltype]
    aucs <- round(ci.auc(roc(response = as.numeric(anocopy$Group == condition), predictor = pred)), 3)
    ##
    RC <- roc(response = as.numeric(anocopy$Group == condition), predictor = pred)
    # pdf(paste0(outFolder,"/sep/","ROC_",paste(condition,batch,tp,celltype,sep="_"),".pdf"))
    # plot(1-RC$specificities, RC$sensitivities, col="black", type="l", lwd=2, xlab="False positive rate", ylab="Sensitivity")
    # abline(0,1)
    # text(0.5, 0.05, paste("AUC: ", aucs[1], aucs[2], aucs[3]))
    # title(paste(condition, batch, tp, celltype))
    # dev.off()
    ##    tibble(auc.l = aucs[1], auc = aucs[2], auc.u = aucs[3],
    ##           celltype = celltype, condition = condition, batch = batch, tp = tp)
    tibble(FPR=rev(1-RC$specificities),Sensitivity=rev(RC$sensitivities),
           celltype = celltype, condition = condition, tp = tp,
           auc.l = aucs[1], auc = aucs[2], auc.u = aucs[3])
  })
  ##
  roc_celltype
})

##write_tsv(rocall,paste0(outFolder,"/sep/","auc_table_all.tsv"))
##mauc %>% filter(auc.l>0.5)
pdf(paste0(outFolder,"ROC_combined_gg.pdf"),height=26,width=7)
rocall %>% ggplot(aes(x=FPR,y=Sensitivity,color=tp)) +
  scale_color_manual(values = c("#1144DD", "#AA5522")) + ## "#E7B800"
  geom_step() + 
  geom_abline(slope=1,intercept=0) + 
  facet_grid(celltype ~ condition) +
  theme_bw()
dev.off()


