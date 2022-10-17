###################################################
### Finding differentially expressed genes ###
###################################################


library(tidyverse)
##library(knitr)
library(DESeq2)
##library(annotables)
library(qqman)

library(Seurat)
library("vsn")

library("BiocParallel")
register(MulticoreParam(16))

parturion_samples<-read_delim("parturition_cv.txt")

# cc <- read_tsv("../Covid19.Samples.txt") %>%
#     mutate(Sample_ID=paste(Pregnancy_ID,Origin,sep="-"))
# head(cc)
# demographic <- read_tsv("Covid19_demographic.txt") %>%
#     mutate(Pregnancy_ID=paste0("HPL",HPL))
# head(demographic)
# cc %>% inner_join(demographic)
# filter_sample<-c("HPL20921","HPL20874")#,"HPL20922")
# Severe_samples<-c("HPL20919","HPL20875")
# filter_sample<-c(filter_sample,Severe_samples)
# new_pg<-c("HPL20919", "HPL20922" ,"HPL20929", "HPL20928")
# filter sample "HPL20874"


#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_noBatchCorrect/")
#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_library_delroute/")
#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_library_delroute_Age/")
#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_library/")
#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res0.8_library/")
outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library_v2/")
#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.5_library/")

#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_SoupX_library/")
#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_SoupX_noBatchCorrect/")
system(paste0("mkdir -p ",outFolder))


#sc <- read_rds("4_harmony_res0.8/sc.NormByLocationRep.Harmony.rds")
#sc <- read_rds("4_harmony_SoupX_res1.0/sc.NormByLocationRep.Harmony.rds")
#sc <- read_rds(paste0("4_harmony_with_covidcontrol_res0.8/sc.NormByLocationRep.Harmony.rds"))
sc <- read_rds(paste0("4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.rds"))
#sc <- read_rds(paste0("4_harmony_with_covidcontrol_res1.5/sc.NormByLocationRep.Harmony.rds"))


#md<-md %>% filter(!is.na(Labor))

#md <- filter(md,!Pregnancy_ID %in% filter_sample)

## Load gene annotations. 
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


gene_symbol <- anno$gene_name
names(gene_symbol) <- anno$kbid

sc<-sc[anno$kbid,]

stopifnot(identical(rownames(sc),anno$kbid))

md <- sc@meta.data
md<-md %>% separate(SNG.BEST.GUESS,into=c("Pcase", "Origin",sep="-",remove=FALSE))

md<-md %>% filter(!is.na(Condition))
md<-md %>% filter (!Pregnancy_ID %in% c("HPL20888", "HPL20874", "HPL20922"))

parturion_samples_data<- parturion_samples %>% select(Pregnancy_ID,Condition_merge=Labor)

md<-md %>% rownames_to_column("barcode") %>% inner_join(parturion_samples_data)
rownames(md)<-md$barcode
#md$Condition [which(md$Condition=="Control")]<-"TIL"

cSize <- md %>% dplyr::count(Location,seurat_clusters,Pregnancy_ID,Origin,Labor)

#sSize <- cSize %>% filter(n>20) %>% dplyr::count(Location,seurat_clusters,Origin,Labor) %>% filter(n>2)
sSize <- cSize %>% filter(n>20) %>% dplyr::count(Location,seurat_clusters,Origin,Labor) %>% filter(n>2)

kSize <- sSize %>% dplyr::count(Location,seurat_clusters,Origin) %>% filter(n>=2) %>%
    mutate(cname=paste(Location,seurat_clusters,Origin,sep="_"))

kSize


md<- md %>% mutate(cname=paste(Location,seurat_clusters,Origin,sep="_"))

cname_cellcount<-table(md$cname)


resList <- lapply(1:nrow(kSize),function(ii){

    ##ii <- 1
    cname <- kSize$cname[ii]
    orig<-unlist(strsplit(cname,"_"))[3]
    cat("#",cname,"\n")
    cat("#",ii,"\n")
    ##
    md_i <- md %>%
        filter(
            Location==kSize$Location[ii],
            seurat_clusters==kSize$seurat_clusters[ii],
            Origin==kSize$Origin[ii]) 
    
    ## samples (control/case)
    #bti <- md_i %>% transmute(bti=paste(Labor,Pregnancy_ID,Library,delRoute,Age, sep="_")) %>% unlist %>% factor
    bti <- md_i %>% transmute(bti=paste(Condition_merge,Pregnancy_ID,Library,delRoute,Age, sep="_")) %>% unlist %>% factor
    #bti <- md_i %>% transmute(bti=paste(Labor,Library,sep="_")) %>% unlist %>% factor
    #bti <- md_i %>% transmute(bti=paste(Labor,sep="_")) %>% unlist %>% factor
    
    ## data selected
    all_i <- sc@assays$RNA@data[,rownames(md_i)]
    ##
    
    ##
    X <- model.matrix( ~ 0 + bti)
    qr.X <- qr(X)
    qr.X$rank
    dim(X)
    YtX <- all_i %*% X
    YtX <- as.matrix(YtX)
    dim(YtX)
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    bti2 <- gsub("bti","",colnames(YtX))
    colnames(YtX) <- bti2
    ##
    cmat<-YtX

    
    anno_i <- tibble(kbid=rownames(cmat),rs=rowSums(cmat),nz=rowSums(cmat>0)) %>%
        inner_join(anno %>% dplyr::select(kbid,Chr,TSS,Strand,gene_name)) %>%
        filter(Chr %in% paste0("chr", 1:22), rs > 10, nz > 3) ## keep only autosomal
    ##
    ##
    table(anno_i$Chr)
    ##
    # genesel <- (unique(anno_i$kbid))
    # cmat <- cmat[genesel,]
    dim(cmat)
    
    
    data<-cmat
 
    ##
    ##Create sample table
    cn<-colnames(cmat)
    x<-strsplit(cn,"_")
    ##
    cvt <- data.frame(matrix(unlist(x), nrow=length(x), byrow=T),stringsAsFactors=FALSE)
    #colnames(cvt)<-c("Group","Indiv","Library")
    
    #colnames(cvt)<-c("Group","Indiv")
    
    #colnames(cvt)<-c("Group","Indiv","Library","delRoute")
    colnames(cvt)<-c("Group","Indiv","Library","delRoute","Age")
    cvt$Age<-as.numeric(cvt$Age)
    #colnames(cvt)<-c("Group","Library")
    #colnames(cvt)<-c("Group")
    ##
    # cvt$Library<-Library
    # cvt$Library<-as.factor(cvt$Library)
    cvt$Group = relevel(factor(cvt$Group),"TNL")
    # cvt$version<-1
    # cvt$version[cvt$Indiv %in% new_pg]<-2
    #    
    # cvt$version<-as.factor(cvt$version)
    # ##
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    dds <- DESeqDataSetFromMatrix(round(cmat),cvt, ~Library+Group)  
    #dds <- DESeqDataSetFromMatrix(round(cmat),cvt, ~ Group)  
    #dds <- DESeqDataSetFromMatrix(round(cmat),cvt, ~version+Group)  
    #dds <- DESeqDataSetFromMatrix(round(cmat),cvt, ~Library+delRoute+Age+Group)  
    dds <- DESeq(dds,parallel=TRUE)
    
    ############################### 
    # PCA
    ###############################
    # vsd <- vst(dds, blind=FALSE)
    # rld <- rlog(dds, blind=FALSE)
    # #head(assay(vsd),3)
    # ntd <- normTransform(dds)
    
    # meanSdPlot(assay(ntd))
    # meanSdPlot(assay(vsd))
    # meanSdPlot(assay(rld))

    #
    # save DEseq object.
    fname <- paste0(outFolder,cname,"_dds_",Sys.Date(),".rds")
    fname
    write_rds(dds,fname)
    #
    # Parse the results
    # ========================================================
    #
    res <- results(dds)
    myres <- as.data.frame(res) %>%
        rownames_to_column("kbid") %>%
        left_join(anno_i)
    myres$cname <- cname
    ##
    write_tsv(myres,paste0(outFolder,cname,".",Sys.Date(),".txt"))
    ##
    nres<-nrow(myres %>% filter(padj<.1,abs(log2FoldChange)>1))
    cat("# Total sig for",cname,": ",nres,"\n")
    myres
 #}
#    cname
 })

 res <- do.call(rbind,resList)
# 
 sum(res$padj<0.1,na.rm=TRUE)
# 
# Sys.Date()
# 
 
 # cell.type.annotation<-read.delim("cell.type.annotation.txt")
 # 
 # 
 # clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
 # clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
 # names(clust2Names)<-as.character(cell.type.annotation$Cluster)
 # 
 # res$Cell_type<-clust2Names[res$Cluster]
 # res$cname<-paste0(res$Cell_type,"_",res$Location,"_",res$Origin)  
 
 
 
system(paste0("mkdir -p ",outFolder))

res %>% filter(padj<0.1,abs(log2FoldChange)>0) %>% dplyr::count(cname) %>%
    write_tsv(paste0(outFolder,"Summary.FDR.",Sys.Date(),".tsv"))

res %>% write_tsv(paste0(outFolder,"ALL.combined.",Sys.Date(),".tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>0.0) %>%
    write_tsv(paste0(outFolder,"SIG.combined.",Sys.Date(),".tsv"))

res %>% filter(padj<0.05,abs(log2FoldChange)>1.0) %>% dplyr::select(kbid,gene_name,cname,log2FoldChange,padj)
# 


# ###
# res1 <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_library/ALL.combined.2022-03-23.tsv" )
# length(unique(res1$gene_name[res1$padj<0.1]))
# res2 <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv" )
# length(unique(res2$gene_name[res2$padj<0.1]))
# res3 <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.5_library/ALL.combined.2022-03-30.tsv" )
# length(unique(res3$gene_name[res3$padj<0.1]))


