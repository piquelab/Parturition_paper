#setwd("//wsu/home/groups/prbgenomics/mouse_pilot/mouse_pilot_analysis/")
library(tidyverse)
##library(knitr)
library(DESeq2)
##library(annotables)
library(qqman)
library(tidyr)
library(Seurat)
library(reshape)

library("BiocParallel")
 register(MulticoreParam(16))




outFolder <- paste0("7_outputs_DESeq_sudobulk/")
system(paste0("mkdir -p ",outFolder))

parturion_samples<-read_delim("parturition_cv.txt")
## Load cells
sc <- read_rds(paste0("4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.rds"))


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



clust2Name<-rep("sudo",33)
names(clust2Name)<-as.character(c(0:32))

md$seurat_clusters<-clust2Name[md$seurat_clusters]
sc$seurat_clusters<-clust2Name[sc$seurat_clusters]


cSize <- md %>% dplyr::count(Location,seurat_clusters,Pregnancy_ID,Origin,Labor)

#sSize <- cSize %>% filter(n>20) %>% dplyr::count(Location,seurat_clusters,Origin,Labor) %>% filter(n>2)
sSize <- cSize %>% filter(n>20) %>% dplyr::count(Location,seurat_clusters,Origin,Labor) %>% filter(n>2)

# kSize <- sSize %>% dplyr::count(Location,seurat_clusters,Origin) %>% filter(n>=2) %>%
#     mutate(cname=paste(Location,seurat_clusters,Origin,sep="_"))

kSize <- sSize %>% dplyr::count(Location,seurat_clusters) %>% filter(n>=2) %>%
    mutate(cname=paste(Location,seurat_clusters,sep="_"))

kSize

kSize

resList <- lapply(1:nrow(kSize),function(ii){
    ##ii <- 1
    cname <- kSize$cname[ii]
    cat("#",cname,"\n")
    cat(ii)
    md_i <- md %>%
        filter(
            Location==kSize$Location[ii],
            seurat_clusters==kSize$seurat_clusters[ii]
            )

    bti <- md_i %>% transmute(bti=paste(Condition_merge,Pregnancy_ID,Library, sep="_")) %>% unlist %>% factor
    
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
    ##
    ##
    anno_i <- tibble(kbid=rownames(cmat),rs=rowSums(cmat),nz=rowSums(cmat>0)) %>%
        inner_join(anno %>% dplyr::select(kbid,Chr,TSS,Strand,gene_name)) %>%
        filter(Chr %in% paste0("chr", 1:22), rs > 10, nz > 3) ## keep only autosomal
    ##
    ##
    #table(anno_i$Chr)
    ##
    data<-cmat
    
    ##
    ##Create sample table
    cn<-colnames(cmat)
    x<-strsplit(cn,"_")
    ##
    cvt <- data.frame(matrix(unlist(x), nrow=length(x), byrow=T),stringsAsFactors=FALSE)
    
    ##
    colnames(cvt)<-c("Group","Indiv","Library")
    ##
    
    cvt$Group = relevel(factor(cvt$Group),"TNL")
    ##
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    
    dds <- DESeqDataSetFromMatrix(round(cmat),cvt, ~Library+Group)  
    dds <- DESeq(dds,parallel=TRUE)
    ##
    ## save DEseq object.    
    fname <- paste0(outFolder,cname,"_dds_",Sys.Date(),".rds")
    fname
    write_rds(dds,fname)
    ##
    ## Parse the results 
    ## ========================================================
    ##
    res <- results(dds)
    myres <- as.data.frame(res) %>%
        rownames_to_column("kbid") %>%
        left_join(anno_i)
    myres$cname <- cname
    ##
    write_tsv(myres,paste0(outFolder,cname,".",Sys.Date(),".txt"))
    ##
    nres<-nrow(myres %>% filter(padj<.1,abs(log2FoldChange)>0.5))
    cat("# Total sig for",cname,": ",nres,"\n")
    myres
})

res <- do.call(rbind,resList)

sum(res$padj<0.1,na.rm=TRUE)



Sys.Date()

system(paste0("mkdir -p ",outFolder))

res %>% filter(padj<0.1,abs(log2FoldChange)>0) %>% dplyr::count(cname) %>%
    write_tsv(paste0(outFolder,"Summary.FDR.",Sys.Date(),".tsv"))

res %>% write_tsv(paste0(outFolder,"ALL.combined.",Sys.Date(),".tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>0.0) %>%
    write_tsv(paste0(outFolder,"SIG.combined",Sys.Date(),".tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>1.0) %>% select(kbid,cname,log2FoldChange,padj) 



res %>% filter(padj<0.1,abs(log2FoldChange)>1) %>% dplyr::count(cname) %>%
    write_tsv(paste0(outFolder,"Summary.FDR.abslogfc_1",Sys.Date(),".tsv"))


