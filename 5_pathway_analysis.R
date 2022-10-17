###################################################
### pathway enrichment analysis

###################################################

library(tidyverse)
library(qqman)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(stringr)
library(magrittr)


#outFolder <- paste0("11_pathway_enrichment/")
outFolder <- paste0("11_pathway_enrichment_with_covidcontrol/")
system(paste0("mkdir -p ",outFolder))

######################################################################################################

pathway_enrich<-function(res_gene=res3,cname_select="CAM_T-cell_M",padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
    cat ("=============== ",cname_select,"=============== ", "\n")
    pathway_enrich_cname_dir<-paste0(outFolder,cname_select,"/")
    system(paste0("mkdir -p ",pathway_enrich_cname_dir))
    result<-list()
    aux <- res_gene %>% filter(cname==cname_select)
    genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    message(".................................")
    message("Number of DE genes: ",length(genes))
    #print(length(genes))
    
    message(".................................")
    message("enrichGO")
    ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize=5)
    print(head(ego))
    result$enrichGO<-ego
    #save(ego,file=paste0(pathway_enrich_cname_dir,"ego.RData"))
    #write.csv(ego,file=paste0(pathway_enrich_cname_dir,"ego.csv"))
    
    print(".................................")
    print("enrichKEGG")
    ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa",minGSSize=5)
    print(head(ekegg)) 
    result$enrichKEGG<-ekegg
    #save(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.RData"))
    #write.csv(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.csv"))
    
    
    message(".................................")
    message("enrichPathway")
    erpath <- enrichPathway(gene=genes,universe=geneUniv,minGSSize=5)
    print(head(erpath))
    result$enrichPathway<-erpath
    #save(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.RData"))
    # write.csv(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.csv"))
    # 
    # message(".................................")
    # message("gseGO")
    # # BP: biological_process, CC: cellular_component, MF: molecular_function
    # gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP",minGSSize=5)
    # print(head(gseGO.res))
    # result$gseGO<-gseGO.res
    # #save(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.RData"))
    # #write.csv(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.csv"))
    
    
    # message(".................................")
    # message("gsePathway")
    # gseRPath.res <- gsePathway(geneList)
    # print(head(gseRPath.res))
    # result$gsePathway<-gseRPath.res
    return (result)
}


# load DE genes
#res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_library/ALL.combined.2022-03-06.tsv")
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv")

# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","Cluster","Origin"),sep="_",remove=FALSE) #Cell_type

# for now
#res$Cell_type<-res$Cluster

#cell.type.annotation<-read.delim("cell.type.annotation.txt")
cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)

res$Cell_type<-clust2Names[res$Cluster]
res$cname<-paste0(res$Cell_type,"_",res$Location,"_",res$Origin)



DEGs_CAM<-res %>% filter(padj<0.1 & Location=="CAM") %>% select(gene_name) %>% unlist() %>% unique() 
DEGs_PVBP<-res %>% filter(padj<0.1 & Location=="PVBP") %>% select(gene_name) %>% unlist() %>% unique()

length(intersect(DEGs_CAM,DEGs_PVBP))



# Removing na pvalues
# Grouping pvalues based on the Location,Cell_type,and Origin
# Adding a column showing the rank of each pvalue devided by the number of pvalues in each group 


#ENTREZID id 
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

# non-na ENTREZID were included
res3 <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

#exploratory analysis
# # counting the number of DE genes per cname   
DE_per_cname<-sapply(unique(res3$cname), function(x,padj_cutoff=0.1,log2FoldChange_cutoff=0.5){
    aux <- res3 %>% filter(cname==x)
    genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    length(genes)
    
})

cname_selected<-names(DE_per_cname)[which(DE_per_cname>0)]
result_pathway_en_list<-lapply(cname_selected, function(x) return(pathway_enrich(res3,x)))
names(result_pathway_en_list)<-cname_selected

save(result_pathway_en_list,file=paste0(outFolder,"pathwayEnrich_result.RData"))
which(DE_per_cname>0)


##########################################################################################  
#####                                  dot plot
##########################################################################################



######################################################
# Pathway enrichment analysis (ORA)
# based on specific cell types 
######################################################

#load(paste0("11_pathway_enrichment_batch_library_corrected/pathwayEnrich_result.RData"))
load(paste0(outFolder,"pathwayEnrich_result.RData"))
cname_selected<-names(result_pathway_en_list)
cname_selected_test<-as.data.frame(cname_selected)
cname_selected_test <- cname_selected_test %>% separate(cname_selected,c("Location","Cluster","Origin"),sep="_",remove=FALSE) #Cell_type
cname_selected_test$Cell_type<-clust2Names[cname_selected_test$Cluster]
cname_selected_test$cname<-paste0(cname_selected_test$Cell_type,"_",cname_selected_test$Location,"_",cname_selected_test$Origin)

names(result_pathway_en_list)<-cname_selected_test$cname
cname_selected<-cname_selected_test$cname
# DE_per_cname_select<-names(DE_per_cname)[which(DE_per_cname>=5)]
# result_pathway_en_list<-result_pathway_en_list [DE_per_cname_select]
#cname_selected<-names(result_pathway_en_list)
#gseGO
res_gseGO_list<-lapply(cname_selected, function(x)
{
    
    rs<-result_pathway_en_list[[x]]$gseGO@result %>% filter(qvalues<=0.05)
    dim1<-dim(rs)[1]
    
    if(min(dim1,5)>0)
    {
        res_en<-rs
        
        # to calculate GeneRatio=count/setSize
        
        #count
        gene_count<- res_en %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
        
        ## merge with the original dataframe
        dot_df<- left_join(res_en, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
        
        dot_df<-dot_df[1:min(dim1,5),c("ID","Description" ,"enrichmentScore","p.adjust","GeneRatio")]
        dot_df$cname<-rep(x,min(dim1,5))
        dot_df
    }
})


res_df_gseGO <- do.call(rbind,res_gseGO_list)

#res_df_gseGO<-res_df_gseGO[1:15,]
pdf(paste0(outFolder,"gseGO_cname_DotPlot.pdf"),width=25,height=30)
ggplot(res_df_gseGO, 
       aes(x = cname, y = Description)) + 
    geom_point(aes(size = enrichmentScore, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
    ylab(NULL)+ 
    xlab(NULL) +
    theme_bw()+
    theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()


#enrichGO
res_enrichGO_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichGO)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichGO@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,5)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,5),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,5))
        }
        res_en }
}
)  

res_df_enrichGO <- do.call(rbind,res_enrichGO_list)

res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichGO<-res_df_enrichGO %>% filter(p.adjust<0.1)

mt<-matrix(nrow=length(unique(res_df_enrichGO$Description)),ncol=length(unique(res_df_enrichGO$cname)),0)
rownames(mt)<-unique(res_df_enrichGO$Description)
colnames(mt)<-unique(res_df_enrichGO$cname)
for ( i in unique(res_df_enrichGO$Description))
{
    print(i)
    inx<-which(res_df_enrichGO$Description==i)
    mt[i,res_df_enrichGO$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]
res_df_enrichGO$orderpathways<-orderpathways[res_df_enrichGO$Description]
res_df_enrichGO <- res_df_enrichGO  %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
res_df_enrichGO$cname <- factor(res_df_enrichGO$cname, levels = unique(res_df_enrichGO$cname))
#res_df_enrichGO$Description <- factor(res_df_enrichGO$Description, levels = unique(res_df_enrichGO$Description))

pdf(paste0(outFolder,"enrichGO_cname_DotPlot.pdf"),width=40,height=30)
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cname, y = reorder(Description,orderpathways))) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+
    xlab(NULL)+
    coord_fixed(ratio = 1)+
    theme_bw()+
    theme(axis.text=element_text(size=25),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=25)) 
dev.off()



#enrichKEGG
res_enrichKEGG_list<-lapply(cname_selected, function(x)
{
    if(length(result_pathway_en_list[[x]]$enrichKEGG)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichKEGG@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,5)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,5),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,5))
        }
        res_en }
}
)  

res_df_enrichKEGG <- do.call(rbind,res_enrichKEGG_list)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})



res_df_enrichKEGG<-res_df_enrichKEGG %>% filter(p.adjust<0.1) 

mt<-matrix(nrow=length(unique(res_df_enrichKEGG$Description)),ncol=length(unique(res_df_enrichKEGG$cname)),0)
rownames(mt)<-unique(res_df_enrichKEGG$Description)
colnames(mt)<-unique(res_df_enrichKEGG$cname)


for ( i in unique(res_df_enrichKEGG$Description))
{
    inx<-which(res_df_enrichKEGG$Description==i)
    mt[i,res_df_enrichKEGG$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichKEGG$orderpathways<-orderpathways[res_df_enrichKEGG$Description]

res_df_enrichKEGG <- res_df_enrichKEGG %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
res_df_enrichKEGG$cname <- factor(res_df_enrichKEGG$cname, levels = unique(res_df_enrichKEGG$cname))



pdf(paste0(outFolder,"enrichKEGG_cname_DotPlot.pdf"),width=20,height=25)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cname, y = reorder(Description,orderpathways))) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    #theme_black()+
    theme_bw()+
    theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()



#enrichPathway

res_enrichPathway_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichPathway)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichPathway@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,5)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,5),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,5))
        }
        res_en }
}
)  

res_df_enrichPathway <- do.call(rbind,res_enrichPathway_list)

res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})


res_df_enrichPathway<-res_df_enrichPathway %>% filter(p.adjust<0.1) 

mt<-matrix(nrow=length(unique(res_df_enrichPathway$Description)),ncol=length(unique(res_df_enrichPathway$cname)),0)
rownames(mt)<-unique(res_df_enrichPathway$Description)
colnames(mt)<-unique(res_df_enrichPathway$cname)


for ( i in unique(res_df_enrichPathway$Description))
{
    inx<-which(res_df_enrichPathway$Description==i)
    mt[i,res_df_enrichPathway$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichPathway$orderpathways<-orderpathways[res_df_enrichPathway$Description]

res_df_enrichPathway <- res_df_enrichPathway %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
res_df_enrichPathway$cname <- factor(res_df_enrichPathway$cname, levels = unique(res_df_enrichPathway$cname))



pdf(paste0(outFolder,"enrichPathway_cname_DotPlot.pdf"),width=30,height=35)
ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cname, y = reorder(Description,orderpathways))) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+
    xlab(NULL)+
    coord_fixed(ratio = 1)+
    theme_bw()+
    theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()

############################################################
### Pathway enrichment analysis (ORA)
### DE genes combined
############################################################


#genes <- filter(res3,padj<0.1,abs(log2FoldChange)>0.5) %>% dplyr::select(ENTREZID) %>% unlist %>% unique
genes <- filter(res3,padj<0.1) %>% dplyr::select(ENTREZID) %>% unlist %>% unique

geneUniv <- res3 %>% dplyr::select(ENTREZID) %>% unlist %>% unique

geneList <- -log10(res3$pvalue)
names(geneList) <- res3$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

message(".................................")
message("enrichGO")
ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
print(head(ego))

res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)

res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})


mt<-matrix(nrow=length(unique(res_df_enrichGO$Description)),ncol=length(unique(res_df_enrichGO$GeneRatio)),0)
rownames(mt)<-unique(res_df_enrichGO$Description)
colnames(mt)<-unique(res_df_enrichGO$GeneRatio)
for ( i in unique(res_df_enrichGO$Description))
{
    print(i)
    inx<-which(res_df_enrichGO$Description==i)
    mt[i,res_df_enrichGO$GeneRatio[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]
res_df_enrichGO$orderpathways<-orderpathways[res_df_enrichGO$Description]
res_df_enrichGO <- res_df_enrichGO  %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
#res_df_enrichGO$GeneRatio <- factor(res_df_enrichGO$GeneRatio, levels = unique(res_df_enrichGO$GeneRatio))



res_df_enrichGO<-res_df_enrichGO[1:min(nrow(res_df_enrichGO),30),]
pdf(paste0(outFolder,"enrichGO_combined_DotPlot.pdf"),width=18,height=20)
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = reorder(Description,orderpathways))) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()



print(".................................")
print("enrichKEGG")
ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
print(head(ekegg)) 

res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.1)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichKEGG<-res_df_enrichKEGG[1:30,]
pdf(paste0(outFolder,"enrichKEGG.combined_DotPlot.pdf"),width=20,height=15)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    theme_bw()+
    theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()

message(".................................")
message("enrichPathway")
erpath <- enrichPathway(gene=genes,universe=geneUniv)
print(head(erpath))

res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.1)

res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})


res_df_enrichPathway<-res_df_enrichPathway[1:30,]
pdf(paste0(outFolder,"enrichPathway.combined_DotPlot.pdf"),width=30,height=25)
ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    theme_bw()+
    theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()






############################################################
### Pathway enrichment analysis (ORA)
### DE genes combined - Location based 
############################################################


sapply (unique(res3$Location), function(x)
    {
    loc<-x
    #genes <- filter(res3,padj<0.1,abs(log2FoldChange)>0.5) %>% dplyr::select(ENTREZID) %>% unlist %>% unique
    genes <- filter(res3,padj<0.1 & Location==loc) %>% dplyr::select(ENTREZID) %>% unlist %>% unique
    
    geneUniv <- filter(res3,Location==loc) %>% dplyr::select(ENTREZID) %>% unlist %>% unique
    
    geneList <- -log10(res3$pvalue)
    names(geneList) <- res3$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    
    message(".................................")
    message("enrichGO")
    ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
    print(head(ego))
    
    res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
    
    res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
        numden<-unlist(strsplit(x,"/"))
        return (as.numeric(numden[1])/as.numeric(numden[2]))
    })
    
    
    mt<-matrix(nrow=length(unique(res_df_enrichGO$Description)),ncol=length(unique(res_df_enrichGO$GeneRatio)),0)
    rownames(mt)<-unique(res_df_enrichGO$Description)
    colnames(mt)<-unique(res_df_enrichGO$GeneRatio)
    for ( i in unique(res_df_enrichGO$Description))
    {
        print(i)
        inx<-which(res_df_enrichGO$Description==i)
        mt[i,res_df_enrichGO$GeneRatio[inx]]<-1
    }
    
    orderpathways<-rowSums(mt)
    orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]
    res_df_enrichGO$orderpathways<-orderpathways[res_df_enrichGO$Description]
    res_df_enrichGO <- res_df_enrichGO  %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
    #res_df_enrichGO$GeneRatio <- factor(res_df_enrichGO$GeneRatio, levels = unique(res_df_enrichGO$GeneRatio))
    
    
    res_df_enrichGO<-res_df_enrichGO %>% filter (p.adjust<0.1)
    
    res_df_enrichGO<-res_df_enrichGO[1:min(nrow(res_df_enrichGO),30),]
    pdf(paste0(outFolder,"enrichGO_combined_",loc,"_DotPlot.pdf"),width=18,height=20)
    ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
           aes(x = GeneRatio, y = reorder(Description,orderpathways))) + 
        geom_point(aes(size = GeneRatio, color = p.adjust)) +
        theme_bw(base_size = 14) +
        #scale_colour_gradient(limits=c(0, 0.10), low="red") +
        scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
        theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
        labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
        ylab(NULL)+ 
        xlab(NULL)+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    dev.off()
    
    
    
    print(".................................")
    print("enrichKEGG")
    ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
    print(head(ekegg)) 
    
    res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.1)
    
    res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
        numden<-unlist(strsplit(x,"/"))
        return (as.numeric(numden[1])/as.numeric(numden[2]))
    })
    
    res_df_enrichKEGG<-res_df_enrichKEGG[1:30,]
    pdf(paste0(outFolder,"enrichKEGG.combined_",loc,"_DotPlot.pdf"),width=15,height=8)
    ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
           aes(x = GeneRatio, y = Description)) + 
        geom_point(aes(size = GeneRatio, color = p.adjust)) +
        theme_bw(base_size = 14) +
        scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
        theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
        labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
        ylab(NULL)+ 
        xlab(NULL)+
        theme_bw()+
        theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    dev.off()
    
    message(".................................")
    message("enrichPathway")
    erpath <- enrichPathway(gene=genes,universe=geneUniv)
    print(head(erpath))
    
    res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.1)
    
    res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
        numden<-unlist(strsplit(x,"/"))
        return (as.numeric(numden[1])/as.numeric(numden[2]))
    })
    
    
    res_df_enrichPathway<-res_df_enrichPathway[1:30,]
    pdf(paste0(outFolder,"enrichPathway.combined_",loc,"_DotPlot.pdf"),width=20,height=20)
    ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
           aes(x = GeneRatio, y = Description)) + 
        geom_point(aes(size = GeneRatio, color = p.adjust)) +
        theme_bw(base_size = 11) +
        scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
        theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
        labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
        ylab(NULL)+ 
        xlab(NULL)+
        theme_bw()+
        theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    dev.off()
    
})




# 
# 
# ################################################## 
# # WikiPathways
# ################################################## 
# 
# 
# res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")
# res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
# ## cluster colors 
# 
# clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# 
# names(clust2Names)<-as.character(c(0:23))
# cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
#                   "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48",
#                   "#AA5485","#4E747A","#C59A89","#C9C76F")   
# names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                          "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                          "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# res$Cell_type<-clust2Names[res$Cell_type]
# 
# #ENTREZID
# eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# names(eg)[1]="gene_name"
# head(eg)
# e2g <- eg$gene_name
# names(e2g) <- eg$ENTREZID
# res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))
# 
# geneList <- -log10(res$pvalue)
# names(geneList) <- res$ENTREZID
# geneList = sort(geneList, decreasing = TRUE)
# 
# 
# res$cname<-paste0(res$Cell_type,"_",res$Origin)
# 
# ########################################################
# # Pathway enrichment analysis (ORA)- cell type specific 
# # WikiPathways
# #qvalueCutoff  = 0.05
# ########################################################
# 
# pathway_enrich<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
# {
#     
#     print(cname_select)
#     
#     aux <- res %>% filter(cname==cname_select)
#     
#     if(nrow(aux)>0)
#     {
#         gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
#         geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
#         ##geneList <- aux$log2FoldChange
#         geneList <- -log10(aux$pvalue)
#         names(geneList) <- aux$ENTREZID
#         geneList = sort(geneList, decreasing = TRUE)
#         #length(genes)
#         
#         wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")
#         
#         wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
#         wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
#         wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
#         
#         
#         print(length(gene))
#         
#         if (length(gene)>0)
#         {
#             
#             ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
#             
#             if (!is.null(ewp))
#             {
#                 res_en<-ewp@result
#                 res_en<-res_en%>%filter(p.adjust<0.1)
#                 dim1<-nrow(res_en)
#                 #print(dim1)
#                 if(min(dim1,15)>0)
#                 {
#                     
#                     
#                     res_en<-res_en[1:min(dim1,15),c("ID","Description" ,"GeneRatio","p.adjust","geneID")]
#                     res_en$Cell_type<-rep(cname_select,min(dim1,15))
#                 }
#                 return (res_en)
#             }
#             
#         }
#         
#     }
#     
# }
# 
# res_enrichwikilist<-lapply(unique(res$cname), function(x) {pathway_enrich(cname_select=x)})  
# res_df_enrichwiki <- do.call(rbind,res_enrichwikilist)
# 
# res_df_enrichwiki$GeneRatio<-sapply(res_df_enrichwiki$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# # WP289: Myometrial Relaxation and Contraction Pathways
# genes_WP289<-unlist(strsplit(res_df_enrichwiki$geneID[which(res_df_enrichwiki$ID=="WP289")],"/"))
# 
# 
# # ORA
# 
# res_df<-res_df_enrichwiki
# pdf(paste0(outFolder,"enrich_wikipathways_cname_DotPlot.pdf"),width=35,height=52)
# ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = Cell_type, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     #theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     coord_fixed(ratio = 1)+
#     #theme_black()+
#     theme_bw()+
#     #theme(axis.text.x = element_text(angle = 45))+
#     theme(axis.text.y = element_text(hjust = 1))+
#     theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# 
# #ggtitle("GO pathway enrichment")
# dev.off()
# 
# 
# ########################################################
# # ORA -all DEGs combined
# # WikiPathways
# ########################################################
# 
# 
# gene<-genes <- filter(res,padj<0.1,abs(log2FoldChange)>0) %>% dplyr::select(ENTREZID) %>% unlist
# geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist
# ##geneList <- aux$log2FoldChange
# geneList <- -log10(res$pvalue)
# names(geneList) <- res$ENTREZID
# geneList = sort(geneList, decreasing = TRUE)
# #length(genes)
# 
# wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")
# 
# wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
# wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
# wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
# 
# ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
# res_df<-ewp@result
# res_df<-res_df%>%filter(p.adjust<0.1)
# 
# res_df$GeneRatio<-sapply(res_df$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df<-res_df[1:15,]
# 
# pdf(paste0(outFolder,"enrich_wikipathways.combined_DotPlot.pdf"),width=12,height=10)
# ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 11) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# ###########################################################
# # Pathway enrichment analysis on subtypes
# ###########################################################
# 
# # load DE genes
# res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")
# 
# # Adding location, cell type, and origin columns 
# res <- res %>% separate(cname,c("Cluster","Origin"),sep="_",remove=FALSE) #Cell_type
# 
# 
# clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
# names(clust2Names)<-c(0:23)
# res$Cell_type<-clust2Names[res$Cluster]
# res$cname<-paste0(res$Cell_type,"_",res$Origin)  
# 
# 
# 
# #ENTREZID id 
# eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# names(eg)[1]="gene_name"
# head(eg)
# 
# e2g <- eg$gene_name
# names(e2g) <- eg$ENTREZID
# 
# # non-na ENTREZID were included
# res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))
# 
# ######################################################################################################################
# # Pathway analysis based on shared DEGs and DEGs specific to each cell type
# # macrophage-1 and macrophage-2 and macrophage-3
# ######################################################################################################################
# 
# res3<-res %>% filter(padj<0.1)
# res_intersect_macrophages<-res3 %>% filter(Cell_type %in% c("Macrophage-3","Macrophage-2","Macrophage-1"))
# 
# Macrophage1<-unique(res3$gene_name[which(res3$Cell_type=="Macrophage-1")])
# Macrophage2<-unique(res3$gene_name[which(res3$Cell_type=="Macrophage-2")])
# Macrophage3<-unique(res3$gene_name[which(res3$Cell_type=="Macrophage-3")])
# x<-list(Macrophage1,Macrophage2,Macrophage3)
# library(VennDiagram)
# 
# library(VennDiagram)
# venn.plot <- draw.triple.venn(
#     area1 = length(Macrophage1),
#     area2 = length(Macrophage2),
#     area3 = length(Macrophage3),
#     n12 = length(intersect(Macrophage1,Macrophage2)),
#     n23 = length(intersect(Macrophage2,Macrophage3)),
#     n13 = length(intersect(Macrophage1,Macrophage3)),
#     n123 = length(intersect(Macrophage1,intersect(Macrophage2,Macrophage3))),
#     category = c("Macrophage-1","Macrophage-2","Macrophage-3"),
#     fill = c(cluster.Colors["Macrophage-1"], cluster.Colors["Macrophage-2"], cluster.Colors["Macrophage-3"]),
#     scaled=FALSE)
# 
# png(filename = paste0(outFolder,"venn_cluster.Macrophages.png"))
# grid.draw(venn.plot)
# dev.off()
# 
# 
# ###########################################################
# ###### shared DEGs (macrophage-1 and macrophage-2 and macrophage-3)
# ###########################################################
# 
# genes <- filter(res_intersect_macrophages,padj<0.1) %>% dplyr::select(ENTREZID) %>% unlist %>% unique
# geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique
# 
# 
# system(paste0("mkdir -p ",outFolder,"Shared_macrophages1-3/"))
# 
# 
# message(".................................")
# message("enrichGO")
# ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize = 5)
# print(head(ego))
# 
# 
# res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
# res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichGO<-res_df_enrichGO %>%filter(p.adjust<0.1)
# 
# res_df_enrichGO<-res_df_enrichGO[1:50,]
# pdf(paste0(outFolder,"Shared_macrophages1-3/","enrichGO_commonDEG_DotPlot.pdf"),width=20,height=30)
# ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# 
# 
# ###########################################################
# print(".................................")
# print("enrichKEGG")
# ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa",minGSSize = 5)
# print(head(ekegg)) 
# 
# 
# res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# res_df_enrichKEGG<-res_df_enrichKEGG %>%filter(p.adjust<0.1)
# res_df_enrichKEGG<-res_df_enrichKEGG[1:50,]
# pdf(paste0(outFolder,"Shared_macrophages1-3/","enrichKEGG_commonDEG_DotPlot.pdf"),width=20,height=30)
# ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# 
# ###########################################################
# message(".................................")
# message("enrichPathway")
# erpath <- enrichPathway(gene=genes,universe=geneUniv,minGSSize = 5)
# print(head(erpath))
# 
# res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichPathway<-res_df_enrichPathway %>%filter(p.adjust<0.1)
# res_df_enrichPathway<-res_df_enrichPathway[1:50,]
# 
# pdf(paste0(outFolder,"Shared_macrophages1-3/","enrichPathway_commonDEG_DotPlot.pdf"),width=20,height=30)
# ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 11) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# ################################################################################################################################################
# # pathway enrichment analysis based on specific cell markers (macrophage-1 and macrophage-2 and macrophage-3)
# ################################################################################################################################################
# 
# 
# ########################################################################
# # DEGs specific to Macrophage-1
# ########################################################################
# 
# genes<-Macrophage1[which(!Macrophage1 %in% c(Macrophage2,Macrophage3))]
# genes<-names(e2g)[which(e2g %in%genes)]
# 
# geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique
# 
# 
# system(paste0("mkdir -p ",outFolder,"Macrophages1-only/"))
# 
# 
# message(".................................")
# message("enrichGO")
# ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize = 5)
# print(head(ego))
# 
# 
# res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
# res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichGO<-res_df_enrichGO %>%filter(p.adjust<0.1)
# 
# res_df_enrichGO<-res_df_enrichGO[1:50,]
# pdf(paste0(outFolder,"Macrophages1-only/","enrichGO_DotPlot.pdf"),width=20,height=30)
# ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# 
# 
# print(".................................")
# print("enrichKEGG")
# ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa",minGSSize = 5)
# print(head(ekegg)) 
# 
# 
# res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# res_df_enrichKEGG<-res_df_enrichKEGG %>%filter(p.adjust<0.1)
# 
# pdf(paste0(outFolder,"Macrophages1-only/","enrichKEGG_DotPlot.pdf"),width=15,height=20)
# ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# 
# message(".................................")
# message("enrichPathway")
# erpath <- enrichPathway(gene=genes,universe=geneUniv,minGSSize = 5)
# print(head(erpath))
# 
# res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichPathway<-res_df_enrichPathway %>%filter(p.adjust<0.1)
# res_df_enrichPathway<-res_df_enrichPathway[1:50,]
# 
# pdf(paste0(outFolder,"Macrophages1-only/","enrichPathway_DotPlot.pdf"),width=20,height=30)
# ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 11) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# ########################################################################
# 
# # DEGs specific to Macrophage-2
# ########################################################################
# 
# genes<-Macrophage2[which(!Macrophage2 %in% c(Macrophage1,Macrophage3))]
# genes<-names(e2g)[which(e2g %in%genes)]
# 
# geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique
# 
# 
# system(paste0("mkdir -p ",outFolder,"Macrophages2-only/"))
# 
# 
# message(".................................")
# message("enrichGO")
# ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize = 5)
# print(head(ego))
# 
# 
# res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
# res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichGO<-res_df_enrichGO %>%filter(p.adjust<0.1)
# 
# 
# pdf(paste0(outFolder,"Macrophages2-only/","enrichGO_DotPlot.pdf"),width=18,height=10)
# ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# 
# print(".................................")
# print("enrichKEGG")
# ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa",minGSSize = 5)
# print(head(ekegg)) 
# 
# 
# res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# res_df_enrichKEGG<-res_df_enrichKEGG %>%filter(p.adjust<0.1)
# 
# pdf(paste0(outFolder,"Macrophages2-only/","enrichKEGG_DotPlot.pdf"),width=15,height=15)
# ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# 
# message(".................................")
# message("enrichPathway")
# erpath <- enrichPathway(gene=genes,universe=geneUniv,minGSSize = 5)
# print(head(erpath))
# 
# res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# # res_df_enrichPathway<-res_df_enrichPathway %>%filter(p.adjust<0.1)
# # res_df_enrichPathway<-res_df_enrichPathway[1:50,]
# # 
# # pdf(paste0(outFolder,"Macrophages2-only/","enrichPathway_DotPlot.pdf"),width=20,height=30)
# # ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
# #        aes(x = GeneRatio, y = Description)) + 
# #     geom_point(aes(size = GeneRatio, color = p.adjust)) +
# #     theme_bw(base_size = 11) +
# #     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
# #     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
# #     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
# #     ylab(NULL)+ 
# #     xlab(NULL)+
# #     theme_bw()+
# #     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# # dev.off()
# 
# 
# ############################################################
# # Pathway enrichment analysis 
# 
# # DEGs specific to Macrophage-3
# ########################################################################
# 
# genes<-Macrophage3[which(!Macrophage3 %in% c(Macrophage1,Macrophage2))]
# genes<-names(e2g)[which(e2g %in%genes)]
# 
# geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique
# 
# 
# system(paste0("mkdir -p ",outFolder,"Macrophages3-only/"))
# 
# 
# message(".................................")
# message("enrichGO")
# ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize = 5)
# print(head(ego))
# 
# 
# res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
# res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichGO<-res_df_enrichGO %>%filter(p.adjust<0.1)
# 
# res_df_enrichGO<-res_df_enrichGO[1:50,]
# pdf(paste0(outFolder,"Macrophages3-only/","enrichGO_DotPlot.pdf"),width=18,height=30)
# ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# 
# print(".................................")
# print("enrichKEGG")
# ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa",minGSSize = 5)
# print(head(ekegg)) 
# 
# 
# res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# res_df_enrichKEGG<-res_df_enrichKEGG %>%filter(p.adjust<0.1)
# 
# pdf(paste0(outFolder,"Macrophages3-only/","enrichKEGG_DotPlot.pdf"),width=18,height=10)
# ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) + 
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+ 
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
# dev.off()
# 
# 
# 
# message(".................................")
# message("enrichPathway")
# erpath <- enrichPathway(gene=genes,universe=geneUniv,minGSSize = 5)
# print(head(erpath))
# 
# res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichPathway<-res_df_enrichPathway %>%filter(p.adjust<0.1)
# # res_df_enrichPathway<-res_df_enrichPathway[1:50,]
# # 
# pdf(paste0(outFolder,"Macrophages3-only/","enrichPathway_DotPlot.pdf"),width=20,height=15)
# ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 11) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# 
# ############################################################
# # pathway enrichment analysis 
# # Stromal cell markers 
# ###########################################################
# 
# 
# allmarkers<-read.csv(file=paste0("5_harmonyClustersDGE/ClusterDEG.csv"),stringsAsFactors = FALSE)
# allmarkers<-allmarkers %>%filter(p_val_adj<0.05)
# allmarkers$Cell_type<-clust2Names[as.character(allmarkers$cluster)]
# 
# 
# 
# outFolder <- paste0("11_pathway_enrichment_batch_library_corrected/")
# subtype<-"Stromal1"
# subtype<-"Stromal2"
# subtype<-"Myofibroblast"
# system(paste0("mkdir -p ",outFolder, subtype,"/"))
# 
# 
# #m2 = read_tsv("5_harmonySubTypesDGE/Stromal/ClusterDEG.tsv")
# m2 = read_tsv("5_harmonySubTypes_v2_DGE/Stromal/ClusterDEG.tsv")
# m2<-m2 %>% filter(p_val_adj<0.1)
# 
# stromal1<-m2$symbol[which(m2$Celltype=="Stromal-1" )][1:100]
# stromal2<-m2$symbol[which(m2$Celltype=="Stromal-2" )][1:100]
# myofibroblast<-m2$symbol[which(m2$Celltype=="Myofibroblast" )][1:100]
# 
# 
# genes<-names(e2g)[which(e2g %in%myofibroblast)]
# if (subtype=="Stromal1") 
#     genes<-names(e2g)[which(e2g %in%stromal1)]
# if (subtype=="Stromal2") 
#     genes<-names(e2g)[which(e2g %in%stromal2)]
# 
# 
# 
# res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")
# res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
# ## cluster colors 
# 
# clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# 
# names(clust2Names)<-as.character(c(0:23))
# cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
#                   "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48",
#                   "#AA5485","#4E747A","#C59A89","#C9C76F")   
# names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                          "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                          "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# res$Cell_type<-clust2Names[res$Cell_type]
# 
# #ENTREZID
# eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# names(eg)[1]="gene_name"
# head(eg)
# e2g <- eg$gene_name
# names(e2g) <- eg$ENTREZID
# res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))
# 
# geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique
# #geneUniv <- unique(m2$symbol)
# #geneUniv<-names(e2g)[which(e2g %in%geneUniv)]
# 
# 
# 
# message(".................................")
# message("enrichGO")
# ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize = 5)
# print(head(ego))
# 
# 
# res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
# res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichGO<-res_df_enrichGO %>%filter(p.adjust<0.1)
# 
# res_df_enrichGO<-res_df_enrichGO[1:50,]
# pdf(paste0(outFolder, subtype,"/","enrichGO_DotPlot.pdf"),width=15,height=30)
# ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# 
# print(".................................")
# print("enrichKEGG")
# ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa",minGSSize = 5)
# print(head(ekegg)) 
# 
# 
# res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# res_df_enrichKEGG<-res_df_enrichKEGG %>%filter(p.adjust<0.1)
# 
# pdf(paste0(outFolder, subtype,"/","enrichKEGG_DotPlot.pdf"),width=20,height=8)
# ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# message(".................................")
# message("enrichPathway")
# erpath <- enrichPathway(gene=genes,universe=geneUniv,minGSSize = 5)
# print(head(erpath))
# 
# res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichPathway<-res_df_enrichPathway %>%filter(p.adjust<0.1)
# # res_df_enrichPathway<-res_df_enrichPathway[1:50,]
# #
# pdf(paste0(outFolder, subtype,"/","enrichPathway_DotPlot.pdf"),width=30,height=22)
# ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 11) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# 
# ################################################################################################################################################
# # pathway enrichment analysis based on specific cell markers (smooth muscle cells )
# # top 100 cell markers
# ################################################################################################################################################
# # SMC-1, SMC-2, SMC-3 
# 
# res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")
# res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
# ## cluster colors 
# 
# clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# 
# names(clust2Names)<-as.character(c(0:23))
# cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
#                   "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
#                   "#AA5485","#4E747A","#C59A89","#C9C76F")   
# names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                          "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                          "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# res$Cell_type<-clust2Names[res$Cell_type]
# 
# #ENTREZID
# eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# names(eg)[1]="gene_name"
# head(eg)
# e2g <- eg$gene_name
# names(e2g) <- eg$ENTREZID
# res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))
# 
# 
# allmarkers<-read.csv(file=paste0("5_harmonyClustersDGE/ClusterDEG.csv"),stringsAsFactors = FALSE)
# allmarkers<-allmarkers %>%filter(p_val_adj<0.05)
# allmarkers$Cell_type<-clust2Names[as.character(allmarkers$cluster)]
# 
# 
# outFolder <- paste0("11_pathway_enrichment_batch_library_corrected/")
# subtype<-"SMC1"
# subtype<-"SMC2"
# subtype<-"SMC3"
# system(paste0("mkdir -p ",outFolder, subtype,"/"))
# 
# 
# #m2 = read_tsv("5_harmonySubTypesDGE/SMC/ClusterDEG.tsv")
# m2 = read_tsv("5_harmonySubTypes_v2_DGE/SMC/ClusterDEG.tsv")
# m2<-m2 %>% filter(p_val_adj<0.1)
# 
# smc1<-m2$symbol[which(m2$Celltype=="Smooth muscle cells-1" )][1:100]
# smc2<-m2$symbol[which(m2$Celltype=="Smooth muscle cells-2" )][1:100]
# smc3<-m2$symbol[which(m2$Celltype=="Smooth muscle cells-3" )][1:100]
# 
# 
# genes<-names(e2g)[which(e2g %in%smc3)]
# if (subtype=="SMC1") 
#     genes<-names(e2g)[which(e2g %in%smc1)]
# if (subtype=="SMC2") 
#     genes<-names(e2g)[which(e2g %in%smc2)]
# 
# 
# # genes<-names(e2g)[which(e2g %in%stromal2)]
# # genes<-names(e2g)[which(e2g %in%myofibroblast)]
# 
# 
# geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique
# #geneUniv <- unique(m2$symbol)
# #geneUniv<-names(e2g)[which(e2g %in%geneUniv)]
# 
# 
# message(".................................")
# message("enrichGO")
# ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize = 5)
# print(head(ego))
# 
# 
# res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
# res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichGO<-res_df_enrichGO %>%filter(p.adjust<0.1)
# 
# res_df_enrichGO<-res_df_enrichGO[1:50,]
# pdf(paste0(outFolder, subtype,"/","enrichGO_DotPlot.pdf"),width=18,height=30)
# ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# 
# print(".................................")
# print("enrichKEGG")
# ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa",minGSSize = 5)
# print(head(ekegg)) 
# 
# 
# res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# res_df_enrichKEGG<-res_df_enrichKEGG %>%filter(p.adjust<0.1)
# 
# pdf(paste0(outFolder, subtype,"/","enrichKEGG_DotPlot.pdf"),width=20,height=25)
# ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# 
# message(".................................")
# message("enrichPathway")
# erpath <- enrichPathway(gene=genes,universe=geneUniv,minGSSize = 5)
# print(head(erpath))
# 
# res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichPathway<-res_df_enrichPathway %>%filter(p.adjust<0.1)
# # res_df_enrichPathway<-res_df_enrichPathway[1:50,]
# #
# pdf(paste0(outFolder, subtype,"/","enrichPathway_DotPlot.pdf"),width=25,height=25)
# ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 11) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# 
# ############################################################
# # pathway enrichment analysis 
# # Macrophages cell markers
# ###########################################################
# 
# res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")
# res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
# ## cluster colors 
# 
# clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# 
# names(clust2Names)<-as.character(c(0:23))
# cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
#                   "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
#                   "#AA5485","#4E747A","#C59A89","#C9C76F")   
# names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                          "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                          "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# res$Cell_type<-clust2Names[res$Cell_type]
# 
# #ENTREZID
# eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# names(eg)[1]="gene_name"
# head(eg)
# e2g <- eg$gene_name
# names(e2g) <- eg$ENTREZID
# res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))
# 
# 
# allmarkers<-read.csv(file=paste0("5_harmonyClustersDGE/ClusterDEG.csv"),stringsAsFactors = FALSE)
# allmarkers<-allmarkers %>%filter(p_val_adj<0.05)
# allmarkers$Cell_type<-clust2Names[as.character(allmarkers$cluster)]
# 
# outFolder <- paste0("11_pathway_enrichment_batch_library_corrected/")
# subtype<-"Macrophage1-cellmarkers"
# subtype<-"Macrophage2-cellmarkers"
# subtype<-"Macrophage3-cellmarkers"
# subtype<-"Macrophage4-cellmarkers"
# system(paste0("mkdir -p ",outFolder, subtype,"/"))
# 
# 
# #m2 = read_tsv("5_harmonySubTypesDGE/Macrophage/ClusterDEG.tsv")
# m2 = read_tsv("5_harmonySubTypes_v2_DGE/Macrophage/ClusterDEG.tsv")
# m2<-m2 %>% filter(p_val_adj<0.1)
# 
# mac1<-m2$symbol[which(m2$Celltype=="Macrophage-1" )][1:100]
# mac2<-m2$symbol[which(m2$Celltype=="Macrophage-2" )][1:100]
# mac3<-m2$symbol[which(m2$Celltype=="Macrophage-3" )][1:100]
# mac4<-m2$symbol[which(m2$Celltype=="Macrophage-4" )][1:100]
# 
# genes<-names(e2g)[which(e2g %in%mac4)]
# if (subtype=="Macrophage1-cellmarkers") 
#     genes<-names(e2g)[which(e2g %in%mac1)]
# if (subtype=="Macrophage2-cellmarkers") 
#     genes<-names(e2g)[which(e2g %in%mac2)]
# if (subtype=="Macrophage3-cellmarkers") 
#     genes<-names(e2g)[which(e2g %in%mac3)]
# 
# 
# geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique
# 
# 
# message(".................................")
# message("enrichGO")
# ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize = 5)
# print(head(ego))
# 
# 
# res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
# res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichGO<-res_df_enrichGO %>%filter(p.adjust<0.1)
# 
# res_df_enrichGO<-res_df_enrichGO[1:50,]
# pdf(paste0(outFolder, subtype,"/","enrichGO_DotPlot.pdf"),width=20,height=20)
# ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# 
# print(".................................")
# print("enrichKEGG")
# ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa",minGSSize = 5)
# print(head(ekegg)) 
# 
# 
# res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# res_df_enrichKEGG<-res_df_enrichKEGG %>%filter(p.adjust<0.1)
# 
# pdf(paste0(outFolder, subtype,"/","enrichKEGG_DotPlot.pdf"),width=15,height=20)
# ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# 
# message(".................................")
# message("enrichPathway")
# erpath <- enrichPathway(gene=genes,universe=geneUniv,minGSSize = 5)
# print(head(erpath))
# 
# res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.1)
# 
# res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichPathway<-res_df_enrichPathway %>%filter(p.adjust<0.1)
#  res_df_enrichPathway<-res_df_enrichPathway[1:50,]
# #
# pdf(paste0(outFolder, subtype,"/","enrichPathway_DotPlot.pdf"),width=15,height=7)
# ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 11) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()
# 
# 
# 
# 
# ########################################################################################
# # Pathway enrichment analysis
# ########################################################################################
# 
# 
# # SMC-1 
# # pathway enrichment analysis DEGs in SMC-1
# 
# res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")
# res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
# ## cluster colors 
# 
# clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# 
# names(clust2Names)<-as.character(c(0:23))
# cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
#                   "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
#                   "#AA5485","#4E747A","#C59A89","#C9C76F")   
# names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
#                          "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
#                          "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
# res$Cell_type<-clust2Names[res$Cell_type]
# 
# #ENTREZID
# eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# names(eg)[1]="gene_name"
# head(eg)
# e2g <- eg$gene_name
# names(e2g) <- eg$ENTREZID
# res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))
# 
# 
# genes<-res %>% filter (padj <0.1 & Cell_type=="Smooth muscle cells-1") %>% dplyr::select (ENTREZID) %>% unlist %>% unique()
# geneUniv<-res %>% filter (Cell_type=="Smooth muscle cells-1") %>% dplyr::select (ENTREZID) %>% unlist %>% unique()
# 
# 
# message(".................................")
# message("enrichGO")
# ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize = 5)
# print(head(ego))
# 
# 
# res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
# res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# 
# res_df_enrichGO<-res_df_enrichGO %>%filter(p.adjust<0.1)
# 
# 
# pdf(paste0(outFolder, "SMC1_enrichGO_DotPlot.pdf"),width=20,height=20)
# ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 14) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30))
# dev.off()



