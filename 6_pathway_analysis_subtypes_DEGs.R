###################################################
### pathway enrichment analysis

###################################################
library(tidyverse)
library(qqman)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(stringr)


outFolder <- paste0("6_pathway_enrichment_compartments/")
system(paste0("mkdir -p ",outFolder))

######################################################################################################





# pathway_enrich<-function(res_gene=res3,cname_select="CAM_T-cell_M",padj_cutoff=0.1,log2FoldChange_cutoff=0)
# {
#     cat ("=============== ",cname_select,"=============== ", "\n")
#     pathway_enrich_cname_dir<-paste0(outFolder,cname_select,"/")
#     system(paste0("mkdir -p ",pathway_enrich_cname_dir))
#     result<-list()
#     aux <- res_gene %>% filter(cname==cname_select)
#     genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
#     geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
#     ##geneList <- aux$log2FoldChange
#     geneList <- -log10(aux$pvalue)
#     names(geneList) <- aux$ENTREZID
#     geneList = sort(geneList, decreasing = TRUE)
#     message(".................................")
#     message("Number of DE genes: ",length(genes))
#     #print(length(genes))
#     
#     message(".................................")
#     message("enrichGO")
#     ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Mm.eg.db,ont="BP",minGSSize=5)
#     print(head(ego))
#     result$enrichGO<-ego
#     #save(ego,file=paste0(pathway_enrich_cname_dir,"ego.RData"))
#     #write.csv(ego,file=paste0(pathway_enrich_cname_dir,"ego.csv"))
#     
#     print(".................................")
#     print("enrichKEGG")
#     ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="mmu",minGSSize=5)
#     print(head(ekegg)) 
#     result$enrichKEGG<-ekegg
#     #save(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.RData"))
#     #write.csv(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.csv"))
#     
#     
#     message(".................................")
#     message("enrichPathway")
#     erpath <- enrichPathway(gene=genes,universe=geneUniv,organism="mouse",minGSSize=5)
#     print(head(erpath))
#     result$enrichPathway<-erpath
#     #save(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.RData"))
#     # write.csv(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.csv"))
#     
#     message(".................................")
#     message("gseGO")
#     # BP: biological_process, CC: cellular_component, MF: molecular_function
#     gseGO.res <- gseGO(geneList,  OrgDb=org.Mm.eg.db,ont="BP",minGSSize=5)
#     print(head(gseGO.res))
#     result$gseGO<-gseGO.res
#     #save(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.RData"))
#     #write.csv(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.csv"))
#     
#     
#     message(".................................")
#     message("gsePathway")
#     gseRPath.res <- gsePathway(geneList,organism="mouse",minGSSize=5)
#     print(head(gseRPath.res))
#     result$gsePathway<-gseRPath.res
#     return (result)
# }
# 
# 
# # load DE genes
# res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
# 
# # Adding location, cell type, and origin columns 
# #res <- res %>% separate(cname,c("Cluster","Origin"),sep="_",remove=FALSE) #Cell_type
# #res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
# res <- res %>% separate(cname,c("Location","Cluster"),sep="_",remove=FALSE)
# 
cell.type.annotation<-read.delim("cell-type.annotation.txt")

clust2Names<-cell.type.annotation$Cell.type #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster #c(0:23)
# res$Cell_type<-clust2Names[res$Cluster]
# res$cname<-paste0(res$Cell_type,"_",res$Location)

# 
# 
# 
# 
# 
# 
# # Removing na pvalues
# # Grouping pvalues based on the Location,Cell_type,and Origin
# # Adding a column showing the rank of each pvalue devided by the number of pvalues in each group 
# # res2 <- res %>% filter(!is.na(pvalue)) %>%
# #     arrange(pvalue) %>%
# #     group_by(Cell_type) %>%
# #     mutate(r=rank(pvalue, ties.method = "random"),pexp=r/length(pvalue))
# 
# 
# #ENTREZID id 
# eg = bitr(res$kbid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
# names(eg)[1]="gene_name"
# head(eg)
# 
# e2g <- eg$gene_name
# names(e2g) <- eg$ENTREZID
# 
# colnames(res)[1]<-"gene_name"
# 
# res3 <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))
# 
# 
# 
# #exploratory analysis
# # # counting the number of DE genes per cname   
# DE_per_cname<-sapply(unique(res3$cname), function(x,padj_cutoff=0.1,log2FoldChange_cutoff=0.5){
#     aux <- res3 %>% filter(cname==x)
#     genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
#     geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
#     geneList <- -log10(aux$pvalue)
#     names(geneList) <- aux$ENTREZID
#     geneList = sort(geneList, decreasing = TRUE)
#     length(genes)
#     
# })
# 
# cname_selected<-names(DE_per_cname)[which(DE_per_cname>5)]
# result_pathway_en_list<-lapply(cname_selected, function(x) return(pathway_enrich(res3,x)))
# names(result_pathway_en_list)<-cname_selected
# 
# save(result_pathway_en_list,file=paste0(outFolder,"pathwayEnrich_result.RData"))
# which(DE_per_cname>0)


##########################################################################################  
#####                                  dot plot
##########################################################################################

cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster #c(0:23)
load(paste0("6_pathway_enrichment/pathwayEnrich_result.RData"))
cname_selected<-names(result_pathway_en_list)

cluster_numbers<-unlist(strsplit(cname_selected,"_"))[seq(1,3*length(cname_selected),by=3)]
cluster_numbers_df<-as.data.frame(cbind(cluster_numbers,cname=cname_selected))
cluster_numbers_df<-cluster_numbers_df %>% group_by (cluster_numbers) %>% mutate(counts=n()) %>% ungroup()



cluster_numbers_df %>% filter(counts>=2 ) %>% select(cluster_numbers) %>% unlist() %>% unique()

clusters1<-as.character(c(3,9))
clusters2<-as.character(c(14,23)) #c(13,14,28)
clusters3<-as.character(c(10,13))
clusters4<-as.character(c(1,21))
clusters5<-as.character(c(0,12))
clusters6<-as.character(c(2,4))
clusters7<-as.character(c(5,7,8))

clusters8<-as.character(c(8,10,14))
clusters9<-as.character(c(17))
clusters10<-as.character(c(8,10,14,28,20,5,7,8,23,11))
clusters11<-as.character(c(6))

clustersset<-list(clusters1,clusters2,clusters3,clusters4,clusters5,clusters6,clusters7,clusters8,clusters9,clusters10,clusters11)

for (i in 1:length(clustersset))
{
  print(clust2Names[clustersset[[i]]])
}



sapply(1:length(clustersset), function(x)
{
  
  
  clusters<-clustersset[[x]]
  subFolder<-paste0(clust2Names[clusters],collapse = "_")
  
  print(clust2Names[clusters])
  outFolder <- paste0("6_pathway_enrichment_compartments/")
  
  subFolder<-gsub("\\(", "",subFolder)
  subFolder<-gsub("\\)", "",subFolder)
  subFolder<-gsub("\\ ", "-",subFolder)
  
  
  system(paste0("mkdir -p ",outFolder,subFolder,"/"))
  thisFolder <- paste0(outFolder,subFolder,"/")
  
  cname_selected<-cluster_numbers_df %>% filter(cluster_numbers %in% clusters) %>% select(cname) %>% unlist()
  
  
  
  #enrichGO
  res_enrichGO_list<-lapply(cname_selected, function(x)
  {
    
    
    if(length(result_pathway_en_list[[x]]$enrichGO)>0)
    {
      rs<-result_pathway_en_list[[x]]$enrichGO@result %>% filter(qvalue<=0.05)
      
      dim1<-dim(rs)[1]
      res_en<-NULL
      if(min(dim1,10)>0)
      {
        res_en<-rs
        res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
        res_en$cname<-rep(x,min(dim1,10))
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
  #res_df_enrichGO<-res_df_enrichGO[1:15,]
  
  
  mt<-matrix(nrow=length(unique(res_df_enrichGO$Description)),ncol=length(unique(res_df_enrichGO$cname)),0)
  rownames(mt)<-unique(res_df_enrichGO$Description)
  colnames(mt)<-unique(res_df_enrichGO$cname)
  
  
  for ( i in unique(res_df_enrichGO$Description))
  {
    inx<-which(res_df_enrichGO$Description==i)
    mt[i,res_df_enrichGO$cname[inx]]<-1
  }
  
  orderpathways<-rowSums(mt)
  orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]
  
  res_df_enrichGO$orderpathways<-orderpathways[res_df_enrichGO$Description]
  
  
  res_df_enrichGO$Location<-sapply(res_df_enrichGO$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
  })
  
  
  res_df_enrichGO <- res_df_enrichGO %>% dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
  
  res_df_enrichGO$cname <- factor(res_df_enrichGO$cname, levels = unique(res_df_enrichGO$cname))
  
  res_df_enrichGO<-res_df_enrichGO %>% arrange(desc(Location))
  res_df_enrichGO$cname <- factor(res_df_enrichGO$cname, levels = unique(res_df_enrichGO$cname))
  
  
  #subFolder<-"allEpithelials"
  pdf(paste0(thisFolder,subFolder,"_","enrichGO.pdf"),width=20,height=12)
  ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
         aes(x = cname, y = reorder(Description,orderpathways))) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    #theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+
    xlab(NULL)+
    coord_fixed(ratio = 1)+
    #theme_black()+
    theme_bw()+
    #facet_grid(.~ Location)+
    #theme(axis.text.x = element_text(angle = 45))+
    theme(text = element_text(size=30)) +
    theme(axis.text.y = element_text(hjust = 1))+
    theme(axis.text.x = element_text(angle = 90, hjust=1)) 
  
  #ggtitle("GO pathway enrichment")
  dev.off()
  
  
  
  ####################
  
  #enrichKEGG
  res_enrichKEGG_list<-lapply(cname_selected, function(x)
  {
    if(length(result_pathway_en_list[[x]]$enrichKEGG)>0)
    {
      rs<-result_pathway_en_list[[x]]$enrichKEGG@result %>% filter(qvalue<=0.05)
      
      dim1<-dim(rs)[1]
      res_en<-NULL
      if(min(dim1,10)>0)
      {
        res_en<-rs
        res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
        res_en$cname<-rep(x,min(dim1,10))
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
  
  res_df_enrichKEGG$Location<-sapply(res_df_enrichKEGG$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
  })
  
  
  
  res_df_enrichKEGG <- res_df_enrichKEGG %>% dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
  res_df_enrichKEGG$cname <- factor(res_df_enrichKEGG$cname, levels = unique(res_df_enrichKEGG$cname))
  
  
  
  #res_df_enrichKEGG<-res_df_enrichKEGG[1:15,]
  pdf(paste0(thisFolder,subFolder,"_","enrichKEGG.pdf"),width=20,height=18)
  ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
         aes(x = cname, y = reorder(Description,orderpathways))) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    #theme_black()+
    theme_bw()+
    #facet_grid(.~Location)+
    #theme(axis.text.x = element_text(angle = 45))+
    theme(axis.text.y = element_text(hjust = 1))+
    #theme(text = element_text(size=40)) +
    theme(text = element_text(size=30)) +
    #theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    theme(axis.text.x = element_text(angle = 90, hjust=1)) 
  dev.off()
  
  
  ####################
  
  #enrichPathway
  
  res_enrichPathway_list<-lapply(cname_selected, function(x)
  {
    
    
    if(length(result_pathway_en_list[[x]]$enrichPathway)>0)
    {
      rs<-result_pathway_en_list[[x]]$enrichPathway@result %>% filter(qvalue<=0.05)
      
      dim1<-dim(rs)[1]
      res_en<-NULL
      if(min(dim1,10)>0)
      {
        res_en<-rs
        res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
        res_en$cname<-rep(x,min(dim1,10))
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
  
  
  
  res_df_enrichPathway$Location<-sapply(res_df_enrichPathway$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
  })
  
  #res_df_enrichPathway<-res_df_enrichPathway[1:15,]
  
  res_df_enrichPathway <- res_df_enrichPathway %>% dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
  
  res_df_enrichPathway$cname <- factor(res_df_enrichPathway$cname, levels = unique(res_df_enrichPathway$cname))
  
  
  pdf(paste0(thisFolder,subFolder,"_","enrichPathway.pdf"),width=25,height=25)
  ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest 
         aes(x = cname, y = reorder(Description,orderpathways))) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    #aes(x = cname, y = Description)) +
    #geom_point(aes(size = GeneRatio, color = p.adjust)) +
    
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+
    xlab(NULL)+
    coord_fixed(ratio = 1)+
    #theme_black()+
    theme_bw()+
    #facet_grid(.~Location)+
    theme(text = element_text(size=30)) +
    #theme(axis.text.x = element_text(angle = 45))+
    theme(axis.text.y = element_text(hjust = 1))+
    #theme(text = element_text(size=40)) +
    #theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    theme(axis.text.x = element_text(angle = 90,  hjust=1)) 
  #scale_x_discrete(expand = c(0,1.9))+
  #scale_y_discrete(expand = c(0,20))
  #ggtitle("GO pathway enrichment")
  dev.off()
  
})