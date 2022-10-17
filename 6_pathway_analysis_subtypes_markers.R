###################################################
### pathway enrichment analysis based on cluster markers

###################################################

library(tidyverse)
library(qqman)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(stringr)
library(magrittr)



#outFolder <- paste0("6_pathway_enrichment_subtypes_cellmarker/")

outFolder <- paste0("6_pathway_enrichment_comparison_subtypes_cellmarker/")
system(paste0("mkdir -p ",outFolder))


################################################################################################################################################
# pathway enrichment analysis based on specific cell markers 
################################################################################################################################################

cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
#cell.type.annotation$color[31]<-"#8B0000"
#write_tsv(cell.type.annotation,file="cell.type.annotation.v2.tsv")
#rownames(cell.type.annotation)<-NULL
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)


cluster.Colors<-cell.type.annotation$color
names(cluster.Colors)<-clust2Names


res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv")
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
res$cname<-paste0(res$Cell_type,"_",res$Location,"_",res$Origin)



############################################################
# pathway enrichment analysis 
# Stromal cell markers 
###########################################################

eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID


#####################################################################
# subtype i
#####################################################################

#m2 = read_tsv("5_harmonyClusters_withcovid19control_DGE/ClusterDEG.tsv")
m2 = read_tsv("5_harmonyClusters_withcovid19control_subtypes_DGE/ClusterDEG.tsv")
m2<-m2 %>%filter(p_val_adj<0.05)
m2 <- m2 %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)
m2$Celltype<-clust2Names[as.character(m2$cluster)]




pathway_enrich<-function(cluster_select)
{
    cat ("=============== ",cluster_select,"=============== ", "\n")
    pathway_enrich_cname_dir<-paste0(outFolder,cluster_select,"/")
    #system(paste0("mkdir -p ",pathway_enrich_cname_dir))
    result<-list()
    
    y<-m2$symbol[which(m2$Celltype==cluster_select )][1:100]
    genes<-names(e2g)[which(e2g %in%y)]
    
    #geneUniv<-m2 %>% filter(Celltype==cluster_select) %>% select (symbol %>% unlist %>% unique()
    geneUniv <- unique(m2$symbol)
    geneUniv<-names(e2g)[which(e2g %in%geneUniv)]
    
    
    #aux <- res_gene %>% filter(cname==cluster_select)
    
    #genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    #geneUniv <- m2 %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
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
    return (result)
}



subtypes<-unique(m2$Celltype)
result_pathway_en_list<-lapply(subtypes, function(x) return(pathway_enrich(x)))
names(result_pathway_en_list)<-subtypes

save(result_pathway_en_list,file=paste0(outFolder,"pathwayEnrich_result.RData"))




clusters1<-as.character(c(0,2,8,20,22,27)) #Stromal
clusters2<-as.character(c(4,6,26,29)) #Macrophage
clusters3<-as.character(c(5,9,15,32))#Tcell
clusters4<-as.character(c(10,12,17)) #Decidual
clusters5<-as.character(c(3,13,16,19,30)) #CTB

clustersset<-list(clusters1,clusters2,clusters3,clusters4,clusters5)
names(clustersset)<-c("Stromal","Macrophage","Tcell","Decidual","CTB")
for (i in 1:length(clustersset))
{
    print(clust2Names[clustersset[[i]]])
}


sapply(1:length(clustersset), function(x)
{
    
    
    clusters<-clustersset[[x]]
    subFolder<-paste0(names(clustersset)[x],"_",paste0(clusters,collapse = "_"))
    
    print(clust2Names[clusters])
   
   
    
    
    system(paste0("mkdir -p ",outFolder,subFolder,"/"))
    thisFolder <- paste0(outFolder,subFolder,"/")
    
    cname_selected<-clust2Names[clusters]
    
    
    
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
    
    
    # res_df_enrichGO$Location<-sapply(res_df_enrichGO$cname, function(x){
    #     x<-unlist(strsplit(x,"_"))
    #     return(x[length(x)])
    #     
    # })
    
    
    res_df_enrichGO <- res_df_enrichGO  %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
    
    res_df_enrichGO$cname <- factor(res_df_enrichGO$cname, levels = unique(res_df_enrichGO$cname))
    
    
    
    
    #subFolder<-"allEpithelials"
    pdf(paste0(thisFolder,subFolder,"_","enrichGO.pdf"),width=20,height=22)
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
   
    
    
    res_df_enrichKEGG <- res_df_enrichKEGG %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
    res_df_enrichKEGG$cname <- factor(res_df_enrichKEGG$cname, levels = unique(res_df_enrichKEGG$cname))
    
    
    
    #res_df_enrichKEGG<-res_df_enrichKEGG[1:15,]
    pdf(paste0(thisFolder,subFolder,"_","enrichKEGG.pdf"),width=14,height=14)
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
    
    
    
    
    
    #res_df_enrichPathway<-res_df_enrichPathway[1:15,]
    
    res_df_enrichPathway <- res_df_enrichPathway  %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
    
    res_df_enrichPathway$cname <- factor(res_df_enrichPathway$cname, levels = unique(res_df_enrichPathway$cname))
    
    
    pdf(paste0(thisFolder,subFolder,"_","enrichPathway.pdf"),width=32,height=20)
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
