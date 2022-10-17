################################################################################################################################################
# pathway enrichment analysis based on specific cell markers
################################################################################################################################################

library(tidyverse)
library(qqman)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(stringr)
library(magrittr)



outFolder <- paste0("6_pathway_enrichment_markers/")
system(paste0("mkdir -p ",outFolder))



allmarkers<-read.csv(file=paste0("5_harmonyClustersDGE/ClusterDEG.csv"),stringsAsFactors = FALSE)



#ENTREZID
eg = bitr(unique(allmarkers$symbol), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

geneUniv<-names(e2g)
allmarkers<-allmarkers %>%filter(p_val_adj<0.05)
allmarkers$cluster<-as.character(allmarkers$cluster)
allmarkers$Cell_type<-allmarkers$cluster#   clust2Names[as.character(allmarkers$cluster)]
m2<-allmarkers


res_enrichPathway_list<-lapply(1:length(unique(allmarkers$cluster)),function(i){
    subtype<-unique(allmarkers$cluster)[i]
    print(subtype)
    y<-m2$symbol[which(m2$cluster==subtype )]
    l=min(length(y),100)
    y<-y[1:l]
    genes<-names(e2g)[which(e2g %in%y)]
    message(".................................")
    message("enrichPathway")
    ego <- enrichPathway(gene=genes,universe=geneUniv,minGSSize = 5)
    print(head(ego))
    res_df_enrich<-ego@result %>% filter(qvalue<=0.1)
    res_df_enrich$GeneRatio<-sapply(res_df_enrich$GeneRatio, function(x){
        numden<-unlist(strsplit(x,"/"))
        return (as.numeric(numden[1])/as.numeric(numden[2]))
    })
    
    
    dot_df<- res_df_enrich<-res_df_enrich %>%filter(p.adjust<0.1)
    if(nrow(dot_df)>0)
    {
        dot_df<-dot_df[1:min(nrow(dot_df),5),c("ID","Description","geneID","p.adjust","GeneRatio")]
        dot_df$cluster<-rep(subtype,min(nrow(dot_df),5))
    }
    
    dot_df
})

res_df_enrichPathway <- do.call(rbind,res_enrichPathway_list)
res_df_enrichPathway<-res_df_enrichPathway %>% filter(p.adjust<0.1) 

mt<-matrix(nrow=length(unique(res_df_enrichPathway$Description)),ncol=length(unique(res_df_enrichPathway$cluster)),0)
rownames(mt)<-unique(res_df_enrichPathway$Description)
colnames(mt)<-unique(res_df_enrichPathway$cluster)


for ( i in unique(res_df_enrichPathway$Description))
{
    inx<-which(res_df_enrichPathway$Description==i)
    mt[i,res_df_enrichPathway$cluster[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichPathway$orderpathways<-orderpathways[res_df_enrichPathway$Description]

res_df_enrichPathway <- res_df_enrichPathway %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
res_df_enrichPathway$cluster <- factor(res_df_enrichPathway$cluster, levels = unique(res_df_enrichPathway$cluster))


pdf(paste0(outFolder,"enrichPathway_DotPlot.pdf"),height=35,width=30)#width=55,height=40) #height=10,width=35)#width=50,height=40)#)#, 
ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cluster, y = reorder(Description,orderpathways))) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.5), low="red") +
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
    #theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.2)) +
    theme(text = element_text(size=30)) 
dev.off()





res_enrichKEGG_list<-lapply(1:length(unique(allmarkers$cluster)),function(i){
    subtype<-unique(allmarkers$cluster)[i]
    y<-m2$symbol[which(m2$cluster==subtype )]
    l=min(length(y),100)
    y<-y[1:l] 
    genes<-names(e2g)[which(e2g %in%y)]
    message(".................................")
    message("enrichKEGG")
    ego <- enrichKEGG(gene=genes,universe=geneUniv,minGSSize = 5,organism="hsa")
    print(head(ego))
    res_df_enrich<-ego@result %>% filter(qvalue<=0.1)
    res_df_enrich$GeneRatio<-sapply(res_df_enrich$GeneRatio, function(x){
        numden<-unlist(strsplit(x,"/"))
        return (as.numeric(numden[1])/as.numeric(numden[2]))
    })
    
    
    dot_df<- res_df_enrich<-res_df_enrich %>%filter(p.adjust<0.1)
    if(nrow(dot_df)>0)
    {
        dot_df<-dot_df[1:min(nrow(dot_df),5),c("ID","Description","geneID","p.adjust","GeneRatio")]
        dot_df$cluster<-rep(subtype,min(nrow(dot_df),5))
    }
    
    dot_df
})

res_df_enrichKEGG <- do.call(rbind,res_enrichKEGG_list)
res_df_enrichKEGG<-res_df_enrichKEGG %>% filter(p.adjust<0.1) 

mt<-matrix(nrow=length(unique(res_df_enrichKEGG$Description)),ncol=length(unique(res_df_enrichKEGG$cluster)),0)
rownames(mt)<-unique(res_df_enrichKEGG$Description)
colnames(mt)<-unique(res_df_enrichKEGG$cluster)


for ( i in unique(res_df_enrichKEGG$Description))
{
    inx<-which(res_df_enrichKEGG$Description==i)
    mt[i,res_df_enrichKEGG$cluster[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichKEGG$orderpathways<-orderpathways[res_df_enrichKEGG$Description]

res_df_enrichKEGG <- res_df_enrichKEGG %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
res_df_enrichKEGG$cluster <- factor(res_df_enrichKEGG$cluster, levels = unique(res_df_enrichKEGG$cluster))

   

pdf(paste0(outFolder,"enrichKEGG_DotPlot.pdf"),height=25,width=20)#width=55,height=40) #height=10,width=35)#width=50,height=40)#)#, 
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cluster, y = reorder(Description,orderpathways))) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.5), low="red") +
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
    #theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.2)) +
    theme(text = element_text(size=30)) 
dev.off()



   
res_enrichGO_list<-lapply(1:length(unique(allmarkers$cluster)),function(i){
    subtype<-unique(allmarkers$cluster)[i]
    print(subtype)
    y<-m2$symbol[which(m2$cluster==subtype )]
    l<-min(100,length(y))
    y<-y[1:l]
    genes<-names(e2g)[which(e2g %in%y)]
    message(".................................")
    message("enrichGO")
    ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP",minGSSize = 5)
    print(head(ego))
    res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)
    res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
        numden<-unlist(strsplit(x,"/"))
        return (as.numeric(numden[1])/as.numeric(numden[2]))
    })
    dot_df<- res_df_enrichGO<-res_df_enrichGO %>%filter(p.adjust<0.1)
    
    
    
    if(nrow(dot_df)>0)
    {
        
        dot_df<-dot_df[1:min(nrow(dot_df),5),c("ID","Description","geneID","p.adjust","GeneRatio")]
        dot_df$cluster<-rep(subtype,min(nrow(dot_df),5))
        
    }
    
    dot_df
})


res_df_enrichGO <- do.call(rbind,res_enrichGO_list)
res_df_enrichGO<-res_df_enrichGO %>% filter(p.adjust<0.1)


mt<-matrix(nrow=length(unique(res_df_enrichGO$Description)),ncol=length(unique(res_df_enrichGO$cluster)),0)
rownames(mt)<-unique(res_df_enrichGO$Description)
colnames(mt)<-unique(res_df_enrichGO$cluster)
for ( i in unique(res_df_enrichGO$Description))
{
    print(i)
    inx<-which(res_df_enrichGO$Description==i)
    mt[i,res_df_enrichGO$cluster[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]
res_df_enrichGO$orderpathways<-orderpathways[res_df_enrichGO$Description]
res_df_enrichGO <- res_df_enrichGO  %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
res_df_enrichGO$cluster <- factor(res_df_enrichGO$cluster, levels = unique(res_df_enrichGO$cluster))



pdf(paste0(outFolder,"enrichGO_DotPlot.pdf"),height=25,width=30)#width=55,height=40) #height=10,width=35)#width=50,height=40)#)#, 
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cluster, y = reorder(Description,orderpathways))) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.5), low="red") +
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
    #theme(text = element_text(size=30)) +
    theme(axis.text.y = element_text(hjust = 1))+
    theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.2)) +
    theme(text = element_text(size=30)) 
#ggtitle("GO pathway enrichment")
dev.off()
