library(tidyverse)
library(dplyr)
library(qqman)
library(org.Hs.eg.db)
library(clusterProfiler)
library(reshape2)
library(ggplot2)


##################################################################
# comparison between between single cell and bulk datasets: 

##################################################################


outFolder="12_comparison_with_bulk/"
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
res <-res %>% filter(!is.na(log2FoldChange) & !is.na(padj))


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

res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))


# bulk datasets
#datasets<-list.files("Bulk_data/  ")
datasets<-c("PTL_IAI-PTL_ALLList.csv","PTL_IAI-PTL_SIAI_ALLList.csv","PTL_SIAI-PTL_ALLList.csv","TLvsTNL_blood_ENTREZ.csv")
names(datasets)<-c("PTL_IAI-PTL","PTL_IAI-PTL_SIAI","PTL_SIAI-PTL","blood")



#############################################
# selecting DEGs from single cell or bulk, or choosing all genes

# set the  threshold and dirfilter parameters:
#threshold<-"Bulk" # selecting DEGs from bulk data
#threshold<-"All" # all genes
threshold<-"single"
#dirfilter<-"Bulk_padj0.1/"
#dirfilter<-"All/"
if(threshold=="single") dirfilter<-"Single_padj0.1/"

for (i in 1:length(datasets))
{
  
  bulk<-names(datasets)[i]
  
  system(paste0("mkdir -p ", outFolder,dirfilter,"_",bulk,"/"))
  
  
  dataset<-datasets[i]
  ref_data<-read_csv(paste0("Bulk_data/",dataset))
  
  if (bulk=="blood")
  {
    
   
    ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ ,Rt=t)
    colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
    ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
    ref_data<-ref_data %>%filter(!is.na(ENTREZID) & !is.na(R.Log2FC) & !is.na(Rpvalue))
    }
  else
  {
    ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC=log2FoldChange,P.Value,adj.P.Val,ENTREZ ,RlfcSE=lfcSE)
    colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","RlfcSE")
    ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
    ref_data<-ref_data %>%filter(!is.na(ENTREZID) & !is.na(R.Log2FC) & !is.na(Rpvalue))
    ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
    ref_data$Rt<-ref_data$R.Log2FC/ref_data$RlfcSE
    }
  
 
  
  cor_locations<-lapply(unique(res$Location), function(loc)
    {
   
    subFolder=loc
    system(paste0("mkdir -p ", outFolder,dirfilter,"/","_",bulk,"/",subFolder,"/"))
    finalFolder<-paste0(outFolder,dirfilter,"_",bulk,"/",subFolder,"/")
    
    cor_celltypes_list<-lapply (unique(res$cname), function(cl ){
      
      res_filter<-res %>% filter(cname ==cl & Location==loc)
      resJoin<-res_filter %>% inner_join(ref_data)
      #resJoin<-resJoin %>% filter(padj< 0.1 )
      #cat(table(resJoin$ref_color),sep=":")
      resJoin$t<-resJoin$log2FoldChange/resJoin$lfcSE
      
      
      
      resJoin_singlecell<-resJoin %>% filter(padj<0.1)
      resJoin_bulk<-resJoin %>% filter(Rpadj<0.1)
      resJoin_all<-resJoin 
      
      #res_filter2<-res %>% filter(padj<0.1)
      #res_filter2<-res_filter2 %>% filter(Cell_type %in% cl)
      # resJoin2<-res_filter2 %>% inner_join(ref_data)
      # #resJoin<-resJoin %>% filter(padj< 0.1 )
      # #cat(table(resJoin$ref_color),sep=":")
      # resJoin2$t<-resJoin2$log2FoldChange/resJoin2$lfcSE
      
      gtitle<-paste0(cl,"-",tolower(loc), " vs ",tolower(bulk))
      
      spearman_all<-c()
      spearman_cor_pval<-c()
      
      if(nrow(resJoin_bulk)>=10 | nrow(resJoin_singlecell)>=10)
      {
        
        
        if(nrow(resJoin_bulk)>=10 & threshold!="single")
        {
          spearman_all<-cor.test(resJoin_bulk$t,resJoin_bulk$Rt,method="spearman",na.rm=TRUE)
        }
          
          
        
        if (nrow(resJoin_singlecell)>=10 & threshold=="single" )
        {
          spearman_all<-cor.test(resJoin_singlecell$t,resJoin_singlecell$Rt,method="spearman",na.rm=TRUE)
        }
          
        if (length(spearman_all)>0)
        {
          spearman_cor_pval<-c(spearman_all$estimate,spearman_all$p.value)
          names(spearman_cor_pval)<-c("spearman_cor","pvalue")
          #gtitle<-paste0(cl,"-",tolower(loc), " vs bulk ",tolower(bulk)," (cor= ",round(spearman_cor_pval[1],2),", p= ",formatC(spearman_cor_pval[2], format = "e", digits = 2),")")
          
        }
        
        if (threshold=="All"  & nrow(resJoin)>=10 )
        {
          spearman_all<-cor.test(resJoin_all$t,resJoin_all$Rt,method="spearman",na.rm=TRUE)
        }
          
        }
      
      res_rest<-res %>% filter(padj<0.1)
      
      
      
      
      
      resJoin$ref_color <- "None"          
      #light blue
      resJoin$ref_color[resJoin$padj<0.1  &resJoin$Rpadj>0.1 ] <- "Only single cell"  
      #purpule  //res_rest union of all DEGs across all cell types
      resJoin$ref_color[(!resJoin$ENTREZID %in% unique(res_rest$ENTREZID)) &  resJoin$Rpadj<0.1  ] <- "Only bulk" 
      #blue
      resJoin$ref_color[resJoin$padj<0.1 & resJoin$Rpadj<0.1]="Single cell and bulk" 
      
      # // resJoin$ref_color<-as.fac
      resJoin$ref_color <- factor(resJoin$ref_color,levels=unique(resJoin$ref_color))
      
      
      table(resJoin$ref_color)
      
      cat(table(resJoin$ref_color),sep=":")
      #resJoin$t<-resJoin$log2FoldChange/resJoin$lfcSE
      if (nrow(resJoin)>0)
      {
        p2 <- resJoin %>% arrange(ref_color) %>%
          ggplot(aes(Rt,t,color=ref_color)) +
          geom_point()+ #aes(colour = ref_color)) +
          #geom_smooth(method=lm, se=FALSE,linetype = "dashed", color="black")+  # filtered way
          xlab(paste0(bulk," (bulk)"," standardized Log2FC"))+
          ylab("Standardized log2FC")+
          scale_color_manual(name="Differentially expressed",values=c("None"="#CCCCCC","Only single cell"="#BDD7EE","Only bulk"="#DD99DD","Single cell and bulk"="#0000EE"))+
          ggtitle(gtitle)+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5))
        
        cl<-gsub("\\(", "",cl)
        cl<-gsub("\\)", "",cl)
        cl<-gsub("\\ ", "-",cl)
        
        
        fname=paste0(finalFolder,paste0(cl,".sc.",loc,"_bulk.",bulk,".png"))
        ggsave(fname,p2,width=9,height=7)
        
      }
        
      if (length(spearman_cor_pval)==0)
        spearman_cor_pval<-c(NA,NA)
      return(spearman_cor_pval)
      
    }) 
    
    cor_celltypes<-do.call(rbind,cor_celltypes_list)
    
    #cor_celltypes<-t(cor_celltypes)
    colnames(cor_celltypes)<-c("spearman_cor","spearman_pvalue")
    #cor_celltypes<-cor_celltypes %>% arrange(desc(spearman_cor,-log(spearman_pvalue)))
    
    cor_celltypes<-as.data.frame(cor_celltypes)
    
    cor_celltypes$celltype<-unique(res$cname)
    cor_celltypes$location<-loc
    
    cor_celltypes<-cor_celltypes %>% filter(!is.na(spearman_cor))
    cor_celltypes<-cor_celltypes %>% arrange(desc(spearman_cor))
    cor_celltypes$spearman_padj<-p.adjust(cor_celltypes$spearman_pvalue,"fdr")
    cor_celltypes2<-cor_celltypes %>% dplyr::select(celltype,spearman_cor,spearman_pvalue,spearman_padj)
    
    fname=paste0(finalFolder,paste0("cor.sc.",loc,"_bulk.",bulk,".csv"))
    write.csv(cor_celltypes2,file=fname)
    
    
    cor_celltypes$sig<-sapply(cor_celltypes$spearman_padj, function(x){
      if (x<=0.0001) x="****"
      else if (x<=0.001)x="***" 
      else if (x<=0.01) x="**"
      else if (x<0.1) x="*"
      else if (x>=0.1) return ("ns")
    })
     
   
    
    #cor_celltypes<- cor_celltypes%>% filter(spearman_cor!=0)
    
    # fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".v.pdf"))
    # pdf(fname,width=10,height=4.5)
    # p2<-ggplot(data=cor_celltypes, aes(x=celltype, y=spearman_cor,fill=celltype)) +
    #   geom_bar(stat="identity",position="stack")+
    #   #geom_text(aes(label=sig,vjust = -sign(spearman_cor)), vjust=1.6, color="black", size=3.5)+
    #   theme_bw()+
    #   scale_fill_manual(values=cluster.Colors) +
    #   theme(axis.text.x = element_text(angle = 45, hjust=1))+
    #   theme(legend.position="none")+
    #   xlab("")+
    #   ylab("Spearman correlation")
    # p2
    # dev.off() 
    # 
    # fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".v.png"))
    # ggsave(fname,p2,width=10,height=4)
    
    
    
    # fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".v.pdf"))
    # pdf(fname,width=10,height=4)
    p2<-ggplot(data=cor_celltypes, aes(x=celltype, y=spearman_cor,fill=celltype)) +
      geom_bar(stat="identity",position="stack")+
      geom_text(aes(label=sig,vjust = -sign(spearman_cor)), vjust=1.6, color="black", size=3.5)+
      theme_bw()+
      scale_fill_manual(values=cluster.Colors) +
      theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position="none")+
      #theme(legend.position="none")+
      xlab("")+
      ylab("Spearman correlation")
    p2
    # dev.off()
    w=10
    if (length(unique(cor_celltypes$celltype))<4) w=5
    fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".withstars.v.png"))
    ggsave(fname,p2,width=w,height=4)
    
    
    # fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".h.pdf"))
    # pdf(fname,width=12,height=6)
    # p2<-ggplot(cor_celltypes, aes(x=reorder(celltype,-spearman_cor), y=spearman_cor,fill=celltype))+#,color=mycolor)+#,fill=DE) +
    #   geom_bar(stat='identity') +
    #   geom_text(aes(label=sig,hjust = -sign(spearman_cor)), vjust=0.7, color="black", size=3.5)+
    #   theme_bw()+
    #   scale_fill_manual(values=cluster.Colors) +
    #   #facet_grid(.~Location,scales="free") + 
    #   coord_flip() +
    #   theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position="none")+
    #   #theme(legend.title=element_blank())+
    #   ylab("Spearman correlation")+
    #   xlab("")
    # p2
    # dev.off() 
    
    # fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_.",bulk,".h.png"))
    # ggsave(fname,p2,width=8,height=10)
    # 
    
    # fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_.",bulk,".withstars.h.pdf"))
    # pdf(fname,width=8,height=10)
    # p2<-ggplot(cor_celltypes, aes(x=reorder(celltype,-spearman_cor), y=spearman_cor,fill=celltype))+#,color=mycolor)+#,fill=DE) +
    #   geom_bar(stat='identity') +
    #   geom_text(aes(label=sig,hjust = -sign(spearman_cor)), vjust=0.7, color="black", size=3.5)+
    #   theme_bw()+
    #   theme(legend.position="none")+
    #   scale_fill_manual(values=cluster.Colors) +
    #   #facet_grid(.~Location,scales="free") + 
    #   coord_flip() +
    #   #theme(legend.title=element_blank())+
    #   ylab("Spearman correlation")+
    #   xlab("")
    # p2
    # dev.off() 
    
    # fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_.",bulk,".withstars.h.png"))
    # ggsave(fname,p2,width=8,height=10)
    # 
    
    return(cor_celltypes)
    
  })
  
  
  total_cor_locations<-do.call(rbind,cor_locations)
  finalFolder<-paste0(outFolder,dirfilter,"/","_",bulk,"/")
  fname=paste0(finalFolder,paste0("cor_all_with_bulk.",bulk,".csv"))
  #total_cor_locations<-total_cor_locations[,-1]
  total_cor_locations$bulk<-bulk
  #total_cor_locations<-total_cor_locations %>% arrange(desc(spearman_cor,-log(spearman_pvalue)))
  write.csv(total_cor_locations,file=fname)
}






# 
# p2 <- resJoin %>% arrange(-padj) %>%
#   ggplot(aes(R.Log2FC,log2FoldChange,color=ref_color)) +
#   #ggplot(aes(Rt,t,color=ref_color)) +
#   geom_point()+ #aes(colour = ref_color)) +
#   geom_smooth(method=lm, se=FALSE,linetype = "dashed", color="black")+
#   xlab(paste0(bulk," (bulk)"," log2FC"))+
#   ylab("Log2FC")+
#   scale_color_manual(name="Differentially expressed",values=c("None"="#CCCCCC","Only single cell"="#BDD7EE","Only bulk"="#DD99DD","Single cell and bulk"="#0000EE"))+
#   ggtitle(gtitle)+
#   theme_bw()
# 
# 
# cl<-gsub("\\(", "",cl)
# cl<-gsub("\\)", "",cl)
# cl<-gsub("\\ ", "-",cl)
# 
# 
# fname=paste0(finalFolder,paste0("test_",cl,".sc.",loc,"_.",bulk,".png"))
# ggsave(fname,p2,width=9,height=7)
# 

# 
# i<-2
# dataset<-datasets[i]
# ref_data<-read_tsv(dataset)
# de<-ref_data
# 
# fname=paste0(outFolder,paste0(names(datasets)[i],"_hist_tstatistic.pdf"))
# pdf(fname,width=8,height=4)
# hist(ref_data$t)
# dev.off()
# 
# fname=paste0(outFolder,paste0(names(datasets)[i],"_hist_logfc.pdf"))
# pdf(fname,width=8,height=4)
# hist(ref_data$logFC)
# dev.off()
# 
# fname=paste0(outFolder,paste0(names(datasets)[i],"_hist_padj.pdf"))
# pdf(fname,width=8,height=4)
# hist(ref_data$adj.P.Val)
# dev.off()
# 
# 
# 
# # Convert directly in the aes()
# p <- ggplot(data=de, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point()
# # Add more simple "theme"
# p <- ggplot(data=de, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()
# # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
# p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
#   geom_hline(yintercept=-log10(0.05), col="red")
# # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
# 
# # add a column of NAs
# de$diffexpressed <- "NO"
# # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
# de$diffexpressed[de$logFC > 0.6 & de$adj.P.Val < 0.05] <- "UP"
# # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
# de$diffexpressed[de$logFC < -0.6 & de$adj.P.Val < 0.05] <- "DOWN"
# # Re-plot but this time color the points with "diffexpressed"
# p <- ggplot(data=de, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + geom_point() + theme_minimal()
# 
# # Add lines as before...
# p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
#   geom_hline(yintercept=-log10(0.05), col="red")
# ## Change point color 
# 
# # 1. by default, it is assigned to the categories in an alphabetical order):
# p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
# 
# # 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
# mycolors <- c("blue", "red", "black")
# names(mycolors) <- c("DOWN", "UP", "NO")
# p3 <- p2 + scale_colour_manual(values = mycolors)
# # Now write down the name of genes beside the points...
# # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
# #de$delabel <- NA
# #de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
# 
# 
# 
# pdf(paste0(outFolder,names(datasets)[i],"_volcano.pdf"),width=8,height=4)
# p<-ggplot(data=de, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + 
#   scale_colour_manual(values = mycolors)+
#   geom_point() + 
#   theme_minimal()
# p
# dev.off()
# 
# 
# pdf(paste0(outFolder,names(datasets)[i],"_standardized_volcano.pdf"),width=8,height=4)
# p<-ggplot(data=de, aes(x=t, y=-log10(adj.P.Val), col=diffexpressed)) + 
#   scale_colour_manual(values = mycolors)+
#   geom_point() + 
#   theme_minimal()
# p
# dev.off()
