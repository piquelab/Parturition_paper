######################################
### plot cell types ###
######################################

## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)

#################
##library(SingleR)

# cell.type.annotation<-read.delim("../../labor_myo/myometrium_analysis/data.colors.mouse.human.txt")
# cluster.Colors.rest<-cell.type.annotation$cluster.Colors_human2
# names(cluster.Colors.rest)<-cell.type.annotation$cluster_human
# cluster.Colors.rest<-cluster.Colors.rest[cluster.Colors.rest!=""]
# 
# 
# 
# cell.type.annotation<-read.delim("cell.type.annotation.txt")
# clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
# clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
# names(clust2Names)<-as.character(cell.type.annotation$Cluster)



sc<-read_rds("4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.final.rds")

 parturion_samples<-read_delim("parturition_cv.txt")
 parturion_samples_data<- parturion_samples %>% select(Pregnancy_ID,Condition=Labor)
#sc<-read_rds("4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.rds")
metadata <- sc@meta.data
metadata<-metadata %>% separate(SNG.BEST.GUESS,into=c("Pcase", "Origin2",sep="-",remove=FALSE))
metadata$Origin[is.na(metadata$Origin)]<-metadata$Origin2[is.na(metadata$Origin)]
metadata<-metadata %>% filter (!Pregnancy_ID %in% c("HPL20888", "HPL20874", "HPL20922"))
metadata<-metadata[,which(colnames(metadata)!="Condition")]
parturion_samples_data<- parturion_samples %>% select(Pregnancy_ID,Condition=Labor)
metadata<-metadata %>% rownames_to_column("barcode") %>% inner_join(parturion_samples_data)
rownames(metadata)<-metadata$barcode





cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
clust2Names<-cell.type.annotation$Potential.final 
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)

cluster.Colors<-cell.type.annotation$color
names(cluster.Colors)<-clust2Names


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


# cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
#                   "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
#                   "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")
# 
# 
# cluster.Colors.rest<-cluster.Colors.rest[which(cluster.Colors.rest %in% cluster.Colors)]
# 
# cluster.Colors<-c(cluster.Colors,cluster.Colors.rest)
# names(cluster.Colors)<-as.character(c(0:43))
# 
# names(cluster.Colors)[1:30]<-clust2Names
# cluster.Colors<-cluster.Colors[1:30]
###########################################

outFolder=paste0("./5_harmony_cellClass_covidcontrol_plots_final_Roger2/")
system(paste0("mkdir -p ", outFolder))
    
#sc$cluster_name<-clust2Names[sc$seurat_clusters]
    
    
   # aa <- FetchData(sc,c("Pregnancy_ID","UMAP_1","UMAP_2","Location","Labor","Condition","Origin","status","FetalSex","seurat_clusters","cluster_name","SNG.BEST.GUESS")) 
   # head(aa)
    
    
    aa <- FetchData(sc,c("UMAP_1","UMAP_2")) 
    aa$barcode<-rownames(aa)
    
    
    metadata$barcode<-rownames(metadata)
    
    
    aa<-aa %>% inner_join(metadata)
    
    #aa<-aa %>% separate(SNG.BEST.GUESS,into=c("Pcase", "Origin",sep="-",remove=FALSE))
    
    aa<-aa %>% filter(!is.na(Origin))
    fname=paste0(outFolder,"UMAP_LocationHarmony.Origin.pdf");
    pdf(fname,width=14,height=5)
    #fname=paste0(outFolder,"UMAP_LocationHarmony.Origin.png");
    #png(fname,width=1600,height=1200)
    p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Origin)) +
        geom_point(size=0.1) +
        ##    scale_color_manual(values=group.colors) +
        guides(colour = guide_legend(override.aes = list(size=10),title="Cell origin")) +
        scale_color_manual("Origin",values=c("M"="#D1D1D1","F"="#A61BB5"))+
        facet_wrap(~Location) +
        theme_bw()
    p2
    ##    theme_black()
    dev.off()
    
    #aa<-aa %>% filter(Labor!="NA" & status!="NA")
    
    
    ## Make a simple plot here:
   
    fname=paste0(outFolder,"UMAP_LocationHarmony.png");
    png(fname,width=1600,height=1200)
    #aa$seurat_clusters <- clust2name[aa$seurat_clusters]
    p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Location)) +
        geom_point(size=0.1) +
        ##    scale_color_manual(values=group.colors) +
        theme_bw()+
        scale_colour_manual(values=c("CAM"="#F4B183","PVBP"="#8AD2CD"))+
        #scale_color_manual(values=c("Control"="#333399","E. coli"="#A50021"))+
        #theme(legend.text=element_text(size=30,face="bold"), axis.text=element_text(size=30,face="bold"), axis.title=element_text(size=20,face="bold"))+
        ##    facet_wrap(~LocTime) +
        guides(colour = guide_legend(override.aes = list(size=10)),title="Location")+
        theme(text = element_text(size=30,face = "bold"),
              plot.title = element_text(size = 25, face = "bold"),
              legend.title=element_text(size=25,face="bold"),
              legend.text=element_text(size=25,face="bold"))
         p1
    ##    theme_black()
    dev.off()
    
    
    
    ## Make a simple plot here:
    fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
    png(fname,width=1600,height=1200)
    p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
        geom_point(size=0.1) +
        scale_color_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
        ##    scale_color_manual(values=group.colors) +
        guides(colour = guide_legend(override.aes = list(size=10)),title="Condition") +
        ##    facet_wrap(~LocTime) +
        theme_bw()+
        theme(text = element_text(size=30,face = "bold"),
              plot.title = element_text(size = 25, face = "bold"),
              legend.title=element_text(size=25,face="bold"), 
              legend.text=element_text(size=25,face="bold"))
    
    p1
    ##    theme_black()
    dev.off()
    
    
    ## Make a simple plot here:
    fname=paste0(outFolder,"UMAP_ConditionHarmony_location.png");
    png(fname,width=1900,height=1200)
    p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
        geom_point(size=0.1) +
        scale_color_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
        ##    scale_color_manual(values=group.colors) +
        guides(colour = guide_legend(override.aes = list(size=10)),title="Condition") +
        ##    facet_wrap(~LocTime) +
        theme_bw()+
        theme(text = element_text(size=30,face = "bold"),
              plot.title = element_text(size = 25, face = "bold"),
              legend.title=element_text(size=25,face="bold"), 
              legend.text=element_text(size=25,face="bold"))+
        facet_wrap(~Location) 
    
    p1
    ##    theme_black()
    dev.off()
    
    
    # ## Make a simple plot here:
    # fname=paste0(outFolder,"UMAP_OriginHarmony.png");
    # png(fname,width=1600,height=1200)
    # p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Origin)) +
    #     geom_point(size=0.1) +
    #     ##    scale_color_manual(values=group.colors) +
    #     guides(colour = guide_legend(override.aes = list(size=5),title="Origin")) +
    #     ##    facet_wrap(~LocTime) +
    #     theme_bw()
    # p1
    # ##    theme_black()
    # dev.off()
    
    
    aa<-aa %>% filter(Labor!="NA" & status!="NA")
    
    ## Make a simple plot here:
    fname=paste0(outFolder,"UMAP_StatusSoC_Harmony.png");
    png(fname,width=1600,height=1200)
    p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=status)) +
        geom_point(size=0.1) +
        ##    scale_color_manual(values=group.colors) +
        guides(colour = guide_legend(override.aes = list(size=10),title="Status")) +
        theme_bw()+
        theme(text = element_text(size=30,face = "bold"),
              plot.title = element_text(size = 25, face = "bold"),
              legend.title=element_text(size=25,face="bold"), 
              legend.text=element_text(size=25,face="bold"))
        ##    facet_wrap(~LocTime) +
        
    p1
    ##    theme_black()
    dev.off()
    
    
    
    
    aa<-aa %>% filter(Condition!="NA" & status!="NA")
    ## Make a simple plot here:
    fname=paste0(outFolder,"UMAP_LocationHarmony.Cell.annotation.png");
    png(fname,width=2200,height=900)
    p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
        geom_point(size=0.1) +
        ##    scale_color_manual(values=group.colors) +
        scale_color_manual(values=cluster.Colors) +
        guides(colour = guide_legend(override.aes = list(size=15),title="Cell type")) +
        #facet_wrap(~Location) +
        facet_grid(Condition ~ Location) +
        theme_bw()+
        theme(text = element_text(size=30,face = "bold"),
              plot.title = element_text(size = 25, face = "bold"),
              legend.title=element_text(size=25,face="bold"), 
              legend.text=element_text(size=25,face="bold"))+
    facet_wrap(~Location) 
    p1
    ##    theme_black()
    dev.off()
    
    
    
    aa<-aa %>% filter(Condition!="NA" & status!="NA")
    ## Make a simple plot here:
    fname=paste0(outFolder,"UMAP_LocationHarmony.Cell_annotation.location.png");
    png(fname,width=2200,height=1400)
    p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
        geom_point(size=0.1) +
        ##    scale_color_manual(values=group.colors) +
        scale_color_manual(values=cluster.Colors) +
        guides(colour = guide_legend(override.aes = list(size=15),title="Cell type")) +
        #facet_wrap(~Location) +
        facet_grid(Condition ~ Location) +
        theme_bw()+
        theme(text = element_text(size=30,face = "bold"),
              plot.title = element_text(size = 25, face = "bold"),
              legend.title=element_text(size=25,face="bold"), 
              legend.text=element_text(size=25,face="bold"))
    p1
    ##    theme_black()
    dev.off()
    
    
    aa<-aa %>% filter(Condition!="NA" & status!="NA")
    ## Make a simple plot here:
    fname=paste0(outFolder,"UMAP.png");
    png(fname,width=2200,height=1400)
    p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
        geom_point(size=0.1) +
        ##    scale_color_manual(values=group.colors) +
        scale_color_manual(values=cluster.Colors) +
        guides(colour = guide_legend(override.aes = list(size=15),title="Cell type")) +
        #facet_wrap(~Location) +
        #facet_grid(Labor ~ Location) +
        theme_bw()+
        theme(text = element_text(size=30,face = "bold"),
              plot.title = element_text(size = 25, face = "bold"),
              legend.title=element_text(size=25,face="bold"), 
              legend.text=element_text(size=25,face="bold"))
    p1
    ##    theme_black()
    dev.off()
    
    
    
    
    mycol=cluster.Colors
    names(mycol)<-as.character(c(0:32))
    
    
    fname=paste0(outFolder,"UMAP_Harmony.png");
    png(fname,width=1000,height=1000)
    DimPlot(sc, cols=mycol,reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) + NoLegend()
    dev.off()
    
    
    ### I think Val asked for this plot. RPR maybe double check.     
    aa<-aa %>% filter(Origin!="NA")
    fname=paste0(outFolder,"UMAP_origin.Barplot.pdf");
    pdf(fname,width=10,height=6)
    p2 <- ggplot(aa,aes(x=cluster_name,fill=Origin)) +
        geom_bar(position="stack") +
        ##    scale_color_manual(values=group.colors) +
        guides(colour = guide_legend(override.aes = list(size=10),title="Cell origin")) +
        scale_fill_manual("Origin",values=c("M"="#D1D1D1","F"="#A61BB5"))+
        facet_grid(.~Location) + coord_flip() +
        theme_bw()
    p2
    ##    theme_black()
    dev.off()
    
    ## Adjust for removing some cells. 
    
    ## Min 30 cells per cluster/origin/compartment/individual
    aa <- aa %>% mutate(gname=paste(Location,cluster_name,Origin,sep="_"))
    
    
    ca2 <- aa %>% group_by(gname) %>% summarize(n=n()) %>% filter(n>600)
    aa2 <- aa %>% filter(gname %in% ca2$gname)
    
    
    fname=paste0(outFolder,"UMAP_origin.Barplot2.pdf");
    pdf(fname,width=10,height=6)
    p2 <- ggplot(aa2,aes(x=cluster_name,fill=Origin)) +
      geom_bar(position="stack") +
      ##    scale_color_manual(values=group.colors) +
      guides(colour = guide_legend(override.aes = list(size=10),title="Cell origin")) +
      scale_fill_manual("Origin",values=c("M"="#D1D1D1","F"="#A61BB5"))+
      facet_grid(.~Location) + coord_flip() +
      theme_bw()
    p2
    ##    theme_black()
    dev.off()
    
    fname=paste0(outFolder,"UMAP_Condition.Barplot2.pdf");
    pdf(fname,width=10,height=6)
    p2 <- ggplot(aa2,aes(x=cluster_name,fill=Condition)) +
      geom_bar(position="stack") +
      ##    scale_color_manual(values=group.colors) +
      guides(colour = guide_legend(override.aes = list(size=10),title="Condition")) +
      scale_fill_manual("Condition", values=c("TNL"="#333399","TIL"="#A50021")) +
      facet_grid(.~Location) + coord_flip() +
      theme_bw()
    p2
    ##    theme_black()
    dev.off()

    ## ENSG00000160654.10  CD3G
    ## ENSG00000167286.9   CD3D
    ## ENSG00000198851.9   CD3E
    ## ENSG00000168685.15  IL7R
    
    a2 <- FetchData(sc,c("ENSG00000160654.10","ENSG00000167286.9","ENSG00000198851.9","ENSG00000168685.15"))
    
    rs.T <- rowSums(a2)[aa$barcode]
    
    aa$rs.T <- rs.T
    
    ## It does not seem that cluster 16 is a T-cell
    aa15 <- aa%>% filter(seurat_clusters==8,FetalSex=="Male",Location=="CAM") 
    
    plot(aa15$percent.Y,aa15$rs.T)
    
    fname=paste0(outFolder,"barplot_pregnancyID.pdf");
    pdf(fname,width=10,height=6)
    p2 <- ggplot(aa,aes(x=cluster_name,fill=Pregnancy_ID)) +
        geom_bar(position="stack") +
        ##    scale_color_manual(values=group.colors) +
        guides(colour = guide_legend(override.aes = list(size=15),title="Pregnancy ID")) +
        facet_grid(.~Location) + coord_flip() +
        theme_bw()
    p2
    ##    theme_black()
    dev.off()
    
    
    
    fname=paste0(outFolder,"UMAP_celltype.origin.location.Barplot.pdf");
    pdf(fname,width=10,height=6)
    p2 <- ggplot(aa,aes(x=Pregnancy_ID,fill=Origin)) +
        geom_bar(position = position_stack(reverse = TRUE)) +
        guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
        facet_grid(.~Location) + coord_flip() +
        scale_fill_manual(values=c("F"="#333399","M"="#A50021"))+
        theme_bw()
    p2
    ##    theme_black()
    dev.off()
    
    
    
    pregnancy.Colors<-cluster.Colors[1:length(unique(aa$Pregnancy_ID))]
    names(pregnancy.Colors)<-unique(aa$Pregnancy_ID)
    
    fname=paste0(outFolder,"UMAP_PregnancyID.png");
    png(fname,width=1600,height=1200)
    p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Pregnancy_ID)) +
        scale_color_manual(values=pregnancy.Colors) +
        geom_point(size=0.1) +
        ##    scale_color_manual(values=group.colors) +
        guides(colour = guide_legend(override.aes = list(size=10),title="Pregnancy ID")) +
        theme_bw()+
        theme(text = element_text(size=30,face = "bold"),
              plot.title = element_text(size = 25, face = "bold"),
              legend.title=element_text(size=25,face="bold"), 
              legend.text=element_text(size=25,face="bold"))
    ##    facet_wrap(~LocTime) +
    
    p1
    ##    theme_black()
    dev.off()
    
    
    #aa<-aa %>% filter(Origin!="NA")
    
    
    aa<-aa %>% filter(Condition!="NA")
    aa<-aa %>% filter(Origin!="NA")
    
for (i in 1:length(unique(aa$cluster_name))){
 #sapply(unique(aa$cluster_name),function(x){
  
     print(i)
     x<-unique(aa$cluster_name)[i]
     print(x)
     aa2<-aa %>% filter (cluster_name==x)#"4:Marophage-2 (Hofbauer)")
     x<-gsub(":","_",x)
     fname=paste0(outFolder,"UMAP_",x,".png")#Marophage-2-Hofbauer.pdf");
     #fname=paste0(outFolder,"UMAP_PregnancyID.png");
     png(fname,width=1800,height=1500)
     #pdf(fname,width=10,height=6)
     p2 <- ggplot(aa2,aes(x=Pregnancy_ID,fill=Condition)) +
         geom_bar(position = position_stack(reverse = TRUE)) +
         guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
         #facet_grid(.~Location) + coord_flip() +
         facet_grid( Location ~ Origin ) + coord_flip() +
         #scale_fill_manual(values=c("F"="#333399","M"="#A50021"))+
         scale_fill_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
         theme_bw()+
         theme(text = element_text(size=30,face = "bold"),
               plot.title = element_text(size = 25, face = "bold"),
               legend.title=element_text(size=25,face="bold"), 
               legend.text=element_text(size=25,face="bold"))
     
     p2
     ##    theme_black()
     dev.off()    
 }  
   

  
    
    
    
    
    
