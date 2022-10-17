library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)
library(NMF)
library(ggalluvial)


##########################################################################
# cellchat plots 1
##########################################################################


outFolder="./10_CellChat_analysis_withcovidcontrol_500_v2_plots/"
system(paste0("mkdir -p ", outFolder))

cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
#rownames(cell.type.annotation)<-NULL
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)



cluster.Colors<-cell.type.annotation$color
names(cluster.Colors)<-clust2Names

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)



########################################################### 
# circle plots per pathways
###########################################################

locations<-c("CAM" ,    "PVBP")
conditions<-c("TIL" ,"TNL")


sapply(locations, function(xlocation){
  
  subFolder<-xlocation
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))

sapply (conditions ,function(xcondition){
  
  print(xcondition)
  
  #cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_200/","cellchat_",xlocation,"_",xcondition,"_","2022-04-07",".rds"))
  #cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_200/","cellchat_",xlocation,"_",xcondition,"_2022-05-05.rds"))
  cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_500/","cellchat_",xlocation,"_",xcondition,"_","2022-05-05",".rds"))
  
  pathways<-cellchat@netP$pathways
  groupSize <- as.numeric(table(cellchat@idents))
  
  pathways<-cellchat@netP$pathways
  
  selected_pathways<-pathways #intersect(pathways,selected_pathways)
  groupSize <- as.numeric(table(cellchat@idents))
  
  
  
  sapply(selected_pathways,function(pathways.show){
    
    
    filename<-"/Pathways_circleplots_top25/"
    # 
    system(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,filename))
    
    
    pdf(paste0(outFolder,subFolder,"/",xcondition,filename,pathways.show,".pdf"),width=12,height=12)
    par(mfrow=c(1,1))
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",color.use=cluster.Colors[rownames(cellchat@net$weight)],vertex.label.cex=1,top=0.25,arrow.width = 2,arrow.size=1)
    dev.off()
    
    pdf(paste0(outFolder,subFolder,"/",xcondition,filename,pathways.show,"_nolabel.pdf"),width=12,height=12)
    
    par(mfrow=c(1,1))
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",color.use=cluster.Colors[rownames(cellchat@net$weight)],vertex.label.cex=0.000001,top=0.25,arrow.width = 2,arrow.size=1)
    dev.off()
    
  })
  
  
})

})

# pathways to focus 
#informationflow_data<-read_rds(paste0("10_CellChat_comparison_conditions_after_filter_200_plots/informationflow_data.rds"))
#selected_pathways_mat<-read.delim("selectedpathways.txt")

# locations
locations<-c("CAM" ,    "PVBP")
conditions<-c("TIL" ,"TNL")


#######################################################################
# alluvial plots
#######################################################################


selected_pathways<-c("IL6","CD80","RESISTIN","TENASCIN","SELE","VEGF","IL1","MHC-II","GALECTIN")

sapply(locations, function(xlocation){
  


lapply (conditions ,function(xcondition){
#for (xcondition in conditions){
    
  #
  subFolder<-xlocation
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))
    
    #cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_200/","cellchat_",xlocation,"_",xcondition,"_","2022-04-07",".rds"))
    cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_500/","cellchat_",xlocation,"_",xcondition,"_","2022-05-05",".rds"))
    
    pathways<-cellchat@netP$pathways
    groupSize <- as.numeric(table(cellchat@idents))
    
    
   
    
    #informationflow_data_specific<-informationflow_data %>% filter(Location==xlocation) 
    
    # top 10 %
    #informationflow_data_specific<-informationflow_data_specific %>% filter(Condition=="Ecoli" ) %>% arrange(pvalues,-contribution.scaled) 
    
    #top10<-round(0.1 *nrow(informationflow_data_specific))
    
    #top10<-15
    #informationflow_data_specific<-informationflow_data_specific %>% top_n(top10,w=c(pvalues))
    #informationflow_data_specific<-informationflow_data_specific %>% head(top10)
    
    
    #informationflow_data_specific<-informationflow_data_specific %>% filter(Condition=="Ecoli" & pvalues<=0.01 ) 
    #selected_pathways<-informationflow_data_specific %>% select (Pathway) %>% unlist %>% unique() %>% as.character()
    
    #selected_pathways<-selected_pathways_mat[,xlocation]
    #selected_pathways<-selected_pathways[selected_pathways!=""]
     system(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,"/"))
    # 
    # # 
    # system(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,"/Pathways_circleplots_top25"))
    # 
    # sapply(pathways,function(pathways.show){
    #   pdf(paste0(outFolder,subFolder,"/",xcondition,"/Pathways_circleplots_top25/",pathways.show,".pdf"),width=20,height=20)
    #   par(mfrow=c(1,1))
    #   #arrow(length = unit(.02, "inches"),type = "closed",angle = 40)
    #   netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",color.use=cluster.Colors[rownames(cellchat@net$weight)],vertex.label.cex=2,top=0.25,arrow.width = 15)
    #   dev.off()
    # })
    # 
    
  ############################################################################
  # alluvial plots
  ############################################################################
  
  
  
  nPatterns<-3 #
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

  cutoff<-0.75

  ############################################################################
  # outgoing data
  ############################################################################
  outgoing_signaling<-cellchat@netP$pattern$outgoing$pattern$signaling
  outgoing_signaling<-outgoing_signaling %>% filter(Contribution>cutoff)

  colnames(outgoing_signaling)[3]<-"Contribution_outgoing_signaling"
  colnames(outgoing_signaling)[1]<-paste0(colnames(outgoing_signaling)[1],"_outgoing")

  outgoing_cell<-cellchat@netP$pattern$outgoing$pattern$cell
  outgoing_cell<-outgoing_cell %>% filter(Contribution>cutoff)

  colnames(outgoing_cell)[3]<-"Contribution_outgoing_cell"
  colnames(outgoing_cell)[2]<-paste0(colnames(outgoing_cell)[2],"_outgoing")

  outgoing_data<-outgoing_signaling %>% inner_join(outgoing_cell)
  outgoing_data<-outgoing_data %>% filter(Signaling %in% selected_pathways)


  ############################################################################
  # incoming data
  ############################################################################
  incoming_signaling<-cellchat@netP$pattern$incoming$pattern$signaling
  incoming_signaling<-incoming_signaling %>% filter(Contribution>cutoff)
  #incoming_signaling$Contribution[incoming_signaling$Contribution < cutoff] <- 0

  colnames(incoming_signaling)[3]<-"Contribution_incoming_signaling"
  colnames(incoming_signaling)[1]<-paste0(colnames(incoming_signaling)[1],"_incoming")

  incoming_cell<-cellchat@netP$pattern$incoming$pattern$cell
  incoming_cell<-incoming_cell %>% filter(Contribution>cutoff)
  #incoming_cell$Contribution[incoming_cell$Contribution < cutoff] <- 0

  colnames(incoming_cell)[3]<-"Contribution_incoming_cell"
  colnames(incoming_cell)[2]<-paste0(colnames(incoming_cell)[2],"_incoming")

  incoming_data<-incoming_signaling %>% inner_join(incoming_cell)

  incoming_data<-incoming_data %>% filter(Signaling %in% selected_pathways)
  
  outgoing_data<-outgoing_data %>% filter(Signaling %in% selected_pathways)
  
  outgoing_data<-outgoing_data %>% mutate(CellGroup_outgoing=CellGroup)#,Signaling,Contribution-outgoing_cell,Contribution_outgoing_signaling)
  incoming_data<-incoming_data %>% mutate(CellGroup_incoming=CellGroup)#,Signaling,Contribution_incoming_signal,Contribution_incoming_cell)
  signalins_incoming<-unique(incoming_data$Signaling)
  outgoing_incoming<-unique(outgoing_data$Signaling)
  overlap_signalings<-intersect(signalins_incoming,outgoing_incoming)


  ############################################################################
  # outgoing alluvial plot
  ############################################################################

  
  # summary 
  majorcells<-rowSums(table(outgoing_data$CellGroup,outgoing_data$Signaling))
  majorcells<-majorcells[order(majorcells,decreasing = TRUE)]
  majorsignalings<-rowSums(table(outgoing_data$Signaling,outgoing_data$Contribution_outgoing_signaling))
  majorsignalings<-majorsignalings[order(majorsignalings,decreasing = TRUE)]
  
  major_cell_signaling<-list(majorcells,majorsignalings)
  names(major_cell_signaling)<-c("majorcells","majorsignalings")
  
  
  outgoing_data$Signaling<-as.character(outgoing_data$Signaling)
  outgoing_data<-outgoing_data[order(outgoing_data$Signaling,decreasing = FALSE),]
  outgoing_data$Signaling <- factor(outgoing_data$Signaling,levels=unique(outgoing_data$Signaling))

  #system(paste0("mkdir -p ", outFolder,subFolder,"/"))

  #fname=paste0(outFolder,subFolder,"/",xcondition,"/",xcondition,"_outgoing_ggalluvial_cuttoff.0.7.pdf");
  fname=paste0(outFolder,subFolder,"/",xcondition,"/",xcondition,"_outgoing_ggalluvial_selectedPathways_cuttoff_",Sys.Date(),"_.pdf");
  
  pdf(fname,width=20,height=20)
  ggplot(data = outgoing_data,
         aes(axis1 = CellGroup_outgoing, axis2 = Signaling, y = Contribution_outgoing_signaling,fill=CellGroup_outgoing)) +
    geom_alluvium(aes(fill=CellGroup_outgoing)) +#aes(fill = Signaling)
    geom_stratum() +
    geom_flow()+
    scale_fill_manual("Cell type",values=cluster.Colors) +
    geom_text(stat = "stratum",aes(label = after_stat(stratum)),size=5) +
    scale_x_discrete(limits = c("Outgoing cell", "Pathway"),
                     expand = c(0.3, 0.1)) +
    theme_bw()
  dev.off()


  #outgoing_data<- outgoing_data %>% filter(Signaling %in% overlap_signalings)
  # fname=paste0(outFolder,subFolder,"/",subFolder,"_outgoing_overlap_ggalluvial.pdf");
  # pdf(fname,width=20,height=28)
  # ggplot(data = outgoing_data,
  #        aes(axis1 = CellGroup_outgoing, axis2 = Signaling, y = Contribution_outgoing_signaling,fill=CellGroup_outgoing)) +
  #   geom_alluvium(aes(fill=CellGroup_outgoing)) +#aes(fill = Signaling)
  #   geom_stratum() +
  #   geom_flow()+
  #   scale_fill_manual("Cell type",values=cluster.Colors) +
  #   geom_text(stat = "stratum",aes(label = after_stat(stratum)),size=5) +
  #   scale_x_discrete(limits = c("Outgoing cell", "Pathway"),
  #                    expand = c(0.3, 0.1)) +
  #   theme_bw()
  # dev.off()

  ############################################################################
  # incoming alluvial plot
  ############################################################################

  incoming_data$Signaling<-as.character(incoming_data$Signaling)
  incoming_data<-incoming_data[order(incoming_data$Signaling,decreasing = FALSE),]
  incoming_data$Signaling <- factor(incoming_data$Signaling,levels=unique(incoming_data$Signaling))


  #fname=paste0(outFolder,subFolder,"/",xcondition,"/",xcondition,"_incoming_ggalluvial_cuttoff.0.7.pdf");
  
  fname=paste0(outFolder,subFolder,"/",xcondition,"/",xcondition,"_incoming_ggalluvial_selectedPathways_cuttoff_",Sys.Date(),"_.pdf")
  pdf(fname,width=20,height=20)
  ggplot(data = incoming_data,
         aes(axis1 = Signaling, axis2 = CellGroup_incoming, y = Contribution_incoming_signaling)) +
    geom_alluvium() +#aes(fill = Signaling)
    geom_stratum(aes(fill=CellGroup_incoming)) +
    geom_flow(aes(fill=CellGroup_incoming))+
    scale_fill_manual("Cell type",values=cluster.Colors) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)),size=5) +
    scale_x_discrete(limits = c("Pathway", "Incoming cell"),
                     expand = c(0.15, 0.05)) +
    theme_bw()
  dev.off()


  }  )

})











########################################################################
# pathways, DEGs , LR pairs 
########################################################################

#outFolder="./10_CellChat_analysis_withcovidcontrol_200_plots/"


locations<-c("CAM" ,    "PVBP")
conditions<-c("TIL" ,"TNL")


cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
#rownames(cell.type.annotation)<-NULL
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)



cluster.Colors<-cell.type.annotation$color
names(cluster.Colors)<-clust2Names

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)


#selected_pathways_mat<-read.delim("selectedpathways.txt")



pathways_genes_celltypes_db_final<-lapply(locations, function(xlocation){
  
  subFolder<-xlocation
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  #pathways_genes_celltypes_db_list<-lapply (conditions ,function(xcondition){
    
  pathways_genes_celltypes_db_list<-list()
   for (i in 1:length(conditions)){
    xcondition<-conditions[i]
    
    #cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_200/","cellchat_",xlocation,"_",xcondition,"_","2022-04-07",".rds"))
    cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_500/","cellchat_",xlocation,"_",xcondition,"_","2022-05-05",".rds"))
    
    pathways<-cellchat@netP$pathways
    groupSize <- as.numeric(table(cellchat@idents))
    
    
    system(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,"/"))
    print(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,"/"))

    pathways_genes_celltypes<-lapply(pathways,function(x){
      
      pairLR <- extractEnrichedLR(cellchat, signaling = x, geneLR.return = FALSE)
      LR_genes<-unique(unlist(str_split(pairLR$interaction_name,"_")))
      #LR_genes<-tolower(LR_genes)
      
      celltypes<-res %>% filter(gene_name %in% LR_genes & Location== xlocation)%>% dplyr::select(Cell_type)%>% unlist %>% unique()
      
      LR_DEGs <-res %>% filter (gene_name %in% LR_genes & Location== xlocation) %>% dplyr::select(gene_name)%>% unlist %>% unique()
      
      LR_DEGs_len<-length(LR_DEGs)
      LR_DEGs<-paste(LR_DEGs,collapse = ", ")
      
      LR_genes<-paste(LR_genes,collapse = ", ")
      
      celltype_len<-length(celltypes)
      celltypes<-paste(celltypes,collapse = ", ")
      if (celltypes=="") celltypes<-"NA"
      return(c(x,LR_genes,LR_DEGs,celltypes,celltype_len,LR_DEGs_len))
    })
    pathways_genes_celltypes_db<-do.call(rbind,pathways_genes_celltypes)
    
    pathways_genes_celltypes_db<-as.data.frame(pathways_genes_celltypes_db)
    
    print(xcondition)
    pathways_genes_celltypes_db$condition<-xcondition
    pathways_genes_celltypes_db$location<-xlocation
    colnames(pathways_genes_celltypes_db)<-c("Pathway","Genes_in_pathway","DEGs_in_pathway","Celltype","celltype_len","LR_DEGs_len","Condition","Location")
    
    pathways_genes_celltypes_db_list[[i]]<-pathways_genes_celltypes_db
    #return(pathways_genes_celltypes_db)
    
  }  #)
  
  pathways_genes_celltypes_db<-do.call(rbind,pathways_genes_celltypes_db_list)
  return(pathways_genes_celltypes_db) 
  
})


pathways_genes_celltypes_db<-do.call(rbind,pathways_genes_celltypes_db_final)

pathways_genes_celltypes_db<-pathways_genes_celltypes_db[order(as.numeric(pathways_genes_celltypes_db$celltype_len),decreasing = TRUE),]

write.csv(pathways_genes_celltypes_db,file=paste0(outFolder,"pathways_genes_celltypes_db.csv"))
write_rds(pathways_genes_celltypes_db,file=paste0(outFolder,"pathways_genes_celltypes_db.rds"))



########################################################################
# pathways, information flow
########################################################################


informationflow_list<-lapply(locations,function(xlocation){
  
  
  subFolder<-xlocation
  
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  
  
  cellchat_TNL<-cellchat_TNL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_500/","cellchat_",xlocation,"_","TNL","_","2022-05-05",".rds"))
  cellchat_TIL<-cellchat_TIL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_500/","cellchat_",xlocation,"_","TIL","_","2022-05-05",".rds"))
  
  
  pathways_TNL<-cellchat_TNL_original@netP$pathways
  controls<-cbind(pathways_TNL, rep("TNL",length(pathways_TNL)))
  controls<-as.data.frame(controls)
  colnames(controls)<-c("pathway","TNL")
  
  
  pathways_TIL<-cellchat_TIL_original@netP$pathways
  TILs<-cbind(pathways_TIL, rep("TIL",length(pathways_TIL)))
  TILs<-as.data.frame(TILs)
  colnames(TILs)<-c("pathway","TIL")
  
  allpathways<-TILs %>% full_join(controls)
  
  allpathways<-allpathways %>% arrange(desc(!is.na(TNL),(!is.na(TIL))))
  
  write.csv(allpathways,file=paste0(outFolder,subFolder,"/","pathways_",xlocation,".csv"))
  
  shared_pathways<-intersect(pathways_TNL,pathways_TIL)
  pathways_TNL_only<-pathways_TNL[!pathways_TNL %in% pathways_TIL]
  pathways_TIL_only<-pathways_TIL[! pathways_TIL %in% pathways_TNL]
  
  
  
  
  # there is an additional population  30_B cell specific to cellchat_TNL compared to TIL
  # we lift up Ecoli by lifting up the cell groups to the same cell labels as cellchat_TNL 
  
  control_cells<-rownames(cellchat_TNL@net$prob)
  TIL_cells<-rownames(cellchat_TIL@net$prob)
  
  TIL_cells_only<-TIL_cells [!TIL_cells %in%control_cells ]
  control_cells_only<-control_cells [! control_cells  %in%TIL_cells ]
  
  if (length(TIL_cells_only)==0)
  {
    group.new = levels(cellchat_TNL@idents)
    cellchat_TIL <- liftCellChat(cellchat_TIL, group.new)
  }else if (length(control_cells_only)==0){
    group.new = levels(cellchat_TIL@idents)
    cellchat_TNL <- liftCellChat(cellchat_TNL, group.new)} else 
    {
      group.new = union(levels( cellchat_TIL@idents), levels(cellchat_TNL@idents))
      cellchat_TNL <- liftCellChat(cellchat_TNL, group.new)
      cellchat_TIL <- liftCellChat(cellchat_TIL, group.new)
    }
  
  
  # now merge
  object.list <- list(TNL = cellchat_TNL, TIL = cellchat_TIL)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  
  
  
  
  #pdf(paste0(outFolder,subFolder,"/","informationflow_",xlocation,".pdf"),width=10,height=18)
  #gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("#333399","#A50021"),font.size=15)
  # gg2 
  # dev.off()
  data<-gg2$data
  colnames(data)<-c("Pathway","contribution","contribution.scaled","Condition","contribution.relative.1","pvalues" )
  data$Location<-xlocation
  #data$Condition[which(data$Condition=="Control")]<-"control"
  data<-data %>% select(Pathway,Location,Condition,contribution,contribution.scaled,pvalues)
  
  return (data)
  })



informationflow_data<-do.call(rbind,informationflow_list)

write.csv(informationflow_data,file=paste0(outFolder,"informationflow_data.csv"))
write_rds(informationflow_data,paste0(outFolder,"informationflow_data.rds"))



#informationflow_data<-read_csv("10_CellChat_analysis_withcovidcontrol_500_v2_plots/informationflow_data.rds")




#informationflow_data$Condition<-tolower(informationflow_data$Condition)
#pathways_genes_celltypes_db$Condition<-tolower(pathways_genes_celltypes_db$Condition)
informationflow_data$Pathway<-as.character(informationflow_data$Pathway)
informationflow_data$Condition<-as.character(informationflow_data$Condition)
informationflow_data$Location<-as.character(informationflow_data$Location)




pathways_genes_celltypes_contribution<-pathways_genes_celltypes_db %>% left_join(informationflow_data)

write.csv(pathways_genes_celltypes_contribution,file=paste0(outFolder,"pathways_genes_celltypes_contribution.csv"))

write_rds(pathways_genes_celltypes_contribution,file=paste0(outFolder,"pathways_genes_celltypes_contribution.rds"))


####################################################################

# selecting top pathways based on information flows
####################################################################
#informationflow_data<-read_rds("10_CellChat_analysis_withcovidcontrol_500_v2_plots/informationflow_data.rds")

pathways<-sapply(locations, function(xlocation){
  
  subFolder<-xlocation
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))

  #for (xcondition in conditions){
  
  xcondition<-"E. coli"
  
  cellchat<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_",xcondition,"_","2021-12-21",".rds"))
  pathways<-cellchat@netP$pathways
  groupSize <- as.numeric(table(cellchat@idents))
  
  
  if (xcondition=="E. coli") xcondition<-"Ecoli" 
  
  informationflow_data_specific<-informationflow_data %>% filter(Location==xlocation) 
  
  # top 10 %
  #informationflow_data_specific<-informationflow_data_specific %>% filter(Condition=="Ecoli" ) %>% arrange(pvalues,-contribution.scaled) 
  
  informationflow_data_specific<-informationflow_data_specific %>% filter(Condition=="Ecoli" ) %>% arrange(pvalues) 
  
  #top10<-round(0.1 *nrow(informationflow_data_specific))
  
  top10<-30
  #informationflow_data_specific<-informationflow_data_specific %>% top_n(top10,w=c(pvalues))
  informationflow_data_specific<-informationflow_data_specific %>% head(top10)
  
  
  #informationflow_data_specific<-informationflow_data_specific %>% filter(Condition=="Ecoli" & pvalues<=0.01 ) 
  selected_pathways<-informationflow_data_specific %>% select (Pathway) %>% unlist %>% unique() %>% as.character()
 
return (selected_pathways)
})

