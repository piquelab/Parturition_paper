################################################################################################################################################
# Cellchat downstream analysis 
# obtaining pathways, cell types, DEGs, sender, reciever, pvalues,contribution degree
########################################################################
# pathways, DEGs , LR pairs 
########################################################################
locations<-c("CAM" ,    "PVBP")
conditions<-c("TIL" ,"TNL")

cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
#rownames(cell.type.annotation)<-NULL
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)


res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)
res$cname<-paste(res$Cell_type,res$Origin,sep="_")


pathways_genes_celltypes_db_final<-lapply(locations, function(xlocation){
  
  subFolder<-xlocation
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  #pathways_genes_celltypes_db_list<-lapply (conditions ,function(xcondition){
  
  pathways_genes_celltypes_db_list<-list()
  for (i in 1:length(conditions)){
    xcondition<-conditions[i]
    
    
    cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_origin_or_500/","cellchat_",xlocation,"_",xcondition,"_","2022-06-01",".rds"))
    
    pathways<-cellchat@netP$pathways
    groupSize <- as.numeric(table(cellchat@idents))
    
    
    system(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,"/"))
    print(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,"/"))
    
    pathways_genes_celltypes<-lapply(pathways,function(x){
      
      pairLR <- extractEnrichedLR(cellchat, signaling = x, geneLR.return = FALSE)
      LR_genes<-unique(unlist(str_split(pairLR$interaction_name,"_")))
      #LR_genes<-tolower(LR_genes)
      
      celltypes<-res %>% filter(gene_name %in% LR_genes & Location== xlocation)%>% dplyr::select(Cell_type)%>% unlist %>% unique()
      
      LR_DEGs <-res %>% filter (gene_name %in% LR_genes& Location== xlocation) %>% dplyr::select(gene_name)%>% unlist %>% unique()
      
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
  
  
  # cellchat_TNL<-cellchat_TNL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_origin_500/","cellchat_",xlocation,"_","TNL","_","2022-05-06",".rds"))
  # cellchat_TIL<-cellchat_TIL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_origin_500/","cellchat_",xlocation,"_","TIL","_","2022-05-06",".rds"))
  # 
  cellchat_TNL<-cellchat_TNL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_origin_or_500/","cellchat_",xlocation,"_","TNL","_","2022-06-01",".rds"))
  cellchat_TIL<-cellchat_TIL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_origin_or_500/","cellchat_",xlocation,"_","TIL","_","2022-06-01",".rds"))
  
  
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
  data<-data %>% dplyr::select(Pathway,Location,Condition,contribution,contribution.scaled,pvalues)
  
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



outFolder="10_CellChat_analysis_withcovidcontrol_origin_500_or_plots/"

#10_CellChat_analysis_withcovidcontrol_origin_500_plots
system(paste0("mkdir -p ", outFolder))


cell.type.annotation<-read_tsv("cell.type.annotation.v2.tsv")
#rownames(cell.type.annotation)<-NULL
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)


# cluster.Colors<-cell.type.annotation$color
# names(cluster.Colors)<-clust2Names


################################
# cells with sending role (from cellchat)
################################

locations<-c("CAM" ,    "PVBP")
conditions<-c("TIL" ,"TNL")




####################################################
####################################################


data_list<-lapply(locations, function(xlocation){
  
 
  
  data_list<-lapply (conditions ,function(xcondition){
    
  #cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_origin_500/","cellchat_",xlocation,"_",xcondition,"_","2022-05-06",".rds"))
  cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_origin_or_500/","cellchat_",xlocation,"_",xcondition,"_2022-06-01.rds"))
  
  pathways<-cellchat@netP$pathways
  groupSize <- as.numeric(table(cellchat@idents))
  nPatterns<-3 #
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  
  cutoff<-0.70
  
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
  #outgoing_data<-outgoing_data %>% filter(Signaling %in% selected_pathways)
  
  
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
  
  #incoming_data<-incoming_data %>% filter(Signaling %in% selected_pathways)
  
  #outgoing_data<-outgoing_data %>% filter(Signaling %in% selected_pathways)
  
  outgoing_data<-outgoing_data %>% mutate(CellGroup_outgoing=CellGroup)#,Signaling,Contribution-outgoing_cell,Contribution_outgoing_signaling)
  incoming_data<-incoming_data %>% mutate(CellGroup_incoming=CellGroup)#,Signaling,Contribution_incoming_signal,Contribution_incoming_cell)
  signalins_incoming<-unique(incoming_data$Signaling)
  outgoing_incoming<-unique(outgoing_data$Signaling)
  overlap_signalings<-intersect(signalins_incoming,outgoing_incoming)
  
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
  
  total_signalings<-unique(c(incoming_data$Signaling, outgoing_data$Signaling))
  
  incoming_data_df<-incoming_data %>% dplyr::select(Signaling,CellGroup_incoming)
  outgoing_data_df<-outgoing_data  %>% dplyr::select(Signaling,CellGroup_outgoing)
  data<-incoming_data_df %>% full_join (outgoing_data_df)
  data$Location<-xlocation
  data$Condition<-xcondition
  
  return(data)
  })
  data<-do.call(rbind,data_list)
return(data)
  })


data<-do.call(rbind,data_list)
  
    

data$signaling_Location_condition<-paste0(data$Signaling,"_",data$Location,"_",data$Condition)

senders<-sapply(unique(data$signaling_Location_condition), function(x)
{
  sender<-data$CellGroup_outgoing[which(data$signaling_Location_condition==x)]
  sender<-as.character(sender)
  sender<-unique(sender)
  sender<-paste(sender,collapse = ", ")
  sender
  
})


receivers<-sapply(unique(data$signaling_Location_condition), function(x)
{
  sender<-data$CellGroup_incoming[which(data$signaling_Location_condition==x)]
  sender<-as.character(sender)
  sender<-unique(sender)
  sender<-paste(sender,collapse = ", ")
  sender
  
})


data<-data.frame(unique(data$signaling_Location_condition),senders,receivers)
data<-data %>% separate(unique.data.signaling_Location_condition.,c("Pathway","Location","Condition"),sep="_",remove=FALSE) 
#data<-data %>% dplyr::select(Pathway=Signaling,Location,Condition,senders,receivers)
  

pathways_genes_celltypes_contribution<-read_rds("10_CellChat_analysis_withcovidcontrol_origin_500_or_plots/pathways_genes_celltypes_contribution.rds")
#pathways_genes_celltypes_contribution<-read_rds("10_CellChat_analysis_withcovidcontrol_200_v2_plots/pathways_genes_celltypes_contribution.rds")

pathways_genes_celltypes_db<-pathways_genes_celltypes_contribution %>% left_join(data)

write.csv(pathways_genes_celltypes_db,file = paste0(outFolder,"pathways_genes_celltypes_sender_receiver_outgoingincoming_cutoff.csv"))
write_rds(pathways_genes_celltypes_db,file = paste0(outFolder,"pathways_genes_celltypes_sender_receiver_outgoingincoming_cutoff.rds"))


########################################################################################################################
# Based on outdegree and indegree cut off 
########################################################################################################################



senders_list<-lapply(locations, function(xlocation){
  #subFolder<-xlocation
  #system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  #pathways_genes_celltypes_db_list<-lapply (conditions ,function(xcondition){
  
  senders_list<-lapply(conditions,function(xcondition){
    
    
    cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_200/","cellchat_",xlocation,"_",xcondition,"_","2022-04-07",".rds"))
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    pathways<-cellchat@netP$pathways 
    senders<-lapply(pathways,function(x){
      
      scores<-cellchat@netP$centr[[x]]$outdeg
      minx<-min(scores)
      maxx<-max(scores)
      scores_scaled<-sapply( scores, function(x){(x-minx)/(maxx-minx)})
      senders<-names(scores_scaled)[which(scores_scaled>=0.75)]
      roles<-rep("sender",length(senders))
      path<-rep(x,length(senders))
      return(cbind(senders,roles,path))
    })
    
    senders_df<-do.call(rbind,senders)
    senders_df<-as.data.frame(senders_df)
    senders_df$Condition<-xcondition
    senders_df
  })
  senders_df<-do.call(rbind,senders_list)
  senders_df<-as.data.frame(senders_df)
  senders_df$Location<-xlocation
  return (senders_df)
})

sender_df<-do.call(rbind,senders_list)


################################
# cells with receiving role (from cellchat)
################################


receivers_list<-lapply(locations, function(xlocation){
  #subFolder<-xlocation
  #system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  #pathways_genes_celltypes_db_list<-lapply (conditions ,function(xcondition){
  
  receivers_list<-lapply(conditions,function(xcondition){
    
    
    cellchat<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_200/","cellchat_",xlocation,"_",xcondition,"_","2022-04-07",".rds"))
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    pathways<-cellchat@netP$pathways 
    receivers<-lapply(pathways,function(x){
      
      scores<-cellchat@netP$centr[[x]]$indeg
      minx<-min(scores)
      maxx<-max(scores)
      scores_scaled<-sapply( scores, function(x){(x-minx)/(maxx-minx)})
      receivers<-names(scores_scaled)[which(scores_scaled>=0.75)]
      roles<-rep("receiver",length(receivers))
      path<-rep(x,length(receivers))
      return(cbind(receivers,roles,path))
    })
    
    receivers_df<-do.call(rbind,receivers)
    receivers_df<-as.data.frame(receivers_df)
    receivers_df$Condition<-xcondition
    receivers_df
  })
  receivers_df<-do.call(rbind,receivers_list)
  receivers_df<-as.data.frame(receivers_df)
  receivers_df$Location<-xlocation
  return (receivers_df)
})

receiver_df<-do.call(rbind,receivers_list)


################################
# merging sender/receiever cells
################################

colnames(sender_df)[1:3]<-c("Celltype","role","Pathway")
colnames(receiver_df)[1:3]<-c("Celltype","role","Pathway")

receiver_senders_df<-rbind(sender_df,receiver_df)
#colnames(receiver_senders_df)<-c("celltype","role","pathway")
receiver_senders_df<-receiver_senders_df[!duplicated(receiver_senders_df), ]
#receiver_senders_df<-receiver_senders_df[order(receiver_senders_df$pathway,decreasing = FALSE),]
#receiver_senders_df<-as.data.frame(receiver_senders_df)
receiver_senders_df<-receiver_senders_df %>% arrange(Pathway)

#write.csv(receiver_senders_df,file=paste0(outFolder,"receiver_senders_df.csv"))
write.csv(receiver_senders_df,file=paste0(outFolder,"receiver_senders_degree_0.75_df.csv"))

sender_receiver_columns<-sapply(1: nrow(pathways_genes_celltypes_db), function(x){
  path<-pathways_genes_celltypes_db$Pathway[x]
  celltypes<-pathways_genes_celltypes_db$Celltype[x]
  celltypes<-unlist(strsplit(celltypes,", ",fixed=TRUE))
  location<-pathways_genes_celltypes_db$Location[x]
  condition<-pathways_genes_celltypes_db$Condition[x]
  senders<- receiver_senders_df %>% filter (Pathway==path &  role=="sender" & Location==location & Condition==condition) %>% select(Celltype)%>% unlist %>% unique
  receivers<-receiver_senders_df %>% filter (Pathway==path & role=="receiver" & Location==location & Condition==condition) %>% select(Celltype)%>% unlist %>% unique
  
  senders<-paste(senders,collapse = ", ")
  receivers<-paste(receivers,collapse = ", ")
  if (senders=="") senders<-"NA"
  if (receivers=="") receivers<-"NA"
  return(c(senders,receivers,location,condition,path))
})



################################
# make the final data frame 
################################


sender_receiver_columns<-t(sender_receiver_columns)
sender_receiver_columns<-as.data.frame(sender_receiver_columns)
colnames(sender_receiver_columns)<-c("senders","receivers","Location","Condition","Pathway")


pathways_genes_celltypes_db<-pathways_genes_celltypes_contribution %>% left_join(sender_receiver_columns)

write.csv(pathways_genes_celltypes_db,file = paste0(outFolder,"pathways_genes_celltypes_sender_receiver_cutoffdegree.csv"))
write_rds(pathways_genes_celltypes_db,file = paste0(outFolder,"pathways_genes_celltypes_db.rds"))



