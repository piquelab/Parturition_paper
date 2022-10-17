library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)

##########################################################################
# cellchat plots 2
##########################################################################

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
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)




#color_tobechanged<-cluster.Colors[clust2Names[which(clust2Names %in% tobechanged)]]

#outFolder="./10_CellChat_analysis_withcovidcontrol_200_plots/"
#outFolder="./10_CellChat_analysis_withcovidcontrol_200_v2_plots/"

outFolder="./10_CellChat_analysis_withcovidcontrol_500_v2_plots/"
system(paste0("mkdir -p ", outFolder))


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


####################################
# Load data
####################################




##################################################################
# circle plots for paper 
##################################################################


# sapply(locations, function(xlocation){
#   
#   
#   
#   subFolder<-xlocation
#   
# 
#   cellchat_TNL<-cellchat_TNL_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_","Control","_","2021-12-21",".rds"))
#   cellchat_TIL<-cellchat_TIL_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_","E. coli","_","2021-12-21",".rds"))
#   # 
# 
#   pathways_TNL<-cellchat_TNL_original@netP$pathways
#   controls<-cbind(pathways_TNL, rep("control",length(pathways_TNL)))
#   controls<-as.data.frame(controls)
#   colnames(controls)<-c("pathway","control")
# 
#   # there is an additional population  30_B cell specific to cellchat_TNL compared to TIL
#   # we lift up TIL by lifting up the cell groups to the same cell labels as cellchat_TNL 
#   
#   TNL_cells<-rownames(cellchat_TNL@net$prob)
#   TIL_cells<-rownames(cellchat_TIL@net$prob)
#   
#   TIL_cells_only<-TIL_cells [!TIL_cells %in%TNL_cells ]
#   TNL_cells_only<-TNL_cells [! TNL_cells  %in%TIL_cells ]
#   
#   if (length(TIL_cells_only)==0)
#   {
#     group.new = levels(cellchat_TNL@idents)
#     cellchat_TIL <- liftCellChat(cellchat_TIL, group.new)
#   }else if (length(TNL_cells_only)==0){
#     group.new = levels(cellchat_TIL@idents)
#     cellchat_TNL <- liftCellChat(cellchat_TNL, group.new)} else 
#     {
#       group.new = union(levels( cellchat_TIL@idents), levels(cellchat_TNL@idents))
#       cellchat_TNL <- liftCellChat(cellchat_TNL, group.new)
#       cellchat_TIL <- liftCellChat(cellchat_TIL, group.new)
#     }
#   
#   
#   # now merge
#   object.list <- list(Control = cellchat_TNL, TIL = cellchat_TIL)
#   cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#   
#   
#     color.use.TIL=cluster.Colors[rownames(cellchat@netP$TIL$prob)]
#     color.use.Control=cluster.Colors[rownames(cellchat@netP$Control$prob)]
#     coloruses<-list(color.use.TIL,color.use.Control)
#   
#   
#     system(paste0("mkdir -p ", outFolder,subFolder,"/nolabel","/"))
#   
#     pdf(paste0(outFolder,subFolder,"/nolabel/","diffInteraction_",xlocation,".pdf"),width=15,height=15)
#     gg2<-netVisual_diffInteraction(cellchat, vertex.weight=6,vertex.label.cex=0.000001, weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)],arrow.width=15)
#     
#     gg2
#     dev.off()
#   
#     pdf(paste0(outFolder,subFolder,"/","diffInteraction_",xlocation,".pdf"),width=15,height=15)
#     gg2<-netVisual_diffInteraction(cellchat, vertex.weight=6, weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)],arrow.width=15)
#     gg2
#     dev.off()
#     
#     
#   
#   
#   
# })


locations<-c("CAM" ,    "PVBP")
conditions<-c("TIL" ,"TNL")

##################################################
#comparison between TIL and pbs plots
######################################################

sapply(locations, function(xlocation){
  
  
  
  subFolder<-xlocation
    
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  

  #cellchat_TNL<-cellchat_TNL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_200/","cellchat_",xlocation,"_","TNL","_","2022-04-07",".rds"))
  #cellchat_TIL<-cellchat_TIL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_200/","cellchat_",xlocation,"_","TIL","_","2022-04-07",".rds"))
 
  # cellchat_TNL<-cellchat_TNL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_200/","cellchat_",xlocation,"_","TNL","_","2022-05-05",".rds"))
  # cellchat_TIL<-cellchat_TIL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_200/","cellchat_",xlocation,"_","TIL","_","2022-05-05",".rds"))
  # 
  cellchat_TNL<-cellchat_TNL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_500/","cellchat_",xlocation,"_","TNL","_","2022-05-05",".rds"))
  cellchat_TIL<-cellchat_TIL_original<-read_rds(paste0("./10_CellChat_analysis_withcovidcontrol_v2_500/","cellchat_",xlocation,"_","TIL","_","2022-05-05",".rds"))
  
  
  pathways_TNL<-cellchat_TNL_original@netP$pathways
  controls<-cbind(pathways_TNL, rep("control",length(pathways_TNL)))
  controls<-as.data.frame(controls)
  colnames(controls)<-c("pathway","control")
  
  
  pathways_TIL<-cellchat_TIL_original@netP$pathways
  TILs<-cbind(pathways_TIL, rep("TIL",length(pathways_TIL)))
  TILs<-as.data.frame(TILs)
  colnames(TILs)<-c("pathway","TIL")
  
  allpathways<-TILs %>% full_join(controls)
  
  allpathways<-allpathways %>% arrange(desc(!is.na(control),(!is.na(TIL))))
  
  write.csv(allpathways,file=paste0(outFolder,subFolder,"/","pathways_",xlocation,".csv"))
  
  shared_pathways<-intersect(pathways_TNL,pathways_TIL)
  pathways_TNL_only<-pathways_TNL[!pathways_TNL %in% pathways_TIL]
  pathways_TIL_only<-pathways_TIL[! pathways_TIL %in% pathways_TNL]
  
  
  
  
  # there is an additional population  30_B cell specific to cellchat_TNL compared to TIL
  # we lift up TIL by lifting up the cell groups to the same cell labels as cellchat_TNL 
  
  TNL_cells<-rownames(cellchat_TNL@net$prob)
  TIL_cells<-rownames(cellchat_TIL@net$prob)
  
  TIL_cells_only<-TIL_cells [!TIL_cells %in%TNL_cells ]
  TNL_cells_only<-TNL_cells [! TNL_cells  %in%TIL_cells ]
  
  if (length(TIL_cells_only)==0)
  {
    group.new = levels(cellchat_TNL@idents)
    cellchat_TIL <- liftCellChat(cellchat_TIL, group.new)
  }else if (length(TNL_cells_only)==0){
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
  
  color.use.TIL=cluster.Colors[rownames(cellchat@netP$TIL$prob)]
  color.use.Control=cluster.Colors[rownames(cellchat@netP$TNL$prob)]
  coloruses<-list(color.use.Control,color.use.TIL)
  
  overalcoloruses<-cluster.Colors[unique(c(rownames(cellchat@netP$TIL$prob), rownames(cellchat@netP$TNL$prob)))]
  
  
#   
  pdf(paste0(outFolder,subFolder,"/","informationflow_",xlocation,".pdf"),width=10,height=18)
  #gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("#333399","#A50021"),font.size=15)
  gg2
  dev.off()
#   
#   
#   
#   
  # Compare the total number of interactions and interaction strength
  # pdf(paste0(outFolder,subFolder,"/","total_number_interactions_",xlocation,".pdf"),width=10,height=4)
  # gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use=c("#A50021","#333399"))
  # gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",color.use=c("#A50021","#333399"))
  # gg1 + gg2
  # dev.off()
#   
#   
#   # Compare the number of interactions and interaction strength among different cell populations
#   
#   # The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, 
#  # where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
#   # 
#   
  
  
  pdf(paste0(outFolder,subFolder,"/","diffInteraction_",xlocation,".pdf"),width=25,height=25)
  #par(mfrow = c(1,2), xpd=TRUE)
  #gg1<-netVisual_diffInteraction(cellchat, weight.scale = T,top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)])
  gg2<-netVisual_diffInteraction(cellchat, arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8,weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)],vertex.label.cex = 2.5,edge.width.max=20)
  #gg1 + gg2
  gg2
  dev.off()
  
  pdf(paste0(outFolder,subFolder,"/","diffInteraction_nolabel_",xlocation,".pdf"),width=25,height=25)
  #par(mfrow = c(1,2), xpd=TRUE)
  #gg1<-netVisual_diffInteraction(cellchat, weight.scale = T,top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)])
  gg2<-netVisual_diffInteraction(cellchat, arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8,weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)],vertex.label.cex = 0.00001,edge.width.max=20)
  #gg1 + gg2
  gg2
  dev.off()
  
  
  
  weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
  pdf(paste0(outFolder,subFolder,"/","interactions_strength_conditions_nolabel_",xlocation,".pdf"),width=26,height=15)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$weight, arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8, vertex.label.cex=0.000001,top=0.25,weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]),color.use = coloruses[[i]])
  }
  dev.off()
  
  
  weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
  pdf(paste0(outFolder,subFolder,"/","interactions_strength_conditions_",xlocation,".pdf"),width=26,height=15)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$weight,arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8,top=0.25,weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]),color.use = coloruses[[i]])
  }
  dev.off()
  
  
  # weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
  # pdf(paste0(outFolder,subFolder,"/","interactions_strength_conditions_",xlocation,".pdf"),width=25,height=20)
  # par(mfrow = c(1,2), xpd=TRUE)
  # for (i in 1:length(object.list)) {
  #   netVisual_circle(object.list[[i]]@net$weight, top=0.25,arrow.width=4,weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]),color.use = coloruses[[i]])
  # }
  # dev.off()
  
  
#   
#   
#   #We can also show differential number of interactions or interaction strength in a greater details using a heatmap. The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). In the colorbar, 
#   # red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.
#   
#   
  pdf(paste0(outFolder,subFolder,"/","diffInteraction_heatmap_",xlocation,".pdf"),width=10,height=10)
  #gg1 <- netVisual_heatmap(cellchat,color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)])
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.use=cluster.Colors[rownames(cellchat@netP$TIL$prob)])
  #> Do heatmap based on a merged object
  #gg1 + gg2
  gg2
  dev.off()
#   
#   
#   
#  # To better control the node size and edge weights of the inferred networks across different datasets, 
#   # we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.
#   
#   
#   
#   
#   
  color.use.TIL=cluster.Colors[rownames(cellchat@netP$TIL$prob)]
  color.use.Control=cluster.Colors[rownames(cellchat@netP$TNL$prob)]
  coloruses<-list(color.use.TIL,color.use.Control)


# vector plot 
  
pdf(paste0(outFolder,subFolder,"/","outgoing_incoming_conditions_",xlocation,".pdf"),width=19,height=10)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.listi <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")

  gg[[i]] <- netAnalysis_signalingRole_scatter(object.listi, title = names(object.list)[i], weight.MinMax = weight.MinMax,color.use = coloruses[[i]])
}
 p<-patchwork::wrap_plots(plots = gg)
 plot(p)
 dev.off()

 
 
 
# vector plot (differential)
label.size = 4
dot.size = c(2, 6)
dot.alpha = 0.6
font.size.title = 13
font.size = 15
df1<-gg[[1]]$data
df1$group<-names(object.list)[1]
df2<-gg[[2]]$data
df2$group<-names(object.list)[2]
newdf<-rbind(df1,df2)
newdf$color_custome<-cluster.Colors[newdf$labels]

newdf$labels<-as.character(newdf$labels)
newdf$color_custome<-cluster.Colors[newdf$labels]

xlabel = "Outgoing interaction strength"
ylabel = "Incoming interaction strength"


overalcoloruses<-cluster.Colors[unique(c(rownames(cellchat@netP$TNL$prob),rownames(cellchat@netP$TIL$prob)))]
overalcoloruses<-overalcoloruses[newdf$labels]

newdf<-newdf %>% arrange(labels)
newdf$labels <- factor(newdf$labels,levels=unique(newdf$labels))
#newdf$group <- factor(newdf$group,levels=unique(newdf$group))




# newdf$color_custome <- factor(newdf$color_custome,levels=unique(newdf$color_custome))
#overalcoloruses<-cluster.Colors[unique(c(rownames(cellchat@netP$TIL$prob), rownames(cellchat@netP$TNL$prob)))]
#overalcoloruses<-overalcoloruses[newdf$labels]

require(grid)
#newdf <- data.frame(x = outgoing.cells, y = incoming.cells, labels = names(incoming.cells), Count = num.link)
#newdf$labels <- factor(newdf$labels, levels = names(incoming.cells))
gg <- ggplot(data = newdf, aes(x, y),show.legend = F) + geom_point(aes(size = Count, colour = labels, fill = labels),show.legend = F)
gg <- gg + CellChat_theme_opts() + theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in")) + labs(title = "TIL + control", x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, face = "plain")) + theme(axis.line.x = element_line(size = 0.25),                                                                                                                                                                               axis.line.y = element_line(size = 0.25))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size)
gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = labels), size = label.size, show.legend = F, segment.size = 0.2, segment.alpha = 0.5)

gg<-gg+geom_path(aes(colour=labels,size = 1.5),arrow = arrow(type = "closed",length=unit(0.2, "inches")))

#gg<-gg+geom_line(aes(group = labels,colour=labels,size = 1.5),arrow = arrow(type = "closed",length=unit(0.2, "inches")))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_colour_manual(values =overalcoloruses, drop = FALSE) + guides(colour = FALSE)+ guides(fill = FALSE)

pdf(paste0(outFolder,subFolder,"/","outgoing_incoming_conditions_both_",xlocation,".pdf"),width=17,height=12)
gg+ theme(legend.position = "none")
dev.off()


require(grid)
#newdf <- data.frame(x = outgoing.cells, y = incoming.cells, labels = names(incoming.cells), Count = num.link)
#newdf$labels <- factor(newdf$labels, levels = names(incoming.cells))

gg <- ggplot(data = newdf, aes(x, y)) + geom_point(aes(size = Count, colour = labels, fill = labels),show.legend = FALSE)
gg <- gg + CellChat_theme_opts() + theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in")) + labs(title = "TIL + control", x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, face = "plain")) + theme(axis.line.x = element_line(size = 0.25),axis.line.y = element_line(size = 0.25))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size)+ guides(fill=FALSE, color=FALSE)
gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = labels), size = 0.0001, show.legend = F, segment.size = 0.2, segment.alpha = 0.5)

gg<-gg+geom_path(aes(colour=labels,size = 1.5),arrow = arrow(type = "closed",length=unit(0.2, "inches")))
#gg<-gg+geom_line(aes(group = labels,colour=labels,size = 1.5),arrow = arrow(type = "closed",length=unit(0.2, "inches")))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_colour_manual(values =overalcoloruses, drop = FALSE) + guides(colour = FALSE)+ guides(fill = FALSE)
#gg + theme(legend.position = "none")
pdf(paste0(outFolder,subFolder,"/","outgoing_incoming_conditions_nolabel_both_",xlocation,".pdf"),width=17,height=12)
gg+ theme(legend.position = "none")

dev.off()


})







