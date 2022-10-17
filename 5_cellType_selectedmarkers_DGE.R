library(Seurat)
library(Matrix)
library(tidyverse)

#############################################################
### cluster marker identification plot (selected markers)
### 
#############################################################

sc <- read_rds("4_harmony_with_covidcontrol_res1.0/sc.NormByLocationRep.Harmony.rds")


# new_names <- read_tsv("ClusterAssignmentSimple.tsv")
# clust2Names <- new_names$scLabor_ID
# names(clust2Names) <- new_names$seurat_clusters
# cc <- new_names %>% dplyr::select(scLabor_ID,color) %>% unique
# cluster.Colors <- cc$color
# names(cluster.Colors) <- cc$scLabor_ID
# sc@meta.data$cluster_name <- clust2Names[sc@meta.data$seurat_clusters]
# Idents(sc) <- "cluster_name"




outFolder="./5_harmonyClusters_withcovid19control_DGE_plots/"
system(paste0("mkdir -p ", outFolder))


fname=paste0(outFolder,"UMAP_Harmony_onlycontrols.png");
png(fname,width=1000,height=1000)
DimPlot(sc, cols = cluster.Colors, label=FALSE , repel=TRUE)
dev.off()


aa <- FetchData(sc,c("UMAP_1","UMAP_2","Location","Condition","Origin","status","FetalSex","seurat_clusters","cluster_name","Library","Pregnancy_ID")) 
head(aa)




myscale = 1/colSums(sc@assays$RNA@counts)*1000000

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


#anno <- read_rds("/nfs/rprdata/scilab/labor2/covid-analysis/2020-10-02/3_MergeDemux_Output/anno.rds")

#genemark <- read_csv("RogerMarkers.csv")

##genesym=c("PTPRC", "ELANE", "CD14", "CD19", "CD3E","CD4","CD8A","KLRC3")
##genesel <- filter(anno,gene_name %in% genesym)


#genesym<-c("FOXP3","IL17A","RORC","TH17")
genesym<-c("CD2", "CD25", "RORC", "FOXP3", "CD4", "CD3G", "CD3D", "CD8A", "IL7R", "CCR7", "S100A4", "GNLY", "NKG7", "PPBP", "FCER1A", "CST3", "CD14", "LYZ", "GH2", "MS4A1")
genesel <- filter(anno,gene_name %in% genesym)


#genesym<-genesym[genesym %in% rownames(sc)]

genexpr <- map_dfr(1:nrow(genesel),function(i){
                aux <- FetchData(sc,genesel$kbid[i])
                        colnames(aux) <- "expr"
                        aux <- cbind(aa,aux*myscale)
                        aux$gene=genesel$kbid[i]
                        aux$symbol=genesel$gene_name[i]
                        #aux$cell_name=genesel$cell_name[i]
                        aux
                    })


# genexpr <- map_dfr(1:length(genesym),function(i){
#     aux <- FetchData(sc,genesym[i])
#     colnames(aux) <- "expr"
#     aux <- cbind(aa,aux*myscale)
#     #aux$gene=genesel$kbid[i]
#     aux$symbol=genesym[i] #genesel$gene_name[i]
#     #aux$cell_name=genesel$cell_name[i]
#     aux
# })



#genexpr$cell_symbol <- paste0(genexpr$cell_name," (",genexpr$symbol,")")

fname=paste0(outFolder,"selectedmarkers_umap2.png");
##png(fname,width=1000,height=800)
p1 <- genexpr %>%  arrange(symbol,expr) %>%
    ggplot(aes(UMAP_1,UMAP_2,color=log10(0.1+expr))) +
    geom_point(size=0.1) +
    ##        scale_fill_viridis_c(option = "plasma") +
    scale_color_gradient(low = "lightblue", high = "darkred") +
    facet_wrap(~symbol) +
    theme_bw()
ggsave(fname,p1,width=11,height=6)
    ##dev.off()                                                                                                                                               cat(fname,"\n")       



sum_rec <- genexpr %>% dplyr::group_by(symbol,seurat_clusters) %>% dplyr::summarize(Prop=mean(expr>0),Expr=mean(expr[expr>0]))

fname=paste0(outFolder,"selectedmarkers_dotplot.png");
p0 <- ggplot(sum_rec,aes(x=seurat_clusters,y=symbol,color=Expr,size=Prop)) +
    geom_point() +
    scale_size_area(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),max_size=6,na.value=0) +
    scale_color_distiller(trans='log10', palette = "RdYlBu",direction = -1,na.value=0,limits=c(10,NA)) +
    ##facet_grid( .~ symbol) +
    theme_classic() +
    labs(y="Gene Symbol",x="Cell type",size="Prop.",color="Expr.(TPM)") + 
    theme(axis.text.x = element_text(angle = 45,hjust=1)) 
ggsave(fname,p0,width=10,height=12)
##dev.off()
cat(fname,"\n")



# fname=paste0(outFolder,"CAM_","_umap3.png");
# ##png(fname,width=1000,height=800)
# p1 <- genexpr %>% filter(Location=='CAM') %>%  arrange(symbol,expr) %>%
#     ggplot(aes(UMAP_1,UMAP_2,color=log10(0.1+expr))) +
#     geom_point(size=0.1) +
#     ##        scale_fill_viridis_c(option = "plasma") +
#     scale_color_gradient(low = "lightblue", high = "darkred") +
#     facet_wrap(~cell_symbol) +
#     theme_bw()
# ggsave(fname,p1,width=15,height=12)
# 
# fname=paste0(outFolder,"PVBP_","_umap3.png");
# ##png(fname,width=1000,height=800)
# p1 <- genexpr %>% filter(Location=='PVBP') %>%  arrange(symbol,expr) %>%
#     ggplot(aes(UMAP_1,UMAP_2,color=log10(0.1+expr))) +
#     geom_point(size=0.1) +
#     ##        scale_fill_viridis_c(option = "plasma") +
#     scale_color_gradient(low = "lightblue", high = "darkred") +
#     facet_wrap(~cell_symbol) +
#     theme_bw()
# ggsave(fname,p1,width=15,height=12)
# 
