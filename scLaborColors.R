newNames = c("B_cells"="B-cell"
            ,"ColumnCyto"="npiCTB"
            ,"Cytotroph"="CTB"
            ,"Decidual"="Decidual"
            ,"Dendr_Macro_A"="Macrophage-1"
            ,"Dendr_Macro_B"="Macrophage-2"
            ,"Endometrial"="Endometrial"
            ,"Endothelial"="LED"
            ,"EVT"="EVT"
            ,"Fibroblasts"="Fibroblast"
            ,"HSC"="HSC"
            ,"Monocytes"="Monocyte"
            ,"(Myeloid)Progenitor"="Stromal-3"
            ,"NK_cells"="NK-cell"
            ,"Stromal_A"="Stromal-1"
            ,"Stromal_B"="Stromal-2"
            ,"Synciotrophoblasts"="STB"
            ,"Tcells_activated"="T-cell-activated"
            ,"Tcells_resting"="T-cell-resting")

group.colors<-c(B_cells="red",ColumnCyto="limegreen",Cytotroph="purple4",Decidual="sienna4",Dendr_Macro_A="hotpink",Dendr_Macro_B="tomato",Endometrial="maroon",Endothelial="yellow4",EVT="plum2",Fibroblasts="navy",HSC="magenta",Monocytes="lightslateblue",Progenitor="peru",NK_cells="royalblue",Stromal_A="turquoise4",Stromal_B="lightpink3",Synciotrophoblasts="yellowgreen",Tcells_activated="seagreen",Tcells_resting="powderblue")


names(group.colors) <- newNames 

