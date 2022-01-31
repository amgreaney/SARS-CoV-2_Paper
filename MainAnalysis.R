# "SARS-CoV-2 leverages airway epithelial protective mechanism for viral infection"

# Original code used to analyze data and generate figures
# for submission to Cell Systems


## login & establish session on cluster

# conda create -n velocytomonocle -c conda-forge \
#                        -c bioconda \
#                        -c eugene_t \
#                        r-ggplot2 \
#                        r-seurat \
#                        r-seuratwrappers \
#                        r-velocyto.r \
#                        r-pagoda2 \
#                        bioconductor-monocle \
#                        r-monocle3 \
#                        bioconductor-complexheatmap \
#                        r-viridis \
#                        r-rcolorbrewer \
#                        r-loomR \
#                        r-hdf5r

# conda activate ##/velocytomonocle

## Enter R
R

library(Seurat)
library(ggplot2)
library(scales)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(circlize)
#library(monocle)
library(monocle3)
#library(loomR)
library(ComplexHeatmap)
library(cowplot)
library(Matrix)
setwd("##")
load('./hbec.integrated2020-10-09.Robj')

# pull and make seurat object w relevant metadata
setwd('##')
hbec <- ReadH5AD('./scv2_200428.h5ad')
metadata <- read.csv('./hbec_metadata_scv2counts.csv')

hbec <- AddMetaData(hbec, metadata = metadata$ctype, col.name = 'ctype')
hbec <- AddMetaData(hbec, metadata = metadata$Condition, col.name = 'Condition')
hbec <- AddMetaData(hbec, metadata = metadata$batch, col.name = 'batch')

hbec <- SplitObject(hbec, split.by = "Condition")

setwd("##")
save(hbec, file = paste("hbec",Sys.Date(),".Robj",sep=""))
load('./hbec2020-05-15.Robj')

mock <- hbec$Mock
DefaultAssay(mock) <- 'RNA'
dpi1 <- hbec$'1dpi'
DefaultAssay(dpi1) <- 'RNA'
dpi2 <- hbec$'2dpi'
DefaultAssay(dpi2) <- 'RNA'
dpi3 <- hbec$'3dpi'
DefaultAssay(dpi3) <- 'RNA'
rm(hbec)

## Perform integration across infection timepoints
# Determine genes to use for CCA, must be highly variable in at least 2 datasets
ob.list <- list(mock, dpi1, dpi2, dpi3)
hbec.anchors <- FindIntegrationAnchors(ob.list, dims = 1:25)
rm(ob.list,mock,dpi1,dpi2,dpi3)
hbec.integrated <- IntegrateData(anchorset = hbec.anchors, dims = 1:25)
DefaultAssay(hbec.integrated) <- "RNA"
rm(hbec.anchors)

save(hbec.integrated, file = paste("hbec.integrated.prescale",Sys.Date(),".Robj",sep=""))
load('./hbec.integrated.prescale2020-05-15.Robj')

# Run the standard workflow for visualization and clustering
ls()
set.seed(2)
hbec.integrated <- NormalizeData(hbec.integrated)
hbec.integrated <- ScaleData(hbec.integrated,vars.to.regress = c('pmito','nCount_RNA'))
hbec.integrated <- RunPCA(hbec.integrated, npcs = 30, verbose = FALSE)
hbec.integrated <- RunUMAP(hbec.integrated, reduction = "pca", dims = 1:30)
hbec.integrated <- FindNeighbors(hbec.integrated,k.param=4,dims = 1:30)
hbec.integrated <- FindClusters(hbec.integrated,resolution = 1)
png('hbec.integrated_umap_1.png', width = 800, height = 800)
UMAPPlot(hbec.integrated,label = T)
dev.off()

png('hbec.integrated_umap_2.png', width = 800, height = 800)
UMAPPlot(hbec.integrated,label = F, group.by = 'Condition')
dev.off()

png('hbec.integrated_umap_3.png', width = 800, height = 800)
UMAPPlot(hbec.integrated,group.by='Condition',split.by='Condition',ncol=2)
dev.off()

png('hbec.integrated_umap_4.png', width = 800, height = 800)
UMAPPlot(hbec.integrated,label = F, group.by = 'ctype')
dev.off()

png('hbec.integrated_fp_1.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('KRT5','SCGB1A1','PIFO','POU2F3'), label = T)
dev.off()

png('hbec.integrated_fp_2.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('TOP2A','TP63','LYPD2','MUC5AC'), label = T)
dev.off()

png('hbec.integrated_fp_3.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('FOXJ1','GNAT3','FOXI1','CFTR'), label = T)
dev.off()

png('hbec.integrated_fp_4.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('ASCL1','CALCA','SOX9','KRT14'), label = T)
dev.off()

png('hbec.integrated_fp_5.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('FN1','TAGLN','VIM','COL1A1'), label = T)
dev.off()

png('hbec.integrated_fp_6.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('MKI67','KRT5','FN1'), label = F)
dev.off()

png('hbec.integrated_fp_7.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('SCGB3A1','MUC5AC','PIFO'), label = F)
dev.off()

png('hbec.integrated_fp_8.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('POU2F3','FOXI1','ASCL1'), label = F)
dev.off()

hbec.integrated.marks <- FindAllMarkers(hbec.integrated,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
hbec.integrated.marks$ratio <- hbec.integrated.marks$pct.1/hbec.integrated.marks$pct.2
save(hbec.integrated.marks,file = paste("hbec.integrated.marks",Sys.Date(),".Robj",sep=""))
write.table(hbec.integrated.marks,file = paste("hbec.integrated.marks",Sys.Date(),".txt",sep=""),sep="\t")

# isolating goblet and tuft cells manually
png('hbec.integrated_fp_c23.png', width = 800, height = 800)
FeaturePlot(subset(hbec.integrated,idents = '23'), c('MUC5AC'))
dev.off()
# subset cluster 23
hbec23 <- subset(hbec.integrated, idents=23)
hbec23.goblet <- subset(hbec23, subset = MUC5AC > 1)
png('hbec23_umap_1.png', width = 800, height = 800)
UMAPPlot(hbec23.goblet,label = T)
dev.off()
# create new column in metadata
hbec.integrated$sub_cluster <- as.character(Idents(hbec.integrated))
# Change the information of cells containing sub-cluster information
hbec.integrated$sub_cluster[Cells(hbec23.goblet)] <- paste("29",Idents(hbec23.goblet))
png('hbec23_umap_2.png', width = 800, height = 800)
DimPlot(hbec.integrated, group.by = "sub_cluster")
dev.off()

# isolating pnec manually
png('hbec.integrated_fp_c20.png', width = 800, height = 800)
FeaturePlot(subset(hbec.integrated,idents = '20'), c('PIFO','ASCL1'))
dev.off()
# subset cluster 20
hbec20 <- subset(hbec.integrated, idents=20)
hbec20.cilia <- subset(hbec20, subset = PIFO > 2)
png('hbec20_umap_1.png', width = 800, height = 800)
UMAPPlot(hbec20.cilia,label = T)
dev.off()
# already have sub_cluster slot in metadata (see above)
# Change the information of cells containing sub-cluster information
hbec.integrated$sub_cluster[Cells(hbec20.cilia)] <- paste("30",Idents(hbec20.cilia))
png('hbec20_umap_2.png', width = 800, height = 800)
DimPlot(hbec.integrated, group.by = "sub_cluster")
dev.off()

hbec.integrated[['initcluster2']] <- Idents(hbec.integrated)

Idents(hbec.integrated) <- hbec.integrated[['initcluster2']]

hbec.integrated <- RenameIdents(hbec.integrated,
                      '0'='Intermediate1',
                      '1'='Club1',
                      '2'='Ciliated1',
                      '3'='Basal1',
                      '4'='Club2',
                      '5'='Basal2',
                      '6'='Ciliated2',
                      '7'='Ciliated3',
                      '8'='Basal3',
                      '9'='Club3',
                      '10'='Intermediate2',
                      '11'='Ciliated4',
                      '12'='Intermediate3',
                      '13'='Club4',
                      '14'='CyclingBasal',
                      '15'='Ciliated5',
                      '16'='Ionocytes',
                      '17'='Ciliated6',
                      '18'='Ciliated7',
                      '19'='Club5',
                      '20'='PNECs',
                      '21'='Ciliated8',
                      '22'='Club6',
                      '23'='Tuft',
                      '24'='BasalEMT',
                      '25'='Intermediate4',
                      '26'='Ciliated9',
                      '27'='Ciliated10',
                      '28'='Basal5',
                      '29 23'='Goblet',
                      '30 20'='Ciliated11')
hbec.integrated[['cell_type1']] <- Idents(hbec.integrated)

png('hbec.integrated_umap_5.png', width = 800, height = 800)
UMAPPlot(hbec.integrated,label = F, group.by = 'cell_type1')
dev.off()

Idents(hbec.integrated) <- hbec.integrated[['initcluster2']]

hbec.integrated <- RenameIdents(hbec.integrated,
                      '0'='Intermediate',
                      '1'='Club',
                      '2'='Ciliated',
                      '3'='Basal',
                      '4'='Club',
                      '5'='Basal',
                      '6'='Ciliated',
                      '7'='Ciliated',
                      '8'='Basal',
                      '9'='Club',
                      '10'='Intermediate',
                      '11'='Ciliated',
                      '12'='Intermediate',
                      '13'='Club',
                      '14'='CyclingBasal',
                      '15'='Ciliated',
                      '16'='Ionocytes',
                      '17'='Ciliated',
                      '18'='Ciliated',
                      '19'='Club',
                      '20'='PNECs',
                      '21'='Ciliated',
                      '22'='Club',
                      '23'='Tuft',
                      '24'='BasalEMT',
                      '25'='Intermediate',
                      '26'='Ciliated',
                      '27'='Ciliated',
                      '28'='Basal',
                      '29 23'='Goblet',
                      '30 20'='Ciliated')
hbec.integrated[['cell_type2']] <- Idents(hbec.integrated)

png('hbec.integrated_umap_7.png', width = 800, height = 800)
UMAPPlot(hbec.integrated,label = F, group.by = 'cell_type2')
dev.off()

Idents(hbec.integrated) <- hbec.integrated[['cell_type4']]
hbec.integrated <- subset(hbec.integrated, idents = 'Remove',invert = T)

# FIGURE PLOTS
# CyclingBasal, Basal, BasalEMT, Intermediate, Club, Goblet, Ciliated, Tuft, Ionocytes, PNECs
celltypecolors <- sample(colorRampPalette(RColorBrewer::brewer.pal(n=8, name='Set1'))(10))
celltypecolors <- c("#F37912", "#D7B32E","#5D6795", "#658E67", "#A35390", "#FFD422", "#E41A1C", "#F781BF", "#43997A", "#B85F49")
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
Idents(hbec.integrated) <- factor(x=Idents(hbec.integrated), levels = c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))
png('hbec.integrated_umap_82.png', width = 600, height = 800)
UMAPPlot(hbec.integrated,label = T,label.size=8,repel=T,cols = celltypecolors) + labs(title = 'Cell Type') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoLegend() + NoAxes()
dev.off()

cellcolorssplit <- c("#D7B32E",'#66C2A5',"#A35390","#658E67",'#E78AC3','#8DA0CB','#FC8D62')
Idents(hbec.integrated.sub) <- hbec.integrated.sub[['cell_type4']]
Idents(hbec.integrated.sub) <- factor(x=Idents(hbec.integrated.sub), levels = c('Basal','Ciliated_Progenitor','Club','Intermediate','Mature_Ciliated1','Mature_Ciliated2','Novel_Infected_Ciliated'))
png('hbec.integrated_umap_83.png', width = 600, height = 800)
UMAPPlot(hbec.integrated.sub,label = F,label.size=8,repel=T,cols = cellcolorssplit) + labs(title = 'Connectomic Cells') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoLegend() + NoAxes()
dev.off()

#Figure 3
png('hbec.integrated_umap_84.png', width = 1200, height = 400)
UMAPPlot(hbec.integrated.sub,cols = cellcolorssplit,split.by='Condition') + theme(plot.title = element_text(size = 38, hjust = 0.5))
dev.off()

cellcolorssplit2 <- c('#D7B32E','#A35390','#658E67','#E41A1C')
Idents(hbec.integrated.sub) <- hbec.integrated.sub[['cell_type2']]
Idents(hbec.integrated.sub) <- factor(x=Idents(hbec.integrated.sub), levels = c('Basal','Club','Intermediate','Ciliated'))
png('hbec.integrated_umap_85.png', width = 1200, height = 400)
UMAPPlot(hbec.integrated.sub,cols = cellcolorssplit2,split.by='Condition') + theme(plot.title = element_text(size = 38, hjust = 0.5))
dev.off()

png('hbec.integrated_umap_9.png', width = 650, height = 800)
UMAPPlot(hbec.integrated,group.by = 'Condition') + labs(title = 'Time Point') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
                                                theme(legend.text = element_text(size = 22),
                                                legend.key.size = unit(1.5, "lines"))
dev.off()

Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
bigcil <- subset(hbec.integrated,idents = c('Ciliated'))
Idents(bigcil) <- bigcil[['ciliated3']]
Idents(bigcil) <- factor(x=Idents(bigcil), levels = c('Mature_Ciliated1','Ciliated_Progenitor','Mature_Ciliated2','Novel_Infected_Ciliated'))
png('hbec.integrated_umap_cil2.png', width = 600, height = 800)
UMAPPlot(bigcil,cols=cilcolors,label = F) + labs(title = 'Ciliated Subtype') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoLegend() + NoAxes()
dev.off()

png('hbec.integrated_fp_bigcil.png', width = 600, height = 800)
FeaturePlot(bigcil,c('SCGB3A1','SCGB1A1','LYPD2','BPIFB1'))
dev.off()

save(hbec.integrated, file = paste("hbec.integrated",Sys.Date(),".Robj",sep=""))
load('./hbec.integrated2020-10-09.Robj')

Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec.integrated.marks <- FindAllMarkers(hbec.integrated,only.pos = T)
hbec.integrated.marks$ratio <- hbec.integrated.marks$pct.1/hbec.integrated.marks$pct.2
save(hbec.integrated.marks,file = paste("hbec.integrated.marks2",Sys.Date(),".Robj",sep=""))
write.table(hbec.integrated.marks,file = paste("hbec.integrated.marks2",Sys.Date(),".txt",sep=""),sep="\t")

png('hbec.integrated_fp_b.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('MKI67','KRT14','KRT5','IGFBP6'), label = T)
dev.off()

png('hbec.integrated_fp_emt.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('MMP9','FN1','VIM'), label = F)
dev.off()

png('hbec.integrated_fp_cl.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('SCGB3A1','SCGB1A1','LYPD2','BPIFB1'), label = T)
dev.off()

png('hbec.integrated_fp_gt.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('AGR2','MUC5AC','POU2F3','ASCL2'), label = F)
dev.off()

png('hbec.integrated_fp_c.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('FOXJ1','DNAH6','PIFO','TUBA1A'), label = T)
dev.off()

png('hbec.integrated_fp_i.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('FOXI1','CFTR','ASCL3','POSTN'), label = F)
dev.off()

png('hbec.integrated_fp_p.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('CHGA','ASCL1','CALCA'), label = F)
dev.off()

png('hbec.integrated_fp_int.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('SPRR1B','LYPD3','CLCA2','NOTCH3'), label = T)
dev.off()




## Cell identity dot plot
# gene list from Ravindra preprint
b.genes <- c('MKI67','TOP2A','KRT5','DAPL1','TP63')  # Basal
bc.genes <- c('CDK1','KRT4','KRT13')  # BC/Club
cl.genes <- c('SCGB3A1','SCGB1A1','LYPD2','BPIFB1')  # Club
g.genes <- c('MUC5AC','MUC5B','GP2')  # Goblet
c.genes <- c('FOXJ1','CCDC153','CCDC113','MLF1','LZTFL1')  # Ciliated
t.genes <- c('POU2F3','ASCL2')  # Tuft
i.genes <- c('FOXI1','CFTR','ASCL3')  # Ionocyte
p.genes <- c('CHGA','ASCL1')  # PNEC

genes.use <- rev(c(b.genes,bc.genes,cl.genes,g.genes,c.genes,t.genes,i.genes,p.genes))

# OR my gene list
b.genes <- c('MKI67','TOP2A','KRT5','KRT14','IGFBP6')  # CyclingBasal + Basal
emt.genes <- c('MMP9','FN1','VIM')  # BasalEMT
cl.genes <- c('SCGB3A1','SCGB1A1','LYPD2','BPIFB1')  # Club
g.genes <- c('AGR2','MUC5AC')  # Goblet
c.genes <- c('FOXJ1','DNAH6','PIFO')  # Ciliated
t.genes <- c('POU2F3','ASCL2')  # Tuft
i.genes <- c('FOXI1','CFTR','ASCL3','POSTN')  # Ionocyte
p.genes <- c('CHGA','ASCL1','CALCA')  # PNEC

genes.use <- rev(c(b.genes,emt.genes,cl.genes,g.genes,c.genes,t.genes,i.genes,p.genes))

# once you load the appropriate gene list...
DefaultAssay(hbec.integrated) <- 'RNA'
hbec.integrated <- ScaleData(hbec.integrated,vars.to.regress = c('pmito','nCount_RNA'),features = genes.use)

# DotPlot from my labels ('cell_type2')
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
Idents(hbec.integrated) <- factor(x=Idents(hbec.integrated), levels = rev(c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs')))

group_colors <- colorRampPalette(RColorBrewer::brewer.pal(n=10, name='Blues'))(2)

dot <- DotPlot(hbec.integrated,features = genes.use,assay='RNA',
            cols = group_colors,scale = T,dot.scale = 14) +
            theme(plot.title = element_text(size=50,color = 'black',hjust = 0.5),
                    axis.text.x=element_text(size=32, color = 'black',angle=90),
                    axis.text.y = element_text(size = 42, color = 'black'),
                    legend.text = element_text(size = 30),
                    legend.key.size = unit(20, "mm"),
                    legend.title = element_text(size=30),
                    strip.text.x = element_text(size=44),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(color = "black"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank())

png('hbec.integrated_dot_1.png', width = 1400, height = 1100)
dot
dev.off()

# DotPlot from their labels ('ctype')
Idents(hbec.integrated) <- hbec.integrated[['ctype']]
Idents(hbec.integrated) <- factor(x=Idents(hbec.integrated), levels = rev(c('Basal cells','BC/Club','Club cells','Goblet cells','Ciliated cells','Tuft cells','Ionocytes','Neuroendocrine cells')))

dot2 <- DotPlot(hbec.integrated,features = genes.use,assay='RNA',
            cols = group_colors,scale = T,dot.scale = 14) +
            theme(plot.title = element_text(size=50,color = 'black',hjust = 0.5),
                    axis.text.x=element_text(size=32, color = 'black',angle=90),
                    axis.text.y = element_text(size = 42, color = 'black'),
                    legend.text = element_text(size = 30),
                    legend.key.size = unit(20, "mm"),
                    legend.title = element_text(size=30),
                    strip.text.x = element_text(size=44),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(color = "black"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank())

png('hbec.integrated_dot_2.png', width = 1500, height = 1100)
dot2
dev.off()





## Generating gene lists for GO
# Bulk differential expression by timepoint w Mock
load('./hbec.integrated2020-05-20.Robj')
DefaultAssay(hbec.integrated) <- 'RNA'

Idents(hbec.integrated) <- hbec.integrated[['Condition']]

hbec.int.Condwmock.marks <- FindAllMarkers(hbec.integrated,only.pos = T)  # min.pct = 0.5,logfc.threshold = 0.5
hbec.int.Condwmock.marks$ratio <- hbec.int.Condwmock.marks$pct.1/hbec.int.Condwmock.marks$pct.2
save(hbec.int.Condwmock.marks,file = paste("hbec.int.Condwmock.marks",Sys.Date(),".Robj",sep=""))
write.table(hbec.int.Condwmock.marks,file = paste("hbec.int.Condwmock.marks",Sys.Date(),".txt",sep=""),sep="\t")



# Cell type over time - CILIATED including mock, 1dpi, 2dpi, 3dpi
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec.ciliated <- subset(hbec.integrated, idents = 'Ciliated')

Idents(hbec.ciliated) <- hbec.ciliated[['Condition']]

hbec.ciliated.marks <- FindAllMarkers(hbec.ciliated,only.pos = T)
hbec.ciliated.marks$ratio <- hbec.ciliated.marks$pct.1/hbec.ciliated.marks$pct.2
save(hbec.ciliated.marks,file = paste("hbec.ciliated.marks2",Sys.Date(),".Robj",sep=""))
write.table(hbec.ciliated.marks,file = paste("hbec.ciliated.marks2",Sys.Date(),".txt",sep=""),sep="\t")


# Cell type over time - CLUB
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec.club <- subset(hbec.integrated, idents = 'Club')

Idents(hbec.club) <- hbec.club[['Condition']]

hbec.club.marks <- FindAllMarkers(hbec.club,only.pos = T)
hbec.club.marks$ratio <- hbec.club.marks$pct.1/hbec.club.marks$pct.2
save(hbec.club.marks,file = paste("hbec.club.marks",Sys.Date(),".Robj",sep=""))
write.table(hbec.club.marks,file = paste("hbec.club.marks",Sys.Date(),".txt",sep=""),sep="\t")


# Cell type over time - BASAL
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec.basal <- subset(hbec.integrated, idents = 'Basal')

Idents(hbec.basal) <- hbec.basal[['Condition']]

hbec.basal.marks <- FindAllMarkers(hbec.basal,only.pos = T)
hbec.basal.marks$ratio <- hbec.basal.marks$pct.1/hbec.basal.marks$pct.2
save(hbec.basal.marks,file = paste("hbec.basal.marks",Sys.Date(),".Robj",sep=""))
write.table(hbec.basal.marks,file = paste("hbec.basal.marks",Sys.Date(),".txt",sep=""),sep="\t")


# Cell type over time - IONOCYTES
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec.ionocytes <- subset(hbec.integrated, idents = 'Ionocytes')

Idents(hbec.ionocytes) <- hbec.ionocytes[['Condition']]

hbec.ionocytes.marks <- FindAllMarkers(hbec.ionocytes,only.pos = T)
hbec.ionocytes.marks$ratio <- hbec.ionocytes.marks$pct.1/hbec.ionocytes.marks$pct.2
save(hbec.ionocytes.marks,file = paste("hbec.ionocytes.marks",Sys.Date(),".Robj",sep=""))
write.table(hbec.ionocytes.marks,file = paste("hbec.ionocytes.marks",Sys.Date(),".txt",sep=""),sep="\t")


# Cell type over time - PNECs
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec.pnecs <- subset(hbec.integrated, idents = 'PNECs')

Idents(hbec.pnecs) <- hbec.pnecs[['Condition']]

hbec.pnecs.marks <- FindAllMarkers(hbec.pnecs,only.pos = T)
hbec.pnecs.marks$ratio <- hbec.pnecs.marks$pct.1/hbec.pnecs.marks$pct.2
save(hbec.pnecs.marks,file = paste("hbec.pnecs.marks",Sys.Date(),".Robj",sep=""))
write.table(hbec.pnecs.marks,file = paste("hbec.pnecs.marks",Sys.Date(),".txt",sep=""),sep="\t")


# Cell type over time - TUFT
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec.tufts <- subset(hbec.integrated, idents = 'Tuft')

Idents(hbec.tufts) <- hbec.tufts[['Condition']]

hbec.tufts.marks <- FindAllMarkers(hbec.tufts,only.pos = T)
hbec.tufts.marks$ratio <- hbec.tufts.marks$pct.1/hbec.tufts.marks$pct.2
save(hbec.tufts.marks,file = paste("hbec.tufts.marks",Sys.Date(),".Robj",sep=""))
write.table(hbec.tufts.marks,file = paste("hbec.tufts.marks",Sys.Date(),".txt",sep=""),sep="\t")


# Cell type over time - INTERMEDIATE
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec.inter <- subset(hbec.integrated, idents = 'Intermediate')

Idents(hbec.inter) <- hbec.inter[['Condition']]

hbec.inter.marks <- FindAllMarkers(hbec.inter,only.pos = T)
hbec.inter.marks$ratio <- hbec.inter.marks$pct.1/hbec.inter.marks$pct.2
save(hbec.inter.marks,file = paste("hbec.inter.marks",Sys.Date(),".Robj",sep=""))
write.table(hbec.inter.marks,file = paste("hbec.inter.marks",Sys.Date(),".txt",sep=""),sep="\t")



## Population Bar Plots # FIGURE PLOT
# Representing cell type breakdown of time points
hbec.integrated$Condition <- factor(hbec.integrated$Condition,levels=c('Mock','1dpi','2dpi','3dpi'))
hbec.integrated$cell_type2 <- factor(hbec.integrated$cell_type2,levels=c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))

data = table(hbec.integrated$cell_type2,hbec.integrated$Condition)
# > data
#               Mock 1dpi 2dpi 3dpi
#  CyclingBasal  258  115  127 1203
#  Basal        5232  917 1257 8202
#  BasalEMT       73   67   66  103
#  Intermediate 3962 1901 2798 4881
#  Club        6134 3750 4529 5859
#  Goblet         12    2    8   27
#  Ciliated     6467 4770 5394 7283
#  Tuft          146   24   69  218
#  Ionocytes     194  152  315  548
#  PNECs         130  131  127  199

prop <- t(t(data)/rowSums(t(data))*100)
prop <- round(prop,1)

# columns should sum to ~100% bc we are breaking down cell type proportions of each time point
sum(prop[,1])
sum(prop[,2])
sum(prop[,3])
sum(prop[,4])

# Organizing data for the bar graph
sample <- c(rep('Mock',10),rep('1dpi',10),rep('2dpi',10),rep('3dpi',10))
cluster <- rep(c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'),4)
value <- as.vector(as.matrix(prop))
dataplot <- data.frame(sample,cluster,value)

# Set the order of your data
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))

# Run/save the ggplot bar graph
png('hbec.integrated_bar_1.png', width = 1200,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = T) + xlab('') + ylab('% of Cell Type per Time Point') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# Bar graph of high proportions - Basal, Intermediate, Club, Ciliated
# seuratcolors <- hue_pal()(10) run in R locally
# > seuratcolors
# [1] "#F8766D" "#D89000" "#A3A500" "#39B600" "#00BF7D" "#00BFC4" "#00B0F6" "#9590FF" "#E76BF3" "#FF62BC"
seuratcolors <- c('#F8766D','#D89000','#A3A500','#39B600','#00BF7D','#00BFC4','#00B0F6','#9590FF','#E76BF3','#FF62BC')

# Selecting just what we want
prop_high <- prop[c('Basal','Intermediate','Club','Ciliated'),]
# > prop_high
#               Mock 1dpi 2dpi 3dpi
#  Basal        23.1  7.8  8.6 28.8
#  Intermediate 17.5 16.1 19.0 17.1
#  Club         27.1 31.7 30.8 20.5
#  Ciliated     28.6 40.3 36.7 25.5

# Organizing data for the bar graph
sample <- c(rep('Mock',4),rep('1dpi',4),rep('2dpi',4),rep('3dpi',4))
cluster <- rep(c('Basal','Intermediate','Club','Ciliated'),4)
value <- as.vector(as.matrix(prop_high))
dataplot <- data.frame(sample,cluster,value)

# Set the order of your data
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('Basal','Intermediate','Club','Ciliated'))
#seuratcolors_high <- rep(c('#D89000','#39B600','#00BF7D','#00B0F6'),4)
celltypecolors_high <- rep(c('#D7B32E','#658E67','#A35390','#E41A1C'),4)

# Run/save the ggplot bar graph
png('hbec.integrated_bar_7.png', width = 900,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = F) +
  scale_fill_manual("legend", values = c("Basal" = "#D7B32E",
                                       "Intermediate" = "#658E67",
                                       "Club" = "#A35390",
                                       'Ciliated' = '#E41A1C')) +
  xlab('') + ylab('% of Cell Type per Time Point') +
  theme(axis.text.x=element_text(size=35,color = 'black'),
        axis.text.y = element_text(size = 35,color = 'black'),
        axis.title.y = element_text(size = 42,color = 'black'),
        legend.title = element_text(size = 0,color = 'black'),
        legend.text = element_text(size = 20,color = 'black'),
        legend.key.size = unit(3, "lines"))
dev.off()



# Bar graph of low proportions - CyclingBasal, BasalEMT, Goblet, Tuft, Ionocytes, PNECs
# Selecting just what we want
prop_low <- prop[c('CyclingBasal','BasalEMT','Goblet','Tuft','Ionocytes','PNECs'),]
# > prop_low
#               Mock 1dpi 2dpi 3dpi
#  CyclingBasal  1.1  1.0  0.9  4.2
#  BasalEMT      0.3  0.6  0.4  0.4
#  Goblet        0.1  0.0  0.1  0.1
#  Tuft          0.6  0.2  0.5  0.8
#  Ionocytes     0.9  1.3  2.1  1.9
#  PNECs         0.6  1.1  0.9  0.7

# Organizing data for the bar graph
sample <- c(rep('Mock',6),rep('1dpi',6),rep('2dpi',6),rep('3dpi',6))
cluster <- rep(c('CyclingBasal','BasalEMT','Goblet','Tuft','Ionocytes','PNECs'),4)
value <- as.vector(as.matrix(prop_low))
dataplot <- data.frame(sample,cluster,value)

# Set the order of your data
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('CyclingBasal','BasalEMT','Goblet','Tuft','Ionocytes','PNECs'))

# Run/save the ggplot bar graph
png('hbec.integrated_bar_7.png', width = 900,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = F) +
  scale_fill_manual("legend", values = c("CyclingBasal" = "#FFD422",
                                         "BasalEMT" = "#B85F49",
                                         "Goblet" = "#D7B32E",
                                         'Tuft' = '#F781BF',
                                         'Ionocytes' = '#43997A',
                                         'PNECs' = '#5D6795')) +
  xlab('') + ylab('% of Cell Type per Time Point') +
  theme(axis.text.x=element_text(size=35,color = 'black'),
        axis.text.y = element_text(size = 35,color = 'black'),
        axis.title.y = element_text(size = 42,color = 'black'),
        legend.title = element_text(size = 0,color = 'black'),
        legend.text = element_text(size = 20,color = 'black'),
        legend.key.size = unit(3, "lines"))
dev.off()

## Evidenced violin plots from GO Terms
Idents(hbec.integrated) <- hbec.integrated[['Condition']]
DefaultAssay(hbec.integrated) <- 'RNA'
png('hbec.integrated_vln_muc2.png', width = 1200,height=700)
VlnPlot(hbec.integrated,c('MUC1','MUC4','MUC16','MUC5B'),ncol = 2, pt.size = 0)
dev.off()

png('hbec.integrated_vln_hla.png', width = 1200,height=700)
VlnPlot(hbec.integrated,c('HLA-C','B2M','PSMB8','PSMA7'),ncol = 2, pt.size = 0)
dev.off()

dpi1 <- subset(hbec.integrated,idents = c('1dpi'))
Idents(dpi1) <- dpi1[['cell_type2']]
png('hbec.integrated_vln_hla1dpi2.png', width = 800,height=400)
VlnPlot(dpi1,c('HLA-C','HLA-A'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_vln_tlr1dpi.png', width = 800,height=400)
VlnPlot(dpi1,c('TLR1','TLR2'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_vln_ace21dpi.png', width = 800,height=400)
VlnPlot(dpi1,c('ACE2','TMPRSS2'),ncol = 2, pt.size = 0.2)
dev.off()

dpi2 <- subset(hbec.integrated,idents = c('2dpi'))
Idents(dpi2) <- dpi2[['cell_type2']]
png('hbec.integrated_vln_hla2dpi.png', width = 800,height=400)
VlnPlot(dpi2,c('HLA-C','HLA-A'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_vln_tlr2dpi.png', width = 800,height=400)
VlnPlot(dpi2,c('TLR1','TLR2'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_vln_ace22dpi.png', width = 800,height=400)
VlnPlot(dpi2,c('ACE2','TMPRSS2'),ncol = 2, pt.size = 0.2)
dev.off()

dpi3 <- subset(hbec.integrated,idents = c('3dpi'))
Idents(dpi3) <- dpi3[['cell_type2']]
png('hbec.integrated_vln_hla3dpi.png', width = 800,height=400)
VlnPlot(dpi3,c('HLA-C','HLA-A'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_vln_tlr3dpi.png', width = 800,height=400)
VlnPlot(dpi3,c('TLR1','TLR2'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_vln_ace23dpi.png', width = 800,height=400)
VlnPlot(dpi3,c('ACE2','TMPRSS2'),ncol = 2, pt.size = 0.2)
dev.off()

Idents(hbec.integrated) <- factor(x=Idents(hbec.integrated),levels=c('Mock','1dpi','2dpi','3dpi'))
png('hbec.integrated_fp_hla.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated,'HLA-C',split.by='Condition')
dev.off()

png('hbec.integrated_vln_tlr.png', width = 1200,height=700)
VlnPlot(hbec.integrated,c('TLR4','TLR1','TLR2','TLR7'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_fp_tlr1.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated,'TLR1',split.by='Condition')
dev.off()

png('hbec.integrated_fp_tlr2.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated,'TLR2',split.by='Condition')
dev.off()

png('hbec.integrated_fp_ace2.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated,'ACE2',split.by='Condition')
dev.off()

png('hbec.integrated_fp_tmprss2.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated,'TMPRSS2',split.by='Condition')
dev.off()

ace2 <- FeaturePlot(hbec.integrated,'ACE2') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
tmp <- FeaturePlot(hbec.integrated,'TMPRSS2') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
p1 <- plot_grid(ace2,tmp,ncol=2)

png('hbec.integrated_fp_covid.png', width = 1300, height = 800)
p1
dev.off()

ace2 <- FeaturePlot(ciliated,'ACE2',pt.size=1) + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
tmp <- FeaturePlot(ciliated,'TMPRSS2',pt.size=1) + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
p1 <- plot_grid(ace2,tmp,ncol=2)

png('ciliated_fp_covid.png', width = 1300, height = 800)
p1
dev.off()



# Connectomic featureplots
png('hbec.integrated_fp_plau2.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated,'PLAU',split.by='Condition',label=T)
dev.off()

png('hbec.integrated_fp_plaur2.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated,'PLAUR',split.by='Condition',label=T)
dev.off()

png('hbec.integrated_fp_serpine12.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated,'SERPINE1',split.by='Condition',label=T)
dev.off()

# by timepoint
Idents(hbec.integrated) <- hbec.integrated[['Condition']]
dpi1 <- subset(hbec.integrated,idents = c('1dpi'))
Idents(dpi1) <- dpi1[['cell_type2']]
dpi2 <- subset(hbec.integrated,idents = c('2dpi'))
Idents(dpi2) <- dpi2[['cell_type2']]
dpi3 <- subset(hbec.integrated,idents = c('3dpi'))
Idents(dpi3) <- dpi3[['cell_type2']]
mock <- subset(hbec.integrated,idents = c('Mock'))
Idents(mock) <- mock[['cell_type2']]

png('hbec.integrated_vln_1dpiplau.png', width = 800, height = 400)
VlnPlot(dpi1,c('PLAU','PLAUR'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_vln_2dpiplau.png', width = 800, height = 400)
VlnPlot(dpi2,c('PLAU','PLAUR'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_vln_3dpiplau.png', width = 800, height = 400)
VlnPlot(dpi3,c('PLAU','PLAUR'),ncol = 2, pt.size = 0.2)
dev.off()

png('hbec.integrated_vln_mockplau.png', width = 800, height = 400)
VlnPlot(mock,c('PLAU','PLAUR'),ncol = 2, pt.size = 0.2)
dev.off()






## Ciliated cell subclustering
set.seed(2)
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
ciliated <- subset(hbec.integrated, idents = 'Ciliated')
DefaultAssay(ciliated) <- 'RNA'
ciliated <- NormalizeData(ciliated)
DefaultAssay(ciliated) <- 'integrated'

ciliated <- FindVariableFeatures(ciliated)
ciliated <- ScaleData(ciliated,vars.to.regress = c('pmito','nCount_RNA'))
ciliated <- RunPCA(ciliated, npcs = 30, verbose = FALSE)

ciliated <- RunUMAP(ciliated, reduction = "pca", dims = 1:22)  # 19 is good, 21
ciliated <- FindNeighbors(ciliated,k.param=10,dims = 1:22)  # 19 is good, 21
ciliated <- FindClusters(ciliated,resolution = 0.5)
DefaultAssay(ciliated) <- 'RNA'
png('ciliated_umap_1.png',width = 800,height = 800)
UMAPPlot(ciliated,label = T)
dev.off()

png('ciliated_umap_2.png', width = 800, height = 800)
UMAPPlot(ciliated,label = F, group.by = 'Condition')
dev.off()

png('ciliated_umap_3.png', width = 800, height = 800)
UMAPPlot(ciliated,group.by='Condition',split.by='Condition',ncol=2)
dev.off()

png('ciliated_umap_4.png', width = 800, height = 800)
UMAPPlot(ciliated,label = F, group.by = 'ctype')
dev.off()

png('ciliated_fp_1.png', width = 800, height = 800)
FeaturePlot(ciliated, c('PIFO','AK9','FOXN4','MKI67'), label = T)
dev.off()

png('ciliated_fp_ace2.png', width = 800, height = 800)
FeaturePlot(ciliated, c('ACE2'))
dev.off()

png('ciliated_fp_ace2_2.png', width = 1600, height = 400)
FeaturePlot(ciliated, c('ACE2'),split.by='Condition')
dev.off()


ciliated[['initcluster']] <- Idents(ciliated)
Idents(ciliated) <- ciliated[['initcluster']]
ciliated.marks <- FindAllMarkers(ciliated,only.pos = T)
ciliated.marks$ratio <- ciliated.marks$pct.1/ciliated.marks$pct.2
save(ciliated.marks,file = paste("ciliated.marks",Sys.Date(),".Robj",sep=""))
write.table(ciliated.marks,file = paste("ciliated.marks",Sys.Date(),".txt",sep=""),sep="\t")

# Marker exploration
png('ciliated_vln_1.png', width = 800, height = 400)
VlnPlot(ciliated, c('nFeature_RNA','nCount_RNA'), pt.size = 0)
dev.off()
# 16 is low-info - remove
ciliated <- subset(ciliated, idents = '16', invert = T)
# 11 is doublets - remove
ciliated <- subset(ciliated, idents = '11', invert = T)
# Will have to re-embed once cell types are identified

png('ciliated_fp_dub.png', width = 800, height = 800)
FeaturePlot(subset(ciliated,idents = c('14','15')),c('PIFO','AK9','ASCL1','FOXI1'))
dev.off()

# Remove PNEC doublets from cluster 15
ciliated15 <- subset(ciliated, idents=15)
ciliated15.pnec <- subset(ciliated15, subset = ASCL1 > 1)
png('ciliated_fp_dub15.png', width = 800, height = 800)
UMAPPlot(ciliated15.pnec)
dev.off()
ciliated$sub_cluster <- as.character(Idents(ciliated))
ciliated$sub_cluster[Cells(ciliated15.pnec)] <- paste("pnecs")
Idents(ciliated) <- ciliated[['sub_cluster']]
ciliated <- subset(ciliated,idents='pnecs',invert=T)
ciliated15 <- subset(ciliated, idents=15)
png('ciliated_fp_wodub15.png', width = 800, height = 800)
UMAPPlot(ciliated15)
dev.off()
Idents(ciliated) <- ciliated[['initcluster']]

# Remove Ionocyte doublets from cluster 14
ciliated14 <- subset(ciliated, idents=14)
ciliated14.ion <- subset(ciliated14, subset = FOXI1 > 1)
png('ciliated_fp_dub14.png', width = 800, height = 800)
UMAPPlot(ciliated14.ion)
dev.off()
ciliated$sub_cluster <- as.character(Idents(ciliated))
ciliated$sub_cluster[Cells(ciliated14.ion)] <- paste("ionocytes")
Idents(ciliated) <- ciliated[['sub_cluster']]
ciliated <- subset(ciliated,idents='ionocytes',invert=T)
ciliated14 <- subset(ciliated, idents=14)
png('ciliated_fp_wodub14.png', width = 800, height = 800)
UMAPPlot(ciliated14)
dev.off()
Idents(ciliated) <- ciliated[['initcluster']]

png('ciliated_fp_cont.png', width = 800, height = 800)
FeaturePlot(subset(ciliated,idents = '15'),c('FOXI1','AK9','ASCL1','CHGA'))  # T-Cell/Neutrophil multiplets
dev.off()

## Possible groupings of clusters
png('ciliated_fp_cl0.png', width = 800, height = 800)
FeaturePlot(ciliated,c('ADH6','DCDC2B','ST3GAL6','CATSPERD'))  # 0
dev.off()
png('ciliated_fp_cl4.png', width = 800, height = 800)
FeaturePlot(ciliated,c('IFITM1','HLA-DPA1','IFI27','SAA2'))  # 4
dev.off()
png('ciliated_fp_cl7.png', width = 800, height = 800)
FeaturePlot(ciliated,c('CLDN9','IL32','ICAM1','SERPINB2'))  # 7
dev.off()

png('ciliated_fp_cl1.png', width = 800, height = 800)
FeaturePlot(ciliated,c('CX3CL1','TCN1','PI3','FCGBP'))  # 1
dev.off()
png('ciliated_fp_cl2.png', width = 800, height = 800)
FeaturePlot(ciliated,c('CXCL14','LGALS7B','FLNA','GJA1'))  # 2
dev.off()

png('ciliated_fp_cl3.png', width = 800, height = 800)
FeaturePlot(ciliated,c('SPRR2D','CLCA4','TGM1','THBD'))  # 3
dev.off()
png('ciliated_fp_cl5.png', width = 800, height = 800)
FeaturePlot(ciliated,c('PLAUR','VEGFA','ADM','ANGPTL4'))  # 5
dev.off()

png('ciliated_fp_cl6.png', width = 800, height = 800)
FeaturePlot(ciliated,c('CFAP70','DNAAF1','DNAH12','DDX17'))  # 6
dev.off()

FeaturePlot(ciliated,c('MTRNR2L8','FTH1'))  # 8 - nonspecific

png('ciliated_fp_cl10.png', width = 800, height = 800)
FeaturePlot(ciliated,c('PLAU','scv2-orf1-10','PLA2G4C','BACH2'))  # infected
dev.off()

png('ciliated_fp_prolif.png', width = 800, height = 800)
FeaturePlot(ciliated,c('MKI67'))
dev.off()

png('ciliated_fp_2.png', width = 800, height = 800)
FeaturePlot(ciliated,c('KRT5','KRT14','SCGB1A1','SCGB3A1'))
dev.off()


Idents(ciliated) <- ciliated[['initcluster']]
ciliated <- RenameIdents(ciliated,
                      '0'='ST3GAL6+',
                      '1'='FCGBP+',
                      '2'='CXCL14+',
                      '3'='ADM+',
                      '4'='IFITM1+',
                      '5'='ADM+',
                      '6'='DNAH12+',
                      '7'='CLDN9+',
                      '8'='junk',
                      '9'='FOXN4+',
                      '10'='Infected',
                      '12'='dead',
                      '13'='Proliferating',
                      '14'='ionocytes',  # ? since ionocytes removed
                      '15'='PNECs')  # ? since PNECs removed
ciliated[['ciliated1']] <- Idents(ciliated)

ciliated <- subset(ciliated,idents=c('junk'),invert=T)

png('ciliated_umap_5.png', width = 800, height = 800)
UMAPPlot(ciliated,label = T, group.by = 'ciliated1')
dev.off()

png('ciliated_umap_6.png', width = 800, height = 800)
UMAPPlot(ciliated,group.by='ciliated1',split.by='Condition',ncol=2)
dev.off()

Idents(ciliated) <- ciliated[['initcluster']]
ciliated <- RenameIdents(ciliated,
                      '0'='Immune_Responsive',
                      '1'='Progenitor',
                      '2'='Progenitor',
                      '3'='Progenitor',
                      '4'='Immune_Responsive',
                      '5'='Progenitor',
                      '6'='Mature_Functional',
                      '7'='Immune_Responsive',
                      '9'='Progenitor',
                      '10'='Novel_Infected',
                      '13'='Progenitor')
ciliated[['ciliated2']] <- Idents(ciliated)

Idents(ciliated) <- ciliated[['initcluster']]
ciliated <- RenameIdents(ciliated,
                      '0'='Mature_Ciliated1',
                      '1'='Ciliated_Progenitor',
                      '2'='Ciliated_Progenitor',
                      '3'='Ciliated_Progenitor',
                      '4'='Mature_Ciliated1',
                      '5'='Ciliated_Progenitor',
                      '6'='Mature_Ciliated2',
                      '7'='Mature_Ciliated1',
                      '9'='Ciliated_Progenitor',
                      '10'='Novel_Infected_Ciliated',
                      '13'='Ciliated_Progenitor')
ciliated[['ciliated3']] <- Idents(ciliated)

png('ciliated_umap_7.png', width = 800, height = 800)
UMAPPlot(ciliated,label = T, group.by = 'ciliated2')
dev.off()

png('ciliated_umap_8.png', width = 800, height = 800)
UMAPPlot(ciliated,group.by='ciliated2',split.by='Condition',ncol=2)
dev.off()

# FIGURE PLOTS
Idents(ciliated) <- ciliated[['ciliated3']]
Idents(ciliated) <- factor(x=Idents(ciliated), levels = c('Mature_Ciliated1','Ciliated_Progenitor','Mature_Ciliated2','Novel_Infected_Ciliated'))
cilcolors <- sample(colorRampPalette(RColorBrewer::brewer.pal(n=4, name='Set2'))(4))
cilcolors <- c("#E78AC3", "#66C2A5", "#8DA0CB", "#FC8D62")
png('ciliated_umap_12.png', width = 600, height = 800)
UMAPPlot(ciliated,cols=cilcolors,label = T,label.size=8,repel=TRUE) + labs(title = 'Ciliated Subtype') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoLegend() + NoAxes()
dev.off()

Idents(ciliated) <- ciliated[['Condition']]
Idents(ciliated) <- factor(x=Idents(ciliated), levels = c('1dpi','2dpi','3dpi','Mock'))
png('ciliated_umap_10.png', width = 650, height = 800)
UMAPPlot(ciliated) + labs(title = 'Time Point') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
                                                theme(legend.text = element_text(size = 22),
                                                legend.key.size = unit(1.5, "lines"))
dev.off()

png('ciliated_vln_ace2.png', width = 800,height=400)
VlnPlot(ciliated,c('ACE2','TMPRSS2'),ncol=2,pt.size = 0.2)
dev.off()

png('ciliated_vln_ace2_2.png', width = 800,height=400)
VlnPlot(ciliated,c('ACE2','TMPRSS2'),ncol=2,pt.size = 0.2)
dev.off()

png('ciliated_vln_krt5.png', width = 800,height=400)
VlnPlot(ciliated,c('KRT5','SCGB3A1'),ncol=2,pt.size = 0.2)
dev.off()

png('ciliated_fp_NI.png', width = 800, height = 800)
FeaturePlot(ciliated,c('scv2-orf1-10','IFIT2','NFKBIA','CXCL2'))
dev.off()

png('ciliated_fp_MF.png', width = 800, height = 800)
FeaturePlot(ciliated,c('CFAP70','DNAAF1','DLEC1','DNAH12'))
dev.off()

png('ciliated_fp_IR.png', width = 800, height = 800)
FeaturePlot(ciliated,c('SAA1','HLA-DRB5','PRDX5','PLA2G10'))
dev.off()

png('ciliated_fp_PR.png', width = 800, height = 800)
FeaturePlot(ciliated,c('KRT5','KRT14','SCGB1A1','SCGB3A1'))
dev.off()

ciliated.marks <- FindAllMarkers(ciliated,only.pos = T)
ciliated.marks$ratio <- ciliated.marks$pct.1/ciliated.marks$pct.2
save(ciliated.marks,file = paste("ciliated.marks",Sys.Date(),".Robj",sep=""))
write.table(ciliated.marks,file = paste("ciliated.marks",Sys.Date(),".txt",sep=""),sep="\t")

save(ciliated, file = paste("ciliated",Sys.Date(),".Robj",sep=""))
load('./ciliated2020-10-09.Robj')



### Ciliated sub-population proportions
seuratcolors <- c('#F8766D','#D89000','#A3A500','#39B600','#00BF7D','#00BFC4','#00B0F6','#9590FF','#E76BF3','#FF62BC')

# ciliated1
Idents(ciliated) <- ciliated[['ciliated1']]
ciliated$Condition <- factor(ciliated$Condition,levels=c('Mock','1dpi','2dpi','3dpi'))
ciliated$ciliated1 <- factor(ciliated$ciliated1,levels=c('ST3GAL6+','FCGBP+','CXCL14+','ADM+','IFITM1+','DNAH12+','CLDN9+','FOXN4+','Infected','Proliferating'))
data = table(ciliated$ciliated1,ciliated$Condition)

prop <- t(t(data)/rowSums(t(data))*100)
prop <- round(prop,1)

# Organizing data for the bar graph
sample <- c(rep('Mock',10),rep('1dpi',10),rep('2dpi',10),rep('3dpi',10))
cluster <- rep(c('ST3GAL6+','FCGBP+','CXCL14+','ADM+','IFITM1+','DNAH12+','CLDN9+','FOXN4+','Infected','Proliferating'),4)
value <- as.vector(as.matrix(prop))
dataplot <- data.frame(sample,cluster,value)

# Set the order of your data
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('ST3GAL6+','FCGBP+','CXCL14+','ADM+','IFITM1+','DNAH12+','CLDN9+','FOXN4+','Infected','Proliferating'))

# Run/save the ggplot bar graph
png('ciliated_bar_1.png', width = 1200,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = T) + xlab('') + ylab('% of Cell Type per Time Point') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 20),
        legend.key.size = unit(3, "lines"))
dev.off()

# just ST3GAL6+
prop_ST3GAL6 <- prop['ST3GAL6+',]
value <- as.vector(as.matrix(prop_ST3GAL6))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_ST3GAL6 <- rep(c('ST3GAL6+'),4)
data_ST3GAL6 <- data.frame(sample_short,cluster_ST3GAL6,value)

data_ST3GAL6$sample_short <- factor(data_ST3GAL6$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_ST3GAL6.png', width = 500,height=700)
ggplot(data_ST3GAL6, aes(fill=cluster_ST3GAL6,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[1], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just FCGBP+
prop_FCGBP <- prop['FCGBP+',]
value <- as.vector(as.matrix(prop_FCGBP))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_FCGBP <- rep(c('FCGBP+'),4)
data_FCGBP <- data.frame(sample_short,cluster_FCGBP,value)

data_FCGBP$sample_short <- factor(data_FCGBP$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_FCGBP.png', width = 500,height=700)
ggplot(data_FCGBP, aes(fill=cluster_FCGBP,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[2], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just CXCL14+
prop_CXCL14 <- prop['CXCL14+',]
value <- as.vector(as.matrix(prop_CXCL14))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_CXCL14 <- rep(c('CXCL14+'),4)
data_CXCL14 <- data.frame(sample_short,cluster_CXCL14,value)

data_CXCL14$sample_short <- factor(data_CXCL14$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_CXCL14.png', width = 500,height=700)
ggplot(data_CXCL14, aes(fill=cluster_CXCL14,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[3], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just ADM+
prop_ADM <- prop['ADM+',]
value <- as.vector(as.matrix(prop_ADM))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_ADM <- rep(c('ADM+'),4)
data_ADM <- data.frame(sample_short,cluster_ADM,value)

data_ADM$sample_short <- factor(data_ADM$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_ADM.png', width = 500,height=700)
ggplot(data_ADM, aes(fill=cluster_ADM,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[4], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just IFITM1+
prop_IFITM1 <- prop['IFITM1+',]
value <- as.vector(as.matrix(prop_IFITM1))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_IFITM1 <- rep(c('IFITM1+'),4)
data_IFITM1 <- data.frame(sample_short,cluster_IFITM1,value)

data_IFITM1$sample_short <- factor(data_IFITM1$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_IFITM1.png', width = 500,height=700)
ggplot(data_IFITM1, aes(fill=cluster_IFITM1,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[5], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just DNAH12+
prop_DNAH12 <- prop['DNAH12+',]
value <- as.vector(as.matrix(prop_DNAH12))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_DNAH12 <- rep(c('DNAH12+'),4)
data_DNAH12 <- data.frame(sample_short,cluster_DNAH12,value)

data_DNAH12$sample_short <- factor(data_DNAH12$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_DNAH12.png', width = 500,height=700)
ggplot(data_DNAH12, aes(fill=cluster_DNAH12,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[6], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just CLDN9+
prop_CLDN9 <- prop['CLDN9+',]
value <- as.vector(as.matrix(prop_CLDN9))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_CLDN9 <- rep(c('CLDN9+'),4)
data_CLDN9 <- data.frame(sample_short,cluster_CLDN9,value)

data_CLDN9$sample_short <- factor(data_CLDN9$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_CLDN9.png', width = 500,height=700)
ggplot(data_CLDN9, aes(fill=cluster_CLDN9,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[7], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just FOXN4+
prop_FOXN4 <- prop['FOXN4+',]
value <- as.vector(as.matrix(prop_FOXN4))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_FOXN4 <- rep(c('FOXN4+'),4)
data_FOXN4 <- data.frame(sample_short,cluster_FOXN4,value)

data_FOXN4$sample_short <- factor(data_FOXN4$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_FOXN4.png', width = 500,height=700)
ggplot(data_FOXN4, aes(fill=cluster_FOXN4,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[8], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Infected
prop_Infected <- prop['Infected',]
value <- as.vector(as.matrix(prop_Infected))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Infected <- rep(c('Infected'),4)
data_Infected <- data.frame(sample_short,cluster_Infected,value)

data_Infected$sample_short <- factor(data_Infected$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_Infected.png', width = 500,height=700)
ggplot(data_Infected, aes(fill=cluster_Infected,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[9], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Proliferating
prop_Proliferating <- prop['Proliferating',]
value <- as.vector(as.matrix(prop_Proliferating))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Proliferating <- rep(c('Proliferating'),4)
data_Proliferating <- data.frame(sample_short,cluster_Proliferating,value)

data_Proliferating$sample_short <- factor(data_Proliferating$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_Proliferating.png', width = 500,height=700)
ggplot(data_Proliferating, aes(fill=cluster_Proliferating,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[10], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()


# ciliated2 # FIGURE PLOT
Idents(ciliated) <- ciliated[['ciliated2']]
ciliated$Condition <- factor(ciliated$Condition,levels=c('Mock','1dpi','2dpi','3dpi'))
ciliated$ciliated2 <- factor(ciliated$ciliated2,levels=c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'))
data = table(ciliated$ciliated2,ciliated$Condition)

prop <- t(t(data)/rowSums(t(data))*100)
prop <- round(prop,1)

# Organizing data for the bar graph
sample <- c(rep('Mock',4),rep('1dpi',4),rep('2dpi',4),rep('3dpi',4))
cluster <- rep(c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'),4)
value <- as.vector(as.matrix(prop))
dataplot <- data.frame(sample,cluster,value)

# Set the order of your data
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'))

# Run/save the ggplot bar graph
png('ciliated_bar_3.png', width = 900,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = F) +
  scale_fill_manual("legend", values = c("Immune_Responsive" = cilcolors[1],
                                       "Progenitor" = cilcolors[2],
                                       "Mature_Functional" = cilcolors[3],
                                       'Novel_Infected' = cilcolors[4])) +
  xlab('') + ylab('% of Subtype per Time Point') +
  theme(axis.text.x=element_text(size=35,color = 'black'),
        axis.text.y = element_text(size = 35,color = 'black'),
        axis.title.y = element_text(size = 42,color = 'black'),
        legend.title = element_text(size = 0,color = 'black'),
        legend.text = element_text(size = 20,color = 'black'),
        legend.key.size = unit(3, "lines"))
dev.off()


# just Immune_Responsive
prop_Immune_Responsive <- prop['Immune_Responsive',]
value <- as.vector(as.matrix(prop_Immune_Responsive))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Immune_Responsive <- rep(c('Immune_Responsive'),4)
data_Immune_Responsive <- data.frame(sample_short,cluster_Immune_Responsive,value)

data_Immune_Responsive$sample_short <- factor(data_Immune_Responsive$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_Immune_Responsive.png', width = 500,height=700)
ggplot(data_Immune_Responsive, aes(fill=cluster_Immune_Responsive,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[1], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Progenitor
prop_Progenitor <- prop['Progenitor',]
value <- as.vector(as.matrix(prop_Progenitor))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Progenitor <- rep(c('Progenitor'),4)
data_Progenitor <- data.frame(sample_short,cluster_Progenitor,value)

data_Progenitor$sample_short <- factor(data_Progenitor$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_Progenitor.png', width = 500,height=700)
ggplot(data_Progenitor, aes(fill=cluster_Progenitor,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[2], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Mature_Functional
prop_Mature_Functional <- prop['Mature_Functional',]
value <- as.vector(as.matrix(prop_Mature_Functional))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Mature_Functional <- rep(c('Mature_Functional'),4)
data_Mature_Functional <- data.frame(sample_short,cluster_Mature_Functional,value)

data_Mature_Functional$sample_short <- factor(data_Mature_Functional$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_Mature_Functional.png', width = 500,height=700)
ggplot(data_Mature_Functional, aes(fill=cluster_Mature_Functional,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[3], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Novel_Infected
prop_Novel_Infected <- prop['Novel_Infected',]
value <- as.vector(as.matrix(prop_Novel_Infected))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Novel_Infected <- rep(c('Novel_Infected'),4)
data_Novel_Infected <- data.frame(sample_short,cluster_Novel_Infected,value)

data_Novel_Infected$sample_short <- factor(data_Novel_Infected$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('ciliated_bar_Novel_Infected.png', width = 500,height=700)
ggplot(data_Novel_Infected, aes(fill=cluster_Novel_Infected,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[4], show.legend = F) + xlab('') + ylab('Percent of Total Ciliated') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()





## Basal cell subclustering
set.seed(2)
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
basal <- subset(hbec.integrated, idents = c('Basal','CyclingBasal','BasalEMT'))
DefaultAssay(basal) <- 'RNA'
basal <- NormalizeData(basal)
DefaultAssay(basal) <- 'integrated'

basal <- FindVariableFeatures(basal)
basal <- ScaleData(basal,vars.to.regress = c('pmito','nCount_RNA'))
basal <- RunPCA(basal, npcs = 30, verbose = FALSE)

basal <- RunUMAP(basal, reduction = "pca", dims = 1:20)
basal <- FindNeighbors(basal,k.param=10,dims = 1:20)
basal <- FindClusters(basal,resolution = 0.5)
DefaultAssay(basal) <- 'RNA'
basal[['initcluster']] <- Idents(basal)
png('basal_umap_1.png',width = 800,height = 800)
UMAPPlot(basal,label = T)
dev.off()

png('basal_umap_2.png', width = 800, height = 800)
UMAPPlot(basal,label = F, group.by = 'Condition')
dev.off()

png('basal_umap_3.png', width = 800, height = 800)
UMAPPlot(basal,group.by='Condition',split.by='Condition',ncol=2)
dev.off()

png('basal_umap_4.png', width = 800, height = 800)
UMAPPlot(basal,label = F, group.by = 'cell_type2')
dev.off()

png('basal_umap_5.png', width = 800, height = 800)
UMAPPlot(basal,label = F, group.by = 'ctype')
dev.off()

png('basal_fp_1.png', width = 800, height = 800)
FeaturePlot(basal, c('KRT5','KRT14','MKI67','FN1'), label = T)
dev.off()

png('basal_fp_2.png', width = 800, height = 800)
FeaturePlot(basal, c('PIFO','POU2F3','FOXI1','ASCL1'), label = T)
dev.off()

# 13, 14 are doublets, 10 is low-info - REMOVE
basal <- subset(basal, idents = c('10','13','14'), invert = T)

Idents(basal) <- basal[['initcluster']]
basal.marks <- FindAllMarkers(basal,only.pos = T)
basal.marks$ratio <- basal.marks$pct.1/basal.marks$pct.2
save(basal.marks,file = paste("basal.marks",Sys.Date(),".Robj",sep=""))
write.table(basal.marks,file = paste("basal.marks",Sys.Date(),".txt",sep=""),sep="\t")

# QC
png('basal_vln_1.png', width = 800, height = 400)
VlnPlot(basal, c('nFeature_RNA','nCount_RNA'), pt.size = 0)
dev.off()

# Marker exploration
png('basal_fp_sec.png', width = 800, height = 800)
FeaturePlot(basal, c('SCGB1A1','SCGB3A1','MUC1','LYPD2'))
dev.off()

png('basal_fp_sec9.png', width = 800, height = 800)
FeaturePlot(subset(basal,idents = 'Pre_Secretory'), c('SCGB1A1','SCGB3A1','MUC1','LYPD2'))
dev.off()

png('basal_fp_emt.png', width = 800, height = 800)
FeaturePlot(basal, c('FN1','MMP9','TAGLN','COL1A1'))  # BasalEMT
dev.off()

png('basal_fp_cyc.png', width = 800, height = 800)
FeaturePlot(basal, c('MKI67','TOP2A','PCNA','TUBA1B'))  # CyclingBasal
dev.off()

# cluster 11? ciliated, remove
png('basal_fp_c11.png', width = 800, height = 800)
FeaturePlot(basal, c('DNAH12','DNAAF1','CFAP70','DNAH11'))  # cluster 11
dev.off()
basal <- subset(basal, idents = c('11'), invert = T)

# cluster 8? novel infected?
png('basal_umap_c8.png', width = 800, height = 800)
UMAPPlot(subset(basal,idents = '8'),group.by='Condition',split.by='Condition',ncol=2)
dev.off()

png('basal_fp_c8.png', width = 800, height = 800)
FeaturePlot(basal, c('ISG15','IFI6','IFITM3','CXCL14'))  # cluster 8
dev.off()

# cluster 5? 0? 3? 4? 5+0+3 are not distinct and may represent one end of an archetype, 4 is somewhat distinct
png('basal_fp_c5.png', width = 800, height = 800)
FeaturePlot(basal, c('ELN','DLL1','DLK2','NGFR'))  # cluster 5
dev.off()
png('basal_fp_c0.png', width = 800, height = 800)
FeaturePlot(basal, c('EDN2','NDUFA4L2','CA2','CYP26A1'))  # cluster 0
dev.off()
png('basal_fp_c3.png', width = 800, height = 800)
FeaturePlot(basal, c('DLL1','IGFBP2','WNT10A','MT2A'))  # cluster 3
dev.off()
png('basal_fp_c4.png', width = 800, height = 800)
FeaturePlot(basal, c('PDGFB','MMP2','GATA3','TNC'))  # cluster 4
dev.off()

# cluster 1 + 2?
png('basal_fp_c2.png', width = 800, height = 800)
FeaturePlot(basal, c('CLCA4','RHCG','IVL','CEACAM1'))  # cluster 2
dev.off()
png('basal_fp_plau.png', width = 800, height = 800)
FeaturePlot(basal, c('PLAU','PLAUR','SERPINE1'))  # plau
dev.off()
png('basal_fp_plaur.png', width = 1600, height = 400)
FeaturePlot(basal, c('PLAUR'),split.by='Condition')  # plaur over time
dev.off()
png('basal_fp_plau2.png', width = 1600, height = 400)
FeaturePlot(basal, c('PLAU'),split.by='Condition')  # plaur over time
dev.off()
png('basal_fp_c1.png', width = 800, height = 800)
FeaturePlot(basal, c('DCLK1','LY6D','SERPINB4','NOTCH3'))  # cluster 1
dev.off()

Idents(basal) <- basal[['initcluster']]
basal <- RenameIdents(basal,
                      '0'='EDN2+',
                      '1'='LY6D+',
                      '2'='CEACAM1+',
                      '3'='WNT10A+',
                      '4'='TNC+',
                      '5'='ELN+',
                      '6'='CyclingBasal',
                      '7'='CyclingBasal',
                      '8'='Novel_Infected',
                      '9'='Pre_Secretory',
                      '12'='BasalEMT')
basal[['basal1']] <- Idents(basal)

png('basal_umap_6.png', width = 800, height = 800)
UMAPPlot(basal,label = T, group.by = 'basal1')
dev.off()

png('basal_umap_7.png', width = 800, height = 800)
UMAPPlot(basal,group.by='basal1',split.by='Condition',ncol=2)
dev.off()

png('basal_vln_ace2.png', width = 800,height=400)
VlnPlot(basal,c('ACE2','TMPRSS2'),ncol=2,pt.size = 0.2)
dev.off()


basal.marks <- FindAllMarkers(basal,only.pos = T)
basal.marks$ratio <- basal.marks$pct.1/basal.marks$pct.2
save(basal.marks,file = paste("basal.marks",Sys.Date(),".Robj",sep=""))
write.table(basal.marks,file = paste("basal.marks",Sys.Date(),".txt",sep=""),sep="\t")

save(basal, file = paste("basal",Sys.Date(),".Robj",sep=""))
load('./basal2020-06-26.Robj')


### Basal sub-population proportions
seuratcolors <- c('#F8766D','#D89000','#A3A500','#39B600','#00BF7D','#00BFC4','#00B0F6','#9590FF','#E76BF3','#FF62BC')

# basal1
Idents(basal) <- basal[['basal1']]
basal$Condition <- factor(basal$Condition,levels=c('Mock','1dpi','2dpi','3dpi'))
basal$basal1 <- factor(basal$basal1,levels=c('EDN2+','LY6D+','CEACAM1+','WNT10A+','TNC+','ELN+','CyclingBasal','Novel_Infected','Pre_Secretory','BasalEMT'))
data = table(basal$basal1,basal$Condition)

prop <- t(t(data)/rowSums(t(data))*100)
prop <- round(prop,1)

# Organizing data for the bar graph
sample <- c(rep('Mock',10),rep('1dpi',10),rep('2dpi',10),rep('3dpi',10))
cluster <- rep(c('EDN2+','LY6D+','CEACAM1+','WNT10A+','TNC+','ELN+','CyclingBasal','Novel_Infected','Pre_Secretory','BasalEMT'),4)
value <- as.vector(as.matrix(prop))
dataplot <- data.frame(sample,cluster,value)

# Set the order of your data
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('EDN2+','LY6D+','CEACAM1+','WNT10A+','TNC+','ELN+','CyclingBasal','Novel_Infected','Pre_Secretory','BasalEMT'))

# Run/save the ggplot bar graph
png('basal_bar_1.png', width = 1200,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = T) + xlab('') + ylab('% of Cell Type per Time Point') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 20),
        legend.key.size = unit(3, "lines"))
dev.off()

# just EDN2+
prop_EDN2 <- prop['EDN2+',]
value <- as.vector(as.matrix(prop_EDN2))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_EDN2 <- rep(c('EDN2+'),4)
data_EDN2 <- data.frame(sample_short,cluster_EDN2,value)

data_EDN2$sample_short <- factor(data_EDN2$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_EDN2.png', width = 500,height=700)
ggplot(data_EDN2, aes(fill=cluster_EDN2,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[1], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just LY6D+
prop_LY6D <- prop['LY6D+',]
value <- as.vector(as.matrix(prop_LY6D))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_LY6D <- rep(c('LY6D+'),4)
data_LY6D <- data.frame(sample_short,cluster_LY6D,value)

data_LY6D$sample_short <- factor(data_LY6D$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_LY6D.png', width = 500,height=700)
ggplot(data_LY6D, aes(fill=cluster_LY6D,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[2], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just CEACAM1+
prop_CEACAM1 <- prop['CEACAM1+',]
value <- as.vector(as.matrix(prop_CEACAM1))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_CEACAM1 <- rep(c('CEACAM1+'),4)
data_CEACAM1 <- data.frame(sample_short,cluster_CEACAM1,value)

data_CEACAM1$sample_short <- factor(data_CEACAM1$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_CEACAM1.png', width = 500,height=700)
ggplot(data_CEACAM1, aes(fill=cluster_CEACAM1,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[3], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just WNT10A+
prop_WNT10A <- prop['WNT10A+',]
value <- as.vector(as.matrix(prop_WNT10A))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_WNT10A <- rep(c('WNT10A+'),4)
data_WNT10A <- data.frame(sample_short,cluster_WNT10A,value)

data_WNT10A$sample_short <- factor(data_WNT10A$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_WNT10A.png', width = 500,height=700)
ggplot(data_WNT10A, aes(fill=cluster_WNT10A,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[4], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just TNC+
prop_TNC <- prop['TNC+',]
value <- as.vector(as.matrix(prop_TNC))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_TNC <- rep(c('TNC+'),4)
data_TNC <- data.frame(sample_short,cluster_TNC,value)

data_TNC$sample_short <- factor(data_TNC$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_TNC.png', width = 500,height=700)
ggplot(data_TNC, aes(fill=cluster_TNC,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[5], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just ELN+
prop_ELN <- prop['ELN+',]
value <- as.vector(as.matrix(prop_ELN))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_ELN <- rep(c('ELN+'),4)
data_ELN <- data.frame(sample_short,cluster_ELN,value)

data_ELN$sample_short <- factor(data_ELN$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_ELN.png', width = 500,height=700)
ggplot(data_ELN, aes(fill=cluster_ELN,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[6], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just CyclingBasal
prop_CyclingBasal <- prop['CyclingBasal',]
value <- as.vector(as.matrix(prop_CyclingBasal))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_CyclingBasal <- rep(c('CyclingBasal'),4)
data_CyclingBasal <- data.frame(sample_short,cluster_CyclingBasal,value)

data_CyclingBasal$sample_short <- factor(data_CyclingBasal$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_CyclingBasal.png', width = 500,height=700)
ggplot(data_CyclingBasal, aes(fill=cluster_CyclingBasal,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[7], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Novel_Infected
prop_Novel_Infected <- prop['Novel_Infected',]
value <- as.vector(as.matrix(prop_Novel_Infected))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Novel_Infected <- rep(c('Novel_Infected'),4)
data_Novel_Infected <- data.frame(sample_short,cluster_Novel_Infected,value)

data_Novel_Infected$sample_short <- factor(data_Novel_Infected$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_Novel_Infected.png', width = 500,height=700)
ggplot(data_Novel_Infected, aes(fill=cluster_Novel_Infected,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[8], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Pre_Secretory
prop_Pre_Secretory <- prop['Pre_Secretory',]
value <- as.vector(as.matrix(prop_Pre_Secretory))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Pre_Secretory <- rep(c('Pre_Secretory'),4)
data_Pre_Secretory <- data.frame(sample_short,cluster_Pre_Secretory,value)

data_Pre_Secretory$sample_short <- factor(data_Pre_Secretory$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_Pre_Secretory.png', width = 500,height=700)
ggplot(data_Pre_Secretory, aes(fill=cluster_Pre_Secretory,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[9], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just BasalEMT
prop_BasalEMT <- prop['BasalEMT',]
value <- as.vector(as.matrix(prop_BasalEMT))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_BasalEMT <- rep(c('BasalEMT'),4)
data_BasalEMT <- data.frame(sample_short,cluster_BasalEMT,value)

data_BasalEMT$sample_short <- factor(data_BasalEMT$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('basal_bar_BasalEMT.png', width = 500,height=700)
ggplot(data_BasalEMT, aes(fill=cluster_BasalEMT,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[10], show.legend = F) + xlab('') + ylab('Percent of Total Basal') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# basal2
Idents(basal) <- basal[['basal2']]
basal$Condition <- factor(basal$Condition,levels=c('Mock','1dpi','2dpi','3dpi'))
basal$basal2 <- factor(basal$basal2,levels=c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'))
data = table(basal$basal2,basal$Condition)

prop <- t(t(data)/rowSums(t(data))*100)
prop <- round(prop,1)

# Organizing data for the bar graph
sample <- c(rep('Mock',4),rep('1dpi',4),rep('2dpi',4),rep('3dpi',4))
cluster <- rep(c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'),4)
value <- as.vector(as.matrix(prop))
dataplot <- data.frame(sample,cluster,value)

# Set the order of your data
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'))

# Run/save the ggplot bar graph
png('basal_bar_2.png', width = 1200,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = T) + xlab('') + ylab('% of Cell Type per Time Point') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 20),
        legend.key.size = unit(3, "lines"))
dev.off()





#### Observing Viral Infection ####
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]

png('hbec.integrated_fp_virus1.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('scv2+','scv2.5+','scv2.10+','scv2-orf1-10'), label = T)
dev.off()

png('hbec.integrated2_fp_virus1.png', width = 800, height = 800)
FeaturePlot(hbec.integrated2, c('scv2+','scv2.5+','scv2.10+','scv2-orf1-10'), label = T)
dev.off()

png('hbec.integrated_fp_virus2.png', width = 800, height = 800)
FeaturePlot(hbec.integrated, c('scv2+'), label = T)
dev.off()

png('hbec.integrated_fp_virus23.png', width = 650, height = 800)
FeaturePlot(hbec.integrated,'scv2+',label = F,label.size=8,repel=T) + labs(title = 'SARS-CoV-2 Infection') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()

png('hbec.integrated_fp_ace2.png', width = 650, height = 800)
FeaturePlot(hbec.integrated,'ACE2',label = F,label.size=8,repel=T) + labs(title = 'ACE2') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()


legend(3,30,legend = c('-','+'),col=c("lightgrey","#00ff00"),title = 'Infection',cex=1,plot=FALSE)

png('hbec.integrated_fp_legend.png', width = 900, height = 700)
draw(lgd)
dev.off()

png('hbec.integrated_fp_virus4.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated, c('scv2+'), split.by = 'Condition',label = F)
dev.off()

png('hbec.integrated_fp_virus42.png', width = 1200, height = 400)
FeaturePlot(hbec.integrated, c('scv2+'), split.by = 'Condition',label = F) +
            theme(legend.text = element_text(size = 22),
                legend.key.size = unit(1.5, "lines"))
dev.off()


hbec.integrated$cell_type2 <- factor(hbec.integrated$cell_type2, levels=c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))
png('hbec.integrated_vln_virus1.png', width = 800,height=400)
VlnPlot(hbec.integrated,c('scv2-orf1-10'),pt.size = 0.2,group.by='cell_type2') + xlab('') + theme(legend.position = 'none')
dev.off()


## Infection Histogram # FIGURE PLOT
hbec.integrated$scv2 <- hbec.integrated$'scv2+'
infected <- subset(hbec.integrated, subset = scv2 > 0.5)

Idents(infected) <- infected[['cell_type2']]
infected$Condition <- factor(infected$Condition,levels=c('Mock','1dpi','2dpi','3dpi'))
infected$cell_type2 <- factor(infected$cell_type2,levels=c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))
data = table(infected$cell_type2,infected$Condition)

# plot infected counts
sample <- c(rep('Mock',10),rep('1dpi',10),rep('2dpi',10),rep('3dpi',10))
cluster <- rep(c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'),4)
value <- as.vector(as.matrix(data))
dataplot <- data.frame(sample,cluster,value)
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))

png('hbec.integrated_bar_virus1.png', width = 1200,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = T) + xlab('') + ylab('# Infected Cells per Time Point') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 20),
        legend.key.size = unit(3, "lines"))
dev.off()

# plot percent infected of cell type population at a time point
hbec.integrated$Condition <- factor(hbec.integrated$Condition,levels=c('Mock','1dpi','2dpi','3dpi'))
hbec.integrated$cell_type2 <- factor(hbec.integrated$cell_type2,levels=c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))
data2 = table(hbec.integrated$cell_type2,hbec.integrated$Condition)

prop <- data/data2*100
prop <- round(prop,1)
sample <- c(rep('Mock',10),rep('1dpi',10),rep('2dpi',10),rep('3dpi',10))
cluster <- rep(c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'),4)
value <- as.vector(as.matrix(prop))
dataplot <- data.frame(sample,cluster,value)
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))

png('hbec.integrated_bar_virus2.png', width = 1200,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = T) + xlab('') + ylab('% of Cell Type Infected') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 20),
        legend.key.size = unit(3, "lines"))
dev.off()


# just CyclingBasal
data_CyclingBasal <- prop['CyclingBasal',]
value <- as.vector(as.matrix(data_CyclingBasal))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_CyclingBasal <- rep(c('CyclingBasal'),4)
data_CyclingBasal <- data.frame(sample_short,cluster_CyclingBasal,value)

data_CyclingBasal$sample_short <- factor(data_CyclingBasal$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('hbec.integrated_bar_virus_CyclingBasal.png', width = 500,height=700)
ggplot(data_CyclingBasal, aes(fill=cluster_CyclingBasal,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[1], show.legend = F) + xlab('') + ylab('% Infected CyclingBasal Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Basal
data_Basal <- prop['Basal',]
value <- as.vector(as.matrix(data_Basal))
sample_short <- c('M','1','2','3')
cluster_Basal <- rep(c('Basal'),4)
data_Basal <- data.frame(sample_short,cluster_Basal,value)

data_Basal$sample_short <- factor(data_Basal$sample_short, levels = c('M','1','2','3'))

png('hbec.integrated_bar_virus_Basal.png', width = 300,height=700)
ggplot(data_Basal, aes(fill=cluster_Basal,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[2], show.legend = F) + xlab('') + ylab('% Infected Basal Cells') +
  theme(axis.text.x=element_text(size=35,color='black'),
        axis.text.y = element_text(size = 35,color='black'),
        axis.title.y = element_text(size = 42,color='black'))
dev.off()

# just BasalEMT
data_BasalEMT <- prop['BasalEMT',]
value <- as.vector(as.matrix(data_BasalEMT))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_BasalEMT <- rep(c('BasalEMT'),4)
data_BasalEMT <- data.frame(sample_short,cluster_BasalEMT,value)

data_BasalEMT$sample_short <- factor(data_BasalEMT$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('hbec.integrated_bar_virus_BasalEMT.png', width = 500,height=700)
ggplot(data_BasalEMT, aes(fill=cluster_BasalEMT,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[3], show.legend = F) + xlab('') + ylab('% Infected BasalEMT Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Intermediate
data_Intermediate <- prop['Intermediate',]
value <- as.vector(as.matrix(data_Intermediate))
sample_short <- c('M','1','2','3')
cluster_Intermediate <- rep(c('Intermediate'),4)
data_Intermediate <- data.frame(sample_short,cluster_Intermediate,value)

data_Intermediate$sample_short <- factor(data_Intermediate$sample_short, levels = c('M','1','2','3'))

png('hbec.integrated_bar_virus_Intermediate.png', width = 500,height=700)
ggplot(data_Intermediate, aes(fill=cluster_Intermediate,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[4], show.legend = F) + xlab('') + ylab('% Infected Intermediate Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Club
data_Club <- prop['Club',]
value <- as.vector(as.matrix(data_Club))
sample_short <- c('M','1','2','3')
cluster_Club <- rep(c('Club'),4)
data_Club <- data.frame(sample_short,cluster_Club,value)

data_Club$sample_short <- factor(data_Club$sample_short, levels = c('M','1','2','3'))

png('hbec.integrated_bar_virus_Club.png', width = 500,height=700)
ggplot(data_Club, aes(fill=cluster_Club,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[5], show.legend = F) + xlab('') + ylab('% Infected Club Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Goblet
data_Goblet <- prop['Goblet',]
value <- as.vector(as.matrix(data_Goblet))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Goblet <- rep(c('Goblet'),4)
data_Goblet <- data.frame(sample_short,cluster_Goblet,value)

data_Goblet$sample_short <- factor(data_Goblet$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('hbec.integrated_bar_virus_Goblet.png', width = 500,height=700)
ggplot(data_Goblet, aes(fill=cluster_Goblet,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[6], show.legend = F) + xlab('') + ylab('% Infected Goblet Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Ciliated
data_Ciliated <- prop['Ciliated',]
value <- as.vector(as.matrix(data_Ciliated))
sample_short <- c('M','1','2','3')
cluster_Ciliated <- rep(c('Ciliated'),4)
data_Ciliated <- data.frame(sample_short,cluster_Ciliated,value)

data_Ciliated$sample_short <- factor(data_Ciliated$sample_short, levels = c('M','1','2','3'))

png('hbec.integrated_bar_virus_Ciliated.png', width = 500,height=700)
ggplot(data_Ciliated, aes(fill=cluster_Ciliated,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[7], show.legend = F) + xlab('') + ylab('% Infected Ciliated Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Tuft
data_Tuft <- prop['Tuft',]
value <- as.vector(as.matrix(data_Tuft))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Tuft <- rep(c('Tuft'),4)
data_Tuft <- data.frame(sample_short,cluster_Tuft,value)

data_Tuft$sample_short <- factor(data_Tuft$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('hbec.integrated_bar_virus_Tuft.png', width = 500,height=700)
ggplot(data_Tuft, aes(fill=cluster_Tuft,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[8], show.legend = F) + xlab('') + ylab('% Infected Tuft Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Ionocytes
data_Ionocytes <- prop['Ionocytes',]
value <- as.vector(as.matrix(data_Ionocytes))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_Ionocytes <- rep(c('Ionocytes'),4)
data_Ionocytes <- data.frame(sample_short,cluster_Ionocytes,value)

data_Ionocytes$sample_short <- factor(data_Ionocytes$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('hbec.integrated_bar_virus_Ionocytes.png', width = 500,height=700)
ggplot(data_Ionocytes, aes(fill=cluster_Ionocytes,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[9], show.legend = F) + xlab('') + ylab('% Infected Ionocytes Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just PNECs
data_PNECs <- prop['PNECs',]
value <- as.vector(as.matrix(data_PNECs))
sample_short <- c('Mock','1dpi','2dpi','3dpi')
cluster_PNECs <- rep(c('PNECs'),4)
data_PNECs <- data.frame(sample_short,cluster_PNECs,value)

data_PNECs$sample_short <- factor(data_PNECs$sample_short, levels = c('Mock','1dpi','2dpi','3dpi'))

png('hbec.integrated_bar_virus_PNECs.png', width = 500,height=700)
ggplot(data_PNECs, aes(fill=cluster_PNECs,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[10], show.legend = F) + xlab('') + ylab('% Infected PNECs Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

## Cowplot infected populations
cil <- ggplot(data_Ciliated, aes(fill=cluster_Ciliated,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", color="blue",size=2,fill = celltypecolors[7], show.legend = F) + xlab('') + ylab('% Cells Infected per Population') +
  labs(title = 'Ciliated') +
  theme(axis.text.x=element_text(size=30,color='black'),
        axis.text.y = element_text(size = 30,color='black'),
        axis.title.y = element_text(size = 30,color='black'),plot.title = element_text(size = 30))
clu <- ggplot(data_Club, aes(fill=cluster_Club,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", color="blue",size=2,fill = celltypecolors[5], show.legend = F) + xlab('') + ylab('') +
  labs(title = 'Club') +
  theme(axis.text.x=element_text(size=30,color='black'),
        axis.text.y = element_text(size = 30,color='black'),plot.title = element_text(size = 30))
int <- ggplot(data_Intermediate, aes(fill=cluster_Intermediate,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", color="blue",size=2,fill = celltypecolors[4], show.legend = F) + xlab('') + ylab('') +
  labs(title = 'Intermediate') +
  theme(axis.text.x=element_text(size=30,color='black'),
        axis.text.y = element_text(size = 30,color='black'),plot.title = element_text(size = 30))
bas <- ggplot(data_Basal, aes(fill=cluster_Basal,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", color="blue",size=2,fill = '#D7B32E', show.legend = F) + xlab('') + ylab('') +
  labs(title = 'Basal') +
  theme(axis.text.x=element_text(size=30,color='black'),
        axis.text.y = element_text(size = 30,color='black'),plot.title = element_text(size = 30))
p1 <- plot_grid(cil, clu, int, bas,ncol=4)

png('hbec.integrated_bar_all3.png', width = 1000,height=500)
p1
dev.off()


### Infection Pie Chart # FIGURE PLOT
Idents(infected) <- infected$Condition2
infected$cell_type2 <- factor(infected$cell_type2, levels = cluster)
inf1 <- subset(infected,idents = '1')
inf2 <- subset(infected,idents = '2')
inf3 <- subset(infected,idents = '3')

pop1 = table(inf1$cell_type2)
pop2 = table(inf2$cell_type2)
pop3 = table(inf3$cell_type2)

cluster <- c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs')
value1 <- as.vector(as.matrix(pop1))
value2 <- as.vector(as.matrix(pop2))
value3 <- as.vector(as.matrix(pop3))
pct1 = value1/sum(value1)*100
pct2 = value2/sum(value2)*100
pct3 = value3/sum(value3)*100
data <- data.frame(cluster,pct1,pct2,pct3)

d1 <- ggplot(data, aes(x='',y=value1,fill=cluster)) + geom_bar(width = 1, stat = "identity",fill=celltypecolors,color='blue',size=2) +
  coord_polar("y", start=0) + theme_minimal() + labs(title = '1dpi') +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=0),
        legend.position = 'none',
        plot.title=element_text(size=30,hjust=0.5))
d2 <- ggplot(data, aes(x='',y=value2,fill=cluster)) + geom_bar(width = 1, stat = "identity",fill=celltypecolors,color='blue',size=2) +
  coord_polar("y", start=0) + theme_minimal() + labs(title = '2dpi') +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=0),
        legend.position = 'none',
        plot.title=element_text(size=30,hjust=0.5))
d3 <- ggplot(data, aes(x='',y=value3,fill=cluster)) + geom_bar(width = 1, stat = "identity",fill=celltypecolors,color='blue',size=2) +
  coord_polar("y", start=0) + theme_minimal() + labs(title = '3dpi') +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=0),
        legend.position = 'none',
        plot.title=element_text(size=30,hjust=0.5))
p1 <- plot_grid(d1, d2, d3, ncol=3)

png('infected_pie_2.png', width = 1500,height=500)
p1
dev.off()





## Ciliated Subpopulation Infection # FIGURE PLOT
ciliated$scv2 <- ciliated$'scv2+'
infected <- subset(ciliated, subset = scv2 > 0.5)

Idents(infected) <- infected[['ciliated2']]
infected$Condition <- factor(infected$Condition,levels=c('Mock','1dpi','2dpi','3dpi'))
infected$ciliated2 <- factor(infected$ciliated2,levels=c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'))
data = table(infected$ciliated2,infected$Condition)

ciliated$Condition <- factor(ciliated$Condition,levels=c('Mock','1dpi','2dpi','3dpi'))
ciliated$ciliated2 <- factor(ciliated$ciliated2,levels=c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'))
data2 = table(ciliated$ciliated2,ciliated$Condition)

prop <- data/data2*100
prop <- round(prop,1)
sample <- c(rep('Mock',4),rep('1dpi',4),rep('2dpi',4),rep('3dpi',4))
cluster <- rep(c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'),4)
value <- as.vector(as.matrix(prop))
dataplot <- data.frame(sample,cluster,value)
dataplot$sample <- factor(dataplot$sample,levels=c('Mock','1dpi','2dpi','3dpi'))
dataplot$cluster <- factor(dataplot$cluster,levels=c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected'))

png('ciliated_bar_virus1.png', width = 1200,height=700)
ggplot(dataplot, aes(fill=cluster,x=sample,y=value)) +
  geom_bar(position="dodge", stat="identity", show.legend = T) + xlab('') + ylab('% of Sub-Type Infected') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 20),
        legend.key.size = unit(3, "lines"))
dev.off()

# just Immune_Responsive
data_Immune_Responsive <- prop['Immune_Responsive',]
value <- as.vector(as.matrix(data_Immune_Responsive))
sample_short <- c('M','1','2','3')
cluster_Immune_Responsive <- rep(c('Immune_Responsive'),4)
data_Immune_Responsive <- data.frame(sample_short,cluster_Immune_Responsive,value)

data_Immune_Responsive$sample_short <- factor(data_Immune_Responsive$sample_short, levels = c('M','1','2','3'))

png('ciliated_bar_virus_Immune_Responsive.png', width = 500,height=700)
ggplot(data_Immune_Responsive, aes(fill=cluster_Immune_Responsive,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[1], show.legend = F) + xlab('') + ylab('% Infected Immune_Responsive Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Progenitor
data_Progenitor <- prop['Progenitor',]
value <- as.vector(as.matrix(data_Progenitor))
sample_short <- c('M','1','2','3')
cluster_Progenitor <- rep(c('Progenitor'),4)
data_Progenitor <- data.frame(sample_short,cluster_Progenitor,value)

data_Progenitor$sample_short <- factor(data_Progenitor$sample_short, levels = c('M','1','2','3'))

png('ciliated_bar_virus_Progenitor.png', width = 500,height=700)
ggplot(data_Progenitor, aes(fill=cluster_Progenitor,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[2], show.legend = F) + xlab('') + ylab('% Infected Progenitor Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Mature_Functional
data_Mature_Functional <- prop['Mature_Functional',]
value <- as.vector(as.matrix(data_Mature_Functional))
sample_short <- c('M','1','2','3')
cluster_Mature_Functional <- rep(c('Mature_Functional'),4)
data_Mature_Functional <- data.frame(sample_short,cluster_Mature_Functional,value)

data_Mature_Functional$sample_short <- factor(data_Mature_Functional$sample_short, levels = c('M','1','2','3'))

png('ciliated_bar_virus_Mature_Functional.png', width = 500,height=700)
ggplot(data_Mature_Functional, aes(fill=cluster_Mature_Functional,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[3], show.legend = F) + xlab('') + ylab('% Infected Mature_Functional Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()

# just Novel_Infected
data_Novel_Infected <- prop['Novel_Infected',]
value <- as.vector(as.matrix(data_Novel_Infected))
sample_short <- c('M','1','2','3')
cluster_Novel_Infected <- rep(c('Novel_Infected'),4)
data_Novel_Infected <- data.frame(sample_short,cluster_Novel_Infected,value)

data_Novel_Infected$sample_short <- factor(data_Novel_Infected$sample_short, levels = c('M','1','2','3'))

png('ciliated_bar_virus_Novel_Infected.png', width = 500,height=700)
ggplot(data_Novel_Infected, aes(fill=cluster_Novel_Infected,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", fill = seuratcolors[4], show.legend = F) + xlab('') + ylab('% Infected Novel_Infected Cells') +
  theme(axis.text.x=element_text(size=35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 42))
dev.off()


## Cowplot infected ciliated subpopulations
ir <- ggplot(data_Immune_Responsive, aes(fill=cluster_Immune_Responsive,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", color="blue",size=2,fill = cilcolors[1], show.legend = F) + xlab('') + ylab('% Cells Infected per Subtype') +
  labs(title = 'MC1') +
  theme(axis.text.x=element_text(size=30,color='black'),
        axis.text.y = element_text(size = 30,color='black'),
        axis.title.y = element_text(size = 30,color='black'),plot.title = element_text(size = 30))
pr <- ggplot(data_Progenitor, aes(fill=cluster_Progenitor,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", color="blue",size=2, fill = cilcolors[2], show.legend = F) + xlab('') + ylab('') +
  labs(title = 'CP') +
  theme(axis.text.x=element_text(size=30,color='black'),
        axis.text.y = element_text(size = 30,color='black'),plot.title = element_text(size = 30))
mf <- ggplot(data_Mature_Functional, aes(fill=cluster_Mature_Functional,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", color="blue",size=2,fill = cilcolors[3], show.legend = F) + xlab('') + ylab('') +
  labs(title = 'MC2') +
  theme(axis.text.x=element_text(size=30,color='black'),
        axis.text.y = element_text(size = 30,color='black'),plot.title = element_text(size = 30))
ni <- ggplot(data_Novel_Infected, aes(fill=cluster_Novel_Infected,x=sample_short,y=value)) +
  geom_bar(position="dodge", stat="identity", color="blue",size=2,fill = cilcolors[4], show.legend = F) + xlab('') + ylab('') +
  labs(title = 'NIC') +
  theme(axis.text.x=element_text(size=30,color='black'),
        axis.text.y = element_text(size = 30,color='black'),plot.title = element_text(size = 30))
p1 <- plot_grid(ir, pr, mf, ni, ncol=4)

png('ciliated_bar_all3.png', width = 1000,height=500)
p1
dev.off()




### Ciliated Infection Pie Chart # FIGURE PLOT
Idents(infected) <- infected$Condition2
inf1 <- subset(infected,idents = '1')
inf2 <- subset(infected,idents = '2')
inf3 <- subset(infected,idents = '3')

pop1 = table(inf1$ciliated2)
pop2 = table(inf2$ciliated2)
pop3 = table(inf3$ciliated2)

cluster <- c('Immune_Responsive','Progenitor','Mature_Functional','Novel_Infected')
value1 <- as.vector(as.matrix(pop1))
value2 <- as.vector(as.matrix(pop2))
value3 <- as.vector(as.matrix(pop3))
pct1 = value1/sum(value1)*100
pct2 = value2/sum(value2)*100
pct3 = value3/sum(value3)*100
data <- data.frame(cluster,pct1,pct2,pct3)

d1 <- ggplot(data, aes(x='',y=value1,fill=cluster)) + geom_bar(width = 1, stat = "identity",fill=cilcolors,color='blue',size=2) +
  coord_polar("y", start=0) + theme_minimal() + labs(title = '1dpi') +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=0),
        legend.position = 'none',
        plot.title=element_text(size=50,hjust=0.5))
d2 <- ggplot(data, aes(x='',y=value2,fill=cluster)) + geom_bar(width = 1, stat = "identity",fill=cilcolors,color='blue',size=2) +
  coord_polar("y", start=0) + theme_minimal() + labs(title = '2dpi') +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=0),
        legend.position = 'none',
        plot.title=element_text(size=50,hjust=0.5))
d3 <- ggplot(data, aes(x='',y=value3,fill=cluster)) + geom_bar(width = 1, stat = "identity",fill=cilcolors,color='blue',size=2) +
  coord_polar("y", start=0) + theme_minimal() + labs(title = '3dpi') +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=0),
        legend.position = 'none',
        plot.title=element_text(size=50,hjust=0.5))
p1 <- plot_grid(d1, d2, d3, ncol=3)

png('infectedcil_pie_3.png', width = 1500,height=500)
p1
dev.off()




#### Infection in Ciliated object
Idents(ciliated) <- ciliated[['ciliated1']]

png('ciliated_fp_virus1.png', width = 800, height = 800)
FeaturePlot(ciliated, c('scv2+','scv2.5+','scv2.10+','scv2-orf1-10'), label = T)
dev.off()

png('ciliated_fp_virus2.png', width = 800, height = 800)
FeaturePlot(ciliated, c('scv2+'), label = T)
dev.off()

png('ciliated_fp_virus5.png', width = 800, height = 800)
FeaturePlot(ciliated, c('scv2-orf1-10'), label = T)
dev.off()

png('ciliated_fp_virus4.png', width = 1600, height = 400)
FeaturePlot(ciliated, c('scv2+'), split.by = 'Condition',label = F)
dev.off()

ciliated$cell_type2 <- factor(ciliated$cell_type2, levels=c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))
png('ciliated_vln_virus1.png', width = 800,height=400)
VlnPlot(ciliated,c('scv2-orf1-10'),pt.size = 0.2,group.by='ciliated1') + xlab('') + theme(legend.position = 'none')
dev.off()

png('ciliated_vln_virus2.png', width = 900,height=600)
VlnPlot(ciliated,c('scv2-orf1-10'),pt.size = 0.2,group.by='ciliated2') + xlab('') + theme(legend.position = 'none')
dev.off()

png('ciliated_fp_virus25.png', width = 650, height = 800)
FeaturePlot(ciliated,'scv2+',label = F,label.size=8,repel=T) + labs(title = 'SARS-CoV-2 Infection') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()

png('ciliated_fp_virus24.png', width = 650, height = 800)
FeaturePlot(ciliated,'scv2+',label = T,label.size=8,repel=T) + labs(title = 'SARS-CoV-2 Infection') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()


png('ciliated_fp_virus43.png', width = 1200, height = 400)
FeaturePlot(bigcil, c('scv2+'), split.by = 'Condition',label = F) +
            theme(legend.text = element_text(size = 22),
                legend.key.size = unit(1.5, "lines"))
dev.off()



#### Infection in Basal object
Idents(basal) <- basal[['basal1']]

png('basal_fp_virus1.png', width = 800, height = 800)
FeaturePlot(basal, c('scv2+','scv2.5+','scv2.10+','scv2-orf1-10'), label = T)
dev.off()

png('basal_fp_virus2.png', width = 800, height = 800)
FeaturePlot(basal, c('scv2+'), label = T)
dev.off()

png('basal_fp_virus3.png', width = 1600, height = 400)
FeaturePlot(basal, c('scv2+'), split.by = 'Condition',label = T)
dev.off()

png('basal_vln_virus1.png', width = 800,height=400)
VlnPlot(basal,c('scv2-orf1-10'),pt.size = 0.2,group.by='basal1') + xlab('') + theme(legend.position = 'none')
dev.off()



#### Adding Ciliated & Basal subclustering metadata back to hbec.integrated object ####
hbec.integrated <- AddMetaData(hbec.integrated, metadata = ciliated$ciliated1, col.name = 'ciliated1')
hbec.integrated <- AddMetaData(hbec.integrated, metadata = ciliated$ciliated2, col.name = 'ciliated2')
hbec.integrated <- AddMetaData(hbec.integrated, metadata = basal$basal1, col.name = 'basal1')
hbec.integrated <- AddMetaData(hbec.integrated, metadata = ciliated$ciliated3, col.name = 'ciliated3')

# Unified sub-clustering
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
# create new column in metadata
hbec.integrated$cell_type3 <- as.character(Idents(hbec.integrated))
# add sub-cluster metadata to cell_type3
Idents(ciliated) <- ciliated[['ciliated1']]
hbec.integrated$cell_type3[Cells(ciliated)] <- paste(Idents(ciliated))
png('hbec.integrated_umap_8.png', width = 800, height = 800)
UMAPPlot(hbec.integrated,group.by = 'cell_type3')
dev.off()

# create new column in metadata
hbec.integrated$cell_type4 <- as.character(Idents(hbec.integrated))
# add sub-cluster metadata to cell_type4
Idents(ciliated) <- ciliated[['ciliated2']]
hbec.integrated$cell_type4[Cells(ciliated)] <- paste(Idents(ciliated))
png('hbec.integrated_umap_9.png', width = 800, height = 800)
UMAPPlot(hbec.integrated,group.by = 'cell_type4')
dev.off()

Idents(hbec.integrated) <- hbec.integrated[['cell_type3']]
hbec.integrated <- RenameIdents(hbec.integrated,'Ciliated'='Remove')
hbec.integrated[['cell_type3']] <- Idents(hbec.integrated)

Idents(hbec.integrated) <- hbec.integrated[['cell_type4']]
hbec.integrated <- RenameIdents(hbec.integrated,'Ciliated'='Remove')
hbec.integrated[['cell_type4']] <- Idents(hbec.integrated)





#### Exploring Ciliated Cells - 7/17/20 ####
Idents(ciliated) <- ciliated[['scv2+']]
ciliated_infected <- subset(ciliated,idents = '1')
ciliated_uninfected <- subset(ciliated,idents = '0')

png('ciliated_umap_infected.png', width = 800, height = 800)
UMAPPlot(ciliated_infected,label = T,group.by='ciliated2')
dev.off()

png('ciliated_umap_uninfected.png', width = 800, height = 800)
UMAPPlot(ciliated_uninfected,label = T,group.by='ciliated2')
dev.off()


## Calculating genes that differentiate infected and uninfected ciliated cells
ciliated.marks <- FindAllMarkers(ciliated,only.pos = T)
ciliated.marks$ratio <- ciliated.marks$pct.1/ciliated.marks$pct.2
save(ciliated.marks,file = paste("ciliated.marks",Sys.Date(),".Robj",sep=""))
write.table(ciliated.marks,file = paste("ciliated.marks",Sys.Date(),".txt",sep=""),sep="\t")

## DotPlot of infection genes
infected.marks <- ciliated.marks[ which(ciliated.marks$cluster=='1'), ]
infected.marks <- infected.marks[order(-infected.marks$avg_logFC),]
genes.use.i <- head(infected.marks,20)
genes.use.i <- genes.use.i$gene

uninfected.marks <- ciliated.marks[ which(ciliated.marks$cluster=='0'), ]
uninfected.marks <- uninfected.marks[order(-uninfected.marks$avg_logFC),]
genes.use.u <- head(uninfected.marks,20)
genes.use.u <- genes.use.u$gene

genes.use <- rev(c(genes.use.i,genes.use.u))

DefaultAssay(ciliated) <- 'RNA'
ciliated <- RenameIdents(ciliated,'0'='Uninfected','1'='Infected')
ciliated[['infect']] <- Idents(ciliated)

ciliated <- ScaleData(ciliated,vars.to.regress = c('pmito','nCount_RNA'),features = genes.use)

# DotPlot
Idents(ciliated) <- ciliated[['infect']]
Idents(ciliated) <- factor(x=Idents(ciliated), levels = rev(c('Infected','Uninfected')))

group_colors <- colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Blues'))(2)

dot <- DotPlot(ciliated,features = genes.use,assay='RNA',
            cols = group_colors,scale = T,dot.scale = 14) +
            theme(plot.title = element_text(size=50,color = 'black',hjust = 0.5),
                    axis.text.x=element_text(size=32, color = 'black',angle=90),
                    axis.text.y = element_text(size = 42, color = 'black'),
                    legend.text = element_text(size = 30),
                    legend.key.size = unit(20, "mm"),
                    legend.title = element_text(size=30),
                    strip.text.x = element_text(size=44),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(color = "black"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank())

png('ciliated_dot_1.png', width = 1500, height = 1000)
dot
dev.off()






#### Overall cell identity heat map #### FIGURE PLOT
b.genes <- c('MKI67','TOP2A','KRT5','KRT14','IGFBP6')  # CyclingBasal + Basal
emt.genes <- c('MMP9','FN1','VIM')  # BasalEMT
cl.genes <- c('SCGB3A1','SCGB1A1','LYPD2','BPIFB1')  # Club
g.genes <- c('AGR2','MUC5AC')  # Goblet
c.genes <- c('FOXJ1','DNAH6','PIFO')  # Ciliated
t.genes <- c('POU2F3','ASCL2')  # Tuft
i.genes <- c('FOXI1','CFTR','ASCL3','POSTN')  # Ionocyte
p.genes <- c('CHGA','ASCL1','CALCA')  # PNEC

genes.use <- c(b.genes,emt.genes,cl.genes,g.genes,c.genes,t.genes,i.genes,p.genes)

DefaultAssay(hbec.integrated) <- 'RNA'
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec.integrated$cell_type2 <- factor(x=hbec.integrated$cell_type2,
                             levels = c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet',
                             'Ciliated','Tuft','Ionocytes','PNECs'))
hbec.integrated <- hbec.integrated[order(hbec.integrated$cell_type2),]
identities <- c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs')

# prepare data
hbec.integrated <- ScaleData(hbec.integrated,vars.to.regress = c('pmito','nCount_RNA'),features = genes.use)
hbec.integrated.avg <- AverageExpression(hbec.integrated,use.scale = T,assays = 'RNA')
hbec.integrated.avg <- t(scale(t(hbec.integrated.avg$RNA)))
unityNormalize <- function(x){
  (x-min(x))/(max(x)-min(x))
}
hbec.integrated.avg <- t(apply(as.matrix(hbec.integrated.avg), 1, unityNormalize))
hbec.integrated.avg <- hbec.integrated.avg[,c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet',
                        'Ciliated','Tuft','Ionocytes','PNECs')]
dim(hbec.integrated.avg)
#[1] 26 10

# Prep & plot
colors.inferno <- colorRamp2(breaks = c(seq(min(hbec.integrated.avg),max(hbec.integrated.avg),length.out=60)), inferno(n=60), space = "RGB")

hbec.integrated.heatmap <- Heatmap(as.matrix(hbec.integrated.avg),
                        col= colors.inferno,
                        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=celltypecolors))),
                        column_split= factor(colnames(hbec.integrated.avg), levels = identities),
                        column_title = identities,
                        column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 18),
                        column_title_rot = 90,
                        cluster_rows=F,
                        cluster_columns=F,
                        cluster_column_slices=F,
                        show_column_names=FALSE,
                        show_row_names = T,
                        row_order = genes.use,
                        show_heatmap_legend = F)

png('hbec.integrated_heatmap_3.png', width = 360, height = 700)
draw(hbec.integrated.heatmap,padding = unit(c(2, 2, 4, 2), "mm"))
dev.off()




#### Ciliated Subtype identity heat map #### FIGURE PLOT
load('./ciliated2020-10-09.Robj')

ir.genes <- c('SAA1','HLA-DRB5','PRDX5','PLA2G10')  # Immune Responsive
p.genes <- c('KRT5','KRT14','SCGB1A1','SCGB3A1')  # Progenitor
mf.genes <- c('CFAP70','DNAAF1','DLEC1','DNAH12')  # Mature Functional
ni.genes <- c('IFIT2','NFKBIA','CXCL2')  # Novel Infected - how to include scv2-orf1-10 ?
entry.genes <- c('HMGB1','NRP1','CTSL','ACE2','TMPRSS2','BSG','FURIN',"ITGB1","SERPINB3","CD81","CAV2","CAV1","RPSA","CLDN1","PPIA","EPHA2","NECTIN1","CTSB")
ifn.genes2 <- unique(c("CYLD","IRF1","AF117829","PLCG2","REL","TNFAIP3","ZC3HAV1","TICAM1","NFKB1","NFKB2","RELB"))

length(rownames(GetAssayData(hbec.integrated, slot = "data")))
length(rownames(GetAssayData(ciliated, slot = "data")))

genes.use <- c(ir.genes,p.genes,mf.genes,ni.genes)
#genes.use <- intersect(rownames(ciliated),genes.use)
genes.use <- entry.genes

DefaultAssay(ciliated) <- 'RNA'
Idents(ciliated) <- ciliated[['ciliated3']]
ciliated$ciliated3 <- factor(x=ciliated$ciliated3,
                             levels = c('Mature_Ciliated1','Ciliated_Progenitor','Mature_Ciliated2','Novel_Infected_Ciliated'))
#ciliated <- ciliated[order(ciliated$ciliated3),]
identities <- c('Mature_Ciliated1','Ciliated_Progenitor','Mature_Ciliated2','Novel_Infected_Ciliated')
identities2 <- c('MC1','CP','MC2','NIC')

# prepare data
ciliated <- ScaleData(ciliated,features = rownames(ciliated))
#genes.int <- intersect(rownames(GetAssayData(ciliated, slot = "scale.data")),genes.use)
ciliated.avg <- AverageExpression(ciliated,assays = 'RNA',slot='scale.data',features=genes.use)
ciliated.avg <- t(scale(t(ciliated.avg$RNA)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
ciliated.avg <- t(apply(as.matrix(ciliated.avg), 1, unityNormalize))
ciliated.avg <- ciliated.avg[,c('Mature_Ciliated1','Ciliated_Progenitor','Mature_Ciliated2','Novel_Infected_Ciliated')]
dim(ciliated.avg)
#[1] 18  4

# Prep & plot
colors.inferno <- colorRamp2(breaks = c(seq(min(ciliated.avg),max(ciliated.avg),length.out=60)), inferno(n=60), space = "RGB")

ciliated.heatmap <- Heatmap(as.matrix(ciliated.avg),
                        col= colors.inferno,
                        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=cilcolors))),
                        column_split= factor(colnames(ciliated.avg), levels = identities),
                        column_title = identities2,
                        column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 18),
                        column_title_rot = 90,
                        cluster_rows=T,
                        row_km = 1,
                        cluster_columns=F,
                        cluster_column_slices=F,
                        show_column_names=FALSE,
                        show_row_names = T,
                        row_order = genes.use,
                        show_heatmap_legend = F)

png('ciliated_heatmap_entry2.png', width = 250, height = 500)
draw(ciliated.heatmap,padding = unit(c(2, 2, 4, 2), "mm"))
dev.off()


myplots <- vector('list', length(entry.genes))
for (i in 1:length(entry.genes)) {
    message(i)
    gene <- entry.genes[i]
    myplots[[i]] <- local({
        i <- i
        gene <- gene
        p1 <- FeaturePlot(ciliated,gene,label = F, pt.size=1) + labs(title = paste(gene)) + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
                    theme(legend.text = element_text(size = 22),
                    legend.key.size = unit(1.5, "lines"))
    })
    png(paste('ciliated_fp_',gene,'.png',sep=""), width = 650, height = 800)
    plot(myplots[[i]])
    dev.off()
}





### Making GO bar plot # FIGURE PLOT
#allgo <- data.frame(read.delim('./hbec.integratedGO.txt'))
#scp /Volumes/10x/hbec.integratedGO3.txt amg99@transfer-farnam.hpc.yale.edu:/home/amg99/project/covid19
#scp /Volumes/10x/ciliatedGO2.txt amg99@transfer-farnam.hpc.yale.edu:/home/amg99/project/covid19
allgo <- data.frame(read.delim('./hbec.integratedGO3.txt'))
allgo$Sample <- factor(allgo$Sample, levels = rev(c('Mock','1dpi','2dpi','3dpi')))
allgo$GO_short <- factor(allgo$GO_short, levels = rev(unique(allgo$GO_short)))

seuratcolors <- hue_pal()(4)
seuratcolorsrep <- rev(c(rep(seuratcolors[1],5),rep(seuratcolors[2],5),rep(seuratcolors[3],5),rep(seuratcolors[4],5)))

#ggsave('hbec.integrated_gobar_2.png',p1,width = 550,height=300,units = c("mm"))
png('hbec.integrated_gobar_4.png', width = 1600,height=1000)
ggplot(allgo,aes(x=neglog10pval,y=GO_short,fill=GO_short)) +
    geom_col(position = position_dodge2(preserve = "single"),show.legend = F) +
    scale_fill_manual(values = seuratcolorsrep) +
    xlab('-log(p_val)') + ylab('') +
    theme(axis.text.x=element_text(size=35,color='black'),
          axis.text.y = element_text(size = 35,color='black'),
          axis.title.x = element_text(size = 40,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
dev.off()

#plot.margin = margin(1, 20, 1, 1, "cm")
png('hbec.integrated_gobar_2.png', width = 1800,height=1500)
ggplot(allgo,aes(x=neglog10pval,y=Sample,fill=GO_short,label=GO_short)) +
    geom_col(position = position_dodge2(preserve = "single"),show.legend = F) +
    scale_fill_manual(values = seuratcolorsrep) +
    xlab('-log(p_val)') + ylab('') + labs(title = 'GO Terms over Time Point') +
    geom_text(position = position_dodge2(0.9,preserve='single'),hjust=-0.02,size=12) +
    scale_x_continuous(limits=c(0, 50)) +
    theme(plot.title = element_text(size = 45, hjust = 0.5,face = 'bold'),
          axis.text.x=element_text(size=35,color='black'),
          axis.text.y = element_text(size = 42,color='black'),
          axis.title.y = element_text(size = 35,color='black'),
          axis.title.x = element_text(size = 35,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
dev.off()

# Making custom legend
png('hbec.integrated_legend_1.png', width = 600,height=600)
ggplot(allgo,aes(x=neglog10pval,y=Sample,fill=Sample)) +
    geom_col(position = position_dodge2(preserve = "single"),show.legend = T) +
    scale_fill_manual(values = rev(seuratcolors),guide = guide_legend(reverse = TRUE)) +
    labs(fill = "Time Point") +
    theme(legend.title = element_text(size = 50,color = 'black',face = 'bold'),
        legend.text = element_text(size = 40,color = 'black'),
        legend.key.size = unit(6, "lines"))
dev.off()




# Ciliated Subtype GO - FIGURE PLOT
cilgo <- data.frame(read.delim('./ciliatedGO2.txt'))
cilgo$Sample <- factor(cilgo$Sample, levels = rev(c('MC1','CP','MC2','NIC')))
cilgo$GO_short <- factor(cilgo$GO_short, levels = rev(unique(cilgo$GO_short)))

cilcolors <- c("#E78AC3", "#66C2A5", "#8DA0CB", "#FC8D62")
cilcolorsrep <- rev(c(rep(cilcolors[2],5),rep(cilcolors[1],5),rep(cilcolors[3],5),rep(cilcolors[4],5)))

png('ciliated_gobar_2.png', width = 1600,height=1000)
ggplot(cilgo,aes(x=neglogpval,y=GO_short,fill=GO_short)) +
    geom_col(position = position_dodge2(preserve = "single"),show.legend = F) +
    scale_fill_manual(values = cilcolorsrep) +
    xlab('-log(p_val)') + ylab('') +
    theme(axis.text.x=element_text(size=35,color='black'),
          axis.text.y = element_text(size = 35,color='black'),
          axis.title.x = element_text(size = 40,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
dev.off()

# Making custom legend
png('ciliated_legend_2.png', width = 600,height=600)
ggplot(cilgo,aes(x=neglogpval,y=Sample,fill=Sample)) +
    geom_col(position = position_dodge2(preserve = "single"),show.legend = T) +
    scale_fill_manual(values = rev(cilcolors),guide = guide_legend(reverse = TRUE)) +
    labs(fill = "Subtype") +
    theme(legend.title = element_text(size = 50,color = 'black',face = 'bold'),
        legend.text = element_text(size = 40,color = 'black'),
        legend.key.size = unit(6, "lines"))
dev.off()







#### Gene variance over time
### Ciliated cells
load('./ciliated2020-10-09.Robj')
load('./hbec.ciliated.marks22020-06-02.Robj')
genes.use <- hbec.ciliated.marks %>% group_by(cluster) %>% top_n(40, avg_logFC)
## Sam's marker list's
# Get genes which have met our statistical thresholds (v1)
# V2
toMatch <- c('VEGF','PDGF','HBEGF','CSF','BMP','ANGPT','IL6','EGF','FGF','HGF','TNF','TGF',
            'CSF','CCL2','IL8','IL6','RETN','VWF','SERPINE1','THBD','SELP','ANGPT2','C1Q','PLAU','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(ciliated), value=TRUE))
matches.exclude <- matches[grep('R',matches)]
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# V3
toMatch <- c('VEGF','PDGF','HBEGF','CYR','BMP','ANGPT','EGF','FGF','HGF','TGF','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(ciliated), value=TRUE))
matches.exclude <- NULL
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# Yap downstream
genes.use <- unique(c('YAP1','AREG','FGF1','CTGF',
                    'SMAD7','AXIN2','SERPINE1','ITGB2','GLI1','BBC3',
                    'AFP','ID1','ID2','NKD1','MYC','CCND1',
                    'SOX2','SNAI2','BIRC2','BIRC5'))

DefaultAssay(ciliated) <- 'RNA'
Idents(ciliated) <- ciliated[['Condition']]
ciliated$Condition <- factor(x=ciliated$Condition,
                             levels = c('Mock','1dpi','2dpi','3dpi'))
ciliated <- ciliated[order(ciliated$Condition),]
identities <- c('Mock','1dpi','2dpi','3dpi')

# prepare data
ciliated <- ScaleData(ciliated,features = genes.use)
ciliated.avg <- AverageExpression(ciliated,assays = 'RNA',slot='scale.data')  # ,features=genes.use
ciliated.avg <- t(scale(t(ciliated.avg$RNA)))
unityNormalize <- function(x){
  (x-min(x))/(max(x)-min(x))
}
ciliated.avg <- t(apply(as.matrix(ciliated.avg), 1, unityNormalize))
ciliated.avg <- ciliated.avg[,c('Mock','1dpi','2dpi','3dpi')]
dim(ciliated.avg)
#[1] 117  4

# Prep & plot
colors.inferno <- colorRamp2(breaks = c(seq(min(ciliated.avg),max(ciliated.avg),length.out=60)), inferno(n=60), space = "RGB")
seuratcolors <- hue_pal()(4)

ciliated.heatmap <- Heatmap(as.matrix(ciliated.avg),
                        col= colors.inferno,
                        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=seuratcolors))),
                        column_split= factor(colnames(ciliated.avg), levels = identities),
                        column_title = identities,
                        column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 18),
                        column_title_rot = 90,
                        cluster_rows=T,
                        row_km = 4,
                        cluster_columns=F,
                        cluster_column_slices=F,
                        show_column_names=FALSE,
                        show_row_names = T,
                        show_heatmap_legend = F)

png('ciliated_timehm_3.png', width = 250, height = 1700)
draw(ciliated.heatmap,padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


### Basal cells
load('./basal2020-06-26.Robj')
load('./hbec.basal.marks2020-06-02.Robj')
genes.use <- hbec.basal.marks %>% group_by(cluster) %>% top_n(40, avg_logFC)
## Sam's marker list's
# Get genes which have met our statistical thresholds (v1)
# V2
toMatch <- c('VEGF','PDGF','HBEGF','CSF','BMP','ANGPT','IL6','EGF','FGF','HGF','TNF','TGF',
            'CSF','CCL2','IL8','IL6','RETN','VWF','SERPINE1','THBD','SELP','ANGPT2','C1Q','PLAU','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(basal), value=TRUE))
matches.exclude <- matches[grep('R',matches)]
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# V3
toMatch <- c('VEGF','PDGF','HBEGF','CYR','BMP','ANGPT','EGF','FGF','HGF','TGF','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(basal), value=TRUE))
matches.exclude <- NULL
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# Yap downstream
genes.use <- unique(c('YAP1','AREG','FGF1','CTGF',
                    'SMAD7','AXIN2','SERPINE1','ITGB2','GLI1','BBC3',
                    'AFP','ID1','ID2','NKD1','MYC','CCND1',
                    'SOX2','SNAI2','BIRC2','BIRC5'))

DefaultAssay(basal) <- 'RNA'
Idents(basal) <- basal[['Condition']]
basal$Condition <- factor(x=basal$Condition,
                             levels = c('Mock','1dpi','2dpi','3dpi'))
basal <- basal[order(basal$Condition),]
identities <- c('Mock','1dpi','2dpi','3dpi')

# prepare data
basal <- ScaleData(basal,features = genes.use)  # $gene
basal.avg <- AverageExpression(basal,assays = 'RNA',slot='scale.data')  # ,features=genes.use
basal.avg <- t(scale(t(basal.avg$RNA)))
unityNormalize <- function(x){
  (x-min(x))/(max(x)-min(x))
}
basal.avg <- t(apply(as.matrix(basal.avg), 1, unityNormalize))
basal.avg <- basal.avg[,c('Mock','1dpi','2dpi','3dpi')]
dim(basal.avg)
#[1] 92  4

#  basal.avg <- basal.avg[-c(3), ]
# Prep & plot
colors.inferno <- colorRamp2(breaks = c(seq(min(basal.avg),max(basal.avg),length.out=60)), inferno(n=60), space = "RGB")
seuratcolors <- hue_pal()(4)

basal.heatmap <- Heatmap(as.matrix(basal.avg),
                        col= colors.inferno,
                        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=seuratcolors))),
                        column_split= factor(colnames(basal.avg), levels = identities),
                        column_title = identities,
                        column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 18),
                        column_title_rot = 90,
                        cluster_rows=T,
                        row_km = 4,
                        cluster_columns=F,
                        cluster_column_slices=F,
                        show_column_names=FALSE,
                        show_row_names = T,
                        show_heatmap_legend = F)

png('basal_timehm_3.png', width = 250, height = 1500)
draw(basal.heatmap,padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


### All cells
load('./hbec.integrated2020-10-09.Robj')
load('./hbec.int.Condwmock.marks2020-06-04.Robj')
genes.use <- hbec.basal.marks %>% group_by(cluster) %>% top_n(40, avg_logFC)
## Sam's marker list's
# Get genes which have met our statistical thresholds (v1)
# V2
toMatch <- c('VEGF','PDGF','HBEGF','CSF','BMP','ANGPT','IL6','EGF','FGF','HGF','TNF','TGF',
            'CSF','CCL2','IL8','IL6','RETN','VWF','SERPINE1','THBD','SELP','ANGPT2','C1Q','PLAU','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(hbec.integrated), value=TRUE))
matches.exclude <- matches[grep('R',matches)]
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# V3
toMatch <- c('VEGF','PDGF','HBEGF','CYR','BMP','ANGPT','EGF','FGF','HGF','TGF','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(hbec.integrated), value=TRUE))
matches.exclude <- NULL
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# Yap downstream
genes.use <- unique(c('YAP1','AREG','FGF1','CTGF',
                    'SMAD7','AXIN2','SERPINE1','ITGB2','GLI1','BBC3',
                    'AFP','ID1','ID2','NKD1','MYC','CCND1',
                    'SOX2','SNAI2','BIRC2','BIRC5'))

DefaultAssay(hbec.integrated) <- 'RNA'
Idents(hbec.integrated) <- hbec.integrated[['Condition']]
hbec.integrated$Condition <- factor(x=hbec.integrated$Condition,
                             levels = c('Mock','1dpi','2dpi','3dpi'))
hbec.integrated <- hbec.integrated[order(hbec.integrated$Condition),]
identities <- c('Mock','1dpi','2dpi','3dpi')

# prepare data
hbec.integrated <- ScaleData(hbec.integrated,features = genes.use)  # $gene
hbec.integrated.avg <- AverageExpression(hbec.integrated,assays = 'RNA',slot='scale.data')  # ,features=genes.use
hbec.integrated.avg <- t(scale(t(hbec.integrated.avg$RNA)))
unityNormalize <- function(x){
  (x-min(x))/(max(x)-min(x))
}
hbec.integrated.avg <- t(apply(as.matrix(hbec.integrated.avg), 1, unityNormalize))
hbec.integrated.avg <- hbec.integrated.avg[,c('Mock','1dpi','2dpi','3dpi')]
dim(hbec.integrated.avg)
#[1] 134   4

#  hbec.integrated.avg <- hbec.integrated.avg[-c(3), ]
# Prep & plot
colors.inferno <- colorRamp2(breaks = c(seq(min(hbec.integrated.avg),max(hbec.integrated.avg),length.out=60)), inferno(n=60), space = "RGB")
seuratcolors <- hue_pal()(4)

hbec.integrated.heatmap <- Heatmap(as.matrix(hbec.integrated.avg),
                        col= colors.inferno,
                        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=seuratcolors))),
                        column_split= factor(colnames(hbec.integrated.avg), levels = identities),
                        column_title = identities,
                        column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 18),
                        column_title_rot = 90,
                        cluster_rows=T,
                        row_km = 4,
                        cluster_columns=F,
                        cluster_column_slices=F,
                        show_column_names=FALSE,
                        show_row_names = T,
                        show_heatmap_legend = F)

png('hbec.integrated_timehm_3.png', width = 250, height = 1800)
draw(hbec.integrated.heatmap,padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


### Intermediate cells
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
intermediate <- subset(hbec.integrated,idents = 'Intermediate')
load('./hbec.inter.marks2020-07-31.Robj')
genes.use <- hbec.inter.marks %>% group_by(cluster) %>% top_n(40, avg_logFC)
## Sam's marker list's
# Get genes which have met our statistical thresholds (v1)
# V2
toMatch <- c('VEGF','PDGF','HBEGF','CSF','BMP','ANGPT','IL6','EGF','FGF','HGF','TNF','TGF',
            'CSF','CCL2','IL8','IL6','RETN','VWF','SERPINE1','THBD','SELP','ANGPT2','C1Q','PLAU','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(intermediate), value=TRUE))
matches.exclude <- matches[grep('R',matches)]
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# V3
toMatch <- c('VEGF','PDGF','HBEGF','CYR','BMP','ANGPT','EGF','FGF','HGF','TGF','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(intermediate), value=TRUE))
matches.exclude <- NULL
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# Yap downstream
genes.use <- unique(c('YAP1','AREG','FGF1','CTGF',
                    'SMAD7','AXIN2','SERPINE1','ITGB2','GLI1','BBC3',
                    'AFP','ID1','ID2','NKD1','MYC','CCND1',
                    'SOX2','SNAI2','BIRC2','BIRC5'))

DefaultAssay(intermediate) <- 'RNA'
Idents(intermediate) <- intermediate[['Condition']]
intermediate$Condition <- factor(x=intermediate$Condition,
                             levels = c('Mock','1dpi','2dpi','3dpi'))
intermediate <- intermediate[order(intermediate$Condition),]
identities <- c('Mock','1dpi','2dpi','3dpi')

# prepare data
intermediate <- ScaleData(intermediate,features = genes.use)  # $gene
intermediate.avg <- AverageExpression(intermediate,assays = 'RNA',slot='scale.data')  # ,features=genes.use
intermediate.avg <- t(scale(t(intermediate.avg$RNA)))
unityNormalize <- function(x){
  (x-min(x))/(max(x)-min(x))
}
intermediate.avg <- t(apply(as.matrix(intermediate.avg), 1, unityNormalize))
intermediate.avg <- intermediate.avg[,c('Mock','1dpi','2dpi','3dpi')]
dim(intermediate.avg)
#[1] 77  4

#  intermediate.avg <- intermediate.avg[-c(3), ]
# Prep & plot
colors.inferno <- colorRamp2(breaks = c(seq(min(intermediate.avg),max(intermediate.avg),length.out=60)), inferno(n=60), space = "RGB")
seuratcolors <- hue_pal()(4)

intermediate.heatmap <- Heatmap(as.matrix(intermediate.avg),
                        col= colors.inferno,
                        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=seuratcolors))),
                        column_split= factor(colnames(intermediate.avg), levels = identities),
                        column_title = identities,
                        column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 18),
                        column_title_rot = 90,
                        cluster_rows=T,
                        row_km = 4,
                        cluster_columns=F,
                        cluster_column_slices=F,
                        show_column_names=FALSE,
                        show_row_names = T,
                        show_heatmap_legend = F)

png('intermediate_timehm_3.png', width = 250, height = 1500)
draw(intermediate.heatmap,padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


### Club cells
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
club <- subset(hbec.integrated,idents = 'Club')
load('./hbec.club.marks2020-06-02.Robj')
genes.use <- hbec.inter.marks %>% group_by(cluster) %>% top_n(40, avg_logFC)
## Sam's marker list's
# Get genes which have met our statistical thresholds (v1)
# V2
toMatch <- c('VEGF','PDGF','HBEGF','CSF','BMP','ANGPT','IL6','EGF','FGF','HGF','TNF','TGF',
            'CSF','CCL2','IL8','IL6','RETN','VWF','SERPINE1','THBD','SELP','ANGPT2','C1Q','PLAU','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(club), value=TRUE))
matches.exclude <- matches[grep('R',matches)]
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# V3
toMatch <- c('VEGF','PDGF','HBEGF','CYR','BMP','ANGPT','EGF','FGF','HGF','TGF','DLL','REG')
matches <- unique(grep(paste(toMatch,collapse="|"),
                       rownames(club), value=TRUE))
matches.exclude <- NULL
matches.include <- matches[!(matches %in% matches.exclude)]
matches.include <- matches.include[matches.include != 'IL6ST']
matches.include <- sort(matches.include,decreasing = T)
genes.use <- matches.include
# Yap downstream
genes.use <- unique(c('YAP1','AREG','FGF1','CTGF',
                    'SMAD7','AXIN2','SERPINE1','ITGB2','GLI1','BBC3',
                    'AFP','ID1','ID2','NKD1','MYC','CCND1',
                    'SOX2','SNAI2','BIRC2','BIRC5'))

DefaultAssay(club) <- 'RNA'
Idents(club) <- club[['Condition']]
club$Condition <- factor(x=club$Condition,
                             levels = c('Mock','1dpi','2dpi','3dpi'))
club <- club[order(club$Condition),]
identities <- c('Mock','1dpi','2dpi','3dpi')

# prepare data
club <- ScaleData(club,features = genes.use)  # $gene
club.avg <- AverageExpression(club,assays = 'RNA',slot='scale.data')  # ,features=genes.use
club.avg <- t(scale(t(club.avg$RNA)))
unityNormalize <- function(x){
  (x-min(x))/(max(x)-min(x))
}
club.avg <- t(apply(as.matrix(club.avg), 1, unityNormalize))
club.avg <- club.avg[,c('Mock','1dpi','2dpi','3dpi')]
dim(club.avg)
# [1] 114   4

#  club.avg <- club.avg[-c(3), ]
# Prep & plot
colors.inferno <- colorRamp2(breaks = c(seq(min(club.avg),max(club.avg),length.out=60)), inferno(n=60), space = "RGB")
seuratcolors <- hue_pal()(4)

club.heatmap <- Heatmap(as.matrix(club.avg),
                        col= colors.inferno,
                        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=seuratcolors))),
                        column_split= factor(colnames(club.avg), levels = identities),
                        column_title = identities,
                        column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 18),
                        column_title_rot = 90,
                        cluster_rows=T,
                        row_km = 4,
                        cluster_columns=F,
                        cluster_column_slices=F,
                        show_column_names=FALSE,
                        show_row_names = T,
                        show_heatmap_legend = F)

png('club_timehm_3.png', width = 250, height = 1700)
draw(club.heatmap,padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()




#### Pulling and observing specific gene sets
Idents(hbec.integrated) <- hbec.integrated[['Condition']]
hbec.integrated$Condition <- factor(hbec.integrated$Condition, levels=c('Mock','1dpi','2dpi','3dpi'))

# manually remove Mock infected contaminates
Idents(hbec.integrated) <- hbec.integrated[['scv2+']]
infected <- subset(hbec.integrated, idents = '1')
Idents(infected) <- infected[['Condition']]
mockinfected <- subset(infected, idents = 'Mock')
hbec.integrated <- subset(hbec.integrated, cells = Cells(mockinfected),invert = T)

## KEGG Notch Pathway Gene List
notch.genes <- c('ADAM17','APH1A','CIR1','CREBBP','CTBP1','CTBP2','DLL1','DLL3','DLL4','DTX1','DTX2','DTX3',
                'DTX3L','DTX4','DVL1','DVL2','DVL3','EP300','HDAC1','HDAC2','HES1','HES5','JAG1','JAG2',
                'KAT2A','KAT2B','LFNG','MAML1','MAML2','MAML3','MFNG','NCOR2','NCSTN','NOTCH1','NOTCH2',
                'NOTCH3','NOTCH4','NUMB','NUMBL','PSEN1','PSEN2','PSENEN','PTCRA','RBPJ','RBPJL','RFNG','SNW1')
notch.genes <- intersect(notch.genes,rownames(hbec.integrated))
notch.genes <- c(notch.genes,'notch.score1')

hippo.genes <- c('AJUBA','AMOT','AMOTL1','AMOTL2','CASP3','DCHS1','DLG5','DVL2','FAT4','LATS1','LATS2',
                'LIMD1','MARK3','MOB1A','MOB1B','MOB3B','NEK8','NF2','NPHP4','PJA2','SAV1','SOX11','STK3',
                'STK4','TEAD1','TEAD2','TEAD3','TEAD4','TJP1','TJP2','WTIP','WWC1','WWC2','WWC3','WWTR1',
                'YAP1','YWHAB','YWHAE')
hippo.genes <- intersect(hippo.genes,rownames(hbec.integrated))
hippo.genes <- c(hippo.genes,'hippo.score1')

hbec.integrated <- ScaleData(hbec.integrated,features = notch.genes)
## ViolinPlot on bulk
Idents(hbec.integrated) <- hbec.integrated[['Condition']]
hbec.integrated <- RenameIdents(hbec.integrated,'Mock'='M','1dpi'='1','2dpi'='2','3dpi'='3')
hbec.integrated[['Condition2']] <- Idents(hbec.integrated)
Idents(hbec.integrated) <- hbec.integrated[['Condition2']]

myplots <- vector('list', length(notch.genes))
for (i in 1:length(notch.genes)) {
    message(i)
    gene <- notch.genes[i]
    myplots[[i]] <- local({
        i <- i
        p1 <- VlnPlot(hbec.integrated,gene,pt.size = 0,group.by='Condition2') + xlab('') + ylab('') + NoLegend() +
            theme(axis.text.x=element_text(size=32, color = 'black',angle = 0,hjust = 0.5),
            axis.text.y = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(size = 38, hjust = 0.5))
    })
}
#,split.plot=T,split.by='scv2+'
p1 <- plot_grid(plotlist=myplots,ncol=6)
png('hbec.integrated_vln_notch.png', width = 1333,height=1778)
p1
dev.off()

#sort( sapply(mget(ls()),object.size) )
#sort( sapply(ls(),function(x){object.size(get(x))}))


### AddModuleScores
hbec.integrated <- AddModuleScore(hbec.integrated,features=list(notch.genes),ctrl = 5,name='notch.score')
hbec.integrated <- AddModuleScore(hbec.integrated,features=list(hippo.genes),ctrl = 5,name='hippo.score')



#### Running on ciliated
load('./ciliated2020-10-09.Robj')
Idents(ciliated) <- ciliated[['Condition']]
ciliated$Condition <- factor(ciliated$Condition, levels=c('Mock','1dpi','2dpi','3dpi'))

### AddModuleScores
ciliated <- AddModuleScore(ciliated,features=list(notch.genes),ctrl = 5,name='notch.score')
ciliated <- AddModuleScore(ciliated,features=list(hippo.genes),ctrl = 5,name='hippo.score')

# manually remove Mock infected contaminates
Idents(ciliated) <- ciliated[['scv2+']]
infected <- subset(ciliated, idents = '1')
Idents(infected) <- infected[['Condition']]
mockinfected <- subset(infected, idents = 'Mock')
ciliated <- subset(ciliated, cells = Cells(mockinfected),invert = T)

## KEGG Notch Pathway Gene List
notch.genes <- c('ADAM17','APH1A','CIR1','CREBBP','CTBP1','CTBP2','DLL1','DLL3','DLL4','DTX1','DTX2','DTX3',
                'DTX3L','DTX4','DVL1','DVL2','DVL3','EP300','HDAC1','HDAC2','HES1','HES5','JAG1','JAG2',
                'KAT2A','KAT2B','LFNG','MAML1','MAML2','MAML3','MFNG','NCOR2','NCSTN','NOTCH1','NOTCH2',
                'NOTCH3','NOTCH4','NUMB','NUMBL','PSEN1','PSEN2','PSENEN','PTCRA','RBPJ','RBPJL','RFNG','SNW1')
notch.genes <- intersect(notch.genes,rownames(ciliated))
notch.genes <- c(notch.genes,'notch.score1')

hippo.genes <- c('AJUBA','AMOT','AMOTL1','AMOTL2','CASP3','DCHS1','DLG5','DVL2','FAT4','LATS1','LATS2',
                'LIMD1','MARK3','MOB1A','MOB1B','MOB3B','NEK8','NF2','NPHP4','PJA2','SAV1','SOX11','STK3',
                'STK4','TEAD1','TEAD2','TEAD3','TEAD4','TJP1','TJP2','WTIP','WWC1','WWC2','WWC3','WWTR1',
                'YAP1','YWHAB','YWHAE')
hippo.genes <- intersect(hippo.genes,rownames(ciliated))
hippo.genes <- c(hippo.genes,'hippo.score1')

ciliated <- ScaleData(ciliated,features = hippo.genes)
## ViolinPlot on bulk
Idents(ciliated) <- ciliated[['Condition']]
ciliated <- RenameIdents(ciliated,'Mock'='M','1dpi'='1','2dpi'='2','3dpi'='3')
ciliated[['Condition2']] <- Idents(ciliated)
Idents(ciliated) <- ciliated[['Condition2']]

myplots <- vector('list', length(hippo.genes))
for (i in 1:length(hippo.genes)) {
    message(i)
    gene <- hippo.genes[i]
    myplots[[i]] <- local({
        i <- i
        p1 <- VlnPlot(ciliated,gene,pt.size = 0,group.by='Condition2') + xlab('') + ylab('') + NoLegend() +
            theme(axis.text.x=element_text(size=32, color = 'black',angle = 0,hjust = 0.5),
            axis.text.y = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(size = 38, hjust = 0.5))
    })
}
#,split.plot=T,split.by='scv2+'
p1 <- plot_grid(plotlist=myplots,ncol=6)
png('ciliated_vln_hippo.png', width = 1333,height=1556)  # 1556, 1778
p1
dev.off()


#### Running on basal
load('./basal2020-06-26.Robj')
Idents(basal) <- basal[['Condition']]
basal$Condition <- factor(basal$Condition, levels=c('Mock','1dpi','2dpi','3dpi'))
### AddModuleScores
basal <- AddModuleScore(basal,features=list(notch.genes),ctrl = 5,name='notch.score')
basal <- AddModuleScore(basal,features=list(hippo.genes),ctrl = 5,name='hippo.score')
# manually remove Mock infected contaminates
Idents(basal) <- basal[['scv2+']]
infected <- subset(basal, idents = '1')
Idents(infected) <- infected[['Condition']]  # no mock infected
## KEGG Notch Pathway Gene List
notch.genes <- c('ADAM17','APH1A','CIR1','CREBBP','CTBP1','CTBP2','DLL1','DLL3','DLL4','DTX1','DTX2','DTX3',
                'DTX3L','DTX4','DVL1','DVL2','DVL3','EP300','HDAC1','HDAC2','HES1','HES5','JAG1','JAG2',
                'KAT2A','KAT2B','LFNG','MAML1','MAML2','MAML3','MFNG','NCOR2','NCSTN','NOTCH1','NOTCH2',
                'NOTCH3','NOTCH4','NUMB','NUMBL','PSEN1','PSEN2','PSENEN','PTCRA','RBPJ','RBPJL','RFNG','SNW1')
notch.genes <- intersect(notch.genes,rownames(basal))
notch.genes <- c(notch.genes,'notch.score1')

hippo.genes <- c('AJUBA','AMOT','AMOTL1','AMOTL2','CASP3','DCHS1','DLG5','DVL2','FAT4','LATS1','LATS2',
                'LIMD1','MARK3','MOB1A','MOB1B','MOB3B','NEK8','NF2','NPHP4','PJA2','SAV1','SOX11','STK3',
                'STK4','TEAD1','TEAD2','TEAD3','TEAD4','TJP1','TJP2','WTIP','WWC1','WWC2','WWC3','WWTR1',
                'YAP1','YWHAB','YWHAE')
hippo.genes <- intersect(hippo.genes,rownames(basal))
hippo.genes <- c(hippo.genes,'hippo.score1')

basal <- ScaleData(basal,features = notch.genes)
## ViolinPlot on bulk
Idents(basal) <- basal[['Condition']]
basal <- RenameIdents(basal,'Mock'='M','1dpi'='1','2dpi'='2','3dpi'='3')
basal[['Condition2']] <- Idents(basal)
Idents(basal) <- basal[['Condition2']]

myplots <- vector('list', length(notch.genes))
for (i in 1:length(notch.genes)) {
    message(i)
    gene <- notch.genes[i]
    myplots[[i]] <- local({
        i <- i
        p1 <- VlnPlot(basal,gene,pt.size = 0,group.by='Condition2') + xlab('') + ylab('') + NoLegend() +
            theme(axis.text.x=element_text(size=32, color = 'black',angle = 0,hjust = 0.5),
            axis.text.y = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(size = 38, hjust = 0.5))
    })
}
#,split.plot=T,split.by='scv2+'
p1 <- plot_grid(plotlist=myplots,ncol=6)
png('basal_vln_notch2.png', width = 1333,height=1778)  # 1556, 1778
p1
dev.off()



#### Running on intermediate
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
intermediate <- subset(hbec.integrated,idents = 'Intermediate')
Idents(intermediate) <- intermediate[['Condition']]
intermediate$Condition <- factor(intermediate$Condition, levels=c('Mock','1dpi','2dpi','3dpi'))
### AddModuleScores
intermediate <- AddModuleScore(intermediate,features=list(notch.genes),ctrl = 5,name='notch.score')
intermediate <- AddModuleScore(intermediate,features=list(hippo.genes),ctrl = 5,name='hippo.score')
# manually remove Mock infected contaminates
Idents(intermediate) <- intermediate[['scv2+']]
infected <- subset(intermediate, idents = '1')
Idents(infected) <- infected[['Condition']]  # no mock infected
## KEGG Notch Pathway Gene List
notch.genes <- c('ADAM17','APH1A','CIR1','CREBBP','CTBP1','CTBP2','DLL1','DLL3','DLL4','DTX1','DTX2','DTX3',
                'DTX3L','DTX4','DVL1','DVL2','DVL3','EP300','HDAC1','HDAC2','HES1','HES5','JAG1','JAG2',
                'KAT2A','KAT2B','LFNG','MAML1','MAML2','MAML3','MFNG','NCOR2','NCSTN','NOTCH1','NOTCH2',
                'NOTCH3','NOTCH4','NUMB','NUMBL','PSEN1','PSEN2','PSENEN','PTCRA','RBPJ','RBPJL','RFNG','SNW1')
notch.genes <- intersect(notch.genes,rownames(intermediate))
notch.genes <- c(notch.genes,'notch.score1')

hippo.genes <- c('AJUBA','AMOT','AMOTL1','AMOTL2','CASP3','DCHS1','DLG5','DVL2','FAT4','LATS1','LATS2',
                'LIMD1','MARK3','MOB1A','MOB1B','MOB3B','NEK8','NF2','NPHP4','PJA2','SAV1','SOX11','STK3',
                'STK4','TEAD1','TEAD2','TEAD3','TEAD4','TJP1','TJP2','WTIP','WWC1','WWC2','WWC3','WWTR1',
                'YAP1','YWHAB','YWHAE')
hippo.genes <- intersect(hippo.genes,rownames(intermediate))
hippo.genes <- c(hippo.genes,'hippo.score1')

intermediate <- ScaleData(intermediate,features = hippo.genes)
## ViolinPlot on bulk
Idents(intermediate) <- intermediate[['Condition']]
intermediate <- RenameIdents(intermediate,'Mock'='M','1dpi'='1','2dpi'='2','3dpi'='3')
intermediate[['Condition2']] <- Idents(intermediate)
Idents(intermediate) <- intermediate[['Condition2']]

myplots <- vector('list', length(hippo.genes))
for (i in 1:length(hippo.genes)) {
    message(i)
    gene <- hippo.genes[i]
    myplots[[i]] <- local({
        i <- i
        p1 <- VlnPlot(intermediate,gene,pt.size = 0,group.by='Condition2') + xlab('') + ylab('') + NoLegend() +
            theme(axis.text.x=element_text(size=32, color = 'black',angle = 0,hjust = 0.5),
            axis.text.y = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(size = 38, hjust = 0.5))
    })
}
# ,split.plot=T,split.by='scv2+'
p1 <- plot_grid(plotlist=myplots,ncol=6)
png('intermediate_vln_hippo2.png', width = 1333,height=1556)  # 1556, 1778
p1
dev.off()


#### Running on club
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
club <- subset(hbec.integrated,idents = 'Club')
Idents(club) <- club[['Condition']]
club$Condition <- factor(club$Condition, levels=c('Mock','1dpi','2dpi','3dpi'))
### AddModuleScores
club <- AddModuleScore(club,features=list(notch.genes),ctrl = 5,name='notch.score')
club <- AddModuleScore(club,features=list(hippo.genes),ctrl = 5,name='hippo.score')
# manually remove Mock infected contaminates
Idents(club) <- club[['scv2+']]
infected <- subset(club, idents = '1')
Idents(infected) <- infected[['Condition']]  # no mock infected
## KEGG Notch Pathway Gene List
notch.genes <- c('ADAM17','APH1A','CIR1','CREBBP','CTBP1','CTBP2','DLL1','DLL3','DLL4','DTX1','DTX2','DTX3',
                'DTX3L','DTX4','DVL1','DVL2','DVL3','EP300','HDAC1','HDAC2','HES1','HES5','JAG1','JAG2',
                'KAT2A','KAT2B','LFNG','MAML1','MAML2','MAML3','MFNG','NCOR2','NCSTN','NOTCH1','NOTCH2',
                'NOTCH3','NOTCH4','NUMB','NUMBL','PSEN1','PSEN2','PSENEN','PTCRA','RBPJ','RBPJL','RFNG','SNW1')
notch.genes <- intersect(notch.genes,rownames(club))
notch.genes <- c(notch.genes,'notch.score1')

hippo.genes <- c('AJUBA','AMOT','AMOTL1','AMOTL2','CASP3','DCHS1','DLG5','DVL2','FAT4','LATS1','LATS2',
                'LIMD1','MARK3','MOB1A','MOB1B','MOB3B','NEK8','NF2','NPHP4','PJA2','SAV1','SOX11','STK3',
                'STK4','TEAD1','TEAD2','TEAD3','TEAD4','TJP1','TJP2','WTIP','WWC1','WWC2','WWC3','WWTR1',
                'YAP1','YWHAB','YWHAE')
hippo.genes <- intersect(hippo.genes,rownames(club))
hippo.genes <- c(hippo.genes,'hippo.score1')

club <- ScaleData(club,features = notch.genes)
## ViolinPlot on bulk
Idents(club) <- club[['Condition']]
club <- RenameIdents(club,'Mock'='M','1dpi'='1','2dpi'='2','3dpi'='3')
club[['Condition2']] <- Idents(club)
Idents(club) <- club[['Condition2']]

myplots <- vector('list', length(notch.genes))
for (i in 1:length(notch.genes)) {
    message(i)
    gene <- notch.genes[i]
    myplots[[i]] <- local({
        i <- i
        p1 <- VlnPlot(club,gene,pt.size = 0,group.by='Condition2') + xlab('') + ylab('') + NoLegend() +
            theme(axis.text.x=element_text(size=32, color = 'black',angle = 0,hjust = 0.5),
            axis.text.y = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(size = 38, hjust = 0.5))
    })
}
#,split.plot=T,split.by='scv2+'
p1 <- plot_grid(plotlist=myplots,ncol=6)
png('club_vln_notch2.png', width = 1333,height=1778)  # 1556, 1778
p1
dev.off()











## ggplot2 violin plot on bulk
DefaultAssay(hbec.integrated) <- 'RNA'
Idents(hbec.integrated) <- hbec.integrated[['Condition']]
hbec.integrated$Condition <- factor(hbec.integrated$Condition, levels=c('Mock','1dpi','2dpi','3dpi'))
hbec.integrated <- hbec.integrated[order(hbec.integrated$Condition),]

gene = c('JAG1')
hbec.integrated <- ScaleData(hbec.integrated,features = gene)
frame1 <- t(data.frame(GetAssayData(hbec.integrated, slot = "scale.data")))
colnames(frame1) <- c('Gene')
rownames(frame1) <- NULL
frame2 <- data.frame(hbec.integrated$Condition)
colnames(frame2) <- c('Condition')
rownames(frame2) <- NULL
frame <- cbind(frame1,frame2)

p <- ggplot(frame, aes(Condition,Gene)) +
    geom_violin(trim=T) + labs(title = gene) +
    theme(axis.text.x=element_text(size=32, color = 'black'),
        axis.text.y = element_text(size = 32, color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 38, hjust = 0.5))

png('hbec.integrated_vln_JAG1.png', width = 400,height=400)
p
dev.off()



#### Heatmaps of my gene lists
DefaultAssay(hbec.integrated) <- 'RNA'
Idents(hbec.integrated) <- hbec.integrated[['Condition']]
hbec.integrated <- RenameIdents(hbec.integrated,'Mock'='M','1dpi'='1','2dpi'='2','3dpi'='3')
hbec.integrated[['Condition2']] <- Idents(hbec.integrated)
Idents(hbec.integrated) <- hbec.integrated[['Condition2']]
hbec.integrated$Condition2 <- factor(x=hbec.integrated$Condition2,
                             levels = c('M','1','2','3'))
hbec.integrated <- hbec.integrated[order(hbec.integrated$Condition2),]
identities <- c('M','1','2','3')

notch.genes <- c('ADAM17','APH1A','CIR1','CREBBP','CTBP1','CTBP2','DLL1','DLL3','DLL4','DTX1','DTX2','DTX3',
                'DTX3L','DTX4','DVL1','DVL2','DVL3','EP300','HDAC1','HDAC2','HES1','HES5','JAG1','JAG2',
                'KAT2A','KAT2B','LFNG','MAML1','MAML2','MAML3','MFNG','NCOR2','NCSTN','NOTCH1','NOTCH2',
                'NOTCH3','NOTCH4','NUMB','NUMBL','PSEN1','PSEN2','PSENEN','PTCRA','RBPJ','RBPJL','RFNG','SNW1')
notch.genes <- intersect(notch.genes,rownames(hbec.integrated))

hippo.genes <- c('AJUBA','AMOT','AMOTL1','AMOTL2','CASP3','DCHS1','DLG5','DVL2','FAT4','LATS1','LATS2',
                'LIMD1','MARK3','MOB1A','MOB1B','MOB3B','NEK8','NF2','NPHP4','PJA2','SAV1','SOX11','STK3',
                'STK4','TEAD1','TEAD2','TEAD3','TEAD4','TJP1','TJP2','WTIP','WWC1','WWC2','WWC3','WWTR1',
                'YAP1','YWHAB','YWHAE')
hippo.genes <- intersect(hippo.genes,rownames(hbec.integrated))

hbec.integrated <- ScaleData(hbec.integrated,features = rownames(hbec.integrated))
hbec.integrated.avg <- AverageExpression(hbec.integrated,assays = 'RNA',slot='scale.data',features=hippo.genes)  # ,features=genes.use
hbec.integrated.avg <- t(scale(t(hbec.integrated.avg$RNA)))
unityNormalize <- function(x){
  (x-min(x))/(max(x)-min(x))
}
hbec.integrated.avg <- t(apply(as.matrix(hbec.integrated.avg), 1, unityNormalize))
hbec.integrated.avg <- hbec.integrated.avg[,c('M','1','2','3')]
dim(hbec.integrated.avg)

# Prep & plot
colors.inferno <- colorRamp2(breaks = c(seq(min(hbec.integrated.avg),max(hbec.integrated.avg),length.out=60)), inferno(n=60), space = "RGB")
seuratcolors <- hue_pal()(4)

hbec.integrated.heatmap <- Heatmap(as.matrix(hbec.integrated.avg),
                        col= colors.inferno,
                        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=seuratcolors))),
                        column_split= factor(colnames(hbec.integrated.avg), levels = identities),
                        column_title = identities,
                        column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 18),
                        cluster_rows=T,
                        row_km = 4,
                        cluster_columns=F,
                        cluster_column_slices=F,
                        show_column_names=FALSE,
                        show_row_names = T,
                        show_heatmap_legend = F)

png('hbec.integrated_hm_hippo2.png', width = 250, height = 1000)
draw(hbec.integrated.heatmap,padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()



#### Making statistically significant plots - TIGHT LOOP SETUP
# load Sam's code 'HBECStatSig.R' to get function HBECStatSig
# load all objects & remove mock infected cells from hbec.integrated and ciliated
#Idents(hbec.integrated) <- hbec.integrated[['scv2+']]
#infected <- subset(hbec.integrated, idents = '1')
#Idents(infected) <- infected[['Condition']]
#mockinfected <- subset(infected, idents = 'Mock')
#hbec.integrated <- subset(hbec.integrated, cells = Cells(mockinfected),invert = T)
load('./ciliated2020-10-09.Robj')
#Idents(ciliated) <- ciliated[['scv2+']]
#infected <- subset(ciliated, idents = '1')
#Idents(infected) <- infected[['Condition']]
#mockinfected <- subset(infected, idents = 'Mock')
#ciliated <- subset(ciliated, cells = Cells(mockinfected),invert = T)
load('./basal2020-06-26.Robj')
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
intermediate <- subset(hbec.integrated,idents = 'Intermediate')
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
club <- subset(hbec.integrated,idents = 'Club')

# Order data & rename short conditions
#Idents(basal) <- basal[['Condition']]
#basal$Condition <- factor(basal$Condition, levels=c('Mock','1dpi','2dpi','3dpi'))
#basal <- RenameIdents(basal,'Mock'='M','1dpi'='1','2dpi'='2','3dpi'='3')
#basal[['Condition2']] <- Idents(basal)
Idents(hbec.integrated) <- hbec.integrated[['Condition2']]
# First pass gene sets
notch.genes <- c('ADAM17','APH1A','CIR1','CREBBP','CTBP1','CTBP2','DLL1','DLL3','DLL4','DTX1','DTX2','DTX3',
                'DTX3L','DTX4','DVL1','DVL2','DVL3','EP300','HDAC1','HDAC2','HES1','HES5','JAG1','JAG2',
                'KAT2A','KAT2B','LFNG','MAML1','MAML2','MAML3','MFNG','NCOR2','NCSTN','NOTCH1','NOTCH2',
                'NOTCH3','NOTCH4','NUMB','NUMBL','PSEN1','PSEN2','PSENEN','PTCRA','RBPJ','RBPJL','RFNG','SNW1')
notch.genes <- intersect(notch.genes,rownames(hbec.integrated))
hippo.genes <- c('AJUBA','AMOT','AMOTL1','AMOTL2','CASP3','DCHS1','DLG5','DVL2','FAT4','LATS1','LATS2',
                'LIMD1','MARK3','MOB1A','MOB1B','MOB3B','NEK8','NF2','NPHP4','PJA2','SAV1','SOX11','STK3',
                'STK4','TEAD1','TEAD2','TEAD3','TEAD4','TJP1','TJP2','WTIP','WWC1','WWC2','WWC3','WWTR1',
                'YAP1','YWHAB','YWHAE')
hippo.genes <- intersect(hippo.genes,rownames(hbec.integrated))
tgfb.genes <- c('ACVR1','ACVR1C','ACVR2A','ACVR2B','ACVRL1','AMH','AMHR2','BMP2','BMP4','BMP5','BMP6','BMP7',
                'BMP8A','BMP8B','BMPR1A','BMPR1B','BMPR2','CDKN2B','CHRD','COMP','CREBBP','CUL1','DCN','E2F4',
                'E2F5','EP300','FST','GDF5','GDF6','GDF7','ID1','ID2','ID3','ID4','IFNG','INHBA','INHBB',
                'INHBC','INHBE','LEFTY1','LEFTY2','LTBP1','MAPK1','MAPK3','MYC','NODAL','NOG','PITX2','PPP2CA',
                'PPP2CB','PPP2R1A','PPP2R1B','RBL1','RBL2','RBX1','RHOA','ROCK1','ROCK2','RPS6KB1','RPS6KB2',
                'SKP1','SMAD1','SMAD2','SMAD3','SMAD4','SMAD5','SMAD6','SMAD7','SMAD9','SMURF1','SMURF2','SP1',
                'TFDP1','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','THBS1','THBS2','THBS3','THBS4','TNF',
                'ZFYVE16','ZFYVE9')
tgfb.genes <- intersect(tgfb.genes,rownames(hbec.integrated))
wnt.genes <- c('APC','APC2','AXIN1','AXIN2','BTRC','CACYBP','CAMK2A','CAMK2B','CAMK2D','CAMK2G','CCND1','CCND2',
                'CCND3','CER1','CHD8','CHP1','CHP2','CREBBP','CSNK1A1','CSNK1A1L','CSNK1E','CSNK2A1',
                'CSNK2A2','CSNK2B','CTBP1','CTBP2','CTNNB1','CTNNBIP1','CUL1','CXXC4','DAAM1','DAAM2','DKK1',
                'DKK2','DKK4','DVL1','DVL2','DVL3','EP300','FBXW11','FOSL1','FRAT1','FRAT2','FZD1','FZD10',
                'FZD2','FZD3','FZD4','FZD5','FZD6','FZD7','FZD8','FZD9','GSK3B','JUN','LEF1','LRP5','LRP6',
                'MAP3K7','MAPK10','MAPK8','MAPK9','MMP7','MYC','NFAT5','NFATC1','NFATC2','NFATC3','NFATC4',
                'NKD1','NKD2','NLK','PLCB1','PLCB2','PLCB3','PLCB4','PORCN','PPARD','PPP2CA','PPP2CB',
                'PPP2R1A','PPP2R1B','PPP2R5A','PPP2R5B','PPP2R5C','PPP2R5D','PPP2R5E','PPP3CA','PPP3CB',
                'PPP3CC','PPP3R1','PPP3R2','PRICKLE1','PRICKLE2','PRKACA','PRKACB','PRKACG','PRKCA','PRKCB',
                'PRKCG','PRKX','PSEN1','RAC1','RAC2','RAC3','RBX1','RHOA','ROCK1','ROCK2','RUVBL1','SENP2',
                'SFRP1','SFRP2','SFRP4','SFRP5','SIAH1','SKP1','SMAD2','SMAD3','SMAD4','SOX17','TBL1X',
                'TBL1XR1','TBL1Y','TCF7','TCF7L1','TCF7L2','TP53','VANGL1','VANGL2','WIF1','WNT1','WNT10A',
                'WNT10B','WNT11','WNT16','WNT2','WNT2B','WNT3','WNT3A','WNT4','WNT5A','WNT5B','WNT6','WNT7A',
                'WNT7B','WNT8A','WNT8B','WNT9A','WNT9B')
wnt.genes <- intersect(wnt.genes,rownames(hbec.integrated))
goi.genes <- c('WNT5A','WNT4','TGFB3','TGFB1','BMP3','BMP4','ID2','ID3','BMPR1A','BMPR1B','BMPR2','RGMB',
                'AREG','EREG','VEGFA','ERBB2','HBEGF')
met.genes <- c('TGFBR2','EGR1','HK2','LDHA','PGK1','PRDX5','GSTA1','ATP5ME','NFE2L2')
genes.use <- c('NFE2L2','GSTA1','GSTP1','ACE2')
# Sam's function to filter genes
notch.use <- HBECStatSig(hbec.integrated,notch.genes,min.pct = 0.1,logfc.threshold = 0)
sig.notch.genes <- sort(unique(notch.use$gene))
hippo.use <- HBECStatSig(hbec.integrated,hippo.genes,min.pct = 0.1,logfc.threshold = 0)
sig.hippo.genes <- sort(unique(hippo.use$gene))
tgfb.use <- HBECStatSig(hbec.integrated,tgfb.genes,min.pct = 0.1,logfc.threshold = 0)
sig.tgfb.genes <- sort(unique(tgfb.use$gene))
wnt.use <- HBECStatSig(hbec.integrated,wnt.genes,min.pct = 0.1,logfc.threshold = 0)
sig.wnt.genes <- sort(unique(wnt.use$gene))
# Calculate and add module scores
hbec.integrated <- AddModuleScore(hbec.integrated,features=list(sig.notch.genes),ctrl = 5,name='notch.score')
hbec.integrated <- AddModuleScore(hbec.integrated,features=list(sig.hippo.genes),ctrl = 5,name='hippo.score')
genes.use <- c(sig.notch.genes,'notch.score1')
genes.use <- c(sig.hippo.genes,'hippo.score1')
hbec.integrated <- AddModuleScore(hbec.integrated,features=list(sig.tgfb.genes),ctrl = 5,name='tgfb.score')
genes.use <- c(sig.tgfb.genes,'tgfb.score1')
hbec.integrated <- AddModuleScore(hbec.integrated,features=list(sig.wnt.genes),ctrl = 5,name='wnt.score')
genes.use <- c(sig.wnt.genes,'wnt.score1')
connectome <- c('TGFB2','BMP4','WNT5A','FST','AREG','EREG','PLAU','SERPINE1','JAG1','TGFBR2','THBS1','CTGF','ITGAV','BMPR1A','ID3','EGFR','BMPR2','NOTCH1','SMAD1','SMAD2','SMAD3','SMAD5','GSTA1','GSTP1','TXN','ACE2','YAP1','WWTR1','SEMA3E')
# Plot violin plots
myplots <- vector('list', length(genes.use))
for (i in 1:length(genes.use)) {
    message(i)
    gene <- genes.use[i]
    myplots[[i]] <- local({
        i <- i
        p1 <- VlnPlot(hbec.integrated,gene,pt.size = 0,group.by='Condition2') + xlab('') + ylab('') + NoLegend() +
            theme(axis.text.x=element_text(size=32, color = 'black',angle = 0,hjust = 0.5),
            axis.text.y = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(size = 38, hjust = 0.5))
    })
}
#,split.plot=T,split.by='scv2+'
ht <- ceiling(length(genes.use)/6)*222.25
p1 <- plot_grid(plotlist=myplots,ncol=6)
png('hbec.integrated_vln_conn2.png', width = 1333,height=ht)  # 1556, 1778
p1
dev.off()

# Plot heatmap
hbec.integrated <- ScaleData(hbec.integrated,features = rownames(hbec.integrated))
hbec.integrated.avg <- AverageExpression(hbec.integrated,assays = 'RNA',slot='scale.data',features=genes.use)
hbec.integrated.avg <- t(scale(t(hbec.integrated.avg$RNA)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
hbec.integrated.avg <- t(apply(as.matrix(hbec.integrated.avg), 1, unityNormalize))
hbec.integrated.avg <- hbec.integrated.avg[,c('M','1','2','3')]
identities <- c('M','1','2','3')
dim(hbec.integrated.avg)
colors.inferno <- colorRamp2(breaks = c(seq(min(hbec.integrated.avg),max(hbec.integrated.avg),length.out=60)), inferno(n=60), space = "RGB")
seuratcolors <- hue_pal()(4)
hbec.integrated.heatmap <- Heatmap(as.matrix(hbec.integrated.avg),
                        col= colors.inferno,
                        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=seuratcolors))),
                        column_split= factor(colnames(hbec.integrated.avg), levels = identities),
                        column_title = identities,
                        column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 18),
                        cluster_rows=T,
                        row_km = 3,
                        cluster_columns=F,
                        cluster_column_slices=F,
                        show_column_names=FALSE,
                        show_row_names = T,
                        show_heatmap_legend = F)
png('hbec.integrated_hm_conn.png', width = 250, height = 850)  # was 1600
draw(hbec.integrated.heatmap,padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()




#### Do IFIT genes correlate with infection?
# Load organize and prepare sample
Idents(ciliated) <- ciliated[['Condition']]
ciliated$Condition <- factor(ciliated$Condition, levels=c('Mock','1dpi','2dpi','3dpi'))
ciliated <- RenameIdents(ciliated,'Mock'='M','1dpi'='1','2dpi'='2','3dpi'='3')
ciliated[['Condition2']] <- Idents(ciliated)
Idents(ciliated) <- ciliated[['Condition2']]
Idents(ciliated) <- ciliated[['ciliated2']]
ciliated$ciliated2 <- factor(ciliated$ciliated2, levels=c('Progenitor','Immune_Responsive','Mature_Functional','Novel_Infected'))
ciliated <- RenameIdents(ciliated,'Progenitor'='P','Immune_Responsive'='IR','Mature_Functional'='MF','Novel_Infected'='NI')
ciliated[['ciliated3']] <- Idents(ciliated)
Idents(ciliated) <- ciliated[['ciliated3']]
# Identify IFIT genes expressed in the sample and filter for significance
ifit.genes <- unique(grep(paste(c('IFIT'),collapse="|"),rownames(ciliated), value=TRUE))
ifit.use <- HBECStatSig(ciliated,ifit.genes)
sig.ifit.genes <- sort(unique(ifit.use$gene))
ciliated <- AddModuleScore(ciliated,features=list(sig.ifit.genes),ctrl = 5,name='ifit.score')
genes.use <- c(sig.ifit.genes,'ifit.score1')
myplots <- vector('list', length(genes.use))
for (i in 1:length(genes.use)) {
    message(i)
    gene <- genes.use[i]
    myplots[[i]] <- local({
        i <- i
        p1 <- VlnPlot(ciliated,gene,pt.size = 0,group.by='ciliated3',split.plot=T,split.by='scv2+') + xlab('') + ylab('') + NoLegend() +
            theme(axis.text.x=element_text(size=32, color = 'black',angle = 0,hjust = 0.5),
            axis.text.y = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(size = 38, hjust = 0.5))
    })
}
ht <- ceiling(length(genes.use)/6)*222.25
p1 <- plot_grid(plotlist=myplots,ncol=6)
png('ciliated_vln_ifit2.png', width = 1333,height=ht)  # 1556, 1778
p1
dev.off()
genes.use <- c(sig.ifit.genes,'scv2+')
myplots <- vector('list', length(genes.use))
for (i in 1:length(genes.use)) {
    message(i)
    gene <- genes.use[i]
    myplots[[i]] <- local({
        i <- i
        p1 <- FeaturePlot(ciliated,gene,label = T,label.size=8,repel=T) + labs(title = paste(gene)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"),
            plot.title = element_text(size = 38, hjust = 0.5))
    })
}
ht <- ceiling(length(genes.use)/6)*222.25*4
p1 <- plot_grid(plotlist=myplots,ncol=3)
png('ciliated_fp_ifit2.png', width = 1333,height=ht)
p1
dev.off()


#### Monocle trajectory gene plots
load('./hbecmono2020-07-30.Robj')
load('./ciliamono2020-07-29.Robj')
# Perform BEAM on ciliated branches
BEAM_res <- BEAM(ciliamono, branch_point = 1, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
png('ciliamono_beam_1.png', width = 800,height=1000)
plot_genes_branched_heatmap(ciliamono[row.names(subset(BEAM_res,
                                          qval < 1e-4)),],
                                          branch_point = 1,
                                          branch_states = c(1,5),
                                          branch_labels = c("Immune_Responsive", "Mature_Functional"),
                                          num_clusters = 4,
                                          cores = 4,
                                          use_gene_short_name = T,
                                          show_rownames = T)
dev.off()

Idents(hbecmono) <- hbecmono[['cell_type2']]
hbecmono[['cell_type2']] <- Idents(hbecmono)



to_be_tested <- row.names(subset(fData(hbecmono),
    gene_short_name %in% c("KRT5", "LYPD2", "DNAH6")))
cds_subset <- hbecmono[to_be_tested,]
colnames(pData(hbecmono))
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores=4)
diff_test_res[,c("gene_short_name", "pval", "qval")]
png('hbecmono_genes_2.png', width = 800,height=1200)
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type2")
dev.off()



#### Monocle3
## Create new monocle object - RNA slot
hbecsubset <- subset(hbec.integrated,idents=c('Basal','Intermediate','Club','Ciliated'))
data <- as(as.matrix(hbecsubset@assays$RNA@data), 'sparseMatrix')
cell_metadata <- as.data.frame(hbecsubset@meta.data, row.names = colnames(data))
gene_metadata <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
hbecmono <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
## Pre-process the data
hbecmono <- preprocess_cds(hbecmono, num_dim = 3,alignment_group='Condition')
hbecmono <- align_cds(hbecmono,alignment_group='Condition')
hbecmono <- reduce_dimension(hbecmono,reduction_method='UMAP',cores=4,umap.fast_sgd=T)
png('hbecmono_umap_1.png', width = 800, height = 800)
plot_cells(hbecmono,color_cells_by='cell_type2', show_trajectory_graph = F,group_label_size = 8)
dev.off()
png('hbecmono_umap_2.png', width = 800, height = 800)
plot_cells(hbecmono,color_cells_by='Condition', show_trajectory_graph = F,group_label_size = 8)
dev.off()
png('hbecmono_umap_3.png', width = 800, height = 800)
plot_cells(hbecmono,color_cells_by='cell_type4', show_trajectory_graph = F,group_label_size = 8)
dev.off()
genes <- c("KRT5","LYPD2","DNAH6",'DCLK1')
png('hbecmono_fp_1.png', width = 800, height = 800)
plot_cells(hbecmono,genes=genes,label_cell_groups=F,show_trajectory_graph=F)
dev.off()

## Cluster your cells
hbecmono <- cluster_cells(hbecmono)
png('hbecmono_umap_4.png', width = 800, height = 800)
plot_cells(hbecmono, color_cells_by = "partition",show_trajectory_graph=F)
dev.off()

levels(hbecmono@clusters$UMAP$partitions)
hbecmono@clusters$UMAP$partitions[hbecmono@clusters$UMAP$partitions == c("2")] <- "1"  #...3
levels(hbecmono@clusters$UMAP$partitions)  # same
length(hbecmono@clusters$UMAP$partitions[hbecmono@clusters$UMAP$partitions == "3"])  # empty
length(hbecmono@clusters$UMAP$partitions[hbecmono@clusters$UMAP$partitions == "1"])  # full

## Learn the trajectory graph & Order the cells in pseudotime
hbecmono <- learn_graph(hbecmono,use_partition=T)
png('hbecmono_umap_5.png', width = 800, height = 800)
plot_cells(hbecmono,color_cells_by = "cell_type2",label_groups_by_cluster=F,label_leaves=F,label_branch_points=F)
dev.off()
png('hbecmono_umap_6.png', width = 800, height = 800)
plot_cells(hbecmono,color_cells_by = "cell_type2",label_cell_groups=F,label_leaves=TRUE,label_branch_points=TRUE,graph_label_size=4)
dev.off()

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(hbecmono, cell_type="Basal"){
  cell_ids <- which(colData(hbecmono)[, "cell_type2"] == cell_type)

  closest_vertex <-
  hbecmono@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(hbecmono), ])
  root_pr_nodes <-
  igraph::V(principal_graph(hbecmono)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
hbecmono <- order_cells(hbecmono, root_pr_nodes=get_earliest_principal_node(hbecmono))

png('hbecmono_umap_7.png', width = 800, height = 800)
plot_cells(hbecmono,color_cells_by = "pseudotime",show_trajectory_graph=F,label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()

save(hbecmono, file = paste("hbecmono",Sys.Date(),".Robj",sep=""))
load('./hbecmono2020-08-10.Robj')

genes <- c("KRT5","KRT14","IGFBP6",'SCGB3A1','LYPD2','BPIFB1','AK9','DNAH12','PIFO')
hbecmono_cds <- hbecmono[rowData(hbecmono)$gene_short_name %in% genes, ]
png('hbecmono_line_2.png', width = 1000, height = 1000)
plot_genes_in_pseudotime(hbecmono_cds,
                         color_cells_by="cell_type2",
                         panel_order=genes,
                         min_expr=0.5,
                         cell_size=0,
                         ncol=3) +
                         theme(plot.title = element_text(size=25,color = 'black',hjust = 0.5),
                         axis.text.x=element_text(size=18, color = 'black',angle=0),
                         axis.text.y = element_text(size = 18, color = 'black'),
                         strip.text.x = element_text(size=25),
                         axis.line = element_line(color = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20)
                         ) + NoLegend()
dev.off()

hbecmono_cds <- hbecmono[rowData(hbecmono)$gene_short_name %in% c('scv2-orf1-10'), ]
png('hbecmono_line_scv2.png', width = 333, height = 333)
plot_genes_in_pseudotime(hbecmono_cds,
                         color_cells_by="cell_type2",
                         min_expr=0.5,
                         cell_size=0) +
                         theme(plot.title = element_text(size=25,color = 'black',hjust = 0.5),
                         axis.text.x=element_text(size=18, color = 'black',angle=0),
                         axis.text.y = element_text(size = 18, color = 'black'),
                         strip.text.x = element_text(size=25),
                         axis.line = element_line(color = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20)
                         ) + NoLegend()
dev.off()

traj.coord <- hbecmono@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
hbec.integrated <- AddMetaData(hbec.integrated, metadata = traj.coord, col.name = 'pseudotime')

png('hbec.integrated_fp_pseudotime.png', width = 600, height = 800)
FeaturePlot(hbecsubset,'pseudotime',cols=viridis(60),label=T,label.size=8) + labs(title = 'Pseudotime') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()
png('hbec.integrated_fp_pseudotime2.png', width = 1200, height = 400)
FeaturePlot(hbecsubset,'pseudotime',cols=viridis(60),label=F,split.by='Condition') +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()

bigcil <- AddMetaData(bigcil, metadata = ciliated[['pseudotime']], col.name = 'pseudotime_cil')
png('hbec.integrated_fp_pseudotimecil.png', width = 600, height = 800)
FeaturePlot(bigcil,'pseudotime_cil',cols=viridis(60),label=T,label.size=8) + labs(title = 'Pseudotime') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()
png('hbec.integrated_fp_pseudotime2cil.png', width = 1200, height = 400)
FeaturePlot(bigcil,'pseudotime_cil',cols=viridis(60),label=F,split.by='Condition') +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()



## Create new monocle object from ciliated subset
load('./ciliated2020-10-09.Robj')
data <- as(as.matrix(ciliated@assays$RNA@data), 'sparseMatrix')
cell_metadata <- as.data.frame(ciliated@meta.data, row.names = colnames(data))
gene_metadata <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
ciliamono <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
## Pre-process the data
ciliamono <- preprocess_cds(ciliamono, num_dim = 3,alignment_group='Condition')
ciliamono <- align_cds(ciliamono,alignment_group='Condition')
ciliamono <- reduce_dimension(ciliamono,reduction_method='UMAP',cores=4,umap.fast_sgd=T)
png('ciliamono_umap_1.png', width = 800, height = 800)
plot_cells(ciliamono,color_cells_by='ciliated2', show_trajectory_graph = F,group_label_size = 8)
dev.off()
png('ciliamono_umap_2.png', width = 800, height = 800)
plot_cells(ciliamono,color_cells_by='Condition', show_trajectory_graph = F,group_label_size = 8)
dev.off()
genes <- c("KRT5","SCGB1A1","DNAH12",'AK9')
png('ciliamono_fp_1.png', width = 800, height = 800)
plot_cells(ciliamono,genes=genes,label_cell_groups=F,show_trajectory_graph=F)
dev.off()

## Cluster your cells
ciliamono <- cluster_cells(ciliamono)
png('ciliamono_umap_3.png', width = 800, height = 800)
plot_cells(ciliamono, color_cells_by = "partition",show_trajectory_graph=F)
dev.off()

levels(ciliamono@clusters$UMAP$partitions)
ciliamono@clusters$UMAP$partitions[ciliamono@clusters$UMAP$partitions == c("2")] <- "1"  #...3
levels(ciliamono@clusters$UMAP$partitions)  # same
length(ciliamono@clusters$UMAP$partitions[ciliamono@clusters$UMAP$partitions == "3"])  # empty
length(ciliamono@clusters$UMAP$partitions[ciliamono@clusters$UMAP$partitions == "1"])  # full

## Learn the trajectory graph & Order the cells in pseudotime
ciliamono <- learn_graph(ciliamono,use_partition=T)
png('ciliamono_umap_4.png', width = 800, height = 800)
plot_cells(ciliamono,color_cells_by = "ciliated2",label_groups_by_cluster=F,label_leaves=F,label_branch_points=F)
dev.off()
png('ciliamono_umap_5.png', width = 800, height = 800)
plot_cells(ciliamono,color_cells_by = "ciliated2",label_cell_groups=F,label_leaves=TRUE,label_branch_points=TRUE,graph_label_size=4)
dev.off()

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(ciliamono, cell_type="Progenitor"){
  cell_ids <- which(colData(ciliamono)[, "ciliated2"] == cell_type)
  closest_vertex <-
  ciliamono@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(ciliamono), ])
  root_pr_nodes <-
  igraph::V(principal_graph(ciliamono)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
ciliamono <- order_cells(ciliamono, root_pr_nodes=get_earliest_principal_node(ciliamono))
png('ciliamono_umap_6.png', width = 800, height = 800)
plot_cells(ciliamono,color_cells_by = "pseudotime",show_trajectory_graph=F,label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()

genes <- c("KRT5","KRT14","SCGB3A1",'SAA1','PLA2G10','PRDX5','CFAP70','DNAAF1','DLEC1','IFITM2','NFKBIA','CXCL3')
ciliamono_cds <- ciliamono[rowData(ciliamono)$gene_short_name %in% genes, ]
png('ciliamono_line_1.png', width = 1000, height = 1067)
plot_genes_in_pseudotime(ciliamono_cds,
                        color_cells_by="ciliated2",
                        panel_order=genes,
                        min_expr=0.5,
                        cell_size=0,
                        ncol=3) +
                        theme(plot.title = element_text(size=25,color = 'black',hjust = 0.5),
                        axis.text.x=element_text(size=18, color = 'black',angle=0),
                        axis.text.y = element_text(size = 18, color = 'black'),
                        strip.text.x = element_text(size=25),
                        axis.line = element_line(color = "black"),
                        axis.title.x = element_text(size=20),
                        axis.title.y = element_text(size=20)
                        ) + NoLegend()
dev.off()

ciliamono_cds <- ciliamono[rowData(ciliamono)$gene_short_name %in% c('scv2-orf1-10'), ]
png('ciliamono_line_scv2.png', width = 333, height = 333)
plot_genes_in_pseudotime(ciliamono_cds,
                         color_cells_by="ciliated2",
                         min_expr=0.5,
                         cell_size=0) +
                         theme(plot.title = element_text(size=25,color = 'black',hjust = 0.5),
                         axis.text.x=element_text(size=18, color = 'black',angle=0),
                         axis.text.y = element_text(size = 18, color = 'black'),
                         strip.text.x = element_text(size=25),
                         axis.line = element_line(color = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20)
                         ) + NoLegend()
dev.off()

traj.coord <- ciliamono@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
ciliated <- AddMetaData(ciliated, metadata = traj.coord, col.name = 'pseudotime')

png('ciliated_fp_pseudotime.png', width = 600, height = 800)
FeaturePlot(ciliated,'pseudotime',cols=viridis(60),label=T,label.size=8) + labs(title = 'Pseudotime') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()
png('ciliated_fp_pseudotime2.png', width = 1200, height = 400)
FeaturePlot(ciliated,'pseudotime',cols=viridis(60),label=F,split.by='Condition') +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()

save(ciliamono, file = paste("ciliamono",Sys.Date(),".Robj",sep=""))
load('./ciliamono2020-08-10.Robj')


### How does Notch vary over differentiation trajectory?
notch.genes <- c('ADAM17','APH1A','CIR1','CREBBP','CTBP1','CTBP2','DLL1','DLL3','DLL4','DTX1','DTX2','DTX3',
                'DTX3L','DTX4','DVL1','DVL2','DVL3','EP300','HDAC1','HDAC2','HES1','HES5','JAG1','JAG2',
                'KAT2A','KAT2B','LFNG','MAML1','MAML2','MAML3','MFNG','NCOR2','NCSTN','NOTCH1','NOTCH2',
                'NOTCH3','NOTCH4','NUMB','NUMBL','PSEN1','PSEN2','PSENEN','PTCRA','RBPJ','RBPJL','RFNG','SNW1')
notch.genes <- intersect(notch.genes,rownames(ciliated))
notch.use <- HBECStatSig(ciliated,notch.genes,min.pct = 0.1,logfc.threshold = 0)
sig.notch.genes <- sort(unique(notch.use$gene))

ciliated_cds <- ciliamono[rowData(ciliamono)$gene_short_name %in% sig.notch.genes, ]
ht <- ceiling(length(sig.notch.genes)/6)*222.25
png('ciliated_line_notch.png', width = 1333, height = ht)  #1067
plot_genes_in_pseudotime(ciliated_cds,
                         color_cells_by="ciliated2",
                         panel_order=sig.notch.genes,
                         min_expr=0.5,
                         cell_size=0,
                         ncol=6) +
                         theme(plot.title = element_text(size=25,color = 'black',hjust = 0.5),
                         axis.text.x=element_text(size=18, color = 'black',angle=0),
                         axis.text.y = element_text(size = 18, color = 'black'),
                         strip.text.x = element_text(size=25),
                         axis.line = element_line(color = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20)
                         ) + NoLegend()
dev.off()





######
load('./ciliamono2020-08-10.Robj')
load('./hbecmono2020-08-10.Robj')

# Ciliated differentiation markers
gene <- c('MCIDAS','FOXJ1','RFX2')
ciliamono_cds <- ciliamono[rowData(ciliamono)$gene_short_name %in% gene, ]
png('ciliamono_line_cil.png', width = 667, height = 333)
plot_genes_in_pseudotime(ciliamono_cds,
                         color_cells_by="ciliated2",
                         min_expr=0.5,
                         cell_size=0,
                         panel_order=gene,
                         ncol=3) +
                         theme(plot.title = element_text(size=25,color = 'black',hjust = 0.5),
                         axis.text.x=element_text(size=18, color = 'black',angle=0),
                         axis.text.y = element_text(size = 18, color = 'black'),
                         strip.text.x = element_text(size=25),
                         axis.line = element_line(color = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20)
                         ) + NoLegend()
dev.off()






### Probing genes of interest
# General Prep
Idents(hbec.integrated) <- hbec.integrated[['Condition2']]
Idents(hbec.integrated) <- factor(x=Idents(hbec.integrated), levels = rev(c('M','1','2','3')))
group_colors <- colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Blues'))(2)

# all-ish
genes <- c('WNT5A','WNT4','TGFB3','EGR1','TGFBR2','BMP3','BMP4','BMPR1A','BMPR1B','BMPR2','RGMB',
            'FGF14','IGFBP1','IGFBP2','IGFBP3','IGFBP4','IGFBP5','IGFBP6','IGFBP7','IGF','AREG',
            'EGF','VEGFA','RBB2','FGF1','JAG1','JAG2','YAP1')
# WNT of interest
genes <- c('WNT5A','WNT4')
# TGFbeta of interest
genes <- c('TGFB3','TGFB1')
# Metabolism switch
genes <- c('TGFBR2','EGR1','HK2', 'LDHA', 'PGK1','PRDX5', 'GSTA1','ATP5ME','NFE2L2')
# BMP of interest
genes <- c('BMP3','BMP4','ID2','ID3','BMPR1A','BMPR1B','BMPR2','RGMB')
# GF of interest
matches <- sort(unique(grep(paste(c('IGFBP','IGF'),collapse="|"),rownames(hbec.integrated), value=TRUE)),decreasing = T)
genes <- c('FGF14','IGFBP2','IGFBP3','IGFBP4','IGFBP5','IGFBP6','IGFBP7')
# EGF of interest
genes <- c('AREG','EREG','VEGFA','ERBB2','HBEGF')
# Notch of interest
genes <- c('JAG1','JAG2','HES1')
# Hippo of interest
genes <- c('YAP1','JAG1','JAG2','HES1')
# IFIT of interest
ifit.genes <- sort(unique(grep(paste(c('IFIT'),collapse="|"),rownames(hbec.integrated), value=TRUE)),decreasing=T)
genes <- c('IFIT1','IFIT2','IFIT3','IFIT5')

dot <- DotPlot(hbec.integrated,features = genes,assay='RNA',
            cols = group_colors,scale = T,dot.scale = 30) +
            theme(plot.title = element_text(size=50,color = 'black',hjust = 0.5),
                    axis.text.x=element_text(size=32, color = 'black',angle=90),
                    axis.text.y = element_text(size = 42, color = 'black'),
                    legend.text = element_text(size = 20),
                    legend.key.size = unit(15, "mm"),
                    legend.title = element_text(size=20),
                    strip.text.x = element_text(size=44),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(color = "black"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank())

png('hbec.integrated_dot_bmp.png', width = (200*length(genes)), height = 700)
dot
dev.off()

Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
gene <- c('GSTP1')
png('hbec.integrated_fp_GSTP1.png', width = 650, height = 800)
FeaturePlot(hbec.integrated,gene,label = T,label.size=8,repel=T) + labs(title = gene) + NoAxes() +
            theme(plot.title = element_text(size = 38, hjust = 0.5),
            legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()

gene <- c('ifn.score1')
png('hbec.integrated_vln_ifn2.png', width = 800,height=700)
VlnPlot(hbec.integrated,gene,pt.size = 0)
dev.off()

Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
Idents(hbec.integrated) <- factor(x=Idents(hbec.integrated), levels = c('CyclingBasal','Basal','BasalEMT','Intermediate','Club','Goblet','Ciliated','Tuft','Ionocytes','PNECs'))
p1 <- VlnPlot(hbec.integrated,'GSTP1',pt.size = 0, cols = celltypecolors) + labs(title = 'GSTP1') + NoLegend() + xlab('') +
            theme(plot.title = element_text(size = 32, hjust = 0.5),
            axis.text.x=element_text(size=20, color = 'black'),
            axis.text.y = element_text(size = 28, color = 'black'))
png('hbec.integrated_vln_GSTP1.png', width = 800,height=500)
p1
dev.off()


## Scoring of interest
glyc.genes <- c("GPI","PKM","TPI1","PGAM1","PGK1","ALDOC","ENO1","ALDOA","ENO2","GAPDH","HK2","PFKP")
hbec.integrated <- AddModuleScore(hbec.integrated,features=list(glyc.genes),ctrl = 5,name='glyc.score')
ox.genes <- c("MT-ND4","MT-ND5","NDUFB3","NDUFB2","UQCR11","COX6A1","MT-ND2","MT-ND1","MT-ATP6","NDUFV2",
        "ATP5MF","SURF1","NDUFA4","NDUFA3","NDUFA2","NDUFC2","NDUFC1","COX6C","ATP5F1E","UQCRQ","NDUFS6",
        "NDUFS5","NDUFAB1","MT-CO3","MT-CYB","MT-CO1","MT-CO2","MT-CYB")
hbec.integrated <- AddModuleScore(hbec.integrated,features=list(ox.genes),ctrl = 5,name='ox.score')
ifn.genes <- unique(c("IFITM3","EGR1","IFITM2","RSAD2","STAT1","MX2","STAT2","MX1","HLA-B","IFI6","ISG15",
        "ADAR","IFIT1","IFIT3","IFIT2","BST2","IFI27","OAS1","OAS2","OAS3","XAF1","IFIH1","IFI16","DDX58",
        "PRKDC","ZC3HAV1","HMGB1","NFKB1","HSPD1","PARP9","PARP14","BST2",'IFIT5'))
hbec.integrated <- AddModuleScore(hbec.integrated,features=list(ifn.genes),ctrl = 5,name='ifn.score')

Idents(hbec.integrated) <- hbec.integrated[['Condition']]
Idents(hbec.integrated) <- factor(x=Idents(hbec.integrated),levels=c('Mock','1dpi','2dpi','3dpi'))
png('hbec.integrated_fp_ifn.png', width = 1600, height = 400)
FeaturePlot(hbec.integrated,'ifn.score1',split.by='Condition',cols = rev(brewer.pal(4,'RdBu'))) + theme(legend.position = 'right')
dev.off()

png('hbec.integrated_fp_ifn3.png', width = 650, height = 800)
FeaturePlot(hbec.integrated,'ifn.score1',label=F,cols = rev(brewer.pal(4,'RdBu'))) + labs(title = 'Ifn Score - Bulk') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()

p1 <- VlnPlot(hbec.integrated,'ifn.score1',pt.size = 0,group.by='Condition2') + xlab('') + labs(title = 'Interferon Signaling') +
            theme(plot.title = element_text(size = 32, hjust = 0.5),
            axis.text.x=element_text(size=28, color = 'black',angle = 0,hjust = 0.5),
            axis.text.y = element_text(size = 28, color = 'black')) + NoLegend()
o1 <- VlnPlot(hbec.integrated,'ox.score1',pt.size = 0,group.by='Condition2') + xlab('') + labs(title = 'Oxidative Phosphorylation') +
            theme(plot.title = element_text(size = 32 , hjust = 0.5),
                axis.text.x=element_text(size=28, color = 'black',angle = 0,hjust = 0.5),
                axis.text.y = element_text(size = 28, color = 'black')) + NoLegend()

p1 <- plot_grid(g1, o1,ncol=2)

png('hbec.integrated_vln_ifn.png', width = 450,height=500)
p1
dev.off()

Idents(ciliated) <- ciliated[['ciliated2']]
ciliated <- AddModuleScore(ciliated,features=list(ifn.genes),ctrl = 5,name='ifn.score')
gene <- c('ifn.score1')
png('ciliated_vln_ifn2.png', width = 800,height=700)
VlnPlot(ciliated,gene,pt.size = 0)
dev.off()


ciliated <- RenameIdents(ciliated, 'Mature_Ciliated1' = 'MC1',
                                    'Ciliated_Progenitor' = 'CP',
                                    'Mature_Ciliated2' = 'MC2',
                                    'Novel_Infected_Ciliated' = 'NIC')
ciliated[['ciliated4']] <- Idents(ciliated)

Idents(ciliated) <- ciliated[['ciliated4']]
Idents(ciliated) <- factor(x=Idents(ciliated), levels = c('MC1','CP','MC2','NIC'))
p1 <- VlnPlot(ciliated,'ifn.score1',pt.size = 0, cols = cilcolors) + labs(title = 'Ifn Score') + NoLegend() + xlab('') +
            theme(plot.title = element_text(size = 32, hjust = 0.5),
            axis.text.x=element_text(size=20, color = 'black'),
            axis.text.y = element_text(size = 28, color = 'black'))
png('ciliated_vln_ifn3.png', width = 400,height=500)
p1
dev.off()


Idents(ciliated) <- ciliated[['Condition']]
Idents(ciliated) <- factor(x=Idents(ciliated),levels=c('Mock','1dpi','2dpi','3dpi'))
png('ciliated_fp_ifn.png', width = 1600, height = 400)
FeaturePlot(ciliated,'ifn.score1',split.by='Condition',cols = rev(brewer.pal(4,'RdBu'))) + theme(legend.position = 'right')
dev.off()

png('ciliated_fp_ifn2.png', width = 650, height = 800)
FeaturePlot(ciliated,'ifn.score1',label=F,pt.size=1,cols = rev(brewer.pal(4,'RdBu'))) + labs(title = 'Ifn Score - Ciliated') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()

png('ciliated_fp_IFIT1.png', width = 650, height = 800)
FeaturePlot(ciliated,'IFIT1',label=F,pt.size=1) + labs(title = 'IFIT1') + theme(plot.title = element_text(size = 38, hjust = 0.5)) + NoAxes() +
            theme(legend.text = element_text(size = 22),
            legend.key.size = unit(1.5, "lines"))
dev.off()


#### Pseudotime heatmap of cells
# subset the object
Idents(hbec.integrated) <- hbec.integrated[['cell_type2']]
hbec_short <- subset(hbec.integrated,idents=c('Basal','Intermediate','Club','Ciliated'))

genes <- c('S100A2','FABP5','KRT5','KRT6A','KRT14','IGFBP6',
        'SCGB1A1','SCGB3A1','CEACAM6','CEACAM5',
        'BPIFB1','LYPD2','NFKBIA','scv2-orf1-10',
        'TUBB4B','TUBA1A','PIFO','TPPP3','CCDC78')
genes <- ifn.genes

hbec_short <- ScaleData(hbec_short,features = rownames(hbec_short))
hbec_store <- hbec_short

hbec_short <- hbec_store
hbec_short.avg <- GetAssayData(hbec_short, slot = "scale.data")[genes, ]
hbec_short.avg <- t(scale(t(hbec_short.avg)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
hbec_short.avg <- t(apply(as.matrix(hbec_short.avg), 1, unityNormalize))

hbec_short.avg <- hbec_short.avg[,order(hbec_short$pseudotime)]
dim(hbec_short.avg)

hbec_short.pseudo <- as.matrix(hbec_short$pseudotime)
hbec_short.pseudo <- as.matrix(hbec_short.pseudo[order(hbec_short$pseudotime),])
hbec_short.celltype2 <- as.matrix(hbec_short$cell_type2)
hbec_short.celltype2 <- as.matrix(hbec_short.celltype2[order(hbec_short$pseudotime),])
hbec_short.inf <- as.matrix(hbec_short$'scv2+')
hbec_short.inf <- as.matrix(hbec_short.inf[order(hbec_short$pseudotime),])

n = ncol(hbec_short.avg)  # running with this rather than 60 takes ages, good for final, not for iteration
colors.inferno <- colorRamp2(breaks = c(seq(min(hbec_short.avg),max(hbec_short.avg),length.out=60)), inferno(n=60), space = "RGB")
n = nrow(hbec_short.pseudo)
colors.viridis <- colorRamp2(breaks = c(seq(min(hbec_short.pseudo),max(hbec_short.pseudo),length.out=n)), viridis(n=n), space = "RGB")

ha = HeatmapAnnotation(df = list(data.frame(hbec_short.pseudo),
                                data.frame(hbec_short.celltype2),
                                data.frame(hbec_short.inf)),
                        col = list(hbec_short.pseudo=colors.viridis,
                                hbec_short.celltype2=c('Basal' = "#F37912","Intermediate" = "#658E67","Club" = "#A35390",'Ciliated' = '#E41A1C'),
                                hbec_short.inf=c('1'='blue','0'='white')))

hbec_short.heatmap <- Heatmap(hbec_short.avg,
                        col= colors.inferno,
                        top_annotation = ha,
                        #column_split= factor(colnames(hbec_short.avg), levels = identities),
                        #column_title = identities,
                        #column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 14),
                        cluster_rows=T,
                        #row_km = 3,
                        cluster_columns=F,
                        show_column_names=F,
                        show_row_names = T,
                        show_heatmap_legend = F)
png('hbec_short_hm_ifn.png', width = 1600, height = 400)  # was 1000
draw(hbec_short.heatmap,padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


#### Pseudotime heatmap of CILIATED cells
genes <- entry.genes

ciliated <- ScaleData(ciliated,features = rownames(ciliated))

ciliated.avg <- GetAssayData(ciliated, slot = "scale.data")[genes, ]
ciliated.avg <- t(scale(t(ciliated.avg)))
unityNormalize <- function(x){(x-min(x))/(max(x)-min(x))}
ciliated.avg <- t(apply(as.matrix(ciliated.avg), 1, unityNormalize))

ciliated.avg <- ciliated.avg[,order(ciliated$pseudotime)]
dim(ciliated.avg)
#[1]    18 21510

ciliated.pseudo <- as.matrix(ciliated$pseudotime)
ciliated.pseudo <- as.matrix(ciliated.pseudo[order(ciliated$pseudotime),])
ciliated.ciliated3 <- as.matrix(ciliated$ciliated3)
ciliated.ciliated3 <- as.matrix(ciliated.ciliated3[order(ciliated$pseudotime),])
ciliated.inf <- as.matrix(ciliated$'scv2+')
ciliated.inf <- as.matrix(ciliated.inf[order(ciliated$pseudotime),])

n = ncol(ciliated.avg)  # running with this rather than 60 takes ages, good for final, not for iteration
colors.inferno <- colorRamp2(breaks = c(seq(min(ciliated.avg),max(ciliated.avg),length.out=60)), inferno(n=60), space = "RGB")
n = nrow(ciliated.pseudo)
colors.viridis <- colorRamp2(breaks = c(seq(min(ciliated.pseudo),max(ciliated.pseudo),length.out=n)), viridis(n=n), space = "RGB")

ha = HeatmapAnnotation(df = list(data.frame(ciliated.pseudo),
                                data.frame(ciliated.ciliated3),
                                data.frame(ciliated.inf)),
                        col = list(ciliated.pseudo=colors.viridis,
                                ciliated.ciliated3=c('Mature_Ciliated1' = "#E78AC3","Ciliated_Progenitor" = "#66C2A5","Mature_Ciliated2" = "#8DA0CB",'Novel_Infected_Ciliated' = '#FC8D62'),
                                ciliated.inf=c('1'='blue','0'='white')))

ciliated.heatmap <- Heatmap(ciliated.avg,
                        col= colors.inferno,
                        top_annotation = ha,
                        #column_split= factor(colnames(ciliated.avg), levels = identities),
                        #column_title = identities,
                        #column_title_gp = gpar(fontsize = 25),
                        row_names_gp = gpar(fontsize = 14),
                        cluster_rows=T,
                        #row_km = 3,
                        cluster_columns=F,
                        show_column_names=F,
                        show_row_names = T,
                        show_heatmap_legend = F)
png('ciliated_cellhm_entry2.png', width = 1200, height = 400)  # was 1000
draw(ciliated.heatmap,padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()





#### Measuring genes in infected cells
Idents(hbec.integrated) <- hbec.integrated[['Condition']]
png('hbec.integrated_vln_ifnb1.png', width = 400,height=400)
VlnPlot(hbec.integrated,c('IFNB1'),pt.size = 0.2,split.plot=T,split.by='scv2+') + xlab('') + ylab('') +
    theme(axis.text.x=element_text(size=32, color = 'black',angle = 0,hjust = 0.5),
    axis.text.y = element_text(size = 28, color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 38, hjust = 0.5))
dev.off()

png('hbec.integrated_vln_il6.png', width = 400,height=400)
VlnPlot(hbec.integrated,c('IL6'),pt.size = 0.2,split.plot=T,split.by='scv2+') + xlab('') + ylab('') +
    theme(axis.text.x=element_text(size=32, color = 'black',angle = 0,hjust = 0.5),
    axis.text.y = element_text(size = 28, color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 38, hjust = 0.5))
dev.off()




