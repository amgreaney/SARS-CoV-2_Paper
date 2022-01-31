# object is a Seurat 3.0 object, function will operate on the active Identity slot
# this is formatted specifically for the HBEC object, which has a 'Condition' metadata slot
# GOI should be a list of genes of interest
# ... arguments to pass to FindMarkers

HBECStatSig <- function(object,GOI,...){
  Idents(object) <- object[['cell_type2']]
  # Stash idents and make new identities which identify each by their 'Condition' metadata
  celltypes <- as.character(unique(Idents(object)))
  celltypes.mock <- paste(celltypes, 'Mock', sep = '_')
  celltypes.1dpi <- paste(celltypes, '1dpi', sep = '_')
  celltypes.2dpi <- paste(celltypes, '2dpi', sep = '_')
  celltypes.3dpi <- paste(celltypes, '3dpi', sep = '_')

  object$celltype.condition <- paste(Idents(object), object$Condition, sep = "_")
  object$celltype <- Idents(object)
  Idents(object) <- "celltype.condition"

  #### Day 1 ####
  # Identify which genes, in each cell population, have an adjusted p-value < 0.05 based on a Wilcoxon rank test when compared to Mock
  genes <- GOI[GOI %in% rownames(object)]
  diff.p.dpi1 <- data.frame()
  for (i in 1:length(celltypes)){
    temp <- NULL
    try(
      temp <- FindMarkers(object,
                        ident.1 = celltypes.1dpi[i],
                        ident.2 = celltypes.mock[i],
                        verbose = FALSE,
                        features = genes,
                        min.cells.group = 0,...)
        )
    if(!is.null(temp)){
    temp2 <- subset(temp, p_val_adj < 0.05)
    if (nrow(temp2)>0){
      temp2$gene <- rownames(temp2)
      temp2$cell <- celltypes[i]
      diff.p.dpi1 <- rbind(diff.p.dpi1, temp2)
    }
    }
  }
  diff.p.dpi1$cell.gene <- paste(diff.p.dpi1$cell,diff.p.dpi1$gene,sep = '.')

  #### Day 2 ####
  diff.p.dpi2 <- data.frame()
  for (i in 1:length(celltypes)){
    temp <- NULL
    try(
      temp <- FindMarkers(object,
                        ident.1 = celltypes.2dpi[i],
                        ident.2 = celltypes.mock[i],
                        verbose = FALSE,
                        features = genes,
                        min.cells.group = 0,...)
       )
    if(!is.null(temp)){
    temp2 <- subset(temp, p_val_adj < 0.05)
    if (nrow(temp2)>0){
      temp2$gene <- rownames(temp2)
      temp2$cell <- celltypes[i]
      diff.p.dpi2 <- rbind(diff.p.dpi2, temp2)
    }
    }
  }
  diff.p.dpi2$cell.gene <- paste(diff.p.dpi2$cell,diff.p.dpi2$gene,sep = '.')

  #### Day 3 ####
  diff.p.dpi3 <- data.frame()
  for (i in 1:length(celltypes)){
    temp <- NULL
    try(
     temp <- FindMarkers(object,
                        ident.1 = celltypes.3dpi[i],
                        ident.2 = celltypes.mock[i],
                        verbose = FALSE,
                        features = genes,
                        min.cells.group = 0,...)
    )
    if(!is.null(temp)){
    temp2 <- subset(temp, p_val_adj < 0.05)
    if (nrow(temp2)>0){
      temp2$gene <- rownames(temp2)
      temp2$cell <- celltypes[i]
      diff.p.dpi3 <- rbind(diff.p.dpi3, temp2)
    }
    }
  }
  diff.p.dpi3$cell.gene <- paste(diff.p.dpi3$cell,diff.p.dpi3$gene,sep = '.')

  # Tag with condition
  diff.p.dpi1$Condition <- 'Dpi1.vs.Mock'
  diff.p.dpi2$Condition <- 'Dpi2.vs.Mock'
  diff.p.dpi3$Condition <- 'Dpi3.vs.Mock'

  # Bind into single output
  differential.genes <- rbind(diff.p.dpi1,diff.p.dpi2,diff.p.dpi3)

  # Erase rownames which are nonsense
  rownames(differential.genes) <- NULL
  return(differential.genes)
}
