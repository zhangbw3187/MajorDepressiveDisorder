
library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)


#######################

ReadBowen <- function(bid,individual,btype){
	bdata <- Read10X(data.dir = bid)
	bsample <- CreateSeuratObject(raw.data = bdata, min.cells = 3, min.genes = 200, project = "10X")
	bsample@meta.data$sampid <- bid
	bsample@meta.data$individual <- individual
	bsample@meta.data$celtype <- btype
	return(bsample)
}

idN6 <- ReadBowen("N6","idN6","N")
idN5 <- ReadBowen("N5","idN5","N")
idN4 <- ReadBowen("N4","idN4","N")
idN3 <- ReadBowen("N3","idN3","N")
idN2 <- ReadBowen("N2","idN2","N")
idN7 <- ReadBowen("N7","idN7","N")
idN8 <- ReadBowen("N8","idN8","N")

idM2.1 <- ReadBowen("M2-1","idM2","M")
idM3.1 <- ReadBowen("M3-1","idM3","M")
idM4.1 <- ReadBowen("M4-1","idM4","M")
idM5.1 <- ReadBowen("M5-1","idM5","M")
idM6.1 <- ReadBowen("M6-1","idM6","M")
idM7.1 <- ReadBowen("M7-1","idM7","M")
idM8.1 <- ReadBowen("M8-1","idM8","M")

idM2.2 <- ReadBowen("M2-2","idM2","M")
idM2.3 <- ReadBowen("M2-3","idM2","M")
idM3.2 <- ReadBowen("M3-2","idM3","M")
idM3.3 <- ReadBowen("M3-3","idM3","M")
idM4.2 <- ReadBowen("M4-2","idM4","M")
idM4.3 <- ReadBowen("M4-3","idM4","M")
idM5.2 <- ReadBowen("M5-2","idM5","M")
idM5.3 <- ReadBowen("M5-3","idM5","M")
idM6.2 <- ReadBowen("M6-2","idM6","M")
idM6.3 <- ReadBowen("M6-3","idM6","M")
idM7.2 <- ReadBowen("M7-2","idM7","M")
idM7.3 <- ReadBowen("M7-3","idM7","M")
idM8.2 <- ReadBowen("M8-2","idM8","M")
idM8.3 <- ReadBowen("M8-3","idM8","M")


###################### 
group = MergeSeurat(obj1, obj2, project = "10X", add.cell.id1 = as.character("M1-1"), add.cell.id2 = as.character("M1-2"))



################

mito.genes <- grep(pattern = "^MT-", x = rownames(x = groupall@data), value = TRUE)
percent.mito <- Matrix::colSums(groupall@raw.data[mito.genes, ])/Matrix::colSums(groupall@raw.data)
groupall <- AddMetaData(object = groupall, metadata = percent.mito, col.name = "percent.mito")
#VlnPlot(object = groupall, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

groupflt <- FilterCells(object = groupall, subset.names = c("nGene", "percent.mito"), low.thresholds = c(100, -Inf), high.thresholds = c(2500, 0.25))
groupflt <- NormalizeData(object = groupflt, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
groupflt <- FindVariableGenes(object = groupflt, mean.function = ExpMean, dispersion.function = LogVMR)
groupflt <- ScaleData(object = groupflt, vars.to.regress = c("nUMI", "percent.mito"))


groupflt <- RunPCA(object = groupflt, pc.genes = groupflt@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
groupflt <- JackStraw(object = groupflt, num.replicate = 100, display.progress = FALSE)

groupcluster <- FindClusters(object = groupflt, reduction.type = "pca", dims.use = 1:15, resolution = 0.6, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)
groupcluster <- RunTSNE(object = groupcluster, dims.use = 1:15, do.fast = TRUE)
TSNEPlot(object = groupcluster,do.label = T, do.return = T, pt.size = 0.5)






