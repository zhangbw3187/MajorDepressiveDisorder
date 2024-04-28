
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(tidyr)

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)



data <- as(as.matrix(cd4@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = cd4@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())


monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)


######
gob<- FindVariableFeatures(cd4)
exp_genes <- VariableFeatures(cd4)
monocle_cds <- setOrderingFilter(monocle_cds, exp_genes)
plot_ordering_genes(monocle_cds)

#########

disp_table <- dispersionTable(monocle_cds)
disp.gene <- subset(disp_table,mean_expression >= 0.15 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
monocle_cds <- setOrderingFilter(monocle_cds,disp.gene)
######

monocle_cds <- reduceDimension(monocle_cds,max_components =2, method="DDRTree")
monocle_cds <- orderCells(monocle_cds)

plot_cell_trajectory(monocle_cds,color_by="seurat_clusters",size=1, show_backbone = T)
plot_cell_trajectory(monocle_cds,color_by="State",size=1, show_backbone = T)


monocle_cds <- orderCells(monocle_cds, root_state=7)

cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a",
		"#66a61e","#e6ab02","#a6761d","#666666")

plot_cell_trajectory(monocle_cds,color_by="seurat_clusters",size=1, show_backbone = T)+scale_color_manual(values=cols)


s.gene <- c("SELL","CCR7","IL7R","CD84","CCL5","S100A4","CD8A","CD3E")
plot_genes_violin(monocle_cds[s.gene,],grouping = "State", color_by="State")


show.gene <- c("PTPRC","IL7R","GPR183","PRKDC","MALT1","CD2","CAMK4","CD46","ANXA1")

show.gene <- c("CCR7","CDC42","CXCR4","CD69")
plot_genes_violin(monocle_cds[show.gene,],grouping = "State", color_by="State")

show.gene <- c("TCF7","IL7R")
plot_genes_violin(monocle_cds[show.gene,],grouping = "cond.res", color_by="cond.res")


show.gene <- c("PTPRC","IL7R","GPR183","PRKDC","MALT1","CD2","CAMK4","CD46","ANXA1")
show.gene <- c("PTPRC","IL7R","GPR183","PRKDC","MALT1","CD4","CD2","CAMK4","CD3G","CD46","ANXA1")
plot_genes_violin(monocle_cds[show.gene,],grouping = "cond.res", color_by="cond.res")
plot_genes_violin(monocle_cds[show.gene,],grouping = "seurat_clusters", color_by="seurat_clusters")

colnames(pData(monocle_cds))

pData(monocle_cds)$PTPRC = log2(exprs(monocle_cds)["PTPRC",] + 1)
pData(monocle_cds)$PRKDC = log2(exprs(monocle_cds)["PRKDC",] + 1)
pData(monocle_cds)$CD2 = log2(exprs(monocle_cds)["CD2",] + 1)
pData(monocle_cds)$IL7R = log2(exprs(monocle_cds)["IL7R",] + 1)

pData(monocle_cds)$CD8A = log2(exprs(monocle_cds)["CD8A",] + 1)
pData(monocle_cds)$CD8B = log2(exprs(monocle_cds)["CD8B",] + 1)

#plot_cell_trajectory(monocle_cds,color_by="PTPRC")

p1 <- plot_cell_trajectory(monocle_cds,color_by="PTPRC")+scale_color_gradient(low="purple",high="yellow")
p2 <- plot_cell_trajectory(monocle_cds,color_by="PRKDC")+scale_color_gradient(low="purple",high="yellow")
p3 <- plot_cell_trajectory(monocle_cds,color_by="CD2")+scale_color_gradient(low="purple",high="yellow")
p4 <- plot_cell_trajectory(monocle_cds,color_by="IL7R")+scale_color_gradient(low="purple",high="yellow")
 plot_cell_trajectory(monocle_cds,color_by="CD8A")+scale_color_gradient(low="purple",high="yellow")
 
 
################
mks.sig2 <- mks.sig[mks.sig$avg_log2FC > 0,]DE <- mks.sig2[,c(7,6)]names(DE) <- c("SYMBOL","Label")table(DE$Label)de <- split(DE, DE$Label)list.de <- NULLfor (i in 1:length(de)){     Label <- unique(de[[i]]$Label)     de[[i]] <- bitr(de[[i]]$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")     de[[i]] <- cbind(de[[i]], Label)     de[[i]]$Label <- as.character(de[[i]]$Label)     list.de <- rbind(list.de, de[[i]])}

uchc_immu_gobp <- compareCluster(ENTREZID~Label, data= list.de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP")dotplot(uchc_immu_gobp,showCategory = 20)





                        
nt1_kk <- enrichGO(gene = de[[2]]$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.1,
               readable      = TRUE)


nt3_kk <- enrichGO(gene = de[[4]]$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.1,
               readable      = TRUE)
               
               
####               
diff_test_res <- differentialGeneTest(monocle_cds,
                    fullModelFormulaStr = "~State")
diff_test_res[,c("gene_short_name", "pval", "qval")]           
               
#######

Time_diff <- differentialGeneTest(monocle_cds[exp_genes,],cores=1,fullModelFormulaStr="~sm.ns(Pseudotime)")

Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()


Time_diff_sig <- Time_diff[Time_diff$pval < 0.05,]

Time_genes <- Time_diff_sig %>% pull(gene_short_name) %>% as.character()


######
count <- table(monocle_cds$State, monocle_cds$cond.res)
count <- count[,2:7]
count <-count/apply(count,1,sum)


####
BEAM_res <- BEAM(monocle_cds[exp_genes,],branch_point = 2,cores=1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res2 <- BEAM_res[,c("gene_short_name","pval","qval")]
BEAM_res3 <- BEAM_res2[!grepl("^RP[SL]",BEAM_res2$gene_short_name),]
BEAM_res3 <- BEAM_res3[!grepl("^MT-",BEAM_res3$gene_short_name),]



plot_genes_branched_heatmap(monocle_cds[row.names(subset(BEAM_res2,qval<1e-4)),],branch_point=2,num_clusters=4,cores=1,use_gene_short_name=T,show_rownames=T)
plot_genes_branched_heatmap(monocle_cds[row.names(BEAM_res3[1:100,]),],branch_point=2,num_clusters=4,cores=1,use_gene_short_name=T,show_rownames=T)



###### send back to seurat


monocle_meta <- data.frame("cell_ID" = row.names(monocle_cds@phenoData@data), "State"=monocle_cds@phenoData@data$State,"Pseudotime"=monocle_cds@phenoData@data$Pseudotime)

row.names(monocle_meta)<- monocle_meta$cell_ID
cd4 <- AddMetaData(cd4,monocle_meta)
Idents(cd4)<- cd4$State
state.mks <- FindAllMarkers(cd4, only.pos=T)
state.mks <- state.mks[state.mks$p_val_adj<0.05,]


##########
################
#state.mks2 <- state.mks[state.mks$avg_log2FC > 0,]

state.mks2 <- state.mks[!grepl("^RP[SL]",state.mks$gene),]
state.mks2 <- state.mks2[!grepl("^MT-",state.mks2$gene),]



top10 <- state.mks2 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

DotPlot(cd4,features = unique(top10$gene),col.min= -1.5, col.max= 1.5)+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("cd4 clusters")+ylab("")


DE <- state.mks2[,c(7,6)]
names(DE) <- c("SYMBOL","Label")
table(DE$Label)

de <- split(DE, DE$Label)
list.de <- NULL
for (i in 1:length(de)){
     Label <- unique(de[[i]]$Label)
     de[[i]] <- bitr(de[[i]]$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
     de[[i]] <- cbind(de[[i]], Label)
     de[[i]]$Label <- as.character(de[[i]]$Label)
     list.de <- rbind(list.de, de[[i]])
}

gobp <- compareCluster(ENTREZID~Label, data= list.de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP")
dotplot(gobp,showCategory = 10)


gomf <- compareCluster(ENTREZID~Label, data= list.de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="MF")
dotplot(gomf,showCategory = 10)


kegg <- compareCluster(ENTREZID~Label, data= list.de, fun="enrichKEGG", organism = "hsa")
dotplot(kegg,showCategory = 10)
dotplot(gobp,showCategory = 20)



                        
nt1_kk <- enrichGO(gene = de[[2]]$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.1,
               readable      = TRUE)


nt3_kk <- enrichGO(gene = de[[4]]$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.1,
               readable      = TRUE)
               
##########
######## volcano plot ####
library(ggrepel)
de.genes <- as.data.frame(de.naive)
#de.genes$gene <- sub("__.*","",rownames(de.genes))
de.genes$sig <- "nonsig"
de.genes$set <- "all"
significance <- function(x){
  x[x$p_val_adj < 0.05,  "sig"] <- "Padj < 0.05"
  return(x)
}
de.genes <- significance(de.genes)
de.genes$set <- "no"
de.genes$set[de.genes$sig == "Padj < 0.05" & abs(de.genes$logFC) > 0.25] <- "yes"

de.genes$color="grey"
de.genes$color[de.genes$sig == "Padj < 0.05"]="black"
de.genes$color[de.genes$sig == "Padj < 0.05" & de.genes$logFC > 0.25] <- "red"
de.genes$color[de.genes$sig == "Padj < 0.05" & de.genes$logFC < -0.25] <- "blue"

#de.genes[abs(de.genes$log2FoldChange) > 1.5, "set"] <- "yes"
ggplot(de.genes, aes(x = logFC, y = -log(p_val_adj, 10), label = SYMBOL)) +
  geom_point(aes(color = color), size = 1) +
  geom_label_repel(aes(label = ifelse(set  == "yes", as.character(SYMBOL), '')),   size = 4) +
  #geom_label_repel(aes(label = ifelse(set  == "yes", as.character(SYMBOL), '')), hjust = 1.25, vjust = 0,   size = 4) +
  scale_colour_manual(name = "Label", 
                      values = c("black" = "black", "grey" = "gray", "red"="red","blue"="blue"),
                      breaks = c("black", "grey", "red","blue"),
                      labels = c("p.adjusted < 0.05", "not significant", "up-regulated in MDD", "down-regulated in MDD")) +
  ggtitle(label = "DEGs of baseline MDD vs HC in Naive T cells") +
  ylab("p-value (-log10)") +
  xlab("log2 fold change") +
  geom_hline(yintercept = -log(0.05,10)) +
  geom_vline(xintercept = 0) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 20))+theme_classic()
