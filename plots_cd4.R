
library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)


#######
anding <- subset(anding, celltype.short %in% c("NaiCD4T"))
anding <- FindVariableFeatures(anding, selection.method = "vst", assay="RNA", nfeatures = 2000)
vft <- VariableFeatures(anding)


DimPlot(anding,label=TRUE)


vft.flt <- vft[grep("^RPS|^RPL|^MT-", vft,invert = T)]
top20 <- head(vft.flt, 20)
top1k <- head(vft.flt, 1000)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(anding)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2

anding <- ScaleData(anding, features=top1k, vars.to.regress = "nCount_RNA")
anding <- RunPCA(anding, npcs=30, features=top1k)
ElbowPlot(anding)
anding <- RunUMAP(anding, dims = 1:15)
anding <- FindNeighbors(anding, dims = 1:15)
anding <- FindClusters(anding, resolution = 0.5)
DimPlot(anding,label=TRUE)
DimPlot(anding,label=TRUE,split.by = "disease")


table(anding$sample)
anding$condition <- paste(anding$stim, anding$time)
table(anding$condition)
DimPlot(anding,label=TRUE,split.by = "condition",ncol=2)


anding$cond.res <- paste(anding$time,anding$responser, sep="_")
anding$cond.res[anding$cond.res == "ctrl_ctrl"]="ctrl"
anding$cond.res <- factor(anding$cond.res , levels=c("t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response","ctrl"))
DimPlot(anding,group.by="cond.res")
DimPlot(anding,split.by="cond.res",ncol=3)

tmp<-data.frame(cell_type= anding@active.ident,sample= anding@meta.data$sample, donor= anding@meta.data$cond.res)
res<-matrix(NA, ncol=4, nrow=length(unique(tmp$cell_type))*length(unique(tmp$sample)))
res<-as.data.frame(res)
names(res)<-c("cell_type","freq","sample","Condition")
n=1;
for (d in unique(tmp$sample)){
  t1<-subset(tmp,sample==d)
  t2<-as.data.frame(table(t1$cell_type))
  res[n:(n+length(unique(tmp$cell_type))-1),1]<-as.character(t2$Var1)
  res[n:(n+length(unique(tmp$cell_type))-1),2]<-t2$Freq
  res[n:(n+length(unique(tmp$cell_type))-1),3]<-rep(as.character(d),length(unique(tmp$cell_type)))
  res[n:(n+length(unique(tmp$cell_type))-1),4]<-rep(as.character(unique(t1$donor)),length(unique(tmp$cell_type)))
  n=n+length(unique(tmp$cell_type));
}
x<-split(res,res$sample)
y<-lapply(x, function(w) {w=transform(w,f=freq/sum(freq))})
y2<-do.call(rbind,y)

y2$Condition <- factor(y2$Condition, levels=c("t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response","ctrl"))
y2$other <- 1-y2$f

p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_bw()+scale_fill_manual(values=c("#fff7bc","#fee391","#fec44f","#c994c7","#df65b0","#ce1256","grey"))+xlab("")+ylab("percent")+theme(axis.text.x=element_text(angle = 70,hjust = 1))

###############


tmp<-data.frame(cell_type= anding@active.ident,sample= anding@meta.data$sample, donor= anding@meta.data$time)
res<-matrix(NA, ncol=4, nrow=length(unique(tmp$cell_type))*length(unique(tmp$sample)))
res<-as.data.frame(res)
names(res)<-c("cell_type","freq","sample","Condition")
n=1;
for (d in unique(tmp$sample)){
  t1<-subset(tmp,sample==d)
  t2<-as.data.frame(table(t1$cell_type))
  res[n:(n+length(unique(tmp$cell_type))-1),1]<-as.character(t2$Var1)
  res[n:(n+length(unique(tmp$cell_type))-1),2]<-t2$Freq
  res[n:(n+length(unique(tmp$cell_type))-1),3]<-rep(as.character(d),length(unique(tmp$cell_type)))
  res[n:(n+length(unique(tmp$cell_type))-1),4]<-rep(as.character(unique(t1$donor)),length(unique(tmp$cell_type)))
  n=n+length(unique(tmp$cell_type));
}
x<-split(res,res$sample)
y<-lapply(x, function(w) {w=transform(w,f=freq/sum(freq))})
y2<-do.call(rbind,y)

#y2$Condition <- factor(y2$Condition, levels=c("t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response","ctrl"))
y2$other <- 1-y2$f

p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_bw()+scale_fill_manual(values=c("#fff7bc","#fee391","#fec44f","#c994c7","#df65b0","#ce1256","grey"))+xlab("")+ylab("percent")+theme(axis.text.x=element_text(angle = 70,hjust = 1))



FeaturePlot(anding, features = c("CD3D","SELL","CD8A","GZMA"),order = T)

Idents(anding)<- "seurat_clusters"
pbmc.mks <- FindAllMarkers(object = anding,
only.pos = TRUE,
min.pct = 0.2,
logfc.threshold = 0.2,
min.diff.pct = 0.1,
test.use = "wilcox"
)


tt.markers <- pbmc.mks[grep("^RPS|^RPL",pbmc.mks$gene,invert = T),]
top10 <- tt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DotPlot(anding,features = unique(top10$gene),group.by="new.id",cols="RdBu")+coord_flip()+xlab("")+ylab("")
tt.markers <- pbmc.mks[grep("^RPS|^RPL|^MT-",pbmc.mks$gene,invert = T),]
top10 <- tt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DotPlot(anding,features = unique(top10$gene),group.by="seurat_clusters",cols="RdBu")+coord_flip()+xlab("")+ylab("")
table(tt.markers$cluster)
tt.markers.sig <- tt.markers[tt.markers$p_val_adj<0.05,]
table(tt.markers.sig$cluster)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
tt.markers[1:4,]
DEall <- tt.markers[,c(6,7)]
names(DEall) <- c("Label", "SYMBOL")
de <- split(DEall, DEall$Label)
list.de <- NULL
for (i in 1:length(de)){
Label <- unique(de[[i]]$Label)
de[[i]] <- bitr(de[[i]]$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
de[[i]] <- cbind(de[[i]], Label)
de[[i]]$Label <- as.character(de[[i]]$Label)
list.de <- rbind(list.de, de[[i]])
}
kegg.mono1 <- compareCluster(ENTREZID~Label, data=list.de, fun="enrichKEGG")
dotplot(kegg.mono1)+theme(axis.text.x = element_text(angle = 50,hjust = 1))
gobp.mono1 <- compareCluster(ENTREZID~Label, data= list.de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP",readable=TRUE)
dotplot(gobp.mono1)+theme(axis.text.x = element_text(angle = 50,hjust = 1))
top10 <- tt.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
DotPlot(anding,features = unique(top10$gene),cols="RdBu")+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+xlab("")+ylab("")
DEall <- tt.markers[,c(6,7)]
names(DEall) <- c("Label", "SYMBOL")
de <- split(DEall, DEall$Label)
list.de <- NULL
for (i in 1:length(de)){
Label <- unique(de[[i]]$Label)
de[[i]] <- bitr(de[[i]]$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
de[[i]] <- cbind(de[[i]], Label)
de[[i]]$Label <- as.character(de[[i]]$Label)
list.de <- rbind(list.de, de[[i]])
}
gobp.mono1 <- compareCluster(ENTREZID~Label, data= list.de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP",readable=TRUE)
dotplot(gobp.mono1)+theme(axis.text.x = element_text(angle = 50,hjust = 1))
DEall <- top10[,c(6,7)]
names(DEall) <- c("Label", "SYMBOL")
de <- split(DEall, DEall$Label)
list.de <- NULL
for (i in 1:length(de)){
Label <- unique(de[[i]]$Label)
de[[i]] <- bitr(de[[i]]$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
de[[i]] <- cbind(de[[i]], Label)
de[[i]]$Label <- as.character(de[[i]]$Label)
list.de <- rbind(list.de, de[[i]])
}
kegg.mono1 <- compareCluster(ENTREZID~Label, data=list.de, fun="enrichKEGG")
dotplot(kegg.mono1)+theme(axis.text.x = element_text(angle = 50,hjust = 1))
gobp.mono1 <- compareCluster(ENTREZID~Label, data= list.de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP",readable=TRUE)
dotplot(gobp.mono1)+theme(axis.text.x = element_text(angle = 50,hjust = 1))
table(anding$responser)

hc <- subset(anding, time == ctrl)
hc <- subset(anding, time == "ctrl")
DimPlot(hc,split.by="sampid",ncol=2)
table(hc$sample)
DimPlot(hc,split.by="sample",ncol=2)


saveRDS(anding,file = "./pbmc.naicd4t.rds")



################################
library(SingleR)
anding <- readRDS("pbmc.naicd4t.rds")
input <- GetAssayData(object = anding, slot = "data", assay = "RNA")

######## optain references ######

hpca.se <- HumanPrimaryCellAtlasData()  #Obtain the HPCA data
blueprint.se <- BlueprintEncodeData(rm.NA="rows") #Obtain human bulk RNA-seq data from Blueprint and ENCODE



singleR.list <- list()
singleR.list$hpca <- SingleR(test = input, 
                             method="single",
                             fine.tune=FALSE,
                             ref = hpca.se, 
                             labels = hpca.se$label.fine)

singleR.list$blueprint <- SingleR(test = input, 
                                  method="single",
                                  fine.tune=FALSE,
                                  ref = blueprint.se, 
                                  labels = blueprint.se$label.fine)


anding$hpca.labels <- singleR.list$hpca$labels
anding$blueprint.labels <- singleR.list$blueprint$labels




####

Idents(anding) <- "blueprint.labels"
anding <- RenameIdents(anding,
`CD4+ T-cells`="CD4+ naive T",
`CD4+ Tcm` = "CD4+ Tcm",
`CD8+ T-cells` = "CD8+ naive T",
`CD8+ Tcm` = "CD8+ Tcm", 
`CLP` = "undefine",  #set to undefine with too weak supports (small number of cells)
`CD8+ Tem`="undefine",#set to undefine with too weak supports (small number of cells)
`Tregs`="undefine",#set to undefine with too weak supports (small number of cells)
`CD4+ Tem`="undefine"#set to undefine with too weak supports (small number of cells)
)
anding$blueprint <- Idents(anding)


#Idents(pbmc.pools)
Idents(anding) <- "blueprint"
#Idents(pbmc.pools) <- "monaco"
tmp<-data.frame(cell_type= anding@active.ident,sample= anding@meta.data$sample, donor= anding@meta.data$cond.res)
res<-matrix(NA, ncol=4, nrow=length(unique(tmp$cell_type))*length(unique(tmp$sample)))
res<-as.data.frame(res)
names(res)<-c("cell_type","freq","sample","Condition")
n=1;
for (d in unique(tmp$sample)){
  t1<-subset(tmp,sample==d)
  t2<-as.data.frame(table(t1$cell_type))
  res[n:(n+length(unique(tmp$cell_type))-1),1]<-as.character(t2$Var1)
  res[n:(n+length(unique(tmp$cell_type))-1),2]<-t2$Freq
  res[n:(n+length(unique(tmp$cell_type))-1),3]<-rep(as.character(d),length(unique(tmp$cell_type)))
  res[n:(n+length(unique(tmp$cell_type))-1),4]<-rep(as.character(unique(t1$donor)),length(unique(tmp$cell_type)))
  n=n+length(unique(tmp$cell_type));
}
x<-split(res,res$sample)
y<-lapply(x, function(w) {w=transform(w,f=freq/sum(freq))})
y2<-do.call(rbind,y)

y2$Condition <- factor(y2$Condition, levels=c("HC","t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response"))
y2$other <- 1-y2$f

p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_bw()+scale_fill_manual(values=c("grey","#fff7bc","#fee391","#fec44f","#c994c7","#df65b0","#ce1256"))+xlab("")+ylab("percent")+theme(axis.text.x=element_text(angle = 70,hjust = 1))




################

tmp<-data.frame(cell= row.names(pbmc.pools@meta.data), blueprint = pbmc.pools@meta.data$blueprint, monaco = pbmc.pools@meta.data$monaco)
row.names(tmp) <- tmp$cell
pbmc.all <- AddMetaData(pbmc.all,metadata = tmp)

library(forcats)
pbmc.all$monaco <- fct_explicit_na(pbmc.all$monaco, "other")
pbmc.all$blueprint <- fct_explicit_na(pbmc.all$blueprint, "other")


pbmc.all$cond.res <- paste(pbmc.all$time, pbmc.all$responser, sep="_")
pbmc.all$cond.res[pbmc.all$cond.res == "ctrl_ctrl"]="ctrl"
pbmc.all$cond.res <- factor(pbmc.all$cond.res , levels=c("t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response","ctrl"))
table(pbmc.all$cond.res)

Idents(pbmc.all) <- "blueprint"
tmp<-data.frame(cell_type= pbmc.all@active.ident,sample= pbmc.all@meta.data$sample, donor= pbmc.all@meta.data$cond.res)
res<-matrix(NA, ncol=4, nrow=length(unique(tmp$cell_type))*length(unique(tmp$sample)))
res<-as.data.frame(res)
names(res)<-c("cell_type","freq","sample","Condition")
n=1;
for (d in unique(tmp$sample)){
  t1<-subset(tmp,sample==d)
  t2<-as.data.frame(table(t1$cell_type))
  res[n:(n+length(unique(tmp$cell_type))-1),1]<-as.character(t2$Var1)
  res[n:(n+length(unique(tmp$cell_type))-1),2]<-t2$Freq
  res[n:(n+length(unique(tmp$cell_type))-1),3]<-rep(as.character(d),length(unique(tmp$cell_type)))
  res[n:(n+length(unique(tmp$cell_type))-1),4]<-rep(as.character(unique(t1$donor)),length(unique(tmp$cell_type)))
  n=n+length(unique(tmp$cell_type));
}
x<-split(res,res$sample)
y<-lapply(x, function(w) {w=transform(w,f=freq/sum(freq))})
y2<-do.call(rbind,y)

y2$Condition <- factor(y2$Condition, levels=c("t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response","ctrl"))
y2$other <- 1-y2$f

y3 <- y2[y2$cell_type != "other",]
p <- ggplot(data=y3, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_bw()+scale_fill_manual(values=c("#fff7bc","#fee391","#fec44f","#c994c7","#df65b0","#ce1256","grey"))+xlab("")+ylab("percent")+theme(axis.text.x=element_text(angle = 70,hjust = 1))




anding <- subset(pbmc.pools, blueprint == "CD4+ naive T")

anding$cond.res <- paste(anding$time, anding$responser, sep="_")
anding$cond.res[anding$cond.res == "ctrl_ctrl"]="HC"
anding$cond.res <- factor(anding$cond.res , levels=c("HC","t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response"))
table(anding$cond.res)
Idents(anding) <- "cond.res"


Idents(anding) <- "cond.res"
allmk <- FindAllMarkers(anding,logfc.threshold = 0.1, only.pos = T)

top10.neu <- allmk %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


###Fig 5B
DotPlot(anding,features = unique(top10.neu$gene))+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("conditions")+ylab("")


library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)


#sig.mk <- allmk[allmk$p_val_adj< 0.05,]

DEall <- sig.mk[,c(6,7)]
names(DEall) <- c("Label", "SYMBOL")

de <- split(DEall, DEall$Label)
list.de <- NULL
for (i in 1:length(de)){
	Label <- unique(de[[i]]$Label)
	de[[i]] <- bitr(de[[i]]$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	de[[i]] <- cbind(de[[i]], Label)
	de[[i]]$Label <- as.character(de[[i]]$Label)
	list.de <- rbind(list.de, de[[i]])
}


head(list.de)

hall <- read.gmt("~/Downloads/h.simple.v7.4.entrez.gmt")
enrich <- compareCluster(ENTREZID~Label, data= list.de, fun="enricher", , TERM2GENE= hall, pvalueCutoff = 0.05)







###############

aa <- levels(as.factor(cd4t$cond.res))[-1]
degs <- NULL
Idents(cd4t) <- "cond.res"
for(i in aa){
	
	tg <- FindMarkers(cd4t, ident.1=i, ident.2="ctrl ctrl", logfc.threshold=0.01)
	tg$gene <- row.names(tg)
	tg$label <- i
	degs <- rbind(degs,tg)
}

head(degs)

####

Idents(cd4t) <- "sample"
#exp <- as.data.frame(covid@assays$RNA@data)
cd4t.sub <- cd4t["MYC",]
dat <- NULL
dnames <- NULL
for(index1 in levels(as.factor(cd4t.sub$sample))){
  #index2 <- gsub("-","_", index1)
  subD <- subset(cd4t.sub, idents = index1)
  #sdat <- 
  sdat <- rowMeans(subD@assays$RNA@data)
  dat <- cbind(dat,sdat)
  dnames <- c(dnames, index1)
 }
 colnames(dat) <- dnames
 
dat2 <- data.frame("MYC"=dat[1,],"label"=c(rep(c("T1","T2","T3"),7),rep("ctrl",7)))
dat2$res = c(rep("_NR",3),rep("_Res",12),rep("_NR",6),rep("s",7))
dat2$cond <- paste(dat2$label,dat2$res,sep="")
dat2$cond <- factor(dat2$cond, levels=c("T1_NR","T2_NR","T3_NR","T1_Res","T2_Res","T3_Res","ctrls"))

ggplot(dat2,aes(x=cond,y=MYC))+geom_boxplot(aes(fill=cond), outlier.shape=NA)+geom_jitter(width = 0.2)+scale_fill_manual(values=c("#fff7bc","#fee391","#fec44f","#c994c7","#df65b0","#ce1256","grey"))+theme_classic()

###############







###### cc.genes


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cd4nt <- CellCycleScoring(cd4nt, s.features = s.genes, g2m.features = g2m.genes)

df <- data.frame(row.names= row.names(cd4nt@meta.data), cluster = cd4nt@meta.data$cond.res, Phase = cd4nt@meta.data$Phase)
df <- df[df$Phase != "Undecided",]
df %>% group_by(cluster, Phase) %>% summarise(n = n()) %>% mutate(freq = (n/sum(n)*100)) -> df
df <- as.data.frame(df)

df$cluster <- factor(df$cluster, levels=c("t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response","ctrl"))
df$Phase <- factor(df$Phase, levels = c("G1", "S", "G2M"))
gg <- ggplot(df, aes(x = cluster, y = freq, group = Phase, fill = Phase)) + geom_bar(stat="identity", position="dodge", color="black", size = 1) + scale_fill_manual(values = c("grey90", "grey50", "black")) + scale_y_continuous(expand = c(0,0), limits = c(0,50)) + ylab("% of cells in cell cylce phase") + ggtitle("Seurat") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))
plot(gg)



#########



s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cd4 <- CellCycleScoring(cd4, s.features = s.genes, g2m.features = g2m.genes)

df <- data.frame(row.names= row.names(cd4@meta.data), cluster = cd4@meta.data$State, Phase = cd4@meta.data$Phase)
df <- df[df$Phase != "Undecided",]
df %>% group_by(cluster, Phase) %>% summarise(n = n()) %>% mutate(freq = (n/sum(n)*100)) -> df
df <- as.data.frame(df)

#df$cluster <- factor(df$cluster, levels=c("t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response","ctrl"))
df$Phase <- factor(df$Phase, levels = c("G1", "S", "G2M"))
gg <- ggplot(df, aes(x = cluster, y = freq, group = Phase, fill = Phase)) + geom_bar(stat="identity", position="dodge", color="black", size = 1) + scale_fill_manual(values = c("grey90", "grey50", "black")) + scale_y_continuous(expand = c(0,0), limits = c(0,50)) + ylab("% of cells in cell cylce phase") + ggtitle("Seurat") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))
plot(gg)




####
Idents(cd4) <- "Phase"
tmp<-data.frame(cell_type= cd4@active.ident,sample= cd4@meta.data$sample, donor= cd4@meta.data$State)
res<-matrix(NA, ncol=4, nrow=length(unique(tmp$cell_type))*length(unique(tmp$sample)))
res<-as.data.frame(res)
names(res)<-c("cell_type","freq","sample","Condition")
n=1;
for (d in unique(tmp$sample)){
  t1<-subset(tmp,sample==d)
  t2<-as.data.frame(table(t1$cell_type))
  res[n:(n+length(unique(tmp$cell_type))-1),1]<-as.character(t2$Var1)
  res[n:(n+length(unique(tmp$cell_type))-1),2]<-t2$Freq
  res[n:(n+length(unique(tmp$cell_type))-1),3]<-rep(as.character(d),length(unique(tmp$cell_type)))
  res[n:(n+length(unique(tmp$cell_type))-1),4]<-rep(as.character(unique(t1$donor)),length(unique(tmp$cell_type)))
  n=n+length(unique(tmp$cell_type));
}
x<-split(res,res$sample)
y<-lapply(x, function(w) {w=transform(w,f=freq/sum(freq))})
y2<-do.call(rbind,y)

#y2$Condition <- factor(y2$Condition, levels=c("t1_nonres", "t2_nonres", "t3_nonres", "t1_response", "t2_response", "t3_response","ctrl"))
#y2$other <- 1-y2$f

p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_bw()+scale_fill_manual(values=c("#fff7bc","#fee391","#fec44f","#c994c7","#df65b0","#ce1256","grey"))+xlab("")+ylab("percent")+theme(axis.text.x=element_text(angle = 70,hjust = 1))

#dev.off()
