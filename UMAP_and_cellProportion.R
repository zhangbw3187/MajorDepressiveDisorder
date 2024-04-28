library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(tidyr)
options(future.globals.maxSize = 10737418240) 

anding <- readRDS("MDD_annotated_PBMC.flt.rds")


DimPlot(anding,split.by="disease")


pdf("UMAP_split.pdf",width=10,height=6)
DimPlot(anding,cols=umap.colors.ct,split.by="disease")
dev.off()


anding$time.res <- paste(anding$time, anding$responser)
 table(anding$time.res)
 
DimPlot(anding,cols=umap.colors.ct,split.by = "time.res",ncol=7)


#################################
FeaturePlot(anding, features = c("CD3D","SELL","CD8A","GZMA"),order = T)


p <- FeaturePlot(anding, features = c("CD3D","CCR7","IL7R","GZMK","GZMA","CD8A","NKG7","CD14","FCGR3A","CST3","CD79A","PPBP"),order = T)


p * (theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),  #remove y axis ticks
        axis.line=element_blank(),
        axis.title=element_blank()
        ))



Idents(anding)<- "celltypeL0"
pbmc.mks <- FindAllMarkers(object = anding,
only.pos = TRUE,
min.pct = 0.2,
logfc.threshold = 0.2,
min.diff.pct = 0.1,
test.use = "wilcox"
)


tt.markers <- pbmc.mks[grep("^RPS|^RPL",pbmc.mks$gene,invert = T),]
top10 <- pbmc.mks %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(anding,features = unique(top10$gene),group.by="new.id",cols="RdBu")+coord_flip()+xlab("")+ylab("")

##################################################
#. Calculate Cell Proportion
###########   T1 HC
library(DirichletReg)
library(ggsignif)

idName <- subset(anding, time %in% c("t1","ctrl"))


Idents(idName) <- "celltypeL0"
tmp<-data.frame(cell_type= idName@active.ident,sample= idName@meta.data$sample,donor= idName@meta.data$time)
res<-matrix(NA,ncol=4,nrow=length(unique(tmp$cell_type))*length(unique(tmp$sample)))
res<-as.data.frame(res)
names(res)<-c("cell_type","freq","sample","donor")
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
y<-lapply(x, function(w) {w=transform(w,freq_per_sample=freq/sum(freq))})
y2<-do.call(rbind,y)

################
p=ggplot(data=y2,aes(x=donor,y=freq_per_sample,fill=cell_type))+geom_boxplot()
p+theme(axis.text.x=element_text(angle = -60,vjust = -0.1,hjust = 0.1))+facet_grid(.~cell_type)

################
p <- ggboxplot(y2, x="donor", y="freq_per_sample", color="donor", palette=c("#00AFBB", "#E7B800"), add="jitter", facet.by="cell_type", short.panel.labs=F)
p + stat_compare_means(label="p.format")


names(y2) <- c("cell_type", "freq",  "sample", "donor", "f")
y2$other <- 1-y2$f

dat.test <- NULL

for(i in levels(as.factor(y2$cell_type))){
#for(i in c("PLA1","PLA2","PLA3")){  
  print(i)
  y2.i <- y2[y2$cell_type==i,]
  y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~donor,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["donort1",3:4]
  dat.test <- rbind(dat.test, c(i,sc,"T1vsHC"))
}
dat.test<- as.data.frame(dat.test)
names(dat.test) <- c("cell_type","Z","Pr","comparison")
dat.test$Pr <- as.numeric(as.character(dat.test$Pr))
dat.test$Z <- as.numeric(as.character(dat.test$Z))

dat.test$sig <- "ns"
dat.test$sig[dat.test$Pr < 0.05]="*"
dat.test$sig[dat.test$Pr < 0.01]="**"
dat.test$sig[dat.test$Pr < 0.001]="***"

dat.sig <- dat.test[dat.test$sig != "ns",]
dat.sig$st="HC"
dat.sig$ed="T1"
maxs <- NULL
for(i in dat.sig$cell_type){maxs <- c(maxs,max(y2$f[y2$cell_type == i])+0.01)}
dat.sig$maxs <- maxs

y2$Condition="HC"
y2$Condition[y2$donor=="t1"]="T1"
p <- ggplot(data=y2, aes(x= Condition,y=f))+geom_boxplot(aes(fill= Condition),outlier.shape=NA)
p2 <- p+theme(axis.text.x=element_text(angle = -60,vjust = -0.1,hjust = 0.1))+
  facet_grid(.~cell_type)+geom_jitter(shape=16, width = 0.2)+theme_bw()
p3 <- p2 + geom_signif(data=dat.sig, aes(xmin=st,xmax=ed, annotations=sig, y_position=maxs), manual = TRUE, vjust=0.8, tip_length=0.005)+#theme(axis.text.x=element_text(angle = 70,hjust = 1))
scale_fill_manual(values=c("#00AFBB", "#E70800"))

pdf("cpp_new.pdf", width=16, height=3)
plot(p3)
dev.off()

##################. all time
anding$condition <- paste(anding$time,anding$responser,sep="_")
anding$condition[anding$condition=="ctrl_ctrl"]="HC"

Idents(anding) <- "celltypeL0"
#tmp<-data.frame(cell_type= anding@active.ident,sample= anding@meta.data$sample,donor= anding@meta.data$condition)
tmp<-data.frame(cell_type= anding@active.ident,sample= anding@meta.data$sample,donor= anding@meta.data$time)
res<-matrix(NA,ncol=4,nrow=length(unique(tmp$cell_type))*length(unique(tmp$sample)))
res<-as.data.frame(res)
names(res)<-c("cell_type","freq","sample","donor")
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
y<-lapply(x, function(w) {w=transform(w,freq_per_sample=freq/sum(freq))})
y2<-do.call(rbind,y)

################
p=ggplot(data=y2,aes(x=donor,y=freq_per_sample,fill=cell_type))+geom_boxplot()
p+theme(axis.text.x=element_text(angle = -60,vjust = -0.1,hjust = 0.1))+facet_grid(.~cell_type)

################
mycomp <- list(c("t1","HC"),c("t2","HC"),c("t3","HC"))

p <- ggboxplot(y2, x="donor", y="freq_per_sample", color="donor", add="jitter", facet.by="cell_type", short.panel.labs=F,ylim = c(0, 0.8),palette = "nejm", font.label = list(size = 6),ylab="cell proportion",xlab="")
p + stat_compare_means(comparisons = mycomp, label="p.format")

##palette=c("#00AFBB", "#E7B800")

y2$Condition <- factor(y2$donor, levels= c( "HC","t1_nonres", "t2_nonres","t3_nonres","t1_response", "t2_response", "t3_response"))
p <- ggboxplot(y2, x="Condition",y="freq_per_sample",facet.by="cell_type", add="jitter",palette="lancet",color="Condition",ylab="cell proportion",xlab="")

mycomp <- list(c("t1_nonres","t1_response"),c("t2_nonres","t2_response"),c("t3_nonres","t3_response"),c("t1_nonres","HC"),c("t1_response","HC"))
p + stat_compare_means(comparisons = mycomp,label="p.format")


p <- ggboxplot(y2, x="donor", y="freq_per_sample", color="donor", palette="jco", add="jitter", facet.by="cell_type", short.panel.labs=F, ,ylab="cell proportion", xlab="stimuli",ylim = c(0, 0.8))
#p + stat_compare_means(label="p.format")
p + stat_compare_means(comparisons = mycomp,label="p.format")




