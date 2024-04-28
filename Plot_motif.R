library(reshape2)

setwd("./ATAC_seq/motif_enrichment/")

motif <- read.table("motif_Res_bg/knownResults.txt",sep="\t", header=T, comment.char = "")
motif <- read.table("motif_NR_bg/knownResults.txt",sep="\t", header=T, comment.char = "")
motif <- read.table("motif_T1_bg/knownResults.txt",sep="\t", header=T, comment.char = "")
motif <- read.table("motif_HC_bg/knownResults.txt",sep="\t", header=T, comment.char = "")



colnames(motif) <- c("name","seqs","Pvalue","LogPvalue","qvalue","Ntarget","Ptarget","NBackground","PBackground")
motif[,c("Mname","Chips","anno")] <- do.call(rbind,strsplit(motif$name, "/"))
motif$p <- exp(motif$LogPvalue)
motif$Padj <- p.adjust(motif$p,method="BH")
motif2 <- motif[1:10,c("Mname","Padj","Ptarget","PBackground")]
motif2$Ptarget <- sapply(motif2$Ptarget, function(x) as.numeric(sub("%","",x))/100)
motif2$PBackground <- sapply(motif2$PBackground, function(x) as.numeric(sub("%","",x))/100)
motif2$P.adjust <- as.character(signif(motif2$Padj,digits=3))
motif3 <- reshape2::melt(motif2,id.vars = c("Mname","Padj","P.adjust"), measure.vars = c("Ptarget","PBackground"), variable.name = "label", value.name = "Proportion")
#motif3$Proportion <- sapply(motif3$Proportion, function(x) as.numeric(sub("%","",x))/100)
aa <- sapply(motif3$Proportion, function(x) as.numeric(sub("%","",x))/100)
p <- ggplot(data=motif3, aes(x=Mname,y=Proportion,label=P.adjust,fill=label))+geom_bar(stat="identity",position="dodge")
p2 <- p+ scale_x_discrete(limits = rev(motif2$Mname))+theme_classic()+theme(axis.text = element_text(size=15))+scale_fill_manual(values=c("Ptarget"="#FF7F00","PBackground"="lightblue"))

p2+geom_text(aes(label=ifelse(label=="Ptarget",paste("P.adj = ",P.adjust,sep=""),''),y=0.45),size=5)+ylim(0,0.62)+coord_flip()



