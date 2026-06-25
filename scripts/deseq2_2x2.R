#!/usr/bin/env Rscript
# DESeq2 2x2 (genotype WT/anac017 x treatment C/AA) on transcript counts (idxstats).
# Key term: genotype:treatment interaction = ANAC017-dependent response to Antimycin A.
suppressMessages({library(DESeq2); library(ggplot2); library(pheatmap); library(RColorBrewer)})
set.seed(42)
IDX <- path.expand("~/deseq2/idxstats"); OUT <- path.expand("~/deseq2/out"); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
samples <- c("WT_C_R1","WT_C_R2","WT_C_R3","WT_AA_R1","WT_AA_R2","WT_AA_R3",
             "anac017-1_C_R1","anac017-1_C_R2","anac017-1_C_R3",
             "anac017-1_AA_R1","anac017-1_AA_R2","anac017-1_AA_R3")
mats <- lapply(samples,function(s){d<-read.delim(file.path(IDX,paste0(s,".txt")),header=FALSE); setNames(d$V3,d$V1)})
tx <- Reduce(intersect, lapply(mats,names))
counts <- sapply(mats,function(m) as.numeric(m[tx])); rownames(counts)<-tx; colnames(counts)<-samples
counts <- counts[rowSums(counts)>=10,]; cat("transcripts kept:",nrow(counts),"\n")
geno <- factor(ifelse(grepl("^WT",samples),"WT","anac017"),levels=c("WT","anac017"))
trt  <- factor(ifelse(grepl("_AA_",samples),"AA","C"),levels=c("C","AA"))
cold <- data.frame(genotype=geno,treatment=trt,row.names=samples)
dds <- DESeqDataSetFromMatrix(round(counts),cold,~genotype+treatment+genotype:treatment)
dds <- DESeq(dds); vsd <- vst(dds,blind=FALSE)
cols <- c("WT.C"="#4C78A8","WT.AA"="#F58518","anac017.C"="#54A24B","anac017.AA"="#E45756")
grp <- paste(cold$genotype,cold$treatment,sep=".")
th <- theme_bw(base_size=16)+theme(plot.title=element_text(face="bold",size=18),axis.title=element_text(face="bold"))

# PCA
pc <- plotPCA(vsd,intgroup=c("genotype","treatment"),returnData=TRUE); pv<-round(100*attr(pc,"percentVar"))
g <- ggplot(pc,aes(PC1,PC2,color=paste(genotype,treatment,sep="."),shape=genotype))+
  geom_point(size=6,alpha=.9)+scale_color_manual(values=cols,name="group")+
  labs(title="DRS transcript-level PCA (2x2 design)",x=paste0("PC1: ",pv[1],"% var"),y=paste0("PC2: ",pv[2],"% var"))+th
ggsave(file.path(OUT,"D1_PCA.png"),g,width=8,height=6,dpi=350); ggsave(file.path(OUT,"D1_PCA.pdf"),g,width=8,height=6)

# dispersion QC
png(file.path(OUT,"D2_dispersion.png"),1000,800,res=130); plotDispEsts(dds); dev.off()

volcano <- function(res,ttl,fn){
  d <- as.data.frame(res); d <- d[!is.na(d$padj),]; d$sig <- d$padj<0.05 & abs(d$log2FoldChange)>1
  g <- ggplot(d,aes(log2FoldChange,-log10(padj),color=sig))+geom_point(size=1.2,alpha=.5)+
    scale_color_manual(values=c("FALSE"="grey75","TRUE"="#E45756"),guide="none")+
    geom_hline(yintercept=-log10(0.05),lty=2)+geom_vline(xintercept=c(-1,1),lty=2)+
    labs(title=ttl,subtitle=paste0(sum(d$sig)," sig (padj<0.05, |LFC|>1)"),x="log2 fold-change",y="-log10(padj)")+th
  ggsave(fn,g,width=7,height=6,dpi=350); ggsave(sub("png$","pdf",fn),g,width=7,height=6)
  write.csv(as.data.frame(res),sub(".png",".csv",fn)); sum(d$sig)
}
n1<-volcano(results(dds,name="treatment_AA_vs_C"),"Antimycin A effect (WT)",file.path(OUT,"D3_volcano_AA_in_WT.png"))
n2<-volcano(results(dds,name="genotype_anac017_vs_WT"),"anac017 vs WT (basal)",file.path(OUT,"D4_volcano_genotype.png"))
n3<-volcano(results(dds,name="genotypeanac017.treatmentAA"),"Genotype x Treatment (ANAC017-dependent AA response)",file.path(OUT,"D5_volcano_interaction.png"))
cat("sig AA-in-WT:",n1," genotype:",n2," interaction:",n3,"\n")

# heatmap top variable transcripts
topv <- head(order(rowVars(assay(vsd)),decreasing=TRUE),40)
mat <- assay(vsd)[topv,]; mat <- mat-rowMeans(mat)
ann <- data.frame(genotype=cold$genotype,treatment=cold$treatment,row.names=colnames(mat))
png(file.path(OUT,"D6_heatmap_top40.png"),1100,1300,res=140)
pheatmap(mat,annotation_col=ann,show_rownames=FALSE,fontsize=12,
         main="Top 40 most variable transcripts (vst, centred)",
         annotation_colors=list(genotype=c(WT="#4C78A8",anac017="#54A24B"),treatment=c(C="grey70",AA="#F58518")))
dev.off()
write.csv(counts,file.path(OUT,"transcript_counts_matrix.csv"))
cat("DESEQ2_DONE\n")
