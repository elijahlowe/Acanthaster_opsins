source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)

############Import and build genome from GFF3 ############
setwd("~/Desktop/Acantaster/data/transcriptome/")

#gff_file <- system.file("extdata", "GFF3_files", "gbr-cotsv1.0.EVM2.edited.gff3", package = "GenomicFeatures")
txdb <- makeTxDbFromGFF("gbr-cotsv1.0.EVM2.edited.gff3", format="gff3")
txdb
gff.file <-"~/Desktop/Acantaster/data/transcriptome/gbr-cotsv1.0.EVM2.edited.gff3"
exonsBy(txdb, by="gene")
dm.genes <- genes(txdb)
head(dm.genes)

biocLite(karyoploteR)
library(karyoploteR)

#### combine gene location with expression ####
mcols(dm.genes) <- c(eyes[names(dm.genes), c("log2FoldChange", "stat", "pvalue", "padj")]) 
head(dm.genes, n=4)

gff.size <- read.table("~/Desktop/Acantaster/data/transcriptome/scaffold65.size", header = FALSE)
gff.size <- cbind(gff.size[1],rep(1,dim(gff.size)[1]),gff.size[2])
gbr.genome <- toGRanges(gff.size[,c(1,2,3)])

features <- import(gff.file)
table(features$type)
genes <- features[features$type=="gene"]

####### filter differentially expressed genes by FDR #########
filtered.dm.genes <- dm.genes[!is.na(dm.genes$padj)]
log.pval <- -log10(filtered.dm.genes$padj)
mcols(filtered.dm.genes)$log.pval <- log.pval
sign.genes <- filtered.dm.genes[filtered.dm.genes$padj < 0.05,]
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval)

####### create data frame for tpm and calucate mean & std #########
eyes_tpm.df <- data.frame(ID=rev(c("c-opsin 1.1a","c-opsin 1.2","c-opsin 1.3")), 
                          EYES1A_S1=all.tpm[c("gbr.65.45","gbr.65.46","gbr.65.47"),]$EYES1A_S11, 
                          EYES1B_S12=all.tpm[c("gbr.65.45","gbr.65.46","gbr.65.47"),]$EYES1B_S12, 
                          EYES1D=all.tpm[c("gbr.65.45","gbr.65.46","gbr.65.47"),]$EYES1D)
eyes_tpm.dfm<-melt(eyes_tpm.df, id.vars="ID")
c_opsin_tpm <- ddply(eyes_tpm.dfm,~ID,summarise,mean=mean(value),sd=sd(value))
rownames(c_opsin_tpm) <- (c_opsin_tpm$ID)
c_opsin_tpm <- c_opsin_tpm[,2:ncol(c_opsin_tpm)]

####### graph the data ######

cex.val <- sqrt(sign.genes$log.pval)/3
fc.ymax <- ceiling(max(abs(range(sign.genes$log2FoldChange))))
fc.ymin <- -fc.ymax
kp <- plotKaryotype(genome=gbr.genome,cytobands = dm.genes, plot.type = 4,chromosomes="gbr_scaffold65")
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r0=0, r1=0.5)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r0=0, r1=0.5,cex=0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r0=0, r1=0.5, cex=0.5)

####### Panel 2: TPM ##########
tpm.ymax <- ceiling(max(max(all.tpm[c("gbr.65.45","gbr.65.46","gbr.65.47"),]$EYES1A_S11)))
tpm.ymin <- 0
kpPoints(kp, dm.genes[c("gbr.65.45","gbr.65.46","gbr.65.47"),], 
         y=c_opsin_tpm[c("gbr.65.45","gbr.65.46","gbr.65.47"),]$mean,
         col="blue", ymax=tpm.ymax,
         ymin=tpm.ymin,pch="<-",cex = 0.75, r0=0.52, r1=1)
kpText(kp, dm.genes[c("gbr.65.45","gbr.65.46","gbr.65.47"),], 
         y=c_opsin_tpm[c("gbr.65.45","gbr.65.46","gbr.65.47"),]$mean,
         x=dm.genes[c("gbr.65.45","gbr.65.46","gbr.65.47"),],
         labels = rev(c("c-opsin 1.1a","c-opsin 1.2","c-opsin 1.3")),
         col="blue", ymax=tpm.ymax, ymin = tpm.ymin,
       cex=0.45, r0=0.52, r1=1, pos=4)

kpPlotMarkers(kp, dm.genes[c("gbr.65.45","gbr.65.46","gbr.65.47"),], labels = names(dm.genes[c("gbr.65.45","gbr.65.46","gbr.65.47"),]), text.orientation = "horizontal",r1=1.1)

