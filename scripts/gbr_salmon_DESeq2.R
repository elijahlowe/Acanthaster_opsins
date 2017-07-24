library("tximportData")
library("DESeq2")

setwd("/Users/elijahlowe/Desktop/Acantaster/data/counts/gbr_mrna_salmon/")
samples <- c("MIX1_S6", "MIX2_S7", "MIX3_S8",
             "EYES1A_S11","EYES1B_S12","EYES1D",
             "RN1_S3", "RN2_S4", "RN3_S5",
             "TF1_S1", "TF2_S9", "TF3_S10"
)
files<-paste(samples, ".quant.sf", sep="")

read.sample <- function(sample.name) {
  file.name <- paste(sample.name, ".quant.sf", sep="")
  result <- read.delim(file.name, head=TRUE,sep="\t", colClasses=c(rep("character",1),rep("numeric",4)))
}

Acan.mix1 <- read.sample(samples[1])
head(Acan.mix1[order(Acan.mix1[,1]),])
nrow(Acan.mix1)

#Read the second sample
Acan.mix2 <- read.sample(samples[2])

#Let's make sure the first and second samples have the same number of rows and the same genes in each row
nrow(Acan.mix1) == nrow(Acan.mix2)
all(Acan.mix1[,1] == Acan.mix2[,1])

all.data <- Acan.mix1
all.data <- cbind(Acan.mix1[, c(1, 5)], Acan.mix2$NumReads)
all.tpm <- cbind(Acan.mix1[, c(1,4)], Acan.mix2$TPM) 
for (c in 3:length(samples)) {
  temp.data <- read.sample(samples[c])
  all.data <- cbind(all.data, temp.data$NumReads)
  all.tpm <- cbind(all.tpm, temp.data$TPM)
}

head(all.data)
head(all.tpm)
colnames(all.data)[2:ncol(all.data)] <- samples
colnames(all.tpm)[2:ncol(all.tpm)] <- samples
rownames(all.tpm) <-all.tpm$Name
all.tpm <- all.tpm[,2:ncol(all.tpm)]

raw.deseq.data <- all.data[,2:ncol(all.data)]
raw.deseq.data <- sapply(raw.deseq.data,as.integer)

#Set row names to the gene names
rownames(raw.deseq.data) <- all.data$Name

head(raw.deseq.data)
Acan.design <- data.frame(
  row.names=samples,
  condition=c("mix","mix","mix","eyes","eyes","eyes",
              "rn","rn","rn","tf","tf","tf"),
  libType=rep("paired-end", 12)
)
#Double check it...
Acan.design

(deseq.data <- DESeqDataSetFromMatrix(countData = raw.deseq.data,
                                      colData = Acan.design,
                                      design = ~condition))

ddsMF1 <- DESeq(deseq.data)
rld_reduced <- rlog(ddsMF1)
plotPCA( rld_reduced, intgroup = c("condition"), title(main="PCA Acantaster tissues"))

eyes <- results( ddsMF1,contrast=c("condition","eyes","mix"))
rn <- results( ddsMF1,contrast=c("condition","rn","mix"))
tf <- results( ddsMF1,contrast=c("condition","tf","mix"))
eyesVrn <- results( ddsMF1,contrast=c("condition","eyes","rn"))
eyesVtf <- results( ddsMF1,contrast=c("condition","eyes","tf"))
rnVtf <- results( ddsMF1,contrast=c("condition","rn","tf"))

summary(eyes)
summary(rn)
summary(tf)
summary(eyesVrn)
summary(eyesVtf)
summary(rnVtf)

ciliary_opsin=c("gbr.65.46.t1","gbr.65.45.t1","gbr.508.2.t1")
melatonin=c("gbr.212.22.t1")
GO_c_opsin=c("gbr.470.6.t1")
echinopsinsB=c("gbr.176.10.t1")
neuropsins=c("gbr.35.57.t1")
peropsins=c("gbr.31.87.t1")
rhabdomeric=c("gbr.176.5.t1")
rgr_opsins=c("gbr.37.118.t1")


png("colored_opsins_29.5.17.png", width = 11.37, height = 7.5, units = 'in', res = 300)
par(mfrow=c(2,3))
plotMA(eyes, main="Eyes vs mixed", ylim=c(-12,12),colNonSig = "azure2",colSig="azure3")
with(eyes[ciliary_opsin, ], {
  points(baseMean, log2FoldChange, col="red", cex=1.5, lwd=2,pch=c(18,9,23))
})
with(eyes[GO_c_opsin, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2,pch=16)
})
#with(eyes[melatonin, ], {
#  points(baseMean, log2FoldChange, col="dimgrey", cex=1.5, lwd=1.5,pch=16)
#})
with(eyes[echinopsinsB, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2,pch=16)
})
with(eyes[neuropsins, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2,pch=15)
})
with(eyes[peropsins, ], {
  points(baseMean, log2FoldChange, col="yellow", cex=1.5, lwd=2,pch=16)
})
with(eyes[rhabdomeric, ], {
  points(baseMean, log2FoldChange, col="blue", cex=1.5, lwd=2,pch=18)
})
with(eyes[rgr_opsins, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2,pch=16)
})
#legend("right", inset=-.05, title="Opsin type",
#       o_type, fill=c(rep("red",3),"green","black","purple","yellow","blue","orange"))

plotMA(rn, main="Radial nerve vs mixed", ylim=c(-12,12),colNonSig = "azure2",colSig="azure3")
with(rn[ciliary_opsin, ], {
  points(baseMean, log2FoldChange, col="red", cex=1.5, lwd=2,pch=8)
})
with(rn[GO_c_opsin, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2,pch=16)
})
with(rn[melatonin, ], {
  points(baseMean, log2FoldChange, col="dimgrey", cex=1.5, lwd=2,pch=17)
})
with(rn[echinopsinsB, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=1.5,pch=17)
})
with(rn[neuropsins, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2,pch=15)
})
with(rn[peropsins, ], {
  points(baseMean, log2FoldChange, col="yellow", cex=1.5, lwd=2,pch=16)
})
with(rn[rhabdomeric, ], {
  points(baseMean, log2FoldChange, col="blue", cex=1.5, lwd=2,pch=18)
})
with(rn[rgr_opsins, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2,pch=16)
})

plotMA(tf, main="Tube feet vs mixed", ylim=c(-12,12),colNonSig = "azure2",colSig="azure3")
with(tf[ciliary_opsin, ], {
  points(baseMean, log2FoldChange, col="red", cex=1.5, lwd=2,pch=8)
})
with(tf[GO_c_opsin, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2,pch=16)
})
with(tf[melatonin, ], {
  points(baseMean, log2FoldChange, col="dimgrey", cex=1.5, lwd=1.5,pch=17)
})
with(tf[echinopsinsB, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2,pch=17)
})
with(tf[neuropsins, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2,pch=15)
})
with(tf[peropsins, ], {
  points(baseMean, log2FoldChange, col="yellow", cex=1.5, lwd=2,pch=16)
})
with(tf[rhabdomeric, ], {
  points(baseMean, log2FoldChange, col="blue", cex=1.5, lwd=2,pch=18)
})
with(tf[rgr_opsins, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2,pch=16)
})

plotMA(eyesVtf, main="Eyes vs Tube feet", ylim=c(-12,12),colNonSig = "azure2",colSig="azure3")
with(eyesVtf[ciliary_opsin, ], {
  points(baseMean, log2FoldChange, col="red", cex=1.5, lwd=2,pch=8)
})
with(eyesVtf[GO_c_opsin, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2,pch=16)
})
with(eyesVtf[melatonin, ], {
  points(baseMean, log2FoldChange, col="dimgrey", cex=1.5, lwd=1.5,pch=16)
})
with(eyesVtf[echinopsinsB, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2,pch=17)
})
with(eyesVtf[neuropsins, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2,pch=15)
})
with(eyesVtf[peropsins, ], {
  points(baseMean, log2FoldChange, col="yellow", cex=1.5, lwd=2,pch=16)
})
with(eyesVtf[rhabdomeric, ], {
  points(baseMean, log2FoldChange, col="blue", cex=1.5, lwd=2,pch=18)
})
with(eyesVtf[rgr_opsins, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2,pch=16)
})

plotMA(eyesVrn, main="Eyes vs Radial nerve", ylim=c(-12,12),colNonSig = "azure2",colSig="azure3")
with(eyesVrn[ciliary_opsin, ], {
  points(baseMean, log2FoldChange, col="red", cex=1.5, lwd=2,pch=8)
})
with(eyesVrn[GO_c_opsin, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2,pch=16)
})
with(eyesVrn[melatonin, ], {
  points(baseMean, log2FoldChange, col="dimgrey", cex=1.5, lwd=1.5,pch=17)
})
with(eyesVrn[echinopsinsB, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2,pch=17)
})
with(eyesVrn[neuropsins, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2,pch=15)
})
with(eyesVrn[peropsins, ], {
  points(baseMean, log2FoldChange, col="yellow", cex=1.5, lwd=2,pch=16)
})
with(eyesVrn[rhabdomeric, ], {
  points(baseMean, log2FoldChange, col="blue", cex=1.5, lwd=2,pch=18)
})
with(eyesVrn[rgr_opsins, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2,pch=16)
})

plotMA(rnVtf, main="Radial nerve vs Tube feet", ylim=c(-12,12),colNonSig = "azure2",colSig="azure3")
with(rnVtf[ciliary_opsin, ], {
  points(baseMean, log2FoldChange, col="red", cex=1.5, lwd=2,pch=8)
})
with(rnVtf[GO_c_opsin, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2,pch=16)
})
with(rnVtf[melatonin, ], {
  points(baseMean, log2FoldChange, col="dimgrey", cex=1.5, lwd=1.5,pch=17)
})
with(rnVtf[echinopsinsB, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2,pch=17)
})
with(rnVtf[neuropsins, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2,pch=15)
})
with(rnVtf[peropsins, ], {
  points(baseMean, log2FoldChange, col="yellow", cex=1.5, lwd=2,pch=16)
})
with(rnVtf[rhabdomeric, ], {
  points(baseMean, log2FoldChange, col="blue", cex=1.5, lwd=2,pch=18)
})
with(rnVtf[rgr_opsins, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2,pch=16)
})
dev.off()

###################################################################
eyes_up<-(rownames(subset(eyes, padj < 0.1 & log2FoldChange >= 1.5)))
rn_up<-(rownames(subset(rn, padj < 0.05 & log2FoldChange >= 1.5)))
tf_up<-(rownames(subset(tf, padj < 0.05 & log2FoldChange >= 1.5)))

eyes_down<-(rownames(subset(tf, padj < 0.05 & log2FoldChange <= -1.5)))
rn_down<-(rownames(subset(tf, padj < 0.05 & log2FoldChange <= -1.5)))
tf_down<-(rownames(subset(tf, padj < 0.05 & log2FoldChange <= -1.5)))

eyes_de<-(rownames(subset(eyes, padj < 0.1)))# & log2FoldChange <= -1.5 | log2FoldChange >= 1.5)))
rn_de<-(rownames(subset(rn, padj < 0.1)))# & log2FoldChange <= -1.5 | log2FoldChange >= 1.5)))
tf_de<-(rownames(subset(tf, padj < 0.1)))# & log2FoldChange <= -1.5 | log2FoldChange >= 1.5)))
write.table(eyes_up, file = "eyes_up_p0.1_1.5.ids", quote=FALSE, row.names = FALSE)
write.table(rn_de, file = "rn.ids", quote=FALSE, row.names = FALSE)
write.table(tf_de, file = "tf.ids", quote=FALSE, row.names = FALSE)

write.csv(eyes, file="eyes_de_salmon_DESeq2.csv", quote=FALSE)
tf_ids<-read.table("../../gbr-tf.id", header=FALSE, colClasses=c(rep("character",3),"numeric"))
eyes_TF<-Reduce(intersect,list(eyes_up,tf_ids$V3))
rn_up_TF<-Reduce(intersect,list(rn_up,tf_ids$V3))
tf_up_TF<-Reduce(intersect,list(tf_up,tf_ids$V3))

library("VennDiagram")
png("venn.png", width = 4.5, height = 4.5, units = 'in', res = 300)
venn.plot <- draw.triple.venn(length(eyes_up), length(rn_up), length(tf_up),
                              sum(eyes_up %in% rn_up), sum(rn_up %in% tf_up), sum(eyes_up %in% tf_up), length(Reduce(intersect,list(eyes_up,rn_up,tf_up))), c("Eyes", "Radial nerve", "Tube feet"))
grid.newpage()
grid.draw(venn.plot)
dev.off()


library(ggplot2)
library(reshape2)
opsins<-all.tpm[c(ciliary_opsin,GO_c_opsin,echinopsinsB,
                  neuropsins,peropsins,rhabdomeric,rgr_opsins),]
#opsins$group <- row.names(opsins)
#opsins.m <- melt(opsins, id.vars = "group")
#ggplot(opsins.m, aes(group, value)) + geom_boxplot()

tissue <- c(rep("mixed",3),rep("eyes",3),rep("rn",3),rep("tf",3))
o_type <- c("c-opsin1.1","c-opsin1.2","c-opsin1.3","go-opsin",
            "echinopsinsB","neuropsins","peropsins","rhabdomeric",
            "RGR opsins")

df<-data.frame(opsins,o_type)
names(df)[1:12]<-tissue
df

#boxplot(df)
#opsins$group <- row.names(o_type)
opsins.m <- melt(df, id.vars = "o_type")
#ggplot(opsins.m, aes(o_type, value, fill = variable)) + geom_boxplot()

opsins.m
library(gridExtra)
library(grid)
library(gtable)
opsins_top <- ggplot(opsins.m, aes(o_type, value, fill = variable,color=variable),log="value") + coord_cartesian(ylim = c(300,4000)) +
  geom_boxplot() + theme(legend.position="none", axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank()) 
opsins_bottom <- ggplot(opsins.m, aes(o_type, value, color=variable),log="value") + coord_cartesian(ylim = c(0,75)) + geom_boxplot() +
  labs(x = "opsin type", color="Tissue", y = "TPM" ) + theme(axis.text.x = element_text(angle=45, vjust=0.75, hjust=1))
gA <- ggplotGrob(opsins_top)
gB <- ggplotGrob(opsins_bottom)
gA$widths <- gB$widths

#grid.arrange(opsins_top, opsins_bottom, nrow = 2)
png("opsin_counts.png", width = 3.5, height = 3, units = 'in', res = 300)
grid.arrange(gA, gB, nrow = 2)
dev.off()


plotCounts(ddsMF1, gene=rhabdomeric, intgroup = c("conditions"))
