---
title: "Salmon_Project1"
author: "Teerapon Sahwangarrom"
date: "27/01/2020"
output: html_document
---

```{r}
library(tximportData)
library(readr)
library(tximport)
tx2gene <- read.table("gencode.v32.metadata.HGNC", sep = "\t",head=F)
class(tx2gene)
names(tx2gene) <- c("TXNAME", "GENE")
tx2gene[1:10,1:2]
```

```{r}
seqDirs <- c("SRR3232052_paired_quant","SRR3232052_unpaired_quant","SRR3232053_paired_quant","SRR3232053_unpaired_quant","SRR5796680_paired_quant","SRR5796680_unpaired_quant","SRR5796681_paired_quant","SRR5796681_unpaired_quant","SRR5796682_paired_quant","SRR5796682_unpaired_quant","SRR5796683_paired_quant","SRR5796683_unpaired_quant")
seqDirs
seqPaths <- paste("/Users/teeraponsahwangarrom/quants",seqDirs,"quant.sf",sep="/")
txi <- tximport(seqPaths,type="salmon",tx2gene=tx2gene,ignoreAfterBar = T)
head(txi$counts)
dim(txi$counts)
class(txi$counts)
```

```{r}
txi.data <- data.frame(SRR5796680=apply(txi$counts[,5:6],1,sum),SRR5796681=apply(txi$counts[,7:8],1,sum),SRR5796682=apply(txi$counts[,9:10],1,sum),SRR5796683=apply(txi$counts[,11:12],1,sum))
prad.data <- as.matrix(txi.data)
```

```{r}
prad.length <- data.frame(SRR5796680=apply(txi$length[,5:6],1,sum),SRR5796681=apply(txi$length[,7:8],1,sum),SRR5796682=apply(txi$length[,9:10],1,sum),SRR5796683=apply(txi$length[,11:12],1,sum))
prad.length <- as.matrix(prad.length)
```

```{r}
library(edgeR)
prad.length <- prad.length/exp(rowMeans(log(prad.length)))
norm.prad.data <- prad.data/prad.length
eff.lib <- calcNormFactors(norm.prad.data) * colSums(norm.prad.data)

prad.length <- sweep(prad.length,2,eff.lib,"*")
prad.length <- log(prad.length)
y <- DGEList(prad.data)
y <- scaleOffset(y,prad.length)
keep <- filterByExpr(y)
y <- y[keep,]
```

```{r}

sampleType <- factor(c("treatment","treatment","control","control"))
##sampleReplicate <- c("rep1","rep2","rep1","rep2","rep1","rep2")
designMat <- model.matrix(~0+sampleType)
designMat

y <- estimateGLMCommonDisp(y,design=designMat)
y <- estimateGLMTrendedDisp(y, design=designMat)
y <- estimateGLMTagwiseDisp(y,design=designMat)
plotBCV(y)

fit <- glmFit(y,designMat)
contrast_treatedAR_vs_treatedIgG <- glmLRT(fit,contrast=makeContrasts(sampleTypetreatment-sampleTypecontrol, levels = designMat))
Salmon_DEG <- topTags(contrast_treatedAR_vs_treatedIgG,n = nrow(y$counts))
dim(Salmon_DEG)
Salmon_DEG["NORAD",]
write.csv(Salmon_DEG,"/Users/teeraponsahwangarrom/Salmon_DEG_RIP-seq_Project1.csv")

```

```{r}
library(ggplot2)
threshold_OE <- DEGgenes$table$FDR < 0.05
length(which(threshold_OE))
DEGgenes$table$threshold <- threshold_OE

DEGgenes_ordered <- DEGgenes$table[order(DEGgenes$table$FDR),]
DEGgenes_ordered$genelabels <- ""
DEGgenes_ordered$genelabels[1:10] <- rownames(DEGgenes_ordered)[1:10]

ggplot(DEGgenes_ordered) + geom_point(aes(x=logFC, y= -log10(FDR), colour= threshold)) + geom_text_repel(aes(x=logFC,y=-log10(FDR),label=ifelse(genelabels==T,rownames(DEGgenes_ordered),""))) + ggtitle("Differential expression genes") + xlab("log2 fold change") + ylab("-log10 FDR") + theme(legend.position= "none", plot.title=element_text(size=rel(1.5), hjust= 0.5), axis.title=element_text(size=rel(1.25)))
```

```{r}
library(DESeq2)
sampleTable <- data.frame(condition = factor(c("control","control","control","control","treatment","treatment","treatment", "treatment", "control", "control","control","control")))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi,sampleTable,~condition)
head(assay(dds))
dds.prad <- as.matrix(assay(dds))
dds.prad.data <- data.frame(SRR3232052= apply(dds.prad[,1:2],1,sum),SRR3232053=apply(dds.prad[,3:4],1,sum),SRR5796680=apply(dds.prad[,5:6],1,sum),SRR5796681=apply(dds.prad[,7:8],1,sum))
head(dds.prad.data)
```
```{r}
sampleTable <- data.frame(condition = c("control","control","treatment","treatment"))
dds.prad.data <- as.matrix(dds.prad.data)
dds.prad.data <- DESeqDataSetFromMatrix(countData=dds.prad.data,colData = sampleTable, design=~ condition)
dds.prad.data

keep <- rowSums(counts(dds.prad.data)) >=10
dds.prad.data <- dds.prad.data[keep,]
```

```{r}
dds <- DESeq(dds.prad.data)
DEG_result <- results(dds)
DEG_result
result_Ordered <- DEG_result[order(DEG_result$pvalue),]
write.csv(as.data.frame(result_Ordered),file="DEG_result_Project1_DESeq2.csv")
```

```{r}
library(dplyr)
library(org.Hs.eg.db)
DEG_result$entrezID <- mapIds(org.Hs.eg.db,keys=rownames(DEG_result),column = ("ENTREZID"),keytype = "SYMBOL")
head(DEG_result)
```

```{r}
geneList <- DEG_result[,2]
names(geneList) <- DEG_result[,7]
geneList <- sort(geneList, decreasing = T)
head(geneList)
gene <- names(geneList)[abs(geneList) > 2]
head(gene)

y <- gsePathway(geneList,nPerm=10000,pvalueCutoff=0.2,pAdjustMethod="BH",verbose=F)
res <- as.data.frame(y)
head(res)
write.csv(res,"/Users/teeraponsahwangarrom/GSEA_data_Salmon_project1.csv")

library(clusterProfiler)
library(enrichplot)
dotplot(y,showCategory=30)
ridgeplot(y)

library(ggupset)
upsetplot(y)

```

```{r}
kk2 <- gseKEGG(geneList = geneList, organism = "hsa", nPerm = 1000, minGSSize = 120, pvalueCutoff = 0.05, verbose = F)
head(kk2)

ego3 <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 1000, minGSSize = 100, maxGSSize = 500, pvalueCutoff = 0.05,verbose = F)
```
```{r}
Salmon_DEG <- as.data.frame(Salmon_DEG)
Salmon_DEG["NORAD",]
library(EnhancedVolcano)
EnhancedVolcano(Salmon_DEG,
    lab = rownames(Salmon_DEG),
    x = 'logFC',
    y = 'FDR',
    selectLab = c("NORAD","MT-ATP6","MT-ND2","MT-ND5","MT-CO2","MT-ND1","CEP78","PAK4","MSH3","RNF216","FAM162A","HMGN1","SRP9","GULP1","H2AFZ"),
    xlim = c(-8,8),
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 2.5,
    labSize = 2.5,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black')
```

