---
title: "Salmon_RNA-seq_Project1"
author: "Teerapon Sahwangarrom"
date: "07/02/2020"
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
seqDirs <- c("SRR5905030_P_quant","SRR5905030_U_quant","SRR5905031_P_quant","SRR5905031_U_quant","SRR5905032_P_quant","SRR5905032_U_quant","SRR5905033_P_quant","SRR5905033_U_quant")
seqDirs
seqPaths <- paste("/Users/teeraponsahwangarrom/quants",seqDirs,"quant.sf",sep="/")
txi <- tximport(seqPaths,type="salmon",tx2gene=tx2gene,ignoreAfterBar = T)
head(txi$counts)
dim(txi$counts)
class(txi$counts)
```

```{r}
txi.data <- data.frame(SRR5905030=apply(txi$counts[,1:2],1,sum),SRR5905031=apply(txi$counts[,3:4],1,sum),SRR5905032=apply(txi$counts[,5:6],1,sum),SRR5905033=apply(txi$counts[,7:8],1,sum))
prad.data <- as.matrix(txi.data)
```

```{r}
prad.length <- data.frame(SRR5905030=apply(txi$length[,1:2],1,sum),SRR5905031=apply(txi$length[,3:4],1,sum),SRR5905032=apply(txi$length[,5:6],1,sum),SRR5905033=apply(txi$length[,7:8],1,sum))
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
group <- factor(c("treated_AR","treated_AR","untreated_AR","untreated_AR"))
design <- model.matrix( ~0+group)
head(design)
design
plotMDS(y, method = "bcv", col=as.numeric(y$samples$group))
legend("topleft", as.character(unique(y$samples$group)),col = 1:3,pch = 20)

y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

plotBCV(y)
fit <- glmFit(y,design)
contrast_untreatedAR_vs_treatedAR <- glmLRT(fit,contrast=makeContrasts(grouptreated_AR-groupuntreated_AR, levels = design))
Salmon_DEG <- topTags(contrast_untreatedAR_vs_treatedAR,n = nrow(y$counts))
Salmon_DEG["NORAD",]
Salmon_DEG[1:10,]
write.csv(Salmon_DEG,file="/Users/teeraponsahwangarrom/RNADEGfor_Project_1_Salmon.csv")
```
```{r}
library(DESeq2)
sampleTable <- data.frame(condition = factor(c("treatment","treatment","treatment","treatment","control","control","control","control")))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi,sampleTable,~condition)
head(assay(dds))
dds.prad <- as.matrix(assay(dds))
dds.prad.data <- data.frame(SRR5905030= apply(dds.prad[,1:2],1,sum),SRR5905031=apply(dds.prad[,3:4],1,sum),SRR5905032=apply(dds.prad[,5:6],1,sum),SRR5905033=apply(dds.prad[,7:8],1,sum))
head(dds.prad.data)
sampleTable <- data.frame(condition = c("treatment","treatment","control","control"))
dds.prad.data <- as.matrix(dds.prad.data)
dds.prad.data <- DESeqDataSetFromMatrix(countData=dds.prad.data,colData = sampleTable, design=~ condition)
dds.prad.data

keep <- rowSums(counts(dds.prad.data)) >=10
dds.prad.data <- dds.prad.data[keep,]
dds <- DESeq(dds.prad.data)
DEG_result <- results(dds)
DEG_result
result_Ordered <- DEG_result[order(DEG_result$pvalue),]
result_Ordered["NORAD",]
write.csv(as.data.frame(result_Ordered),file="DEG_result_Project1_DESeq2_Salmon.csv")
```

```{r}
Salmon_DEG <- as.data.frame(Salmon_DEG)
library(EnhancedVolcano)
EnhancedVolcano(Salmon_DEG,
    lab = rownames(Salmon_DEG),
    x = 'logFC',
    y = 'FDR',
    selectLab = c("NORAD","TTN","HPGD","SLC2A3","CCDC141","UGT2B28","ST6GALNAC1","EAF2","STEAP4","PDE3A","NNMT","PPWD1","DSN1","SLC2A1"),
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

