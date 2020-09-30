---
title: "STAR_RIPseq_Analysis"
author: "Teerapon Sahwangarrom"
date: "16/02/2020"
output: html_document
---

```{r}
counts.data <- read.table("/Users/teeraponsahwangarrom/Desktop/Teerapon_project1/project_1_featurecounts.Rmatrix.txt",sep="\t",skip = 1, head =T,row.names = 1)
class(counts.data)
dim(counts.data)
counts.data[1:10,1:6]
colnames(counts.data) <- c("SRR3232052_1U", "SRR3232052_2U", "SRR3232052_P", "SRR3232053_1U", "SRR3232053_2U", "SRR3232053_P","SRR5796680_1U", "SRR5796680_2U", "SRR5796680_P", "SRR5796681_1U","SRR5796681_2U","SRR5796681_P","SRR5796682_1U","SRR5796682_2U","SRR5796682_P","SRR5796683_1U","SRR5796683_2U","SRR5796683_P")
counts.data[1:5,1:5]
counts.data1 <- data.frame(SRR5796680=apply(counts.data[,7:9],1,sum),SRR5796681=apply(counts.data[,10:12],1,sum),SRR5796682=apply(counts.data[,13:15],1,sum),SRR5796683=apply(counts.data[,16:18],1,sum))
class(counts.data1)
dim(counts.data1)
counts.data1[1:4,1:4]

```

```{r}
library(edgeR)
counts.data2 <- as.matrix(counts.data1)
counts.data2[1:4,1:4]
group <- factor(c("treated_AR","treated_AR","treated_IgG","treated_IgG"))
design <- model.matrix( ~0+group)
head(design)
design
dim(counts.data2)
counts.data2 <- counts.data2[ rowSums(cpm(counts.data2)>100) >= 2, ]
dim(counts.data2)

dgList_prad <- DGEList(counts=counts.data2, group = group, lib.size = colSums(counts.data2))
dgList_prad <- calcNormFactors(dgList_prad,method="TMM")
##dgList_prad$samples$group <- c(1,1,2,2,3,3)
plotMDS(dgList_prad, method = "bcv", col=as.numeric(dgList_prad$samples$group))
##dgList_prad$samples$group <- c("treated_AR","treated_AR","treated_IgG","treated_IgG")
legend("topleft", as.character(unique(dgList_prad$samples$group)),col = 1:2,pch = 20)
dgList_prad
```

```{r}
dgList_prad <- estimateGLMCommonDisp(dgList_prad,design)
dgList_prad <- estimateGLMTrendedDisp(dgList_prad,design)
dgList_prad <- estimateGLMTagwiseDisp(dgList_prad,design)
plotBCV(dgList_prad)

fit <- glmFit(dgList_prad,design)
contrast_treatedAR_vs_treatedIgG <- glmLRT(fit,contrast=makeContrasts(grouptreated_AR-grouptreated_IgG, levels = design))
STAR_DEG <- topTags(contrast_treatedAR_vs_treatedIgG,n = nrow(dgList_prad$counts))
write.csv(STAR_DEG,file="/Users/teeraponsahwangarrom/DEGfor_Project_1_STAR.csv")
dim(STAR_DEG)
class(STAR_DEG)
STAR_DEG <- as.data.frame(STAR_DEG)
rownames(STAR_DEG) <- gsub("\\.[0-9][0-9]","",rownames(STAR_DEG))
head(STAR_DEG)
rownames(STAR_DEG) <- gsub("\\.[0-9]","",rownames(STAR_DEG))
head(STAR_DEG)

library(org.Hs.eg.db)
STAR_DEG$symbols <- mapIds(org.Hs.eg.db,keys=rownames(STAR_DEG),column="SYMBOL",keytype="ENSEMBL")
head(STAR_DEG)
STAR_DEG$entrezID <- mapIds(org.Hs.eg.db,keys=rownames(STAR_DEG),column="ENTREZID",keytype="ENSEMBL")
head(STAR_DEG)
STAR_DEG["ENSG00000260032",]
write.csv(STAR_DEG,"/Users/teeraponsahwangarrom/RIPSTAR_DEG_edgeR.csv")

library(ggplot2)
library(ggrepel)
threshold_OE <- STAR_DEG$FDR < 0.05
STAR_DEG$threshold <- threshold_OE
DEGgenes_ordered <- STAR_DEG[order(STAR_DEG$FDR),]
DEGgenes_ordered$genelabels <- ""
DEGgenes_ordered$genelabels[2011] <- DEGgenes_ordered$symbols[which(rownames(DEGgenes_ordered)== "ENSG00000260032")]
DEGgenes_ordered
ggplot(DEGgenes_ordered) + geom_point(aes(x=logFC, y= -log10(FDR), colour= threshold)) + geom_text_repel(aes(x=logFC,y=-log10(FDR),label=ifelse(DEGgenes_ordered$genelabels==T,as.character(DEGgenes_ordered$symbols),""))) + ggtitle("Differential expression genes") + xlab("log2 fold change") + ylab("-log10 FDR") + theme(legend.position= "none", plot.title=element_text(size=rel(1.5), hjust= 0.5), axis.title=element_text(size=rel(1.25)))




```

```{r}
library(EnhancedVolcano)
EnhancedVolcano(DEGgenes_ordered,
    lab = DEGgenes_ordered$symbols,
    x = 'logFC',
    y = 'FDR',
    selectLab = c("ATP6","ND2","ND5","COX2","ND1","WWP1","EIF3B","ENTPD6","IMPDH1","PRKCD","NORAD","C6orf62","RTL8A","TUBB","KDELR2"),
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

