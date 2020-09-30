---
title: "RNAseq_STAR_Project1_hg38"
author: "Teerapon Sahwangarrom"
date: "11/02/2020"
---

```{r}
counts.data <- read.table("/Users/teeraponsahwangarrom/Desktop/project_1_2featurecounts.Rmatrix.txt",sep="\t",skip = 1, head =T, row.names = 1)
class(counts.data)
dim(counts.data)
counts.data[1:10,1:6]
```

```{r}
colnames(counts.data) <- c("SRR5905030_1U", "SRR5905030_2U", "SRR5905030_P", "SRR5905031_1U", "SRR5905031_2U", "SRR5905031_P","SRR5905032_1U", "SRR5905032_2U", "SRR5905032_P", "SRR5905033_1U","SRR5905033_2U","SRR5905033_P")
colnames(counts.data)
counts.data[1:5,1:5]
dim(counts.data)
```

```{r}
counts.data1 <- data.frame(SRR5905030=apply(counts.data[,1:3],1,sum),SRR5905031=apply(counts.data[,4:6],1,sum),SRR5905032=apply(counts.data[,7:9],1,sum),SRR5905033=apply(counts.data[,10:12],1,sum))
class(counts.data1)
dim(counts.data1)
counts.data1[1:4,1:4]
NORAD_data <- data.frame(Treated_AR=apply(counts.data1[,1:2],1,sum),Untreated_AR=apply(counts.data1[,3:4],1,sum))
NORAD_data
NORAD_data1 <- NORAD_data["ENSG00000260032.2",]
NORAD_data1
#write.csv(NORAD_data,"/Users/teeraponsahwangarrom/Desktop/NORAD_data.csv")

```

```{r}
library(edgeR)
counts.data2 <- as.matrix(counts.data1)
counts.data2[1:4,1:4]
group <- factor(c("treated_AR","treated_AR","untreated_AR","untreated_AR"))
design <- model.matrix( ~0+group)
head(design)
design


dim(counts.data2)
counts.data2 <- counts.data2[ rowSums(cpm(counts.data2)>100) >= 2, ]
dim(counts.data2)

dgList_prad <- DGEList(counts=counts.data2, group = group, lib.size = colSums(counts.data2))
dgList_prad <- calcNormFactors(dgList_prad,method="TMM")
plotMDS(dgList_prad, method = "bcv", col=as.numeric(dgList_prad$samples$group))
legend("topleft", as.character(unique(dgList_prad$samples$group)),col = 1:3,pch = 20)

dgList_prad <- estimateGLMCommonDisp(dgList_prad,design)
dgList_prad <- estimateGLMTrendedDisp(dgList_prad,design)
dgList_prad <- estimateGLMTagwiseDisp(dgList_prad,design)

plotBCV(dgList_prad)

fit <- glmFit(dgList_prad,design)
contrast_untrestedAR_vs_treatedIgG <- glmLRT(fit,contrast=makeContrasts(grouptreated_AR-groupuntreated_AR, levels = design))
STAR_DEG <- topTags(contrast_untrestedAR_vs_treatedIgG,n = nrow(dgList_prad$counts))
write.csv(STAR_DEG,file="/Users/teeraponsahwangarrom/DEGfor_Project_1_STAR.csv")

dim(STAR_DEG)

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
STAR_DEG[1:10,]
write.csv(STAR_DEG,"/Users/teeraponsahwangarrom/RNASTAR_DEG_edgeR.csv")

```

```{r}


library(EnhancedVolcano)
EnhancedVolcano(STAR_DEG,
    lab = STAR_DEG$symbols,
    x = 'logFC',
    y = 'FDR',
    selectLab = c("NORAD","WSB1","HPGD","NDRG1","FKBP5","UGT2B15","OTULINL","MBOAT2","IGF1R","LIFR","PHC3","EIF3E"),
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

