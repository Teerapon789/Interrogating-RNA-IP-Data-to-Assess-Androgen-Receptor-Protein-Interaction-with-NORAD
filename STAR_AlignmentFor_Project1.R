---
title: "STAR_AlignmentCode"
author: "Teerapon Sahwangarrom"
date: "06/12/2019"

---
### Generating genome indexes, 
###--runThreadN option defines the number of threads to be used for genome generation
###--runMode genomeGenerate to run genome indices generation job
###--genomeDir specifies path to the directory, where the genome indices are stored
###--genomeFastaFiles specifies FASTA files with the genome reference sequences
###--sjdbGTFfile specifies the path to files with annotated transcripts in the standard GTF format
###--sjdbOverhang specifies the length of the genomic sequence around the annotated junction to be used in contructing the splice junctions database. For Illumina paired-end reads, the ideal value is 99 (ReadLength - 1)
```{bash}
./STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /rds/general/user/ts2513/home/rnaseq/refs --genomeFastaFiles /rds/general/user/ts2513/home/rnaseq/refs/GRCh38.primary_assembly.genome.fa --sjdbGTFfile /rds/general/user/ts2513/home/rnaseq/refs/gencode.v32.primary_assembly.annotation.gtf --sjdbOverhang 99
```


### Running mapping jobs
###--genomeDir specifies path to the genome directory, where genome indices are generated
###--readFilesIn names with path of the files containing the sequences to be mappled (RNA-seq FASTQ files)
###--outFileNamePrefix all output files are written in the specific directory
###--outSAMtype BAM SortedByCoordinate specifies output alignments directly in binary BAM format and it can also sort BAM files by coordinates
###--outSAMunmapped Within specifies output unmapped reads within the main SAM file
###--outSAMattributes a string of desired SAM attributes, in the order desired for the output SAM

####For trimmed-paired sequences 
```{bash}
./STAR --genomeDir /rds/general/user/ts2513/home/rnaseq/refs/ --runThreadN 6 --readFilesIn /rds/general/user/ts2513/home/rnaseq/SRR3232052.trimmed_1P.fq /rds/general/user/ts2513/home/rnaseq/SRR3232052.trimmed_2P.fq --outFileNamePrefix /rds/general/user/ts2513/home/rnaseq/results/SRR3232052_P_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
```

####For trimmed-unpaired sequence (Forward)
```{bash}
./STAR --genomeDir /rds/general/user/ts2513/home/rnaseq/refs/ --runThreadN 6 --readFilesIn /rds/general/user/ts2513/home/rnaseq/SRR5796681.trimmed_1U.fq --outFileNamePrefix /rds/general/user/ts2513/home/rnaseq/results/SRR5796681_1U_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
```

####For trimmed-unpaired sequence (Reverse)
```{bash}
./STAR --genomeDir /rds/general/user/ts2513/home/rnaseq/refs/ --runThreadN 6 --readFilesIn /rds/general/user/ts2513/home/rnaseq/SRR5796680.trimmed_2U.fq --outFileNamePrefix /rds/general/user/ts2513/home/rnaseq/results/SRR5796680_2U_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
```


##### Counting reads associated with genes
### Running featureCounts
###-a option specifies path to annotation file, GTF format
###-o option specifies path to and name of the text output (count matrix)
### /rds/general/user/ts2513/home/rnaseq/results/bams/*.out.bam, the list of all the bam files, which are collected count information for
```{bash}
./featureCounts -T 4 -a /rds/general/user/ts2513/home/rnaseq/refs/gencode.v32.primary_assembly.annotation.gtf -o /rds/general/user/ts2513/home/rnaseq/results/project_1_featurecounts.txt /rds/general/user/ts2513/home/rnaseq/results/bams/*.out.bam
```

#### Cleaning up the featureCounts matrix
### To split some columns, which are not useful to analyse the differential expression gene (DEG)
```{bash}
cut -f1,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 /rds/general/user/ts2513/home/rnaseq/results/project_1_featurecounts.txt > /rds/general/user/ts2513/home/rnaseq/results/project_1_featurecounts.Rmatrix.txt
```

### Preparing dataset before performing differential expression gene analysis (DEG)
### Merging unpaired- and paired alignment data together 
```{r}
counts.data <- read.table("/Users/teeraponsahwangarrom/Desktop/project_1_featurecounts.Rmatrix.txt",sep="\t",skip = 1, head =T)
class(counts.data)
dim(counts.data)
counts.data[1:10,1:6]

colnames(counts.data) <- c("Geneid","SRR3232052_1U", "SRR3232052_2U", "SRR3232052_P", "SRR3232053_1U", "SRR3232053_2U", "SRR3232053_P","SRR5796680_1U", "SRR5796680_2U", "SRR5796680_P", "SRR5796681_1U","SRR5796681_2U","SRR5796681_P","SRR5796682_1U","SRR5796682_2U","SRR5796682_P","SRR5796683_1U","SRR5796683_2U","SRR5796683_P")
colnames(counts.data)
counts.data[1:5,1:5]

counts.data1 <- data.frame(SRR3232052=apply(counts.data[,1:3],1,sum),SRR3232053=apply(counts.data[,4:6],1,sum),SRR5796680=apply(counts.data[,7:9],1,sum),SRR5796681=apply(counts.data[,10:12],1,sum),SRR5796682=apply(counts.data[,13:15],1,sum),SRR5796683=apply(counts.data[,16:18],1,sum))
class(counts.data1)
dim(counts.data1)
counts.data1[1:6,1:6]

```

#### Differential gene expression analysis (edgeR)
```{r}

### Creating a DGEList object, which contains the counts and also all the associated metadata
library(edgeR)
counts.data2 <- as.matrix(counts.data1)
counts.data2[1:6,1:6]
dgList_prad <- DGEList(counts=counts.data2, genes=rownames(counts.data2))
dgList_prad

### Filtering
#### removing genes that have lower than 1 cpm (counts per million) reads at least 2 samples 
countCheck <- countsPerMillion >1
head(countCheck)
keep <- which(rowSums(countCheck) >=2)
dgList_prad <- dgList_prad[keep,]
summary(cpm(dgList_prad))

### Normalisation
## To normalise RNA-seq both within and between samples
dgList_prad <- calcNormFactors(dgList_prad,method="TMM")

### Data exploration 
## To examine inter-sample relationships by creating plot based on multidimensional scaling
plotMDS(dgList_prad)

### Setting up the model
sampleType <- c("control","control","treatment","treatment","treatment","treatment")
designMat <- model.matrix(~sampleType)
designMat

### Estimating dispersions
 ## To estimate the dispersion parameter for negative binomial model 
dgList_prad <- estimateGLMCommonDisp(dgList_prad,design=designMat)
dgList_prad <- estimateGLMTrendedDisp(dgList_prad,design=designMat)
dgList_prad <- estimateGLMTagwiseDisp(dgList_prad,design=designMat)
plotBCV(dgList_prad)

### Differential expression analysis, using topTags function
fit <- glmFit(dgList_prad,designMat)
lrt <- glmLRT(fit)
write.csv(topTags(lrt,n=nrow(lrt)),file="/Users/teeraponsahwangarrom/DEGfor_Project_1_STAR.csv")

STAR_edgeR_DEG <- topTags(lrt,n=nrow(lrt))
dim(STAR_edgeR_DEG)
class(STAR_edgeR_DEG)
STAR_DEG <- as.data.frame(STAR_edgeR_DEG)
rownames(STAR_DEG) <- gsub("\\.[0-9][0-9]","",rownames(STAR_DEG))
head(STAR_DEG)
rownames(STAR_DEG) <- gsub("\\.[0-9]","",rownames(STAR_DEG))
head(STAR_DEG)

### Mapping ensembl gene ID to gene symbols
library(org.Hs.eg.db)
STAR_DEG$symbols <- mapIds(org.Hs.eg.db,keys=rownames(STAR_DEG),column="SYMBOL",keytype="ENSEMBL")
head(STAR_DEG)
### Mapping ensembl gene ID to entrez ID
STAR_DEG$entrezID <- mapIds(org.Hs.eg.db,keys=rownames(STAR_DEG),column="ENTREZID",keytype="ENSEMBL")
head(STAR_DEG)
write.csv(STAR_DEG,"/Users/teeraponsahwangarrom/STAR_DEG_edgeR.csv")
#STAR_DEG[na.position,"symbols"] <- STAR_DEG[na.position,"genes"]

```


### Volcano plot
## To visualise the result of differential expression gene analysis
```{r}
library(EnhancedVolcano)
EnhancedVolcano(STAR_DEG,lab = STAR_DEG$symbols,x="logFC",y="PValue",xlim= c(-15,15), ylim=c(0,40))
```


#### Functional enrichment analysis (GSEA)
```{r}
d <- read.csv("/Users/teeraponsahwangarrom/STAR_DEG_edgeR.csv",header=T,sep=",")
DEG.data <- filter(d,!is.na(entrezID))
geneList <- DEG.data$logFC
names(geneList) <- as.character(DEG.data$entrezID)
geneList <- sort(geneList,decreasing=T)

y <- gsePathway(geneList,nPerm=10000,pvalueCutoff=0.2,pAdjustMethod="BH",verbose=F)
res <- as.data.frame(y)
head(res)
write.csv(res,"/Users/teeraponsahwangarrom/GSEA_data_STAR_project1.csv")

library(clusterProfiler)
library(enrichplot)
dotplot(y,showCategory=30)
ridgeplot(y)

library(ggupset)
upsetplot(y)

```

