---
title: "ChIP-seq Analysis"
author: "Teerapon Sahwangarrom"
date: "15/12/2019"
---
```{bash}
### Make folders for this analysis
mkdir -p ChIp-seq/{data,results,scripts,software}
mkdir -p ChIP-seq/results/{bams,peaks,fastqc,motifs}
cd ChIP-seq/data

fastq-dump SRR2646261
fastq-dump SRR5238073
fastq-dump SRR653104
fastq-dump SRR3728821
fastq-dump SRR3728866

find *fastq | parallel -j 4 fastqc {} -o ../results/fastqc
REF=~/GRCh38.primary_assembly.genome.fa
bowtie2-build $REF $REF

bowtie -p 10 --best --chunkmbs 320 $REF -q anti-AR_DHT_LNCaP.fastq -S ../results/bams/anti-AR_DHT_LNCaP.sam

bowtie -p 10 --best --chunkmbs 320 $REF -q anti-AR_LNCaP.fastq -S ../results/bams/anti-AR_LNCaP.sam

bowtie -p 10 --best --chunkmbs 320 $REF -q AR_LNCaP.fastq -S ../results/bams/AR_LNCaP.sam

bowtie -p 10 --best --chunkmbs 320 $REF -q MYC_LNCaP.fastq -S ../results/bams/MYC_LNCaP.sam

bowtie -p 10 --best --chunkmbs 320 $REF -q input_LNCaP.fastq -S ../results/bams/input_LNCaP.sam
```

```{bash}
cd ../results/bams
samtools view -bS anti-AR_DHT_LNCaP.sam | samtools sort -@ 4 - -T anti-AR_DHT_LNCaP -o anti-AR_DHT_LNCaP.sorted.bam
samtools index anti-AR_DHT_LNCaP.sorted.bam

samtools view -bS anti-AR_LNCaP.sam | samtools sort -@ 4 - -T anti-AR_LNCaP -o anti-AR_LNCaP.sorted.bam
samtools index anti-AR_LNCaP.sorted.bam

samtools view -bS AR_LNCaP.sam | samtools sort -@ 4 - -T AR_LNCaP -o AR_LNCaP.sorted.bam
samtools index AR_LNCaP.sorted.bam

samtools view -bS MYC_LNCaP.sam | samtools sort -@ 4 - -T MYC_LNCaP -o MYC_LNCaP.sorted.bam
samtools index MYC_LNCaP.sorted.bam

samtools view -bS input_LNCaP.sam | samtools sort -@ 4 - -T input_LNCaP -o input_LNCaP.sorted.bam
samtools index input_LNCaP.sorted.bam

rm *sam
```

```{bash}
conda install -c bioconda macs2

macs2 callpeak -t anti-AR_DHT_LNCaP.sorted.bam -c input_LNCaP.sorted.bam -f BAM -g hs -n anti-AR_DHT_LNCaP -B -q 0.05

macs2 callpeak -t anti-AR_LNCaP.sorted.bam -c input_LNCaP.sorted.bam -f BAM -g hs -n anti-AR_LNCaP -B -q 0.05

macs2 callpeak -t AR_LNCaP.sorted.bam -c input_LNCaP.sorted.bam -f BAM -g hs -n AR_LNCaP -B -q 0.05

macs2 callpeak -t MYC_LNCaP.sorted.bam -c input_LNCaP.sorted.bam -f BAM -g hs -n MYC_LNCaP -B -q 0.05


```

```{bash}
conda install -c bioconda deeptools
bamCoverage -b anti-AR_DHT_LNCaP.sorted.bam --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 10 --extendReads 200 -o anti-AR_DHT_LNCaP.bw

bamCoverage -b anti-AR_LNCaP.sorted.bam --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 10 --extendReads 200 -o anti-AR_LNCaP.bw

bamCoverage -b AR_LNCaP.sorted.bam --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 10 --extendReads 200 -o AR_LNCaP.bw

bamCoverage -b input_LNCaP.sorted.bam --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 10 --extendReads 200 -o input_LNCaP.bw

bamCoverage -b MYC_LNCaP.sorted.bam --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 10 --extendReads 200 -o MYC_LNCaP.bw

```



