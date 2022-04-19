library(DESeq2)
library(dplyr)

cutoff=0.1
metadatafile="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/scripts/atacbammeta.csv"
coveragetablefile="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/counts/ATAC_Peaks_and_genes_122100.coverage.csv"
indir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/counts/"
annot_file = paste(indir, "hg38_refseq_genenames_and_2021_updated_ATACseqFullFiles_MUMERGE.saf", sep="")
outdir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/deseq2/"

#read the metadata
metadata <- read.csv(metadatafile, header=TRUE, sep=",")
head(metadata)
metadata$D21orT21 <- relevel(metadata$D21orT21 , "D21")
head(metadata)
metadata$tempature <- relevel(metadata$tempature , "thirtyseven")
metadata$group <- paste(metadata$D21orT21, metadata$tempature , sep="_")
metadata$genotype <- relevel(metadata$D21orT21 , "D21")
metadata$treatment <- relevel(metadata$tempature , "thirtyseven")
head(metadata)

safdf = read.table(annot_file, sep=",", header=TRUE)
head(safdf)
tail(safdf)

coveragetable <- read.csv(coveragetablefile, row.names=1)
# this rearagnes the rows so they are in the same order as the metadata
countdat <- coveragetable %>% dplyr::select(as.vector(metadata$label))
countdat_geneonly <- countdat %>% rownames_to_column() %>% filter(!str_detect(rowname, "chr")) %>% column_to_rownames()
countdat_peakonly <- countdat %>% rownames_to_column() %>% filter(str_detect(rowname, "chr")) %>% column_to_rownames()

dim(countdat)
dim(countdat_geneonly)
dim(countdat_peakonly)


sample2="thirtyseven"
sample1="fourtytwo"
contrastcol = "treatment"

#run Deseq on ATAC over genes
ddsgene <- DESeqDataSetFromMatrix(countData =  countdat_geneonly, colData = metadata, design = ~genotype + treatment+genotype:treatment) #use this for all samples
DEddsgene <- DESeq(ddsgene)

resgeneD21=results(DEddsgene, contrast=c(contrastcol,sample1,sample2))
resgeneT21=results(DEddsgene, contrast=list( c("treatment_fourtytwo_vs_thirtyseven", "genotypeT21.treatmentfourtytwo" ))) #what changes in T21 after hs 
resgenegenotype=results(DEddsgene, contrast=c("genotype","T21","D21")) #what changes becuase of genotype
resgenediffhs=results(DEddsgene, contrast=list("genotype_T21_vs_D21", "genotypeT21.treatmentfourtytwo" )) #what changes in T21 after hs that doesn't change in D21 or 
summary(resgeneD21)
DESeq2::plotMA(resgeneD21, main="just genes")
summary(resgeneT21)
DESeq2::plotMA(resgeneT21, main="just genes")
summary(resgenegenotype)
DESeq2::plotMA(resgenegenotype, main="just genes")
summary(resgenediffhs)
DESeq2::plotMA(resgenediffhs, main="just genes")
resSig_gene_anddiff <- subset(resgenediffhs, padj < cutoff)
head(resgeneD21[ order( resgeneD21$padj ), ])
resSig_gene <- subset(resgeneD21, padj < cutoff)
resSig_gene_up <- subset(resSig_gene, log2FoldChange > 0)
resSig_gene_down <- subset(resSig_gene, log2FoldChange < 0)
#merge with saf and get bed file out
head(resSig_gene_up[ order( resSig_gene_up$padj ), ])
head(resSig_gene_down[ order( resSig_gene_down$padj ), ])

#run Deseq on ATAC over peaks
ddspeaks <- DESeqDataSetFromMatrix(countData =  countdat_peakonly, colData = metadata, design = ~genotype + treatment+genotype:treatment) #use this for all samples
DEddspeaks <- DESeq(ddspeaks)

respeaks=results(DEddspeaks, contrast=c(contrastcol,sample1,sample2))
summary(respeaks)
DESeq2::plotMA(respeaks, main="just peaks")
head(respeaks[ order( respeaks$padj ), ])
resSig_peak <- subset(respeaks, padj < cutoff)
resSig_peak$GeneID <- rownames(resSig_peak)
resSig_peak_up <- subset(resSig_peak, log2FoldChange > 0)
resSig_peak_down <- subset(resSig_peak, log2FoldChange < 0)

head(resSig_peak_up[ order( resSig_peak_up$padj ), ])
head(resSig_peak_down[ order( resSig_peak_down$padj ), ])



beddf <- merge(as.data.frame(resSig_peak_up),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "uppeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

beddf <- merge(as.data.frame(resSig_peak_down),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "downpeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)



#run Deseq on RNA-seq with 21
ddsboth <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~genotype + treatment+genotype:treatment) #use this for all samples
DEddsboth <- DESeq(ddsboth)

resboth=results(DEddsboth, contrast=c(contrastcol,sample1,sample2))
summary(resboth)
DESeq2::plotMA(resboth,  main="genes and peaks")
head(resboth[ order( resboth$padj ), ])
resSig_both <- subset(resboth, padj < cutoff)
resSig_both_up <- subset(resSig_both, log2FoldChange > 0)
resSig_both_down <- subset(resSig_both, log2FoldChange < 0)


#make genes control genes that will be used for size factors and the peaks will not be used for size factors
controlgenes <- ifelse(grepl('^chr', rownames(DEddsboth)), FALSE, TRUE)

ddsboth <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~genotype + treatment+genotype:treatment) #use this for all samples
dds_controlgenes <-estimateSizeFactors(ddsboth,controlGenes=controlgenes)
dds_controlgenes <- estimateDispersionsGeneEst(dds_controlgenes) #adds mu to dds_controlgenes <-estimateDispersionsFit(dds_controlgenesr)
dds_controlgenes<-estimateDispersionsFit(dds_controlgenes) 
dds_controlgenes<- estimateDispersionsMAP(dds_controlgenes)
dds_controlgenes <- nbinomWaldTest(dds_controlgenes)

res_controled=results(dds_controlgenes, contrast=c(contrastcol,sample1,sample2))
summary(res_controled)
DESeq2::plotMA(res_controled,  main="genes and peaks with genes controling sizefactors")
resSig_controled <- subset(res_controled, padj < cutoff)
resSig_controled$GeneID <- rownames(resSig_controled)
resSig_controled_up<- subset(resSig_controled, log2FoldChange > 0)
resSig_controled_down <- subset(resSig_controled, log2FoldChange < 0)
head(resSig_controled_up[ order( resSig_controled_up$padj ), ])
resSig_controled_up_nogenes <- as.data.frame(resSig_controled_up) %>% rownames_to_column() %>% filter(str_detect(rowname, "chr")) %>% column_to_rownames()
resSig_controled_down_nogenes <- as.data.frame(resSig_controled_down) %>% rownames_to_column() %>% filter(str_detect(rowname, "chr")) %>% column_to_rownames()



beddf <- merge(as.data.frame(resSig_controled_up_nogenes),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "uppeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

beddf <- merge(as.data.frame(resSig_controled_down_nogenes),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "downpeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)



df <- rbind(sizeFactors(DEddsgene), sizeFactors(DEddspeaks), sizeFactors(DEddsboth), sizeFactors(dds_controlgenes))
rownames(df) <- c("genes_no_control","peaks_no_control","geneandpeaks_no_control", "geneandpeaks_genes_as_sizefacto_control")
df = as.data.frame(df)
df$analysis <- rownames(df)
longdf = gather(df, sample, sizefactor, -analysis)


p<-ggplot(data=longdf, aes(x=sample, y=sizefactor, fill=analysis)) +
  geom_bar(stat="identity", position=position_dodge())+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

plotCounts(dds_controlgenes, gene="chr19:14512865-14514017", intgroup=c("tempature", "D21orT21"))


