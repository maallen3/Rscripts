library(DESeq2)
library(dplyr)

cutoff=0.1
metadatafile="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/scripts/atacbammeta.csv"
coveragetablefile="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/counts/ATAC_Peaks_and_genes_122100.coverage.csv"
indir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/counts/"
annot_file = paste(indir, "hg38_refseq_genenames_and_2021_updated_ATACseqFullFiles_MUMERGE.saf", sep="")
outdir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/deseq2/peakonly/"

#read the metadata
metadata <- read.csv(metadatafile, header=TRUE, sep=",")
head(metadata)
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




#run Deseq on ATAC over peaks
ddspeaks <- DESeqDataSetFromMatrix(countData =  countdat_peakonly, colData = metadata, design = ~genotype + treatment+genotype:treatment) #use this for all samples
keep <- rowSums(counts(ddspeaks)) >= 30
ddspeaks <- ddspeaks[keep,]
DEddspeaks <- DESeq(ddspeaks)
normcounts <- as.data.frame(counts(DEddspeaks, normalize=TRUE))
write.csv(as.data.frame(normcounts), file=paste(outdir,"normalizedcounts.csv",sep=""))




sample2="thirtyseven"
sample1="fourtytwo"
contrastcol = "treatment"

respeaksD21=results(DEddspeaks, contrast=c(contrastcol,sample1,sample2))
summary(respeaksD21, alpha=cutoff)
DESeq2::plotMA(respeaksD21, main="D21 differential peaks", alpha=cutoff)
write.csv(as.data.frame(respeaksD21), file=paste(outdir,"D21res.csv",sep=""))

respeaksT21=results(DEddspeaks, contrast=list( c("treatment_fourtytwo_vs_thirtyseven", "genotypeT21.treatmentfourtytwo" )))
summary(respeaksT21)
DESeq2::plotMA(respeaksT21, main="T21 differential peaks")
write.csv(as.data.frame(respeaksT21), file=paste(outdir,"T21res.csv",sep=""))

D21resSig_peak <- subset(respeaksD21, padj < cutoff)
D21resSig_peak$GeneID <- rownames(D21resSig_peak)
D21resSig_peak_up <- subset(D21resSig_peak, log2FoldChange > 0)
D21resSig_peak_down <- subset(D21resSig_peak, log2FoldChange < 0)

head(D21resSig_peak_up[ order( D21resSig_peak_up$padj ), ])
head(D21resSig_peak_down[ order( D21resSig_peak_down$padj ), ])

beddf <- merge(as.data.frame(D21resSig_peak_up),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "D21uppeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

beddf <- merge(as.data.frame(D21resSig_peak_down),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "D21downpeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)



T21resSig_peak <- subset(respeaksT21, padj < cutoff)
T21resSig_peak$GeneID <- rownames(T21resSig_peak)
T21resSig_peak_up <- subset(T21resSig_peak, log2FoldChange > 0)
T21resSig_peak_down <- subset(T21resSig_peak, log2FoldChange < 0)

head(T21resSig_peak_up[ order( T21resSig_peak_up$padj ), ])
head(T21resSig_peak_down[ order( T21resSig_peak_down$padj ), ])

beddf <- merge(as.data.frame(T21resSig_peak_up),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "T21uppeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

beddf <- merge(as.data.frame(T21resSig_peak_down),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "T21downpeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)


