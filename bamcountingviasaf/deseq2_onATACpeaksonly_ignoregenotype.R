library(DESeq2)
library(dplyr)

cutoff=0.1
metadatafile="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/scripts/atacbammeta.csv"
coveragetablefile="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/counts/ATAC_Peaks_and_genes_122100.coverage.csv"
indir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/counts/"
annot_file = paste(indir, "hg38_refseq_genenames_and_2021_updated_ATACseqFullFiles_MUMERGE.saf", sep="")
outdir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/deseq2/peakonly_ignoregenotype/"

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
ddspeaks <- DESeqDataSetFromMatrix(countData =  countdat_peakonly, colData = metadata, design = ~genotype + treatment) #use this for all samples
keep <- rowSums(counts(ddspeaks)) >= 30
ddspeaks <- ddspeaks[keep,]
DEddspeaks <- DESeq(ddspeaks)
normcounts <- as.data.frame(counts(DEddspeaks, normalize=TRUE))
write.csv(as.data.frame(normcounts), file=paste(outdir,"normalizedcounts.csv",sep=""))




sample2="thirtyseven"
sample1="fourtytwo"
contrastcol = "treatment"

respeaks=results(DEddspeaks, contrast=c(contrastcol,sample1,sample2))
summary(respeaks, alpha=cutoff)
DESeq2::plotMA(respeaks, main="heat shock peaks", alpha=cutoff)
write.csv(as.data.frame(respeaks), file=paste(outdir,"res.csv",sep=""))

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

