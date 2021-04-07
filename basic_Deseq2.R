#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")


#set indirectory outdirectory and gtf
indir <- "/Shares/down/heatshock/analysis/PROseq/genedeseq2/"
outdir <- "/Shares/down/heatshock/analysis/PROseq/genedeseq2/"
gtf="/scratch/Shares/public/genomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"

#read the counts csv
coveragetable <- read.csv("/Shares/down/heatshock/analysis/PROseq/genedeseq2/featureCounts_attr_gene_id_feature_exon_125256.coverage.csv", row.names=1)
head(coveragetable)

#read the metadata
metadata <- read.csv("/Users/allenma/scripts/Rscripts/Rscripts/PROmeta.txt", header=TRUE, sep="\t")
head(metadata)

#make sure they are in the same order and have the same samples!!!!!
countdat <- coveragetable[,colnames(coveragetable) %in% metadata$label]

#set up the deseq object
dds <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~group)

#run deseq
DEdds <- DESeq(dds)

#Check the size factors
sizeFactors(DEdds) #this should be pretty simlilar to ratios in the millions mapped. CHECK IT!


#these two graphs show you before and after normalization of samples. In this data set that doesn't matter. In some it does.
allcounts <- as.data.frame(counts(DEdds))
allcountslong <- allcounts %>% gather(key = "sample", value = "signal")
ggplot(allcountslong, aes(x = sample, y = signal))+ geom_violin(trim = FALSE) + theme(axis.text.x = element_text(angle = 90))+scale_y_continuous(trans='log2')

normcounts <- as.data.frame(counts(DEdds, normalize=TRUE))
normcountslong <- normcounts %>% gather(key = "sample", value = "signal")
ggplot(normcountslong, aes(x = sample, y = signal))+ geom_violin(trim = FALSE) + theme(axis.text.x = element_text(angle = 90))+scale_y_continuous(trans='log2')

#plot the dispersion of the data
plotDispEsts(DEdds)

#get the results dataframes from comparing two samples 
#my samples are D21_fourtytwo D21_thirtyseven T21_fourtytwo T21_thirtyseven
#Then save the results file
sample1="D21_thirtyseven"
sample2="D21_fourtytwo"
res=results(DEdds, contrast=c("group",sample1,sample2))
fileroot = "featureCounts_attr_gene_id_feature_exon_125256"
write.csv(res, paste(outdir,fileroot,"_",sample1, "_",sample2,".results.csv", sep=""))


#look at the results
summary(res)
plotMA(res)
resSig <- subset(res, padj < 0.01)
head(res[ order( res$padj ), ])

#plot one gene
plotCounts(DEdds, gene="HSPH1", intgroup=c( "group"))
