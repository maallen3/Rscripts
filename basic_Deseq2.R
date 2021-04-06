
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

allcounts <- as.data.frame(counts(DEdds))
allcountslong <- allcounts %>% gather(key = "sample", value = "signal")
ggplot(allcountslong, aes(x = sample, y = signal))+ geom_violin(trim = FALSE) + theme(axis.text.x = element_text(angle = 90))+scale_y_continuous(trans='log2')

normcounts <- as.data.frame(counts(DEdds, normalize=TRUE))
normcountslong <- normcounts %>% gather(key = "sample", value = "signal")
ggplot(normcountslong, aes(x = sample, y = signal))+ geom_violin(trim = FALSE) + theme(axis.text.x = element_text(angle = 90))+scale_y_continuous(trans='log2')

#plot the dispersion of the data
plotDispEsts(DEdds)

#get the results dataframes (which was not )
resD21_37vs42 <- results(DEdds, contrast=c("group","D21_thirtyseven","D21_fourtytwo"))
resT21_37vs42 <- results(DEdds, contrast=c("group","T21_thirtyseven","T21_fourtytwo"))
res37_D21vsT21 <- results(DEdds, contrast=c("group","D21_thirtyseven","T21_thirtyseven"))
res42_D21vsT21 <- results(DEdds, contrast=c("group","D21_fourtytwo","T21_fourtytwo"))


#look at and output the results
res = resD21_37vs42
summary(res)
plotMA(res)
resSig <- subset(res, padj < 0.01)
head(res[ order( res$padj ), ])

#plot one gene
plotCounts(DEdds, gene="HSPH1", intgroup=c( "group"))
