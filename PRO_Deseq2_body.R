#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")


#set indirectory outdirectory and gtf
indir <- "/Shares/down/heatshock/analysis/PROseq/genedeseq2/"
outdir <- "/Shares/down/heatshock/analysis/PROseq/genedeseq2/"
cutoff = 0.001

#read the counts csv
coveragetable <- read.csv("/Shares/down/heatshock/analysis/PROseq/genedeseq2/featureCounts_body_122529.coverage.csv", row.names=1)
fileroot = "featureCounts_body_122529"
head(coveragetable)

#read the metadata
metadata <- read.csv("/Users/allenma/scripts/Rscripts/Rscripts/PROmeta.txt", header=TRUE, sep="\t")
head(metadata)
metadata$genotype <- relevel(metadata$genotype , "D21")
head(metadata)
metadata$tempature <- relevel(metadata$tempature , "thirtyseven")
head(metadata)

#make sure they are in the same order and have the same samples!!!!!
countdat <- coveragetable %>% select(as.vector(metadata$label))

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



sample1="T21_thirtyseven"
sample2="T21_fourtytwo"
res_T21=results(DEdds, contrast=c("group",sample1,sample2))
write.csv(res_T21, paste(outdir,fileroot,"_",sample1, "_",sample2,".results.csv", sep=""))

sample1="D21_thirtyseven"
sample2="D21_fourtytwo"
res_D21=results(DEdds, contrast=c("group",sample1,sample2))
write.csv(res_D21, paste(outdir,fileroot,"_",sample1, "_",sample2,".results.csv", sep=""))


#look at the results
summary(res)
DESeq2::plotMA(res)
resSig <- subset(res, padj < cutoff)
head(res[ order( res$padj ), ])
dim(resSig)

#look at the results
summary(res_T21)
DESeq2::plotMA(res_T21)
resSigT21 <- subset(res_T21, padj < cutoff)
head(res_T21[ order( res_T21$padj ), ])
dim(resSigT21)

#look at the results
summary(res_D21)
DESeq2::plotMA(res_D21)
resSigD21 <- subset(res_D21, padj < cutoff)
head(res_D21[ order( res_D21$padj ), ])
dim(resSigD21)

#plot one gene
plotCounts(DEdds, gene="NM_020161", intgroup=c( "group"))

plotlogfc<- function(res1, res2, name1, name2){
  res1 = as.data.frame(res1)
  res2 = as.data.frame(res2)
  df <- merge(res1, res2, by.x=0, by.y=0, suffixes=c(paste(".",name1, sep=""), paste(".",name2, sep="")))
  }

df <- plotlogfc(res_D21, res_T21, "D21", "T21")
df_sig <- df %>% filter(padj.T21<cutoff | padj.D21<cutoff)
df_sig <- df_sig %>%  mutate(max_padj = pmax(padj.T21,padj.D21))
df_sig <- df_sig %>%  mutate(min_padj = pmin(padj.T21,padj.D21))
df_sig <- df_sig %>% mutate(cat = case_when(max_padj <cutoff~"both", padj.T21 <cutoff ~"T21", padj.D21 <cutoff ~"D21"))

df <- df %>%  mutate(max_padj = pmax(padj.T21,padj.D21))
df <- df %>%  mutate(min_padj = pmin(padj.T21,padj.D21))
df <- df %>% mutate(cat = case_when(max_padj <cutoff~"both", padj.T21 <cutoff ~"T21", padj.D21 <cutoff ~"D21", max_padj >cutoff~"N.S"))


df_sig %>% group_by(cat) %>% summarize(count=n()) 

df %>% group_by(cat) %>% summarize(count=n()) 


ggplot(df_sig, aes(x=log2FoldChange.D21, y=log2FoldChange.T21, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)

ggplot(df_sig, aes(x=log2FoldChange.D21, y=log2FoldChange.T21, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)+facet_wrap(~cat)


hist(df_sig$padj.T21, breaks=100)
hist(df_sig$padj.D21, breaks=100)

ggplot(df, aes(x=baseMean.D21, y=log2FoldChange.D21, color=cat))+geom_point()+scale_x_continuous(trans="log10")+facet_wrap(~cat)
ggplot(df, aes(x=baseMean.T21, y=log2FoldChange.T21, color=cat))+geom_point()+scale_x_continuous(trans="log10")+facet_wrap(~cat)


normcounts_group <- normcounts[rownames(normcounts) %in% noquote(row.names(resSig)),]
normcounts_group <- as.matrix(normcounts_group)


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

normcounts_group_rowscaled <- t(apply(normcounts_group, 1, cal_z_score))

pheatmap(normcounts_group_rowscaled, cluster_cols = F)

