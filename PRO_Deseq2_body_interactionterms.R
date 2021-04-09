#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")


#set indirectory outdirectory and gtf
indir <- "/Shares/down/heatshock/analysis/PROseq/genedeseq3/"
outdir <- "/Shares/down/heatshock/analysis/PROseq/genedeseq3/"
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
dds <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~genotype+tempature+genotype:tempature)


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

resultsNames(DEdds)

# the tempature effect for D21 (the main effect of tempature)
resD21hs <- results(DEdds, contrast=c("tempature","thirtyseven","fourtytwo"))
head(resD21hs[ order( resD21hs$padj ), ])

# the tempature effect for T21
resT21hs <- results(DEdds, contrast=list( c("tempature_fourtytwo_vs_thirtyseven","genotypeT21.tempaturefourtytwo") ))
head(resT21hs[ order( resT21hs$padj ), ])

# the interaction term for tempature effect in T21 vs D21.
# this tests if the tempature effect is different in T21 compared to D21.
resdiff <- results(DEdds, name="genotypeT21.tempaturefourtytwo")
head(resdiff[ order( resdiff$padj ), ])
DESeq2::plotMA(resdiff)





#look at the results
summary(resT21hs)
DESeq2::plotMA(resT21hs)
resSigT21 <- subset(resT21hs, padj < cutoff)
head(res_T21[ order( resT21hs$padj ), ])
dim(resSigT21)

#look at the results
summary(resD21hs)
DESeq2::plotMA(resD21hs)
resSigD21 <- subset(resD21hs, padj < cutoff)
head(res_D21[ order( resD21hs$padj ), ])
dim(resD21hs)

#plot one gene
plotCounts(DEdds, gene="NR_110849", intgroup=c( "group"))

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
both <- df_sig %>% filter(cat=="both")
lm(both$log2FoldChange.D21~ both$log2FoldChange.T21)

df <- df %>%  mutate(max_padj = pmax(padj.T21,padj.D21))
df <- df %>%  mutate(min_padj = pmin(padj.T21,padj.D21))
df <- df %>% mutate(cat = case_when(max_padj <cutoff~"both", padj.T21 <cutoff ~"T21", padj.D21 <cutoff ~"D21", max_padj >cutoff~"N.S"))
df$lfcvalsdiff = df$log2FoldChange.D21-df$log2FoldChange.T21


df_sig %>% group_by(cat) %>% summarize(count=n()) 

df %>% group_by(cat) %>% summarize(count=n()) 


ggplot(both, aes(x=log2FoldChange.D21, y=log2FoldChange.T21, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)+stat_smooth(method = "lm", col = "red")

ggplot(df_sig, aes(x=log2FoldChange.D21, y=log2FoldChange.T21, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)+stat_smooth(method = "lm", col = "red")

hist(df_sig$padj.T21, breaks=100)
hist(df_sig$padj.D21, breaks=100)

ggplot(df, aes(x=baseMean.D21, y=log2FoldChange.D21, color=cat))+geom_point()+scale_x_continuous(trans="log10")+facet_wrap(~cat)
ggplot(df, aes(x=baseMean.T21, y=log2FoldChange.T21, color=cat))+geom_point()+scale_x_continuous(trans="log10")+facet_wrap(~cat)


normcounts_group <- normcounts[rownames(normcounts) %in% noquote(both$Row.names),]
normcounts_group <- as.matrix(normcounts_group)
normcounts_grouplog <- na.omit(normcounts_group)
normcounts_grouplog <-log(normcounts_grouplog)
normcounts_grouplog <- na.omit(normcounts_grouplog)

#pheatmap(normcounts_grouplog, cluster_cols = F)
#I can't get this to work!


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

normcounts_group_rowscaled <- t(apply(normcounts_group, 1, cal_z_score))

pheatmap(normcounts_group_rowscaled, cluster_cols = F)

lfcvals <- both %>% select("Row.names", "log2FoldChange.D21", "log2FoldChange.T21")
lfcvals$diff = lfcvals$log2FoldChange.D21-lfcvals$log2FoldChange.T21
head(lfcvals[ order( lfcvals$diff ), ])
head(lfcvals[ order( lfcvals$diff ,decreasing = TRUE), ])
lfcvalsdown <- lfcvals %>% filter(lfcvals$log2FoldChange.D21<0) %>% arrange(desc(diff))
lfcvalsup <- lfcvals %>% filter(lfcvals$log2FoldChange.D21>0) %>% arrange(diff)

hist(lfcvals$diff, breaks=100)  
hist(lfcvalsup$diff, breaks=100)
hist(lfcvalsdown$diff, breaks=100)



onegene="NM_145202"
plotCounts(DEdds, gene=onegene, intgroup=c( "group"))
