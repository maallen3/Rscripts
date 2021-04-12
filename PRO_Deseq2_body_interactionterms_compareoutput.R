#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")


#set indirectory outdirectory and gtf
indir <- "/Shares/down/heatshock/analysis/PROseq/genedeseq3/"
outdir <- "/Shares/down/heatshock/analysis/PROseq/genedeseq3/"
cutoff = 0.001

#read the counts csv
fileroot = "featureCounts_body_153616"
coveragetable <- read.csv(paste(indir,fileroot,".coverage.csv", sep=""), row.names=1)

head(coveragetable)

#read the metadata
metadata <- read.csv("/Users/allenma/scripts/Rscripts/Rscripts/PROmeta.txt", header=TRUE, sep="\t")
head(metadata)
#if you are using the genotype:tempature interaction term its really imporant that you make sure the reference level is correct!
metadata$genotype <- relevel(metadata$genotype , "D21")
head(metadata)
metadata$tempature <- relevel(metadata$tempature , "thirtyseven")
head(metadata)
metadata$group <- relevel(metadata$group , "D21_thirtyseven")
head(metadata)

#make sure they are in the same order and have the same samples!!!!!
countdat <- coveragetable %>% select(as.vector(metadata$label))

#set up the deseq object
dds0 <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~genotype+tempature+genotype:tempature)
dds1 <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~genotype+tempature)
dds2 <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~tempature)
dds3 <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~group)
dds0_r <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~replicate+genotype+tempature+genotype:tempature)
dds1_r <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~replicate+genotype+tempature)
dds2_r <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~replicate+tempature)
dds3_r <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~replicate+group)



#run deseq
DEdds0 <- DESeq(dds0)
DEdds1 <- DESeq(dds1)
DEdds2 <- DESeq(dds2)
DEdds3 <- DESeq(dds3)
DEdds0_r <- DESeq(dds0_r)
DEdds1_r <- DESeq(dds1_r)
DEdds2_r <- DESeq(dds2_r)
DEdds3_r <- DESeq(dds3_r)

#plot the dispersion of the data
plotDispEsts(DEdds0)
plotDispEsts(DEdds1)
plotDispEsts(DEdds2)
plotDispEsts(DEdds3)
plotDispEsts(DEdds0_r)
plotDispEsts(DEdds1_r)
plotDispEsts(DEdds2_r)
plotDispEsts(DEdds3_r)

#get the results dataframes from comparing two samples 
#my samples are D21_fourtytwo D21_thirtyseven T21_fourtytwo T21_thirtyseven
#Then save the results file

resultsNames(DEdds0)
resultsNames(DEdds1)
resultsNames(DEdds2)
resultsNames(DEdds3)
resultsNames(DEdds0_r)
resultsNames(DEdds1_r)
resultsNames(DEdds2_r)
resultsNames(DEdds3_r)

# the tempature effect for D21 (the main effect of tempature)
resD21hs0 <- results(DEdds0, contrast=c("tempature","thirtyseven","fourtytwo"))
resD21hs0_r <- results(DEdds0_r, contrast=c("tempature","thirtyseven","fourtytwo"))


# the tempature effect for T21
resT21hs0 <- results(DEdds0, contrast=list( c("tempature_fourtytwo_vs_thirtyseven","genotypeT21.tempaturefourtytwo") ))
resT21hs0_r <- results(DEdds0_r, contrast=list( c("tempature_fourtytwo_vs_thirtyseven","genotypeT21.tempaturefourtytwo") ))


# the interaction term for tempature effect in T21 vs D21.
# this tests if the tempature effect is different in T21 compared to D21.
resdiff_D21_T210 <- results(DEdds0, name="genotypeT21.tempaturefourtytwo")
resdiff_D21_T210_r <- results(DEdds0_r, name="genotypeT21.tempaturefourtytwo")

# the tempature effect, with knowledge of genotype """but not genotype specific""" (the main effect of tempature)
reshs1 <- results(DEdds1, contrast=c("tempature","thirtyseven","fourtytwo"))
reshs1_r <- results(DEdds1_r, contrast=c("tempature","thirtyseven","fourtytwo"))


# the tempature effect for tempature, without knowledge of genotype (the main effect of tempature)
reshs2 <- results(DEdds2, contrast=c("tempature","thirtyseven","fourtytwo"))
reshs2_r <- results(DEdds2_r, contrast=c("tempature","thirtyseven","fourtytwo"))


# the tempature effect for D21 
resD21hs3 <- results(DEdds3, contrast=c("group","D21_thirtyseven","D21_fourtytwo"))
resD21hs3_r <- results(DEdds3_r, contrast=c("group","D21_thirtyseven","D21_fourtytwo"))


# the tempature effect for T21
resT21hs3 <- results(DEdds3, contrast=c("group","T21_thirtyseven", "T21_fourtytwo") )
resT21hs3_r <- results(DEdds3_r, contrast=c("group","T21_thirtyseven", "T21_fourtytwo") )



#These are summarys of all the effects we have created:

#the tempature effect 
DESeq2::plotMA(reshs1)
summary(reshs1)
DESeq2::plotMA(reshs1_r)
summary(reshs1_r)
#the tempature effect  even if genotype changes the start level
DESeq2::plotMA(reshs2)
summary(reshs2)
DESeq2::plotMA(reshs2_r)
summary(reshs2_r)
# the tempature effect for D21
DESeq2::plotMA(resD21hs0)
summary(resD21hs0)
DESeq2::plotMA(resD21hs3)
summary(resD21hs3)
# the tempature effect for T21
DESeq2::plotMA(resT21hs0)
summary(resT21hs0)
DESeq2::plotMA(resT21hs3)
summary(resT21hs3)
#major differences in how D21 and T21 change the temp effect
DESeq2::plotMA(resdiff_D21_T210)
summary(resdiff_D21_T210)

#combine the calls from D21 and T21 heat shock differences to determine
#the minor difference betten how D21 and T21 change the temp effect
mergeres<- function(res1, res2, name1, name2){
  res1 = as.data.frame(res1)
  res2 = as.data.frame(res2)
  df <- merge(res1, res2, by.x=0, by.y=0, suffixes=c(paste(".",name1, sep=""), paste(".",name2, sep="")))
}

transpadj <- function(x){
  (-log10(x))^(1/6)}

df <- mergeres(resD21hs3, resT21hs3, "D21", "T21")
df <- df %>% mutate(max_padj = pmax(padj.T21,padj.D21)) %>%  mutate(min_padj = pmin(padj.T21,padj.D21)) %>% mutate(cat = case_when(max_padj <cutoff~"both", padj.T21 <cutoff ~"T21", padj.D21 <cutoff ~"D21", max_padj >cutoff~"N.S")) %>% mutate(padj.D21.trans = transpadj(padj.D21)) %>%mutate(padj.T21.trans = transpadj(padj.T21))
df$lfcvalsdiff = df$log2FoldChange.D21-df$log2FoldChange.T21
df_sig_either <- df %>% filter(padj.T21<cutoff | padj.D21<cutoff) 
both <- df_sig_either %>% filter(cat=="both")


#MApolt with D21 fc 
ggplot(df, aes(x=baseMean.D21, y=log2FoldChange.D21, color=cat))+geom_point()+scale_x_continuous(trans="log10")
#MApolt with T21 fc 
ggplot(df, aes(x=baseMean.T21, y=log2FoldChange.T21, color=cat))+geom_point()+scale_x_continuous(trans="log10")


#How many genes change via HS for D21 or T21 or both?
df %>% group_by(cat) %>% summarize(count=n()) 

#Why?
hist(df_sig$padj.T21, breaks=100)
hist(df_sig$padj.D21, breaks=100)


#plot of the fold change in D21 vs. T21 for genes that change in both, or just one of the two
ggplot(df_sig_either, aes(x=log2FoldChange.D21, y=log2FoldChange.T21, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)


#are the changes seen in only one sample just below the cuttoff in the other sample
#dotted line is a padj value of 0.05
#our cut off is 0.01
ggplot(df_sig_either, aes(x=padj.D21.trans, y=padj.T21.trans, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)+ggtitle("adj-pval transformed by -log10(x)^(1/6)") +
  xlab("padj D21") + ylab("padj T21")+ geom_vline(xintercept = transpadj(0.05), linetype="dotted", color = "blue", size=1.5)+ geom_hline(yintercept = transpadj(0.05), linetype="dotted", color = "green", size=1.5)


#plot of the fold change in D21 vs. T21 for genes that change in both, or just one of the two
ggplot(df_sig_either, aes(x=log2FoldChange.D21, y=log2FoldChange.T21, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)

#Does it look like the Genes going up in both samples are going up more in T21?
lm(both$log2FoldChange.D21~ both$log2FoldChange.T21)
ggplot(both, aes(x=log2FoldChange.D21, y=log2FoldChange.T21, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)+stat_smooth(method = "lm", col = "red")
#YEP!

#Lets check another way. Histogram of the the differnce in log fold changes
lfcvals <- both %>% select("Row.names", "log2FoldChange.D21", "log2FoldChange.T21")
lfcvals$diff = lfcvals$log2FoldChange.D21-lfcvals$log2FoldChange.T21
head(lfcvals[ order( lfcvals$diff ), ])
head(lfcvals[ order( lfcvals$diff ,decreasing = TRUE), ])
lfcvalsdown <- lfcvals %>% filter(lfcvals$log2FoldChange.D21<0) %>% arrange(desc(diff))
lfcvalsup <- lfcvals %>% filter(lfcvals$log2FoldChange.D21>0) %>% arrange(diff)
hist(lfcvals$diff, breaks=100)  
hist(lfcvalsup$diff, breaks=100)
hist(lfcvalsdown$diff, breaks=100)
#T21 increase reaction to heat shock


#Lets exlore how using interation terms changed the calls---

# notice that includeing the replicate helps Deseq to call more genes as differental
DESeq2::plotMA(reshs1)
summary(reshs1)
DESeq2::plotMA(reshs1_r)
summary(reshs1_r)
# notice that inlcudeing the genotype information decreased the total number of genes called as an effect of tempature. This is because there are less samples per group. 
#example 1
DESeq2::plotMA(reshs1)
summary(reshs1)
DESeq2::plotMA(reshs2)
summary(reshs2)
#example 2 (with rep info)
DESeq2::plotMA(reshs1_r)
summary(reshs1_r)
DESeq2::plotMA(reshs2_r)
summary(reshs2_r)

#example of a gene differenal in one but not in other
df_diffinteractions <- mergeres(reshs1, reshs2, "temp", "temp_know_geno")
df_diffinteractions <- df_diffinteractions %>% mutate(max_padj = pmax(padj.temp,padj.temp_know_geno)) %>%
  mutate(min_padj = pmin(padj.temp,padj.temp_know_geno)) %>% 
  mutate(cat = case_when(max_padj <cutoff~"both", padj.temp <cutoff ~"temp", padj.temp_know_geno <cutoff ~"temp_know_geno", max_padj >cutoff~"N.S")) %>%
  mutate(padj.temp.trans = transpadj(padj.temp)) %>%
  mutate(padj.temp_know_geno.trans = transpadj(padj.temp_know_geno))
df_diffinteractions$lfcvalsdiff = df_diffinteractions$log2FoldChange.temp-df_diffinteractions$log2FoldChange.temp_know_geno
df_diffinteractions_sig_either <- df_diffinteractions %>% filter(padj.temp<cutoff | padj.temp_know_geno<cutoff) 
df_diffinteractions_both <- df_diffinteractions_sig_either %>% filter(cat=="both")

#number of genes in each group
df_diffinteractions %>% group_by(cat) %>% summarize(count=n()) 

#genes only seen when Deseq2 is aware of the genotype
df_temp_know_geno_only <- df_diffinteractions %>% filter(cat=="temp_know_geno") %>% arrange(desc(padj.temp_know_geno))
df_temp <- df_diffinteractions %>% filter(cat=="temp") %>% arrange(desc(padj.temp_know_geno)) 

#plot all genes only in temp alone
rld <- rlog(DEdds2, blind=TRUE)
rld_mat <- assay(rld)
cluster_rlog <- rld_mat[df_temp$Row.names, ]
normcounts <- as.data.frame(counts(DEdds2, normalize=TRUE))
meta <-as.data.frame(cbind(paste(DEdds2$label), paste(DEdds2$genotype),paste(DEdds2$group),paste(DEdds2$replicate), paste(DEdds2$tempature)))
colnames(meta)<-c("label","genotype", "group", "replicate", "tempature")
row.names(meta)<-meta$label
meta$tempature <- relevel(meta$tempature, ref="thirtyseven")
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "tempature", col="genotype")
cluster_groups <- clusters$df
for (i in unique(cluster_groups$cluster)){
  group <- clusters$df %>% filter(cluster == i)
  normcounts_group <- normcounts[rownames(normcounts) %in% noquote(group$genes),]
  normcounts_group <- as.matrix(normcounts_group)
  title = paste0("Design only with temp cluster", i, sep=" ")
  heatmap(normcounts_group,Colv = NA, main=title)}




#plot all genes only in temp knows genotype
rld <- rlog(DEdds1, blind=TRUE)
rld_mat <- assay(rld)
cluster_rlog <- rld_mat[df_temp_know_geno_only$Row.names, ]
normcounts <- as.data.frame(counts(DEdds1, normalize=TRUE))
meta <-as.data.frame(cbind(paste(DEdds1$label), paste(DEdds1$genotype),paste(DEdds1$group),paste(DEdds1$replicate), paste(DEdds1$tempature)))
colnames(meta)<-c("label","genotype", "group", "replicate", "tempature")
row.names(meta)<-meta$label
meta$tempature <- relevel(meta$tempature, ref="thirtyseven")
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "tempature", col="genotype")
cluster_groups <- clusters$df
for (i in unique(cluster_groups$cluster)){
  group <- clusters$df %>% filter(cluster == i)
  normcounts_group <- normcounts[rownames(normcounts) %in% noquote(group$genes),]
  normcounts_group <- as.matrix(normcounts_group)
  title = paste0("Design with temp and genotype cluster", i, sep=" ")
  heatmap(normcounts_group,Colv = NA, main=title)}




