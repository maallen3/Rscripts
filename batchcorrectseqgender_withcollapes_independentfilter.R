#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("sva")
#BiocManager::install("DEGreport")
#BiocManager::install("vsn")

library('dplyr')
library('tidyverse')
library('pheatmap')
library("vsn")
library('sva')
library('limma')
library("DESeq2")
library("DEGreport")

#function for making outdirectories
mkdirs <- function(fp) {
  if(!file.exists(fp)) {
    mkdirs(dirname(fp))
    dir.create(fp)
  }
}

#set indirectory outdirectory 
masterdir="/Users/maryallen/Desktop/Quests/CUBoulder/projects/2021/hoeffer/Hoeffer_results_aug17_2021/"
scriptdir=paste0(masterdir, "scripts/", sep="")
indir <- paste0(masterdir, "input/", sep="")
outdir=paste0(masterdir, "August2021/", sep="")
mkdirs(outdir)
subdir="batchcorrectseqgender_collapse_independentFiltering/"
outdir=paste0(outdir, subdir, sep="")
mkdirs(outdir)

#read in annotation files
genes <- read.csv(paste(indir, "featureCounts_attr_gene_name_feature_exon_125500.annotation.csv", sep=""), row.names=1)
#chromosome  of gene
genes<- genes %>% separate(Chr, c("Chr_first", "Chr_rest"), extra = "merge", fill = "right")
#first position of gene
genes<- genes %>% separate(Start, c("Start_first", "Start_rest"), extra = "merge", fill = "right")
genes<- genes %>% separate(Strand, c("Strand_first", "Strand_rest"), extra = "merge", fill = "right", sep=";")

#last position of gene
get.last.val <- function(name) {
  lapply(ifelse(grepl(";",as.character(name)),strsplit(as.character(name),";"),as.character(name)),`[[`,1)
}
genes <- genes %>% mutate("last_END"=get.last.val(End))

#make a list of the triplicated genes
genes_onlyDp16tri <- genes %>% filter(Chr_first=="chr16")  %>% filter(Start_first>=75540513) %>% filter(last_END<=97962622)
genes_onlyDp16trinotRcan <-genes_onlyDp16tri %>% filter(GeneID!="Rcan1")
genes$score = "."

#Make a bed file likst file from the annotation.
generegion <-select(genes,Chr_first, Start_first, last_END, GeneID,score, Strand_first)


#read the counts csv
#paired end and single end samples had to be counted separate becuase feature counts has a paired end flag
#paired end samples
coveragetable1 <- read.csv(paste(indir, "featureCounts_attr_gene_name_feature_exon_125500.coverage.csv", sep=""), row.names=1)
#single end samples
coveragetable2 <- read.csv(paste(indir, "featureCounts_attr_gene_name_feature_exon_142156.coverage.csv", sep=""), row.names=1)
coveragetable <-cbind(coveragetable1, coveragetable2)
#this is all the samples counts
head(coveragetable)
dim(coveragetable)
#read the metadata file
metadata <- read.csv(paste0(indir, "metadata.csv"), header=TRUE)
filtermetadata = metadata[metadata$failed=="n",] #filter out the samples that failed becuase of adapter content or ribosomal content
#remove samples from countdata and put them in order
countdat <- coveragetable %>% select(as.vector(filtermetadata$specificlabel))

counts <- as.matrix(countdat)
batch1 <-filtermetadata$singleorpaired
batch2 <-filtermetadata$gender



#remove the genes with 0 variation in counts, they are not useful
m_countdat=counts[apply(counts,1,var)>0,]
## And here we use the Bayesian approach for batch correction included in the sva library
sva_corrected <- ComBat_seq(m_countdat, batch=filtermetadata$singleorpaired)
sva_corrected <- ComBat_seq(sva_corrected, batch=filtermetadata$gender)

## Basic PCA plot, this is too messy to interpret with random data but
## shows how to plot the results of a PCA on your corrected data.
batches <- filtermetadata$singleorpaired
diagnostic <- TRUE
if(diagnostic) {
  counts_pca <- as_tibble(prcomp(t(vst(counts)))$x) %>% mutate(group=as_factor(batches))
  sva_pca <- as_tibble(prcomp(t(vst(sva_corrected)))$x) %>% mutate(group=as_factor(batches))
  ggplot() +
    geom_point(data=counts_pca, aes(x=PC1, y=PC2, color="Original", shape=group)) +
    geom_point(data=sva_pca, aes(x=PC1, y=PC2, color="SVA", shape=group))
}

batches <- filtermetadata$gender
diagnostic <- TRUE
if(diagnostic) {
  counts_pca <- as_tibble(prcomp(t(vst(counts)))$x) %>% mutate(group=as_factor(batches))
  sva_pca <- as_tibble(prcomp(t(vst(sva_corrected)))$x) %>% mutate(group=as_factor(batches))
  ggplot() +
    geom_point(data=counts_pca, aes(x=PC1, y=PC2, color="Original", shape=group)) +
    geom_point(data=sva_pca, aes(x=PC1, y=PC2, color="SVA", shape=group))
}

dds <- DESeqDataSetFromMatrix(countData = sva_corrected, colData = filtermetadata, design = ~singleorpaired+gender+type_tested)
dds <- collapseReplicates( dds,groupby = dds$label,run = dds$label )
DEdds <- DESeq(dds)

filename = paste(outdir, "normcounts.txt", sep="")
normcounts <- as.data.frame(counts(DEdds, normalize=TRUE))
write.table(normcounts, file = filename, append = FALSE, sep = "\t")


resWTDp16 <- results(DEdds, contrast=c("type_tested","Dp16","WT"))
resWTRcan1 <- results(DEdds, contrast=c("type_tested","Rcan1het","WT"))
resWTDp16Rcan1 <- results(DEdds, contrast=c("type_tested","Dp16dsRcan1","WT"))
resDp16Dp16Rcan1 <- results(DEdds, contrast=c("type_tested","Dp16dsRcan1","Dp16"))
resRcan1hetDp16Rcan1 <- results(DEdds, contrast=c("type_tested","Dp16dsRcan1","Rcan1het"))

graphsandsinks<-function(res, resname){
  filename = paste(outdir, resname, "_summary.txt", sep="")
  sink(filename)
  summary(res)
  sink()
  resexpressed <- as.data.frame(res)
  filename = paste(outdir, resname, "_res.txt", sep="")
  write.table(resexpressed, file=filename, quote = FALSE, sep="\t")
  resexpressed <-resexpressed[!is.na(resexpressed$padj),]
  filename = paste(outdir, resname, "_expressedgenenames.txt", sep="")
  write.table(rownames(resexpressed), file=filename, row.names=FALSE, col.names = FALSE, quote = FALSE)
  resSig <-as.data.frame(subset(res, padj < 0.1))
  resSig<- resSig[ order( resSig$padj ), ]
  filename = paste(outdir, resname, "_siggenenames.txt", sep="")
  write.table(rownames(resSig), file=filename, row.names=FALSE, col.names = FALSE, quote = FALSE)
  rnkdf <- tibble(gene = rownames(res),rnk = -log(res$pvalue) * sign(res$log2FoldChange)) %>%arrange(desc(rnk)) %>% drop_na()
  filename = paste(outdir, resname, "_deseq_res_for_gsea.rnk", sep="")
  write.table(rnkdf, file = filename, append = FALSE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  jpeg(paste0(outdir,resname,'_plotMA.jpg', sep=""))
  plotMA(res)
  dev.off()
}


res <-resDp16Dp16Rcan1
resname <- "resDp16Dp16Rcan1"
graphsandsinks(res, resname)
res <-resWTDp16
resname <- "resWTDp16"
graphsandsinks(res, resname)
res <-resWTRcan1
resname <- "resWTRcan1"
graphsandsinks(res, resname)
res <-resWTDp16Rcan1
resname <- "resWTDp16Rcan1"
graphsandsinks(res, resname)
res <-resRcan1hetDp16Rcan1
resname <- "resRcan1hetDp16Rcan1"
graphsandsinks(res, resname)
