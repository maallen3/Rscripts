#install.packages("remotes")
#remotes::install_github("lpantano/DEGreport")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

library("DESeq2")
library("DEGreport")

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
countdat <- coveragetable %>% select(as.vector(metadata$label))

#set up the deseq object
dds <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~replicate+group)

#run deseq
DEdds <- DESeq(dds, test="LRT", reduced = ~ 1)


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
res_LRT <- results(DEdds)

#look at and output the results of one comparison
res = res_LRT
summary(res)
plotMA(res)
resSig <- subset(res, padj < 0.01)
head(res[ order( res$padj ), ])

#plot one gene
plotCounts(DEdds, gene="FER1L5", intgroup=c( "group"))


padj.cutoff <-0.01
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)

length(sigLRT_genes)

#clustering takes a while so we are only going to cluster on the top 1000 genes in class
#make sure to use all genes when running this for real
clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) %>%
  head(n=1000)

# Obtain rlog values for those significant genes
rld <- rlog(DEdds, blind=TRUE)
rld_mat <- assay(rld)
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

#make a metadata frame for this anaysis... different becuase the index must be the column names for the counts
meta <-as.data.frame(cbind(paste(DEdds$label), paste(DEdds$genotype),paste(DEdds$group),paste(DEdds$replicate)))
colnames(meta)<-c("label","genotype", "group", "replicate")
row.names(meta)<-meta$label


#dev.off()
#jpeg(paste0(outdir,outfilename,"/",'clusters.jpg', sep=""))
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "group", col=NULL)
#dev.off()

#get more info on clusters
cluster_groups <- clusters$df
for (i in unique(cluster_groups$cluster)){
  group <- clusters$df %>% filter(cluster == i)
  normcounts_group <- normcounts[rownames(normcounts) %in% noquote(group$genes),]
  normcounts_group <- as.matrix(normcounts_group)
  title = paste0("cluster", i, sep=" ")
  heatmap(normcounts_group,Colv = NA, main=title)}
  

#jpeg(paste0(outdir,outfilename,"/",'cluster',i,'.jpg', sep=""))
#dev.off()
  #write.table(normcounts_group,paste0(outdir,outfilename,"/",'cluster',i,'.txt', sep=""), sep="\t",append = FALSE, quote = FALSE)







