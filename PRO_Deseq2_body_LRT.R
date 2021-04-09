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
dds <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~group)
#run deseq
DEdds <- DESeq(dds, test="LRT", reduced = ~1)
#DEdds <- DESeq(dds, test="LRT", reduced = ~replicate)
#plot the dispersion of the data
plotDispEsts(DEdds)

#get the results dataframes (which was not )
res_LRT <- results(DEdds)

res = res_LRT
summary(res)
plotMA(res)
resSig <- subset(res, padj < 0.01)
head(res[ order( res$padj ), ])



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
cluster_rlog <- rld_mat[sig_res_LRT$gene, ]

#make a metadata frame for this anaysis... different becuase the index must be the column names for the counts
meta <-as.data.frame(cbind(paste(DEdds$label), paste(DEdds$genotype),paste(DEdds$group),paste(DEdds$replicate), paste(DEdds$tempature)))
colnames(meta)<-c("label","genotype", "group", "replicate", "tempature")
row.names(meta)<-meta$label
meta$tempature <-relevel(meta$tempature , "thirtyseven")

#dev.off()
#jpeg(paste0(outdir,outfilename,"/",'clusters.jpg', sep=""))
#clusters <- degPatterns(cluster_rlog, metadata = meta, time = "group", col = NULL)
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "tempature", col = "genotype", reduce=TRUE)

#dev.off()

#get more info on clusters
cluster_groups <- clusters$df
for (i in unique(cluster_groups$cluster)){
  group <- clusters$df %>% filter(cluster == i)
  normcounts_group <- normcounts[rownames(normcounts) %in% noquote(group$genes),]
  normcounts_group <- as.matrix(normcounts_group)
  title = paste0("cluster", i, sep=" ")
  heatmap(normcounts_group,Colv = NA, main=title)}
  

maindf <- merge(as.data.frame(res), cluster_groups, by=0)

cluster1 = maindf %>% filter(cluster==1) %>% arrange(padj)

onegene="NM_030930"
plotCounts(DEdds, gene=onegene, intgroup=c( "group"))
