#install.packages("/Users/allenma/lasso2_1.2-20.tar.gz", repos = NULL, type="source")
#BiocManager::install("DEGreport")

library(gplots)
library(dplyr)
library("tibble")
library(tidyr)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(RColorBrewer)
library(corrgram)

theme_set(
  theme_classic(base_size = 18)
)


RNAbed <-"/scratch/Shares/public/genomes/allfilesmadefromigenomesfiles/bed12from_gene_idgtf/mm10.bed"
RNAindir="/Shares/down/RNA/hoeffer/Oct2020/"
outdir="/Shares/down/RNA/hoeffer/Oct2020/"
filetable <- read.csv(paste0(outdir, "fileinfo.tech.csv"), header=TRUE)
RNAcoveragedat="nextflowhoeffermm10_072920a_and_b_res_featureCounts_gene_idfull_213620.RData"
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19")
masterannotationdf <-read.table(RNAbed, sep="\t", col.names=c("chr", "start", "stop", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"))
masterannotationdf <- masterannotationdf[masterannotationdf$chr %in% minichrs,]
masterannotationdf_onlyDp16tri <- masterannotationdf %>% filter(chr=="chr16")  %>% filter(start>=75540513) %>% filter(stop<=97962622)
masterannotationdf_onlyDp16tri_removercan <-masterannotationdf_onlyDp16tri %>% filter(name!="Rcan1") 
anuplodygenes <-masterannotationdf_onlyDp16tri[["name"]] #Dp16
anuplodygenesDp16 <-anuplodygenes
length(anuplodygenes)
anuplodygenesnotRcan1 <- masterannotationdf_onlyDp16tri_removercan[["name"]] 
length(anuplodygenesnotRcan1)
#masterannotationdf_onlyTs1Cjetri <- masterannotationdf %>% filter(chr=="chr16")  %>% filter(start>=90220762) %>% filter(stop<=97462906) #why does this end after Dp16?
#masterannotationdf_onlyTs1Cjetri[masterannotationdf_onlyTs1Cjetri$name=="NM_001081549.2",]
#Ts1Cje has Rcan1
anuplodygenes <-masterannotationdf_onlyDp16tri[["name"]] #Dp16
baseploidy <- 2
alt_ploidy <-3


#this is loading the data from feature counts
load(paste0(RNAindir, RNAcoveragedat))
#this is loading the metadata
RNAfiletable <- read.csv(paste0(RNAindir, "fileinfo.tech.csv"), header=TRUE)
RNAsamples<-RNAfiletable
RNAsamples$samplegroup <-paste(RNAsamples$type, "_", RNAsamples$RNA.prep, sep="")
#RNAsamples <- RNAsamples %>% mutate(ploidy = ifelse(startsWith(as.vector(type), "Dp"), "ploidy_Dp16", "ploidy_typical"))
RNAsamples <- RNAsamples %>% mutate(ploidy = ifelse(type=="Dp16", "ploidy_Dp16", ifelse(type=="Dp16dsRcan1", "ploidy_Dp16notRcan","ploidy_typical" )))

mkdirs <- function(fp) {
  if(!file.exists(fp)) {
    mkdirs(dirname(fp))
    dir.create(fp)
  }
} 


outfilename="allsamples"
mkdirs(paste0(outdir, outfilename, sep="/"))
#make sure that ddsFull is set right for this sample. It includes the date variable-but the selections below can't incude the date variable. 




#only use 1/21/19
#outfilename="012119samples"
#mkdirs(paste0(outdir, outfilename, sep="/"))
#RNAsamples <- RNAsamples %>%filter(RNA.prep_date=="1/21/19")

#outfilename="012119samples_drop10and11"
#mkdirs(paste0(outdir, outfilename, sep="/"))
#RNAsamples <- RNAsamples %>%filter(RNA.prep_date=="1/21/19")
#RNAsamples <- RNAsamples  %>%filter(samplegroup!="Rcan1het_11") %>%filter(samplegroup!="Dp16dsRcan1_10")

#outfilename="RNA02072020allsamples"
#mkdirs(paste0(outdir, outfilename, sep="/"))
#RNAsamples <- RNAsamples %>%filter(RNA.prep_date!="1/21/19")

#outfilename="RNA02072020drop_22and28"
#mkdirs(paste0(outdir, outfilename, sep="/"))
#RNAsamples <- RNAsamples %>%filter(RNA.prep_date!="1/21/19")
#RNAsamples <- RNAsamples %>%filter(samplegroup!="Dp16dsRcan1_22") %>%filter(samplegroup!="Dp16_23")%>%filter(samplegroup!="WT_28")

#RNAsamples <- RNAsamples %>%filter(samplegroup!="Dp16dsRcan1_22") %>%filter(samplegroup!="Dp16_23")%>%filter(samplegroup!="WT_28") %>%filter(samplegroup!="Rcan1het_11") %>%filter(samplegroup!="Dp16dsRcan1_10") %>%filter(samplegroup!="Dp16dsRcan1_6")

#outfilename="temp"
#mkdirs(paste0(outdir, outfilename, sep="/"))
#RNAsamples <- RNAsamples %>%filter(RNA.prep_date=="1/21/19")
#RNAsamples <- RNAsamples  %>%filter(samplegroup!="Rcan1het_11") %>%filter(samplegroup!="Dp16dsRcan1_10")


dim(coverage$counts)

#this is where you would remove samples if they are bad
countdat <- coverage$counts[,colnames(coverage$counts) %in% RNAsamples$label]
dim(countdat)

#remove genes that are not in the masterchrlist
countdat<-countdat[noquote(row.names(countdat)) %in% unlist(masterannotationdf$name),]
dim(countdat)



#run Deseq on RNA-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = countdat, colData = RNAsamples, design = ~ RNA.prep_date+mother+type) #use this for all samples
#ddsFull <- DESeqDataSetFromMatrix(countData = countdat, colData = RNAsamples, design = ~ mother+type)
ddsFull$type <- relevel(ddsFull$type, "WT")
dds <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )
ddsFull$samplegroup
dds$samplegroup

controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) 
samplepinfo<-as.data.frame(colData(dds))
ploidy_Dp16 = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_Dp16notRcan = ifelse(rownames(dds) %in% anuplodygenesnotRcan1, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))
normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy
unique(normFactors)


ddsCollapsed<-DESeq(dds)

ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
#and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.
#Adds a ton of columns to elementMetadata(ddsCollapsed) 
#including Intercept and all values they expect you to compare 
#SE_Intercept (SE stands for standard error)  and SE for every thing they expect you to compare. 
#Wald static for what they expect you to compare. 
#betaConv  betaIter         deviance  maxCooks

jpeg(paste0(outdir,outfilename,"/",'plotestdisp.jpg', sep=""))
plotDispEsts(ddsCollapsed)
dev.off()

jpeg(paste0(outdir,outfilename,"/",'plotestdisp_correctforTri.jpg', sep=""))
plotDispEsts(ddsCollapsed_normfactor)
dev.off()


resWTDp16 <- results(ddsCollapsed_normfactor, contrast=c("type","Dp16","WT"))
resWTRcan1 <- results(ddsCollapsed_normfactor, contrast=c("type","Rcan1het","WT"))
resWTDp16Rcan1 <- results(ddsCollapsed_normfactor, contrast=c("type","Dp16dsRcan1","WT"))
resDp16Dp16Rcan1 <- results(ddsCollapsed_normfactor, contrast=c("type","Dp16dsRcan1","Dp16"))
resRcan1hetDp16Rcan1 <- results(ddsCollapsed_normfactor, contrast=c("type","Dp16dsRcan1","Rcan1het"))

thepvalsWTDp16 <- as.data.frame(resWTDp16)
colnames(thepvalsWTDp16)<-paste(colnames(thepvalsWTDp16), "WTDp16", sep=".")
thepvalsWTDp16$Row.names<-row.names(thepvalsWTDp16)

thepvalsWTRcan1 <- as.data.frame(resWTRcan1)
colnames(thepvalsWTRcan1)<-paste(colnames(thepvalsWTRcan1), "WTRcan1", sep=".")
thepvalsWTRcan1$Row.names<-row.names(thepvalsWTRcan1)

thepvalsWTDp16Rcan1 <- as.data.frame(resWTDp16Rcan1)
colnames(thepvalsWTDp16Rcan1)<-paste(colnames(thepvalsWTDp16Rcan1), "WTDp16Rcan1", sep=".")
thepvalsWTDp16Rcan1$Row.names<-row.names(thepvalsWTDp16Rcan1)

thepvalsDp16Dp16Rcan1 <- as.data.frame(resDp16Dp16Rcan1)
colnames(thepvalsDp16Dp16Rcan1)<-paste(colnames(thepvalsDp16Dp16Rcan1), "Dp16Dp16Rcan1", sep=".")
thepvalsDp16Dp16Rcan1$Row.names<-row.names(thepvalsDp16Dp16Rcan1)

thepvalsRcan1hetDp16Rcan1<- as.data.frame(resRcan1hetDp16Rcan1)
colnames(thepvalsRcan1hetDp16Rcan1)<-paste(colnames(thepvalsRcan1hetDp16Rcan1), "Rcan1hetDp16", sep=".")
thepvalsRcan1hetDp16Rcan1$Row.names<-row.names(thepvalsRcan1hetDp16Rcan1)


thepvals <- merge(thepvalsWTDp16, thepvalsWTRcan1,by="Row.names", all=TRUE)
thepvals <- merge(thepvals, thepvalsWTDp16Rcan1,by="Row.names", all=TRUE)
thepvals <- merge(thepvals, thepvalsDp16Dp16Rcan1,by="Row.names", all=TRUE)
thepvals <- merge(thepvals, thepvalsRcan1hetDp16Rcan1,by="Row.names", all=TRUE)

thepvals_highex <- thepvals[thepvals$baseMean.WTDp16>1000,]


thepvalsonly <- thepvals %>% select(contains("pvalue"))
thepvalsonly_highex <- thepvals_highex %>% select(contains("pvalue"))
cor(thepvalsonly, use="pairwise.complete.obs", method="spearman")
cor(thepvalsonly_highex, use="pairwise.complete.obs", method="spearman")

corrgram(thepvalsonly, order=NULL, lower.panel=panel.shade,
         upper.panel=NULL, text.panel=panel.txt,
         main="p values",cor.method="spearman")

corrgram(thepvalsonly_highex, order=NULL, lower.panel=panel.shade,
         upper.panel=NULL, text.panel=panel.txt,
         main="p values", cor.method="spearman")

jpeg(paste0(outdir,outfilename,"/",'plotMAresWTDp16.jpg', sep=""))
plotMA(resWTDp16)
dev.off()
jpeg(paste0(outdir,outfilename,"/",'plotMAresWTRcan1.jpg.jpg', sep=""))
plotMA(resWTRcan1)
dev.off()
jpeg(paste0(outdir,outfilename,"/",'plotMAresWTDp16Rcan1.jpg', sep=""))
plotMA(resWTDp16Rcan1)
dev.off()
jpeg(paste0(outdir,outfilename,"/",'plotMAresDp16Dp16Rcan1.jpg', sep=""))
plotMA(resDp16Dp16Rcan1)
dev.off()
jpeg(paste0(outdir,outfilename,"/",'plotMAresRcan1hetDp16Rcan1.jpg', sep=""))
plotMA(resRcan1hetDp16Rcan1)
dev.off()

summary(resWTDp16)
summary(resWTRcan1)
summary(resWTDp16Rcan1)
summary(resDp16Dp16Rcan1)
summary(resRcan1hetDp16Rcan1)


res <-resDp16Dp16Rcan1
rnkdf <- tibble(gene = rownames(res),
                rnk = -log(res$pvalue) * sign(res$log2FoldChange)) %>%
  arrange(desc(rnk)) %>% drop_na()

## Write out the table without any additional information
write.table(rnkdf, paste0(outdir,outfilename,"/","deseq_res_for_gsea.rnk",sep=""),append = FALSE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

resSig <- as.data.frame(subset(resDp16Dp16Rcan1, padj < 0.1))
resSig$name <- as.factor(noquote(row.names(resSig)))
resSig <- left_join(resSig, masterannotationdf, by.x = name, by.y = name)
bedfile <- resSig %>% select(chr, start, stop,name, padj, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts)
dim(bedfile)
write.table(resSig,paste0(outdir,outfilename,"/",'sigDp16Dp16Rcan1.txt', sep=""), sep="\t",append = FALSE, quote = FALSE)
write.table(bedfile,paste0(outdir,outfilename,"/",'sigDp16Dp16Rcan1.bed', sep=""), sep="\t",row.names=FALSE, col.names=FALSE,append = FALSE, quote = FALSE)

normcounts0 <- as.data.frame(counts(ddsCollapsed, normalize=TRUE))
normcounts <- as.data.frame(counts(ddsCollapsed_normfactor, normalize=TRUE))
unnorm_counts=normcounts[,1]*normFactors[,1]
for (i in 2:dim(normcounts)[2]){unnorm_counts<-cbind(unnorm_counts,normcounts[,i]*normFactors[,i])}
unnorm_counts <- as.data.frame(unnorm_counts)
colnames(unnorm_counts)<-colnames(normcounts)
rownames(unnorm_counts)<-rownames(normcounts)
orderofcols <-c("WT_1", "WT_4", "WT_8", "WT_12", "Rcan1het_3", "Rcan1het_7", "Rcan1het_11" ,"Rcan1het_16", "Dp16_5", "Dp16_9","Dp16_13", "Dp16_14","Dp16dsRcan1_2", "Dp16dsRcan1_6","Dp16dsRcan1_10","Dp16dsRcan1_15","WT_17", "WT_19", "WT_25", "WT_28", "Dp16_18", "Dp16_20", "Dp16_23", "Dp16_26", "Dp16dsRcan1_21", "Dp16dsRcan1_22","Dp16dsRcan1_24", "Dp16dsRcan1_27")
#normcounts0 <- normcounts0 %>% select(any_of(orderofcols))
#normcounts <- normcounts %>% select(any_of(orderofcols))
orderofcols <- orderofcols[orderofcols %in% colnames(unnorm_counts)]
unnorm_counts <- unnorm_counts %>% select(all_of(orderofcols))

corrgram(unnorm_counts, order=NULL, lower.panel=panel.shade,
         upper.panel=NULL, text.panel=panel.txt,
         main="counts", cor.method="pearson")

#plot a heat map of all sig genes

normcounts_select <- unnorm_counts[rownames(unnorm_counts) %in% noquote(resSig$name),]
normcounts_select <-as.matrix(normcounts_select)
jpeg(paste0(outdir,outfilename,"/",'sigDp16Dp16Rcan1.jpg', sep=""))
heatmap(normcounts_select,Colv = NA)
write.table(normcounts_select,paste0(outdir,outfilename,"/",'sigDp16Dp16Rcan1.counts.txt', sep=""), sep="\t",append = FALSE, quote = FALSE)
dev.off()





#plot Rcan1 remember I've corrected Dp16 and over corrected Dp16Rcan1het
normcounts_select <- unnorm_counts[rownames(unnorm_counts) %in% noquote(anuplodygenes),]
normcounts_select <-as.matrix(normcounts_select)
jpeg(paste0(outdir,outfilename,"/",'aneuploid_genes.jpg', sep=""))
heatmap(normcounts_select,Colv = NA)
write.table(normcounts_select,paste0(outdir,outfilename,"/",'aneuploid_genes.txt', sep=""), sep="\t",append = FALSE, quote = FALSE)
dev.off()


#plot Rcan1 remember I've corrected Dp16 and over corrected Dp16Rcan1het
Rcan1genes <- c("Rcan1")
normcounts_select <- unnorm_counts[rownames(unnorm_counts) %in% noquote(Rcan1genes),]
normcounts_select <-as.matrix(normcounts_select)
jpeg(paste0(outdir,outfilename,"/",'Rcan1.jpg', sep=""))
heatmap(normcounts_select,Colv = NA)
write.table(normcounts_select,paste0(outdir,outfilename,"/",'Rcan1.txt', sep=""), sep="\t",append = FALSE, quote = FALSE)
dev.off()

#plot Rcan1
plotCounts(ddsCollapsed_normfactor, intgroup="type","Rcan1")




#normcounts["NM_001081549.2",]


#plot any gene you want
plotCounts(ddsCollapsed_normfactor, intgroup="type","NM_178195.2")


#for (genename in bedfile$name){
#message(genename)
#plotCounts(ddsCollapsed_normfactor, intgroup="type", genename)}
#normcounts["NM_001039934.1",]

#Dp16dsRcan1_27, Dp16_23, Dp16_26


rld <- rlog(ddsCollapsed_normfactor, blind=TRUE)
jpeg(paste0(outdir,outfilename,"/",'plotPCAtype.jpg', sep=""))
plotPCA(rld, intgroup="type")
dev.off()
jpeg(paste0(outdir,outfilename,"/",'plotRNA.prep_date.jpg', sep=""))
plotPCA(rld, intgroup="RNA.prep_date")
dev.off()
jpeg(paste0(outdir,outfilename,"/",'plotPCAsamplegroup.jpg', sep=""))
plotPCA(rld, intgroup="samplegroup")
dev.off()

rld_mat <- assay(rld)  
rld_cor <- cor(rld_mat)   
jpeg(paste0(outdir,outfilename,"/",'pheatmap_rld_cor.jpg', sep=""))
pheatmap(rld_cor)
dev.off()
dev.off()
dev.off()
dev.off()

#masterannotationdf_onlyDp16tr_mRNA <- masterannotationdf_onlyDp16tri %>% filter(grepl("NM",name))
#for (genename in masterannotationdf_onlyDp16tr_mRNA$name){
#  message(genename)
#  plotCounts(ddsCollapsed_normfactor, intgroup="type", genename)}
#c[genename,]

#normcounts["NM_173443.3",]

dds_lrt <- DESeq(ddsCollapsed_normfactor, test="LRT", reduced = ~ 1)
#normcounts <- as.data.frame(counts(dds_lrt, normalize=TRUE))
#normcountslong <- normcounts %>% gather(key = "sample", value = "signal")
#ggplot(normcountslong, aes(x = sample, y = signal))+ geom_violin(trim = FALSE) + theme(axis.text.x = element_text(angle = 90))+scale_y_continuous(trans='log2')
#genename="NM_001083319.1"
#normcounts[genename,]
#plotCounts(dds_lrt, intgroup="type",genename)
res_LRT <- results(dds_lrt)
# Subset the LRT results to return genes with padj < 0.05
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
# Subset results for faster cluster finding (for classroom demo purposes)

clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) %>%
  head(n=1000)

#when you do alll genes you should do this instead
#clustering_sig_genes <- sig_res_LRT

# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

#make meta
meta <-as.data.frame(cbind(paste(ddsCollapsed_normfactor$type),paste(ddsCollapsed_normfactor$samplegroup),paste(ddsCollapsed_normfactor$RNA.prep_date), paste(ddsCollapsed_normfactor$mother)))
colnames(meta)<-c("type","samplegroup", "RNA.prep_date", "mother")
row.names(meta)<-meta$samplegroup

dev.off()
dev.off()
dev.off()
dev.off()
jpeg(paste0(outdir,outfilename,"/",'clusters.jpg', sep=""))
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "type", col=NULL)
dev.off()

#get more info on clusters
cluster_groups <- clusters$df
for (i in unique(cluster_groups$cluster)){
  group <- clusters$df %>% filter(cluster == i)
  normcounts_group <- unnorm_counts[rownames(unnorm_counts) %in% noquote(group$genes),]
  normcounts_group <- as.matrix(normcounts_group)
  title = paste0("cluster", i, sep=" ")
  jpeg(paste0(outdir,outfilename,"/",'cluster',i,'.jpg', sep=""))
  heatmap(normcounts_group,Colv = NA, main=title)
  dev.off()
  write.table(normcounts_group,paste0(outdir,outfilename,"/",'cluster',i,'.txt', sep=""), sep="\t",append = FALSE, quote = FALSE)}




