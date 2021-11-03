library(EnrichedHeatmap)
library(GenomicRanges)
library(plyranges)
library(circlize)


outdir = "/scratch/Users/allenma/cov_aroundp53/"
metadatafile=paste(outdir,"metadata.txt", sep="")
line = "This is the metadata for this directory."
write(line,file=metadatafile)
motiffile = "/scratch/Shares/dowell/motifs/HOCOMOCO_HUMAN_v11_p1e-5_grch38/P53_HUMAN.H11MO.0.A.bed"
chipfile <- "/scratch/Shares/dowell/Tim_ChIP/hg38/Peaks/broadpeaks/motif_p10-5_intersects/hct116.nutlin1hr.hg38_peaks.broadPeak"
line = paste("motif file=", motiffile, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("chip file=", chipfile, sep="")
write(line,file=metadatafile,append=TRUE)
usecentermotif=TRUE
line = paste("usecentermotif=", usecentermotif, sep="")
write(line,file=metadatafile,append=TRUE)
windowaroundmotif=1500
line = paste("windowaroundmu=", windowaroundmotif, sep="")
write(line,file=metadatafile,append=TRUE)

#import the motifs
motifs = read.table(motiffile)
motifs["loc"] = paste(motifs$V1,":",  motifs$V2,"-",motifs$V3, sep="")
grmotifs <- GRanges(seqnames = motifs[[1]], ranges = IRanges(motifs[[2]], motifs[[3]]), score = motifs[[5]], name= motifs[[4]],  strand=motifs[[6]], loc=motifs[,"loc"])
#grab the center of the motifs and Tfit calls and expand the motif regions by 1500 on both sides
if(usecentermotif==TRUE){
  grmotifs=  mutate(anchor_center(grmotifs), width = 1)}
grmotifs1500 = mutate(anchor_center(grmotifs), width = 3000)
n_occur <- data.frame(table(motifs$loc))
bothstrands <- motifs[motifs$loc %in% n_occur$Var1[n_occur$Freq > 1],]
onestrand <- motifs[motifs$loc %in% n_occur$Var1[n_occur$Freq == 1],]
grmotifsbothstrands <- GRanges(seqnames = bothstrands[[1]], ranges = IRanges(bothstrands[[2]], bothstrands[[3]]), score = bothstrands[[5]], name= bothstrands[[4]],  strand=bothstrands[[6]], loc=bothstrands[,"loc"])
grmotifsonestrand <- GRanges(seqnames = onestrand[[1]], ranges = IRanges(onestrand[[2]], onestrand[[3]]), score = onestrand[[5]], name=onestrand[[4]],  strand=onestrand[[6]], loc=onestrand[,"loc"])
line = paste("motif file, n=", length(grmotifs))
write(line,file=metadatafile,append=TRUE)
line = paste("There are , n=", length(grmotifsbothstrands)/2, " motifs on both strands. Therefore they are n=",length(grmotifsbothstrands), "of the total motifs.")
write(line,file=metadatafile,append=TRUE)

#I'm going to require my motifs to be in a chip-seq bound region

chips = read.table(chipfile)
grchip <- GRanges(seqnames = chips[[1]], ranges = IRanges(chips[[2]], chips[[3]]), score = chips[[5]])

motifs_in_chip = join_overlap_inner(grmotifs1500, grchip)
motifs_in_chip = granges(motifs_in_chip)


###### list the samples
#sampledir="/scratch/Shares/dowell/dbnascent/out/Allen2014global/bedgraphs/"
sampledir="/scratch/Users/allenma/runallen2014/Allen2014global/mapped/rcc_bedgraphs/"
samples <-c("SRR1105736", "SRR1105738", "SRR1105737", "SRR1105739", "SRR1105740", "SRR1105741")
#samplepostfix <- ".sorted.bedGraph"
samplepostfix <- ".rcc.bedGraph"
samplenames <- c("DMSO_rep1", "Nutlin_rep1", "DMSO_rep2", "Nutlin_rep2", "p53null_DMSO", "p53null_Nutlin")

#prefilter grmotifs for motifs at least 1 count in one sample

#I'm going to require my motifs to have at least some signal in any one of my pro-seq samples
#I do this by collecting the names of motifs that have singal in each sample and then grabing the unique sample names

grmotifs_with_signal <- vector(mode="list", length=length(samplenames)*2)

for (i in seq(1, length(grmotifs_with_signal), by=2)){
  message(i)
  samplenum = 1
  sampleval=samples[samplenum]
  samplename = samplenames[samplenum]
  samplenum = samplenum+1
  proseqfile<-paste(sampledir,sampleval,samplepostfix, sep="")
  line = paste("bedgraph file=", proseqfile, sep="")
  write(line,file=metadatafile,append=TRUE)
  #read the bedgraph and split to postive and negative strand counts
  df = read.table(proseqfile)
  gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), coverage = df[[4]])
  #gr <- gr %>% filter(seqnames=="chr21")#heret
  pos = gr %>% filter(coverage>0)
  neg = gr %>% filter(coverage<0)
  mcols(neg)$coverage <-mcols(neg)$coverage*-1
  line = paste("proseq bedgraph file, lines=", length(gr))
  write(line,file=metadatafile,append=TRUE)
  gr_motifs_with_pos_signal <- join_overlap_inner(motifs_in_chip, pos)
  gr_motifs_with_neg_signal <- join_overlap_inner(motifs_in_chip, neg)
  gr_motifs_with_pos_signal <- granges(gr_motifs_with_pos_signal)
  gr_motifs_with_neg_signal <- granges(gr_motifs_with_neg_signal)
  gr_motifs_with_pos_signal <- mutate(anchor_center(gr_motifs_with_pos_signal), width = 1)
  gr_motifs_with_neg_signal <- mutate(anchor_center(gr_motifs_with_neg_signal), width = 1)
  grmotifs_with_signal[[i]] <- gr_motifs_with_pos_signal
  grmotifs_with_signal[[i+1]] <- gr_motifs_with_neg_signal
}

gr_motifs_with_signal <- do.call(c, grmotifs_with_signal)
gr_motifs_with_signal <- unique(gr_motifs_with_signal)


#now that I've filtered my motifs requiring them ot be in in a chip-seq site and have some pro-seq reads in some sample
#Im' going to make my normalize matrixes

matposlist <- list()
matneglist <- list()
quants <- list()

for (i in 1:length(samplenames)){
  message(i)
  sampleval=samples[i]
  samplename = samplenames[i]
  proseqfile<-paste(sampledir,sampleval,samplepostfix, sep="")
  line = paste("bedgraph file=", proseqfile, sep="")
  write(line,file=metadatafile,append=TRUE)
  #read the bedgraph and split to postive and negative strand counts
  df = read.table(proseqfile)
  gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), coverage = df[[4]])
  #gr <- gr %>% filter(seqnames=="chr21")#heret
  pos = gr %>% filter(coverage>0)
  neg = gr %>% filter(coverage<0)
  mcols(neg)$coverage <-mcols(neg)$coverage*-1
  line = paste("proseq bedgraph file, lines=", length(gr))
  write(line,file=metadatafile,append=TRUE)
  numberofregions = as.character(length(gr_motifs_with_signal))
  #pos <= pos %>% filter(seqnames=="chr21")
  mat0pos <- normalizeToMatrix(pos, gr_motifs_with_signal, value_column = "coverage",extend = 1500, w=5)
  matposlist[[i]]<-mat0pos
  #neg <= neg %>% filter(seqnames=="chr21")
  mat0neg <- normalizeToMatrix(neg, gr_motifs_with_signal, value_column = "coverage",extend = 1500,  w=5)
  matneglist[[i]]<-mat0neg
  max99 <- quantile(mat0pos, c(0.99))
  quants[i]<-max99
}



maxquant <- max(unlist(quants))
#when max quant is 0 it has to be set manualy
#maxquant<-10
col_fun = colorRamp2(c(0, maxquant), c("white", "blue"))

# Maybe I should add a heat map that is distance to nearest other motif in list?

#run the first sample and make a list
samplename = samplenames[1]
mat0pos <- matposlist[[1]]
mat0neg <- matneglist[[1]]
p<-EnrichedHeatmap(mat0pos, col = col_fun, use_raster = TRUE, name=paste(samplename, " Pro-seq positive strand", sep=""), column_title = paste(samplename," +"))
n<-EnrichedHeatmap(mat0neg, col = col_fun, use_raster = TRUE,  name=paste(samplename, " Pro-seq negative strand", sep=""), column_title = paste(samplename," -"))
ht_listpos <- p
ht_listneg <- n

#loop through the other samples
for (i in 2:length(samplenames)){
  samplename = samplenames[i]
  message(i)
  mat0pos <- matposlist[[i]]
  mat0neg <- matneglist[[i]]
  p<-EnrichedHeatmap(mat0pos, col = col_fun, use_raster = TRUE, name=paste(samplename, " Pro-seq positive strand", sep=""), column_title = paste(samplename," +"))
  n<-EnrichedHeatmap(mat0neg, col = col_fun, use_raster = TRUE,  name=paste(samplename, " Pro-seq negative strand", sep=""), column_title = paste(samplename," -"))
  ht_listpos <- ht_listpos+p
  ht_listneg <- ht_listneg+n}


pos_diff_rep1 <- matposlist[[2]] - matposlist[[1]]
pos_diff_rep2 <- matposlist[[4]] - matposlist[[3]]
pos_diff_null <- matposlist[[6]] - matposlist[[5]]

col_fun2 <-colorRamp2(c(-0.5,  0, 0.5), c("red", "white", "green"))
  #

ht_listpos = EnrichedHeatmap(pos_diff_rep2, name="pos diff rep 2", col = col_fun2, column_title = "pos diff rep 2")+EnrichedHeatmap(pos_diff_rep1, name="pos diff rep 1", col = col_fun2,  column_title = "pos diff rep 1")+EnrichedHeatmap(pos_diff_null, name="diff null", col = col_fun2, column_title = "pos diff p53 null")+ht_listpos

draw(ht_listpos)

png(paste(outdir,"Pos_pro_near_motifs_EnrichedHeatmap_of_pro-seq_for_signalsetto1.png", sep=""), width = 1500, height = 500)
draw(ht_listpos)
dev.off()

neg_diff_rep1 <- matneglist[[2]] - matneglist[[1]]
neg_diff_rep2 <- matneglist[[4]] - matneglist[[3]]
neg_diff_null <- matneglist[[6]] - matneglist[[5]]

#

ht_listneg = EnrichedHeatmap(neg_diff_rep2, name="neg diff rep 2", col = col_fun2, column_title = "neg diff rep 2")+EnrichedHeatmap(neg_diff_rep1, name="neg diff rep 1", col = col_fun2,  column_title = "neg diff rep 1")+EnrichedHeatmap(neg_diff_null, name="neg diff null", col = col_fun2, column_title = "neg diff p53 null")+ht_listneg

draw(ht_listneg)

png(paste(outdir,"Neg_pro_near_motifs_EnrichedHeatmap_of_pro-seq_for_signalsetto1.png", sep=""), width = 1500, height = 500)
draw(ht_listneg)
dev.off()


png(paste(outdir,"both_strands_pro_near_motifs_EnrichedHeatmap_of_pro-seq_for_signalsetto1.png", sep=""), width = 1500, height = 500)
draw(ht_listpos+ht_listneg)
dev.off()

