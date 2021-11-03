library(EnrichedHeatmap)
library(plyranges)
library(circlize)

#outdir = "/Users/allenma/motif_strand/onesamplep53/"
outdir = "/Users/allenma/motif_strand/manysamples/"
metadatafile=paste(outdir,"metadata.txt", sep="")
line = "This is the metadata for this directory."
write(line,file=metadatafile)
#run date
#number of mu
#number of motif
#all stuff below

#suff to put in metadata
#mu for one sample
#mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105739.sorted-1_bidir_predictions.bed"
mufile="/Users/allenma/mumerge/p53results/merge1hour__MUMERGE.bed"
line = paste("mu file=", mufile, sep="")
write(line,file=metadatafile,append=TRUE)
#or mumerge file for a group of samples
#mufile="/Users/allenma/mumerge/p53results/merge1hour__MUMERGE.bed"
motiffile = "/scratch/Shares/dowell/motifs/HOCOMOCO_HUMAN_v11_p1e-5_grch38/P53_HUMAN.H11MO.0.A.bed"
line = paste("motif file=", motiffile, sep="")
write(line,file=metadatafile,append=TRUE)
usecentermotif=TRUE
line = paste("usecentermotif=", usecentermotif, sep="")
write(line,file=metadatafile,append=TRUE)
windowaroundmu=1500
line = paste("windowaroundmu=", windowaroundmu, sep="")
write(line,file=metadatafile,append=TRUE)

#import the motifs
motifs = read.table(motiffile)
motifs["loc"] = paste(motifs$V1,":",  motifs$V2,"-",motifs$V3, sep="")
grmotifs <- GRanges(seqnames = motifs[[1]], ranges = IRanges(motifs[[2]], motifs[[3]]), score = motifs[[5]], name= motifs[[4]],  strand=motifs[[6]], loc=motifs[,"loc"])
#grab the center of the motifs and Tfit calls and expand the motif regions by 1500 on both sides
if(usecentermotif==TRUE){
  grmotifs=  mutate(anchor_center(grmotifs), width = 1)}
n_occur <- data.frame(table(motifs$loc))
bothstrands <- motifs[motifs$loc %in% n_occur$Var1[n_occur$Freq > 1],]
onestrand <- motifs[motifs$loc %in% n_occur$Var1[n_occur$Freq == 1],]
grmotifsbothstrands <- GRanges(seqnames = bothstrands[[1]], ranges = IRanges(bothstrands[[2]], bothstrands[[3]]), score = bothstrands[[5]], name= bothstrands[[4]],  strand=bothstrands[[6]], loc=bothstrands[,"loc"])
grmotifsonestrand <- GRanges(seqnames = onestrand[[1]], ranges = IRanges(onestrand[[2]], onestrand[[3]]), score = onestrand[[5]], name=onestrand[[4]],  strand=onestrand[[6]], loc=onestrand[,"loc"])

line = paste("motif file, n=", length(grmotifs))
write(line,file=metadatafile,append=TRUE)

line = paste("There are , n=", length(grmotifsbothstrands)/2, " motifs on both strands. Therefore they are n=",length(grmotifsbothstrands), "of the total motifs.")
write(line,file=metadatafile,append=TRUE)



#read the mu file
mus = read.table(mufile)
#if I'm using the mu-merge I have to make up my own names for the regions
if (length(colnames(mus))==3){mus$name = paste(mus$V1,":",  mus$V2,"-",mus$V3, sep="")}
grmus = GRanges(seqnames = mus[[1]], ranges = IRanges(mus[[2]], mus[[3]]), names=mus[[4]])
#grab the center of the mu file
grmus<-mutate(anchor_center(grmus), width = 1)

line = paste("mu file, n=", length(grmus))
write(line,file=metadatafile,append=TRUE)

#mu near motifs
widegrmotifs=  mutate(anchor_center(grmotifs), width = windowaroundmu*2)
munearmotifswide <- join_overlap_inner(grmus,widegrmotifs) #This will give repeat the mu if the mu is by two motifs
grmunearmotifs <-granges(munearmotifswide) #removes the subject so you can use the queary ranges again
mcols(grmunearmotifs)$name = munearmotifswide$names
grmunearmotifsnorepeats<- unique(grmunearmotifs)

line = paste("There are ", length(munearmotifswide), " mu near motif (with repeats). This number includes repeats when a single mu is near multiple motifs. ", sep="")
write(line,file=metadatafile,append=TRUE)

line = paste("There are ", length(grmunearmotifsnorepeats), " mu near motif. This number does not include repeats when a single mu is near multiple motifs. ")
write(line,file=metadatafile,append=TRUE)


distsfrommu=distanceToNearest(grmunearmotifs,grmotifs,select="arbitrary")
distsfrommu_norepeats=distanceToNearest(grmunearmotifsnorepeats,grmotifs,select="arbitrary")
numberofregions <- length(grmunearmotifs)
numberofregionsnorepeats <- length(grmunearmotifsnorepeats)

nofbins=windowaroundmu
png(paste(outdir,"dis_from_mu_to_motif_warn_repeatmus.bins",nofbins,".png", sep=""))
hist(mcols(distsfrommu)$distance, breaks=nofbins, main=paste("center on mu n=",numberofregions, " mu can repeat if more than one motif nearby.",sep=""))
dev.off()
png(paste(outdir,"dis_from_mu_to_motif.bins",nofbins,".png", sep=""))
hist(mcols(distsfrommu_norepeats)$distance, breaks=nofbins, main=paste("center on mu n=",numberofregionsnorepeats, sep=""))
dev.off()

windowaroundmudiv10=windowaroundmu/10
nofbins=windowaroundmudiv10
png(paste(outdir,"dis_from_mu_to_motif_warn_repeatmus.bins",nofbins,".png", sep=""))
hist(mcols(distsfrommu)$distance, breaks=nofbins, main=paste("center on mu n=",numberofregions, " mu can repeat if more than one motif nearby.",sep=""))
dev.off()
png(paste(outdir,"dis_from_mu_to_motif.bins",nofbins,".png", sep=""))
hist(mcols(distsfrommu_norepeats)$distance, breaks=nofbins, main=paste("center on mu n=",numberofregionsnorepeats, sep=""))
dev.off()

#mu near motifs
widegrmu=  mutate(anchor_center(grmus), width = windowaroundmu*2)
motifnearmuswide <- join_overlap_inner(grmotifs,widegrmu) #This will give repeat the motif if the motif is by two mu
grmotifnearmuswide <-granges(motifnearmuswide) #removes the subject so you can use the queary ranges again
mcols(grmotifnearmuswide)$name = motifnearmuswide$names
grmotifnearmuswidenorepeats<- unique(grmotifnearmuswide)


distsfrommotif=distanceToNearest(grmotifnearmuswide,grmus,select="arbitrary")
distsfrommotif_norepeats=distanceToNearest(grmotifnearmuswidenorepeats,grmus,select="arbitrary")
numberofregions <- length(distsfrommotif)
numberofregionsnorepeats <- length(distsfrommotif_norepeats)

nofbins=windowaroundmu
png(paste(outdir,"dis_from_motif_to_mu_warn_repeatmus.bins",nofbins,".png", sep=""))
hist(mcols(distsfrommotif)$distance, breaks=nofbins, main=paste("center on mu n=",numberofregions, " mu can repeat if more than one motif nearby.",sep=""))
dev.off()
png(paste(outdir,"dis_from_motif_to_mu.bins",nofbins,".png", sep=""))
hist(mcols(distsfrommotif_norepeats)$distance, breaks=nofbins, main=paste("center on motif n=",numberofregionsnorepeats, sep=""))
dev.off()

windowaroundmudiv10=windowaroundmu/10
nofbins=windowaroundmudiv10
png(paste(outdir,"dis_from_motif_to_mu_warn_repeatmus.bins",nofbins,".png", sep=""))
hist(mcols(distsfrommotif)$distance, breaks=nofbins, main=paste("center on mu n=",numberofregions, " mu can repeat if more than one motif nearby.",sep=""))
dev.off()
png(paste(outdir,"dis_from_motif_to_mu.bins",nofbins,".png", sep=""))
hist(mcols(distsfrommotif_norepeats)$distance, breaks=nofbins, main=paste("center on motif n=",numberofregionsnorepeats, sep=""))
dev.off()

df <- as.data.frame(distsfrommu_norepeats)
q <- as.data.frame(grmunearmotifsnorepeats[df$queryHits]) #seqnames     start       end width strand
colnames(q) <- paste("mu", colnames(q), sep = "_")
s<- as.data.frame(grmotifs[df$subjectHits]) #start       end width strand   score "name"        "loc"         "df$distance"
colnames(s) <- paste("motif", colnames(s), sep = "_")
qs <- cbind(q,s, df$distance)
qs$mu_strand <- as.character(qs$mu_strand)
qs$motif_strand <- as.character(qs$motif_strand)
qs$strand_cat <- ifelse(qs$motif_name %in% mcols(grmotifsbothstrands)$name, paste("both",qs$motif_strand, sep="_"), qs$motif_strand)
qs$motif_middle <- round((qs$motif_end-qs$motif_start)/2+qs$motif_start)
qs$motif_relativetomu <- ifelse(qs$mu_start<qs$motif_middle, "left", "right")
qs$group <- paste("motifstrand_", qs$strand_cat, qs$motif_relativetomu, sep="_")
qsmunearmotif <- GRanges(seqnames = qs$mu_seqnames, ranges = IRanges(qs$mu_start, qs$mu_end, distomotif=qs$`df$distance`))
write.csv(qs, paste(outdir,"munearmotifs.csv", sep=""))

qsmunearmotif <- qsmunearmotif %>%  filter(seqnames=="chr21")

line = paste("There are ", length(qsmunearmotif), "unique mu near motifs. This number does not include repeats when a single mu is near multiple motifs. ", sep="")
write(line,file=metadatafile,append=TRUE)


###### list the samples
#sampledir="/scratch/Shares/dowell/dbnascent/out/Allen2014global/bedgraphs/"
sampledir="/scratch/Users/allenma/runallen2014/Allen2014global/mapped/rcc_bedgraphs/"
samples <-c("SRR1105736", "SRR1105738", "SRR1105737", "SRR1105739", "SRR1105740", "SRR1105741")
#samplepostfix <- ".sorted.bedGraph"
samplepostfix <- ".rcc.bedGraph"
samplenames <- c("DMSO_rep1", "Nutlin_rep1", "DMSO_rep2", "Nutlin_rep2", "p53null_DMSO", "p53null_Nutlin")

quants <- list()
col_dis_fun = colorRamp2(c(0,150, 1500), c("yellow", "lightgreen","green"))
disfrom_mu = Heatmap(qsmunearmotif$distomotif, name = "distance to motif", col=col_dis_fun,width = unit(5, "mm"), use_raster = TRUE)
ht_list <- disfrom_mu
ht_listpos <- ht_list
ht_listneg <- ht_list

matposlist <- list()
matneglist <- list()

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
  gr <- gr %>% filter(seqnames=="chr21")#heret
  pos = gr %>% filter(coverage>0)
  neg = gr %>% filter(coverage<0)
  mcols(neg)$coverage <-mcols(neg)$coverage*-1
  line = paste("proseq bedgraph file, lines=", length(gr))
  write(line,file=metadatafile,append=TRUE)
  numberofregions = as.character(length(qsmunearmotif))
  pos <= pos %>% filter(seqnames=="chr21")
  mat0pos <- normalizeToMatrix(pos, qsmunearmotif, value_column = "coverage",extend = 1500)
  matposlist[[i]]<-mat0pos
  neg <= neg %>% filter(seqnames=="chr21")
  mat0neg <- normalizeToMatrix(neg, qsmunearmotif, value_column = "coverage",extend = 1500)
  matneglist[[i]]<-mat0neg
  max99 <- quantile(mat0pos, c(0.99))
  quants[i]<-max99
}

maxquant <- max(unlist(quants))
#maxquant<-10
col_fun = colorRamp2(c(0, maxquant), c("white", "blue"))


for (i in 1:length(samplenames)){
  samplename = samplenames[i]
  message(i)
  mat0pos <- matposlist[[i]]
  mat0neg <- matneglist[[i]]
  p<-EnrichedHeatmap(mat0pos, col = col_fun, use_raster = TRUE, name=paste(samplename, " Pro-seq positive strand", sep=""), column_title = samplename)
  n<-EnrichedHeatmap(mat0neg, col = col_fun, use_raster = TRUE,  name=paste(samplename, " Pro-seq negative strand", sep="", column_title = samplename))
  ht_listpos <- ht_listpos+p
  ht_listneg <- ht_listneg+n}


#hm_list_1<-ht_listpos[]

#hm_mat_mean_1 = getSignalsFromList(hm_list_1, mean)
#hm_mat_mean_2 = getSignalsFromList(hm_list_2, mean)
#hm_mat_diff = hm_mat_mean_1 - hm_mat_mean_2

png(paste(outdir,"mu_near_motifs_centermu_EnrichedHeatmap_of_pro-seq_pos_for_coverage_wdistance.png", sep=""))
draw(ht_listpos, column_title = paste("PRO-seq over mu", "_", numberofregions , sep=""))
dev.off()
row_order = row_order(ht_listpos)
png(paste(outdir,"mu_near_motifs_centermu_EnrichedHeatmap_of_pro-seq_neg_for_coverage_wdistance.png", sep=""))
draw(ht_listneg, row_order = row_order,column_title = paste("PRO-seq over mu", "_", numberofregions , sep=""))
dev.off()

#needs to be fixed below this
png(paste(outdir,"mu_near_motifs_centermu_EnrichedHeatmap_of_pro-seq_for_coverage.png", sep=""))
draw(, column_title = paste("PRO-seq over mu", "_", numberofregions , sep=""))
dev.off()

numberofregions = as.character(length(qsmunearmotif))
mat0pos <- normalizeToMatrix(pos, qsmunearmotif,extend = 1500)
mat0neg <- normalizeToMatrix(neg, qsmunearmotif,extend = 1500)
col_fun = colorRamp2(quantile(mat0pos, c(0, 1)), c("white", "blue"))
p<-EnrichedHeatmap(mat0pos, col = col_fun, use_raster = TRUE, name="pos strand pro-seq")
n<-EnrichedHeatmap(mat0neg, col = col_fun, use_raster = TRUE, name="neg strand pro-seq")
col_dis_fun = colorRamp2(c(0,150, 1500), c("yellow", "lightgreen","green"))
disfrom_mu = Heatmap(qsmunearmotif$distomotif, name = "distance to motif", col=col_dis_fun,width = unit(5, "mm"), use_raster = TRUE)
ht_list <- disfrom_mu+n+p
png(paste(outdir,"mu_near_motifs_centermu_EnrichedHeatmap_of_pro-seq_for_signalsetto1_wdistance.png", sep=""))
draw(ht_list, column_title = paste("PRO-seq over mu", "_", numberofregions , sep=""))
dev.off()
ht_list <- n+p
png(paste(outdir,"mu_near_motifs_centermu_EnrichedHeatmap_of_pro-seq_for_signalsetto1.png", sep=""))
draw(ht_list, column_title = paste("PRO-seq over mu", "_", numberofregions , sep=""))
dev.off()

