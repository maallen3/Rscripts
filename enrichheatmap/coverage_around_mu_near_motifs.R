library(EnrichedHeatmap)
library(plyranges)
library(circlize)
library(tidyr)

outdir = "/scratch/Users/allenma/cov_around_mu_nearp53/"
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
windowaroundmu=1500
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

mufile="/Users/allenma/mumerge/p53results/merge1hour__MUMERGE.bed"
line = paste("mu file=", mufile, sep="")
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


widegrmotifs=  mutate(anchor_center(grmotifs), width = windowaroundmu*2)
munearmotifswide <- join_overlap_inner(grmus,widegrmotifs) #This will give repeat the mu if the mu is by two motifs

overlapdf <- as.data.frame(munearmotifswide)
overlapdf$loc2 <- overlapdf$loc
overlapdf = subset(overlapdf, select = -c(name) )
overlapdf <- distinct(overlapdf)
overlapdf<-overlapdf %>% tidyr::separate(loc2, c("motif_seqname","motif_loc"), sep = "([.?:])")
overlapdf<-overlapdf %>% tidyr::separate(motif_loc, c("motif_start","motif_stop"), sep = "([-])")
overlapdf$motif_start <- as.integer(overlapdf$motif_start)
overlapdf$motif_stop <- as.integer(overlapdf$motif_stop)
overlapdf$dis_from_start = overlapdf$start-overlapdf$motif_start
overlapdf$dis_from_stop = overlapdf$start-overlapdf$motif_stop
overlapdf$dis <- if_else(abs(overlapdf$dis_from_start) <= abs(overlapdf$dis_from_stop), overlapdf$dis_from_start,  overlapdf$dis_from_stop)
overlapdf$abs_dis <- abs(overlapdf$dis)
overlapdf$one_mu <- paste(as.character(overlapdf$seqnames),as.character(overlapdf$start), sep="_")
  
minidf <- overlapdf[1:30,]

minidf$rank = 0


minidf<- minidf %>% dplyr::group_by(one_mu) %>% mutate(my_ranks = order(abs_dis))



grmunearmotifs <-granges(munearmotifswide) #removes the subject so you can use the queary ranges again
mcols(grmunearmotifs)$name = munearmotifswide$names
grmunearmotifsnorepeats<- unique(grmunearmotifs)

