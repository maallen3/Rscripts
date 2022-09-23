library(plyranges)
library(tidyr)

outDir = "/scratch/Users/allenma/testdistable/"
motiffile = "/scratch/Shares/dowell/motifs/HOCOMOCO_HUMAN_v11_p1e-5_grch38/P53_HUMAN.H11MO.0.A.bed"
mufile="/Users/allenma/mumerge/p53results/merge1hour__MUMERGE.bed"



dir.create(file.path(outDir), showWarnings = FALSE)
setwd(file.path(outDir))
metadatafile=paste(outDir,"metadata.txt", sep="")
line = "This is the metadata for this directory. This program will make a distance table between a motif bed file and a mu file. This program makes several assumptions. Distances are measured from center to center of regions. If a motif is found on both strands, then the motif.csv will have both strands, but the distance table will only have this site in the table once. The motif strand with the higher motif score is kept. In this file it is possible for one mu to have several motifs. In the column dis_rank the lowest number is the motif closest to mu. In the column score_rank the lowest number is the motif with the higher score has a rank of 1. The column motif_counts is the number of motifs at that mu. ."
write(line,file=metadatafile)
line = "This program makes several assumptions. Distances are measured from center to center of regions. If a motif is found on both strands, then the motif.csv will have both strands, but the distance table will only have this site in the table once. The motif strand with the higher motif score is kept." 
write(line,file=metadatafile)
line = "In this file it expected that some of the time a signle mu can have several motifs. The column motif_counts is the number of motifs at that mu (within the windowaroundmotif from the center mu). In the column dis_rank the lowest number is the motif closest to mu. In the column score_rank the lowest number is the motif with the higher score has a rank of 1.  "
write(line,file=metadatafile)
line = paste("motif file=", motiffile, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("mu file=", mufile, sep="")
write(line,file=metadatafile,append=TRUE)
usecentermotif=TRUE
line = paste("usecentermotif=", usecentermotif, sep="")
write(line,file=metadatafile,append=TRUE)
windowaroundmotif=1500
line = paste("windowaroundmotif=", windowaroundmotif, sep="")
write(line,file=metadatafile,append=TRUE)


#make the motif csv
motifs = read.table(motiffile, col.names=c("seqnames", "start", "stop", "name", "score", "strand"))
motifs["loc"] = paste(motifs$seqnames,":",  motifs$start,"-",motifs$stop, sep="")
n_occur <- data.frame(table(motifs$loc))
bothstrands <- motifs[motifs$loc %in% n_occur$Var1[n_occur$Freq > 1],]
posstrand <- motifs[motifs$strand=="+",]
motifs<-motifs %>% mutate(strand_group=if_else(motifs$loc %in% bothstrands$loc, "both", if_else(motifs$loc %in% posstrand$loc, "+", "-")))
motifs$id <- rownames(motifs)
motifs$motif_id <- paste(motifs$loc,motifs$id, sep=";motif_")
write.csv(motifs, paste(outDir, "motif.csv", sep=""))
line = paste("motif file, n=", nrow(motifs))
write(line,file=metadatafile,append=TRUE)

#make the mu csv
mus = read.table(mufile)
#if I'm using the mu-merge I have to make up my own names for the regions but if I use a tfit file then the name is there
if (length(colnames(mus))==3){mus$name = paste(mus$V1,":",  mus$V2,"-",mus$V3, sep="")}
colnames(mus) <- c("seqnames", "start", "stop", "name")
mus["loc"] = paste(mus$seqnames,":",  mus$start,"-",mus$stop, sep="")
mus$id <- rownames(mus)
mus$region_id <-  paste(mus$loc,mus$id, sep=";region_")
write.csv(mus, paste(outDir, "mus.csv", sep=""))
line = paste("mu file, n=", nrow(mus))
write(line,file=metadatafile,append=TRUE)

#turn them into gr_ranges and if they need to be centered, center them
motifs_onestrand <- motifs %>% arrange(score)%>%distinct(loc, .keep_all = TRUE)
line = paste("motif file after removeing duplicate motif hits when found on both strands, n=", nrow(motifs_onestrand))
write(line,file=metadatafile,append=TRUE)
grmotifs <- GRanges(seqnames = motifs_onestrand$seqnames, ranges = IRanges(motifs_onestrand$start, motifs_onestrand$stop),motif_id=motifs_onestrand$motif_id, score=motifs_onestrand$score)
#grab the center of the motifs and Tfit calls and expand the motif regions by 1500 on both sides
grmotifs=  mutate(anchor_center(grmotifs), width = 1)
mcols(grmotifs)$motif_start <- start(grmotifs)


#read the mu file
grmus = GRanges(seqnames = mus$seqnames, ranges = IRanges(mus$start, mus$stop), region_id=mus$region_id)
#grab the center of the mu file
grmus<-mutate(anchor_center(grmus), width = 1)

widegrmotifs=  mutate(anchor_center(grmotifs), width = windowaroundmotif*2)
munearmotifswide <- join_overlap_inner(grmus,widegrmotifs) #This will give repeat the mu if the mu is by two motifs

overlapdf <- as.data.frame(munearmotifswide)
overlapdf <- distinct(overlapdf)
overlapdf$dis = overlapdf$start-overlapdf$motif_start
overlapdf$abs_dis <- abs(overlapdf$dis)
overlapdf$one_mu <- paste(as.character(overlapdf$seqnames),as.character(overlapdf$start), sep="_")
overlapdf <-overlapdf %>% group_by(one_mu) %>% arrange(abs_dis) %>% mutate(dis_rank = order(abs_dis))
overlapdf <-overlapdf %>% group_by(one_mu) %>% arrange(score) %>% mutate(score_rank = order(score, decreasing = TRUE))
overlapdf <-overlapdf  %>% group_by(one_mu) %>% mutate(motif_counts=n_distinct(motif_start))
overlapdf<- overlapdf  %>% mutate(barcode_bin=ceiling(dis/10)*10)
write.csv(overlapdf, paste(outDir, "distance_table.csv", sep=""))


#raw distance as barcode
barcodedf <- as.data.frame(table(overlapdf$dis))
write.csv(barcodedf, paste(outDir, "raw_barcode_vals.csv", sep=""))
#ten base pair bins as barcode
barcodedf2 <- as.data.frame(table(overlapdf$barcode_bin))
write.csv(barcodedf2, paste(outDir, "barcode_vals_10nt_bins.csv", sep=""))




