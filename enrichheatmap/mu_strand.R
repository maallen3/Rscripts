library(EnrichedHeatmap)
library(plyranges)
library(circlize)

#set file paths
#motif path
motifdir="/scratch/Shares/dowell/motifs/HOCOMOCO_HUMAN_v11_p1e-5_grch38/"
motifname = "P53_HUMAN.H11MO.0.A.bed"
motiffile = paste(motifdir,motifname, sep="")
#bedgrpah path
begraphdir="/scratch/Shares/dowell/dbnascent/out/Allen2014global/bedgraphs/"
SRR = "SRR1105739"
proseqfile <- paste(begraphdir,SRR,".sorted.bedGraph", sep="")


# tfit file
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105739.sorted-1_bidir_predictions.bed"
mus = read.table(mufile)
mus$strand = "+"
grmus = GRanges(seqnames = mus[[1]], ranges = IRanges(mus[[2]], mus[[3]]), name= mus[[4]], strands=mus[,"strand"])
#chip-seq regions
#chip-seq Tim read
chippeakdir = "/projects/dowellLab/Taatjes/Tim/170420_K00262_0091_AHJTMHBBXX_TIMNUTLIN/MACS2/"
chipsample="Nutlin1Hr_1_peaks_hg19tohg38liftover.bed"
chip = read.table(paste(chippeakdir, chipsample, sep=""))
grchip = GRanges(seqnames = chip[[1]], ranges = IRanges(chip[[2]], chip[[3]]), name= chip[[4]])


#read the bedgraph and split to postive and negative strand counts
df = read.table(proseqfile)
gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), coverage = df[[4]])
gr <- gr %>% plyranges::filter(seqnames=="chr1") #for ease of drawing just used chr1
pos = gr %>% plyranges::filter(coverage>0)
neg = gr %>% plyranges::filter(coverage<0)


#read the motifs
motifs = read.table(motiffile)
motifs["loc"] = paste(motifs$V1,":",  motifs$V2,"-",motifs$V3, sep="")
n_occur <- data.frame(table(motifs$loc))
motifs <- merge(motifs, n_occur, by.x = "loc", by.y = "Var1")
grmotifs <- GRanges(seqnames = motifs[,"V1"], ranges = IRanges(motifs[,"V2"], motifs[,"V3"]), score = motifs[,"V5"], name= motifs[,"V4"],  strand=motifs[,"V6"], loc=motifs[,"loc"], freq=motifs[,"Freq"])
grmotifs <- grmotifs %>% plyranges::filter(seqnames=="chr1")
#center of the motifs
grmotifs1=  mutate(anchor_center(grmotifs), width = 1)

#mu near motifs
grmotifs1500=  mutate(anchor_center(grmotifs1), width = 3000)
grmotifs1500$extrastrand = strand(grmotifs1500)
grmu1<-mutate(anchor_center(grmus), width = 1)
munearmotifs1500 <- join_overlap_inner(grmu1,grmotifs1500) #might not be right...
grmunearmotifs <-granges(munearmotifs1500) #removes what ever you overlapped so you can use the ranges again



#distance from mu to motifs within 1500 (of mu to nearest near motifs 1500)
#center on mu, where is motif
d = nearest(munearmotifs1500,grmotifs,select="arbitrary",ignore.strand=TRUE)
distsfrommu=distanceToNearest(munearmotifs1500,grmotifs1,select="arbitrary",ignore.strand=FALSE)

hist(mcols(distsfrommu)$distance, breaks=300, main="center on mu, where is motif")

#this should make a new orinetation column. I need to use that column to graph things on positive and negative separate


