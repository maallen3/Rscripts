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
grmus = GRanges(seqnames = mus[[1]], ranges = IRanges(mus[[2]], mus[[3]]), name= mus[[4]])
#chip-seq regions
#chip-seq Tim read
chippeakdir = "/projects/dowellLab/Taatjes/Tim/170420_K00262_0091_AHJTMHBBXX_TIMNUTLIN/MACS2/"
chipsample="Nutlin1Hr_1_peaks_hg19tohg38liftover.bed"
chip = read.table(paste(chippeakdir, chipsample, sep=""))
grchip = GRanges(seqnames = chip[[1]], ranges = IRanges(chip[[2]], chip[[3]]), name= chip[[4]])

#chip-seq bedgraph



#read the bedgraph and split to postive and negative strand counts
df = read.table(proseqfile)
gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), coverage = df[[4]])
gr <- gr %>% plyranges::filter(seqnames=="chr1") #for ease of drawing just used chr1
pos = gr %>% plyranges::filter(coverage>0)
neg = gr %>% plyranges::filter(coverage<0)


# Goal: draw proseq over all p53 motifs (ignore strand)

#read the motifs
motifs = read.table(motiffile)
grmotifs <- GRanges(seqnames = motifs[[1]], ranges = IRanges(motifs[[2]], motifs[[3]]), score = motifs[[5]], name= motifs[[4]],  strand=motifs[[6]])
grmotifs <- grmotifs %>% plyranges::filter(seqnames=="chr1")
#center of the motifs
grmotifs1=  mutate(anchor_center(grmotifs), width = 1)

#create a matrix and plot it
mat1 = normalizeToMatrix(pos, grmotifs)
EnrichedHeatmap(mat1, name = "pro-seq pos strand coverage over p53 motifs")
#there are 15002 motifs on chr1
#So some motifs have transcription and some do not....

#just plot motifs that are near mu
#motifs near mu
mus = read.table(mufile)
grmus = GRanges(seqnames = mus[[1]], ranges = IRanges(mus[[2]], mus[[3]]), name= mus[[4]])
grmu1<-mutate(anchor_center(grmus), width = 1)
grmuswith1500 <-mutate(anchor_center(grmu1), width = 3000)
motifsnearmu1500 <- join_overlap_inner(grmotifs1, grmuswith1500)
grmotifsnearmu <-granges(motifsnearmu1500) #removes what ever you overlapped so you can use the ranges again
mcols(grmotifsnearmu)$name = motifsnearmu1500$name.x


#some places in the genome have a motif on both the postive and negative strands, how many?
#split the motifs to both, neg and pos dfs, but some are on both strands
#turn them into granges
motifs = read.table(motiffile)
motifs["loc"] = paste(motifs$V1,":",  motifs$V2,"-",motifs$V3, sep="")
n_occur <- data.frame(table(motifs$loc))
bothstrands <- motifs[motifs$loc %in% n_occur$Var1[n_occur$Freq > 1],]
onestrand <- motifs[motifs$loc %in% n_occur$Var1[n_occur$Freq == 1],]
grmotifsbothstrands <- GRanges(seqnames = bothstrands[[1]], ranges = IRanges(bothstrands[[2]], bothstrands[[3]]), score = bothstrands[[5]], name= bothstrands[[4]],  strand=bothstrands[[6]], loc=bothstrands[,"loc"])
grmotifsonestrand <- GRanges(seqnames = onestrand[[1]], ranges = IRanges(onestrand[[2]], onestrand[[3]]), score = onestrand[[5]], name=onestrand[[4]],  strand=onestrand[[6]], loc=onestrand[,"loc"])
negstrandmotifs =  grmotifsonestrand %>% filter(strand=="-")
posstrandmotifs =  grmotifsonestrand %>% filter(strand=="+")
negstrandmotifs_pres_bothstrand =  grmotifsbothstrands %>% filter(strand=="-")
posstrandmotifs_pres_bothstrand =  grmotifsbothstrands %>% filter(strand=="+")
negstrandmotifs$nstrand ="one"
negstrandmotifs_pres_bothstrand$nstrand ="both"
allnegstrandmotif <- c(negstrandmotifs, negstrandmotifs_pres_bothstrand)
allnegstrandmotifinchip <- join_overlap_inner(allnegstrandmotif, grchip)
posstrandmotifs$nstrand ="one"
posstrandmotifs_pres_bothstrand$nstrand ="both"
allposstrandmotif <- c(posstrandmotifs, posstrandmotifs_pres_bothstrand)
allposstrandmotifinchip <- join_overlap_inner(allposstrandmotif, grchip)


#create a matrix and plot it
mat1 = normalizeToMatrix(neg, allnegstrandmotifinchip)
partition = paste0("cluster", kmeans(mat1, centers = 3)$cluster)
lgd = Legend(at = c("cluster1", "cluster2", "cluster3"), title = "Clusters", 
             type = "lines", legend_gp = gpar(col = 2:4))
#col_fun = colorRamp2(quantile(mat1, c(0,0.5,0.9)), c("blue", "white", "red"))
ht_list = Heatmap(partition, col = structure(2:4, names = paste0("cluster", 1:3)), name = "partition",
                  show_row_names = FALSE, width = unit(3, "mm")) +
  EnrichedHeatmap(mat1,name = "pro-seq pos strand coverage over p53 motifs in chip-peaks", top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))))+
  Heatmap(allnegstrandmotifinchip$nstrand, name="Number of strands",show_row_names = FALSE, width = unit(15, "mm"),
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:4),
                                                                    outline = FALSE, axis_param = list(side = "right"))))

draw(ht_list, split = partition, annotation_legend_list = list(lgd), 
     ht_gap = unit(c(2, 8, 8), "mm"))


#create a matrix and plot it
mat1 = normalizeToMatrix(pos, allposstrandmotifinchip)
partition = paste0("cluster", kmeans(mat1, centers = 3)$cluster)
lgd = Legend(at = c("cluster1", "cluster2", "cluster3"), title = "Clusters", 
             type = "lines", legend_gp = gpar(col = 2:4))
#col_fun = colorRamp2(quantile(mat1, c(0,0.5,0.9)), c("blue", "white", "red"))
ht_list = Heatmap(partition, col = structure(2:4, names = paste0("cluster", 1:3)), name = "partition",
                  show_row_names = FALSE, width = unit(3, "mm")) +
  EnrichedHeatmap(mat1,name = "pro-seq pos strand coverage over p53 motifs in chip-peaks", top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))))+
  Heatmap(allposstrandmotifinchip$nstrand, name="Number of strands",show_row_names = FALSE, width = unit(15, "mm"),
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:4),
                                                                    outline = FALSE, axis_param = list(side = "right"))))

draw(ht_list, split = partition, annotation_legend_list = list(lgd), 
     ht_gap = unit(c(2, 8, 8), "mm"))


#ignore strand of motif, plot pro-seq on the postive strand for all p53 motifs in all samples
