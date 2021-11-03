library(ComplexHeatmap)

motiffile = "/scratch/Shares/dowell/motifs/HOCOMOCO_HUMAN_v11_p1e-5_grch38/P53_HUMAN.H11MO.0.A.bed"
windowaroundmu=1500


fahrenheit_to_celsius <- function(temp_F) {
  temp_C <- (temp_F - 32) * 5 / 9
  return(temp_C)
}

createHeatmapIQR <- function(mufile, motiffile, thissamplename){
  usecentermotif==TRUE
  #read the mu file
  mus = read.table(mufile)
  message("number of mus")
  message(dim(mus))
  #if I'm using the mu-merge I have to make up my own names for the regions
  if (length(colnames(mus))==3){mus$name = paste(mus$V1,":",  mus$V2,"-",mus$V3, sep="")}
  grmus = GRanges(seqnames = mus[[1]], ranges = IRanges(mus[[2]], mus[[3]]), names=mus[[4]])
  #grab the center of the mu file
  grmus<-mutate(anchor_center(grmus), width = 1)
  length(grmus)
  #import the motifs
  motifs = read.table(motiffile)
  message("number of motifs")
  message(dim(motifs))
  motifs["loc"] = paste(motifs$V1,":",  motifs$V2,"-",motifs$V3, sep="")
  grmotifs <- GRanges(seqnames = motifs[[1]], ranges = IRanges(motifs[[2]], motifs[[3]]), score = motifs[[5]], name= motifs[[4]],  strand=motifs[[6]], loc=motifs[,"loc"])
  #grab the center of the motifs and Tfit calls and expand the motif regions by 1500 on both sides
  if(usecentermotif==TRUE){
    grmotifs=  mutate(anchor_center(grmotifs), width = 1)}
  #mu near motifs
  widegrmotifs=  mutate(anchor_center(grmotifs), width = windowaroundmu*2)
  munearmotifswide <- join_overlap_inner(grmus,widegrmotifs) #This will give repeat the mu if the mu is by two motifs
  grmunearmotifs <-granges(munearmotifswide) #removes the subject so you can use the queary ranges again
  mcols(grmunearmotifs)$name = munearmotifswide$names
  grmunearmotifsnorepeats<- unique(grmunearmotifs)
  distsfrommu_norepeats=distanceToNearest(grmunearmotifsnorepeats,grmotifs,select="arbitrary")
  d<-mcols(distsfrommu_norepeats)$distance
  message("distances")
  message(length(d))
  thebins <-cut(d,seq(-1,1500, 10))
  bindf <- as.data.frame(cbind(d,thebins))
  n_occur_dis <- data.frame(table(bindf$thebins))
  n_occur_dis$Var1 <- as.integer(n_occur_dis$Var1)
  n_occur_dis <- n_occur_dis[order(n_occur_dis$Freq),] 
  q1=n_occur_dis[37,]$Freq
  q3=n_occur_dis[112,]$Freq
  q2median = n_occur_dis[75,]$Freq
  iqr = q3-q1
  altop= q2median+(iqr*2)
  albottom= q2median-(iqr*2)
  top= q2median+(iqr*10)
  bottom= q2median-(iqr*10)
  message("q1 ", q1)
  message("q2 ", q2median) 
  message("q3 ", q3)
  message("iqr ", iqr)
  n_occur_dis_bypos <-n_occur_dis[order(n_occur_dis$Var1),]
  freqplotdf <- t(n_occur_dis_bypos)
  freqplotdf = freqplotdf["Freq",]
  col_fun = colorRamp2(c(bottom,albottom, q1,q2median,q3,altop, top), c("red", "#feeeed","#fff8f6", "white", "#f6fff7","#edfee7" ,"green"))
  a <- ComplexHeatmap::Heatmap(t(freqplotdf), col=col_fun, cluster_columns = FALSE,row_title=thissamplename, name=paste(thissamplename, " Color based on IQR", sep=""))
  return(a)}

#mufile="/Shares/dbnascent/Allen2014global/tfit/SRR1105736.sorted-1_bidir_predictions.bed"

mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105736.sorted-1_bidir_predictions.bed"
p53_DMSO_rep1 <- createHeatmapIQR(mufile, motiffile, "p53_DMSO_rep1")
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105737.sorted-1_bidir_predictions.bed"
p53_DMSO_rep2 <- createHeatmapIQR(mufile, motiffile, "p53_DMSO_rep2")
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105738.sorted-1_bidir_predictions.bed"
p53_Nutlin_rep1 <- createHeatmapIQR(mufile, motiffile, "p53_Nutlin_rep1")
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105739.sorted-1_bidir_predictions.bed"
p53_Nutlin_rep2 <- createHeatmapIQR(mufile, motiffile, "p53_Nutlin_rep2")
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105740.sorted-1_bidir_predictions.bed"
p53null_DMSO_rep1 <- createHeatmapIQR(mufile, motiffile, "p53null_DMSO_rep1")
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105741.sorted-1_bidir_predictions.bed"
p53null_Nutlin_rep1 <- createHeatmapIQR(mufile, motiffile, "p53null_Nutlin_rep1")
ht_list = p53_DMSO_rep1%v% p53_Nutlin_rep1 %v% p53_DMSO_rep2 %v% p53_Nutlin_rep2 %v% p53null_DMSO_rep1 %v% p53null_Nutlin_rep1
draw(ht_list)


createHeatmapzscoreonIQR <- function(mufile, motiffile, thissamplename){
  usecentermotif==TRUE
  #read the mu file
  mus = read.table(mufile)
  message("number of mus")
  message(dim(mus))
  #if I'm using the mu-merge I have to make up my own names for the regions
  if (length(colnames(mus))==3){mus$name = paste(mus$V1,":",  mus$V2,"-",mus$V3, sep="")}
  grmus = GRanges(seqnames = mus[[1]], ranges = IRanges(mus[[2]], mus[[3]]), names=mus[[4]])
  #grab the center of the mu file
  grmus<-mutate(anchor_center(grmus), width = 1)
  length(grmus)
  #import the motifs
  motifs = read.table(motiffile)
  message("number of motifs")
  message(dim(motifs))
  motifs["loc"] = paste(motifs$V1,":",  motifs$V2,"-",motifs$V3, sep="")
  grmotifs <- GRanges(seqnames = motifs[[1]], ranges = IRanges(motifs[[2]], motifs[[3]]), score = motifs[[5]], name= motifs[[4]],  strand=motifs[[6]], loc=motifs[,"loc"])
  #grab the center of the motifs and Tfit calls and expand the motif regions by 1500 on both sides
  if(usecentermotif==TRUE){
    grmotifs=  mutate(anchor_center(grmotifs), width = 1)}
  #mu near motifs
  widegrmotifs=  mutate(anchor_center(grmotifs), width = windowaroundmu*2)
  munearmotifswide <- join_overlap_inner(grmus,widegrmotifs) #This will give repeat the mu if the mu is by two motifs
  grmunearmotifs <-granges(munearmotifswide) #removes the subject so you can use the queary ranges again
  mcols(grmunearmotifs)$name = munearmotifswide$names
  grmunearmotifsnorepeats<- unique(grmunearmotifs)
  distsfrommu_norepeats=distanceToNearest(grmunearmotifsnorepeats,grmotifs,select="arbitrary")
  d<-mcols(distsfrommu_norepeats)$distance
  message("distances")
  message(length(d))
  thebins <-cut(d,seq(-1,1500, 10))
  bindf <- as.data.frame(cbind(d,thebins))
  n_occur_dis <- data.frame(table(bindf$thebins))
  n_occur_dis$Var1 <- as.integer(n_occur_dis$Var1)
  n_occur_dis <- n_occur_dis[order(n_occur_dis$Freq),] 
  q1=n_occur_dis[37,]$Freq
  q3=n_occur_dis[112,]$Freq
  q2median = n_occur_dis[75,]$Freq
  iqr = q3-q1
  altop= q2median+(iqr*2)
  albottom= q2median-(iqr*2)
  top= q2median+(iqr*10)
  bottom= q2median-(iqr*10)
  message("q1 ", q1)
  message("q2 ", q2median) 
  message("q3 ", q3)
  message("IQR ", iqr)
  n_occur_dis$Freq_IQR<-ifelse((n_occur_dis$Freq>q3|n_occur_dis$Freq<q1),NA,n_occur_dis$Freq)
  n_occur_dis <- n_occur_dis %>% mutate(zscore_basedonmid = (Freq - mean(Freq_IQR,na.rm=TRUE))/sd(Freq_IQR,na.rm=TRUE))
  n_occur_dis_bypos <-n_occur_dis[order(n_occur_dis$Var1),]
  freqplotdf <- t(n_occur_dis_bypos)
  freqplotdf = freqplotdf["zscore_basedonmid",]
  col_fun = colorRamp2(c(-30,-6,-3,0,3,6,30), c("red", "#feeeed","#fff8f6", "white", "#f6fff7","#edfee7" ,"green"))
  a<- ComplexHeatmap::Heatmap(t(freqplotdf), cluster_columns = FALSE,col=col_fun, name="Zscore based on IQR remove background",left_annotation = rowAnnotation(text = anno_text(thissamplename, just = "right", location=1)))
  #a <- ComplexHeatmap::Heatmap(t(freqplotdf), , cluster_columns = FALSE)
  return(a)}


createHeatmapzscore_afteroutlierremoval <- function(mufile, motiffile, thissamplename){
  usecentermotif==TRUE
  #read the mu file
  mus = read.table(mufile)
  message("number of mus")
  message(dim(mus))
  #if I'm using the mu-merge I have to make up my own names for the regions
  if (length(colnames(mus))==3){mus$name = paste(mus$V1,":",  mus$V2,"-",mus$V3, sep="")}
  grmus = GRanges(seqnames = mus[[1]], ranges = IRanges(mus[[2]], mus[[3]]), names=mus[[4]])
  #grab the center of the mu file
  grmus<-mutate(anchor_center(grmus), width = 1)
  length(grmus)
  #import the motifs
  motifs = read.table(motiffile)
  message("number of motifs")
  message(dim(motifs))
  motifs["loc"] = paste(motifs$V1,":",  motifs$V2,"-",motifs$V3, sep="")
  grmotifs <- GRanges(seqnames = motifs[[1]], ranges = IRanges(motifs[[2]], motifs[[3]]), score = motifs[[5]], name= motifs[[4]],  strand=motifs[[6]], loc=motifs[,"loc"])
  #grab the center of the motifs and Tfit calls and expand the motif regions by 1500 on both sides
  if(usecentermotif==TRUE){
    grmotifs=  mutate(anchor_center(grmotifs), width = 1)}
  #mu near motifs
  widegrmotifs=  mutate(anchor_center(grmotifs), width = windowaroundmu*2)
  munearmotifswide <- join_overlap_inner(grmus,widegrmotifs) #This will give repeat the mu if the mu is by two motifs
  grmunearmotifs <-granges(munearmotifswide) #removes the subject so you can use the queary ranges again
  mcols(grmunearmotifs)$name = munearmotifswide$names
  grmunearmotifsnorepeats<- unique(grmunearmotifs)
  distsfrommu_norepeats=distanceToNearest(grmunearmotifsnorepeats,grmotifs,select="arbitrary")
  d<-mcols(distsfrommu_norepeats)$distance
  message("distances")
  message(length(d))
  thebins <-cut(d,seq(-1,1500, 10))
  bindf <- as.data.frame(cbind(d,thebins))
  n_occur_dis <- data.frame(table(bindf$thebins))
  n_occur_dis$Var1 <- as.integer(n_occur_dis$Var1)
  n_occur_dis <- n_occur_dis[order(n_occur_dis$Freq),] 
  q1=n_occur_dis[37,]$Freq
  q3=n_occur_dis[112,]$Freq
  q2median = n_occur_dis[75,]$Freq
  iqr = q3-q1
  altop= q2median+(iqr*2)
  albottom= q2median-(iqr*2)
  top= q2median+(iqr*10)
  bottom= q2median-(iqr*10)
  message("q1 ", q1)
  message("q2 ", q2median) 
  message("q3 ", q3)
  message("IQR ", iqr)
  message("removing below",q1-(iqr*1.5))
  message("removing above",q3+(iqr*1.5))
  n_occur_dis$Freq_IQR<-ifelse((n_occur_dis$Freq>q3+(iqr*1.5)|n_occur_dis$Freq<q1-(iqr*1.5)),NA,n_occur_dis$Freq)
  message("calculaing zscores")
  n_occur_dis <- n_occur_dis %>% mutate(zscore_basedonmid = (Freq - mean(Freq_IQR,na.rm=TRUE))/sd(Freq_IQR,na.rm=TRUE))
  n_occur_dis_bypos <-n_occur_dis[order(n_occur_dis$Var1),]
  freqplotdf <- t(n_occur_dis_bypos)
  freqplotdf = freqplotdf["zscore_basedonmid",]
  col_fun = colorRamp2(c(-30,-6,-3,0,3,6,30), c("red", "#feeeed","#fff8f6", "white", "#f6fff7","#edfee7" ,"green"))
  a<- ComplexHeatmap::Heatmap(t(freqplotdf), cluster_columns = FALSE,col=col_fun, name="Zscore based on background (outliers removed)",left_annotation = rowAnnotation(text = anno_text(thissamplename, just = "right", location=1)))
  #a <- ComplexHeatmap::Heatmap(t(freqplotdf), , cluster_columns = FALSE)
  return(a)}

mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105736.sorted-1_bidir_predictions.bed"
p53_DMSO_rep1 <- createHeatmapzscore_afteroutlierremoval(mufile, motiffile, "p53_DMSO_rep1")
p53_DMSO_rep1
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105737.sorted-1_bidir_predictions.bed"
p53_DMSO_rep2 <- createHeatmapzscore_afteroutlierremoval(mufile, motiffile, "p53_DMSO_rep2")
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105738.sorted-1_bidir_predictions.bed"
p53_Nutlin_rep1 <- createHeatmapzscore_afteroutlierremoval(mufile, motiffile, "p53_Nutlin_rep1")
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105739.sorted-1_bidir_predictions.bed"
p53_Nutlin_rep2 <- createHeatmapzscore_afteroutlierremoval(mufile, motiffile, "p53_Nutlin_rep2")
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105740.sorted-1_bidir_predictions.bed"
p53null_DMSO_rep1 <- createHeatmapzscore_afteroutlierremoval(mufile, motiffile, "p53null_DMSO_rep1")
mufile = "/Shares/dbnascent/Allen2014global/tfit/SRR1105741.sorted-1_bidir_predictions.bed"
p53null_Nutlin_rep1 <- createHeatmapzscore_afteroutlierremoval(mufile, motiffile, "p53null_Nutlin_rep1")
ht_list = p53_DMSO_rep1%v% p53_Nutlin_rep1 %v% p53_DMSO_rep2 %v% p53_Nutlin_rep2 %v% p53null_DMSO_rep1 %v% p53null_Nutlin_rep1
draw(ht_list)

#left_annotation
#cn = colnames(mat)
#Heatmap(mat, show_column_names = FALSE, 
#        bottom_annotation = HeatmapAnnotation(
#          text = anno_text(thissamplename, rot = 90, offset = unit(1, "npc"), just = "right"),
#
#        )
#)

# left_annotation= HeatmapAnnotation(
#          text = anno_text(cn, rot = 90, offset = unit(1, "npc"), just = "right"),
#          annotation_height = max_text_width(cn)


n_occur_dis <- n_occur_dis %>% mutate(Freq_removeback = (Freq - mean(Freq_IQR,na.rm=TRUE)))
n_occur_dis_bypos <-n_occur_dis[order(n_occur_dis$Var1),]



