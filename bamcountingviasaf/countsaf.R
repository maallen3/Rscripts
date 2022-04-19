library("Rsubread")

metadatadir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/scripts/"
metadatafile="atacbammeta.csv"
metadata= paste(metadatadir, metadatafile, sep="")
outdir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/counts/"
indir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/counts/"
annot_file = paste(indir, "hg38_refseq_genenames_and_2021_updated_ATACseqFullFiles_MUMERGE.saf", sep="")

safdf = read.table(annot_file, sep=",", header=TRUE)
head(df)
tail(df)



filetable <- read.csv(metadata, header=TRUE, sep=",")
filetable$bamfile <- paste(filetable$directory, filetable$filename, sep="")
filetable$label <- paste(filetable$protocol, filetable$D21orT21, filetable$replicate, filetable$tempature, sep="_")
#write.csv(filetable, metadata, quote=FALSE, row.names = FALSE)
filelist <-as.vector(filetable$bamfile)

#check if all the bam files exist
if (!all(file.exists(filelist))) {
  print("WARNING: Not all specified files exist")
}

coverage <- featureCounts(files=filelist,
                          annot.ext=safdf,
                          isGTFAnnotationFile=FALSE,
                          isPairedEnd=TRUE,
                          nthreads=32)

colnames(coverage$counts) <- filetable$label

time <- strsplit(as.character(Sys.time()), split = " ")[[1]][2]
time <- paste(strsplit(time, split = ":")[[1]], collapse = '')
time

#you can save the whole session as a R image. I would not suggest it. 
#save.image(paste0(outdir, "nextflowhoeffermm10_072920a_and_b_res_featureCounts_", GTFattrType, "full", "_", time, ".RData"))


fileroot<-paste0(outdir, "ATAC_Peaks_and_genes","_", time)
#fileroot<-paste0(outdir, "featureCounts_tss_", time)

save.image(paste0(fileroot, ".RData"))

write.csv(coverage$counts, paste(fileroot,".coverage.csv", sep=""))
write.csv(coverage$stat, paste(fileroot,".stat.csv", sep=""))
write.csv(coverage$annotation, paste(fileroot,".annotation.csv", sep=""))
write.csv(coverage$targets, paste(fileroot,".targets.csv", sep=""))
