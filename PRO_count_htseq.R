
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rsubread")
library("Rsubread")
part="tss"

#set your out directory and the gtf for your genome
outdir="/Shares/down/heatshock/analysis/PROseq/genedeseq3/"
#annot_file <- "/scratch/Shares/public/genomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
annot_file <-paste("/scratch/Shares/dowell/genomes/hg38/PROparts_genenames/",part,".saf", sep="")


#Make a list of your bam files. 
#You can make the list in a text file with all the metadata (you will need this for deseq2 anyway), or you can do that in the R script.
filetable <- read.csv("/Users/allenma/scripts/Rscripts/Rscripts/PROmeta.txt", header=TRUE, sep="\t")
filelist <-as.vector(filetable$bamfile)

#check if all the bam files exist
if (!all(file.exists(filelist))) {
  print("WARNING: Not all specified files exist")
}

coverage <- featureCounts(files=filelist,
                          annot.ext=annot_file,
                          isGTFAnnotationFile=FALSE,
                          isPairedEnd=FALSE,
                          nthreads=32)

colnames(coverage$counts) <- filetable$label

time <- strsplit(as.character(Sys.time()), split = " ")[[1]][2]
time <- paste(strsplit(time, split = ":")[[1]], collapse = '')
time

#you can save the whole session as a R image. I would not suggest it. 
#save.image(paste0(outdir, "nextflowhoeffermm10_072920a_and_b_res_featureCounts_", GTFattrType, "full", "_", time, ".RData"))


fileroot<-paste0(outdir, "featureCounts_",part,"_", time)
#fileroot<-paste0(outdir, "featureCounts_tss_", time)

save.image(paste0(fileroot, ".RData"))

write.csv(coverage$counts, paste(fileroot,".coverage.csv", sep=""))
write.csv(coverage$stat, paste(fileroot,".stat.csv", sep=""))
write.csv(coverage$annotation, paste(fileroot,".annotation.csv", sep=""))
write.csv(coverage$targets, paste(fileroot,".targets.csv", sep=""))



