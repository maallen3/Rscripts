
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rsubread")
library("Rsubread")


fastqdir1="/Shares/down/RNA/hoeffer/nextflowhoeffermm10_072920a/mapped/bams/"
fastqdir2="/Shares/down/RNA/hoeffer/nextflowhoeffermm10_072920b/mapped/bams/"
outdir="/Shares/down/RNA/hoeffer/Oct2020/"
#gtf <- "/scratch/Shares/dowell/genomes/mm10/mm10_refseq.gtf"
gtf <- "/scratch/Shares/public/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"


annot_file <- gtf
filetable <- read.csv(paste0(outdir, "fileinfo.tech.csv"), header=TRUE)
filelist <-as.vector(filetable$filelist)

if (!all(file.exists(filelist))) {
  print("WARNING: Not all specified files exist")
}


#GTFattrType="gene"
#GTFattrType="Dbxref"
GTFattrType="gene_id"

coverage <- featureCounts(files=filelist,
                          annot.ext=annot_file,
                          isGTFAnnotationFile=TRUE,
                          useMetaFeatures=TRUE,
                          GTF.featureType="exon",
                          GTF.attrType=GTFattrType,
                          allowMultiOverlap=TRUE,
                          largestOverlap=TRUE,
                          isPairedEnd=TRUE,
                          requireBothEndsMapped=FALSE,
                          nthreads=32)

colnames(coverage$counts) <- filetable$label

time <- strsplit(as.character(Sys.time()), split = " ")[[1]][2]
time <- paste(strsplit(time, split = ":")[[1]], collapse = '')
time

save.image(paste0(outdir, "nextflowhoeffermm10_072920a_and_b_res_featureCounts_", GTFattrType, "full", "_", time, ".RData"))

