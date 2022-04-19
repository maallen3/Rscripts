library(dplyr)

outdir="/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/counts/"

#read in the bed file which has three columns
bedfile="/scratch/Users/joca4543/191120_Cardiello_Dowell-992/results/2020_01_15/MuMerge/2021_updated_ATACseqFullFiles_MUMERGE.bed"
df = read.table(bedfile, sep="\t", col.name=c("Chr", "Start", "End"))

#add a geneID and strand column. In this case strand doesn't matter so I use a .
df$GeneID = paste(df$Chr, ":", df$Start, "-", df$End, sep="")
df$Strand = "."


#reanage the columns to standard saf format
df <- df %>% dplyr::select(GeneID, Chr, Start, End, Strand)

#output the saf file
write.csv(df, paste(outdir, "2021_updated_ATACseqFullFiles_MUMERGE.saf", sep=""), quote = FALSE, row.names = FALSE)
