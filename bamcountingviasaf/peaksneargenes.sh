#!/bin/bash 
#SBATCH --job-name=peaksneargenes # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=23:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/e_and_o/peaksneargenes.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/e_and_o/peaksneargenes.%j.err # Standard error log



module load bedtools/2.28.0 



peakfile=/scratch/Users/joca4543/191120_Cardiello_Dowell-992/results/2020_01_15/MuMerge/2021_updated_ATACseqFullFiles_MUMERGE.bed
genefile=/scratch/Shares/dowell/genomes/hg38/hg38_refseq_genenames.bed
outfile1=/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/peaksneargenes/2021_updated_ATACseqFullFiles_MUMERGE_neargenes.bed
outfile2=/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/peaksneargenes/genesnear_2021_updated_ATACseqFullFiles_MUMERGE.bed

bedtools window -w 25000 -u -a $peakfile -b $genefile >$outfile1
bedtools window -w 25000 -u -a $genefile -b $peakfile >$outfile2 

wc -l $peakfile
wc -l $genefile
wc -l $outfile1
wc -l $outfile2
