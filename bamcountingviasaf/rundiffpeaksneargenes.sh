



#need unsortedbedfile, peakfile, genesnearpeaksfile
beddir=/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/deseq2/peakandgene_controlsfwithgenes_ignoregenotype/
#beddir=/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/deseq2/peakonly_ignoregenotype/
#beddir=/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/deseq2/peakonly/
#beddir=/Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/deseq2/peakandgene_controlsfwithgenes/
outdir=$beddir

for pathandfilename in `ls ${beddir}*.unsorted.bed`; do
rootname=`basename $pathandfilename .unsorted.bed`
peakfile=${outdir}${rootname}.sorted.bed
genesnearpeaksfile=${outdir}genesnear_${rootname}.bed
echo $rootname
sbatch --export=unsortedbedfile=$pathandfilename,peakfile=$peakfile,genesnearpeaksfile=$genesnearpeaksfile /Shares/down/heatshock/analysis/ATAC_peaks_diff_exp/scripts/diffpeaksneargenes.sh
done

