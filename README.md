# Rscripts
These are small R scripts that show one how to use Deseq. 

There are two examples of how to count using feature counts:
basic_count_htseq.R
count_htseq_oct2020.R

Then there are three examples of runing Deseq2 with counts files. 

basic_Deseq2.R #goes over a basic run
basic_Deseq2_LRT.R #goes over a basic LRT run
basic_Deseq2_withbatcheffects.R #adds batch effects to the design



# Running the scripts

In order to run the scripts you will need the Deseq2 library. 

You can get the file of counts used for the basic deseq scripts by 
example:
wget https://www.colorado.edu/lab/allen/sites/default/files/attached-files/featurecounts_attr_gene_id_feature_exon_125256.coverage.csv_.txt

This file is a csv, so you should change the name removing the _.txt

example:
mv featurecounts_attr_gene_id_feature_exon_125256.coverage.csv_.txt featurecounts_attr_gene_id_feature_exon_125256.coverage.csv


