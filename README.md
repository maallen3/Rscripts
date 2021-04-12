# Rscripts
These are small R scripts that show one how to use Deseq. 

There are two examples of how to count using feature counts:
```
basic_count_htseq.R
count_htseq_oct2020.R
```
Then there are three examples of runing Deseq2 with counts files. 
```
basic_Deseq2.R #goes over a basic run
basic_Deseq2_LRT.R #goes over a basic LRT run
basic_Deseq2_withbatcheffects.R #adds batch effects to the design
```


# Running the scripts

In order to run the scripts you will need the Deseq2 library. 

You can get the file of counts used for the basic deseq scripts by 


example:
```
wget https://www.colorado.edu/lab/allen/sites/default/files/attached-files/featurecounts_attr_gene_id_feature_exon_125256.coverage.csv_.txt
wget https://www.colorado.edu/lab/allen/sites/default/files/attached-files/featurecounts_body_153616.coverage.csv_.txt
```

These files are csv, so you should change the name removing the _.txt

example:

```
mv featurecounts_attr_gene_id_feature_exon_125256.coverage.csv_.txt featurecounts_attr_gene_id_feature_exon_125256.coverage.csv
mv featurecounts_body_153616.coverage.csv_.txt featurecounts_body_153616.coverage.csv
```
# Looking under the hood
What is 
```
DEdds <- DESeq(dds)
```
really doing?

These are the steps underhood and the slots they add to the dds object
```
dds <-estimateSizeFactors(dds) #adds a column to colData(dds) called sizeFactor if you use the parmater normmatrix then normalizationFactors(dds) instead.
dds <- estimateDispersionsGeneEst(dds) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(dds) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
dds<-estimateDispersionsFit(dds) # and fills in dispersionFunction(dds) and adds a column in  elementMetadata(dds) called dispFit
dds <- estimateDispersionsMAP(dds) #adds a attr(,"dispPriorVar") to dispersionFunction(dds) 
#and adds 4 columns to elementMetadata(dds): dispersion, dispIter, dispOutlier, dispMAP
dds <- nbinomWaldTest(dds) #adds both H and cooks to assays(dds)  both are for each sample you have.
#Adds a ton of columns to elementMetadata(dds) 
#including Intercept and all values they expect you to compare 
#SE_Intercept (SE stands for standard error)  and SE for every thing they expect you to compare. 
#Wald static for what they expect you to compare. 
#betaConv  betaIter         deviance  maxCooks
```


# design

Lets say you have two genotypes (D21 and T21)
And one treatment (tempature)


to get effect of treatment in each group seperately

dds$group <- factor(paste0(dds$genotype, "_", dds$treatment))

if design = ~group
and
res = (dds, contrast("group", "D21_thirtyseven", "D21_fourtytwo")
you get the difference caused by the treatment in the D21 
res = (dds, contrast("group", "T21_thirtyseven", "T21_fourtytwo")
you get the difference caused by the treatment in the T21


If you only want the common features of the treatment (regardless of where the gene started in genotype 1 and 2)
design = ~genotype+treatment

if you want the effect of the treatment on the first level of genotype with the realtive responce of the other levels of genotype
~genotype + treatment+genotype:treatment

![Comparison of design](https://github.com/maallen3/Rscript/blob/main/compare_design.png)



