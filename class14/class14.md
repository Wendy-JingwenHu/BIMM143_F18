Class 14: Genomic Analysis
================
11/15/2018

Asthma SNPs
===========

Examine the asthma SNPs in the MXL data
---------------------------------------

``` r
SNP <- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
```

Lets focus on the second column of the genotype info to see the percentage of each SNP

``` r
genotypes <- round (table(SNP[,2])/nrow(SNP)*100, 2)
```

There are 34.38% AA genotypes in the MXL population

### Interpret base qualities in R

``` r
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred
```

    ##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
    ## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31

Population genetics analysis on SNP
===================================

``` r
expression <- read.table("https://bioboot.github.io/bimm143_F18/class-material/rs8067378_ENSG00000172057.6.txt")
inds.aa <- expression$geno == "A/A"
summary(expression$exp[inds.aa])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   11.40   27.02   31.25   31.82   35.92   51.52

``` r
inds.ag <- expression$geno == "A/G"
summary(expression$exp[inds.ag])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   7.075  20.626  25.065  25.397  30.552  48.034

``` r
inds.gg <- expression$geno == "G/G"
summary(expression$exp[inds.gg])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   6.675  16.903  20.074  20.594  24.457  33.956

``` r
boxplot(exp ~ geno, data = expression)
```

![](class14_files/figure-markdown_github/unnamed-chunk-4-1.png)

Start from fastq file (the input), do quality control (FASTQC), alignment (mapping) using TopHat, then do Cufflink to count the entries, get the table in R, make a boxplot to show the gene expression is different
