Class 11: Structural Bioinformatics
================
11/6/2018

1.1 PDB Statistics
==================

Download PDB statictics as CSV file from <http://www.rcsb.org/stats/summary>

``` r
pdbstats <- read.csv("Data Export Summary.csv", row.names = 1)
```

lets look at the table

``` r
library(knitr)
kable(pdbstats)
```

|                     |  Proteins|  Nucleic.Acids|  Protein.NA.Complex|  Other|   Total|
|---------------------|---------:|--------------:|-------------------:|------:|-------:|
| X-Ray               |    122263|           1960|                6333|     10|  130566|
| NMR                 |     10898|           1263|                 253|      8|   12422|
| Electron Microscopy |      1822|             31|                 657|      0|    2510|
| Other               |       244|              4|                   6|     13|     267|
| Multi Method        |       119|              5|                   2|      1|     127|

Q1: Determine the percentage of structures solved by X-Ray and Electron Microscopy. What proportion of structures are protein?

``` r
sumTotal <- sum(pdbstats$Total)
percentage <- round(pdbstats$Total / sumTotal*100, 2)
```

The percentage of structures solved by X-ray is 89.49% and the perentage of structures solved by electron microscopy is 1.72% in the protein database as 2018-11-06

adding the percentage of each method as a new column of the table

``` r
nstats <- pdbstats
nstats$Percent <- percentage
kable(nstats)
```

|                     |  Proteins|  Nucleic.Acids|  Protein.NA.Complex|  Other|   Total|  Percent|
|---------------------|---------:|--------------:|-------------------:|------:|-------:|--------:|
| X-Ray               |    122263|           1960|                6333|     10|  130566|    89.49|
| NMR                 |     10898|           1263|                 253|      8|   12422|     8.51|
| Electron Microscopy |      1822|             31|                 657|      0|    2510|     1.72|
| Other               |       244|              4|                   6|     13|     267|     0.18|
| Multi Method        |       119|              5|                   2|      1|     127|     0.09|

Calculate the proportion of protein entries

``` r
protein <- sum(pdbstats$Proteins)
proportionP <- round(protein / sumTotal*100, 2)
```

The proportion of protein entries from the summary table is 92.77%
