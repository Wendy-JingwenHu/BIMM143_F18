Class 8
================
10/25/2018

K-means clustering
------------------

Our first example with **kmeans()** function.

``` r
# Example plot to see how Rmarkdown works
# Shortcut for insert code chunks: option+command+i
plot(1:10, typ = "l")
```

![](class08_files/figure-markdown_github/unnamed-chunk-1-1.png)

Back to kmeans

``` r
# Generate some example data for clustering
tmp <- c(rnorm(30,-3), rnorm(30,3))  # this is a vector, c is used to concatenate
x <- cbind(x=tmp, y=rev(tmp))
plot(x)
```

![](class08_files/figure-markdown_github/unnamed-chunk-2-1.png)

Use the kmeans() function setting k to 2 and nstart=20

``` r
k <- kmeans(x, centers = 2, nstart = 20)
```

Inspect/print the results

``` r
print(k)
```

    ## K-means clustering with 2 clusters of sizes 30, 30
    ## 
    ## Cluster means:
    ##           x         y
    ## 1  3.124985 -2.990950
    ## 2 -2.990950  3.124985
    ## 
    ## Clustering vector:
    ##  [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1
    ## [36] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## 
    ## Within cluster sum of squares by cluster:
    ## [1] 48.93355 48.93355
    ##  (between_SS / total_SS =  92.0 %)
    ## 
    ## Available components:
    ## 
    ## [1] "cluster"      "centers"      "totss"        "withinss"    
    ## [5] "tot.withinss" "betweenss"    "size"         "iter"        
    ## [9] "ifault"

Q. How many points are in each cluster?

``` r
k$size
```

    ## [1] 30 30

-   there are 30 points in each cluster

Q. What ‘component’ of your result object details - cluster size: `k$size` - cluster assignment/membership: clustering vector

``` r
table(k$cluster)
```

    ## 
    ##  1  2 
    ## 30 30

    - cluster center: values in cluster means

``` r
k$centers
```

    ##           x         y
    ## 1  3.124985 -2.990950
    ## 2 -2.990950  3.124985

    - Within cluster sum of squares by cluster: square distance to the center

Plot x colored by the kmeans cluster assignment and add cluster centers as blue points

``` r
palette(c("blue", "black"))
plot(x, col = k$cluster)
points(k$centers, col = "red", pch = 20, cex = 3)
```

![](class08_files/figure-markdown_github/unnamed-chunk-8-1.png)

Q. Repeat for k=3, which has the lower tot.withinss?

``` r
k3 <- kmeans(x, centers = 3, nstart = 20)
k4 <- kmeans(x, centers = 4, nstart = 20)
k10 <- kmeans(x, centers = 10, nstart = 20)
k10$tot.withinss
```

    ## [1] 21.27397

``` r
k4$tot.withinss
```

    ## [1] 62.23547

``` r
k3$tot.withinss
```

    ## [1] 80.05128

``` r
k$tot.withinss
```

    ## [1] 97.8671

``` r
plot(c(k$tot.withinss, k3$tot.withinss, k4$tot.withinss, k10$tot.withinss))
```

![](class08_files/figure-markdown_github/unnamed-chunk-9-1.png)

-   k=3 has a lower tot.withinss (larger k will result lower tot.withinss) Major disadvantage is that we have to give a k value each time.

Hierarchical clustering in R
============================

Let's try out the **hclust()** function for Hierarchical clustering in R. This function need a **distance matrix** as input.

``` r
d <- dist(x)
hc <- hclust(d)
#this plots a dendrogram 
plot(hc)
# Draw a line on the dendrogram
abline(h = 8, col = "red")
```

![](class08_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# Cut the tree to yeild cluster membership vector
cutree(hc, h = 8)  #cut by height
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
cutree(hc, k = 2)  #cut by k groups
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

Using different hierarchical clustering methods

``` r
hc.complete <- hclust(d, method="complete")
hc.average <- hclust(d, method="average")
hc.single <- hclust(d, method="single")
```

### A more real example of hierarchical clustering

``` r
### Step 1. Generate some example data for clustering
x <- rbind(
matrix(rnorm(100, mean=0, sd = 0.3), ncol = 2), # c1
matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2), # c2
matrix(c(rnorm(50, mean = 1, sd = 0.3), # c3
rnorm(50, mean = 0, sd = 0.3)), ncol = 2))
colnames(x) <- c("x", "y")
# Step 2. Plot the data without clustering
plot(x)
```

![](class08_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
# Step 3. Generate colors for known clusters
# (just so we can compare to hclust results)
palette(c("red","blue","black"))
col <- as.factor( rep(c("c1","c2","c3"), each=50) )
plot(x, col=col, pch = 20)
```

![](class08_files/figure-markdown_github/unnamed-chunk-12-2.png)

Q. Use the dist(), hclust(), plot() and cutree() functions to return 2 and 3 clusters

``` r
disX <- dist(x)
hcx <- hclust(disX)
plot(hcx)
```

![](class08_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
cluster2 <- cutree(hcx, k =2)  #return 2 clusters
cluster3 <- cutree(hcx, k =3)  #return 3 clusters
palette(c("red","blue", "yellow"))
plot(x, col = cluster2, pch = 20)
```

![](class08_files/figure-markdown_github/unnamed-chunk-13-2.png)

``` r
plot(x, col = cluster3, pch = 20)
```

![](class08_files/figure-markdown_github/unnamed-chunk-13-3.png)

Q. How does this compare to your known 'col' groups?

``` r
table(cluster3, col)
```

    ##         col
    ## cluster3 c1 c2 c3
    ##        1 40  0  0
    ##        2  6 48  4
    ##        3  4  2 46

Using principal component analysis (PCA) in R
=============================================

Generate a gene expression table

``` r
# Initialize a blank 100 row by 10 column matrix
mydata <- matrix(nrow=100, ncol=10)
# Lets label the rows gene1, gene2 etc. to gene100
rownames(mydata) <- paste("gene", 1:100, sep="")
# Lets label the first 5 columns from wt1 to wt5 and the last 5 fromko1 to ko5
colnames(mydata) <- c( paste("wt", 1:5, sep=""),
paste("ko", 1:5, sep="") )
# Fill in some fake read counts
for(i in 1:nrow(mydata)) {
wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
mydata[i,] <- c(wt.values, ko.values)
}
head(mydata)
```

    ##       wt1 wt2 wt3 wt4 wt5  ko1 ko2 ko3  ko4  ko5
    ## gene1 181 155 177 170 218  656 676 650  651  708
    ## gene2 660 637 655 705 666 1006 984 973 1021 1029
    ## gene3 651 725 642 656 661  220 209 239  245  218
    ## gene4 714 620 659 644 683  742 731 727  777  763
    ## gene5 730 724 762 745 725   60  54  57   65   63
    ## gene6 742 742 737 742 723  559 568 566  549  587

Note the prcomp() function wants us to take the transposed data of our matrix with column and rows flipped

``` r
head(t(mydata))
```

    ##     gene1 gene2 gene3 gene4 gene5 gene6 gene7 gene8 gene9 gene10 gene11
    ## wt1   181   660   651   714   730   742   337   200   912    225    140
    ## wt2   155   637   725   620   724   742   334   214   850    224    146
    ## wt3   177   655   642   659   762   737   340   200   893    220    155
    ## wt4   170   705   656   644   745   742   317   185   853    243    147
    ## wt5   218   666   661   683   725   723   318   230   809    204    129
    ## ko1   656  1006   220   742    60   559   571   306   665    945    646
    ##     gene12 gene13 gene14 gene15 gene16 gene17 gene18 gene19 gene20 gene21
    ## wt1    739    212    964    577    895    135    806    442    156    557
    ## wt2    642    201   1007    594    891    146    792    432    143    533
    ## wt3    717    224   1016    543    927    165    777    399    166    553
    ## wt4    667    218    987    572    904    165    813    416    155    541
    ## wt5    669    216    968    557    927    147    774    441    172    590
    ## ko1    203    172    621    112    536    359    234    404    863    543
    ##     gene22 gene23 gene24 gene25 gene26 gene27 gene28 gene29 gene30 gene31
    ## wt1    327    272    834    715     43     66    389    593    827    611
    ## wt2    324    244    817    651     62     92    420    596    812    601
    ## wt3    311    261    779    657     48     71    411    544    788    650
    ## wt4    332    272    793    673     50     83    420    574    788    593
    ## wt5    309    260    823    708     52     74    439    580    827    666
    ## ko1    409    864    708    391    728    518    334    438     82    780
    ##     gene32 gene33 gene34 gene35 gene36 gene37 gene38 gene39 gene40 gene41
    ## wt1    303    360   1005    951    807    883    677    965    344    698
    ## wt2    337    346    986    948    802    890    638   1002    331    643
    ## wt3    310    370   1005    968    810    919    672    971    350    676
    ## wt4    334    336   1039    962    760    885    652    935    338    654
    ## wt5    363    351    997    973    828    948    644    928    327    688
    ## ko1     61    672    236    882    590    491    930    841     83    732
    ##     gene42 gene43 gene44 gene45 gene46 gene47 gene48 gene49 gene50 gene51
    ## wt1    237    611    633    334    762    953    973    471    210    624
    ## wt2    243    570    601    363    803    961    935    437    208    621
    ## wt3    254    616    632    334    833    899    949    419    200    569
    ## wt4    272    607    613    330    808    850    960    483    201    603
    ## wt5    241    584    614    321    748    834    978    472    184    585
    ## ko1    479    889    559    341    129    486    336    769    608    488
    ##     gene52 gene53 gene54 gene55 gene56 gene57 gene58 gene59 gene60 gene61
    ## wt1    997    854    462    813    697    172    330    617    907    286
    ## wt2    908    851    449    786    750    158    315    692    861    243
    ## wt3    892    849    489    799    696    184    332    689    791    274
    ## wt4    968    879    454    802    733    186    306    673    870    245
    ## wt5    878    867    491    841    670    176    337    638    893    286
    ## ko1    579    563    243    240    969    193    308    534    501     91
    ##     gene62 gene63 gene64 gene65 gene66 gene67 gene68 gene69 gene70 gene71
    ## wt1    653    718    691    657    226    737    451    433    161    469
    ## wt2    662    740    669    705    235    695    451    515    126    464
    ## wt3    650    844    659    685    260    674    465    511    137    455
    ## wt4    656    787    714    680    248    649    427    460    138    422
    ## wt5    649    727    712    686    228    709    381    433    139    425
    ## ko1    489     14     83    695    956    792    165    121    526    181
    ##     gene72 gene73 gene74 gene75 gene76 gene77 gene78 gene79 gene80 gene81
    ## wt1    598    577    254    159    609    200    439    436    424    472
    ## wt2    600    623    229    144    726    219    442    431    450    449
    ## wt3    572    622    216    130    668    187    460    412    433    429
    ## wt4    607    688    217    137    708    212    468    420    474    444
    ## wt5    585    623    207    135    669    217    431    457    456    427
    ## ko1    820    691    976    833    342    710    256    378    847     79
    ##     gene82 gene83 gene84 gene85 gene86 gene87 gene88 gene89 gene90 gene91
    ## wt1    408     46    579    371    777    248    843    679    572    804
    ## wt2    400     38    600    361    792    249    826    700    585    777
    ## wt3    425     36    572    379    790    239    780    694    593    807
    ## wt4    404     45    611    400    838    238    837    687    547    810
    ## wt5    420     47    594    372    807    248    806    667    545    829
    ## ko1    731    355    193    102    353    171    445    336    646    493
    ##     gene92 gene93 gene94 gene95 gene96 gene97 gene98 gene99 gene100
    ## wt1    500    833    423    396    620    692    718    729     550
    ## wt2    502    834    415    377    628    680    705    713     548
    ## wt3    455    848    371    412    633    651    710    716     529
    ## wt4    464    853    404    364    591    642    650    708     523
    ## wt5    490    855    410    397    631    633    695    701     529
    ## ko1    809     84    477    131    233    569    408     82     423

Now let's try find structure in this data with **prcomp()**

``` r
pca <- prcomp(t(mydata), scale = TRUE)
#see what's returned by prcomp using attributes()
attributes(pca)
```

    ## $names
    ## [1] "sdev"     "rotation" "center"   "scale"    "x"       
    ## 
    ## $class
    ## [1] "prcomp"

Make a PC plot of PC1 vs PC2, using `$x` component of output (pca$x)

``` r
plot(pca$x[,1], pca$x[,2])
```

![](class08_files/figure-markdown_github/unnamed-chunk-18-1.png)

Looking the variance from the original data that we are capturing

``` r
#variance: squared stdev
pca.var <- pca$sdev^2
pca.var
```

    ##  [1] 9.019482e+01 2.304601e+00 2.029525e+00 1.489292e+00 1.195600e+00
    ##  [6] 9.497404e-01 8.713764e-01 6.223631e-01 3.426834e-01 7.323252e-30

``` r
#getting the percentage variance 
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
pca.var.per
```

    ##  [1] 90.2  2.3  2.0  1.5  1.2  0.9  0.9  0.6  0.3  0.0

Making the basic PCA plot more useful

``` r
## A vector of colors for wt and ko samples
colvec <- colnames(mydata)
colvec[grep("wt", colvec)] <- "red"
colvec[grep("ko", colvec)] <- "blue"
plot(pca$x[,1], pca$x[,2], col=colvec, pch=16,
xlab=paste0("PC1 (", pca.var.per[1], "%)"),
ylab=paste0("PC2 (", pca.var.per[2], "%)"))
```

![](class08_files/figure-markdown_github/unnamed-chunk-20-1.png)

Making a "Scree-plot"

``` r
barplot(pca.var.per, main="Scree Plot",
xlab="Principal Component", ylab="Percent Variation")
```

![](class08_files/figure-markdown_github/unnamed-chunk-21-1.png)

Looking at the loading scores

``` r
# Lets focus on PC1 as it accounts for > 90% of variance
loading_scores <- pca$rotation[,1]
summary(loading_scores)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -0.10527 -0.10453 -0.09779 -0.02257  0.10282  0.10514

``` r
# We are interested in the magnitudes of both plus and minus contributing genes
gene_scores <- abs(loading_scores)
# Sort by magnitudes from high to low
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
# Find the names of the top 5 genes
top_5_genes <- names(gene_score_ranked[1:5])
# Show the scores (with +/- sign)
pca$rotation[top_5_genes,1]
```

    ##     gene99     gene93     gene30      gene5     gene34 
    ## -0.1052712 -0.1052250 -0.1052243 -0.1051922 -0.1051883

Try using what learned today to analyze UK food data
====================================================

First read the food data and its dimension

``` r
x <- read.csv("UK_foods.csv")
#reformat the row names
rownames(x) <- x[,1]
x <- x[,-1]
head(x)
```

    ##                England Wales Scotland N.Ireland
    ## Cheese             105   103      103        66
    ## Carcass_meat       245   227      242       267
    ## Other_meat         685   803      750       586
    ## Fish               147   160      122        93
    ## Fats_and_oils      193   235      184       209
    ## Sugars             156   175      147       139

``` r
#check the dimension of x 
dim(x)
```

    ## [1] 17  4

Read the whole data

``` r
knitr::kable(x, caption="The full UK foods data table")
```

|                     |  England|  Wales|  Scotland|  N.Ireland|
|---------------------|--------:|------:|---------:|----------:|
| Cheese              |      105|    103|       103|         66|
| Carcass\_meat       |      245|    227|       242|        267|
| Other\_meat         |      685|    803|       750|        586|
| Fish                |      147|    160|       122|         93|
| Fats\_and\_oils     |      193|    235|       184|        209|
| Sugars              |      156|    175|       147|        139|
| Fresh\_potatoes     |      720|    874|       566|       1033|
| Fresh\_Veg          |      253|    265|       171|        143|
| Other\_Veg          |      488|    570|       418|        355|
| Processed\_potatoes |      198|    203|       220|        187|
| Processed\_Veg      |      360|    365|       337|        334|
| Fresh\_fruit        |     1102|   1137|       957|        674|
| Cereals             |     1472|   1582|      1462|       1494|
| Beverages           |       57|     73|        53|         47|
| Soft\_drinks        |     1374|   1256|      1572|       1506|
| Alcoholic\_drinks   |      375|    475|       458|        135|
| Confectionery       |       54|     64|        62|         41|

Plot a heatmap from this data set

``` r
#par(mar=c(20, 4, 4, 2))
heatmap(as.matrix(x))
```

![](class08_files/figure-markdown_github/unnamed-chunk-25-1.png)

Functions covered today
=======================

kmeans(x) hclust(dist(x)) cutree(hc) prcomp(x, scale = TRUE) table()
