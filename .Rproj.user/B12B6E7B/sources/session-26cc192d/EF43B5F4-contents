Lupus Data Wrangling
================

``` r
knitr::opts_chunk$set(echo = TRUE,message = FALSE, warning = FALSE)

library(Seurat)
```

    ## Attaching SeuratObject

``` r
library(tidyverse)
```

    ## Warning in system("timedatectl", intern = TRUE): running command 'timedatectl'
    ## had status 1

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.4.4     ✔ purrr   1.0.2
    ## ✔ tibble  3.2.1     ✔ dplyr   1.1.3
    ## ✔ tidyr   1.3.0     ✔ stringr 1.5.0
    ## ✔ readr   2.1.4     ✔ forcats 1.0.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(edgeR)
```

    ## Loading required package: limma

``` r
library(scater)
```

    ## Loading required package: SingleCellExperiment

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'SummarizedExperiment'

    ## The following object is masked from 'package:SeuratObject':
    ## 
    ##     Assays

    ## The following object is masked from 'package:Seurat':
    ## 
    ##     Assays

    ## 
    ## Attaching package: 'SingleCellExperiment'

    ## The following object is masked from 'package:edgeR':
    ## 
    ##     cpm

    ## Loading required package: scuttle

    ## 
    ## Attaching package: 'scater'

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMDS

``` r
source("/mnt/auxf.R")
```

## Read file

The dataset has been obtained from GEO accession number GSE174188

``` r
lupus <- readRDS("/domino/datasets/local/RUV/230126_lupusdataSCE_raw.rds")
lupus <- as.Seurat(lupus, counts = 'X', data=NULL)
dim(lupus)
```

    ## [1]   32738 1263676

## Filter out low quality genes and cells

``` r
#genes expressed in the cells

nF<-lupus@meta.data[["nFeature_originalexp"]]
min(nF) # all cells with counts in more than 100 genes
```

    ## [1] 114

``` r
# transcript  in the cells

nC <-lupus@meta.data[["nCount_originalexp"]]
min(nC) # all cells with more than 400 counts
```

    ## [1] 404

``` r
# mitochondrial DNA
lupus$mitoPercent <- PercentageFeatureSet(lupus, pattern = '^MT-')
mP <- lupus@meta.data[["mitoPercent"]]
fivenum(mP)# cellls with high percentage of mitochondrial genes
```

    ## [1]  0.000000  2.356021  3.167421  4.255319 80.192235

``` r
lupus <- subset(lupus, subset= mitoPercent<15) #remove cells with more than 15% of mitocondrial DNA
dim(lupus)
```

    ## [1]   32738 1260840

## Edit metadata

``` r
lupus@meta.data$Age <- as.numeric(as.character(lupus@meta.data$Age))

lupus$lab <-case_when(grepl('HC', lupus$ind_cov) ~'lab1', grepl('IGTB',lupus$ind_cov)~'lab2', grepl('control',lupus$ind_cov)~'lab2', .default = 'lab3') # Adds laboratory info

lupus$cell_id <- Cells(lupus) # Include the cell names in the metadata

lupus$sample_cell<-paste0(lupus$ind_cov_batch_cov,'_',lupus$cg_cov) # adds identifier for pseudobulk
lupus <- SetIdent(lupus, value = "sample_cell")
# write_rds(lupus,"/domino/datasets/local/RUV/lupusfiltcor.rds")

DimPlot(lupus, group.by = 'cg_cov')
```

![](data-wrangling_files/figure-gfm/Additional%20variables-1.png)<!-- -->

## Simulation subset

``` r
subsample <- dplyr::filter(lupus@meta.data, SLE_status =='Healthy', between(Age,24,28), Sex=='Female') #Filter healthy control females between 24 and 28 years old

total_cell_cg <- subsample%>% group_by(cg_cov) %>% summarise(ncell_cg=n()) # count cells per cell type
subsample <- left_join(subsample,total_cell_cg)%>% dplyr::filter(ncell_cg>500) # remove cells from cell types with less than 200 cells

total_cell_s <- subsample%>% group_by(sample_cell) %>% summarise(ncell_s=n()) # count cells per pseudobulk sample
subsample <- left_join(subsample,total_cell_s) %>% dplyr::filter(ncell_s>3) # remove cells from pseudobluk samples with less than 4 cells

subsample_bool<- lupus$cell_id %in% subsample$cell_id # create identifier
lupus$orig.index <- subsample_bool
Idents(object = lupus) <- "orig.index"
ds2 <- subset(x = lupus, idents = TRUE) # filter the selected cells 
ds2@meta.data <- droplevels(ds2@meta.data) # remove empty levels

dim(ds2)
```

    ## [1]  32738 131862

``` r
#write_rds(ds2,"/domino/datasets/local/RUV/sclupusds2cor.rds")

celltypes <- levels(ds2$cg_cov)
```

### Pseudobulk object

``` r
ds2 <- SetIdent(ds2, value = "sample_cell")
ds2 <- FindVariableFeatures(ds2) #### adds metafeatures on ds2@assays[["originalexp"]]@meta.features
ds2 <- ScaleData(ds2) #### adds scaledata on ds2@assays[["originalexp"]]@data

PBC2 <- AggregateExpression(ds2, group_by = "sample_cell", assays = 'originalexp',slot='counts', return.seurat = F)

metacovs <- colnames(ds2@meta.data)[c(4:6,8,10:16, 18,20)]

PBC2 <- DGEList(counts = PBC2$originalexp,samples = unique(ds2@meta.data[,metacovs]))
colnames(PBC2) <- unique(ds2$sample_cell)
PBC2 <- PBC2[rowSums(PBC2$counts >= 5) >= 5,]
dim(PBC2)
```

    ## [1] 13877   296

``` r
#write_rds(PBC2,"/domino/datasets/local/RUV/pblupusds2.rds")
```

``` r
PCAplot_ct <- function (PBC,ct){
  
  y <- PBC[rowSums(PBC$counts >= 5) >= 5,PBC$samples$cg_cov==ct]
  pmin <- find_p(y)
  y <- calcNormFactors(y, method="upperquartile",p=pmin)
  logy <- edgeR::cpm(y,log = T)
  
  pca <- calculatePCA(logy,ncomponents=2)
  
  ppca <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(y$samples, by = "sample_cell") %>% mutate(control= ind_cov=='ICC_control') %>%
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort)) + geom_point(aes(shape=pop_cov),size=3) + 
    geom_text(aes(label = ifelse(control, 'control', "")),color='black',vjust = 0, nudge_y = 0.5) +
      labs(color = "Processing cohort", shape = "Ethnicity") + 
    theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1')
    
  return(ppca)
  
}


pcaplots <- lapply(celltypes, function (x) PCAplot_ct(PBC2,x))
names(pcaplots) <- celltypes
```

### PCA plots

``` r
for (i in celltypes){
  
  cat(knitr::knit_expand(text=paste0('\n#### ',i, ' {-}\n')))
  
  print(pcaplots[[i]])

  cat(knitr::knit_expand(text=paste0("\n\n")))

}
```

#### B

![](data-wrangling_files/figure-gfm/PCA%20ds2%20-1.png)<!-- -->

#### NK

![](data-wrangling_files/figure-gfm/PCA%20ds2%20-2.png)<!-- -->

#### T4

![](data-wrangling_files/figure-gfm/PCA%20ds2%20-3.png)<!-- -->

#### T8

![](data-wrangling_files/figure-gfm/PCA%20ds2%20-4.png)<!-- -->

#### cDC

![](data-wrangling_files/figure-gfm/PCA%20ds2%20-5.png)<!-- -->

#### cM

![](data-wrangling_files/figure-gfm/PCA%20ds2%20-6.png)<!-- -->

#### ncM

![](data-wrangling_files/figure-gfm/PCA%20ds2%20-7.png)<!-- -->

#### pDC

![](data-wrangling_files/figure-gfm/PCA%20ds2%20-8.png)<!-- -->

## Study case subset Asian

``` r
#write_rds(lupus.filt,"/domino/datasets/local/RUV/lupusfiltcor.rds")
celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC') # select same cell types as in simulations

subsample <- dplyr::filter(lupus@meta.data, Status %in% c('Healthy','Managed'),Sex=="Female", between(Age,30,34),pop_cov == 'Asian', Sex=='Female', cg_cov %in% celltypes) #Filter asian females between 30 and 34 years old

total_cell_cg <- subsample%>% group_by(cg_cov) %>% summarise(ncell_cg=n()) # count cells per cell type
subsample <- left_join(subsample,total_cell_cg)%>% dplyr::filter(ncell_cg>200) # remove cells from cell types with less than 200 cells

total_cell_s <- subsample%>% group_by(sample_cell) %>% summarise(ncell_s=n()) # count cells per pseudobulk sample
subsample <- left_join(subsample,total_cell_s) %>% dplyr::filter(ncell_s>3) # remove cells from pseudobluk samples with less than 4 cells

subsample_bool<- lupus$cell_id %in% subsample$cell_id # create identifier
lupus$orig.index <- subsample_bool
Idents(object = lupus) <- "orig.index"
ds3 <- subset(x = lupus, idents = TRUE) # filter the selected cells 
ds3@meta.data <- droplevels(ds3@meta.data) # remove empty levels

dim(ds3)
```

    ## [1] 32738 80190

``` r
#write_rds(ds3,"/domino/datasets/local/RUV/sclupusds3cor.rds")
```

### Pseudobulk object

``` r
ds3 <- SetIdent(ds3, value = "sample_cell")
ds3 <- FindVariableFeatures(ds3) #### adds metafeatures on ds3@assays[["originalexp"]]@meta.features
ds3 <- ScaleData(ds3) #### adds scaledata on ds3@assays[["originalexp"]]@data

PBC3 <- AggregateExpression(ds3, group_by = "sample_cell", assays = 'originalexp',slot='counts', return.seurat = F)

metacovs <- colnames(ds3@meta.data)[c(4:6,8,10:16, 18,20)]

PBC3 <- DGEList(counts = PBC3$originalexp,samples = unique(ds3@meta.data[,metacovs]))
colnames(PBC3) <- unique(ds3$sample_cell)
PBC3 <- PBC3[rowSums(PBC3$counts >= 5) >= 5,]
dim(PBC3)
```

    ## [1] 12538   160

``` r
#write_rds(PBC3,"/domino/datasets/local/RUV/pblupusds3.rds")
```

``` r
PCAplot3_ct <- function (PBC,ct){
  
  y <- PBC[rowSums(PBC$counts >= 5) >= 5,PBC$samples$cg_cov==ct]
  pmin <- find_p(y)
  y <- calcNormFactors(y, method="upperquartile",p=pmin)
  logy <- edgeR::cpm(y,log = T)
  
  pca <- calculatePCA(logy,ncomponents=2)
  
  ppca <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(y$samples, by = "sample_cell") %>% 
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort)) + geom_point(aes(shape= SLE_status),size=3) + 
      labs(color = "Processing cohort", shape = "Status") + 
    theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1')
    
  return(ppca)
  
}


pcaplots3 <- lapply(celltypes, function (x) PCAplot3_ct(PBC3,x))
names(pcaplots3) <- celltypes
```

### PCA plots

``` r
for (i in celltypes){
  
  cat(knitr::knit_expand(text=paste0('\n#### ',i, ' {-}\n')))
  
  print(pcaplots3[[i]])

  cat(knitr::knit_expand(text=paste0("\n\n")))

}
```

#### B

![](data-wrangling_files/figure-gfm/PCA%20ds3%20-1.png)<!-- -->

#### NK

![](data-wrangling_files/figure-gfm/PCA%20ds3%20-2.png)<!-- -->

#### T4

![](data-wrangling_files/figure-gfm/PCA%20ds3%20-3.png)<!-- -->

#### T8

![](data-wrangling_files/figure-gfm/PCA%20ds3%20-4.png)<!-- -->

#### cDC

![](data-wrangling_files/figure-gfm/PCA%20ds3%20-5.png)<!-- -->

#### cM

![](data-wrangling_files/figure-gfm/PCA%20ds3%20-6.png)<!-- -->

#### ncM

![](data-wrangling_files/figure-gfm/PCA%20ds3%20-7.png)<!-- -->

#### pDC

![](data-wrangling_files/figure-gfm/PCA%20ds3%20-8.png)<!-- -->

## Study case subset mix

``` r
#lupus <- readRDS("/domino/datasets/local/RUV/lupusfiltcor.rds")

celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC') # select same cell types as in simulations

subsample <- dplyr::filter(lupus@meta.data, Status %in% c('Healthy','Managed'), between(Age,29,34), Sex=='Female', cg_cov %in% celltypes) #Filter females between 29 and 34 years old

total_cell_cg <- subsample%>% group_by(cg_cov) %>% summarise(ncell_cg=n()) # count cells per cell type
subsample <- left_join(subsample,total_cell_cg)%>% dplyr::filter(ncell_cg>200) # remove cells from cell types with less than 200 cells

total_cell_s <- subsample%>% group_by(sample_cell) %>% summarise(ncell_s=n()) # count cells per pseudobulk sample
subsample <- left_join(subsample,total_cell_s) %>% dplyr::filter(ncell_s>3) # remove cells from pseudobluk samples with less than 4 cells

subsample_bool<- lupus$cell_id %in% subsample$cell_id # create identifier
lupus$orig.index <- subsample_bool
Idents(object = lupus) <- "orig.index"
ds4 <- subset(x = lupus, idents = TRUE) # filter the selected cells 
ds4@meta.data <- droplevels(ds4@meta.data) # remove empty levels

dim(ds4)
```

    ## [1]  32738 284781

``` r
# write_rds(ds4,"/domino/datasets/local/RUV/sclupusds4.rds")
```

### Pseudobulk object

``` r
ds4 <- SetIdent(ds4, value = "sample_cell")
ds4 <- FindVariableFeatures(ds4) #### adds metafeatures on ds4@assays[["originalexp"]]@meta.features
ds4 <- ScaleData(ds4) #### adds scaledata on ds4@assays[["originalexp"]]@data

PBC4 <- AggregateExpression(ds4, group_by = "sample_cell", assays = 'originalexp',slot='counts', return.seurat = F)

metacovs <- colnames(ds4@meta.data)[c(4:6,8,10:16, 18,20)]

PBC4 <- DGEList(counts = PBC4$originalexp,samples = unique(ds4@meta.data[,metacovs]))
colnames(PBC4) <- unique(ds4$sample_cell)
PBC4 <- PBC4[rowSums(PBC4$counts >= 5) >= 5,]
dim(PBC4)
```

    ## [1] 14009   680

``` r
# write_rds(PBC4,"/domino/datasets/local/RUV/pblupusds4.rds")
```

``` r
PCAplot4_ct <- function (PBC,ct){
  
  y <- PBC[rowSums(PBC$counts >= 5) >= 5,PBC$samples$cg_cov==ct]
  pmin <- find_p(y)
  y <- calcNormFactors(y, method="upperquartile",p=pmin)
  logy <- edgeR::cpm(y,log = T)
  
  pca <- calculatePCA(logy,ncomponents=2)
  
  ppca <- data.frame(pca) %>%
  tibble::rownames_to_column("sample_cell") %>% # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(y$samples, by = "sample_cell") %>% 
  ggplot(aes(x = PC1,y = PC2,color= Processing_Cohort)) + geom_point(aes(shape= SLE_status),size=3) + 
      labs(color = "Processing cohort", shape = "Status") + 
    theme_minimal()+ theme(legend.position = "bottom") + scale_color_brewer(palette='Set1')
    
  return(ppca)
  
}


pcaplots4 <- lapply(celltypes, function (x) PCAplot4_ct(PBC4,x))
names(pcaplots4) <- celltypes
```

### PCA plots

``` r
for (i in celltypes){
  
  cat(knitr::knit_expand(text=paste0('\n#### ',i, ' {-}\n')))
  
  print(pcaplots4[[i]])

  cat(knitr::knit_expand(text=paste0("\n\n")))

}
```

#### B

![](data-wrangling_files/figure-gfm/PCA%20ds4%20-1.png)<!-- -->

#### NK

![](data-wrangling_files/figure-gfm/PCA%20ds4%20-2.png)<!-- -->

#### T4

![](data-wrangling_files/figure-gfm/PCA%20ds4%20-3.png)<!-- -->

#### T8

![](data-wrangling_files/figure-gfm/PCA%20ds4%20-4.png)<!-- -->

#### cDC

![](data-wrangling_files/figure-gfm/PCA%20ds4%20-5.png)<!-- -->

#### cM

![](data-wrangling_files/figure-gfm/PCA%20ds4%20-6.png)<!-- -->

#### ncM

![](data-wrangling_files/figure-gfm/PCA%20ds4%20-7.png)<!-- -->

#### pDC

![](data-wrangling_files/figure-gfm/PCA%20ds4%20-8.png)<!-- -->
