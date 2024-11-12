RUV Normalization
================

# Problem description

The single cell RNA Sequencing dataset used in the
*[paper](https://www.science.org/doi/10.1126/science.abf1970)* from
*Perez et. al.*, available at GEO accession number GSE174188 has a
strong technical effect linked to the cohorts where the samples were
processed. This processing cohort effect introduces unwanted variation
that can obscure the biological signals of interest.To address this
issue, we select a subset of healthy controls and evaluate the
effectiveness of RUV methods (RUV2, RUVIII, RUV4) to normalize the data
by removing the unwanted factors *W*, source of the unwanted variation.
Additionally, we included a fourth approach for comparison: the direct
Removal of the processing cohort effect, estimated using a linear model
of the log-transformed data.

# Benchmarking tools

- PCA (Principal Component Analysis) plots: To visually inspect
  clustering patterns and check if samples separate by cohort or by
  biological conditions after normalization. The PCA is performed over
  the logCPM values of the 500 most variable genes.

- Silhouette coefficients: To measure how well-defined the processing
  cohort clusters are before and after normalization. In our case, lower
  or negative silhouette coefficients are desirable, as they indicate
  reduced clustering driven by a known technical factor (the processing
  cohort). The silhouette coefficients are computed using the first 10
  PCs and the Manhattan distance.

- RLE (Relative Log Expression) plots: To observe differences in the
  gene expression of samples from different processing cohorts.

- p-value histograms (under no expected differential expression): To
  check if the normalization method induces a higher proportion of false
  positives. The p-values are computed using the Limma-Voom method to
  test for differential expression between 2 groups (mock treatment A vs
  mock treatment B) where no true differences are expected.

We demonstrate that the inclusion of pseudobulk pseudosamples (PBPS) in
the RUVIII method outperforms other methods (including the direct
removal of the Processing cohort effect) to reduce the unwanted
variation in the normalized datasets.

# Analysis

## Files

The following files were generated in the Data wrangling vignette and
contain the pseudobulk matrix and the pseudobulk pseudosamples.

``` r
PBC <- readRDS("~/PBPS/Data/pblupusds2.rds")
PBPSC  <- readRDS("~/PBPS/Data/pbps10rep_ds2.rds")
```

## Negative control genes and high variable genes

The following analysis uses the *[negative control
genes](https://www.worldscientific.com/doi/abs/10.1142/S0219720020400041)*
proposed by *Deeke and Gagnon-Bartsch*. We remove from the negative
control genes set the genes with a standardized variance higher than
1.8.

``` r
hkGagnon <- read.csv(paste0(path,"/Genes_Gagnon.txt"), sep="")
hkGagnon <- unlist(hkGagnon)
filter.vst <- FindVariableFeatures(PBC$counts, selection.method = 'vst')
top_high_regvarg<-rownames(arrange(filter.vst,-variance))[1:500]
high_varg <- which(filter.vst[hkGagnon,]$variance.standardized>1.8)

gene.D <- data.frame(gene=rownames(PBC))
rownames(gene.D) <- gene.D$gene
hk.ind <- rownames(PBC) %in% hkGagnon[-high_varg]
```

## Assigning the mock treatment

We randomly assign a mock treatment with 2 levels to later test for
differential expression.

``` r
subs <- unique(PBC$samples$ind_cov)
set.seed(1)
treatment <- sample(c("A","B"),length(subs), replace=T)

PBC$samples <- PBC$samples %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')
```

## Subsample

We have N = 37 assays from J = 29 subjects, the samples from processing
cohort 4 come from a different laboratory, with the exception of the
sample labelled as *control*, which was collected in the same laboratory
than the samples from the first 3 processing cohorts. The mock
treatments were assigned as follows:

![](results-norm_files/figure-gfm/dataset%20-1.png)<!-- -->

It is a standard in pseudobulk studies to analyze the samples from
different cell types separately. Therefore we will review the results
for 8 cell types present in the data.

## Unwanted variation

Before using any normalization method, we observe a strong effect
associated with the processing cohort in all cell types. The differences
between the RLE across samples is not high, but the silhouette
coefficients are, as well as the clustering effects in the first 2
principal components. Most of the p-values histograms are roughly
homogeneous, indicating no confounding effect between the processing
cohorts and the mock treatment.

### B

![](results-norm_files/figure-gfm/rawplots%20-1.png)<!-- -->

### NK

![](results-norm_files/figure-gfm/rawplots%20-2.png)<!-- -->

### T4

![](results-norm_files/figure-gfm/rawplots%20-3.png)<!-- -->

### T8

![](results-norm_files/figure-gfm/rawplots%20-4.png)<!-- -->

### cDC

![](results-norm_files/figure-gfm/rawplots%20-5.png)<!-- -->

### cM

![](results-norm_files/figure-gfm/rawplots%20-6.png)<!-- -->

### ncM

![](results-norm_files/figure-gfm/rawplots%20-7.png)<!-- -->

### pDC

![](results-norm_files/figure-gfm/rawplots%20-8.png)<!-- -->

## Removing the Processing Cohort effect

we attempted to estimate the processing cohort effect directly and
subtract it from the data in order to reduce technical variation.
However, this approach alone can actually introduce technical noise in
the genes where there was no strong processing cohort effect. As a
result we have worse PCA and RLE plots for most cell types.

### B

![](results-norm_files/figure-gfm/PCplots%20-1.png)<!-- -->

### NK

![](results-norm_files/figure-gfm/PCplots%20-2.png)<!-- -->

### T4

![](results-norm_files/figure-gfm/PCplots%20-3.png)<!-- -->

### T8

![](results-norm_files/figure-gfm/PCplots%20-4.png)<!-- -->

### cDC

![](results-norm_files/figure-gfm/PCplots%20-5.png)<!-- -->

### cM

![](results-norm_files/figure-gfm/PCplots%20-6.png)<!-- -->

### ncM

![](results-norm_files/figure-gfm/PCplots%20-7.png)<!-- -->

### pDC

![](results-norm_files/figure-gfm/PCplots%20-8.png)<!-- -->

## RUV methods

All RUV methods estimate 5 unwanted factors to capture the unwanted
variation in the data, associated for instance with the processing
cohort.

### RUV2

The RUV2 method uses solely the negative control genes to estimate the
unwanted factors via EFA, for a more detail explanation check the
*[technical report](https://statistics.berkeley.edu/tech-reports/820)*
from *Gagnon-Bartsch et. al.*.

Here we observe a lower average silhouette coefficient (but not bellow
0), and similar distributions in the RLE plots, indicating the removal
of the technical noise associated with the processing cohorts. One issue
observed with RUV2 was its inability to properly normalize the control
sample, which resulted in it appearing as an outlier in the PCA plots of
some cell types. This suggests that the method still has room for
improvement.

``` r
for (i in 1:length(celltypes)){
  
  cat(knitr::knit_expand(text=paste0('\n#### ',celltypes[i], ' {-}\n')))
  
  plots <- list(SilplotsT22[[i]], PCAplotsT22[[i]], RLEplotsT22[[i]], HistsT22[[i]])
 
  print(ggpubr::ggarrange(plotlist=plots,nrow=2,ncol=2, common.legend = T, legend="bottom"))

  cat(knitr::knit_expand(text=paste0("\n\n")))

}
```

#### B

![](results-norm_files/figure-gfm/ruv2T2plots%20-1.png)<!-- -->

#### NK

![](results-norm_files/figure-gfm/ruv2T2plots%20-2.png)<!-- -->

#### T4

![](results-norm_files/figure-gfm/ruv2T2plots%20-3.png)<!-- -->

#### T8

![](results-norm_files/figure-gfm/ruv2T2plots%20-4.png)<!-- -->

#### cDC

![](results-norm_files/figure-gfm/ruv2T2plots%20-5.png)<!-- -->

#### cM

![](results-norm_files/figure-gfm/ruv2T2plots%20-6.png)<!-- -->

#### ncM

![](results-norm_files/figure-gfm/ruv2T2plots%20-7.png)<!-- -->

#### pDC

![](results-norm_files/figure-gfm/ruv2T2plots%20-8.png)<!-- -->

### RUVIII

The RUVIII method uses both negative control genes and negative control
samples to estimate the unwanted factors via EFA, for a more detail
explanation check the *[original
paper](https://academic.oup.com/nar/article/47/12/6073/5494770?login=false)*
from *Molania et. al.*.

The dataset has technical replicates mostly between the processing
cohorts 2 and 3, therefore it was expected to observe the removal of
technical noise mainly between those two sets. Therefore the inclusion
of PBPS is relevant to capture the technical noise from all processing
cohorts.

#### B

![](results-norm_files/figure-gfm/ruv3T2plots%20-1.png)<!-- -->

#### NK

![](results-norm_files/figure-gfm/ruv3T2plots%20-2.png)<!-- -->

#### T4

![](results-norm_files/figure-gfm/ruv3T2plots%20-3.png)<!-- -->

#### T8

![](results-norm_files/figure-gfm/ruv3T2plots%20-4.png)<!-- -->

#### cDC

![](results-norm_files/figure-gfm/ruv3T2plots%20-5.png)<!-- -->

#### cM

![](results-norm_files/figure-gfm/ruv3T2plots%20-6.png)<!-- -->

#### ncM

![](results-norm_files/figure-gfm/ruv3T2plots%20-7.png)<!-- -->

#### pDC

![](results-norm_files/figure-gfm/ruv3T2plots%20-8.png)<!-- -->

### RUV4

The RUV4 method uses the negative control genes and the information from
the factor of interest (Mock treatment) to estimate the unwanted factors
via EFA, for a more detail explanation check the *[technical
report](https://statistics.berkeley.edu/tech-reports/820)* from
*Gagnon-Bartsch et. al.*.

The RUV4 method increases the false discovery rate in the samples from
T4 and NK cells, this is observed in the skewed p-values histrograms.
Therefore its use to normalize datasets is not recommended.

#### Code

#### B

![](results-norm_files/figure-gfm/ruv4T2plots%20-1.png)<!-- -->

#### NK

![](results-norm_files/figure-gfm/ruv4T2plots%20-2.png)<!-- -->

#### T4

![](results-norm_files/figure-gfm/ruv4T2plots%20-3.png)<!-- -->

#### T8

![](results-norm_files/figure-gfm/ruv4T2plots%20-4.png)<!-- -->

#### cDC

![](results-norm_files/figure-gfm/ruv4T2plots%20-5.png)<!-- -->

#### cM

![](results-norm_files/figure-gfm/ruv4T2plots%20-6.png)<!-- -->

#### ncM

![](results-norm_files/figure-gfm/ruv4T2plots%20-7.png)<!-- -->

#### pDC

![](results-norm_files/figure-gfm/ruv4T2plots%20-8.png)<!-- -->

### RUVIII with PBPs

To improve the removal of unwanted technical variation, we extended our
analysis by incorporating pseudobulk pseudosamples (PBPS) as described
in the PBPS_ISCB vignette. These PBPS were generated and included in the
dataset specifically for estimating unwanted factors using the RUVIII
method. However, it’s important to note that the pseudosamples were
exclusively utilized for calculating the unwanted factors and are not
meant to be used in any downstream analyses. They were included in the
first PCA plot only for informative purposes. The use of PBPS improves
the previous results. Using RUVIII with PBPS led to the lowest average
silhouette coefficients, indicating minimal clustering due to batch
effects (i.e., processing cohorts), it also resulted in homogeneous
distributions in the RLE plots, uniform p-values histogram and a better
batch correction in the control sample.

#### B

![](results-norm_files/figure-gfm/ruvpbps%20plots%20-1.png)<!-- -->

#### NK

![](results-norm_files/figure-gfm/ruvpbps%20plots%20-2.png)<!-- -->

#### T4

![](results-norm_files/figure-gfm/ruvpbps%20plots%20-3.png)<!-- -->

#### T8

![](results-norm_files/figure-gfm/ruvpbps%20plots%20-4.png)<!-- -->

#### cDC

![](results-norm_files/figure-gfm/ruvpbps%20plots%20-5.png)<!-- -->

#### cM

![](results-norm_files/figure-gfm/ruvpbps%20plots%20-6.png)<!-- -->

#### ncM

![](results-norm_files/figure-gfm/ruvpbps%20plots%20-7.png)<!-- -->

#### pDC

![](results-norm_files/figure-gfm/ruvpbps%20plots%20-8.png)<!-- -->
