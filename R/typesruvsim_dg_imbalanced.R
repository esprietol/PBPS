rm(list=ls())
setwd("/mnt")

BiocManager::install(c('EDASeq','iCOBRA','RUVSeq'),update = F)
install.packages(c('rlist','ruv','factoextra'))
system("R CMD INSTALL /domino/datasets/local/RUV/swapper-master/")

library(tidyverse)
library(BiocManager)
library(DESeq2)
library(EDASeq)
library(scuttle)
library(scater)
library(uwot)
library(edgeR)
library(ruv)
library(Seurat)
library(swapper)
library(rlist)
library(iCOBRA)
library(ggpubr)
library(cluster)
library(factoextra)


# sc data -----------------------------------------------------------------

source("/mnt/auxf.R")
ds2 <- readRDS("/domino/datasets/local/RUV/sclupusds2.rds")
ds2 <- ds2[,ds2$Sex=='Female']
Idents(object = ds2) <- "cg_cov"

celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')
ds2s <- subset(x = ds2, idents = celltypes)
rm(list = 'ds2')

pdata <- droplevels(ds2s@meta.data)
pdata$Age <- as.numeric(as.character(pdata$Age))

ds2s <- SetIdent(ds2s, value = "sample_cell") # change ident to group later
ds2s <- FindVariableFeatures(ds2s) #### adds metafeatures on ds2s@assays[["originalexp"]]@meta.features
ds2s <- ScaleData(ds2s) #### adds scaledata on ds2s@assays[["originalexp"]]@data

# pseudobulk --------------------------------------------------------------

PBC <- AggregateExpression(ds2s, group_by = "sample_cell", assays = 'originalexp',slot='counts', return.seurat = F)
PBC <- PBC$originalexp # same result than manual procedure
colnames(PBC) <- unique(pdata$sample_cell)
metacovs <- colnames(pdata)[c(4:6,8,10:16, 18:19)] # variables at type cell-sample level
pbc.metaD <- unique(pdata[,metacovs])
bc.metaD <- unique(pbc.metaD[,!colnames(pbc.metaD) %in% c("cg_cov" ,"sample_cell" )])
rownames(pbc.metaD) <- pbc.metaD$sample_cell 
pbc.metaD$Age <- as.numeric(as.character(pbc.metaD$Age))


# NCG ---------------------------------------------------------------------

path <- ('/domino/datasets/local/RUV')
hkGagnon <- read.csv(paste0(path,"/Genes_Gagnon.txt"), sep="")
hkGagnon <- unlist(hkGagnon)
filter.vst<-FindVariableFeatures(PBC,selection.method = 'vst')
high_varg <- which(filter.vst[hkGagnon,]$vst.variance.standardized>1.8)


# Fake treatment ----------------------------------------------------------

subs <- unique(pbc.metaD$ind_cov)
# set.seed(1)
# treatment <- sample(c("A","B"),length(subs), replace=T)

iter <- subs %in% pbc.metaD[pbc.metaD$Processing_Cohort=='4.0','ind_cov']
set.seed(1)
treatment <- NULL
treatment[iter] <- sample(c("A","B"),sum(iter), replace=T,prob=c(0.9,0.1))
treatment[!iter] <- sample(c("A","B"),sum(!iter), replace=T,prob=c(0.5,0.5))

pbc.metaD <- pbc.metaD %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')
pdata <- pdata %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')
bc.metaD <- bc.metaD %>% left_join(bind_cols(ind_cov = subs, fk.tr = treatment), by='ind_cov')




# Global filter and DGEList -----------------------------------------------------------

PBC <- PBC[rowSums(PBC >= 5) >= 5,]
gene.D <- data.frame(gene=rownames(PBC))
rownames(gene.D) <- gene.D$gene
hk.ind <- rownames(PBC) %in% hkGagnon[-high_varg]


# hyperparameters -----------------------------------------------------------

# ct = 'T4'
# i=1
seeds <- 1:100*1000
k=5
sw2its <- list()
swpits <- list()
pvalitsruv2 <- list()
pvalitsruv3 <- list()
pvalitsruv4 <- list()
truthits <- list()
padjits <- list()
for (i in 1:length(seeds)) {
  
  
  
  # In silico DEG -----------------------------------------------------------
  hk.ind <- rownames(PBC) %in% hkGagnon[-high_varg]
  
  samp_to_swap <- pbc.metaD$fk.tr == "A"
  
  sim <- PBC
  
  for(j in 1:length(celltypes)){
    set.seed(seeds[i]+j-1)
    sim.a <- simulateDE(PBC[!hk.ind,pbc.metaD$cg_cov==celltypes[j] ], which_cols = samp_to_swap[pbc.metaD$cg_cov==celltypes[j]], prop_DE = 0.1)
    gene.D[, paste0('trueDE_',celltypes[j])] <- FALSE
    gene.D[!hk.ind,  paste0('trueDE_',celltypes[j])] <- sim.a@elementMetadata@listData[["is_DE"]]
    sim.a <- rbind(assays(sim.a)$counts,PBC[hk.ind,pbc.metaD$cg_cov==celltypes[j] ])
    sim.a <- sim.a[gene.D$gene,]
    sim[,pbc.metaD$cg_cov==celltypes[j] ] <- sim.a    
  }
  
  
  sim.t3 <- split(pbc.metaD$sample_cell,pbc.metaD$ind_cov_batch_cov)
  sim.b <- sapply(sim.t3, function(x) rowSums(sim[,x])  )
  
  UQDE <- list()
  UQ_sw2 <- list()
  UQ_swp<- list()
  UQBDE <- list()
  UQB_sw2 <- list()
  UQB_swp<- list()
  
  t1ruv2DE <- list()
  t1ruv3DE <- list()
  t1ruv4DE <- list()
  t1ruv2_sw2 <- list()
  t1ruv2_swp<- list()
  t1ruv3_sw2 <- list()
  t1ruv3_swp<- list()
  t1ruv4_sw2 <- list()
  t1ruv4_swp<- list()
  
  t2ruv2DE <- list()
  t2ruv3DE <- list()
  t2ruv4DE <- list()
  t2ruv2_sw2 <- list()
  t2ruv2_swp<- list()
  t2ruv3_sw2 <- list()
  t2ruv3_swp<- list()
  t2ruv4_sw2 <- list()
  t2ruv4_swp<- list()
  
  t3ruv2DE <- list()
  t3ruv3DE <- list()
  t3ruv4DE <- list()
  t3ruv2_sw2 <- list()
  t3ruv2_swp<- list()
  t3ruv3_sw2 <- list()
  t3ruv3_swp<- list()
  t3ruv4_sw2 <- list()
  t3ruv4_swp<- list()
  
  truth <- list()
  
  
  # W Type 1 ------------------------------------------------------------------
  
  pbc.all <- sim
  pbc.all <-  pbc.all[rowSums(pbc.all >= 5) >= 5,]
  gene.D.all <- gene.D[ rownames(pbc.all),]
  
  hk.ind <- rownames(pbc.all) %in% hkGagnon[-high_varg]
  y.all <- DGEList(counts=pbc.all )
  pmin.all <- find_p(pbc.all )
  y.all <- calcNormFactors(y.all, method="upperquartile", p=pmin.all)
  nf.all <- y.all$samples$norm.factors
  logy.all <- edgeR::cpm(y.all,log = T)
  
  pbc.metaD$ind_covT1 <- factor(paste(pbc.metaD$ind_cov, pbc.metaD$cg_cov,sep="."))
  
  M <- replicate.matrix(pbc.metaD$ind_covT1)
  rownames(M) <- pbc.metaD$ind_covT1
  
  trT1 <- factor(paste(pbc.metaD$fk.tr, pbc.metaD$cg_cov,sep="."))
  
  ruv2T1 <- RUV2mod(Y = t(logy.all), X = trT1, ctl = hk.ind, k=k)
  ruv3T1 <- RUVIIIW(Y = t(logy.all), M = M, ctl = hk.ind,  k = k, return.info = T)
  ruv4T1 <- RUV4mod(Y = t(logy.all), X = trT1, ctl = hk.ind, k = k, Z = NULL)
  
  
  # W Type 3 ------------------------------------------------------------------
  
  bc <- sim.b
  bc <-  bc[rowSums(bc >= 5) >= 5,]
  gene.D.bc <- gene.D[ rownames(bc),]
  
  hk.ind <- rownames(bc) %in% hkGagnon[-high_varg]
  y.bc <- DGEList(counts=bc )
  pmin.bc <- find_p(bc )
  y.bc <- calcNormFactors(y.bc, method="upperquartile", p=pmin.bc)
  nf.bc <- y.bc$samples$norm.factors
  logy.bc <- edgeR::cpm(y.bc,log = T)
  
  M.bc <- replicate.matrix(bc.metaD$ind_cov)
  rownames(M.bc) <- bc.metaD$ind_cov
  
  ruv2T3 <- RUV2mod(Y = t(logy.bc), X = bc.metaD$fk.tr, ctl = hk.ind, k=k)
  ruv3T3 <- RUVIIIW(Y = t(logy.bc), M = M.bc, ctl = hk.ind,  k = k, return.info = T)
  ruv4T3 <- RUV4mod(Y = t(logy.bc), X = bc.metaD$fk.tr, ctl = hk.ind, k=k)
  
  
  # W Type 2 --------------------------------------------------------------
  
  for (ct in celltypes) {
    indsamp.ct <- bc.metaD$ind_cov_batch_cov %in% pbc.metaD$ind_cov_batch_cov[pbc.metaD$cg_cov==ct]
    pbc.ct <- sim[,pbc.metaD$cg_cov==ct]
    pbc.ct <-  pbc.ct[rowSums(pbc.ct >= 5) >= 5,]
    gene.D.ct <- gene.D[ rownames(pbc.ct),]
    gene.D.ct$ind <- paste0(ct,'_',gene.D.ct$gene,'_',i)
    truth[[ct]] <- gene.D.ct
    
    hk.ind <- rownames(pbc.ct) %in% hkGagnon[-high_varg]
    y.ct <- DGEList(counts=pbc.ct )
    pmin.ct <- find_p(pbc.ct )
    y.ct <- calcNormFactors(y.ct, method="upperquartile", p=pmin.ct)
    nf.ct <- y.ct$samples$norm.factors
    logy.ct <- edgeR::cpm(y.ct,log = T)
    
    indg.ct <- rownames(bc) %in% rownames(logy.ct)
    indg.ct2 <- rownames(logy.ct) %in% rownames(bc) 
    
    # UQ ----------------------------------------------------------------------
    
    pca <- calculatePCA(logy.ct,ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    UQ_sw2[[ct]] <- mean(sil2[,3])
    UQ_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct] ) 
    colnames(design) <- c('(Intercept)','tr')
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct$trueDE
    
    UQDE[[ct]] <- etable
    
    # UQB ----------------------------------------------------------------------
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('Proc', 2:4))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    alpha <- vfit$coefficients[,3:5] # same coefficients as efit and returned by toptable
    newY <- t(logy.ct) - design[,3:5]%*%t(alpha)
    
    pca <- calculatePCA(t(newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    UQB_sw2[[ct]] <- mean(sil2[,3])
    UQB_swp[[ct]] <- mean(silp[,3])
    
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct$trueDE
    
    UQBDE[[ct]] <- etable
    
    
    # RUV2 --------------------------------------------------------------
    
    pca <- calculatePCA(t(ruv2T1$newY[pbc.metaD$cg_cov==ct,]),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t1ruv2_sw2[[ct]] <- mean(sil2[,3])
    t1ruv2_swp[[ct]] <- mean(silp[,3])
    
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv2T1$W[pbc.metaD$cg_cov==ct,] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t1ruv2DE[[ct]] <- etable
    
    
    ruv2 <- RUV2mod(Y = t(logy.ct), X = pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct], ctl = hk.ind, k=k)
    
    pca <- calculatePCA(t(ruv2$newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv2_sw2[[ct]] <- mean(sil2[,3])
    t2ruv2_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv2$W ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t2ruv2DE[[ct]] <- etable
    
    
    newY <- t(logy.ct[indg.ct2,]) - ruv2T3$W[indsamp.ct,] %*% ruv2T3$fullalpha[,indg.ct]
    
    pca <- calculatePCA(t(newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t3ruv2_sw2[[ct]] <- mean(sil2[,3])
    t3ruv2_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv2T3$W[indsamp.ct,] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t3ruv2DE[[ct]] <- etable
    
    
    
    
    # RUV3 --------------------------------------------------------------------
    
    pca <- calculatePCA(t(ruv3T1$newY[pbc.metaD$cg_cov==ct,]),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t1ruv3_sw2[[ct]] <- mean(sil2[,3])
    t1ruv3_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv3T1$W[pbc.metaD$cg_cov==ct,] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t1ruv3DE[[ct]] <- etable
    
    
    Mct <- replicate.matrix(pbc.metaD$ind_cov[pbc.metaD$cg_cov==ct])
    rownames(Mct) <- pbc.metaD$ind_cov[pbc.metaD$cg_cov==ct]
    
    ruv3 <- RUVIIIW(Y = t(logy.ct), M = Mct, ctl=hk.ind,  k = k,
                    return.info = T)
    
    pca <- calculatePCA(t(ruv3$newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv3_sw2[[ct]] <- mean(sil2[,3])
    t2ruv3_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv3$W ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t2ruv3DE[[ct]] <- etable
    
    newY <- t(logy.ct[indg.ct2,]) - ruv3T3$W[indsamp.ct,] %*% ruv3T3$fullalpha[,indg.ct]
    
    pca <- calculatePCA(t(newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t3ruv3_sw2[[ct]] <- mean(sil2[,3])
    t3ruv3_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv3T3$W[indsamp.ct,] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t3ruv3DE[[ct]] <- etable
    
    # RUV4 --------------------------------------------------------------
    
    pca <- calculatePCA(t(ruv4T1$newY[pbc.metaD$cg_cov==ct,]),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t1ruv4_sw2[[ct]] <- mean(sil2[,3])
    t1ruv4_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv4T1$W[pbc.metaD$cg_cov==ct,] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t1ruv4DE[[ct]] <- etable
    
    
    ruv4 <- RUV4mod(Y = t(logy.ct), X = pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct], ctl = hk.ind, k=k)
    
    pca <- calculatePCA(t(ruv4$newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t2ruv4_sw2[[ct]] <- mean(sil2[,3])
    t2ruv4_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv4$W ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t2ruv4DE[[ct]] <- etable
    
    newY <- t(logy.ct[indg.ct2,]) - ruv4T3$W[indsamp.ct,] %*% ruv4T3$fullalpha[,indg.ct]
    
    pca <- calculatePCA(t(newY),ncomponents=10)
    pco <- as.integer(pbc.metaD$Processing_Cohort[pbc.metaD$cg_cov==ct])
    sil2 <- silhouette(x=pco,dist= dist(pca[,1:2], "euclidean"))
    silp <- silhouette(x=pco,dist= dist(pca, "manhattan"))
    
    t3ruv4_sw2[[ct]] <- mean(sil2[,3])
    t3ruv4_swp[[ct]] <- mean(silp[,3])
    
    design <- model.matrix(~ pbc.metaD$fk.tr[pbc.metaD$cg_cov==ct]+ ruv4T3$W[indsamp.ct,] ) 
    colnames(design) <- c('(Intercept)','tr',paste0('W', 1:k))
    v <- voom(y.ct, design, plot=F) 
    vfit <- lmFit(v, design)
    efit <- eBayes(vfit)
    
    etable <- topTable(efit,number=dim(y.ct)[1] , coef="tr", adjust="BH")
    etable <- etable[gene.D.ct$gene,]
    etable$DEG <- ifelse(etable$adj.P.Val<0.01,T,F)
    etable$trueDEG <- gene.D.ct[rownames(etable),2]
    
    t3ruv4DE[[ct]] <- etable
    
    
  }
  
  
  pval2 <- NULL 
  pval3 <- NULL
  pval4 <- NULL
  sw2 <- NULL
  swp <- NULL
  
  
  
  pval2[['ind']] <- mapply(function(x,ct) paste0(ct,'_',row.names(x),'_',i),UQDE,names(UQDE),SIMPLIFY = F)
  pval2[['UQ']] <- sapply(UQDE, function(x) x$P.Value, simplify = F)
  pval2[['UQ_Batch']] <- sapply(UQBDE, function(x) x$P.Value, simplify = F)
  pval2[['T1']] <- sapply(t1ruv2DE, function(x) x$P.Value, simplify = F)
  pval2[['T2']] <- sapply(t2ruv2DE, function(x) x$P.Value, simplify = F)
  pval2[['T3']] <- sapply(t3ruv2DE, function(x) x$P.Value, simplify = F)
  
  pval3[['ind']] <- mapply(function(x,ct) paste0(ct,'_',row.names(x),'_',i),UQDE,names(UQDE),SIMPLIFY = F)
  pval3[['UQ']] <- sapply(UQDE, function(x) x$P.Value, simplify = F)
  pval3[['UQ_Batch']] <- sapply(UQBDE, function(x) x$P.Value, simplify = F)
  pval3[['T1']] <- sapply(t1ruv3DE, function(x) x$P.Value, simplify = F)
  pval3[['T2']] <- sapply(t2ruv3DE, function(x) x$P.Value, simplify = F)
  pval3[['T3']] <- sapply(t3ruv3DE, function(x) x$P.Value, simplify = F)
  
  
  pval4[['ind']] <- mapply(function(x,ct) paste0(ct,'_',row.names(x),'_',i),UQDE,names(UQDE),SIMPLIFY = F)
  pval4[['UQ']] <- sapply(UQDE, function(x) x$P.Value, simplify = F)
  pval4[['UQ_Batch']] <- sapply(UQBDE, function(x) x$P.Value, simplify = F)
  pval4[['T1']] <- sapply(t1ruv4DE, function(x) x$P.Value, simplify = F)  
  pval4[['T2']] <- sapply(t2ruv4DE, function(x) x$P.Value, simplify = F)
  pval4[['T3']] <- sapply(t3ruv4DE, function(x) x$P.Value, simplify = F)
  
  sw2[['UQ']] <- UQ_sw2
  sw2[['UQ_Batch']] <- UQB_sw2
  sw2[['ruv2T1']] <- t1ruv2_sw2
  sw2[['ruv2T2']] <- t2ruv2_sw2
  sw2[['ruv2T3']] <- t3ruv2_sw2
  sw2[['ruv3T1']] <- t1ruv3_sw2
  sw2[['ruv3T2']] <- t2ruv3_sw2
  sw2[['ruv3T3']] <- t3ruv3_sw2
  sw2[['ruv4T1']] <- t1ruv4_sw2
  sw2[['ruv4T2']] <- t2ruv4_sw2
  sw2[['ruv4T3']] <- t3ruv4_sw2
  
  swp[['UQ']] <- UQ_swp
  swp[['UQ_Batch']] <- UQB_swp
  swp[['ruv2T1']] <- t1ruv2_swp
  swp[['ruv2T2']] <- t2ruv2_swp
  swp[['ruv2T3']] <- t3ruv2_swp
  swp[['ruv3T1']] <- t1ruv3_swp
  swp[['ruv3T2']] <- t2ruv3_swp
  swp[['ruv3T3']] <- t3ruv3_swp
  swp[['ruv4T1']] <- t1ruv4_swp
  swp[['ruv4T2']] <- t2ruv4_swp
  swp[['ruv4T3']] <- t3ruv4_swp
  
  sw2its[[i]] <- list_transpose(sw2) #%>% bind_cols()
  swpits[[i]] <- list_transpose(swp)
  
  pvalitsruv2[[i]] <- list_transpose(pval2) %>% lapply(bind_cols)
  pvalitsruv3[[i]] <- list_transpose(pval3) %>% lapply(bind_cols)
  pvalitsruv4[[i]] <- list_transpose(pval4) %>% lapply(bind_cols)
  
  truthits[[i]] <- truth
  
  print(i)
  
}


pvals <- list(pvalitsruv2,pvalitsruv3,pvalitsruv4,truthits)
saveRDS(pvals,paste0(path,'/pvals_ruvdg_tr_imbalance.rds'))

saveRDS(sw2its,paste0(path,'/asw2_ruvdg_tr_imbalance.rds'))
saveRDS(swpits,paste0(path,'/aswp_ruvdg_tr_imbalance.rds'))

# sw2its<- list_transpose(sw2its)
# 
# meanits2<- sapply(sw2its,function (x) bind_rows(x) %>% colMeans)
# 
# con.sil2 <- purrr::imap(lapply(sw2its,bind_rows), ~mutate(.x, cg_cov = .y)) %>% bind_rows()
# 
# con.sil2 <- pivot_longer(con.sil2,-cg_cov, names_to ='Method', values_to="ASW")
# 
# ggplot(filter(con.sil2, Method %in% c('UQ', 'UQ_Batch', 'ruv2T1','ruv2T2','ruv2T3')), aes(y= ASW, x = Method, color=Method)) + geom_boxplot() + facet_wrap('cg_cov')
# 
# 
# swpits<- list_transpose(swpits)
# 
# meanitsp<- sapply(swpits,function (x) bind_rows(x) %>% colMeans)
# 
# con.silp <- purrr::imap(lapply(swpits,bind_rows), ~mutate(.x, cg_cov = .y)) %>% bind_rows()
# 
# con.silp <- pivot_longer(con.silp,-cg_cov, names_to ='Method', values_to="ASW")
# 
# ggplot(filter(con.silp, Method %in% c('UQ', 'UQBatch', 'ruv2T1','ruv2T2','ruv2T3')), aes(y= ASW, x = Method, color=Method)) + geom_boxplot() + facet_wrap('cg_cov')
# 
# ggplot(filter(con.silp, Method %in% c('UQ', 'UQBatch', 'ruv3T1','ruv3T2','ruv3T3')), aes(y= ASW, x = Method, color=Method)) + geom_boxplot() + facet_wrap('cg_cov')
# 
# ggplot(filter(con.silp, Method %in% c('UQ', 'UQBatch', 'ruv4T1','ruv4T2','ruv4T3')), aes(y= ASW, x = Method, color=Method)) + geom_boxplot() + facet_wrap('cg_cov')
# 
# 
# 

# pvalits2 <- list_transpose(pvalitsruv2)
# pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )
# 
# truthits2 <- list_transpose(truthits)
# truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )
# 
# 
# COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
# COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
# cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
#                         binary_truth = paste0("trueDE_",celltypes))
# cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2", 
#                       facetted = TRUE)
# 
# plotsC <- list()
# 
# plotsC [['Fruv2']] <- lapply(cobratoplot,plot_fdrtprcurve,plottype='points',pointsize=0.5)
# plotsC [['Truv2']] <- lapply(cobratoplot,plot_tpr)
# 
# pvalits2 <- list_transpose(pvalitsruv3)
# pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )
# 
# COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
# COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
# cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
#                         binary_truth = paste0("trueDE_",celltypes))
# cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2", 
#                       facetted = TRUE)
# 
# 
# plotsC [['Fruv3']] <- lapply(cobratoplot,plot_fdrtprcurve,plottype='points',pointsize=0.5)
# plotsC [['Truv3']] <- lapply(cobratoplot,plot_tpr)
# 
# 
# pvalits2 <- list_transpose(pvalitsruv4)
# pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )
# 
# COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
# COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
# cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
#                         binary_truth = paste0("trueDE_",celltypes))
# cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2", 
#                       facetted = TRUE)
# 
# 
# plotsC [['Fruv4']] <- lapply(cobratoplot,plot_fdrtprcurve,plottype='points',pointsize=0.5)
# plotsC [['Truv4']] <- lapply(cobratoplot,plot_tpr)
# 
# 
# 
# saveRDS(plotsC, paste0(path,'/cobra_ruvdg_tr_imbalance.rds'))
# 


