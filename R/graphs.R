load('/domino/datasets/local/RUV//simsD1.RData')
#library(egg)
library(tidyverse)
library(rlist)
library(ggpubr)
library(iCOBRA)
source("/mnt/auxf.R")

path <- ('/domino/datasets/local/RUV')
plotscobra <- readRDS(paste0(path,'/cobra_limma_ruvps.rds')) 

plotscobra[[1]]+ 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.7) + ylim(0.82,0.95) +theme_minimal()
  theme(strip.background = element_blank(),strip.text =element_blank(),
         axis.title.x = element_text(size=rel(1)), axis.title.y = element_text(size=rel(1))
         )


  
plotscobra1 <- lapply(plotscobra[[1]], function(x) x + 
                        geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                        xlim(0,0.5) +
                        #ylim(0.75,0.92) +
                        theme_minimal()+
                        theme(strip.background = element_blank(),strip.text =element_blank()#,
                              # axis.title.x = element_text(size=rel(1)),
                              # axis.title.y = element_text(size=rel(1))
                        ))  
  

plotscobra2 <- lapply(plotscobra[[1]], function(x) x + 
                        geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                        xlim(0,0.15) + ylim(0.75,0.92) + theme_minimal()+
                        theme(strip.background = element_blank(),strip.text =element_blank()#,
                              # axis.title.x = element_text(size=rel(1)),
                              # axis.title.y = element_text(size=rel(1))
                              ))

plotscobra3 <- lapply(plotscobra[[3]], function(x) x + 
                        #geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                        xlim(0.78,0.92) + 
                        #ylim(0.75,0.92) +
                        theme_minimal()+
                        theme(strip.background = element_blank(),strip.text =element_blank()#,
                              # axis.title.x = element_text(size=rel(1)),
                              # axis.title.y = element_text(size=rel(1))
                        ))

path2 <- ('/domino/datasets/local/RUV/jul24')

 ggpubr::ggarrange(plotlist=plotscobra1,nrow=2,ncol=4,labels=names(plotscobra2),
                                   common.legend = T, legend="bottom")

ggsave(paste0(path2,'fdr_nonc.png'))

ggpubr::ggarrange(plotlist=plotscobra2,nrow=2,ncol=4,labels=names(plotscobra2),
                                   common.legend = T, legend="bottom")

ggsave(paste0(path2,'fdr_zoom.png'))

ggpubr::ggarrange(plotlist=plotscobra3,nrow=2,ncol=4,labels=names(plotscobra2),
                                   common.legend = T, legend="bottom")

ggsave(paste0(path2,'tpr_nonc_zoom.png'))


plotscobra3 <- lapply(plotscobra2, function(x) x + xlim(0,0.2) + ylim(0.82,0.91) )



plotscobraC2 <- lapply(plotscobraC, function(x) x  +
                      xlim(0,0.4) + ylim(0.6,0.95) + theme_minimal()+
                      theme(strip.background = element_blank(),strip.text =element_blank()))






ggpubr::ggarrange(plotlist=plotscobraC2,nrow=2,ncol=4,labels=names(plotscobra2),
                  common.legend = T, legend="bottom")

ggpubr::ggarrange(plotlist=plotscobra3,nrow=2,ncol=4,labels=names(plotscobra3),
                  common.legend = T, legend="bottom")

ggpubr::ggarrange(plotlist=lapply(plotscobra, function(x) x+ xlim(0,0.1)+ylim(0.5,0.9)),
                  nrow=2,ncol=4,labels=names(plotscobra),common.legend = T, legend="bottom")

ggpubr::ggarrange(plotlist=lapply(plotscobratpr, function(x) x+ xlim(0.75,1) +theme_minimal()+
                    theme(strip.background = element_blank(),strip.text =element_blank())),
                                    nrow=2,ncol=4, labels=names(plotscobratpr),
                    common.legend = T,legend="bottom")



# Pvals no PBPS -------------------------------------------------------------------
pvals_ruv_diffG <- readRDS("/domino/datasets/local/RUV/pvals_ruv_diffG.rds")
pvalitsruv2 <- pvals_ruv_diffG[[1]]
pvalitsruv3 <- pvals_ruv_diffG[[2]]
pvalitsruv4 <- pvals_ruv_diffG[[3]]
truthits <- pvals_ruv_diffG[[4]]


pvalits2 <- list_transpose(pvalitsruv2)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                               xlim(0,0.5) +
                               ylim(0.75,0.9) + 
                               theme_minimal()+
                               scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                               theme(strip.background = element_blank(),strip.text =element_blank())
                               )

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV2FDRDiffG.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv3)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.5) +
                     ylim(0.75,0.9) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV3FDRDiffG.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv4)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.5) +
                     ylim(0.75,0.9) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV4FDRDiffG.png'),height = 4, width = 12)

#pvals_ruvdg_tr_imbalance <- readRDS("/domino/datasets/local/RUV/pvals_ruvdg_tr_imbalance.rds")

pvalitsruv2 <- pvals_ruvdg_tr_imbalance[[1]]
pvalitsruv3 <- pvals_ruvdg_tr_imbalance[[2]]
pvalitsruv4 <- pvals_ruvdg_tr_imbalance[[3]]
truthits <- pvals_ruvdg_tr_imbalance[[4]]


pvalits2 <- list_transpose(pvalitsruv2)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.5) +
                     ylim(0.75,0.9) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV2FDRDimb.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv3)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.5) +
                     ylim(0.75,0.9) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV3FDRDimb.png'),height = 4, width = 12)

pvalits2 <- list_transpose(pvalitsruv4)
pvalits2 <- lapply(pvalits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthits)
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot, colorscheme = "Dark2",
                      facetted = TRUE)

fdrplots <- lapply(cobratoplot,function (x) plot_fdrtprcurve(x,plottype='points',pointsize=0.5) + geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                     xlim(0,0.5) +
                     ylim(0.75,0.9) + 
                     theme_minimal()+
                     scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
                     theme(strip.background = element_blank(),strip.text =element_blank())
)

ggpubr::ggarrange(plotlist=fdrplots,nrow=2,ncol=4,labels=names(fdrplots),
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/allcellsRUV4FDRDimb.png'),height = 4, width = 12)




# pvals types -------------------------------------------------------------

pvals_ruvdg <- readRDS("/domino/datasets/local/RUV/pvals_ruv_diffG.rds")

pvalitsruv2 <- pvals_ruvdg[[1]]
pvalitsruv3 <- pvals_ruvdg[[2]]
pvalitsruv4 <- pvals_ruvdg[[3]]
truthits <- pvals_ruvdg[[4]]

pvalits2 <- list_transpose(pvalitsruv2)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 
truthits2 <- list_transpose(truthits)$T4
truthits2 <- bind_rows(truthits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV2PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.3) +
  ylim(0.75,0.9) + 
  theme_minimal()+ ggtitle("RUV2") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

RUV2PLOT

pvalits2 <- list_transpose(pvalitsruv3)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV3PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.3) +
  ylim(0.75,0.9) +
  theme_minimal()+ ggtitle("RUVIII") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

pvalits2 <- list_transpose(pvalitsruv4)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 


COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV4PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.3) +
  ylim(0.75,0.9) + 
  theme_minimal()+ ggtitle("RUV4") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

RUV2PLOT <- RUV2PLOT + xlim(0,0.2)+ ylim(0.82,0.9)
RUV3PLOT <- RUV3PLOT + xlim(0,0.2)+ ylim(0.82,0.9)
RUV4PLOT <- RUV4PLOT + xlim(0,0.5) + ylim(0.82,0.9)

ggpubr::ggarrange(plotlist=list(RUV2PLOT,RUV3PLOT,RUV4PLOT),nrow=1,ncol=3,
                  common.legend = T, legend="bottom")
ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRDiffGT4.png'),height = 4, width = 12)


pvals_ruvdg_tr_imbalance <- readRDS("/domino/datasets/local/RUV/pvals_ruvdg_tr_imbalance.rds")

pvalitsruv2 <- pvals_ruvdg_tr_imbalance[[1]]
pvalitsruv3 <- pvals_ruvdg_tr_imbalance[[2]]
pvalitsruv4 <- pvals_ruvdg_tr_imbalance[[3]]
truthits <- pvals_ruvdg_tr_imbalance[[4]]

pvalits2 <- list_transpose(pvalitsruv2)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 
truthits2 <- list_transpose(truthits)$T4
truthits2 <- bind_rows(truthits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)

RUV2PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.5) +
  ylim(0.82,0.9) +
  theme_minimal()+ ggtitle("RUV2") + 
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())


pvalits2 <- list_transpose(pvalitsruv3)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 

COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV3PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.5) +
  ylim(0.82,0.9) + 
  theme_minimal()+ ggtitle("RUVIII") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

pvalits2 <- list_transpose(pvalitsruv4)$T4
pvalits2 <- bind_rows(pvalits2) %>% remove_rownames()%>%column_to_rownames('ind') 


COBRA <- COBRAData(pval=pvalits2,truth = truthits2)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUV4PLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.5) +
  ylim(0.82,0.9) + 
  theme_minimal()+ ggtitle("RUV4") +
  scale_color_brewer(palette='Set1',labels = c('T1', 'T2', 'T3', 'UQ', 'UQ Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())


ggpubr::ggarrange(plotlist=list(RUV2PLOT,RUV3PLOT,RUV4PLOT),nrow=1,ncol=3,
                  common.legend = T, legend="bottom")
ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRDimbT4.png'),height = 4, width = 12)


# new code ----------------------------------------------------------------

pvalspbps10_dgim <- readRDS("/domino/datasets/local/RUV/pvalspbps10_dgim.rds")
truthpbps10_dgim <- readRDS("/domino/datasets/local/RUV/truthpbps10_dgim.rds")

# celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')
celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')

pvalits2 <- list_transpose(pvalspbps10_dgim)
pvalits2 <- pvalits2[celltypes]
pvalits2 <- lapply(pvalits2 , function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthpbps10_dgim )
truthits2 <- truthits2[celltypes]
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot,
                      facetted = TRUE)

plotscobra <- lapply(cobratoplot,plot_fdrtprcurve,plottype='points',pointsize=0.5)


plotscobra1 <- mapply( function(x,y) x + 
                         geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                         xlim(0,0.4) +
                         ylim(0.81,0.93) +
                         theme_minimal()+ ggtitle(y) +
                         scale_color_brewer(palette='Set1',labels = c('RUVIII', 'RUVIII+Batch', 'RUVIII PBPS', 'UQ Batch')) +
                         theme(strip.background = element_blank(),strip.text =element_blank())
                         ,x=plotscobra,y=names(plotscobra),SIMPLIFY = F) 


ggpubr::ggarrange(plotlist=plotscobra1,nrow=2,ncol=4,
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRruviii.png'),height = 8, width = 12)

ggpubr::ggarrange(plotlist=list(plotscobra1$B,plotscobra1$T4,plotscobra1$ncM),nrow=1,ncol=3,
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRruv3sub.png'),height = 4, width = 12)

plotscobraC <- lapply(cobratoplot,plot_fdrtprcurve)


ggpubr::ggarrange(plotlist=plotscobraC2,nrow=2,ncol=4,labels=names(plotscobra2),
                  common.legend = T, legend="bottom")


 pvalspbps10_dgim <- readRDS("/domino/datasets/local/RUV/pvalspbps10_dgim5p.rds")
truthpbps10_dgim <- readRDS("/domino/datasets/local/RUV/truthpbps10_dgim5p.rds")

# celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')
celltypes <- c('B', 'NK', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC')

pvalits2 <- list_transpose(pvalspbps10_dgim)
pvalits2 <- pvalits2[celltypes]
pvalits2 <- lapply(pvalits2 , function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )

truthits2 <- list_transpose(truthpbps10_dgim )
truthits2 <- truthits2[celltypes]
truthits2 <- lapply(truthits2, function (x) bind_rows(x) %>% remove_rownames()%>%column_to_rownames('ind') )


COBRADatalist <- mapply(COBRAData,pval=pvalits2,truth = truthits2, SIMPLIFY = F)
COBRADatalist <- lapply(COBRADatalist,calculate_adjp)
cobraperflist <- mapply(calculate_performance, cobradata = COBRADatalist,
                        binary_truth = paste0("trueDE_",celltypes))
cobratoplot <- lapply(cobraperflist, prepare_data_for_plot,
                      facetted = TRUE)

plotscobra <- lapply(cobratoplot,plot_fdrtprcurve,plottype='points',pointsize=0.5)


plotscobra1 <- mapply( function(x,y) x + 
                         geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
                         xlim(0,0.2) +
                         ylim(0.75,0.9) +
                         theme_minimal()+ ggtitle(y) +
                         theme(strip.background = element_blank(),strip.text =element_blank()#,
                               # axis.title.x = element_text(size=rel(1)),
                               # axis.title.y = element_text(size=rel(1))
                         ),x=plotscobra,y=names(plotscobra),SIMPLIFY = F) 


ggpubr::ggarrange(plotlist=plotscobra1,nrow=2,ncol=4,
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRruv35p.png'),height = 8, width = 12)

ggpubr::ggarrange(plotlist=list(plotscobra1$B,plotscobra1$T4,plotscobra1$ncM),nrow=1,ncol=3,
                  common.legend = T, legend="bottom")

ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/DEA/FDRruv3sub.png'),height = 4, width = 12)

plotscobraC <- lapply(cobratoplot,plot_fdrtprcurve)


ggpubr::ggarrange(plotlist=plotscobraC2,nrow=2,ncol=4,labels=names(plotscobra2),
                  common.legend = T, legend="bottom")





# 
# 
# #awspDG <- readRDS("/domino/datasets/local/RUV/aswp_ruv_diffG.rds")
# #awspDG <- readRDS("/domino/datasets/local/RUV/aswp_ruv_diffG_fktr.rds")
# awspDG <- readRDS("/domino/datasets/local/RUV/aswp_ruv_diffG_et.rds")
# 
# con.silp <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/aswp_ruv_diffG_et.rds"))
# 
# for (ct in celltypes){
#   
#   aws_plot(con.silp, ct=ct) 
#   #ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/normalisation/',ct,'batch_silhouettep.png'),height = 5, width = 6)  
#   #ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/normalisation/',ct,'fktr_silhouettep.png'),height = 5, width = 6)
#   ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/normalisation/',ct,'et_silhouettep.png'),height = 5, width = 6)
#   
#    
# }




awsT4fktr <- asw_plot_prep( readRDS("/domino/datasets/local/RUV/aswp_ruv_diffG_fktr.rds"))%>%filter(Approach %in% c('T3','UQ','UQ_Batch')&RUV%in%c('no_ruv','ruv2'))|> aws_plot( ct='T4') 



awsT4PC <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/aswp_ruv_diffG.rds")) |> aws_plot( ct='T4') 
awsT4fktr <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/aswp_ruv_diffG_fktr.rds"))|> aws_plot( ct='T4') 
awsT4et <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/aswp_ruv_diffG_et.rds"))|> aws_plot( ct='T4') 



# ggpubr::ggarrange(plotlist=list(awsT4PC, awsT4fktr, awsT4et),nrow=1,ncol=3,labels=c('A', 'B', 'C'),
#                   common.legend = T, legend="bottom")
# 
# ggsave('/domino/datasets/local/RUV/paper/graphs/normalisation/allT4_silhouettep.png',height = 5, width = 7)




awsT4PC <- awsT4PC+facet_wrap("RUV", ncol=4,labeller = as_labeller(c('no_ruv'='No RUV', 'ruv2'='RUV2', 'ruv3'='RUVIII', 'ruv4'='RUV4')))
awsT4fktr <- awsT4fktr+facet_wrap("RUV", ncol=4,labeller = as_labeller(c('no_ruv'='No RUV', 'ruv2'='RUV2', 'ruv3'='RUVIII', 'ruv4'='RUV4')))

ggpubr::ggarrange(plotlist=list(awsT4PC, awsT4fktr),nrow=2,ncol=1,labels=c('A', 'B'),
                  common.legend = T, legend="bottom")

ggsave('/domino/datasets/local/RUV/paper/graphs/normalisation/allT4_silhouettep.png',height = 5, width = 7)



awsT4PC <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG.rds")) |> aws_plot( ct='T4') 
awsT4fktr <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG_fktr.rds"))|> aws_plot( ct='T4') 
awsT4et <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG_et.rds"))|> aws_plot( ct='T4') 


ggpubr::ggarrange(plotlist=list(awsT4PC, awsT4fktr, awsT4et),nrow=1,ncol=3,labels=c('A', 'B', 'C'),
                  common.legend = T, legend="bottom")

ggsave('/domino/datasets/local/RUV/paper/graphs/normalisation/allT4_silhouette2.png',height = 5, width = 7)



#aws2DG <- readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG.rds")
aws2DG <- readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG_fktr.rds")

sw2its<- list_transpose(aws2DG)

meanits2<- sapply(sw2its,function (x) bind_rows(x) %>% colMeans)

con.sil2 <- purrr::imap(lapply(sw2its,bind_rows), ~mutate(.x, cg_cov = .y)) %>% bind_rows()

con.sil2 <- pivot_longer(con.sil2,-cg_cov, names_to ='Method', values_to="ASW")

con.sil2 <- con.sil2 %>%mutate(Approach =sub(".*ruv.", "", Method), RUV= str_extract(Method, "ruv.")) %>%
  mutate(Approach = trimws(Approach)) 

con.sil2$RUV[is.na(con.sil2$RUV)] <- 'no_ruv'

for (ct in celltypes){
  
  ggplot(filter(con.sil2,cg_cov==ct), aes(y= ASW,  color=Approach, group=Approach)) + 
    geom_boxplot() + facet_wrap('RUV') + theme_minimal() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +
    ggtitle(paste0(ct," Processing Cohort's Average Silhouette Width"))
  
  ggsave(paste0('/domino/datasets/local/RUV/paper/graphs/normalisation/',ct,'batch_silhouette2.png'),height = 5, width = 6)
  
  
}



asw2fktr <- readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG_fktr.rds")





# pbps asw ----------------------------------------------------------------

aswpbps10_dgim <- readRDS("/domino/datasets/local/RUV/aswpbps10_dgim.rds")

test <- map_depth(aswpbps10_dgim,3, function(x) x[1,2])


swpits<- list_transpose(test)

con.silp <- purrr::imap(lapply(swpits,bind_rows), ~mutate(.x, Method = .y)) %>% bind_rows()

con.silp <- pivot_longer(con.silp,-Method, names_to ='cg_cov', values_to="ASW")

con.silp <- con.silp %>%mutate(RUV= "", Approach= Method) 

p1 <- aws_plot(con.silp, ct="T4")+ labs('no_ruv'='No RUV', 'ruv2'='RUV2', 'ruv3'='RUVIII', 'ruv4'='RUV4')

test <- map_depth(aswpbps10_dgim,3, function(x) x[2,2])


swpits<- list_transpose(test)

con.silp <- purrr::imap(lapply(swpits,bind_rows), ~mutate(.x, Method = .y)) %>% bind_rows()

con.silp <- pivot_longer(con.silp,-Method, names_to ='cg_cov', values_to="ASW")

con.silp <- con.silp %>%mutate(RUV= "", Approach= Method) 

p2 <- aws_plot(con.silp, ct="T4")+



ggpubr::ggarrange(plotlist=list(p1, p2),nrow=1,ncol=2,labels=c('A', 'B'),
                  common.legend = T, legend="bottom")

ggsave('/domino/datasets/local/RUV/paper/graphs/normalisation/allpbpsT4_silhouettep.png',height = 5, width = 7)


awsT4PC <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG.rds")) |> aws_plot( ct='T4') 
awsT4fktr <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG_fktr.rds"))|> aws_plot( ct='T4') 
awsT4et <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG_et.rds"))|> aws_plot( ct='T4') 


ggpubr::ggarrange(plotlist=list(awsT4PC, awsT4fktr, awsT4et),nrow=1,ncol=3,labels=c('A', 'B', 'C'),
                  common.legend = T, legend="bottom")

ggsave('/domino/datasets/local/RUV/paper/graphs/normalisation/allT4_silhouette2.png',height = 5, width = 7)





# pbps balanced design asw ----------------------------------------------------------------

aswpbps10_dgbal <- readRDS("/domino/datasets/local/RUV/aswpbps10_dgbal.rds")


swpits<- list_transpose(map_depth(aswpbps10_dgbal,3, function(x) x[1,2]))

con.silp <- purrr::imap(lapply(swpits,bind_rows), ~mutate(.x, Method = .y)) %>% bind_rows()

con.silp <- pivot_longer(con.silp,-Method, names_to ='cg_cov', values_to="ASW")

con.silp <- con.silp %>%mutate(RUV= "", Approach= Method) 

awsT4PC <- aws_plot(con.silp, ct="T4")

swpits<- list_transpose(map_depth(aswpbps10_dgbal,3, function(x) x[1,2]))

con.silp <- purrr::imap(lapply(swpits,bind_rows), ~mutate(.x, Method = .y)) %>% bind_rows()

con.silp <- pivot_longer(con.silp,-Method, names_to ='cg_cov', values_to="ASW")

con.silp <- con.silp %>%mutate(RUV= "", Approach= Method) 







unique(con.silp$Approach)



ggpubr::ggarrange(plotlist=list(awsT4PC, awsT4fktr, awsT4et),nrow=1,ncol=3,labels=c('A', 'B', 'C'),
                  common.legend = T, legend="bottom")

ggsave('/domino/datasets/local/RUV/paper/graphs/normalisation/allT4_silhouettep.png',height = 5, width = 7)


awsT4PC <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG.rds")) |> aws_plot( ct='T4') 
awsT4fktr <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG_fktr.rds"))|> aws_plot( ct='T4') 
awsT4et <- asw_plot_prep(readRDS("/domino/datasets/local/RUV/asw2_ruv_diffG_et.rds"))|> aws_plot( ct='T4') 


ggpubr::ggarrange(plotlist=list(awsT4PC, awsT4fktr, awsT4et),nrow=1,ncol=3,labels=c('A', 'B', 'C'),
                  common.legend = T, legend="bottom")

ggsave('/domino/datasets/local/RUV/paper/graphs/normalisation/allT4_silhouette2.png',height = 5, width = 7)






# poster FDR --------------------------------------------------------------
pvalspbps10_dgim <- readRDS("/domino/datasets/local/RUV/pvalspbps10_dgim.rds")
truthpbps10_dgim <- readRDS("/domino/datasets/local/RUV/truthpbps10_dgim.rds")

pvals_ruvdg_tr_imbalance <- readRDS("/domino/datasets/local/RUV/pvals_ruvdg_tr_imbalance.rds")
pvalitsruv2 <- pvals_ruvdg_tr_imbalance[[1]]
truthits <- pvals_ruvdg_tr_imbalance[[4]]

pvalits2 <- list_transpose(pvalitsruv2)$T4
pvalits2 <- bind_rows(pvalits2) 
pvalitspb <- list_transpose(pvalspbps10_dgim)$T4
pvalitspb <- bind_rows(pvalitspb) 
colnames(pvalitspb)<- c("ind", "UQ_Batchpb", "RUV3", "RUV3B", "RUV3PS")
pvalits <- inner_join(pvalits2,pvalitspb, by='ind')%>% remove_rownames()%>%column_to_rownames('ind')



truthits2 <- list_transpose(truthits)$T4
truthits2 <- bind_rows(truthits2) #%>% remove_rownames()%>%column_to_rownames('ind') 
truthitspb <- list_transpose(truthpbps10_dgim)$T4
truthitspb <- bind_rows(truthitspb)
colnames(truthitspb)<- c(paste0(colnames(truthitspb)[1:9],'2'),'ind') 
truthit <- inner_join(truthits2,truthitspb, by = 'ind')%>% remove_rownames()%>%column_to_rownames('ind')
check <- sum(truthit$trueDE_T4!=truthit$trueDE_T42)

pvalits <- pvalits[,c("UQ","UQ_Batchpb", "RUV3",  "RUV3B", "RUV3PS")] #both batch are equal
COBRA <- COBRAData(pval=pvalits,truth = truthit)
COBRA <- calculate_adjp(COBRA)
COBRA <- calculate_performance(cobradata = COBRA, binary_truth = "trueDE_T4")
cobratoplot <- prepare_data_for_plot(COBRA, facetted = TRUE)
RUVPLOT <- plot_fdrtprcurve(cobratoplot,plottype='points',pointsize=0.5) + 
  geom_point(size=5,shape=17,aes(colour= method,fill=method)) +
  xlim(0,0.6) +
  ylim(0.8,0.95) + 
  theme_minimal()+ ggtitle("In silico differential expression") + #theme(legend.position = "bottom")+
  scale_color_brewer(palette='Set1',labels = c('RUVIII', 'RUVIII + Batch', 'RUVIII PBPS','Raw','Batch')) +
  theme(strip.background = element_blank(),strip.text =element_blank())

RUVPLOT

ggsave('/domino/datasets/local/RUV/paper/graphs/FDRT4poster.png',height = 5, width = 5)
