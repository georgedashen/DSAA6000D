######################################
## Author: Chenzhuoyang
## Date: 2022-11-19
## Content: build network
#####################################

setwd("D:/WorkFile/HKUST/graph/network/")
dataDir <- "D:/WorkFile/Cornell_file/Thesis/TCGA/Validation/"
source("D:/WorkFile/Cornell_file/Thesis/Functions.R")
CancerList <- list.files(dataDir)
CancerList <- CancerList[-c(5,24,30,31,38)]

##pipeline for cancer in CancerList
load(paste0(dataDir,cancer,"/m_Seg_",tolower(cancer),"_firehose.Rdata"))
load("D:/WorkFile/Cornell_file/Thesis/TCGA/Validation/ACC/m_Seg_acc_firehose.Rdata")
expr <- m_Seg$RSEM
rm(m_Seg)
gc()
##filter low expression
##filter low variation
cat("\nFiltering sdRSEM:\n")
sdRSEM <- pbapply(expr, 1, sd, cl=7)
meanRSEM <- pbapply(expr, 1, mean, cl=7)
ptsRSEM <- pbapply(expr, 1, function(x){sum(x>0)/length(x)}, cl=7)
geneFilter <- sdRSEM > 1 & meanRSEM > 1 & ptsRSEM > 0.5
cat(paste0(length(geneFilter)," genes in total\n"))
cat(paste0(sum(geneFilter)," genes left (", 
           sum(geneFilter)/length(geneFilter),"%)\n"))
expr <- expr[geneFilter,]
expr <- log(expr+0.1)

##correlation matrix
mu <- pbapply(expr, 2, mean, cl=7)
sigma <- pbapply(expr, 2, sd, cl=7)
Mu <- rep(mu, times=nrow(expr)) %>% 
  matrix(.,nrow=nrow(expr),ncol=length(mu),byrow=TRUE)
Sigma <- rep(sigma, times=nrow(expr)) %>% 
  matrix(.,nrow=nrow(expr),ncol=length(sigma),byrow=TRUE)
expr <- (expr - Mu) / Sigma
#cor_mat <- cor(t(expr))
library(WGCNA)
cor_mat <- WGCNA::cor(t(expr), nThreads = 7, method="pearson")

##adjacency matrix

##save(processed, correlation, adjacency)
network <- list(expr, cor_mat, adj_mat)
names(network) <- c("sample_Z, cor_pearson, adj")
save(network, file = paste0(cancer,"_coexpression_network.Rdata"))


