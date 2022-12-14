######################################
## Author: Chenzhuoyang
## Date: 2022-11-23
## Content: generate network for cancers
#####################################

setwd("D:/WorkFile/HKUST/graph/network/")
dataDir <- "D:/WorkFile/Cornell_file/Thesis/TCGA/Validation/"
source("simFunctions.R")
CancerList <- list.files(dataDir)
CancerList <- CancerList[-c(5,15,24,30,31,38)]
start <- "LGG"
stop <- "LGG"
#LAML separately processed due to missing m_Seg .Rdata

##pipeline for cancer in CancerList
unsign_unweight = TRUE # binary 0 1
unsign_weight = TRUE # [0,1]
sign_unweight = TRUE # -1, 0, 1
sign_weight = TRUE # [-1,1]
cor_cutoff = 0.5
tom_cutoff = 0.1
add_mode = FALSE

for(cancer in CancerList){
  if(match(cancer, CancerList) < match(start, CancerList)) next
  start = CancerList[match(start, CancerList)+1]
  if(match(cancer, CancerList) > match(stop, CancerList)) break
  if(!add_mode){
    cat(paste0("\n",cancer," analysis begins:\n"))
    load(paste0(dataDir,cancer,"/m_Seg_",tolower(cancer),"_firehose.Rdata"))
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
    rm(sdRSEM, meanRSEM, ptsRSEM, geneFilter);gc()
    
    ##correlation matrix
    cat("\nGenerating Z-scores:\n")
    mu <- pbapply(expr, 2, mean, cl=7)
    sigma <- pbapply(expr, 2, sd, cl=7)
    Mu <- rep(mu, times=nrow(expr)) %>% 
      matrix(.,nrow=nrow(expr),ncol=length(mu),byrow=TRUE)
    Sigma <- rep(sigma, times=nrow(expr)) %>% 
      matrix(.,nrow=nrow(expr),ncol=length(sigma),byrow=TRUE)
    expr <- (expr - Mu) / Sigma
    rm(Mu, Sigma)
    gc()
    
    cat("\nCalculating correlation:\n")
    R <- rcorr(t(expr),type="pearson")
    R$n <- NULL;gc()
    cat("\nGet cor_mat:\n")
    cor_mat <- R$r
    R$r <- NULL;gc()
    cat("\nGet p_mat:\n")
    p_adj <- R$P
    rm(R);gc()
    diag(p_adj) <- 0
    p_adj <- p_adj * (ncol(expr)-1) * (ncol(expr)-2) / 2 #bonferroni
    cor_mat[p_adj>0.05] = 0
    cor_mat <- as(cor_mat, "sparseMatrix")
    
    save(cor_mat, file = paste0(cancer,"_cor_mat.Rdata"))
    rm(expr, p_adj);gc()
    
    cat("\nCalculating adjacency matrix:\n")
    diag(cor_mat) <- 0
    network <- list()
    if(sign_weight){
      cor_mat_sw <- cor_mat
      cor_mat_sw[abs(cor_mat)<cor_cutoff] = 0
      cor_mat_sw <- as(cor_mat_sw, "sparseMatrix")
      network <- c(network,list(cor_mat_sw))
      names(network)[length(network)] <- "sign_weight"
    }
    if(sign_unweight){
      cor_mat_su <- cor_mat
      cor_mat_su[abs(cor_mat)<cor_cutoff] = 0
      cor_mat_su[cor_mat_su>0] = 1
      cor_mat_su[cor_mat_su<0] = -1
      cor_mat_su <- as(cor_mat_su, "sparseMatrix")
      network <- c(network,list(cor_mat_su))
      names(network)[length(network)] <- "sign_unweight"
    }
    if(unsign_unweight){
      cor_mat_uu <- cor_mat
      cor_mat_uu[abs(cor_mat)<cor_cutoff] = 0
      cor_mat_uu[cor_mat_uu!=0] = 1
      cor_mat_uu <- as(cor_mat_uu, "sparseMatrix")
      network <- c(network,list(cor_mat_uu))
      names(network)[length(network)] <- "unsign_unweight"
    }
    if(unsign_weight){
      cor_mat_uw <- cor_mat
      cor_mat_uw[abs(cor_mat)<cor_cutoff] = 0
      cor_mat_uw[cor_mat_uw<0] = abs(cor_mat_uw[cor_mat_uw<0])
      cor_mat_uw <- as(cor_mat_uw, "sparseMatrix")
      network <- c(network,list(cor_mat_uw))
      names(network)[length(network)] <- "unsign_weight"
    }
    
    save(network, file = paste0(cancer,"_cor_threshold",cor_cutoff,"_network.Rdata"))
    rm(cor_mat, network, cor_mat_su, cor_mat_sw, cor_mat_uu, cor_mat_uw)
    gc()
    
  } else{
    #retrieve cor_mats from the pre-defined cor_cutoff
    load(file = paste0(cancer,"_cor_threshold",cor_cutoff,"_network.Rdata"))
    
    #define the additional matrix here, sparseMatrix
    new_mat <- NULL
    new_mat <- as(new_mat, "sparseMatrix")
    #... other matrix
    #... sparse matrix
    
    network <- c(network,list(new_mat)) #append new matrix
    names(network)[length(network)] <- "newName" #assign a new name
    # ...append
    # ...assign name
    
    
    save(network, file = paste0(cancer,"_cor_threshold",cor_cutoff,"_network.Rdata"))
    gc()
  }  
  
  
}
