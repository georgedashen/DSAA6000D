#########################
## Author: Chenzhuoyang
## Date: 2022-11-21
## Content: functions for investigate graph similarity
########################

library(MASS) #matrix manipulation
library(data.table) #fread csv
library(magrittr) #pipe operation
library(WGCNA) #co-expression network
library(Matrix) #sparse matrix
library(RSpectra) #high efficiency svd and eigens

#### matrixConsistent =====================
checkMatrixConsistent <- function(mat1, mat2){
  c <- nrow(mat1)==nrow(mat2) && ncol(mat1)==ncol(mat2) && nrow(mat1)==ncol(mat1) && nrow(mat2)==ncol(mat2)
  return(c)
}




#### makeConsistent =======================
makeConsistent <- function(mat1, mat2){
  cor_mat1 <- as(mat1, "sparseMatrix")
  cor_mat2 <- as(mat2, "sparseMatrix")
  rm(mat1,mat2)
  
  aGene <- setdiff(colnames(cor_mat2), colnames(cor_mat1))
  con_cor_mat1 <- matrix(0, nrow = ncol(cor_mat1) + length(aGene), 
                         ncol = ncol(cor_mat1) + length(aGene))
  colnames(con_cor_mat1) = rownames(con_cor_mat1) = c(rownames(cor_mat1), aGene)
  con_cor_mat1 <- as(con_cor_mat1, "sparseMatrix")
  con_cor_mat1[1:nrow(cor_mat1),1:ncol(cor_mat1)] = cor_mat1
  diag(con_cor_mat1) <- 1
  
  bGene <- setdiff(colnames(cor_mat1), colnames(cor_mat2))
  con_cor_mat2 <- matrix(0, nrow = ncol(cor_mat2) + length(bGene), 
                         ncol = ncol(cor_mat2) + length(bGene))
  colnames(con_cor_mat2) = rownames(con_cor_mat2) = c(rownames(cor_mat2), bGene)
  con_cor_mat2 <- as(con_cor_mat2, "sparseMatrix")
  con_cor_mat2[1:nrow(cor_mat2),1:ncol(cor_mat2)] = cor_mat2
  diag(con_cor_mat2) <- 1
  
  con_cor_mat2[match(rownames(con_cor_mat1),rownames(con_cor_mat2)),
               match(rownames(con_cor_mat1),rownames(con_cor_mat2))] -> con_cor_mat2
  rownames(con_cor_mat2) = colnames(con_cor_mat2) = rownames(con_cor_mat1)
  
  con_mat <- list(con_cor_mat1, con_cor_mat2)
  names(con_mat) <- c("con_cor_mat1","con_cor_mat2")
  
  return(con_mat)
}



#### Hamming ==============================
Hamming_sim <- function(mat1, mat2){
  if(!checkMatrixConsistent(mat1,mat2)) 
    stop("Dimension inconsistent!")
  
  Hamming <- sum(abs(mat1 - mat2)) / (nrow(mat1) * (nrow(mat1)-1))
  return(1-Hamming)
}




#### Jaccard ==============================
Jaccard_sim <- function(mat1, mat2){
  if(!checkMatrixConsistent(mat1, mat2)) 
    stop("Dimension inconsistent!")
  
  Jaccard <- sum(abs(mat1 - mat2)) / 
    (sum((abs(mat1) + abs(mat2)) !=0)-nrow(mat1))
  return(1-Jaccard)
}




#### spectral ==============================
spectral_dist <- function(mat1, mat2, p=1, k=min(nrow(mat1),nrow(mat2))){
  if(nrow(mat1)==nrow(mat2) || ncol(mat1)==ncol(mat2)) 
    warning("Dimension not need to be consistent!")
  
  if(k > min(nrow(mat1),nrow(mat2))) stop("Please decrease k!")
  
  
  poly1 <- eigs_sym(mat1, k = k)$values
  poly2 <- eigs_sym(mat2, k = k)$values
  
  if(p == 1){
    return((sum(abs(poly1-poly2)))^(1/p))
  } else{
    return((sum((poly1-poly2)^p))^(1/p))
  }
}




#### TOMsimilarity =====================
TOMsimilarity <- function(cor_mat, type="unsigned"){
  if(!type%in%c("signed","unsigned"))
    stop("Please select correct type!")
  
  if(type == "unsigned"){
    cor_mat <- abs(cor_mat)
    connectivity=colSums(cor_mat)
    Dhelp1 = matrix(connectivity,ncol=length(connectivity),nrow=length(connectivity))
    denomTOM = pmin(as.dist(Dhelp1),as.dist(t(Dhelp1))) + as.dist(1-cor_mat);
    rm(Dhelp1); gc()
    numTOM=as.dist(cor_mat %*% cor_mat +cor_mat);
    TOM = as.matrix(numTOM/denomTOM)
    rm(denomTOM, numTOM); gc()
    TOM = as(TOM, "sparseMatrix")
    return(TOM)
  } else{
    connectivity=colSums(abs(cor_mat))
    Dhelp1 = matrix(connectivity,ncol=length(connectivity),nrow=length(connectivity))
    denomTOM = pmin(as.dist(Dhelp1),as.dist(t(Dhelp1))) + as.dist(1-abs(cor_mat));
    rm(Dhelp1); gc()
    numTOM=as.dist(cor_mat %*% cor_mat +cor_mat);
    TOM = as.matrix(numTOM/denomTOM)
    rm(denomTOM, numTOM); gc()
    TOM = as(TOM, "sparseMatrix")
    return(TOM)
  }
}




#### inverse_lbp =======================
inverse_lbp <- function(mat){
  #LBP = linearized belief propagation
  mat <- as(mat, "sparseMatrix")
  N = nrow(mat)
  I = diag(1, N, N)
  I = as(I, "sparseMatrix")
  
  # Sparse degree-diagonal matrix, D[i,i] = sum(graph[i,])
  x <- rowSums(mat, sparseResult = TRUE)
  D <- sparseMatrix(c(1:N), c(1:N), x = x, dims = c(N, N))
  
  # Compute about-half homophily factor to guarantee covergence
  c1 = sum(D) + 2
  c2 = sum(D^2) - 1
  h_h = sqrt((-c1 + sqrt(c1^2 + 4 * c2)) / (8 * c2))
  
  # Compute constant ah and ch
  ah = 4 * h_h^2 / (1 - 4 * h_h^2)
  ch = 2 * h_h / (1 - 4 * h_h^2)
  
  # Invert matrices M1 and M2
  M = ch * graph  - ah * D
  
  # Calculate inverse of M
  inv_ = I
  mat_ = M
  pow = 1
  while(max(mat_) > 1e-09 && pow < 10) {
    inv_ = inv_ + mat_
    mat_ = mat_ %*% M
    pow = pow + 1
  }

  return(inv_)
}




#### DeltaCon =======================
DeltaCon <- function(mat1, mat2){
  if(!checkMatrixConsistent(mat1, mat2)) 
    stop("Dimension inconsistent!")
  inv1 <- inverse_lbp(mat1)
  inv2 <- inverse_lbp(mat2)
  delta_con <- 1 / (1 + sqrt(sum( (sqrt(inv1) - sqrt(inv2))^2 )))
  
  return(delta_con)
}


