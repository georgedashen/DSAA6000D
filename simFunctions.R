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
library(pbapply) #faster apply and progress bar
library(Hmisc) #cor.test for matrix
library(dplyr) #dataframe manipulation

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
  
  Hamming <- sum(abs(mat1 - mat2)) / ((nrow(mat1)-1) * (nrow(mat1)-2))
  return(1-Hamming)
}




#### Jaccard ==============================
Jaccard_sim <- function(mat1, mat2){
  if(!checkMatrixConsistent(mat1, mat2)) 
    stop("Dimension inconsistent!")
  
  Jaccard <- sum(abs(mat1 - mat2)) / 
    sum((abs(mat1) + abs(mat2)) !=0)
  return(1-Jaccard)
}




#### spectral ==============================
spectral_dist <- function(mat1, mat2, p=1, k=min(nrow(mat1),nrow(mat2))){
  # if(nrow(mat1)==nrow(mat2) || ncol(mat1)==ncol(mat2)) 
  #   warning("Dimension not need to be consistent!")
  # 
  if(k > min(nrow(mat1),nrow(mat2))) stop("Please decrease k!")
  
  
  poly1 <- eigs_sym(mat1, k = k)$values
  poly2 <- eigs_sym(mat2, k = k)$values
  
  if(p == 1){
    return((sum(abs(poly1-poly2))))
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
inverse_lbp <- function(mat, power = 10){
  #LBP = linearized belief propagation
  N = nrow(mat)
  I = diag(1, N, N)
  I = as(I, "sparseMatrix")
  
  # Sparse degree-diagonal matrix, D[i,i] = sum(mat[i,])
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
  M = ch * mat  - ah * D
  
  # Calculate inverse of M
  inv_ = I
  mat_ = M
  pow = 1
  while(max(mat_) > 1e-09 && pow <= power) {
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




#### removeLowDegree ================
removeLowDegree <- function(mat1, mat2, degree = 1000){
  geneList <- rownames(mat1) #assume that mat1 and mat2 has same rownames
  g1 <- geneList[rowSums(mat1!=0)<degree]
  g2 <- geneList[rowSums(mat2!=0)<degree]
  return(intersect(g1, g2))
}




#### simList2Mat ==========================================
simList2Mat <- function(s){
  mat <- matrix(0, nrow = 33, ncol = 33)
  rownames(mat) = colnames(mat) = CancerTerm$Cohort[match(CancerList, CancerTerm$Cohort)]
  for(i in 1:nrow(drug_sim)){
    mat[drug_sim$x[i], drug_sim$y[i]] = s[i]
  }
  mat = mat + t(mat)
  diag(mat) = 1
  rownames(mat) = colnames(mat) = CancerTerm$Tissue[match(CancerList, CancerTerm$Cohort)]
  return(mat)
}


#### simHeatmap =====================
simHeatmap <- function(mat, Ncluster = 10, legend = NULL, ylab = NULL, titleFont = 20, columnFont = 15,
                       showRow = FALSE, gridBorder = 0){
  set.seed(667)
  Heatmap(mat,
          column_km = Ncluster, column_km_repeats = 500,
          row_km = Ncluster, row_km_repeats = 500,
          col = colorRamp2(c(0, 1), c("white", "red")),
          #cluster_rows = cluster_within_group(expr[,rownames(df)], as.factor(sorted)), 
          #cluster_columns = T, cluster_rows = T,
          show_column_names = T,show_row_names = showRow, show_row_dend = showRow,
          heatmap_legend_param = list(title = legend, direction = "vertical",
                                      title_position = "leftcenter-rot",at=c(0,1),legend_height = unit(3, "cm")),
          #top_annotation = top_anno,left_annotation = left_anno,
          row_title = NULL, column_title = ylab, column_title_side = "top",
          column_title_gp = gpar(fontsize = titleFont, fontface = "bold"),
          column_names_gp = grid::gpar(fontsize = columnFont),
          column_names_rot = 45,
          rect_gp = gpar(col = "white", lwd = gridBorder)) -> ht
  draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))
}




#### getVEnum ====================
getVEnum <- function(mat){
  return(c(nrow(mat),sum(mat!=0)))
}




#### getSignEnum =================
getSignEnum <- function(mat){
  return(c(sum(mat>0),sum(mat<0)))
}




#### getSignEnum =================
getWeigthDiff <- function(mat, upper = 1.0, lower = 0.9){
  gc()
  L <- mat>=lower
  U <- mat<=upper
  U[!L]=FALSE
  return(sum(U))
}




#### getDegree =================
getDegree <- function(mat){
  return(colSums(mat!=0))
}
