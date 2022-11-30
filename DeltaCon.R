####### Please Check github.com/bxshi/rdsg for new updates, this gist will not be change in the future #######
## However, you can still use this since it is correct so far ##

# Port of original DeltaCon to R
# Baoxu(Dash) Shi
# Data Sciense Group
# University of Notre Dame
#
# Citation
#   D. Koutra, J. T. Vogelstein, and C. Faloutsos:
#   DeltaCon: A Principled Massive-Graph Similarity Function.
#   SIAM 2013: 162–170.
#   D. Koutra, T. Ke, U. Kang, D. H. Chau, H. K. Pao, C. Faloutsos:
#   Unifying Guilt-by-Association Approaches: Theorems and Fast Algorithms.
#   ECML/PKDD (2) 2011: 245-260

library(Matrix)
library(SparseM)
library(pracma)

.MAX_POWER = 10
.p = 0.51

output_time <- function(debug, tim, s) {
  if(debug) {
    print(paste(s, "user time:", tim[1], "system time:", tim[2], "elapsed time:", tim[3]))
  }
}

inverse_lbp <- function(graph, priors=NULL, power = 10, times = 10, debug = FALSE) {
  
  # Sparse identity matrix
  nnodes = nrow(graph)
  I = diag(1, nnodes, nnodes)
  I = as(I, "sparseMatrix")
  
  # Sparse degree-diagonal matrix, D[i,i] = sum(mat[i,])
  x <- Matrix::rowSums(graph)
  D <- sparseMatrix(c(1:nnodes), c(1:nnodes), x = x, dims = c(nnodes, nnodes))
  
  # Compute about-half homophily factor to guarantee covergence
  c1 = sum(D) + 2
  c2 = sum(D^2) - 1
  h_h = sqrt((-c1 + sqrt(c1^2 + 4 * c2)) / (8 * c2))
  
  # Compute constant ah and ch
  ah = 4 * h_h^2 / (1 - 4 * h_h^2)
  ch = 2 * h_h / (1 - 4 * h_h^2)
  
  # Invert matrices M1 and M2
  M = NULL
  M = ch * graph  - ah * D
  
  # Calculate inverse of M
  if (is.null(priors)) {
    inv_ = I
    mat_ = M
    pow = 1
    tim <- system.time({
      while(max(mat_) > 1e-09 && pow < power) {
        inv_ = inv_ + mat_
        mat_ = mat_ %*% M
        pow = pow + 1
      }
    })
    output_time(debug, tim, "Invert of matrix")
    return(inv_)
  } else {
    final_mat <- NULL
    tim <- system.time({
      for(i in c(1:times)) {
        inv_ = matrix(priors[, i], nnodes, 1)
        mat_ = matrix(priors[, i], nnodes, 1)
        pow = 1
        while(max(mat_) > 1e-09 && pow < power) {
          mat_ = M %*% mat_
          inv_ = inv_ + mat_
          pow = pow + 1
        }
        if (i == 1) {
          final_mat <- matrix(inv_)
        } else {
          final_mat <- cbind(final_mat, matrix(inv_))
        }
      }
    })
    output_time(debug, tim, "Approximate invert of matrix")
    return(final_mat)
  }
  
}

init_priors_percent <- function(percent, nnodes) {
  times <- ceiling(1 / percent)
  init_nodes <- floor(percent * nnodes)
  
  rand_vector <- runif(nnodes) 
  
  rand_mat <- repmat(matrix(rand_vector, nnodes, 1), 1, times)
  
  for(i in c(1:times)) {
    rand_mat[rand_mat[,i] >= (i-1)*percent & rand_mat[,i] < i *percent, i] <- 1
    rand_mat[rand_mat[,i] != 1, i] <- 0
  }
  
  return(rand_mat)
}

delta_con <- function(g1, g2, method = "naive", power = 10,
                      percent = 0.1, removeDegree = 1000, debug = FALSE) {
  
  removeGene <- removeLowDegree(g1, g2, degree = removeDegree)
  removeGene <- match(removeGene, geneList) #reduce memory
  g1 <- g1[-removeGene, -removeGene]
  g1 <- g1*0.5+0.5 #positive
  g2 <- g2[-removeGene, -removeGene]
  g2 <- g2*0.5+0.5 #positive
  
  nnodes = nrow(g1)
  if (method == "fast") {
    # Number of groups to be initialized
    ngroups <- ceiling(1 / percent)
    
    repetitions <- 5
    t_all <- rep(0, repetitions)
    sim <- rep(0, repetitions)
    for(i in c(1:repetitions)) {
      priors <- NULL
      tim <- system.time({
        priors <- init_priors_percent(percent, nnodes)
      })
      output_time(debug, tim, "Calculate priors")
      
      inv1 <- inverse_lbp(g1, priors, power = power, ngroups, debug = debug)# * (.p - 0.5)
      inv2 <- inverse_lbp(g2, priors, power = power, ngroups, debug = debug)# * (.p - 0.5)
      sim[i] <- 1 / (1 + sqrt(sum( (sqrt(inv1) - sqrt(inv2))^2 )))
    }
    delta_con <- mean(sim)
    return(delta_con)
  } else {
    # Naive FaBP
    inv1 <- inverse_lbp(g1, debug = debug) * (.p - 0.5)
    inv2 <- inverse_lbp(g2, debug = debug) * (.p - 0.5)
    
    # Compute DeltaCon similarity score
    delta_con <- 1 / (1 + sqrt(sum( (sqrt(inv1) - sqrt(inv2))^2 )))
    return(delta_con)
  }
}

delta_con_example <- function() {
  g1 <- as.data.frame(cbind(c(1,1,2,2,3,3,4,8,5),
                            c(2,3,4,5,6,7,8,9,10)))
  g2 <- as.data.frame(cbind(c(1,1,2,2,3,3,4,8,5,9,10,5,6),
                            c(2,3,4,5,6,7,8,9,10,11,11,12,12)))
  
  print(delta_con(g1,g1,max(g1,g2), debug = TRUE))
  print(delta_con(g1,g2,max(g1,g2), debug = TRUE))
  print(delta_con(g2,g1,max(g1,g2), debug = TRUE))
  print(delta_con(g2,g2,max(g1,g2), debug = TRUE))
  
  print(delta_con(g1,g1,max(g1,g2), node_list = c(1,2,3,4,5,6,7,8,9,10), method="fast", debug = TRUE))
  print(delta_con(g1,g2,max(g1,g2), node_list = c(1,2,3,4,5,6,7,8,9,10), method="fast", debug = TRUE))
  print(delta_con(g2,g1,max(g1,g2), method="fast", debug = TRUE))
  print(delta_con(g2,g2,max(g1,g2), method="fast", debug = TRUE))
  
  print(delta_con(g1,g1,max(g1,g2), method="fast", symmetrical = FALSE))
  print(delta_con(g1,g2,max(g1,g2), method="fast", symmetrical = FALSE))
  print(delta_con(g2,g1,max(g1,g2), method="fast", symmetrical = FALSE))
  print(delta_con(g2,g2,max(g1,g2), method="fast", symmetrical = FALSE))
  
  print(delta_con(g1,g1,max(g1,g2), symmetrical = FALSE))
  print(delta_con(g1,g2,max(g1,g2), symmetrical = FALSE))
  print(delta_con(g2,g1,max(g1,g2), symmetrical = FALSE))
  print(delta_con(g2,g2,max(g1,g2), symmetrical = FALSE))
  
  g3 <- as.data.frame(cbind(c(2,3),
                            c(1,1)))
  g4 <- as.data.frame(cbind(c(2,3),
                            c(4,4)))
  
  g5 <- as.data.frame(cbind(c(2,3,3,1,5,6,6),
                            c(1,1,2,4,4,4,5)))
  g6 <- as.data.frame(cbind(c(2,3,3,5,6,6),
                            c(1,1,2,4,4,5)))
}