######################################
## Author: Chenzhuoyang
## Date: 2022-11-24
## Content: generate similarity for cancers
#####################################

setwd("D:/WorkFile/HKUST/graph/network/")
dataDir <- "D:/WorkFile/Cornell_file/Thesis/TCGA/Validation/"
source("simFunctions.R")
source("DeltaCon.R")
load("drug_similarity.Rdata")
cor_cutoff = 0.5

#### load all
networkList = list()
for (i in CancerList) {
  load(file = paste0(i,"_cor_threshold",cor_cutoff,"_network.Rdata"))
  networkList <- c(list(network), networkList)
}
names(networkList) <- CancerList
rm(network); gc()

save(networkList, file = "networkList.Rdata")

geneList=c()
for(n in networkList){
  geneList <- c(rownames(n$unsign_unweight), geneList)
}
geneList <- unique(geneList)
geneList <- sort(geneList)

matrixUnify <- function(test){
  addName <- setdiff(geneList, rownames(test))
  Nold <- nrow(test)
  N <- length(geneList)
  test <- rbind(test, Matrix(0,sparse = T, nrow = N-Nold,ncol = Nold))
  test <- cbind(test, Matrix(0,sparse = T, nrow = N,ncol = N-Nold))
  rownames(test)[(Nold+1):N] = colnames(test)[(Nold+1):N] = addName
  test <- test[order(rownames(test)),order(colnames(test))]
  return(test)
}


## unify
matList <- c("unsign_unweight","unsign_weight","sign_unweight","sign_weight")
for(m in matList[-1]){
  load("networkList.Rdata")
  for (c in CancerList) {
    for(adj in matList[-match(m,matList)]){
      networkList[[c]][[adj]] = NULL
    }
    
  }
  gc()
  for (c in CancerList) {
    networkList[[c]] <- matrixUnify(networkList[[c]][[m]])
  }
  
  saveRDS(networkList, file = paste0(m,"_unified.rds"))
  rm(networkList);gc()
}


## Jaccard
#load("drug_similarity.Rdata")
unsign_unweight <- readRDS("unsign_unweight_unified.rds")
s = drug_sim$overlap
start <- Sys.time()
for(i in 1:length(s)){
    s[i] = Jaccard_sim(unsign_unweight[[drug_sim$x[i]]], unsign_unweight[[drug_sim$y[i]]])
    print(paste0(format(i/length(s)*100,digits=2),"%"))
}
paste0(Sys.time() - start) #30.045 min
saveRDS(s, file = "unsign_unweight_Jaccard.rds")
plot(order(s), order(drug_sim$overlap))

##
sign_unweight <- readRDS("sign_unweight_unified.rds")
s = drug_sim$overlap
start <- Sys.time()
for(i in 1:length(s)){
  s[i] = Jaccard_sim(sign_unweight[[drug_sim$x[i]]], sign_unweight[[drug_sim$y[i]]])
  print(paste0(format(i/length(s)*100,digits=3),"%"))
}
paste0(Sys.time() - start) #31.53 min
saveRDS(s, file = "sign_unweight_Jaccard.rds")
plot(order(s), order(drug_sim$overlap))

##
sign_weight <- readRDS("sign_weight_unified.rds")
s = drug_sim$overlap
start <- Sys.time()
for(i in 1:length(s)){
  s[i] = Jaccard_sim(sign_weight[[drug_sim$x[i]]], sign_weight[[drug_sim$y[i]]])
  print(paste0(format(i/length(s)*100,digits=3),"%"))
}
paste0(Sys.time() - start) #32.69 min
saveRDS(s, file = "sign_weight_Jaccard.rds")
plot(order(s), order(drug_sim$overlap))

##
unsign_weight <- readRDS("unsign_weight_unified.rds")
s = drug_sim$overlap
start <- Sys.time()
for(i in 1:length(s)){
  s[i] = Jaccard_sim(unsign_weight[[drug_sim$x[i]]], unsign_weight[[drug_sim$y[i]]])
  print(paste0(format(i/length(s)*100,digits=3),"%"))
}
paste0(Sys.time() - start) #33.008 min
saveRDS(s, file = "unsign_weight_Jaccard.rds")
plot(order(s), order(drug_sim$overlap))

## Spectral
analyzeSpectral <- function(matrixType, p = 2, k = 10){
  if(!matrixType %in% c("unsign_unweight","unsign_weight","sign_unweight","sign_weight"))
    stop("Unknown matrix type!")
  matList <- readRDS(paste0(matrixType, "_unified.rds"))
  s = drug_sim$overlap
  start <- Sys.time()
  for(i in 1:length(s)){
    s[i] = Jaccard_sim(matList[[drug_sim$x[i]]], matList[[drug_sim$y[i]]])
    print(paste0(format(i/length(s)*100,digits=3),"%"))
  }
  print(paste0(Sys.time() - start))
  saveRDS(s, file = paste0(matrixType, "_Spectral.rds"))
  plot(order(s), order(drug_sim$overlap))
}

analyzeSpectral("unsign_unweight") #33.069
analyzeSpectral("sign_unweight") #32.49
analyzeSpectral("unsign_weight") #33.66
analyzeSpectral("sign_weight") #33.30


## DeltaCon
analyzeDeltaCon <- function(matrixType){
  if(!matrixType %in% c("unsign_unweight","unsign_weight","sign_unweight","sign_weight"))
    stop("Unknown matrix type!")
  matList <- readRDS(paste0(matrixType, "_unified.rds"))
  s = drug_sim$overlap
  start <- Sys.time()
  for(i in 1:length(s)){
    s[i] = delta_con(matList[[drug_sim$x[i]]], matList[[drug_sim$y[i]]], 
                     method="fast", power = 5, 
                     percent = 1, removeDegree = 1000)
    print(paste0(format(i/length(s)*100,digits=3),"%"))
  }
  print(paste0(Sys.time() - start))
  saveRDS(s, file = paste0(matrixType, "_DeltaCon.rds"))
  plot(order(s), order(drug_sim$overlap))
}

analyzeDeltaCon("unsign_unweight") # 2h 11m
analyzeDeltaCon("sign_unweight") # 2h 17m
analyzeDeltaCon("unsign_weight") # 2h 23m
analyzeDeltaCon("sign_weight") # 2h 6m