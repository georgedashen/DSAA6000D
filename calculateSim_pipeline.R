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
    s[i] = spectral_dist(matList[[drug_sim$x[i]]], matList[[drug_sim$y[i]]], p = p, k = k)
    print(paste0(format(i/length(s)*100,digits=3),"%"))
  }
  print(paste0(Sys.time() - start))
  s = (max(s) - s) + 1
  saveRDS(s, file = paste0(matrixType, "_Spectral.rds"))
  plot(order(s), order(drug_sim$overlap))
}

analyzeSpectral("unsign_unweight") #11.89 min
analyzeSpectral("sign_unweight") #11.25
analyzeSpectral("unsign_weight") #11.69
analyzeSpectral("sign_weight") #10.73


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

analyzeDeltaCon("unsign_unweight") # 2h 6m
analyzeDeltaCon("sign_unweight") # 2h 10m
analyzeDeltaCon("unsign_weight") # 2h 14m
analyzeDeltaCon("sign_weight") # 2h 4m


##size of mat
unsign_unweight <- readRDS("unsign_unweight_unified.rds")
Size <- pblapply(unsign_unweight, function(m) getVEnum(m)) %>% bind_rows()
Size <- t(Size) %>% as.data.frame()
colnames(Size) <- c("Num_nodes","Num_edges")
summary(Size$Num_edges)
#      Min.  1st Qu.   Median     Mean      3rd Qu.     Max. 
# 0.776168  2.764940  4.976554  7.751560  8.024516 32.474904 (Million)
saveRDS(Size, file = "processed_size.rds")


getConsistSize <- function(g1, g2) {
  geneList <- rownames(g1)
  removeGene <- removeLowDegree(g1, g2)
  removeGene <- match(removeGene, geneList) #reduce memory
  g1 <- g1[-removeGene, -removeGene]
  g2 <- g2[-removeGene, -removeGene]
  
  return(c(nrow(g1),sum(g1!=0),nrow(g2),sum(g2!=0)))
}

index <- 1:nrow(drug_sim)
names(index) <- as.character(1:nrow(drug_sim))
SizeDeltaCon <- pblapply(index, function(i){
  getConsistSize(unsign_unweight[[drug_sim$x[i]]],unsign_unweight[[drug_sim$y[i]]])}) %>% bind_rows()
SizeDeltaCon <- t(SizeDeltaCon) %>% as.data.frame()
colnames(SizeDeltaCon) <- c("Num_nodes_g1","Num_edges_g1","Num_nodes_g2","Num_edges_g2")
summary(SizeDeltaCon)
summary(c(SizeDeltaCon$Num_edges_g1,SizeDeltaCon$Num_edges_g2))
# node range from 60 to 13374, mean 4285
# edge range from 482 to 30M, mean 4M, median 1M
# It's reported that 1M edge graphs takes 160sec for DeltaCon
# One could reduce the accuracy of approximation by decrease repetition, decrease power and tolerance
