######################################
## Author: Chenzhuoyang
## Date: 2022-11-19
## Content: build network
#####################################

setwd("D:/WorkFile/HKUST/graph/network/")
dataDir <- "D:/WorkFile/Cornell_file/Thesis/TCGA/Validation/"
source("simFunctions.R")
CancerList <- list.files(dataDir)
CancerList <- CancerList[-c(5,24,30,31,38)]
saveRDS(CancerList, file = "cancerList.rds")

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
rm(Mu, Sigma)
gc()
cor_mat <- cor(t(expr))
#cor_mat <- WGCNA::cor(t(expr), nThreads = 7, method="pearson")

## adjacency matrix
#k <- softConnectivity(t(expr), type = "signed", power = 15, indent = 1)
#hist(k)
#scaleFreePlot(k)
#sampleTree = hclust(dist(t(expr)), method = "average") %>% plot()
#my_pca <- prcomp(t(expr), scale = FALSE, center = FALSE, retx = T)
#plot(my_pca$x[,1],my_pca$x[,2])

# Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
# sft = pickSoftThreshold.fromSimilarity(cor_mat, powerVector = powers, verbose = 2)

##output plot 1
# plot(x=sft$fitIndices[,1], 
#      y=-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1]+0.5, -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=0.9,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")

##output plot 2
hist(cor_mat[upper.tri(cor_mat)])

# {0.5*|1+cor|}^6
# A <- adjacency.fromSimilarity(cor_mat, type = "signed", power = 6)
# adjacency <- A
# adjacency[adjacency < 0] = 0
# adjacency[adjacency > 1] = 1

#should start from unthresholding one
connectivity=colSums(cor_mat) #abs
maxADJconst=1
Dhelp1 = matrix(connectivity,ncol=length(connectivity),nrow=length(connectivity))
denomTOM = pmin(as.dist(Dhelp1),as.dist(t(Dhelp1))) + as.dist(maxADJconst-cor_mat);
rm(Dhelp1); gc()
numTOM=as.dist(cor_mat %*% cor_mat +cor_mat);
TOM = as.matrix(numTOM/denomTOM)
rm(denomTOM, numTOM); gc()
TOM = as(TOM, "sparseMatrix")
sum(TOM > 0.1)
#cor_mat["ANK2|287","?|10357"]
#TOM["ANK2|287","?|10357"]
#cor_mat["AGPS|8540","?|10357"]
#TOM["AGPS|8540","?|10357"]


adj[adj > 0.1] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj)
network <- simplify(network)  # removes self-loops
results <- blockwiseModules(data, power=6, TOMType="signed", networkType="signed")
V(network)$color <- results$colors
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
seed(12345)
plot(network, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)



##save(processed, correlation, adjacency)
network <- list(expr, cor_mat2)
names(network) <- c("sample_Z", "cor_pearson")
save(network, file = paste0(cancer,"_coexpression_network.Rdata"))


##similarity
cor_mat1 <- cor_mat
cor_mat1[abs(cor_mat1)<0.5] = 0
cor_mat1 <- as(cor_mat1, "sparseMatrix")

cor_mat2 <- network$cor_pearson
cor_mat2[abs(cor_mat2)<0.5] = 0
cor_mat2 <- as(cor_mat2, "sparseMatrix")

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

#local distance
Hamming = sum(abs(con_cor_mat1 - con_cor_mat2)) / (nrow(con_cor_mat1) * (nrow(con_cor_mat1)-1))
Jaccard = sum(abs(con_cor_mat1 - con_cor_mat2)) / (sum((abs(con_cor_mat1) + abs(con_cor_mat2)) !=0)-nrow(con_cor_mat1))

#global distance
decomp1 <- eigen(cor_mat1, only.values = TRUE, symmetric = TRUE)
poly1 <- sort(decomp1$values)
decomp2 <- eigen(cor_mat2, only.values = TRUE, symmetric = TRUE)
poly2 <- sort(decomp2$values)

if(p == 1){
  return((sum(abs(poly1-poly2)))^(1/p))
} else{
  return((sum((poly1-poly2)^p))^(1/p))
}

library(RSpectra)
poly1 <- eigs_sym(cor_mat, k = 1000)$values
poly2 <- eigs_sym(cor_mat, k = 1000)$values

#mesoscale distance
Deltacon <- DeltaCon(con_cor_mat1, con_cor_mat2)


## ground truth
drug_overlap <- function(cancer1, cancer2){
  list1 <- filter(drug_aggr, `TCGA Classification` == cancer1) %>%
    .$`Drug Name`
  list2 <- filter(drug_aggr, `TCGA Classification` == cancer2) %>%
    .$`Drug Name`
  overlap <- length(intersect(list1,list2)) / ( sqrt(length(list1)) * sqrt(length(list2)) )
}

drug <- fread("C:/Users/Chenzhuoyang/Downloads/PANCANCER_IC_GDSC1.csv")
drug <- filter(drug, `TCGA Classification` %in% CancerList)
drug1 <- drug %>% group_by(`Drug Name`,`TCGA Classification`) %>%
  summarise(mean_IC50 = mean(IC50))
drug_aggr1 <- filter(drug1, mean_IC50 < 3)

drug <- fread("C:/Users/Chenzhuoyang/Downloads/PANCANCER_IC_GDSC2.csv")
drug <- filter(drug, `TCGA Classification` %in% CancerList)
drug2 <- drug %>% group_by(`Drug Name`,`TCGA Classification`) %>%
  summarise(mean_IC50 = mean(IC50))
drug_aggr2 <- filter(drug2, mean_IC50 < 3)

drug_aggr <- rbind(drug_aggr1, drug_aggr2)
drug_aggr <- drug_aggr[,c(2,1)]
drug_aggr <- arrange(drug_aggr, `TCGA Classification`)

drug_sim <- data.frame(x=rep(CancerList, each = length(CancerList)),
                       y=rep(CancerList, times = length(CancerList)))
drug_sim <- filter(drug_sim, match(x, CancerList)<match(y, CancerList))

drug_sim$overlap <- pbsapply(1:length(drug_sim), function(i) drug_overlap(drug_sim$x[i],drug_sim$y[i]))
unique(drug_aggr$`Drug Name`)
table(drug_aggr$`TCGA Classification`)
hist(drug_sim$overlap)

save(drug_sim, file = "drug_similarity.Rdata")

## LAML data
load("D:/WorkFile/Cornell_file/Thesis/TCGA/Validation/LAML/RSEM.Rdata")
m_Seg <- list(expr)
names(m_Seg) <- "RSEM"
save(m_Seg, file = "D:/WorkFile/Cornell_file/Thesis/TCGA/Validation/LAML/m_Seg_laml_firehose.Rdata")
