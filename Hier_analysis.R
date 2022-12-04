######################################
## Author: Chenzhuoyang
## Date: 2022-11-30
## Content: Hierarchical relationship analysis
#####################################

setwd("D:/WorkFile/HKUST/graph/network/")
dataDir <- "D:/WorkFile/Cornell_file/Thesis/TCGA/Validation/"
source("simFunctions.R")
load("drug_similarity.Rdata")
library(PerformanceAnalytics)
library(ComplexHeatmap)
library(circlize)

# CancerList
CancerTerm <- fread("CancerTerm.csv", data.table = F)
CancerList <- readRDS("cancerList.rds")

# geneList
for(cancer in CancerList){
  load(paste0(dataDir,cancer,"/m_Seg_",tolower(cancer),"_firehose.Rdata"))
  if(cancer == "ACC"){
    geneList <- rownames(m_Seg$RSEM)
  }else{
    geneList <- intersect(rownames(m_Seg$RSEM), geneList)
  }
}
saveRDS(geneList, file = "geneList.rds")



#baseline
dataDir <- "D:/WorkFile/Cornell_file/Thesis/TCGA/Validation/"
geneList <- readRDS("geneList.rds")

mean_expr_mat <- data.frame(gene = geneList)
for (cancer in CancerList) {
  load(paste0(dataDir,cancer,"/m_Seg_",tolower(cancer),"_firehose.Rdata"))
  expr <- m_Seg$RSEM[geneList,]
  expr <- log(expr+0.1)
  sizefactor <- pbapply(expr, 2, sum, cl=7)
  SizeFactor <- rep(sizefactor, times=nrow(expr)) %>% 
    matrix(.,nrow=nrow(expr),ncol=length(sizefactor),byrow=TRUE)
  expr <- expr / SizeFactor
  mean_expr_mat$new <- rowMeans(expr)
  colnames(mean_expr_mat)[ncol(mean_expr_mat)] <- cancer
}

rownames(mean_expr_mat) <- mean_expr_mat$gene
mean_expr_mat <- mean_expr_mat[,-1]
save(mean_expr_mat, file = "mean_expr_mat.Rdata")

# hist(mean_expr_mat$ACC)
#colnames(mean_expr_mat) <- CancerTerm$Cohort[match(colnames(mean_expr_mat), CancerTerm$Tissue)]
colnames(mean_expr_mat) <- CancerTerm$Tissue[match(colnames(mean_expr_mat), CancerTerm$Cohort)]
png("meanExpr_TCGAsim.png",width = 10, height = 8, units = "in", res = 300)
simHeatmap(cor(mean_expr_mat, method = "spearman"), Ncluster = 5, legend = "Expression Similarity", ylab = "Correlation of Mean Expression\nacross Cancers")
dev.off()


## test
# Heatmap
library(ComplexHeatmap)
library(circlize)
# top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
# color_v=ColorPalette[1:N_group] #subcluster
# names(color_v)=sort(unique(df$leiden_subcluster))
# color_h=ColorPalette[1:N_class]
# names(color_h)=sort(unique(df$class))
set.seed(670)
Heatmap(mat, #绘图数据的CB顺序和注释CB顺序保持一致
        #column_km = 10, column_km_repeats = 1000,
        split = 5,
        #row_km = 5, row_km_repeats = 1000,
        col = colorRamp2(c(0, 1), c("white", "red")),
        #cluster_rows = cluster_within_group(expr[,rownames(df)], as.factor(sorted)), 
        cluster_columns = T, cluster_rows = T,
        show_column_names = T,show_row_names = F,
        heatmap_legend_param = list(title = "Jaccard Similarity", direction = "vertical",
                                    title_position = "leftcenter-rot",at=c(0,1),legend_height = unit(3, "cm")),
        #top_annotation = top_anno,left_annotation = left_anno, #添加注释
        row_title = NULL, column_title = "Cancer")


#Jaccard
s <- readRDS("unsign_unweight_Jaccard.rds")
drug_sim$uu_Jaccard <- s
uu <- simList2Mat(s)

s <- readRDS("unsign_weight_Jaccard.rds")
drug_sim$uw_Jaccard <- s
uw <- simList2Mat(s)

s <- readRDS("sign_unweight_Jaccard.rds")
drug_sim$su_Jaccard <- s
su <- simList2Mat(s)

s <- readRDS("sign_weight_Jaccard.rds")
drug_sim$sw_Jaccard <- s
sw <- simList2Mat(s)

png("unsign_unweigt_Jaccard_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(uu, Ncluster = 5, legend = "Jaccard Similarity", ylab = "Jaccard: Unsign Unweight Adj", columnFont = 18)
dev.off()

png("unsign_weigt_Jaccard_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(uw, Ncluster = 5, legend = "Jaccard Similarity", ylab = "Jaccard: Unsign Weight Adj", columnFont = 18)
dev.off()

png("sign_unweigt_Jaccard_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(su, Ncluster = 5, legend = "Jaccard Similarity", ylab = "Jaccard: Sign Unweight Adj", columnFont = 18)
dev.off()

png("sign_weigt_Jaccard_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(sw, Ncluster = 5, legend = "Jaccard Similarity", ylab = "Jaccard: Sign Weight Adj", columnFont = 18)
dev.off()

#
png("Jaccard_Adj_heatmap.png", width = 6, height = 5, units = "in", res = 300)
simHeatmap(cor(drug_sim[c(4:7)], method = "spearman"), Ncluster = 2, showRow = T, 
           legend = "Spearman rho", ylab = "Correlation of Jaccard Similarity \nacorss AdjMat Types")
dev.off()

png("Jaccard_Adj_corchat.png", width = 10, height = 8, units = "in", res = 300)
chart.Correlation(drug_sim[,c(4:7)], method = "spearman")
dev.off()

png("Jaccard_Adj_corchat_pearson.png", width = 10, height = 8, units = "in", res = 300)
chart.Correlation(drug_sim[,c(4:7)], method = "pearson")
dev.off()

#Spectral
s <- readRDS("unsign_unweight_Spectral.rds")
drug_sim$uu_Spectral <- s
uu <- simList2Mat(s)

s <- readRDS("unsign_weight_Spectral.rds")
drug_sim$uw_Spectral <- s
uw <- simList2Mat(s)

s <- readRDS("sign_unweight_Spectral.rds")
drug_sim$su_Spectral <- s
su <- simList2Mat(s)

s <- readRDS("sign_weight_Spectral.rds")
drug_sim$sw_Spectral <- s
sw <- simList2Mat(s)

png("unsign_unweigt_Spectral_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(uu, Ncluster = 5, legend = "Spectral Similarity", ylab = "Spectral: Unsign Unweight Adj", columnFont = 18)
dev.off()

png("unsign_weigt_Spectral_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(uw, Ncluster = 5, legend = "Spectral Similarity", ylab = "Spectral: Unsign Weight Adj", columnFont = 18)
dev.off()

png("sign_unweigt_Spectral_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(su, Ncluster = 5, legend = "Spectral Similarity", ylab = "Spectral: Sign Unweight Adj", columnFont = 18)
dev.off()

png("sign_weigt_Spectral_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(sw, Ncluster = 5, legend = "Spectral Similarity", ylab = "Spectral: Sign Weight Adj", columnFont = 18)
dev.off()

#
png("Spectral_Adj_heatmap.png", width = 6, height = 5, units = "in", res = 300)
simHeatmap(cor(drug_sim[c(8:11)], method = "spearman"), Ncluster = 2, showRow = T, 
           legend = "Spearman rho", ylab = "Correlation of Spectral-dist Similarity\nacorss AdjMat Types")
dev.off()

png("Spectral_Adj_corchat.png", width = 10, height = 8, units = "in", res = 300)
chart.Correlation(drug_sim[,c(8:11)], method = "spearman")
dev.off()

png("Spectral_Adj_corchat_pearson.png", width = 10, height = 8, units = "in", res = 300)
chart.Correlation(drug_sim[,c(8:11)], method = "pearson")
dev.off()

# DeltaCon
s <- readRDS("unsign_unweight_DeltaCon.rds")
drug_sim$uu_DeltaCon <- s
uu <- simList2Mat(s)

s <- readRDS("unsign_weight_DeltaCon.rds")
drug_sim$uw_DeltaCon <- s
uw <- simList2Mat(s)

s <- readRDS("sign_unweight_DeltaCon.rds")
drug_sim$su_DeltaCon <- s
su <- simList2Mat(s)

s <- readRDS("sign_weight_DeltaCon.rds")
drug_sim$sw_DeltaCon <- s
sw <- simList2Mat(s)

png("unsign_unweigt_DeltaCon_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(uu, Ncluster = 5, legend = "DeltaCon Similarity", ylab = "DeltaCon: Unsign Unweight Adj", columnFont = 18)
dev.off()

png("unsign_weigt_DeltaCon_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(uw, Ncluster = 5, legend = "DeltaCon Similarity", ylab = "DeltaCon: Unsign Weight Adj", columnFont = 18)
dev.off()

png("sign_unweigt_DeltaCon_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(su, Ncluster = 5, legend = "DeltaCon Similarity", ylab = "DeltaCon: Sign Unweight Adj", columnFont = 18)
dev.off()

png("sign_weigt_DeltaCon_TCGAsim.png", width = 10, height = 8, units = "in", res = 300)
simHeatmap(sw, Ncluster = 5, legend = "DeltaCon Similarity", ylab = "DeltaCon: Sign Weight Adj", columnFont = 18)
dev.off()

#
png("DeltaCon_Adj_heatmap.png", width = 6, height = 5, units = "in", res = 300)
simHeatmap(cor(drug_sim[c(12:15)], method = "spearman"), Ncluster = 2, showRow = T, 
           legend = "Spearman rho", ylab = "Correlation of DeltaCon Similarity \nacorss AdjMat Types")
dev.off()

png("DeltaCon_Adj_corchat.png", width = 10, height = 8, units = "in", res = 300)
chart.Correlation(drug_sim[,c(12:15)], method = "spearman")
dev.off()

#
png("All_Adj_heatmap.png", width = 6, height = 5, units = "in", res = 300)
simHeatmap(cor(drug_sim[c(4:15)], method = "spearman"), Ncluster = 2, showRow = T, 
           legend = "Spearman rho", ylab = "Correlation of All Similarity \nacorss AdjMat Types")
dev.off()

png("All_Adj_corchat.png", width = 10, height = 8, units = "in", res = 300)
chart.Correlation(drug_sim[,c(4:15)], method = "spearman")
dev.off()

#
corMeanExpr <- cor(mean_expr_mat, method = "spearman")
colnames(corMeanExpr) = rownames(corMeanExpr) = CancerTerm$Cohort[match(colnames(corMeanExpr), CancerTerm$Tissue)]
drug_sim$meanExpr <- pbsapply(1:nrow(drug_sim), function(i) corMeanExpr[drug_sim$x[i], drug_sim$y[i]])

png("Allexpr_Adj_heatmap.png", width = 6, height = 5, units = "in", res = 300)
simHeatmap(cor(drug_sim[c(4:16)], method = "spearman"), Ncluster = 2, showRow = T, 
           legend = "Spearman rho", ylab = "Correlation of All Similarity \nacorss AdjMat Types")
dev.off()

png("Allexpr_Adj_corchat.png", width = 10, height = 8, units = "in", res = 300)
chart.Correlation(drug_sim[,c(4:16)], method = "spearman")
dev.off()

#conclusion 1:
### For Jaccard, sign/unsigned doesn't matter
### For Spectral, sign/unsigned doesn't matter
### For DeltaCon, both sign/unsigned matter, and sign is the major difference
### This is natural, because deltacon measures information flow in the vicinity
### the sign of deltacon determine the effect of change and weight determine the amount of change

#conclusion 2:
### DeltaCon is the most different methods compared to others, but what these similarities
### represent still not clear, the expression matrix similarity itself already well represents
### the cell origin and tissue characteristic.

timeConsume <- data.frame(method = rep(c("Jaccard", "Spectral(eigen50)", "DeltaCon(power4)"), each = 4),
                          adjmat = rep(c("unsign_unweight", "sign_unweight", "unsign_weight", "sign_weight"), times = 3))
timeConsume$method <- factor(timeConsume$method, levels = c("Jaccard","Spectral(eigen50)","DeltaCon(power4)"))
timeConsume$adjmat <- factor(timeConsume$adjmat, levels = c("sign_unweight", "sign_weight", "unsign_unweight", "unsign_weight"))
timeConsume$elapsed <- c(30.1, 31.5, 33.0, 32.7, 40.6, 41.6, 39.4, 40.2, 126, 130, 134, 124)

library(ggplot2)
library(Seurat)
ggplot(timeConsume, aes(adjmat, elapsed, fill = adjmat)) + geom_bar(stat="identity") + 
  geom_text(aes(label = elapsed), vjust = -0.5, hjust = 0.5, ) +
  facet_wrap(~method) + ylab("elapsed(min)") + 
  scale_y_continuous(breaks = c(0,25,50,75,100,125,150),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_bw(base_size = 20) + theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) + 
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"cm"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,0,0,0)) +
  RotatedAxis() + scale_fill_manual(values=c("#fdae61","#d7191c","#abd9e9","#2c7bb6"))
ggsave("elapsed_time.png", width = 9, height = 6)
