# DSAA6000D
This is for the final project for DSAA6000D at https://github.com/georgedashen/DSAA6000D

######################################################################################
Project title: Comparing different graph similarites on gene co-expression network

Author: Zhuoyang Chen
######################################################################################

Since the processing of raw data is not included in the code because preprocessed data already existed from previous work, plus all these matrix objects are very large and it's not suitable for direct download. The missing part is only related to format the matrix but not change any values.

Raw data are available and downloaded directly from https://gdac.broadinstitute.org/

#######################################################################################
firehose.R: Not used in main analysis. Gives examples for data preprocessing, correlation calculation, adjacency matrix generation and similarity calculation.

DeltaCon.R: A revised version to fix some bugs and errors and make it suitable for direct sparseMatrix input, originated from github.com/bxshi

simFunctions.R: Source this. Self-defined functions for preprocessing, distance and similarity calculation algorithms.

generateAdj_pipeline.R: Generates signed-weighted, signed-unweighted, unsigned_weighted and unsigned_unweighted adjacency matrix for all 33 cancers expression data.

calculateSim_pipeline.R: Calculates Jaccard distance, Spectral distance by top 10 highest eigen values using Euclidean, DeltaCon similarity up to power of 5 for all four categories of adjacency matrix, for all combinations of cancer pairs, only keep nodes more than 1000 degrees. Size of graphs are summarized.

Hier_analysis.R: Generates results .png of heatmaps and correlation chats, including similarity between cancers by different methods for different adj_mats, hierachical structure of cancers, correlation between different adj_mats, correlation between different methods. Elapsed time and graph size visualization are shown.

#######################################################################################
Firstly use generateAdj_pipeline.R to generate the four adjacency matrices for futher similarity computations. Then use calculateSim_pipeline.R to calculate all 528 pair-wise graph similarity for four different matrix representations under three methods. Finally, use Hier_analysis.R generate the heatmap and correlation chart results.

Conclusion:
All three methods: Jaccard, spectral-based and DeltaCon similarity can reflex some extent of cancer similarity regarding cell/tissue origin and pathways. All shared medium similarity with DeltaCon on signed graphs the most different. When using Jaccard similarity, signs do not need to be considered and weight/unweight graphs yield similar results; when using Spectral similarity, neither signs or weights need to be considered and all generate similar results; for DeltaCon, weights are not considered and sign/unsign graphs produce different results. Jaccard similarity takes least time to compute and DeltaCon similarity takes longest, while spectral-based similarity depends on how many eigen values chosen to calculate. Elapsed time for all three methods scale with number of edges linear.
