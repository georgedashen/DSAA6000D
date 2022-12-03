# DSAA6000D
This is for the final project for DSAA6000D

Project title: Comparing different graph similarites on gene co-expression network
Author: Zhuoyang Chen

Raw data are available and downloaded directly from https://gdac.broadinstitute.org/

firehose.R: Not used in main analysis. Gives examples for data preprocessing, correlation calculation, adjacency matrix generation and similarity calculation.

DeltaCon.R: A revised version to fix some bugs and errors and make it suitable for direct sparseMatrix input, originated from github.com/bxshi

simFunctions.R: Source this. Self-defined functions for preprocessing, distance and similarity calculation algorithms.

generateAdj_pipeline.R: Generates signed-weighted, signed-unweighted, unsigned_weighted and unsigned_unweighted adjacency matrix for all 33 cancers expression data.

calculateSim_pipeline.R: Calculates Jaccard distance, Spectral distance by top 10 highest eigen values using Euclidean, DeltaCon similarity up to power of 5 for all four categories of adjacency matrix, for all combinations of cancer pairs, only keep nodes more than 1000 degrees. Size of graphs are summarized.

Hier_analysis.R: Generates results .png of heatmaps and correlation chats, including similarity between cancers by different methods for different adj_mats, hierachical structure of cancers, correlation between different adj_mats, correlation between different methods. Elapsed time and graph size visualization are shown.
