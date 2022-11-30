# DSAA6000D
This is for the final project for DSAA6000D

Project title: Calculating Graph Similarity for Signed Weighted Graphs from Gene Co-expression Analysis

Author: Zhuoyang Chen

firehose.R: Gives examples for data preprocessing, correlation calculation, adjacency matrix generation and similarity calculation

simFunctions.R: Self-defined functions for preprocessing, distance and similarity calculation algorithms

generateAdj_pipeline.R: Generate signed-weighted, signed-unweighted, unsigned_weighted and unsigned_unweighted adjacency matrix for all 33 cancers expression data.

calculateSim_pipeline.R: Calculate Jaccard distance, Spectral distance up to 10 eigen values using Euclidean, DeltaCon similarity up to power of 5 for all four categories of adjacency matrix, for all combinations of cancer pairs
