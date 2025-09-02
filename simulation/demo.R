source('functions_updated.r')
source('London_model_simulations_updated.R')
pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)
library(ggrepel)

pacman::p_load("doParallel")
registerDoParallel(detectCores()-1)



#   n: number of nodes
#   m: number of time points
#   p: transition probability before the changepoint
#   q: transition probability after the changepoint
#   nmc: number of Monte Carlo simulations
#   mm: vector of different m values to test
#   d: vector of shuffling ratios
#   final_errors: optional, pre-existing list to store results
#   set_of_n: a vector of network sizes (number of nodes) to be tested 

n = 300
tmax <- m <- 40
p <- 0.3
q <- 0.4
delta <- (1-0.1)/tmax
tstar <- tmax / 2 
nmc = 500
mm <- c(12, 20, 30, 40)
d <- seq(0, 1, by = 0.05)
final_errors <- NULL
set_of_n <- c(100, 300, 800)

# Output:
# 1st-dim MDS overtime with true changepoint and estimated changepoint
# Original data, fitted linear model, estimated changepoint
analyze_network_changepoint(n, m, p, q) 

# Output: 
# 1st-dim MDS of dMV distance over time
# returns a list of L-inf errors for average edges, W1 distance, and dMV distance
Linf_errors(n, m, p, q) 

# Output:
# Given m values, plot showing mean MSE and conf interval for each error across 
# difference m; print final summary dataframe with corresponding values
simulate_network_changepoint(n, p, q, nmc, mm)

# Output:
# heatmap analyzing errors in shuffled network simulations: original vs. shuffled paired
analyze_shuffling_errors(m, tstar, delta, p, q, d, n, nmc, final_errors = NULL)

# Output:
# heatmap analyzing the impact of network size (number of nodes) and shuffling on change point detection in dynamic networks.
analyze_network_shuffling_wrt_n(m, tstar, delta, p, q, d, set_of_n, nmc)