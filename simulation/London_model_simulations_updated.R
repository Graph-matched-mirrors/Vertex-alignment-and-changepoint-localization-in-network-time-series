pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)
library(ggrepel)

pacman::p_load("doParallel")
registerDoParallel(detectCores()-1)

## mds_nochange and df_true_nochange generation

#   n: number of nodes
#   m: number of time points
#   p: transition probability before the changepoint
#   q: transition probability after the changepoint

analyze_network_changepoint <- function(n, m, p, q) {
  
  tmax = m 
  set.seed(2)

  delta <- (1 - 0.1) / tmax
  tstar <- tmax / 2
  
  df <- doSim_London(n, tmax, delta, p, q, tstar)
  
  # Calculate true and estimated distance matrices
  D_True_dMV <- true_London_dMV(tmax, tstar, p, q)
  D_dMV <- getD(df$xhat)
  
  # Multidimensional Scaling
  df.mds_dMV <- doMDS(D_dMV, doplot = FALSE)

  slope_cp <- find_slope_changepoint_with_plot(df.mds_dMV$mds[,1], doplot = TRUE)
  
  # L-infinity error
  ecp <- linf_error(df.mds_dMV$mds[,1])[2]
  
  # Plot MDS results
  plot(1:tmax, df.mds_dMV$mds[,1])
  abline(v = tstar)
  points(x = ecp, y = df.mds_dMV$mds[ecp, 1], pch = 2, col = "red")

  return(linf_error(sqrt(unlist(df$avg_edges))) * m)
}

#Linf erros for avg_edges, 1st-MDS-dim of W1/dMV distance
Linf_errors <- function(n, m, p, q) {
  tmax <- m
  delta <- (1-0.1)/tmax
  tstar <- tmax / 2 
  
  df <- doSim_London(n, tmax, delta, p, q, tstar)
  Dhat_W1 <- getD_W1(df$xhat)
  D_dMV <- getD(df$xhat)
  
  df.mds_W1 <- doMDS(Dhat_W1, doplot = FALSE)
  df.mds_dMV <- doMDS(D_dMV, doplot = FALSE)
  
  error_avg_edges <- linf_error(sqrt(unlist(df$avg_edges)))*m
  error_W1 <- linf_error(df.mds_W1$mds[,1])*m
  error_dMV <- linf_error(df.mds_dMV$mds[,1])*m
  
  plot(1:tmax, df.mds_dMV$mds[,1])
  
  return(list(
    error_avg_edges = error_avg_edges,
    error_W1 = error_W1,
    error_dMV = error_dMV
  ))
}

#   n: number of nodes
#   p: transition probability before changepoint
#   q: transition probability after changepoint
#   nmc: number of Monte Carlo simulations
#   mm: vector of different m values to test

simulate_network_changepoint <- function(n, p, q, nmc, mm) {

  set.seed(2)
  
  # Initialize a list to store results
  results_list <- list()
  
  # Error types
  error_types <- c("error_avg_degree", "error_W1", "error_W2", "error_dMV")
  
  for (i in 1:length(mm)) {
    tmax <- m <- mm[i]
    tstar <- tmax / 2
    delta <- (1-0.1)/tmax
    
    matrix_name <- paste0("error_matrix_m_", m, "_n_", n)
    assign(matrix_name, matrix(nrow = 4, ncol = nmc))
    
    # Monte Carlo simulations
    for (mc in 1:nmc) {
      df <- doSim_London(n, tmax, delta, p, q, tstar)
      
      # Compute square root of avg_edges
      sqrt_edges <- sqrt(unlist(df$avg_edges))
      
      # Compute distance matrices
      Dhat_W1 <- getD_W1(df$xhat)
      Dhat_W2 <- getD_W2(df$xhat)
      D_dMV  <- getD(df$xhat)
      
      # Perform MDS
      df.mds_W1 <- doMDS(Dhat_W1, doplot = FALSE)
      df.mds_W2 <- doMDS(Dhat_W2, doplot = FALSE)
      df.mds_dMV <- doMDS(D_dMV, doplot = FALSE)
      
      # Calculate errors
      error_avg_degree <- linf_error(sqrt_edges)
      error_W1 <- linf_error(df.mds_W1$mds[,1])
      error_W2 <- linf_error(df.mds_W2$mds[,1])
      error_dMV <- linf_error(df.mds_dMV$mds[,1])
      
      # Update error matrix
      temp_matrix <- get(matrix_name)
      temp_matrix[, mc] <- c(error_avg_degree, error_W1, error_W2, error_dMV)
      assign(matrix_name, temp_matrix)
      
      if (mc %% 100 == 0) {
        print(c(i, mc))
      }
    }
    
    # Compute mean squared error and confidence interval
    final_matrix <- get(matrix_name)
    mse_mean <- apply(abs(final_matrix)^2, 1, mean)
    mse_sd <- apply(abs(final_matrix)^2, 1, sd)
    conf_interval <- 1.96 * mse_sd / sqrt(nmc)
    
    temp_df <- data.frame(
      m = rep(m, 4),
      error_type = error_types,
      Lower_CI = mse_mean - conf_interval,
      Mean_MSE = mse_mean,
      Upper_CI = mse_mean + conf_interval
    )
    
    results_list[[i]] <- temp_df
  }
  
  summary_df <- do.call(rbind, results_list)
  
  plottt1 <- ggplot(summary_df, aes(x = m, y = Mean_MSE, color = error_type)) +
    geom_line() +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.5) +
    scale_x_continuous(breaks = mm) +
    labs(y = 'MSE', x = 'm', title =  paste('n =', n, 
                                            'nmc =', nmc, 
                                            'p =', p, 
                                            'q =', q)) +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15, face = "bold")
    )
  
  # Print the plot
  print(plottt1)
  
  # Display the final summary table
  print(summary_df)
  
  return(summary_df)
}


# Function to analyze errors in shuffled network simulations
#   m: number of time points
#   tstar: changepoint time
#   delta: step size for latent positions
#   p: transition probability before changepoint
#   q: vector of transition probabilities after changepoint
#   d: vector of shuffling ratios
#   n: number of nodes
#   nmc: number of Monte Carlo simulations
#   final_errors: optional, pre-existing list to store results
analyze_shuffling_errors <- function(m, tstar, delta, p, q, d, n, nmc, final_errors = NULL) {
  
  # Initialize if not provided
  if (is.null(final_errors)) {
    final_errors <- list()
  }
  
  for(i in 1:length(q)){
    temp_errors <- paired_error_in_shuffling(nmc = nmc, n = n, p = p, q = q[i], m = m, delta = delta, tstar = tstar, del = d)
    
    # Calculate MSE and standard deviations
    final_errors[[paste0("row_mse_", q[i])]] <- temp_errors[,1]
    final_errors[[paste0("row_sds_", q[i])]] <- temp_errors[,2]
    
    print(paste(q[i],"is done!"))
    print(final_errors[[paste0("row_mse_", q[i])]])
    print(final_errors[[paste0("row_sds_", q[i])]])
  }
  
  mse_keys <- grep("^row_mse_", names(final_errors), value = TRUE)
  mse_values <- final_errors[mse_keys]
  mse_matrix <- do.call(cbind, mse_values)
  colnames(mse_matrix) <- q 
  rownames(mse_matrix) <- c(0,d)
  
  sd_keys <- grep("^row_sds_", names(final_errors), value = TRUE)
  sd_values <- final_errors[sd_keys]
  sd_matrix <- do.call(cbind, sd_values)
  colnames(sd_matrix) <- q 
  rownames(sd_matrix) <- c(0,d)
  
  # Heatmap of MSE results
  heatmap(mse_matrix, Rowv = NA, Colv = NA,
          main = paste("n =",n, ", p =", p, ", nmc =", nmc),
          xlab = "q", ylab = "shuffling ratio", scale = "none") 
  
  # Save results to CSV
  write.csv(final_errors, "~/final_errors.csv", row.names = FALSE)
  
  return(list(final_errors = final_errors, mse_matrix = mse_matrix, sd_matrix = sd_matrix))
}

#   m: number of time points
#   tstar: changepoint time
#   delta: step size for latent positions
#   p: transition probability before changepoint
#   q: transition probability after changepoint
#   d: vector of shuffling ratios
#   set_of_n: vector of network sizes to test
#   nmc: number of Monte Carlo simulations

# analyze the impact of network size (number of nodes) and shuffling on change point detection in dynamic networks.
analyze_network_shuffling_wrt_n <- function(m = 50, tstar = 25, delta = 0.1, p = 0.4, q = 0.25, 
                                      d = seq(0.05, 1, by = 0.05), set_of_n = c(100,300,800), 
                                      nmc = 100) {
  
  final_errors_diff_n <- list()
  
  for(i in 1:length(set_of_n)){
    # Perform paired error in shuffling simulation
    temp_errors <- paired_error_in_shuffling(nmc = nmc, n = set_of_n[i], p = p, q = q, 
                                             m = m, delta = delta, tstar = tstar, del = d)
    
    # Store MSE and standard deviations
    final_errors_diff_n[[paste0("row_mse_", set_of_n[i])]] <- temp_errors[,1]
    final_errors_diff_n[[paste0("row_sds_", set_of_n[i])]] <- temp_errors[,2]
    
    print(paste(set_of_n[i],"is done!"))
    print(final_errors_diff_n[[paste0("row_mse_", set_of_n[i])]])
    print(final_errors_diff_n[[paste0("row_sds_", set_of_n[i])]])
  }
  
  mse_keys <- grep("^row_mse_", names(final_errors_diff_n), value = TRUE)
  mse_values <- final_errors_diff_n[mse_keys]
  mse_matrix <- do.call(cbind, mse_values)
  colnames(mse_matrix) <- set_of_n
  rownames(mse_matrix) <- c(0,d)
  
  # Extract and format standard deviation results
  sd_keys <- grep("^row_sds_", names(final_errors_diff_n), value = TRUE)
  sd_values <- final_errors_diff_n[sd_keys]
  sd_matrix <- do.call(cbind, sd_values)
  colnames(sd_matrix) <- set_of_n
  rownames(sd_matrix) <- c(0,d)
  
  # Create heatmap of MSE results
  heatmap(mse_matrix, Rowv = NA, Colv = NA,
          main = paste("q =",q, ", p =", p, ", nmc =", nmc),
          xlab = "n", ylab = "shuffling ratio", scale = "none")
  
  write.csv(final_errors_diff_n, "~/final_errors_diff_n.csv", row.names = FALSE)
  
  return(list(final_errors_diff_n = final_errors_diff_n, 
              mse_matrix = mse_matrix, 
              sd_matrix = sd_matrix))
}







