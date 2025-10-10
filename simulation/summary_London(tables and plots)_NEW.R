## This is the file used to summarized the results generated from MCSim_London.R


## for local load 
setwd('.../Vertex-alignment-and-changepoint-localization-in-network-time-series/simulation/Simulation_results')
result_name = 'out_dd_Londonn200_m20_p0.4_q0.3_max_iter100_20250328_0038.RData' #This is left pannel Table 1 
result_name = 'out_dd_Londonn200_m20_p0.3_q0.4_max_iter100_20250328_0114.RData'

# Extract simulation parameters from result_name
pattern <- "n(\\d+)_m(\\d+)_p([0-9.]+)_q([0-9.]+)_max_iter(\\d+)"
matches <- regmatches(result_name, regexec(pattern, result_name))[[1]]
matches

n         <- as.numeric(matches[2])
m         <- as.numeric(matches[3])
p         <- as.numeric(matches[4])
q         <- as.numeric(matches[5])
max_iter  <- as.numeric(matches[6])



load(result_name)



#out_dd

# Create a summary data frame from out_dd (ignoring the example_df)
summary_df <- do.call(rbind, lapply(seq_along(out_dd), function(i) {
  out <- out_dd[[i]]
  data.frame(
    mc = i,
    true_mds1_iso1 = out$tmp1[1],
    true_mds2_iso1 = out$tmp1[2],
    true_mds3_iso1 = out$tmp1[3],
    shuffled_mds1_iso1 = out$tmp2[1],
    shuffled_mds2_iso1 = out$tmp2[2],
    shuffled_mds3_iso1 = out$tmp2[3],
    gm_alltoone_mds1_iso1 = out$tmp3[1],
    gm_alltoone_mds2_iso1 = out$tmp3[2],
    gm_alltoone_mds3_iso1 = out$tmp3[3],
    gm_pairwise_mds1_iso1 = out$tmp4[1],
    gm_pairwise_mds2_iso1 = out$tmp4[2],
    gm_pairwise_mds3_iso1 = out$tmp4[3],
    tmp_W1 = out$tmp_W1,
    tmp_avg_edges = out$tmp_avg_edges
  )
}))



num_df <- summary_df[, !names(summary_df) %in% "mc"]

squared_df <- num_df^2

nmc <- nrow(summary_df)

# Compute the mean and standard deviation for each metric (each column)
col_means <- colMeans(squared_df)
col_sds   <- apply(squared_df, 2, sd)

lower_bound <- col_means - 1.96 * col_sds / sqrt(nmc)
upper_bound <- col_means + 1.96 * col_sds / sqrt(nmc)

ci_df <- data.frame(lower_bound = lower_bound, mean = col_means, upper_bound = upper_bound)
ci_df$metric <- rownames(ci_df)


ci_df$iso <- ifelse(grepl("mds1", ci_df$metric), "iso1",
                    ifelse(grepl("mds2", ci_df$metric), "iso2",
                           ifelse(grepl("mds3", ci_df$metric), "iso3", NA)))


ci_df$category <- ifelse(grepl("true", ci_df$metric), "true alignment",
             ifelse(grepl("shuffled", ci_df$metric), "shuffled",
                ifelse(grepl("alltoone", ci_df$metric), "alltoone",
                     ifelse(grepl("pairwise", ci_df$metric), "pairwise", NA))))

ci_df$iso[grepl("tmp_W1|tmp_avg_edges", ci_df$metric)] <- "iso1"

ci_df$category[grepl("tmp_W1", ci_df$metric)] <- "W1"
ci_df$category[grepl("tmp_avg_edges", ci_df$metric)] <- "avg_edges"

ci_df[ci_df$iso == 'iso1',] # This is Table 1 in the paper, left or right depends on which file you use.


support = seq(2/m-0.5, (m-1)/m-0.5,by=1/m)
chance_level = sum(support^2/length(support))
chance_level


