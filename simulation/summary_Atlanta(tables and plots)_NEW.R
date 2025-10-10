#This is the file used to summarized the results generated from MCSim_Atlanta.R

setwd('.../Vertex-alignment-and-changepoint-localization-in-network-time-series/simulation/Simulation_results')



result_name = "out_dd_n500_m20_p0.4_q0.2_nmc300_num_state50_max_iter100_20250908_1644.RData"
#result_name = "out_dd_n800_m20_p0.4_q0.2_nmc300_num_state50_max_iter100_20250912_1716.RData"


# Extract simulation parameters from result_name
pattern <- "n(\\d+)_m(\\d+)_p([0-9.]+)_q([0-9.]+)_nmc(\\d+)_num_state(\\d+)_max_iter(\\d+)"
matches <- regmatches(result_name, regexec(pattern, result_name))[[1]]
n         <- as.numeric(matches[2])
m         <- as.numeric(matches[3])
p         <- as.numeric(matches[4])
q         <- as.numeric(matches[5])
nmc       <- as.numeric(matches[6])
num_state <- as.numeric(matches[7])
max_iter  <- as.numeric(matches[8])

load(result_name)

#df = out_dd[[100]]$example_df
res_matrix <- do.call(rbind, lapply(out_dd, function(element) {
  # Select only the list elements that are NOT 'example_df'
  numeric_elements <- element[names(element) != "example_df"]
  
  # Now unlist, which will only operate on the numeric vectors
  row <- unlist(numeric_elements)
  
  return(row)
}))


# Convert to a data frame for easier use
res_matrix  <- as.data.frame(res_matrix)

# Assign the column names to the final data frame
colnames(res_matrix ) <- c(
  "true1", "true_iso_d1", "true_iso_d4", "true_iso_d8",
  "shuffle1", "shuffle_iso_d1", "shuffle_iso_d4", "shuffle_iso_d8",
  "gm_alltoone1", "gm_alltoone_iso_d1", "gm_alltoone_iso_d4", "gm_alltoone_iso_d8",
  "gm_pairwise1", "gm_pairwise_iso_d1", "gm_pairwise_iso_d4", "gm_pairwise_iso_d8",
  "W1", "Avg_degree"
)


# Number of simulations used in the analysis
nmc <- length(out_dd)

# Calculate summary statistics for each column in res_matrix
col_means <- apply(abs(res_matrix)^2, 2, mean)
col_sds   <- apply(abs(res_matrix)^2, 2, sd)

# Organize the summary into a data frame
summary_df <- data.frame(
  metric = names(col_means),
  mean   = col_means,
  sd     = col_sds
)

# Compute the lower and upper bounds of the 95% Confidence Interval
summary_df$lower_bound <- summary_df$mean - 1.96 * summary_df$sd / sqrt(nmc)
summary_df$upper_bound <- summary_df$mean + 1.96 * summary_df$sd / sqrt(nmc)


library(dplyr)

# Create two new variables: 'type' and 'iso' from the metric names
summary_df <- summary_df %>%
  mutate(type = case_when(
           grepl("^true", metric)       ~ "dMV true alignment",
           grepl("^shuffle", metric)      ~ "dMV shuffled",
           grepl("^gm_alltoone", metric)  ~ "GM all to one",
           grepl("^gm_pairwise", metric)  ~ "GM consecutive"
         ),
         iso = case_when(
           grepl("iso_d1", metric) ~ "iso_d1+D",
           grepl("iso_d4", metric) ~ "iso_d4+D",
           grepl("iso_d8", metric) ~ "iso_d8+D",
           TRUE                  ~ "MDS1+D^2"   # Default case: plain "1"
         ))


summary_df_filtered = summary_df[summary_df$iso == 'MDS1+D^2',]
formatted_df <- summary_df_filtered %>%
  mutate(across(where(is.numeric), ~ formatC(., digits = 4, format = "g", flag = "#")))

formatted_df[,c(1,4,2,5)] ## This is table 2

