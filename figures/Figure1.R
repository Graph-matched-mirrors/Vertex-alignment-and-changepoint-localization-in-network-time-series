source('~/Desktop/Research /PhDresearch/London model with GM/Vertex-alignment-and-changepoint-localization-in-network-time-series/simulation/utility_functions.r')
# --- 2. Generate Data for London Model ---
set.seed(2)
n = 200
tmax <- m <- 30
p_lon <- 0.3
q_lon <- 0.9
delta <- (1-0.1)/tmax
tstar <- m/2

# True dMV
True_D_lon = true_London_dMV(tmax,tstar,p_lon,q_lon,delta)
MDS_True_D_lon = doMDS(True_D_lon, doplot = F)

# Simulation
df_lon <- doSim_London(n, tmax, delta, p_lon, q_lon, tstar)

# Estimated dMV
HatdMV_D_lon <- getD(df_lon$xhat)
mds_hat_dMV_lon <- doMDS(HatdMV_D_lon,doplot = F)

# 100% Shuffled Estimated dMV
df_lon <- df_lon %>%
  mutate(shuffle_Xhat_1 = map( xhat ,~shuffle_X_optimized(., 1)) )
D2_shuffle_1_lon <- getD(df_lon$shuffle_Xhat_1)
mds_hat_shuffled_dMV_lon = doMDS(D2_shuffle_1_lon, doplot= F)

# Theoretical Line
psi_z_lon = (psi_Z(1:tmax/tmax,p_lon,q_lon)- mean(psi_Z(1:tmax/tmax,p_lon,q_lon)))*0.9

# --- 3. Generate Data for Atlanta Model ---
set.seed(12)
p_atl = 0.1
q_atl = 0.45
tmax = m = 30
tstar = 15
n = 1000
c=(.9-0.1)
num_state = 50
delta = delta_atl = c/(num_state-1)

# True dMV (Analytically derived d^2MV, then sqrt)
True_dmv_square_atl=true_Atlanta_dmv(p_atl,q_atl,num_state,m,tstar,delta_atl)
True_dmv_atl = sqrt(True_dmv_square_atl)
MDS_True_D_atl =doMDS(True_dmv_atl, doplot = F)

# Simulation
xt=matrix(0,nrow = n, ncol = m+1)
initila_state_all_nodes=sample(seq(.1,.9,by=delta_atl),n,replace = TRUE)
xt[,1]=initila_state_all_nodes
for (i in 1:n) {
  for (j in 2:(tstar)) {
    xt[i,j]=update_function(xt[i,j-1],p_atl,delta)
  }
  for (j in ((tstar+1):(m+1)) ) {
    xt[i,j]=update_function(xt[i,j-1],q_atl,delta)
  }
}

df_atl <- tibble(time=1:m) %>%
  mutate(Xt = map(time, function(x) matrix(xt[,x],n,1) ))%>%
  mutate(g = map(Xt, ~rdpg.sample(.))) %>%
  mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

# Estimated dMV
D2_atl <- getD(df_atl$xhat)
MDS_hat_dMV_atl = doMDS(D2_atl, doplot= F)

# 100% Shuffled Estimated dMV
df_atl <- df_atl %>%
  mutate(shuffle_Xhat_1 = map( xhat ,~shuffle_X_optimized(., 1)) )
D2_shuffle_1_atl <- getD(df_atl$shuffle_Xhat_1)
MDS_hat_shuffle_dMV_atl = doMDS(D2_shuffle_1_atl, doplot= F)

# --- 4. Create Data Frame for Plotting ---

create_plot_df <- function(mds_true, mds_est, mds_shuff, model_name) {
  bind_rows(
    # True Panel
    tibble(x = 1:tmax/tmax, y = mds_true, metric = "True", type = "point", model = model_name),
    # Estimated Panel
    tibble(x = 1:tmax/tmax, y = mds_est, metric = "Estimated (No Shuffling)", type = "point", model = model_name),
    # Shuffled Panel
    tibble(x = 1:tmax/tmax, y = mds_shuff, metric = "Estimated (100% Shuffling)", type = "point", model = model_name),
  )
}

df_london_plot <- create_plot_df(MDS_True_D_lon$mds[,1], mds_hat_dMV_lon$mds[,1], mds_hat_shuffled_dMV_lon$mds[,1],"London")
df_atlanta_plot <- create_plot_df(MDS_True_D_atl$mds[,1], MDS_hat_dMV_atl$mds[,1], MDS_hat_shuffle_dMV_atl$mds[,1], "Atlanta")


df_final <- bind_rows(df_london_plot, df_atlanta_plot) %>%
  mutate(model = factor(model, levels = c("London", "Atlanta"))) %>%
  mutate(metric = factor(metric, levels = c("True", "Estimated (No Shuffling)", "Estimated (100% Shuffling)"))) %>%
  filter(type == "point") 

# --- 5. Generate the Plot ---

gg <- ggplot(df_final, aes(x = x, y = y)) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", alpha = 1.2) +
  
  geom_line(data = df_final %>% filter(metric == "True"),
            aes(color = metric), size = 1.2, alpha = 0.8) +
  
  geom_point(aes(color = metric), size = 2.3, alpha = 0.6) +
  
  facet_grid(model ~ metric, scales = "free_y", switch = "y") +
  
  scale_color_manual(values = c("True" = "red", 
                                "Estimated (No Shuffling)" = "blue", 
                                "Estimated (100% Shuffling)" = "#00ba38")) + 
  
  # 6. Theme and Labels
  theme_bw() +
  labs(x = "time", y = expression(d[MV] ~ mirror)) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 15),
    strip.placement = "outside" 
  )

# Display the plot
print(gg)



ggsave(
  filename = "AL2.pdf",
  plot = gg,
  width = 15, 
  height = 10, 
  units = "in",
  device = cairo_pdf # Matches the "Use cairo_pdf device" option
)
