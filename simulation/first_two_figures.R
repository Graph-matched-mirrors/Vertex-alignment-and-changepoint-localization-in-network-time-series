pacman::p_load(segmented, igraph, irlba ,RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)
library(ggrepel)
library(ggplot2)
library(dplyr)

## del is the set of shuffle percetage you want
doSim_London <- function(n=300, tmax=40, delta=0.1, p=0.4, q=0.9, tstar=20 , del = c(1))
{    
  glist <- NULL
  Xlist <- NULL
  Xt <- matrix(0,n,(tmax+1))
  
  for (t in 2:(tstar+1)) {
    tmp <- runif(n) < p
    Xt[,t] <- Xt[,t] + Xt[,t-1]
    Xt[tmp,t] <- Xt[tmp,t] + delta
  }
  
  for (t in (tstar+2):(tmax+1)) {
    tmp <- runif(n) < q
    Xt[,t] <- Xt[,t] + Xt[,t-1]
    Xt[tmp,t] <- Xt[tmp,t] + delta
  }
  Xt <- Xt[,-1]
  Xt <- Xt+0.1
  
  df <- tibble(time=1:tmax) %>%
    mutate(Xt = map(time, function(x) matrix(Xt[,x],n,1)  )) %>%
    mutate(g = map(Xt, ~rdpg.sample(.))) %>%
    mutate(avg_edges = map( g ,  ~sum(as.matrix(.))    )   ) %>%
    mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
  
  for(perc in del){
    df <- df %>%
      mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~shuffle_X(., perc)) )
  }
  
  df
}


true_distance_London=function(m,p){
  return( (m*(m-1)*p^2+m*p)*delta^2 )
}

##CMDS on the true distance matrix for model with change point                  
true_London_dMV = function(tt,t0,p,q){
  
  D=matrix(0,tt,tt)
  for (i in 1:t0) {
    for (j in (i):t0) {
      D[i,j]=sqrt(true_distance(abs(i-j),p))
    }
  }
  
  for (i in t0:tt) {
    for (j in (i):tt) {
      D[i,j]=sqrt(true_distance(abs(i-j),q))
    }
  }
  
  for (i in 1:(t0-1)) {
    for (j in (t0+1):tt) {
      m=t0-i
      n=j-t0
      D[i,j]=sqrt(true_distance(m,p)+true_distance(n,q)+2*m*n*p*q*delta^2)
    }
  }
  
  D=D+t(D)
  
  return(D)
  #D2=D^2
  
  #P=diag(rep(1,tt))-rep(1,tt)%*%t(rep(1,tt))/tt
  #B=(-1/(2))*P%*%D2%*%P
  
  #svdb= irlba(B,cdim)
  #df_true=as.data.frame(svdb$u[,1:cdim]%*%sqrt(diag(svdb$d[1:cdim])))
  
  #if(df_true$V1[1]>0){
  #  df_true$V1=-df_true$V1
  #}
  
  #df_true
}

true_shuffled_dMV_square_London = function(tt,t0,p,q){
  
  D=matrix(0,tt,tt)
  for (i in 1:t0) {
    for (j in (i):t0) {
      D[i,j]=sqrt(delta^2*(p^2*(i-j)^2 + (p-p^2)*(i+j)))
    }
  }
  
  for (i in t0:tt) {
    for (j in (i):tt) {
      D[i,j]=sqrt(delta^2*(q^2*(i-j)^2 + (q-q^2)*(i+j) +2*t0*(p-p^2-q+q^2) ))
    }
  }
  
  for (i in 1:(t0-1)) {
    for (j in (t0+1):tt) {
      m=t0-i
      n=j-t0
      D[i,j]=sqrt( delta^2*( ( (j-t0)*q+(t0-i)*p  )^2 + (p-p^2)*(i+t0) +(q-q^2)*(j-t0)    )  )
    }
  }
  D=D+t(D)
  D2=D^2
  diag(D2) = 0
  return(D2)
}

true_W1_square_London = function(tt,t0,p,q){
  D=matrix(0,tt,tt)
  for (i in 1:t0) {
    for (j in (i):t0) {
      D[i,j]=sqrt(delta^2*(p^2*(i-j)^2 ))
    }
  }
  
  for (i in t0:tt) {
    for (j in (i):tt) {
      D[i,j]=sqrt(delta^2*(q^2*(i-j)^2 ))
    }
  }
  
  for (i in 1:(t0-1)) {
    for (j in (t0+1):tt) {
      m=t0-i
      n=j-t0
      D[i,j]=sqrt( delta^2*( ( (j-t0)*q+(t0-i)*p  )^2    )  )
    }
  }
  D=D+t(D)
  D2=D^2
  return(D2)
}

true_avg_degree = function(x){
  
  if (x < tstar +1 ){
    y =  0.1+p*delta*x 
  } else{
    y = 0.1 + p*delta*tstar + (x - tstar)*q*delta
  }
  
  return( y^2*(n-1))
  
}

true_avg_degree = Vectorize(true_avg_degree)

psi_Z = function(x, p, q){
  if(x < 0.5){
    y = p*x
  } else {
    y = p*0.5+(x - 0.5)* q
  }
  return(y)
}
psi_Z = Vectorize(psi_Z)

psi_Z_center = function(x, p, q){
  c0 = 0.5*(p-q)*(0.25-1)-q/2
  if(x < 0.5){
    y = p*x +c0
  } else {
    y = p*0.5+(x - 0.5)* q +c0
  }
  return(y)
}
psi_Z_center = Vectorize(psi_Z_center)

#tmax = 10
#x = 1:tmax/tmax
#mean(psi_Z_center(x,p,q))
#sum(psi_Z_center(1:10000/10000,p,q))/10000
#psi_Z_center(1:1000/1000,p,q)



pacman::p_load("doParallel")
registerDoParallel(detectCores()-1)

## mds_nochange and df_true_nochange generation
n = 100
tmax <- m <- 30
p <- 0.3
q <- 0.9
delta <- (1-0.1)/tmax
tstar <- m/2


#True_dmv = true_London_dMV(tmax, tstar, p ,q )
#True_dmv_square = True_dmv^2
#True_shuffle_dmv_square = true_shuffled_dMV_square_London(tmax, tstar, p ,q )


True_D = true_London_dMV(tmax,tstar,p,q) 
MDS_True_D =doMDS(True_D, doplot = F)
True_shuffled_D = sqrt(true_shuffled_dMV_square_London(tmax,tstar,p,q))

MDS_True_shuffled_D = doMDS(True_shuffled_D , doplot = F)

set.seed(5)
df <- doSim_London(n, tmax, delta, p, q, tstar)

#df <- df %>%
#  mutate(shuffle_g = map(g, ~optimized_shuffle_graph(.))) %>%
#  mutate(Xhat_shuffle = map(shuffle_g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

del = c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,1)

for (alpha in del) {
  nm <- paste0("True_shuffle_dmv_", alpha)
  assign(nm,
         sqrt((1 - alpha) * True_D^2 + alpha * True_shuffled_D^2),
         envir = .GlobalEnv)
}

MDS_True_shuffle_0.2 <- doMDS(True_shuffle_dmv_0.2, doplot = F)

for(perc in del){
  df <- df %>%
    mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~ shuffle_X_optimized(., perc)) ) %>%
    mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat ,~shuffle_X_optimized(., perc)) )
}


D2_shuffle_0.1_hat <- getD(df$shuffle_Xhat_0.1)
MDS_shuffle_0.1_hat <- doMDS(D2_shuffle_0.1_hat, doplot = F)

D2_shuffle_0.2_hat <- getD(df$shuffle_Xhat_0.2)
MDS_shuffle_0.2_hat <- doMDS(D2_shuffle_0.2_hat, doplot = F)

D2_shuffle_1 <- getD(df$shuffle_Xhat_1)
mds_hat_shuffled_dMV = doMDS(D2_shuffle_1, doplot= F)

HatdMV_D <- getD(df$xhat) 
mds_hat_dMV <- doMDS(HatdMV_D,doplot = F)

D_london_W1 = getD_W1(df$xhat)
London_W1_mds = doMDS(D_london_W1, doplot = T)
true_London_W1_mds = doMDS(sqrt(true_W1_square_London(tmax, tstar , p , q )), doplot=F)
#par(mfrow= c(1,3))

#plot(1:tmax, mds_hat_dMV$mds[,1], col = 'blue',xlab = '', ylab = 'MDS1 ')
#points(1:tmax, MDS_True_D$mds[,1])
#(psi_Z(1:tmax/tmax,p,q)- mean(psi_Z(1:tmax/tmax,p,q)))*0.9
#MDS_True_D$mds[,1]

#True_shuffled_D = sqrt(true_shuffled_dMV_square_London(tmax,tstar,p,q))
#MDS_True_shuffled_D =doMDS(True_shuffled_D, doplot = F)


#Hat_shuffled_dMV_D <- getD(df$Xhat_shuffle) 
#mds_hat_shuffled_dMV <- doMDS(Hat_shuffled_dMV_D,doplot = F)

#plot(1:tmax, mds_hat_shuffled_dMV$mds[,1], col = 'blue' ,xlab = '', ylab = 'MDS1 ')
#points(1:tmax, MDS_True_shuffled_D$mds[,1])

#plot(1:tmax, sqrt(true_avg_degree(1:tmax)) , ylab = 'Avg degree' , xlab = '')
#points(1:tmax, sqrt(unlist(df$avg_edges)/n), col = 'blue')



library(tibble)
library(dplyr)
library(ggplot2)

df_plot <- bind_rows(
  # dMV
  tibble(x      = 1:tmax/tmax,
         y      = MDS_True_D$mds[,1],
         metric = "dMV",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = mds_hat_dMV$mds[,1],
         metric = "dMV",
         type   = "estimated"),
  tibble (x = 1:tmax/tmax,
          y = (psi_Z(1:tmax/tmax,p,q)- mean(psi_Z(1:tmax/tmax,p,q)))*0.9,
          metric = "dMV",
          type = 'psi z'
  ),
  
  # 20% shuffled-dMV
  tibble(x      = 1:tmax/tmax,
         y      = MDS_True_shuffle_0.2$mds[,1],
         metric = "20%-shuffled-dMV",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = MDS_shuffle_0.2_hat$mds[,1],
         metric = "20%-shuffled-dMV",
         type   = "estimated"),
  tibble (x = 1:tmax/tmax,
          y = (psi_Z(1:tmax/tmax,p,q)- mean(psi_Z(1:tmax/tmax,p,q)))*0.9,
          metric = "20%-shuffled-dMV",
          type = 'psi z'
  ),
  
  # 100% shuffled-dMV
  tibble(x      = 1:tmax/tmax,
         y      = MDS_True_shuffled_D$mds[,1],
         metric = "100%-shuffled-dMV",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = mds_hat_shuffled_dMV$mds[,1],
         metric = "100%-shuffled-dMV",
         type   = "estimated"),
  tibble (x = 1:tmax/tmax,
          y = (psi_Z(1:tmax/tmax,p,q)- mean(psi_Z(1:tmax/tmax,p,q)))*0.9,
          metric = "100%-shuffled-dMV",
          type = 'psi z'
  ),
  # W1 distance
  tibble(x      = 1:tmax/tmax,
         y      = (psi_Z(1:tmax/tmax, p, q) - mean(psi_Z(1:tmax/tmax, p, q)))*0.9,
         metric='W1',
         type   = "psi z"),
  tibble(x      = 1:tmax/tmax,
         y      = true_London_W1_mds$mds[,1],
         metric='W1',
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = London_W1_mds$mds[,1],
         metric='W1',
         type   = "estimated"
  ),
  # avg degree
  tibble(x      = 1:tmax/tmax,
         y      = sqrt(true_avg_degree(1:tmax))/sqrt(n-1)-mean(sqrt(true_avg_degree(1:tmax))/sqrt(n-1)),
         metric = "avg degree",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = sqrt(unlist(df$avg_edges) / (n) )/sqrt(n-1)- mean(sqrt(unlist(df$avg_edges) / (n) )/sqrt(n-1)),
         metric = "avg degree",
         type   = "estimated"),
  tibble (x = 1:tmax/tmax,
          y = (psi_Z(1:tmax/tmax,p,q)- mean(psi_Z(1:tmax/tmax,p,q)))*0.9,
          metric = "avg degree",
          type = 'psi z')
) %>%
  mutate(metric = factor(
    metric,
    levels = c("dMV" ,"20%-shuffled-dMV","100%-shuffled-dMV",  "W1","avg degree"),
    labels = c(
      "d[MV]",
      "20*\"%\"~shuffled~d[MV]",
      "100*\"%\"~shuffled~d[MV]",
      "W[1]",
      "sqrt(plain('avg degree')/(n-1))"
    )
  ))

ggplot(df_plot, aes(x = x, y = y, color = type)) +
  geom_point(
    data = df_plot %>% filter(type %in% c("true","estimated")),
    aes(x = x, y = y, color = type),
    size  = 2,
    alpha = 0.5
  ) +
  geom_line(
    data = df_plot %>% filter(type == "psi z"),
    aes(x = x, y = y),
    color = "red",
    size  = 1
  )+
  facet_wrap(~ metric,
             nrow     = 1,
             scales   = "free_y",
             labeller = label_parsed) +
  labs(x = "time",
       y = "graph dynamics representation") +
  scale_color_manual(values = c(true = "black",
                                estimated = "blue")) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12)
  )

ggsave("London-5-metrics.pdf", width=15, height=3, units="in")



ggplot(df_plot, aes(x = x, y = y, color = type)) +
  geom_point(
    data = df_plot %>% filter(type %in% c("true")) %>% filter(metric %in% c("d[MV]") ),
    aes(x = x, y = y, color = type),
    size  = 2,
    alpha = 0.5
  ) +
  labs(x = "time",
       y = "mirror") +
  scale_color_manual(values = c(true = "black",
                                estimated = "blue")) +
  geom_vline(xintercept = 0.5, color = "black", linetype = "dashed") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12)
  )

#ggsave("London-1-metric-pre.pdf", width=3, height=3, units="in")

ggplot(df_plot %>% filter(metric %in% c("d[MV]", "100*\"%\"~shuffled~d[MV]")), aes(x = x, y = y, color = type)) +
  geom_point(
    data = df_plot %>% filter(type %in% c("true", "estimated") & metric %in% c("d[MV]", "100*\"%\"~shuffled~d[MV]")),
    aes(x = x, y = y, color = type),
    size  = 2,
    alpha = 0.5
  ) +
  facet_wrap(~ metric,
             nrow     = 1,
             scales   = "free_y",
             labeller = label_parsed) +
  labs(x = "time",
       y = "mirror") +
  scale_color_manual(values = c(true = "black",
                                estimated = "blue")) +
  geom_vline(xintercept = 0.5, color = "black", linetype = "dashed") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12)
  )
ggsave("London-2-metric-pre.pdf", width=6, height=3, units="in")


## Atlanta
p = 0.05
q = 0.45
tmax = m = 30
tstar = 15

c=(.9-0.1)
num_state = 50
delta = c/(num_state-1)

True_dmv_square=true_Atlanta_dmv(p,q,num_state,m,tstar,delta) ## This is the analytically d^2MV result.
True_dmv = sqrt(True_dmv_square)
MDS_True_D_square =doMDS(True_dmv_square, doplot = F)
MDS_True_D =doMDS(sqrt(True_dmv_square), doplot = F)


yy = MDS_True_D_square$mds[,1]*num_state*(num_state-1)/(2*c^2*m)

if( yy[1] > 0){ yy = -yy} else {
  yy = yy
}


plot(1:tmax/tmax, psi_Z(1:tmax/tmax, p ,q) - mean(psi_Z(1:tmax/tmax, p ,q)) )
points(1:tmax/tmax , yy)

max(yy -( psi_Z(1:tmax/tmax, p ,q) - mean(psi_Z(1:tmax/tmax, p ,q)))) ## as N increase this will be smaller with order 1/N

#(psi_Z(1:tmax/tmax, p ,q) - mean(psi_Z(1:tmax/tmax, p ,q)))[1]
#yy[1] 


True_shuffle_dmv_square=true_shuffle_Atlanta_dmv(c, num_state , m )
True_shuffle_dmv = sqrt(True_shuffle_dmv_square)
MDS_True_shuffle_D = doMDS(True_shuffle_dmv, doplot = F)

set.seed(2)
n = 1000
xt=matrix(0,nrow = n, ncol = m+1)
initila_state_all_nodes=sample(seq(.1,.9,by=delta),n,replace = TRUE)

mean(initila_state_all_nodes)
the_var = c^2/12*( (num_state+1)/(num_state-1) )
print(paste('theretical variacne = ', the_var , 'sample variance = ',var(initila_state_all_nodes) )  )
sqrt(2*the_var) #is the entry of true 100% shuffling dMV
True_shuffle_dmv[1,2]

xt[,1]=initila_state_all_nodes
for (i in 1:n) {
  for (j in 2:(tstar)) {
    xt[i,j]=update_function(xt[i,j-1],p)
  }
  for (j in ((tstar+1):(m+1)) ) {
    xt[i,j]=update_function(xt[i,j-1],q)
  }
}
df <- tibble(time=1:m) %>%
  mutate(Xt = map(time, function(x) matrix(xt[,x],n,1)  ))  %>%
  mutate(g = map(Xt, ~rdpg.sample(.))) %>%
  mutate(avg_edges = map( g ,  ~sum(as.matrix(.))    )   ) %>%
  mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

D2 <- getD(df$xhat)
MDS_hat_dMV_square = doMDS(D2^2, doplot= F)

D_atlanta_W1 = getD_W1(df$xhat)
atlanta_W1_mds = doMDS(D_atlanta_W1, doplot = T)
#True_dMV square is harder to estimate
#summary(as.vector(abs(D2^2 - True_dmv_square)/True_dmv_square))
#summary(as.vector(abs(D2 - True_dmv)/True_dmv))



MDS_hat_dMV = doMDS(D2, doplot= F)


del = c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,1)

for (alpha in del) {
  nm <- paste0("True_shuffle_dmv_", alpha)
  assign(nm,
         sqrt((1 - alpha) * True_dmv_square + alpha * True_shuffle_dmv_square),
         envir = .GlobalEnv)
}

MDS_True_shuffle_0.2 <- doMDS(True_shuffle_dmv_0.2, doplot = F)

P=diag(rep(1,tmax))-rep(1,tmax)%*%t(rep(1,tmax))/tmax
B=(-1/(2))*P%*%True_dmv_square%*%P

sqrt(0.8+0.2*True_shuffle_dmv_square[1,2]/(2*eigen(B)$values[1]))
MDS_True_shuffle_0.2$mds[,1]/MDS_True_D$mds[,1]

for(perc in del){
  df <- df %>%
    mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~ shuffle_X_optimized(., perc)) ) %>%
    mutate(!!paste0("shuffle_Xhat_", perc) := map( xhat ,~shuffle_X_optimized(., perc)) )
}


D2_shuffle_0.1_hat <- getD(df$shuffle_Xhat_0.1)
MDS_shuffle_0.1_hat <- doMDS(D2_shuffle_0.1_hat, doplot = F)

D2_shuffle_0.2_hat <- getD(df$shuffle_Xhat_0.2)
MDS_shuffle_0.2_hat <- doMDS(D2_shuffle_0.2_hat, doplot = F)




D2_shuffle_1 <- getD(df$shuffle_Xhat_1)
MDS_hat_shuffle_dMV = doMDS(D2_shuffle_1, doplot= F)


df_plot <- bind_rows(
  # d^2MV
  tibble(x      = 1:tmax/tmax,
         y      = yy,
         metric = "d^2MV",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = MDS_hat_dMV_square$mds[,1]*num_state*(num_state-1)/(2*c^2*m),
         metric = "d^2MV",
         type   = "estimated"),
  tibble (x = 1:tmax/tmax,
          y = psi_Z(1:tmax/tmax,p,q)- mean(psi_Z(1:tmax/tmax,p,q)),
          metric = "d^2MV",
          type = 'psi z'
  ),
  #dMV
  tibble(x      = 1:tmax/tmax,
         y      = MDS_True_D$mds[,1],
         metric = "dMV",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = MDS_hat_dMV$mds[,1],
         metric = "dMV",
         type   = "estimated"),
  # 20% shuffled-dMV
  tibble(x      = 1:tmax/tmax,
         y      = MDS_True_shuffle_0.2$mds[,1],
         metric = "20%-shuffled-dMV",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = MDS_shuffle_0.2_hat$mds[,1],
         metric = "20%-shuffled-dMV",
         type   = "estimated"),
  # 100%-shuffled-dMV
  tibble(x      = 1:tmax/tmax,
         y      = MDS_True_shuffle_D$mds[,1],
         metric = "100%-shuffled-dMV",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = MDS_hat_shuffle_dMV$mds[,1],
         metric = "100%-shuffled-dMV",
         type   = "estimated"),
  tibble(x      = 1:tmax/tmax,
         y      = rep(0,tmax) ,
         metric = "W1",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = atlanta_W1_mds$mds[,1],
         metric = "W1",
         type   = "estimated"),
  # avg deg
  tibble(x      = 1:tmax/tmax,
         y      = rep(0.5^2*(n-1),tmax),
         metric = "avg degree",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = unlist(df$avg_edges) / n,
         metric = "avg degree",
         type   = "estimated")
) %>%
  mutate(metric = factor(
    metric,
    levels = c("d^2MV","dMV","20%-shuffled-dMV","100%-shuffled-dMV","W1", "avg degree"),
    labels = c(
      "frac(N*(N-1), 2*c[A]*m) ~ d[MV]^2",
      "d[MV]",
      "20*\"%\"~shuffled~d[MV]",
      "100*\"%\"~shuffled~d[MV]",
      "W[1]",
      "'avg degree'"
    )
  ))

ggplot() +
  # 1) your black/blue points for “true” & “estimated”
  geom_point(
    data = df_plot %>% filter(type %in% c("true","estimated")),
    aes(x = x, y = y, color = type),
    size  = 2,
    alpha = 0.5
  ) +
  scale_color_manual(values = c(true = "black", estimated = "blue")) +
  
  # 2) overlay psi_Z as a red line in the d^2MV panel
  geom_line(
    data = df_plot %>% filter(type == "psi z"),
    aes(x = x, y = y),
    color = "red",
    size  = 1
  ) +
  
  # 3) faceting and labels
  facet_wrap(~ metric,
             nrow     = 2,
             scales   = "free_y",
             labeller = label_parsed) +
  labs(x = "time",
       y = "graph dynamics representation") +
  
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12)
  )
ggsave("Atlanta-6-metrics_23.pdf", width=9, height=6, units="in")






ggplot(df_plot %>% filter(metric %in% c("d[MV]", "100*\"%\"~shuffled~d[MV]")), aes(x = x, y = y, color = type)) +
  geom_point(
    data = df_plot %>% filter(type %in% c("true", "estimated") & metric %in% c("d[MV]", "100*\"%\"~shuffled~d[MV]")),
    aes(x = x, y = y, color = type),
    size  = 2,
    alpha = 0.5
  ) +
  facet_wrap(~ metric,
             nrow     = 1,
             scales   = "free_y",
             labeller = label_parsed) +
  labs(x = "time",
       y = "mirror") +
  scale_color_manual(values = c(true = "black",
                                estimated = "blue")) +
  geom_vline(xintercept = 0.5, color = "black", linetype = "dashed") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12)
  )
ggsave("Atlanta-2-metric-pre.pdf", width=6, height=3, units="in")





#plot(1:tmax, MDS_shuffle_0.2_hat$mds[,1], col = 'blue' )
#points(1:tmax, MDS_True_shuffle_0.2$mds[,1], col = 'black' )


#D2_shuffle_0.2 <- getD(df$shuffle_Xt_0.2)
#summary(as.vector(abs(True_shuffle_dmv_0.2 - D2_shuffle_0.2)/True_shuffle_dmv_0.2))


#D2_shuffle_0.5 <- getD(df$shuffle_Xt_0.5)
#summary(as.vector(abs(True_shuffle_dmv_0.5 - D2_shuffle_0.5)/True_shuffle_dmv_0.5))

#D2_shuffle_1 <- getD(df$shuffle_Xt_1)
#summary(as.vector(abs(True_shuffle_dmv_1 - D2_shuffle_1)/True_shuffle_dmv_1))



#doMDS(D2_shuffle_0.01,doplot = T)

#D2_shuffle_0.05 <- getD(df$shuffle_Xhat_0.05)
#doMDS(D2_shuffle_0.05,doplot = T)

#D2_shuffle_0.1 <- getD(df$shuffle_Xhat_0.1)
#doMDS(D2_shuffle_0.1,doplot = T)

#D2_shuffle_0.2 <- getD(df$shuffle_Xhat_0.2)
#doMDS(D2_shuffle_0.2,doplot = T)

#D2 <- getD(df$xhat)
#MDS_hat_dMV_square = doMDS(D2^2, doplot= T)
#summary(as.vector(abs(True_dmv - D2)/True_dmv))





#mean(df$Xt[[4]])
#plot(1:tmax, rep(0.5^2*(n-1),tmax) , ylab = 'Avg degree' , xlab = '')
#points(1:tmax, unlist(df$avg_edges)/n, col = 'blue')


## W1 distance

set.seed(3)

n=1000
c=(.9-0.1)
p <- 0.05
q <- 0.45
num_state = 50
delta = c/(num_state-1)
xt=matrix(0,nrow = n, ncol = m+1)
initila_state_all_nodes=sample(seq(.1,.9,by=delta),n,replace = TRUE)

xt[,1]=initila_state_all_nodes
for (i in 1:n) {
  for (j in 2:(tstar)) {
    xt[i,j]=update_function(xt[i,j-1],p)
  }
  for (j in ((tstar+1):(m+1)) ) {
    xt[i,j]=update_function(xt[i,j-1],q)
  }
}
df <- tibble(time=1:m) %>%
  mutate(Xt = map(time, function(x) matrix(xt[,x],n,1)  ))  %>%
  mutate(g = map(Xt, ~rdpg.sample(.))) %>%
  mutate(avg_edges = map( g ,  ~sum(as.matrix(.))    )   ) %>%
  mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

D_atlanta_W1 = getD_W1(df$xhat)
atlanta_W1_mds = doMDS(D_atlanta_W1, doplot = T)

n = 100
tmax <- m <- 30
p <- 0.3
q <- 0.9
delta <- (1-0.1)/tmax
tstar <- m/2

df <- doSim_London(n, tmax, delta, p, q, tstar)
df <- df %>%
  mutate(shuffle_g = map(g, ~optimized_shuffle_graph(.))) %>%
  mutate(Xhat_shuffle = map(shuffle_g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

D_london_W1 = getD_W1(df$xhat)
London_W1_mds = doMDS(D_london_W1, doplot = F)
true_London_W1_mds = doMDS(sqrt(true_W1_square_London(tmax, tstar , p , q )), doplot=F)
#plot(1:tmax/tmax, (psi_Z(1:tmax/tmax, p, q) - mean(psi_Z(1:tmax/tmax, p, q)))*0.9 )
#points(1:tmax/tmax , true_London_W1_mds$mds[,1])
#points(1:tmax/tmax, London_W1_mds$mds[,1])





df_plot <- bind_rows(
  # London
  tibble(x      = 1:tmax/tmax,
         y      = (psi_Z(1:tmax/tmax, p, q) - mean(psi_Z(1:tmax/tmax, p, q)))*0.9,
         metric='W1',
         model = "london",
         type   = "psi_Z"),
  tibble(x      = 1:tmax/tmax,
         y      = true_London_W1_mds$mds[,1],
         metric='W1',
         model = "london",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = London_W1_mds$mds[,1],
         model = "london",
         metric='W1',
         type   = "estimated"),
  # Atlanta
  tibble(x      = 1:tmax/tmax,
         y      = rep(0,tmax) ,
         model = "Atlanta",
         metric = "W1",
         type   = "true"),
  tibble(x      = 1:tmax/tmax,
         y      = atlanta_W1_mds$mds[,1],
         model = "Atlanta",
         metric = "W1",
         type   = "estimated")
) 


# 1) ensure 'model' is a factor in the order you want
df_plot <- df_plot %>%
  mutate(
    model = factor(
      model,
      levels = c("london", "Atlanta"),
      labels = c("London", "Atlanta")
    )
  )

# 2) draw the two‐panel plot, overriding each strip label to "W1"
ggplot() +
  # black & blue points
  geom_point(
    data = df_plot %>% filter(type %in% c("true", "estimated")),
    aes(x = x, y = y, color = type),
    size  = 3,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(true = "black", estimated = "blue")) +
  
  # red psi_Z line (only for London, since only there you have type=="psi_Z")
  geom_line(
    data = df_plot %>% filter(type == "psi_Z"),
    aes(x = x, y = y),
    color = "red",
    size  = 1
  ) +
  
  # facet into two side‐by‐side panels, but force both strip titles to "W1"
  facet_wrap(
    ~ model,
    nrow     = 1,
    scales   = "free_y",
    labeller = labeller(model = c(London = "W1", Atlanta = "W1"))
  ) +
  
  # shared axis labels
  labs(
    x = "time",
    y = "graph dynamics representation"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12)
  )



ggsave("W1-both.pdf", width=8, height=4, units="in")
