## processing procedure(DO NOT RUN, data is not available)

#datanum = 1
#dat = data.table::fread(paste0("~/Dropbox/SofA/NRL/MillingDataExamples/MillTimeSeries_",datanum,".csv"), header=T) |> as_tibble()
#tind = as.numeric(names(dat)) |> round(2)

# build Xt
# 2-101 for x, 102-201 for y, 202-301 for x-vel, 302-401 for y-vel
#(n = 100)
#names(dat) = as.character(tind)
#dat = dat |> mutate(i=rep(1:n, times=4)) |>
#  mutate(type = rep(c("x","y","xv","yv"), each=n))
# rearrange dat as {i,x,y,xv,yv,time}, `time` is the name of the columns
#df = dat |> pivot_longer(cols=-c(i,type), names_to="time", values_to="val") |>
#  mutate(time=as.numeric(time)) |>
#  pivot_wider(names_from=type, values_from=val) |>
#  dplyr::select(i,x,y,xv,yv,time) |> arrange(i,time)

# subset t
#use_t = seq(780, 820, by=0.1)
#tstart = use_t[1]
#tstep = use_t[2] - use_t[1]
#tstar = c(800, 810)

#df2 = df |> filter(time %in% use_t)

# build glist
# from t=1...,t*_1
#tind2 = use_t
#(ts = length(tind2))
#if (FALSE) {
#  Xt <- glist <- NULL
#  set.seed(12345)
#  for (j in 1:ts) {
#    Xt[[j]] <- df2 %>% filter(time==tind2[j]) %>% dplyr::select(x,y) %>% as.matrix()
#    glist[[j]] <- rdpg.sample(Xt[[j]])
#  }
#  save(glist, paste0(file="~/Dropbox/SofA/glist-NRL-",tstart,"-by",round(tstep,2),".RData"))
#} else {
#  load(paste0("~/Dropbox/SofA/glist-NRL-",tstart,"-by",round(tstep,2),".RData"))
#}


###################(RUN below)


## Functions note getDW1 and getDW2 functions are different from other R files functions, because in this case we need to calculate for 2 dimension. 

pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)
library(ggrepel)
library(transport)



add_noise_g = function(g, shuffled.v=0) {
  n = vcount(g)
  vid = 1:n
  cand = sample(n, round(n*shuffled.v))
  vid2 = replace(vid, cand, sample(vid[cand]))
  g = igraph::permute(g, vid2)
  return(g)
}

# this is the same as yours except `center_graph`, which was advised by the author of the package!
graph_matching_YP <- function(g1, g2) {
  gm = gm(center_graph(g1), center_graph(g2), start = "rds")
  perm = diag(length(gm[,2]))[gm[,2],]
  new_graph = as.matrix(perm%*%as.matrix(g2)%*% t(perm))
  #print(c("Pre_Frob=",sqrt(sum((G1-G2)^2)),"Post_Frob=", sqrt(sum((G1-new_graph)^2)) ) )
  return( graph_from_adjacency_matrix(new_graph, mode = 'undirected') )
}

##ASE for a network A with embedding dimension d
full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  require(irlba)
  
  # doptr
  if (doptr) {
    g <- ptr(A)
    A <- g[]
  } else {
    A <- A[]
  }
  
  # diagaug
  if (diagaug) {
    diag(A) <- rowSums(A) / (nrow(A)-1)
  }
  
  A.svd <- svds(A,k=d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  
  if (!isSymmetric(A)) {
    Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d))
  }
  
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), Xhat.R=Xhat.R))
}



graph_mathing <- function(stand,mess,max_it){
  G1=as.matrix( stand )
  G2=as.matrix( mess )
  
  gm=gm(A=stand,B=mess,start = "rds", max_iter = max_it)
  perm=diag(length(gm[,2]))[gm[,2],]
  
  new_graph=as.matrix(perm%*%as.matrix( mess )%*% t(perm))
  
  #print(c("Pre_Frob=",sqrt(sum((G1-G2)^2)),"Post_Frob=", sqrt(sum((G1-new_graph)^2)) ) )
  
  return( graph_from_adjacency_matrix(new_graph,mode = 'undirected') )
}


#graph_mathing(df$g[[2]],df$g[[3]])

## Procruste/i.e gain W in the estimated d_MV distance. Note in our case latent positions are 1d so this is trivial with w=1or -1
procrustes2 <- function(X, Y) {
  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  newX <- X %*% W
  return(list(newX = newX, error = norm(newX-Y, type="F"), W = W))
}



## Get distance matrix
getD <- function(Xlist, k=0, etype="proc") {
  m <- length(Xlist)
  n <- nrow(Xlist[[1]])
  
  if (k==0) {
    ind <- 1:n
  } else {
    ind <- which(Yhat==k)
  }
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    #cat("i = ", i, ", j = ", j, "\n")
    
    if (etype == "proc") {
      Xhati <- Xlist[[i]][ind,] # 32277 x Khat
      Xhatj <- Xlist[[j]][ind,]
      proc <- procrustes2(as.matrix(Xhati), as.matrix(Xhatj))
      Xhati <- Xhati %*% proc$W
    } else {
      Xhati <- Xlist[[i]][ind,] # 32277 x Khat
      Xhatj <- Xlist[[j]][ind,]
    }
    
    D <- norm(Xhati - Xhatj, type="2")^2/n
    tibble(i=i, j=j, D=D)
  }
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- (D2 + t(D2)) / 1
  #as.dist(D2)
  D2 <- sqrt(D2)
  D2
}

## Apply CMDS on distance matrix 
doMDS <- function(D, doplot=TRUE)
{
  tmax <- m <- nrow(D)
  mds <- cmdscale(D, m-1)
  
  df.mds <- tibble(
    ind = 1:tmax,
    time = sprintf("%2d", 1:tmax),
    x = mds[, 1],  # Always take the first column
    
    # Check if additional columns exist before assigning
    y = if (ncol(mds) >= 2) mds[, 2] else NA, 
    z = if (ncol(mds) >= 3) mds[, 3] else NA,
    w = if (ncol(mds) >= 4) mds[, 4] else NA
  )
  
  if (doplot) {
    
    library('gridExtra')
    library('ggplot2')
    library('ggrepel')
    library('grid')  # Needed for textGrob()
    
    # Get the matrix name for the big title
    matrix_name <- deparse(substitute(D))
    big_title <- textGrob(paste("CMDS on", matrix_name), gp=gpar(fontsize=18, fontface="bold"))
    
    # Convert the result of apply(mds, 2, sd) into a data frame for ggplot
    sd_data <- data.frame(
      dimension = 1:ncol(mds),  # Assumes mds is a matrix or data frame
      column_sd = apply(mds, 2, sd)
    )
    
    # First plot: standard deviation by dimension (converted to ggplot)
    p0 <- ggplot(sd_data, aes(x = dimension, y = column_sd)) +
      geom_point() +
      geom_line() +
      labs(x = "dimension", y = "column stdev") +
      theme_minimal()
    
    # Second plot: MDS1 vs time
    p1 <- df.mds %>%
      ggplot(aes(x=ind, y=x, color=time, group=1)) +
      geom_point(size=3) +
      geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") +
      labs(x="time", y="mds1")
    
    # Third plot: MDS2 vs time
    p2 <- df.mds %>%
      ggplot(aes(x=ind, y=y, color=time, group=1)) +
      geom_point(size=3) +
      geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") +
      labs(x="time", y="mds2")
    
    # Fourth plot: MDS3 vs time
    p3 <- df.mds %>%
      ggplot(aes(x=ind, y=z, color=time, group=1)) +
      geom_point(size=3) +
      geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") +
      labs(x="time", y="mds3")
    
    # Fifth plot: MDS1 vs MDS2
    p4 <- df.mds %>%
      ggplot(aes(x=x, y=y, color=time)) +
      geom_point(size=3) +
      geom_label_repel(aes(label=time), size=2) +
      theme(legend.position = "none") +
      labs(x="mds1", y="mds2")
    
    # Arrange all four ggplots in a 1x4 grid
    grid.arrange(
      big_title,  # The title at the top
      arrangeGrob(p0, p1, p2, p3, ncol=4),  # Arrange plots in a 2x2 layout
      ncol=1, heights=c(0.3, 4)  # Adjust title height relative to plots
    )    
    
  }
  
  return(list(mds=mds, df.mds=df.mds))
}

#apply ISOMAP on the CMDS result with chosen dimension mdsd from CMDS step defaultly it always embeds to 1 
doIso <- function(mds, mdsd=2, isod=1, doplot=F)
{
  df.iso <- NULL
  dis <- vegdist(mds[,1:mdsd,drop=F], "euclidean")
  knn <- 1
  success <- FALSE
  while(!success) {
    tryCatch({
      iso = isomap(dis, k=knn, ndim=isod, path="shortest")$points
      success <- TRUE
    },
    error = function(e) {
      knn <<- knn + 1
    })
  }
  iso2 <- tibble(iso=iso[,1]) %>% mutate(i=1:nrow(mds), ind=1:nrow(mds), time=1:nrow(mds), knn=knn)
  df.iso <- rbind(df.iso, cbind(iso2, mdsd=mdsd))
  df.iso <- df.iso %>% group_by(mdsd) %>% mutate(iso = if(iso[1] > 0) {-iso} else {iso}) %>% ungroup()
  
  if (doplot) {
    
    p <- df.iso %>% filter(mdsd==mdsd) %>%
      ggplot(aes(x = 1: length(iso) , y=iso, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      theme(legend.position = "none") +
      labs(x="time", y="isomap embedding") +
      # scale_x_date(breaks = scales::breaks_pretty(8), labels=label_date_short()) +
      theme(axis.text.x=element_text(hjust=0.7))
    # theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3),
    #       axis.text.y = element_text(size = 12),
    #       axis.title = element_text(size = 14, face="bold"))
    #    p <- p + scale_x_date(breaks = scales::breaks_pretty(8), labels=label_date_short())
    print(p)
    
    df.isok <- df.iso %>% filter(mdsd==mdsd) #%>% mutate(date2 = format(ymd(paste0(date,"-01")),"%m/%y"))
    row.names(df.isok) <- df.isok$time
    fit <- lm(iso ~ i, data=df.isok)
    # print(tidy(fit))
    # print(glance(fit))
    myfor <- augment(fit)
    myfor2 <- myfor %>% mutate(date=.rownames,
                               ranks = rank(.sigma),
                               mycol=sprintf("%2d",rank(.fitted)))
    p <- myfor2 %>%
      ggplot(aes(.fitted, .resid)) +
      geom_point(aes(color=mycol)) +
      geom_hline(yintercept = 0, linetype="dashed", color="grey") +
      geom_smooth(method="loess", se=FALSE) +
      labs(x="Fitted Values", y="Residuals") +
      theme(legend.position = "none",
            axis.title = element_text(size=14, face="bold"))
    p <- p + geom_label_repel(aes(label=date), data=myfor2 %>% filter(ranks %in% 1:3))
    print(p)
  }
  
  return(df.iso)
}

## This is another slope change point algorithm called segmented that is not used in the paper 
break_point_dection=function(D,k){
  tmax <- nrow(D)
  df.mds <- doMDS(D,doplot = F)
  mds <- df.mds$mds
  df.iso <- doIso(mds, mdsd=k)
  x=1:tmax/tmax
  y1=df.iso$iso
  os1<-segmented(lm(y1~x),psi=c(0.2))
  result=as.data.frame(matrix(0,1,3))
  result[1,1:4]=c((os1$psi[1,2]-1.95*os1$psi[1,3])*tmax,(os1$psi[1,2])*tmax,(os1$psi[1,2]+1.95*os1$psi[1,3])*tmax,os1$psi[1,3]*tmax)
  return(result)
  ## this result returns you the confidence interval of point estimation and the standard deviation
}


find_slope_changepoint_with_plot <- function(y, doplot = TRUE) {
  tmax <- length(y)
  x <- 1:tmax
  best_cp <- NULL
  min_loss <- Inf
  best_coefs <- NULL
  
  for (cp in 2:(tmax - 1)) {
    # Construct design matrix based on the given model
    X <- cbind(1, (x - cp), (x > cp) * (x - cp))
    
    # Fit the model using least squares
    fit <- lm(y ~ X - 1)  # "-1" removes intercept as it's already in design matrix
    
    # Get the coefficients and ensure beta_L != beta_R
    coef_fit <- coef(fit)
    alpha_hat <- coef_fit[1]
    beta_L_hat <- coef_fit[2]
    beta_R_hat <- coef_fit[3] + beta_L_hat  # Adjust for slope change
    
    if (beta_L_hat != beta_R_hat) {
      # Calculate sum of squared residuals
      fitted_values <- alpha_hat + beta_L_hat * (x - cp) + (beta_R_hat - beta_L_hat) * (x - cp) * (x > cp)
      loss <- sum((y - fitted_values)^2)
      
      # Update best changepoint if this loss is the minimum
      if (loss < min_loss) {
        min_loss <- loss
        best_cp <- cp
        best_coefs <- c(alpha = alpha_hat, beta_L = beta_L_hat, beta_R = beta_R_hat)
      }
    }
  }
  
  best_coefs
  # Generate fitted values using the best coefficients
  if (!is.null(best_coefs)) {
    fitted_y <- best_coefs["alpha.X1"] + best_coefs["beta_L.X2"] * (x - best_cp) + 
      (best_coefs["beta_R.X3"] - best_coefs["beta_L.X2"]) * (x - best_cp) * (x > best_cp)
  } else {
    fitted_y <- rep(NA, tmax)  # Return NA values if no changepoint was found
  }
  
  # Plot the results if doplot is TRUE
  if (doplot && !is.null(best_cp)) {
    plot(x, y, pch = 16, col = "black", main = "", xlab = "time", ylab = "mirror")
    lines(x, fitted_y, col = "red", lwd = 2)
    abline(v = best_cp, col = "red", lwd = 2, lty = 2)
    abline(v = tstar, col = "black", lwd = 2, lty = 2)
    legend("topleft", legend = c("Data", "Fitted Line l2", "Estimated_CP l2"),
           col = c("black", "red", "red"), pch = c(16, NA, NA), lty = c(NA, 1, 2), lwd = 2 , cex =0.5)
  }
  
  # Return results
  return(list(changepoint = best_cp, coefficients = best_coefs , error = (best_cp- tstar)/tmax   ))
}


## Implementation of the 3rd step in Algorithm 2 by recasting it as a linear programming problem.
## This function will return the objective function value Sk in the paper for a given change point t 
linf_cp=function(t,y,cp){
  n=length(t)
  nl=sum(t<cp)+1
  XL=matrix(1,nrow = nl,ncol=4)
  XL[,4]=0
  XL[,3]=t[1:nl]-cp
  
  #XL; y[1:nl]
  
  XL2=XL
  XL2[,2:3]=-XL[,2:3]
  
  #rbind(XL,XL2); c(y[1:nl],-y[1:nl])
  
  XR=matrix(1,nrow = n-nl,ncol = 4)
  XR[,3]=0
  XR[,4]=t[(nl+1):n]-cp
  
  XR2=XR
  XR2[,c(2,4)]=-XR[,c(2,4)]
  
  
  X=rbind(XL,XR,XL2,XR2)
  Y=c(y,-y)
  
  library(lpSolveAPI)
  lprec <- make.lp(0,4)
  set.objfn(lprec,c(1,0,0,0))
  for (i in 1:(nrow(X)) ) {
    add.constraint(lprec, X[i,], ">=", Y[i])
  }
  
  set.bounds(lprec, lower = c(0,-Inf,-Inf,-Inf), columns = c(1,2,3,4))
  ColNames <- c('Z', "alpha", "bl","br")
  dimnames(lprec)[[2]] <- ColNames
  solve(lprec)
  return(get.variables(lprec))
}


linf_error=function(x){
  obf=NULL
  for (nk in 2:(tmax-1)) { ## find the point which minimize the obj func Sk, that is the change point 
    obf[nk]=linf_cp(1:tmax,x,nk)[1]
  }
  ecp=min(which(obf==min(obf[-1])))
  return( c((ecp-tstar)/tmax ,  ecp) )
}


shuffle_X <- function(X,del){
  n=nrow(X)
  dn=floor(del*n)
  permu_vec=sample(1:dn)
  random_perm=diag(dn)[permu_vec,]
  a=bdiag(diag(n-dn), random_perm)
  random_perm=as.matrix(a)
  return(random_perm %*% X)
}

shuffle_graph_delta_perc <- function(A,del=0.2){
  G=as.matrix(A)
  n=nrow(G)
  dn=floor(del*n)
  permu_vec=sample(1:dn)
  random_perm=diag(dn)[permu_vec,]
  a=bdiag(diag(n-dn), random_perm)
  random_perm=as.matrix(a)
  permu_G=as.matrix(random_perm%*%G%*%t(random_perm))
  G_graph=graph_from_adjacency_matrix(permu_G,mode ="undirected")
  return(G_graph)
}



getIso = function(glist, dmax=2, dhat=2, mdsdmax=100, mdsd=2, isod=1, shuffled.v=0, doaug=TRUE, doptr=TRUE)
{
  
  pacman::p_load(vegan, doParallel)
  n = vcount(glist[[1]])
  m = length(glist)
  registerDoParallel(min(50, detectCores()-2))
  
  Xhat = lapply(glist, function(x) full.ase(x, d=dmax, diagaug=doaug, doptr=doptr)$Xhat)
  cat("** ASE is done! **\n")
  
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    #        cat("i = ", i, ", j = ", j, "\n")
    
    Xhat1 = Xhat[[i]][,1:dhat,drop=F]
    Xhat2 = Xhat[[j]][,1:dhat,drop=F]
    proc <- procrustes2(as.matrix(Xhat1), as.matrix(Xhat2))
    Xhat1 <- Xhat1 %*% proc$W
    D <- norm(Xhat1 - Xhat2, type="2")^2/n
    tibble(i=i, j=j, D=D)
  }
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- (D2 + t(D2)) / 1
  D2 <- sqrt(D2)
  cat("** D is done! **\n")
  
  # do mds
  mdsdmax = min(mdsdmax, m-1)
  mds <- cmdscale(D2, mdsdmax)
  cat("** CMDS is done! **\n")
  
  # do isomap
  dis <- vegdist(mds[,1:mdsd,drop=F], "euclidean")
  knn <- 1
  success <- FALSE
  while(!success) {
    tryCatch({
      iso = isomap(dis, k=knn, ndim=isod, path="shortest")$points
      success <- TRUE
    },
    error = function(e) {
      knn <<- knn + 1
    })
    # warning = function(e) {print(paste("warning : ", e))})
  }
  cat("** Isomap is done! **\n")
  
  iso2 <- tibble(time=1:nrow(mds), iso=iso[,1]) %>% mutate(i=1:nrow(mds), knn=knn)
  df.iso <- iso2 %>% mutate(mdsd=mdsd)
  df.iso <- df.iso %>% group_by(mdsd) %>%
    mutate(iso = if(iso[1] > 0) {-iso} else {iso}) %>% ungroup()
  
  return(list(D=D2, mds=mds, iso=df.iso))
}


calculate_w1_distance <- function(matrix1, matrix2) {
  if (ncol(matrix1) != ncol(matrix2)) {
    stop("The two matrices must have the same number of columns.")
  }
  
  n1 <- nrow(matrix1)
  n2 <- nrow(matrix2)
  
  weights1 <- rep(1/n1, n1)
  weights2 <- rep(1/n2, n2)
  
  cost_matrix <- as.matrix(dist(rbind(matrix1, matrix2)))[1:n1, (n1+1):(n1+n2)]
  
  transport_plan <- transport(weights1, weights2, costm = cost_matrix)
  
  indices <- cbind(transport_plan$from, transport_plan$to)
  costs_for_plan <- cost_matrix[indices]
  
  w1_distance <- sum(transport_plan$mass * costs_for_plan)
  
  return(w1_distance)
}

calculate_w2_distance <- function(matrix1, matrix2) {
  if (ncol(matrix1) != ncol(matrix2)) {
    stop("The two matrices must have the same number of columns.")
  }
  
  n1 <- nrow(matrix1)
  n2 <- nrow(matrix2)
  
  weights1 <- rep(1/n1, n1)
  weights2 <- rep(1/n2, n2)
  
  cost_matrix <- as.matrix(dist(rbind(matrix1, matrix2))^2)[1:n1, (n1+1):(n1+n2)]
  
  transport_plan <- transport(weights1, weights2, costm = cost_matrix)
  
  indices <- cbind(transport_plan$from, transport_plan$to)
  costs_for_plan <- cost_matrix[indices]
  
  w2_distance <- sqrt(sum(transport_plan$mass * costs_for_plan))
  
  return(w2_distance)
}

getDW1 <- function(Xlist, k=0, etype="proc") {
  m <- length(Xlist)
  n <- nrow(Xlist[[1]])
  if (k==0) {
    ind <- 1:n
  } else {
    # 'Yhat' must exist in the global environment if k != 0
    ind <- which(Yhat==k)
  }
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    
    if (etype == "proc") {
      Xhati <- Xlist[[i]][ind,]
      Xhatj <- Xlist[[j]][ind,]
      proc <- procrustes2(as.matrix(Xhati), as.matrix(Xhatj))
      Xhati <- Xhati %*% proc$W
    } else {
      Xhati <- Xlist[[i]][ind,]
      Xhatj <- Xlist[[j]][ind,]
    }
    
    # --- REPLACED CALCULATION ---
    # Calculate the W1 distance between the two matrices
    D <- calculate_w1_distance(Xhati, Xhatj)
    
    tibble(i=i, j=j, D=D)
  }
  
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- D2 + t(D2)
  
  D2
}



getDW2 <- function(Xlist, k=0, etype="proc") {
  m <- length(Xlist)
  n <- nrow(Xlist[[1]])
  if (k==0) {
    ind <- 1:n
  } else {
    # 'Yhat' must exist in the global environment if k != 0
    ind <- which(Yhat==k)
  }
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    
    if (etype == "proc") {
      Xhati <- Xlist[[i]][ind,]
      Xhatj <- Xlist[[j]][ind,]
      proc <- procrustes2(as.matrix(Xhati), as.matrix(Xhatj))
      Xhati <- Xhati %*% proc$W
    } else {
      Xhati <- Xlist[[i]][ind,]
      Xhatj <- Xlist[[j]][ind,]
    }
    
    # --- REPLACED CALCULATION ---
    # Calculate the W2 distance between the two matrices
    D <- calculate_w2_distance(Xhati, Xhatj)
    
    tibble(i=i, j=j, D=D)
  }
  
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- D2 + t(D2)
  
  D2
}


#####################



set.seed(1)
load("/.../glist-NRL-780-by0.1.RData")


del = c(0.5,1)

max_iter = 200

df_new <- tibble(g = glist) %>%
  mutate(Xhat = map(g, function(x) full.ase(as.matrix(as_adjacency_matrix(x)), 2)$Xhat[, 1:2, drop = FALSE]))

for(perc in del){
  df_new <- df_new %>%
    mutate(!!paste0("g_shuffle_", perc) := map(g, ~shuffle_graph_delta_perc(., perc))) 
  #%>%
    #mutate(!!paste0("g_shuffle_YP_", perc) := map(g, ~add_noise_g(., perc)))
  
  df_new <- df_new %>%
    mutate(!!paste0("shuffle_g_GM_alltoone", perc) := map(!!sym(paste0("g_shuffle_", perc)),
                                                          ~graph_mathing(df_new[[paste0("g_shuffle_", perc)]][[1]], .x, max_iter))) 
  #%>%
  #  mutate(!!paste0("shuffle_g_YP_GM_alltoone", perc) := map(!!sym(paste0("g_shuffle_YP_", perc)), 
  #                                                        ~graph_matching_YP(df_new[[paste0("g_shuffle_YP_", perc)]][[1]], .x)))
}


tstar = c(800,810)

df.iso_raw = getIso(df_new$g, dmax=2, dhat=2, mdsd=2, isod=1, doptr=FALSE) # graphs are symmetric, binary, hollow.
df.iso_shuffle1 = getIso(df_new$g_shuffle_1, dmax=2, dhat=2, mdsd=2, isod=1, doptr=FALSE) # graphs are symmetric, binary, hollow.
df.iso_shuffle_g_GM_alltoone1 = getIso(df_new$shuffle_g_GM_alltoone1, dmax=2, dhat=2, mdsd = 3, isod=1, doptr=FALSE) # graphs are symmetric, binary, hollow.


x_values <- seq(780, 820, by = 0.1)
data <- data.frame(x = x_values, y = df.iso_raw$iso$iso)
lm_model <- lm(y ~ x, data = data)
tmp = selgmented(lm_model, seg.Z = ~ x, type = "bic", Kmax = 20, msg = TRUE, plot.ic = TRUE)
summary(tmp)
p0 <- df.iso_raw$iso|> 
  ggplot(aes(x=seq(780,820,by=0.1), y=iso)) + 
  geom_line(color="grey") +
  geom_point(alpha=0.5, color="grey") + 
  geom_vline(xintercept = tstar, color="black", linetype="dashed", linewidth=1) +
  #geom_smooth(method = "loess", se=T, span=0.2, fill="blue", alpha=0.2) +
  #    geom_line(data=df.iso3, aes(x=x,y=y), color="red", linewidth=1) +
  geom_line(data=tibble(x=seq(780,820,by=0.1), y=broken.line(tmp)$fit), aes(x=x,y=y), color="red", linewidth=1) +
  labs(x="time", y="iso-mirror", title = expression(d['MV']~"with true alignment") ) +
  guides(color=guide_legend(title="")) +
  theme(legend.position = c(0.15,0.95),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height=unit(0.8,"line"),
        legend.background = element_rect(fill=NA, color=NA))

print(p0)

df.iso_shuffle1$iso$time = seq(780, 820, by = 0.1)

x_values <- seq(780, 820, by = 0.1)
data <- data.frame(x = x_values, y = df.iso_shuffle1$iso$iso)
lm_model <- lm(y ~ x, data = data)
tmp = selgmented(lm_model, seg.Z = ~ x, type = "bic", Kmax = 20, msg = TRUE, plot.ic = TRUE)
summary(tmp)
linear_fit <- lm(iso ~ time, data = df.iso_shuffle1$iso)

print(linear_fit)
plot <- plot + 
  geom_line(aes(y = predict(linear_fit)), color = "red", linewidth = 1)

p1 <- df.iso_shuffle1$iso|> 
  ggplot(aes(x=seq(780,820,by=0.1), y=iso)) + 
  geom_line(color="grey") +
  geom_point(alpha=0.5, color="grey") + 
  geom_vline(xintercept = tstar, color="black", linetype="dashed", linewidth=1) +
  #geom_smooth(method = "loess", se=T, span=0.2, fill="blue", alpha=0.2) +
  #    geom_line(data=df.iso3, aes(x=x,y=y), color="red", linewidth=1) +
  #geom_line(data=tibble(x=df.iso_shuffle_g_YP_GM_alltoone0.5$iso$time, y=broken.line(tmp)$fit), aes(x=x,y=y), color="red", linewidth=1) +
  geom_line(aes(y = predict(linear_fit)), color = "red", linewidth = 1)+
  labs(x="time", y="iso-mirror", title =  expression(d['MV']~"with 100% shuffling")) +
  guides(color=guide_legend(title="")) +
  theme(legend.position = c(0.15,0.95),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height=unit(0.8,"line"),
        legend.background = element_rect(fill=NA, color=NA))
print(p1)

x_values <- seq(780, 820, by = 0.1)
data <- data.frame(x = x_values, y = df.iso_shuffle_g_GM_alltoone1$iso$iso)
lm_model <- lm(y ~ x, data = data)
tmp = selgmented(lm_model, seg.Z = ~ x, type = "bic", Kmax = 20, msg = TRUE, plot.ic = TRUE)
summary(tmp)
p2 <-df.iso_shuffle_g_GM_alltoone1$iso|> 
  ggplot(aes(x=seq(780,820,by=0.1), y=iso)) + 
  geom_line(color="grey") +
  geom_point(alpha=0.5, color="grey") + 
  geom_vline(xintercept = tstar, color="black", linetype="dashed", linewidth=1) +
  #geom_smooth(method = "loess", se=T, span=0.2, fill="blue", alpha=0.2) +
  #    geom_line(data=df.iso3, aes(x=x,y=y), color="red", linewidth=1) +
  geom_line(data=tibble(x=seq(780,820,by=0.1), y=broken.line(tmp)$fit), aes(x=x,y=y), color="red", linewidth=1) +
  labs(x="time", y="iso-mirror", title =  expression(d['MV']~"with all to one graph matching on 100% shuffling")  ) +
  guides(color=guide_legend(title="")) +
  theme(legend.position = c(0.15,0.95),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height=unit(0.8,"line"),
        legend.background = element_rect(fill=NA, color=NA))
print(p2)



D_W1 = getDW1(df_new$Xhat)
W1mds <-  doMDS(D_W1, doplot = F)
W1_iso <- doIso(W1mds$mds, mdsd = 1, doplot = F)

x_values <- seq(780, 820, by = 0.1)
data <- data.frame(x = x_values, y = W1_iso$iso)
lm_model <- lm(y ~ x, data = data)
tmp = selgmented(lm_model, seg.Z = ~ x, type = "bic", Kmax = 20, msg = TRUE, plot.ic = TRUE)
summary(tmp)

tstar = c(800,810)

p_W1 <- W1_iso|> 
  ggplot(aes(x=seq(780,820,by=0.1), y=iso)) + 
  geom_line(color="grey") +
  geom_point(alpha=0.5, color="grey") + 
  geom_vline(xintercept = tstar, color="black", linetype="dashed", linewidth=1) +
  #geom_smooth(method = "loess", se=T, span=0.2, fill="blue", alpha=0.2) +
  #    geom_line(data=df.iso3, aes(x=x,y=y), color="red", linewidth=1) +
  geom_line(data=tibble(x=seq(780,820,by=0.1), y=broken.line(tmp)$fit), aes(x=x,y=y), color="red", linewidth=1) +
  labs(x="time", y="iso-mirror", title = expression(W[1])) +
  guides(color=guide_legend(title="")) +
  theme(legend.position = c(0.15,0.95),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height=unit(0.8,"line"),
        legend.background = element_rect(fill=NA, color=NA))

print(p_W1)


D_W2 = getDW2(df_new$Xhat)
W2mds <-  doMDS(D_W2, doplot = F)

W2_iso <- doIso(W2mds$mds, mdsd = 1, doplot = F)
x_values <- seq(780, 820, by = 0.1)
data <- data.frame(x = x_values, y = W2_iso$iso)
lm_model <- lm(y ~ x, data = data)
tmp = selgmented(lm_model, seg.Z = ~ x, type = "bic", Kmax = 20, msg = TRUE, plot.ic = TRUE)
summary(tmp)

tstar = c(800,810)

p_W2 <- W2_iso|> 
  ggplot(aes(x=seq(780,820,by=0.1), y=iso)) + 
  geom_line(color="grey") +
  geom_point(alpha=0.5, color="grey") + 
  geom_vline(xintercept = tstar, color="black", linetype="dashed", linewidth=1) +
  #geom_smooth(method = "loess", se=T, span=0.2, fill="blue", alpha=0.2) +
  #    geom_line(data=df.iso3, aes(x=x,y=y), color="red", linewidth=1) +
  geom_line(data=tibble(x=seq(780,820,by=0.1), y=broken.line(tmp)$fit), aes(x=x,y=y), color="red", linewidth=1) +
  labs(x="time", y="iso-mirror", title = expression(W[2])) +
  guides(color=guide_legend(title="")) +
  theme(legend.position = c(0.15,0.95),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height=unit(0.8,"line"),
        legend.background = element_rect(fill=NA, color=NA))

print(p_W2)




df_new <- df_new %>%
  mutate(avg_deg = map(g, function(x) sum(as.matrix(as_adjacency_matrix(x)))/vcount(x))   )

plot(seq(780,820,by=0.1),unlist(df_new$avg_deg))
x_values <- seq(780, 820, by = 0.1)
data <- data.frame(x = x_values, y = unlist(df_new$avg_deg))
lm_model <- lm(y ~ x, data = data)
tmp = selgmented(lm_model, seg.Z = ~ x, type = "bic", Kmax = 20, msg = TRUE, plot.ic = TRUE)
summary(tmp)


p_avg_deg <- data|> 
  ggplot(aes(x=seq(780,820,by=0.1), y=y)) + 
  geom_line(color="grey") +
  geom_point(alpha=0.5, color="grey") + 
  geom_vline(xintercept = tstar, color="black", linetype="dashed", linewidth=1) +
  #geom_smooth(method = "loess", se=T, span=0.2, fill="blue", alpha=0.2) +
  #    geom_line(data=df.iso3, aes(x=x,y=y), color="red", linewidth=1) +
  geom_line(data=tibble(x=seq(780,820,by=0.1), y=broken.line(tmp)$fit), aes(x=x,y=y), color="red", linewidth=1) +
  labs(x="time", y="average deg", title = 'average degree') +
  guides(color=guide_legend(title="")) +
  theme(legend.position = c(0.15,0.95),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height=unit(0.8,"line"),
        legend.background = element_rect(fill=NA, color=NA))

print(p_avg_deg)

library(gridExtra)
grid.arrange(p0, p1, p2,p_W1,p_W2,p_avg_deg,nrow = 6)

grid.arrange(p0, p1, p2, nrow = 3)

g <- arrangeGrob(p0, p1, p2,p_W1,p_W2,p_avg_deg, nrow = 6)

ggsave(
  filename = "NRL.pdf",
  plot     = g,
  width    = 6,
  height   = 11,
  units    = "in",
  device   = cairo_pdf  # or just "pdf"
)

#####################

df.iso_shuffle_g_YP_GM_alltoone1 = getIso(df_new$shuffle_g_YP_GM_alltoone1, dmax=2, dhat=2, mdsd=2, isod=1, doptr=FALSE) # graphs are symmetric, binary, hollow.
tmp = selgmented(df.iso_shuffle_g_YP_GM_alltoone1$iso$iso, type="bic", Kmax=20, msg=TRUE, plot.ic=T)
p2 <-df.iso_shuffle_g_YP_GM_alltoone1$iso|> 
  ggplot(aes(x=time, y=iso)) + 
  geom_line(color="grey") +
  geom_point(alpha=0.5, color="grey") + 
  geom_vline(xintercept = tstar, color="black", linetype="dashed", linewidth=1) +
  geom_smooth(method = "loess", se=T, span=0.2, fill="blue", alpha=0.2) +
  #    geom_line(data=df.iso3, aes(x=x,y=y), color="red", linewidth=1) +
  geom_line(data=tibble(x=df.iso_shuffle_g_YP_GM_alltoone0.5$iso$time, y=broken.line(tmp)$fit), aes(x=x,y=y), color="red", linewidth=1) +
  labs(x="time", y="iso-mirror", title = 'YP GM after shuffling 100% percent') +
  guides(color=guide_legend(title="")) +
  theme(legend.position = c(0.15,0.95),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height=unit(0.8,"line"),
        legend.background = element_rect(fill=NA, color=NA))
print(p2)

df.iso_shuffle_YP_1 = getIso(df_new$g_shuffle_YP_1, dmax=2, dhat=2, mdsd=2, isod=1, doptr=FALSE) # graphs are symmetric, binary, hollow.
tmp = selgmented(df.iso_shuffle_YP_1$iso$iso, type="bic", Kmax=20, msg=TRUE, plot.ic=T)
#broken.line(tmp)
linear_fit <- lm(iso ~ time, data = df.iso_shuffle_YP_1$iso)

# Add the fitted linear regression line to the plot
plot <- plot + 
  geom_line(aes(y = predict(linear_fit)), color = "red", linewidth = 1)

p1 <- df.iso_shuffle_YP_1$iso|> 
  ggplot(aes(x=time, y=iso)) + 
  geom_line(color="grey") +
  geom_point(alpha=0.5, color="grey") + 
  geom_vline(xintercept = tstar, color="black", linetype="dashed", linewidth=1) +
  geom_smooth(method = "loess", se=T, span=0.2, fill="blue", alpha=0.2) +
  #    geom_line(data=df.iso3, aes(x=x,y=y), color="red", linewidth=1) +
  #geom_line(data=tibble(x=df.iso_shuffle_g_YP_GM_alltoone0.5$iso$time, y=broken.line(tmp)$fit), aes(x=x,y=y), color="red", linewidth=1) +
  geom_line(aes(y = predict(linear_fit)), color = "red", linewidth = 1)+
  labs(x="time", y="iso-mirror", title = 'YP shuffle 100% percent') +
  guides(color=guide_legend(title="")) +
  theme(legend.position = c(0.15,0.95),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height=unit(0.8,"line"),
        legend.background = element_rect(fill=NA, color=NA))

