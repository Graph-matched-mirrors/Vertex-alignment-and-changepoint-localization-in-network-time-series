pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)
library(ggrepel)

##generate RDPG from given latent position X
rdpg.sample <- function(X, rdpg_scale=FALSE) {
  P <- X %*% t(X)
  if (rdpg_scale) {
    P <- scales::rescale(P)
  }
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  diag(U) <- runif(n)
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(graph_from_adjacency_matrix(A,"undirected"))
}

##ASE for a network A with embedding dimension d
full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  #    require(irlba)
  
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
    mutate( avg_edges = map( g ,  ~sum(as.matrix(.))    )   ) %>%
    mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
  
  for(perc in del){
    df <- df %>%
      mutate(!!paste0("shuffle_Xt_", perc) := map( Xt ,~shuffle_X(., perc)) )
  }
  
  df
}


## true dmv distance
true_distance=function(m,p){
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



getD_W1 <- function(Xlist) {
  Xlist= df$Xt
  m <- length(Xlist)
  ind <- 1:n
  
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    
    Xhati <- Xlist[[i]][ind,] 
    Xhatj <- Xlist[[j]][ind,]
    
    D <- mean( abs( sort(Xhatj) - sort(Xhati) ) )
    tibble(i=i, j=j, D=D)
  }
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- (D2 + t(D2)) / 1
  D2
}


getD_W2 <- function(Xlist) {
  Xlist= df$Xt
  m <- length(Xlist)
  ind <- 1:n
  
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    
    Xhati <- Xlist[[i]][ind,] 
    Xhatj <- Xlist[[j]][ind,]
    
    D <- sqrt(mean( ( sort(Xhatj) - sort(Xhati) )^2 ))
    tibble(i=i, j=j, D=D)
  }
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- (D2 + t(D2)) / 1
  D2
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
    grid.arrange(p0, p1, p2, p3, ncol=4)
    
    
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

## find the point in time series with significant change in the slope
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


# changepoint detection using L-infinity norm
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

## true distance matrix for Atlanta model
true_Atlanta_dmv <- function(p,q,num_state,m,tstar,delta){
  library(expm) 
  markov_p = matrix(0, num_state, num_state)
  markov_q = matrix(0, num_state, num_state)
  True_dmv_square = matrix(0, m, m)
  
  diag(markov_p) = 1 - 2 * p
  diag(markov_p[-1, ]) = p
  diag(markov_p[, -1]) = p
  markov_p[1, 1] = 1 - p
  markov_p[1, 2] = p
  markov_p[num_state, num_state - 1] = p
  markov_p[num_state, num_state] = 1 - p
  
  diag(markov_q) = 1 - 2 * q
  diag(markov_q[-1, ]) = q
  diag(markov_q[, -1]) = q
  markov_q[1, 1] = 1 - q
  markov_q[1, 2] = q
  markov_q[num_state, num_state - 1] = q
  markov_q[num_state, num_state] = 1 - q
  
  # Precompute the diff matrix
  diff = outer(1:num_state, 1:num_state, function(i, j) (i - j)^2)
  
  # Precompute Markov powers for all steps up to m
  markov_powers = lapply(0:(m - 1), function(x) markov_p %^% x)
  markov_q_powers = lapply(0:(m - 1), function(x) markov_q %^% x)
  
  # Calculate True_dmv_square
  for (s in 1:(m - 1)) {
    for (k in (s + 1):m) {
      if ((s < tstar) & (k > tstar)) {
        True_dmv_square[s, k] = sum((markov_q_powers[[abs(k - tstar) + 1]] %*% markov_powers[[abs(tstar - s) + 1]]) * diff) * delta^2 / num_state
      } else if ((s < tstar) & (k < tstar)) {
        True_dmv_square[s, k] = sum(markov_powers[[abs(s - k) + 1]] * diff) * delta^2 / num_state
      } else if (s == tstar) {
        True_dmv_square[s, k] = sum(markov_q_powers[[abs(s - k) + 1]] * diff) * delta^2 / num_state
      } else if (k == tstar) {
        True_dmv_square[s, k] = sum(markov_powers[[abs(s - k) + 1]] * diff) * delta^2 / num_state
      } else if ((s > tstar) & (k > tstar)) {
        True_dmv_square[s, k] = sum(markov_q_powers[[abs(s - k) + 1]] * diff) * delta^2 / num_state
      }
    }
  }
  
  
  True_dmv_square=True_dmv_square+t(True_dmv_square)
  return(True_dmv_square)
  
}

# true squared distance matrix for shufflted Atlanta 
true_shuffle_Atlanta_dmv <- function(c, num_state , m){
  the_var = c^2/12*( (num_state+1)/(num_state-1) )
  True_shuffle_Dmv_square = matrix( 2*the_var ,m,m)
  diag(True_shuffle_Dmv_square)=0
  return(True_shuffle_Dmv_square)
}


jump_Lbd=function(cp,p){
  if(runif(1)<1-p){
    np=cp
  } else {
    np=0.1+delta
  }
  return(np)
}


jump_Rbd=function(cp,p){
  if(runif(1)<1-p){
    np=cp
  } else {
    np=0.9-delta
  }
}

jump_middle=function(cp,p){
  u=runif(1)
  if(u<p){
    
    np=cp+delta
    
  }
  
  if(p<u & u<2*p){
    np=cp-delta
  }
  
  if(u>2*p){
    np=cp
  }
  return(np)
}


# Use jump functions to determine next position
update_function=function(current_position,p){
  if( abs(current_position-0.1)<10^(-10) ){
    next_position=jump_Lbd(current_position,p)
  }
  if(current_position> (0.1+10^(-10)) & current_position <(0.9-10^(-10)) ){
    next_position=jump_middle(current_position,p)
  }
  if(abs(current_position-0.9)<10^(-10)){
    next_position=jump_Rbd(current_position,p)
  }
  return(next_position)
}

shuffle_perc_graph <- function(A, del, n){
  G=as.matrix(A)
  print
  dn = floor(n*del)
  if(dn == 0)
    return(A)
  else{
    permu_vec=sample(1:dn)
    random_perm=diag(dn)[permu_vec,]
    a=bdiag(diag(n-dn), random_perm)
    random_perm = as.matrix(a)
    permu_G=as.matrix(random_perm%*%G%*%t(random_perm))
    G_graph=graph_from_adjacency_matrix(permu_G,mode ="undirected")
    return(G_graph)
  }
}

paired_error_in_shuffling_once <- function(n = 1000, p = 0.4, q = 0.15, m = 50, delta = 0.1, tstar = 25, del = c(0.1,0.2)){
  
  xt=matrix(0,nrow = n, ncol = m+1)
  initial_state_all_nodes=sample(seq(0,1,by=delta),n,replace = TRUE)
  xt[,1]=initial_state_all_nodes
  for (i in 1:n) {
    for (j in 2:(tstar)) {
      xt[i,j]=update_function(xt[i,j-1],p,delta)
    }
    for (j in ((tstar+1):(m+1)) ) {
      xt[i,j]=update_function(xt[i,j-1],q,delta)
    }
  }
  
  
  
  
  df <- tibble(time=1:m) %>%
    mutate(Xt = map(time, function(x) matrix(xt[,x],n,1)  )) %>%
    mutate(g = map(Xt, ~rdpg.sample(.)))%>%
    mutate(xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
  
  
  for(perc in del){
    df <- df %>%
      mutate(!!paste0("shuffle_", perc) := map(g,~shuffle_perc_graph(., perc, n))) %>%
      mutate(!!paste0("xhat_", perc) := map(!!sym(paste0("shuffle_", perc)), function(x) full.ase(x,2)$Xhat[,1,drop=F]))
  }
  
  ### make graph
  
  D2=getD(df$xhat)
  df.mds <- doMDS(D2,doplot = FALSE)
  df.iso <- doIso(df.mds, mdsd=10, AtlantaFlag=T)
  errors <- NULL
  errors[1] <- linf_error(df.iso$iso, m)
  i <- 2
  for(perc in del){
    D2_shuffle=getD(df[[paste0("xhat_", perc)]])
    df.mds_shuffle <- doMDS(D2_shuffle,doplot = FALSE)
    df.iso_shuffle <- doIso(df.mds_shuffle, mdsd=10, AtlantaFlag=T)
    errors[i] <- linf_error(df.iso_shuffle$iso, m)
    i <- i + 1
  }
  print(paste(n,q,Sys.time()))
  errors
}
paired_error_in_shuffling <- function(nmc = 50, n = 1000, p = 0.4, q = 0.15, m = 50, delta = 0.1, tstar = 25, del = 0.1){
  mc_errors <- sapply(1:nmc, function(i) paired_error_in_shuffling_once(n, p, q, m, delta, tstar, del))
  row_mse <- apply(mc_errors, 1, function(row) mean(abs(row)^2)) # mean(abs(errors[row,])^2)
  row_sds <- apply(mc_errors, 1, function(row) sd(abs(row)^2) / sqrt(nmc)) # sd(abs(errors[row,])^2) / sqrt(nmc)
  cbind(row_mse, row_sds)
}