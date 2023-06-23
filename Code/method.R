# may be redundant 

library(RColorBrewer)
library(colourvalues)
library(grDevices)
library(SDMTools)
library(network)
library(ggraph)
library(tidygraph)
library(dplyr)
library(gasper)
library(readxl)
library(forecast)
library(ggfortify)
library(Metrics)
library(GNAR)
library(DTWBI)
library(vars)
library(geosphere)
library(xlsx)
library(scales)
library(igraph)


#######################################
## LOCAAT 
#######################################

integrate_on_graph <- function(f, G, S){
  res <- 0
  for(i in S){
    nbrs <- which(G$A[i,]!=0)
    wt <- sum(G$dist[i,nbrs])
    res <- res + wt*f[i]
  }
  
  return(res)
}


update_graph <- function(G, i, nbrs, given){
  # given = TRUE: spatial data without prespecified graph structure
  
  if(given==FALSE){
    # remove link connected to node i
    G$A[i,] <- 0
    G$A[,i] <- 0
    if(is.vector(G$sA[((G$sA[,1]!=i) & (G$sA[,2]!=i)),])){
      G$sA <- matrix(G$sA[((G$sA[,1]!=i) & (G$sA[,2]!=i)),], nrow=1)
    } else{
      G$sA <- G$sA[((G$sA[,1]!=i) & (G$sA[,2]!=i)),]
    }
    
    A <- G$sdist[((G$sdist[,1] %in% nbrs) & (G$sdist[,2] %in% nbrs)),]
    if(is.vector(A)){
      A <- matrix(A, nrow=1)
      G1 <- graph.data.frame(matrix(A[,1:2], nrow=1), directed=FALSE)
    } else{
      G1 <- graph.data.frame(A[,1:2], directed=FALSE)
    }
    
    E(G1)$weight <- A[,3]
    mst <- minimum.spanning.tree(G1)
    
    edge.wt <- igraph::as_data_frame(mst, what="edges")
    tmp <- nrow(edge.wt)
    edge.wt <- sapply(edge.wt, as.numeric)
    
    # remove all linked structure within neighbors
    G$A[nbrs, nbrs] <- 0
    G$sA <- G$sA[(!((G$sA[,1] %in% nbrs) & (G$sA[,2] %in% nbrs))),]
    if(tmp!=0){
      edge.wt <- matrix(edge.wt, ncol=3)
      edge.wt[,3] <- exp(-edge.wt[,3]/20)
      G$sA <- rbind(G$sA, edge.wt)
      for(j in 1:tmp){
        G$A[edge.wt[j,1], edge.wt[j,2]] <- edge.wt[j,3]
        G$A[edge.wt[j,2], edge.wt[j,1]] <- edge.wt[j,3]
      }
    }
  } else{
    # remove link connected to node i
    G$A[i,] <- 0
    G$A[,i] <- 0
    if(is.vector(G$sA[((G$sA[,1]!=i) & (G$sA[,2]!=i)),])){
      G$sA <- matrix(G$sA[((G$sA[,1]!=i) & (G$sA[,2]!=i)),], nrow=1)
    } else{
      G$sA <- G$sA[((G$sA[,1]!=i) & (G$sA[,2]!=i)),]
    }
    
    if(length(nbrs)!=1){
      for(j in 1:(length(nbrs)-1)){
        for(k in (j+1):length(nbrs)){
          if(G$A[nbrs[j], nbrs[k]]==0){
            G$A[nbrs[j], nbrs[k]] <- exp(-(G$dist[nbrs[j], i]+G$dist[nbrs[k], i]))
            G$A[nbrs[k], nbrs[j]] <- exp(-(G$dist[nbrs[j], i]+G$dist[nbrs[k], i]))
            G$sA <- rbind(G$sA, c(nbrs[j], nbrs[k], 
                                  exp(-(G$dist[nbrs[j], i]+G$dist[nbrs[k], i]))))
            G$sA <- rbind(G$sA, c(nbrs[k], nbrs[j], 
                                  exp(-(G$dist[nbrs[j], i]+G$dist[nbrs[k], i]))))
            G$dist[nbrs[j], nbrs[k]] <- G$dist[nbrs[j], i]+G$dist[nbrs[k], i]
            G$dist[nbrs[k], nbrs[j]] <- G$dist[nbrs[j], i]+G$dist[nbrs[k], i]
          }
        }
      }
      tmp <- as.data.frame(G$sA)
      tmp$V1 <- as.character(tmp$V1)
      tmp$V2 <- as.character(tmp$V2)
      tmp$V3 <- as.numeric(as.character(tmp$V3))
      G$sA <- tmp
    }
  }
  return(G)
}


# TO BE FIXED!!!!! doesn't work when stop <= 2
LOCAAT <- function(f, G, stop, given=FALSE){
  # edge weight : inv_dist / sum of inv_dist of neighbors
  n <- length(f)
  S = c(1:n) # indices of the scaling coef
  D = c() # indices of the detail coef
  c <- f # vector for scaling coef, initialized with signal values
  d <- rep(0,n) # vector for detail coef
  phi <- diag(1,n) # each row is phi_r,k(x)
  integrals <- rep(0,n)
  NBRS <- list()
  WT <- list()
  V <- list()
  for(p in 1:n){
    integrals[p] <- integrate_on_graph(phi[p,], G, S)
  }
  for(r in n:stop){
    # stage r
    
    ## split ##
    i <- S[which.min(integrals[S])] # node to be lifted
    nbrs <- which(G$A[i,]!=0)
    
    D <- c(i, D) # add node i to the set of detail coef indices
    S <- setdiff(S, i) # remove node i from the set of scaling coef indices
    NBRS[[r]] <- nbrs 
    
    ## predict ##
    # get wavelet coef
    wt <- 1/G$dist[i,nbrs]/sum(1/G$dist[i,nbrs]) # inverse distance prediction weights
    d[i] <- c[i] - sum(wt*c[nbrs])
    WT[[r]] <- wt
    
    ## update ##
    # update scaling func
    phi[nbrs,] <- phi[nbrs,] + matrix(wt, ncol=1) %*% phi[i,]
    integrals[nbrs] <- integrals[nbrs] + wt*integrals[i]
    
    denom <- sum(integrals[nbrs]^2)
    v <- integrals[i] / denom * integrals[nbrs]
    c[nbrs] <- c[nbrs] + v*d[i]
    V[[r]] <- v
    
    # update neighborhood structure
    G <- update_graph(G, i, nbrs, given)
    
    for(p in nbrs){
      integrals[p] <- integrate_on_graph(phi[p,], G, S) # integral update considering network structure change
    }
  }
  res <- list()
  res$S <- S
  res$D <- D
  res$d <- d[D]
  res$c <- c[S]
  res$phi <- phi
  res$nbrs <- NBRS
  res$wt <- WT
  res$v <- V
  
  return(res)
}


inverse_LOCAAT <- function(coef, info, stop){
  # info: result of LOCAAT for the given graph
  n <- length(coef)
  S <- info$S
  D <- info$D
  nbrs <- info$nbrs
  wt <- info$wt
  V <- info$v
  
  for(k in stop:n){
    i <- D[k+1-stop]
    j <- nbrs[[k]]
    w <- wt[[k]]
    v <- V[[k]]
    coef[j] <- coef[j] - v*coef[i]
    coef[i] <- coef[i] + sum(w*coef[j])
  }
  return(coef)
}



#######################################
## GNARI with intercept
#######################################

simulate.GNARfit0 <- function(object, nsim=object$frbic$time.in, seed=NULL,
                              future=TRUE, set.noise=NULL, allcoefs=FALSE, ...){
  stopifnot(is.GNARfit(object))
  stopifnot(floor(nsim)==nsim)
  stopifnot(nsim>0)
  if(!is.null(seed)){
    stopifnot(floor(seed)==seed)
    set.seed(seed)
  }
  if(!is.null(object$frbic$fact.var)){
    stop("fact.var not currently supported with simulate")
  }
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  if(!is.null(set.noise)){
    sig <- set.noise
  }else{
    sig <- sigma(object$mod)
  }
  if(!allcoefs){
    nas <- is.na(object$mod$coefficients)
    pvs <- summary(object$mod)$coefficients[,4] < 0.05
    vals <- rep(0, length(pvs))
    vals[pvs] <- summary(object$mod)$coefficients[pvs,1]
    coefvec <- rep(0, length(nas))
    coefvec[(!nas)] <- vals
    
  }else{
    coefvec <- object$mod$coefficients
    coefvec[is.na(coefvec)] <- 0
  }
  
  #use GNARsim
  if(object$frbic$globalalpha){
    #global alpha has one alpha per time lag
    alphaout <-  vector(mode="list", length=object$frbic$alphas.in)
    betaout <- as.list(rep(0,length=object$frbic$alphas.in))
    count <- 1
    for(jj in 1:object$frbic$alphas.in){
      alphaout[[jj]] <- rep(coefvec[count], object$frbic$nnodes)
      if(object$frbic$betas.in[jj]>0){
        betaout[[jj]] <- coefvec[(count+1):(count+object$frbic$betas.in[jj])]
      }
      count <- count + object$frbic$betas.in[jj] + 1
    }
    
  }else{
    #multiple alphas per time lag
    alphaout <-  vector(mode="list", length=object$frbic$alphas.in)
    betaout <- as.list(rep(0,length=object$frbic$alphas.in))
    count <- 1
    for(jj in 1:object$frbic$alphas.in){
      alphaout[[jj]] <- coefvec[count:(count+object$frbic$nnodes-1)]
      if(object$frbic$betas.in[jj]>0){
        betaout[[jj]] <- coefvec[(count+object$frbic$nnodes):(count+
                                                                object$frbic$nnodes+object$frbic$betas.in[jj]-1)]
      }
      count <- count + object$frbic$nnodes + object$frbic$betas.in[jj]
    }
  }
  if(!future){
    newser <- GNARsim(n=nsim, net = object$frbic$net.in,
                      alphaParams = alphaout, betaParams = betaout,
                      sigma=sig)
  }else{
    nnodes <- object$frbic$nnodes
    max.nei <- max(unlist(lapply(betaout, length)))
    nei.mats <- vector(mode="list", length=max.nei)
    #create weight matrices for neighbours
    #flip network so that NofNeighbours gives into node information
    netmat <- as.matrix(object$frbic$net.in, normalise=FALSE)
    if(!isSymmetric(netmat)){
      net <- as.GNARnet(t(netmat))
    }else{
      net <- object$frbic$net.in
    }
    for(ii in 1:max.nei){
      nei.mats[[ii]] <- as.matrix(x=net, stage=ii, normalise=TRUE)
      if(sum(nei.mats[[ii]])==0){
        warning("beta order too large for network, neighbour set ",ii," is empty")
      }
    }
    
    xx.init <- object$frbic$final.in
    ntimedep <- object$frbic$alphas.in
    stopifnot(nrow(xx.init)==ntimedep)
    
    xx.gen <- matrix(NA, nrow=nsim+ntimedep, ncol=nnodes)
    xx.gen[1:ntimedep,] <- xx.init
    
    for(tt in (ntimedep+1):(nsim+ntimedep)){
      
      for(ii in 1:ntimedep){
        if(ii==1){
          time.forecast <- alphaout[[ii]]*xx.gen[(tt-ii),]
        }else{
          tmpi <- alphaout[[ii]]*xx.gen[(tt-ii),]
          time.forecast <- time.forecast + tmpi
        }
        
        
      }
      
      nei.forecast <- 0
      beta.pos <- NULL
      for(aa in 1:ntimedep){
        bb <- length(betaout[[aa]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast <- nei.forecast + betaout[[aa]][dd]*xx.gen[tt-aa,]%*%t(nei.mats[[dd]])
          }
        }
      }
      xx.gen[tt,] <- time.forecast+nei.forecast+rnorm(n=object$frbic$nnodes, mean=0, sd=sig) + 
        object$mod$coefficients[(length(object$mod$coefficients)-object$frbic$nnodes+1):length(object$mod$coefficients)]
    }
    if(nsim==1){
      newser <- as.ts(t(xx.gen[(ntimedep+1):(nsim+ntimedep),]), start=1, end=nsim)
    }else{
      newser <- as.ts(xx.gen[(ntimedep+1):(nsim+ntimedep),], start=1, end=nsim)
    }
  }
  
  return(newser)
}

predict.GNARfit0 <- function(object, n.ahead=1, ...){
  stopifnot(is.GNARfit(object))
  if(!is.null(object$frbic$fact.var)){
    stop("fact.var not currently supported with predict")
  }
  return(simulate.GNARfit0(object, nsim=n.ahead, future=TRUE, set.noise=0, ...))
}



GNARfit0 <- function (vts = GNAR::fiveVTS, net = GNAR::fiveNet, alphaOrder = 2, 
                      betaOrder = c(1, 1), fact.var = NULL, globalalpha = TRUE, 
                      tvnets = NULL, netsstart = NULL, ErrorIfNoNei = TRUE) 
{
  stopifnot(is.GNARnet(net))
  stopifnot(ncol(vts) == length(net$edges))
  stopifnot(alphaOrder > 0)
  stopifnot(floor(alphaOrder) == alphaOrder)
  stopifnot(length(betaOrder) == alphaOrder)
  stopifnot(floor(betaOrder) == betaOrder)
  if (!is.null(fact.var)) {
    stopifnot(length(fact.var) == length(net$edges))
  }
  stopifnot(is.matrix(vts))
  stopifnot(is.logical(globalalpha))
  if (!is.null(tvnets)) {
    cat("Time-varying networks not yet supported")
  }
  stopifnot(is.null(tvnets))
  useNofNei <- 1
  frbic <- list(nnodes = length(net$edges), alphas.in = alphaOrder, 
                betas.in = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                xtsp = tsp(vts), time.in = nrow(vts), net.in = net, 
                final.in = vts[(nrow(vts) - alphaOrder + 1):nrow(vts), 
                ])
  dmat <- GNARdesign(vts = vts, net = net, alphaOrder = alphaOrder, 
                     betaOrder = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                     tvnets = tvnets, netsstart = netsstart)
  if (ErrorIfNoNei) {
    if (any(apply(dmat == 0, 2, all))) {
      parname <- strsplit(names(which(apply(dmat == 0, 
                                            2, all)))[1], split = NULL)[[1]]
      betastage <- parname[(which(parname == ".") + 1):(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", 
           betastage)
    }
  }
  predt <- nrow(vts) - alphaOrder
  yvec <- NULL
  for (ii in 1:length(net$edges)) {
    yvec <- c(yvec, vts[((alphaOrder + 1):(predt + alphaOrder)), 
                        ii])
  }
  tmp <- NULL
  tmp2 <- nrow(dmat) / length(net$edges)
  for(l in 1:length(net$edges)){
    tmp3 <- rep(0, nrow(dmat))
    tmp3[(tmp2*(l-1)+1) : (tmp2*l)] <- 1
    tmp <- cbind(tmp, tmp3)
  }
  colnames(tmp) <- paste("intercept", 1:length(net$edges), sep="")
  
  if (sum(is.na(yvec)) > 0) {
    yvec2 <- yvec[!is.na(yvec)]
    dmat2 <- dmat[!is.na(yvec), ]
    modNoIntercept <- lm(yvec2 ~ dmat2 + tmp + 0)
  }
  else {
    modNoIntercept <- lm(yvec ~ dmat + tmp + 0)
  }
  out <- list(mod = modNoIntercept, y = yvec, dd = dmat, frbic = frbic)
  class(out) <- "GNARfit"
  return(out)
}

NARIMAselect0 <- function(vts, net, max.alpha, max.beta, globalalpha = TRUE){
  minbic <- 10^5
  minalpha <- max.alpha
  for(i in 1:max.alpha){
    tryCatch({
      skip_to_next <<- FALSE
      modelbic <- BIC(GNARfit0(vts = vts, net = net, alphaOrder = i,
                               betaOrder = rep(0,i), globalalpha = globalalpha))
    }, error=function(e){skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    
    if(!(is.nan(modelbic))){
      if(is.nan(minbic)){
        minbic <- modelbic
        minalpha <- i
      } else{
        if(minbic > modelbic){
          minbic <- modelbic
          minalpha <- i
        }
      }
    }
  }
  
  blist <- expand.grid(rep(list(0:max.beta), minalpha))
  minbic <- 10^5
  minbeta <- as.numeric(blist[nrow(blist),])
  
  for(i in 1:nrow(blist)){
    tryCatch({
      skip_to_next <<- FALSE
      modelbic <- BIC(GNARfit0(vts = vts,
                               net = net, alphaOrder = minalpha, betaOrder = as.numeric(blist[i,]), globalalpha = globalalpha))
    }, error=function(e){skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    
    if(!(is.nan(modelbic))){
      if(is.nan(minbic)){
        minbic <- modelbic
        minbeta <- as.numeric(blist[i,])
      } else{
        if(minbic > modelbic){
          minbic <- modelbic
          minbeta <- as.numeric(blist[i,])
        }
      }
    }
  }
  
  res <- list()
  res$alphaOrder <- minalpha
  res$betaOrder <- minbeta
  
  return(res)
}



forecast_narima0 <- function(vts, h, N, net, max.alpha, max.beta, globalalpha = TRUE, centering=TRUE){
  if(centering){
    vmean <- rowMeans(t(vts)) 
    vts <- t(t(vts)-vmean)
  } 
  
  narimaorder.gnar <- NARIMAselect0(vts = vts, net = net, max.alpha = max.alpha, max.beta = max.beta, globalalpha = globalalpha)
  print(narimaorder.gnar)
  model.fit <- GNARfit0(vts = vts, net = net, alphaOrder = narimaorder.gnar$alphaOrder, 
                        betaOrder = narimaorder.gnar$betaOrder, globalalpha = globalalpha)
  print(model.fit)
  pred.gnar <- predict.GNARfit0(model.fit, n.ahead=h)
  if(centering){
    res <- t(pred.gnar) + vmean
  } else{
    res <- t(pred.gnar)
  }
  return(res)
}




###################################
## GFT-based
###################################

# GFT function
compute_gft <- function(f, evectors){
  res <- c()
  N <- nrow(f)
  for(i in 1:N){
    res<- rbind(res, colSums(evectors[,i]*f))
  }
  return(res)
}


forecast_gft <- function(fc, h, period){
  # fc: fourier coef matrix where columns are date, rows are evalues
  # h: number of periods for forecasting
  res <- matrix(0, nrow(fc), ncol=h)
  
  for(i in 1:nrow(fc)){
    data <- fc[i,]
    model.fit <- auto.arima(ts(data, frequency=period))
    res[i,] <- as.numeric(forecast(model.fit, h=h)$mean)
  }
  
  return(res)
}


# inverse GFT
inverse_gft <- function(fc, evectors){
  return(evectors %*% fc)
}




###################################
## SGWT-based
###################################

# for each node, scalewise ARMA forecasting (independency check theoretically)
forecast_sgwt.arima <- function(wcf, h, period){
  # wcf: wcf.data matrix where columns are date, rows are nodes, scales, values = wavelet coefficients
  # h: number of periods for forecasting
  res <- matrix(0, nrow(wcf), ncol=h)
  
  for(i in 1:nrow(wcf)){
    data <- wcf[i,]
    model.fit <- auto.arima(ts(data, frequency=period))
    res[i,] <- as.numeric(forecast(model.fit, h=h)$mean)
  }
  
  return(res)
}




####################################
## simulation
####################################

GNARsim.custom <- function(n=200, net=GNAR::fiveNet, initsignal, alphaParams=list(c(rep(0.2,5))), betaParams=list(c(0.5)), sigma=1,
                           tvnets=NULL, netsstart=NULL, drift=FALSE){
  #use 0s in betaParams to discount dependencies
  stopifnot(is.GNARnet(net))
  stopifnot(length(alphaParams)==length(betaParams))
  stopifnot(!is.null(betaParams))
  stopifnot(length(alphaParams[[1]])==length(net$edges))
  stopifnot(length(initsignal)==length(net$edges))
  
  if(!is.null(tvnets)){
    cat("Time-varying nets not currently supported")
  }
  stopifnot(is.null(tvnets))
  
  nnodes <- length(net$edges)
  max.nei <- max(unlist(lapply(betaParams, length)))
  nei.mats <- vector(mode="list", length=max.nei)
  #create weight matrices for neighbours
  #flip network so that NofNeighbours gives into node information
  netmat <- as.matrix(net, normalise=FALSE)
  if(!isSymmetric(netmat)){
    net <- as.GNARnet(t(netmat))
  }
  for(ii in 1:max.nei){
    nei.mats[[ii]] <- as.matrix(x=net, stage=ii, normalise=TRUE)
    if(sum(nei.mats[[ii]])==0){
      warning("beta order too large for network, neighbour set ",ii," is empty")
    }
  }
  #print(length(alphaParams))
  
  #seed the process from normal dist with mean 0 and sigma as given
  #do this for as many time points as needed for alpha
  ntimedep <- length(alphaParams)
  #print(length(alphaParams))
  
  xx.init <- matrix(initsignal, nrow=ntimedep, ncol=nnodes)
  
  xx.gen <- matrix(NA, nrow=n+50, ncol=nnodes)
  xx.gen[1:ntimedep,] <- xx.init
  # print(length(alphaParams))
  # print("just before the tt indexed loop")
  for(tt in (ntimedep+1):(n+50)){
    # print(tt)
    # print(alphaParams)
    # print(rev(alphaParams))
    for(ii in 1:ntimedep){
      if(ii==1){
        time.forecast <- alphaParams[[ii]]*xx.gen[(tt-ii),]
      }else{
        tmpi <- alphaParams[[ii]]*xx.gen[(tt-ii),]
        time.forecast <- time.forecast + tmpi
      }
      
      
    }
    
    nei.forecast <- 0
    beta.pos <- NULL
    for(aa in 1:ntimedep){
      bb <- length(betaParams[[aa]])
      if(bb>0){
        for(dd in 1:bb){
          nei.forecast <- nei.forecast + betaParams[[aa]][dd]*xx.gen[tt-aa,]%*%t(nei.mats[[dd]])
        }
      }
    }
    xx.gen[tt,] <- time.forecast+nei.forecast+rnorm(nnodes, mean=0, sd=sigma)
    if(!identical(drift,FALSE)){
      xx.gen[tt,] <- time.forecast+nei.forecast+rnorm(nnodes, mean=0, sd=sigma) + drift
    }
    
  }
  return(as.ts(xx.gen[51:(n+50),]))
}
