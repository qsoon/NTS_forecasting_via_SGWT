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

# load functions
source("/home/kyu9510/NTS_forecasting_via_SGWT/Code/method.R", chdir=TRUE)

###################################################
## Scenario1 : Random Network Time Series
###################################################

set.seed(1)
er <- erdos.renyi.game(n=30, p=0.15, type="gnp")
# plot(er, vertex.size=6, vertex.label=NA)  
is_connected(er)
sig_list <- c(0.5, 1, 1.5)
for(sig in sig_list){
  plt0.scenario1.mat <- c()
  plt1.scenario1.mat <- c()
  plt2.scenario1.mat <- c()
  plt3.scenario1.mat <- c()
  
  plt0.scenario1.list <- list()
  plt1.scenario1.list <- list()
  plt2.scenario1.list <- list()
  plt3.scenario1.list <- list()
  
  set.seed(100*which(sig_list==sig))
  for(k in 1:20){
    print(paste("iteration:", k))
    scenario1.net <- igraphtoGNAR(er)
    
    scenario1.nts <- matrix(rnorm(220*30, 0, sig), nrow=220, ncol=30)
    colnames(scenario1.nts) <- 1:30
    
    # plot(scenario3.nts[1,], type="l")
    
    h1 <- 20
    N1 <- 30
    
    scenario1 <- list()
    A <- matrix(0, nrow=N1, ncol=N1)
    
    for(i in 1:length(scenario1.net$edges)){
      A[i,scenario1.net$edges[[i]]] <- scenario1.net$dist[[i]]
    }
    scenario1$A <- A
    
    train_scenario1 <- t(scenario1.nts[1:200,])
    true_scenario1 <- t(scenario1.nts[201:220,])
    
    pred.gnar.scenario1 <- forecast_narima0(vts = ts(t(train_scenario1)), h = h1, N = N1, 
                                            net = scenario1.net, max.alpha = 5, 
                                            max.beta = 3, globalalpha = TRUE, centering=FALSE)
    
    
    rownames(pred.gnar.scenario1) <- paste("Node",1:N1, sep="")
    colnames(pred.gnar.scenario1) <- paste("t", 1:h1, sep="")
    
    
    L1 <- diag(rowSums(scenario1$A)) - scenario1$A 
    val1 <- eigensort(L1)
    evalues.scenario1 <- val1$evalues
    evectors.scenario1 <- val1$evectors
    #- largest eigenvalue
    lmax.scenario1 <- max(evalues.scenario1)
    #- parameter that controls the scale number b <- 2
    b.scenario1 <- 3.5
    tf.scenario1 <- tight_frame(evalues.scenario1, evectors.scenario1, b=b.scenario1)
    
    m1 <- nrow(tf.scenario1)/N1
    
    wcf.scenario1 <- c() 
    
    # b: hyperparameter
    for(i in 1:ncol(train_scenario1)){
      f <- train_scenario1[,i]
      wcf <- forward_sgwt(f, evalues.scenario1, evectors.scenario1, b=b.scenario1)
      wcf.scenario1 <- cbind(wcf.scenario1, wcf)
    }
    colnames(wcf.scenario1) <- paste("t", 1:ncol(train_scenario1))
    rownames(wcf.scenario1) <- paste(paste("Node",1:N1, sep=""), rep(1:m1, each=N1), sep="")
    
    
    pred.wc.arima.scenario1 <- forecast_sgwt.arima(wcf.scenario1, h=h1, period=1)
    
    pred.sgwt.arima.scenario1 <- matrix(0, nrow=N1, ncol=h1) 

    rownames(pred.sgwt.arima.scenario1) <- paste("Node",1:N1, sep="")
    colnames(pred.sgwt.arima.scenario1) <- paste("t", 1:h1, sep="")
    
    for(i in 1:h1){
      pred.sgwt.arima.scenario1[,i] <- inverse_sgwt(pred.wc.arima.scenario1[,i], evalues.scenario1, evectors.scenario1,  b=b.scenario1)
    }
    
    
    fc.data.scenario1 <- compute_gft(train_scenario1, evectors.scenario1)
    
    true_fc.scenario1 <- compute_gft(true_scenario1, evectors.scenario1)
    
    # fc.data: N  x  train date   (for train data)
    # predicted fourier coefficients
    pred.fc.scenario1 <- forecast_gft(fc.data.scenario1, h=h1, period=1) # N by h (predicted data for h time points)
    # takes about 10 minutes
    
    # predict
    pred.gft.scenario1 <- inverse_gft(pred.fc.scenario1, evectors.scenario1) # predicted people num, row: stations, col: pred_date
    # rownames(pred.gft.covid.log.diff) <- district.name
    rownames(pred.gft.scenario1) <- paste("Node",1:N1, sep="")
    colnames(pred.gft.scenario1) <- paste("t", 1:h1, sep="")
    
    # nodewise ARIMA
    res <- matrix(0, nrow=N1, ncol=h1)
    for(i in 1:N1){
      model.fit <- auto.arima(ts(train_scenario1[i,], frequency=1))
      res[i,] <- as.numeric(forecast(model.fit, h=h1)$mean)
    }
    pred.nodewise.arima.scenario1 <- res
    
    
    plt0.scenario1 <- vector(length=20)
    for(i in 1:20){
      plt0.scenario1[i] <- rmse(c(true_scenario1[,i]), c(pred.nodewise.arima.scenario1[,i]))
    }
    
    plt1.scenario1 <- vector(length=20)
    for(i in 1:20){
      plt1.scenario1[i] <- rmse(c(true_scenario1[,i]), c(pred.gnar.scenario1[,i]))
    }
    
    plt2.scenario1 <- vector(length=20)
    for(i in 1:20){
      plt2.scenario1[i] <- rmse(c(true_scenario1[,i]), c(pred.sgwt.arima.scenario1[,i]))
    }
    
    plt3.scenario1 <- vector(length=20)
    for(i in 1:20){
      plt3.scenario1[i] <- rmse(c(true_scenario1[,i]), c(pred.gft.scenario1[,i]))
    }
    
    plt0.scenario1.list[[k]] <- pred.nodewise.arima.scenario1
    plt1.scenario1.list[[k]] <- pred.gnar.scenario1
    plt2.scenario1.list[[k]] <- pred.sgwt.arima.scenario1
    plt3.scenario1.list[[k]] <- pred.gft.scenario1
    
    plt0.scenario1.mat <- rbind(plt0.scenario1.mat, plt0.scenario1)
    plt1.scenario1.mat <- rbind(plt1.scenario1.mat, plt1.scenario1)
    plt2.scenario1.mat <- rbind(plt2.scenario1.mat, plt2.scenario1)
    plt3.scenario1.mat <- rbind(plt3.scenario1.mat, plt3.scenario1)
  }
  
  if(sig == sig_list[1]){
    plt0.scenario1.mat0 <- plt0.scenario1.mat
    plt1.scenario1.mat0 <- plt1.scenario1.mat
    plt2.scenario1.mat0 <- plt2.scenario1.mat
    plt3.scenario1.mat0 <- plt3.scenario1.mat
    
    plt0.scenario1.list0 <- plt0.scenario1.list
    plt1.scenario1.list0 <- plt1.scenario1.list
    plt2.scenario1.list0 <- plt2.scenario1.list
    plt3.scenario1.list0 <- plt3.scenario1.list
  }
  if(sig == sig_list[2]){
    plt0.scenario1.mat1 <- plt0.scenario1.mat
    plt1.scenario1.mat1 <- plt1.scenario1.mat
    plt2.scenario1.mat1 <- plt2.scenario1.mat
    plt3.scenario1.mat1 <- plt3.scenario1.mat
    
    plt0.scenario1.list1 <- plt0.scenario1.list
    plt1.scenario1.list1 <- plt1.scenario1.list
    plt2.scenario1.list1 <- plt2.scenario1.list
    plt3.scenario1.list1 <- plt3.scenario1.list
  }
  if(sig == sig_list[3]){
    plt0.scenario1.mat2 <- plt0.scenario1.mat
    plt1.scenario1.mat2 <- plt1.scenario1.mat
    plt2.scenario1.mat2 <- plt2.scenario1.mat
    plt3.scenario1.mat2 <- plt3.scenario1.mat
    
    plt0.scenario1.list2 <- plt0.scenario1.list
    plt1.scenario1.list2 <- plt1.scenario1.list
    plt2.scenario1.list2 <- plt2.scenario1.list
    plt3.scenario1.list2 <- plt3.scenario1.list
  }
}

# paper 
par(mfrow=c(1,3))
train_scenario1.list <- list()
true_scenario1.list <- list()
sig_list <- c(0.5, 1, 1.5)
k <- 0
for(sig in sig_list){
  k <- k + 1
  set.seed(100*which(sig_list==sig))
  scenario1.nts <- matrix(rnorm(220*30, 0, sig), nrow=220, ncol=30)
  colnames(scenario1.nts) <- 1:30
  
  train_scenario1.list[[k]] <- t(scenario1.nts[1:200,])
  true_scenario1.list[[k]] <- t(scenario1.nts[201:220,])
  plot(scenario1.nts[100,], type="l", main=bquote(sigma ==.(sig)), cex.main =2, ylab="", xlab="Node index", cex.lab=1.5, ylim=c(-4, 4))
}

par(mfrow=c(3,1))
train_scenario1.list <- list()
true_scenario1.list <- list()
sig_list <- c(0.5, 1, 1.5)
k <- 0
for(sig in sig_list){
  k <- k + 1
  set.seed(100*which(sig_list==sig))
  scenario1.nts <- matrix(rnorm(220*30, 0, sig), nrow=220, ncol=30)
  colnames(scenario1.nts) <- 1:30
  
  train_scenario1.list[[k]] <- t(scenario1.nts[1:200,])
  true_scenario1.list[[k]] <- t(scenario1.nts[201:220,])
  plot(scenario1.nts[,1], type="l", main=bquote(sigma ==.(sig)), cex.main =2, ylab="", xlab="time", cex.lab=1.8, ylim=c(-4, 4))
}

############ sigma = 0.5 #################
MM <-list()
plotdf <- data.frame(time = as.factor(rep(1:h1, each=20)), 
                     val = c(plt2.scenario1.mat0 - plt0.scenario1.mat0))

MM[[1]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="black", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 1: Proposed method RMSE - Nodewise ARIMA RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h1, each=20)), 
                     val = c(plt2.scenario1.mat0 - plt1.scenario1.mat0))

MM[[2]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="green", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 1: Proposed method RMSE - GNARI RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())

plotdf <- data.frame(time = as.factor(rep(1:h1, each=20)), 
                     val = c(plt2.scenario1.mat0 - plt3.scenario1.mat0))

MM[[3]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="deepskyblue3", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 1: Proposed method RMSE - GFT method RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


grid.arrange(MM[[1]], MM[[2]], MM[[3]],  nrow=3)



############ sigma = 1 #################
MM <-list()
plotdf <- data.frame(time = as.factor(rep(1:h1, each=20)), 
                     val = c(plt2.scenario1.mat1 - plt0.scenario1.mat1))

MM[[1]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="black", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 1: Proposed method RMSE - Nodewise ARIMA RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h1, each=20)), 
                     val = c(plt2.scenario1.mat1 - plt1.scenario1.mat1))

MM[[2]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="green", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 1: Proposed method RMSE - GNARI RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h1, each=20)), 
                     val = c(plt2.scenario1.mat1 - plt3.scenario1.mat1))

MM[[3]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="deepskyblue3", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 1: Proposed method RMSE - GFT method RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


grid.arrange(MM[[1]], MM[[2]], MM[[3]], nrow=3)


############ sigma = 1.5 #################
MM <-list()
plotdf <- data.frame(time = as.factor(rep(1:h1, each=20)), 
                     val = c(plt2.scenario1.mat2 - plt0.scenario1.mat2))

MM[[1]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="black", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 1: Proposed method RMSE - Nodewise ARIMA RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h1, each=20)), 
                     val = c(plt2.scenario1.mat2 - plt1.scenario1.mat2))

MM[[2]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="green", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 1: Proposed method RMSE - GNARI RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h1, each=20)), 
                     val = c(plt2.scenario1.mat2 - plt3.scenario1.mat2))

MM[[3]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="deepskyblue3", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 1: Proposed method RMSE - GFT method RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())



grid.arrange(MM[[1]], MM[[2]], MM[[3]], nrow=3)





###################################################################
## Scenario2 : Network Time Series with Inhomogeneous Structure
###################################################################

sig_list <- c(0.5, 1, 1.5)
for(sig in sig_list){
  plt0.scenario2.mat <- c()
  plt1.scenario2.mat <- c()
  plt2.scenario2.mat <- c()
  plt3.scenario2.mat <- c()
  plt3_0.scenario2.mat <- c()
  
  plt0.scenario2.list <- list()
  plt1.scenario2.list <- list()
  plt2.scenario2.list <- list()
  plt3.scenario2.list <- list()
  plt3_0.scenario2.list <- list()
  
  set.seed(100*which(sig_list==sig))
  for(k in 1:20){
    print(paste("iteration:", k))
    scenario2A <- list()
    N <- 30
    loc <- cbind(1:N, 0)
    scenario2A$xy <- loc
    
    # weight matrix
    A <- matrix(0, nrow=N, ncol=N)
    for(i in 1:(N-1)){
      A[i,i+1] <- 1
      A[i+1,i] <- 1
    }
    
    scenario2A$A <- A
    
    
    scenario2A.net <- matrixtoGNAR(scenario2A$A)
    
    scenario2A.nts <- GNARsim.custom(n=220, net=scenario2A.net, initsignal = rnorm(N,0,1), 
                                      alphaParams = list(rep(0.55,N)), betaParams = list(c(0.4)), sigma=sig,
                                      drift = rep(1,N))
    
    scenario2B <- list()
    N <- 10
    loc <- cbind((31):(30+N), 0)
    scenario2B$xy <- loc
    
    # weight matrix
    A <- matrix(0, nrow=N, ncol=N)
    for(i in 1:(N-1)){
      A[i,i+1] <- 1
      A[i+1,i] <- 1
    }
    
    scenario2B$A <- A
    
    
    scenario2B.net <- matrixtoGNAR(scenario2B$A)
    
    scenario2B.nts <- GNARsim.custom(n=220, net=scenario2B.net, initsignal = rnorm(N,0,1), alphaParams = list(rep(0.55,N)), 
                                      betaParams = list(c(0.1)), sigma=sig,
                                      drift = rep(3,N))
    
    scenario2C <- list()
    N <- 2
    loc <- cbind((41):(40+N), 0)
    scenario2C$xy <- loc
    
    # weight matrix
    A <- matrix(0, nrow=N, ncol=N)
    for(i in 1:(N-1)){
      A[i,i+1] <- 1
      A[i+1,i] <- 1
    }
    
    scenario2C$A <- A
    
    
    scenario2C.net <- matrixtoGNAR(scenario2C$A)
    
    scenario2C.nts <- GNARsim.custom(n=220, net=scenario2C.net, initsignal = rnorm(N,0,1), 
                                      alphaParams = list(rep(0.55,N)), betaParams = list(c(0.4)), sigma=sig,
                                      drift = rep(-0.5,N))
    
    
    scenario2D <- list()
    N <- 18
    loc <- cbind((43):(42+N), 0)
    scenario2D$xy <- loc
    
    # weight matrix
    A <- matrix(0, nrow=N, ncol=N)
    for(i in 1:(N-1)){
      A[i,i+1] <- 1
      A[i+1,i] <- 1
    }
    
    scenario2D$A <- A
    
    
    scenario2D.net <- matrixtoGNAR(scenario2D$A)
    
    scenario2D.nts <- GNARsim.custom(n=220, net=scenario2D.net, initsignal = rnorm(N,0,1), alphaParams = list(rep(0.55,N)), 
                                      betaParams = list(c(0.1)), sigma=sig,
                                      drift = rep(3,N))
    
    
    scenario2 <- list()
    N <- 60
    loc <- cbind(1:60, 0)
    scenario2$xy <- loc
    
    # weight matrix
    A <- matrix(0, nrow=N, ncol=N)
    for(i in 1:(N-1)){
      A[i,i+1] <- 1
      A[i+1,i] <- 1
    }
    A[30,31] <- 0.2
    A[31,30] <- 0.2
    A[40,41] <- 0.2
    A[41,40] <- 0.2
    A[42,43] <- 0.2
    A[43,42] <- 0.2
    
    scenario2$A <- A
    
    scenario2.net <- matrixtoGNAR(scenario2$A)
    
    scenario2.nts <- cbind(scenario2A.nts, scenario2B.nts , scenario2C.nts , scenario2D.nts)
    colnames(scenario2.nts) <- 1:60
    
    h2 <- 20
    N2 <- 60
    
    train_scenario2 <- t(scenario2.nts[1:200,])
    true_scenario2 <- t(scenario2.nts[201:220,])
    
    pred.gnar.scenario2 <- forecast_narima0(vts = ts(t(train_scenario2)), h = h2, N = N2, 
                                             net = scenario2.net, max.alpha = 5, 
                                             max.beta = 3, globalalpha = TRUE, centering=FALSE)
    rownames(pred.gnar.scenario2) <- paste("Node",1:N2, sep="")
    colnames(pred.gnar.scenario2) <- paste("t", 1:h2, sep="")
    
    
    L2 <- diag(rowSums(scenario2$A)) - scenario2$A 
    val2 <- eigensort(L2)
    evalues.scenario2 <- val2$evalues
    evectors.scenario2 <- val2$evectors
    #- largest eigenvalue
    lmax.scenario2 <- max(evalues.scenario2)
    #- parameter that controls the scale number b <- 2
    b.scenario2 <- 2.2
    tf.scenario2 <- tight_frame(evalues.scenario2, evectors.scenario2, b=b.scenario2)
    
    m2 <- nrow(tf.scenario2)/N2
    
    wcf.scenario2 <- c() 
    
    # b: hyperparameter
    for(i in 1:ncol(train_scenario2)){
      f <- train_scenario2[,i]
      wcf <- forward_sgwt(f, evalues.scenario2, evectors.scenario2, b=b.scenario2)
      wcf.scenario2 <- cbind(wcf.scenario2, wcf)
    }
    colnames(wcf.scenario2) <- paste("t", 1:ncol(train_scenario2))
    rownames(wcf.scenario2) <- paste(paste("Node",1:N2, sep=""), rep(1:m2, each=N2), sep="")
    
    
    pred.wc.arima.scenario2 <- forecast_sgwt.arima(wcf.scenario2, h=h2, period=1)
    
    pred.sgwt.arima.scenario2 <- matrix(0, nrow=N2, ncol=h2) 

    rownames(pred.sgwt.arima.scenario2) <- paste("Node",1:N2, sep="")
    colnames(pred.sgwt.arima.scenario2) <- paste("t", 1:h2, sep="")
    
    for(i in 1:h2){
      pred.sgwt.arima.scenario2[,i] <- inverse_sgwt(pred.wc.arima.scenario2[,i], evalues.scenario2, evectors.scenario2, b=b.scenario2)
    }
    
    
    fc.data.scenario2 <- compute_gft(train_scenario2, evectors.scenario2)

    # fc.data: N  x  train date   (for train data)
    # predicted fourier coefficients
    pred.fc.scenario2 <- forecast_gft(fc.data.scenario2, h=h2, period=1) # N by h (predicted data for h time points)
    
    # predict
    pred.gft.scenario2 <- inverse_gft(pred.fc.scenario2, evectors.scenario2) 
    rownames(pred.gft.scenario2) <- paste("Node",1:N2, sep="")
    colnames(pred.gft.scenario2) <- paste("t", 1:h2, sep="")
    
    
    L2_1 <- diag(rowSums(scenario2$A[1:30, 1:30])) - scenario2$A[1:30, 1:30] 
    val2_1<- eigensort(L2_1)
    evectors.scenario2_1 <- val2_1$evectors
    
    fc.data.scenario2_1 <- compute_gft(train_scenario2[1:30,], evectors.scenario2_1)
    
    # fc.data: N  x  train date   (for train data)
    # predicted fourier coefficients
    pred.fc.scenario2_1 <- forecast_gft(fc.data.scenario2_1, h=h2, period=1) # N by h (predicted data for h time points)
    # takes about 2 minutes
    
    # predict
    pred.gft.scenario2_1 <- inverse_gft(pred.fc.scenario2_1, evectors.scenario2_1) 
    rownames(pred.gft.scenario2_1) <- paste("Node",1:30, sep="")
    colnames(pred.gft.scenario2_1) <- paste("t", 1:h2, sep="")
    
    
    L2_2 <- diag(rowSums(scenario2$A[31:40, 31:40])) - scenario2$A[31:40, 31:40] 
    val2_2<- eigensort(L2_2)
    evectors.scenario2_2 <- val2_2$evectors
    
    fc.data.scenario2_2 <- compute_gft(train_scenario2[31:40,], evectors.scenario2_2)
    
    # fc.data: N  x  train date   (for train data)
    # predicted fourier coefficients
    pred.fc.scenario2_2 <- forecast_gft(fc.data.scenario2_2, h=h2, period=1) # N by h (predicted data for h time points)
    
    # predict
    pred.gft.scenario2_2 <- inverse_gft(pred.fc.scenario2_2, evectors.scenario2_2) 
    rownames(pred.gft.scenario2_2) <- paste("Node",31:40, sep="")
    colnames(pred.gft.scenario2_2) <- paste("t", 1:h2, sep="")
    
    L2_3 <- diag(rowSums(scenario2$A[41:42, 41:42])) - scenario2$A[41:42, 41:42] 
    val2_3<- eigensort(L2_3)
    evectors.scenario2_3 <- val2_3$evectors
    
    fc.data.scenario2_3 <- compute_gft(train_scenario2[41:42,], evectors.scenario2_3)
    
    # fc.data: N  x  train date   (for train data)
    # predicted fourier coefficients
    pred.fc.scenario2_3 <- forecast_gft(fc.data.scenario2_3, h=h2, period=1) # N by h (predicted data for h time points)
    # takes about 2 minutes
    
    # predict
    pred.gft.scenario2_3 <- inverse_gft(pred.fc.scenario2_3, evectors.scenario2_3) 
    rownames(pred.gft.scenario2_3) <- paste("Node",41:42, sep="")
    colnames(pred.gft.scenario2_3) <- paste("t", 1:h2, sep="")
    
    L2_4 <- diag(rowSums(scenario2$A[43:60, 43:60])) - scenario2$A[43:60, 43:60] 
    val2_4<- eigensort(L2_4)
    evectors.scenario2_4 <- val2_4$evectors
    
    fc.data.scenario2_4 <- compute_gft(train_scenario2[43:60,], evectors.scenario2_4)
    
    # fc.data: N  x  train date   (for train data)
    # predicted fourier coefficients
    pred.fc.scenario2_4 <- forecast_gft(fc.data.scenario2_4, h=h2, period=1) # N by h (predicted data for h time points)
    # takes about 2 minutes
    
    # predict
    pred.gft.scenario2_4 <- inverse_gft(pred.fc.scenario2_4, evectors.scenario2_4) 
    rownames(pred.gft.scenario2_4) <- paste("Node",43:60, sep="")
    colnames(pred.gft.scenario2_4) <- paste("t", 1:h2, sep="")
    
    
    
    pred.gft.scenario2_0 <- rbind(pred.gft.scenario2_1, pred.gft.scenario2_2, pred.gft.scenario2_3, pred.gft.scenario2_4)
    
    
    # nodewise ARIMA
    res <- matrix(0, nrow=N2, ncol=h2)
    for(i in 1:N2){
      model.fit <- auto.arima(ts(train_scenario2[i,], frequency=1))
      res[i,] <- as.numeric(forecast(model.fit, h=h2)$mean)
    }
    pred.nodewise.arima.scenario2 <- res
    
    
    plt0.scenario2 <- vector(length=20)
    for(i in 1:20){
      plt0.scenario2[i] <- rmse(c(true_scenario2[,i]), c(pred.nodewise.arima.scenario2[,i]))
    }
    
    plt1.scenario2 <- vector(length=20)
    for(i in 1:20){
      plt1.scenario2[i] <- rmse(c(true_scenario2[,i]), c(pred.gnar.scenario2[,i]))
    }
    
    plt2.scenario2 <- vector(length=20)
    for(i in 1:20){
      plt2.scenario2[i] <- rmse(c(true_scenario2[,i]), c(pred.sgwt.arima.scenario2[,i]))
    }
    
    plt3.scenario2 <- vector(length=20)
    for(i in 1:20){
      plt3.scenario2[i] <- rmse(c(true_scenario2[,i]), c(pred.gft.scenario2[,i]))
    }
    
    plt3_0.scenario2 <- vector(length=20)
    for(i in 1:20){
      plt3_0.scenario2[i] <- rmse(c(true_scenario2[,i]), c(pred.gft.scenario2_0[,i]))
    }
    
    plt0.scenario2.list[[k]] <- pred.nodewise.arima.scenario2
    plt1.scenario2.list[[k]] <- pred.gnar.scenario2
    plt2.scenario2.list[[k]] <- pred.sgwt.arima.scenario2
    plt3.scenario2.list[[k]] <- pred.gft.scenario2
    plt3_0.scenario2.list[[k]] <- pred.gft.scenario2_0
    
    plt0.scenario2.mat <- rbind(plt0.scenario2.mat, plt0.scenario2)
    plt1.scenario2.mat <- rbind(plt1.scenario2.mat, plt1.scenario2)
    plt2.scenario2.mat <- rbind(plt2.scenario2.mat, plt2.scenario2)
    plt3.scenario2.mat <- rbind(plt3.scenario2.mat, plt3.scenario2)
    plt3_0.scenario2.mat <- rbind(plt3_0.scenario2.mat, plt3_0.scenario2)
  }
  
  if(sig == sig_list[1]){
    plt0.scenario2.mat0 <- plt0.scenario2.mat
    plt1.scenario2.mat0 <- plt1.scenario2.mat
    plt2.scenario2.mat0 <- plt2.scenario2.mat
    plt3.scenario2.mat0 <- plt3.scenario2.mat
    plt3_0.scenario2.mat0 <- plt3_0.scenario2.mat
    
    plt0.scenario2.list0 <- plt0.scenario2.list
    plt1.scenario2.list0 <- plt1.scenario2.list
    plt2.scenario2.list0 <- plt2.scenario2.list
    plt3.scenario2.list0 <- plt3.scenario2.list
    plt3_0.scenario2.list0 <- plt3_0.scenario2.list
  }
  if(sig == sig_list[2]){
    plt0.scenario2.mat1 <- plt0.scenario2.mat
    plt1.scenario2.mat1 <- plt1.scenario2.mat
    plt2.scenario2.mat1 <- plt2.scenario2.mat
    plt3.scenario2.mat1 <- plt3.scenario2.mat
    plt3_0.scenario2.mat1 <- plt3_0.scenario2.mat
    
    plt0.scenario2.list1 <- plt0.scenario2.list
    plt1.scenario2.list1 <- plt1.scenario2.list
    plt2.scenario2.list1 <- plt2.scenario2.list
    plt3.scenario2.list1 <- plt3.scenario2.list
    plt3_0.scenario2.list1 <- plt3_0.scenario2.list
  }
  if(sig == sig_list[3]){
    plt0.scenario2.mat2 <- plt0.scenario2.mat
    plt1.scenario2.mat2 <- plt1.scenario2.mat
    plt2.scenario2.mat2 <- plt2.scenario2.mat
    plt3.scenario2.mat2 <- plt3.scenario2.mat
    plt3_0.scenario2.mat2 <- plt3_0.scenario2.mat
    
    plt0.scenario2.list2 <- plt0.scenario2.list
    plt1.scenario2.list2 <- plt1.scenario2.list
    plt2.scenario2.list2 <- plt2.scenario2.list
    plt3.scenario2.list2 <- plt3.scenario2.list
    plt3_0.scenario2.list2 <- plt3_0.scenario2.list
  }
}


par(mfrow=c(1,3))
sig_list <- c(0.5, 1, 1.5)
k <- 0
train_scenario2.list <- list()
true_scenario2.list <- list()
for(sig in sig_list){
  k <- k +1
  set.seed(100*which(sig_list==sig))
  scenario2A <- list()
  N <- 30
  loc <- cbind(1:N, 0)
  scenario2A$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  
  scenario2A$A <- A
  
  
  scenario2A.net <- matrixtoGNAR(scenario2A$A)
  
  scenario2A.nts <- GNARsim.custom(n=220, net=scenario2A.net, initsignal = rnorm(N,0,1), 
                                    alphaParams = list(rep(0.55,N)), betaParams = list(c(0.4)), sigma=sig,
                                    drift = rep(1,N))
  
  scenario2B <- list()
  N <- 10
  loc <- cbind((31):(30+N), 0)
  scenario2B$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  
  scenario2B$A <- A
  
  
  scenario2B.net <- matrixtoGNAR(scenario2B$A)
  
  scenario2B.nts <- GNARsim.custom(n=220, net=scenario2B.net, initsignal = rnorm(N,0,1), alphaParams = list(rep(0.55,N)), 
                                    betaParams = list(c(0.1)), sigma=sig,
                                    drift = rep(3,N))
  
  scenario2C <- list()
  N <- 2
  loc <- cbind((41):(40+N), 0)
  scenario2C$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  
  scenario2C$A <- A
  
  
  scenario2C.net <- matrixtoGNAR(scenario2C$A)
  
  scenario2C.nts <- GNARsim.custom(n=220, net=scenario2C.net, initsignal = rnorm(N,0,1), 
                                    alphaParams = list(rep(0.55,N)), betaParams = list(c(0.4)), sigma=sig,
                                    drift = rep(-0.5,N))
  
  
  scenario2D <- list()
  N <- 18
  loc <- cbind((43):(42+N), 0)
  scenario2D$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  
  scenario2D$A <- A
  
  
  scenario2D.net <- matrixtoGNAR(scenario2D$A)
  
  scenario2D.nts <- GNARsim.custom(n=220, net=scenario2D.net, initsignal = rnorm(N,0,1), alphaParams = list(rep(0.55,N)), 
                                    betaParams = list(c(0.1)), sigma=sig,
                                    drift = rep(3,N))
  
  
  
  scenario2 <- list()
  N <- 60
  loc <- cbind(1:60, 0)
  scenario2$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  A[30,31] <- 0.2
  A[31,30] <- 0.2
  A[40,41] <- 0.2
  A[41,40] <- 0.2
  A[42,43] <- 0.2
  A[43,42] <- 0.2
  
  scenario2$A <- A
  
  scenario2.net <- matrixtoGNAR(scenario2$A)
  
  scenario2.nts <- cbind(scenario2A.nts, scenario2B.nts , scenario2C.nts , scenario2D.nts)
  colnames(scenario2.nts) <- 1:60
  train_scenario2.list[[k]] <- t(scenario2.nts[1:200,])
  true_scenario2.list[[k]] <- t(scenario2.nts[201:220,])
  plot(scenario2.nts[100,], type="l", main=bquote(sigma ==.(sig)), cex.main =2, ylab="", xlab="Node index", cex.lab=1.5, ylim=c(-15, 25))
}


par(mfrow=c(3,1))
sig_list <- c(0.5, 1, 1.5)
k <- 0
train_scenario2.list <- list()
true_scenario2.list <- list()
for(sig in sig_list){
  k <- k +1
  set.seed(100*which(sig_list==sig))
  scenario2A <- list()
  N <- 30
  loc <- cbind(1:N, 0)
  scenario2A$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  
  scenario2A$A <- A
  
  
  scenario2A.net <- matrixtoGNAR(scenario2A$A)
  
  scenario2A.nts <- GNARsim.custom(n=220, net=scenario2A.net, initsignal = rnorm(N,0,1), 
                                    alphaParams = list(rep(0.55,N)), betaParams = list(c(0.4)), sigma=sig,
                                    drift = rep(1,N))
  
  scenario2B <- list()
  N <- 10
  loc <- cbind((31):(30+N), 0)
  scenario2B$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  
  scenario2B$A <- A
  
  
  scenario2B.net <- matrixtoGNAR(scenario2B$A)
  
  scenario2B.nts <- GNARsim.custom(n=220, net=scenario2B.net, initsignal = rnorm(N,0,1), alphaParams = list(rep(0.55,N)), 
                                    betaParams = list(c(0.1)), sigma=sig,
                                    drift = rep(3,N))
  
  scenario2C <- list()
  N <- 2
  loc <- cbind((41):(40+N), 0)
  scenario2C$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  
  scenario2C$A <- A
  
  
  scenario2C.net <- matrixtoGNAR(scenario2C$A)
  
  scenario2C.nts <- GNARsim.custom(n=220, net=scenario2C.net, initsignal = rnorm(N,0,1), 
                                    alphaParams = list(rep(0.55,N)), betaParams = list(c(0.4)), sigma=sig,
                                    drift = rep(-0.5,N))
  
  
  scenario2D <- list()
  N <- 18
  loc <- cbind((43):(42+N), 0)
  scenario2D$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  
  scenario2D$A <- A
  
  
  scenario2D.net <- matrixtoGNAR(scenario2D$A)
  
  scenario2D.nts <- GNARsim.custom(n=220, net=scenario2D.net, initsignal = rnorm(N,0,1), alphaParams = list(rep(0.55,N)), 
                                    betaParams = list(c(0.1)), sigma=sig,
                                    drift = rep(3,N))
  
  
  scenario2 <- list()
  N <- 60
  loc <- cbind(1:60, 0)
  scenario2$xy <- loc
  
  # weight matrix
  A <- matrix(0, nrow=N, ncol=N)
  for(i in 1:(N-1)){
    A[i,i+1] <- 1
    A[i+1,i] <- 1
  }
  A[30,31] <- 0.2
  A[31,30] <- 0.2
  A[40,41] <- 0.2
  A[41,40] <- 0.2
  A[42,43] <- 0.2
  A[43,42] <- 0.2
  
  scenario2$A <- A
  
  scenario2.net <- matrixtoGNAR(scenario2$A)
  
  scenario2.nts <- cbind(scenario2A.nts, scenario2B.nts , scenario2C.nts , scenario2D.nts)
  colnames(scenario2.nts) <- 1:60
  train_scenario2.list[[k]] <- t(scenario2.nts[1:200,])
  true_scenario2.list[[k]] <- t(scenario2.nts[201:220,])
  plot(scenario2.nts[,1], type="l", main=bquote(sigma ==.(sig)), cex.main =2, ylab="", xlab="time", cex.lab=1.8, ylim=c(10, 30))
}



############ sigma = 0.5 #################
MM <-list()
plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat0 - plt0.scenario2.mat0))

MM[[1]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="black", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - Nodewise ARIMA RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat0 - plt1.scenario2.mat0))

MM[[2]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="green", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - GNARI RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat0 - plt3.scenario2.mat0))

MM[[3]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="deepskyblue3", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - GFT method RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat0 - plt3_0.scenario2.mat0))

MM[[4]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="cyan", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - GFT method with clustering RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


grid.arrange(MM[[1]], MM[[2]], MM[[3]], MM[[4]], nrow=4)



############ sigma = 1 #################
MM <-list()
plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat1 - plt0.scenario2.mat1))

MM[[1]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="black", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - Nodewise ARIMA RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat1 - plt1.scenario2.mat1))

MM[[2]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="green", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - GNARI RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat1 - plt3.scenario2.mat1))

MM[[3]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="deepskyblue3", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - GFT method RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())

plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat1 - plt3_0.scenario2.mat1))

MM[[4]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="cyan", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - GFT method with clustering RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


grid.arrange(MM[[1]], MM[[2]], MM[[3]], MM[[4]], nrow=4)


############ sigma = 2 #################
MM <-list()
plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat2 - plt0.scenario2.mat2))

MM[[1]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="black", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - Nodewise ARIMA RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat2 - plt1.scenario2.mat2))

MM[[2]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="green", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - GNARI RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat2 - plt3.scenario2.mat2))

MM[[3]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="deepskyblue3", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - GFT method RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())

plotdf <- data.frame(time = as.factor(rep(1:h2, each=20)), 
                     val = c(plt2.scenario2.mat2 - plt3_0.scenario2.mat2))

MM[[4]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="cyan", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 2: Proposed method RMSE - GFT method with clustering RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


grid.arrange(MM[[1]], MM[[2]], MM[[3]], MM[[4]], nrow=4)




#################################################
## Scenario3 : Network Time Series with Bumps
#################################################

sig_list <- c(0.5, 1, 1.5)
for(sig in sig_list){
  plt0.scenario3.mat <- c()
  plt1.scenario3.mat <- c()
  plt2.scenario3.mat <- c()
  plt3.scenario3.mat <- c()
  
  plt0.scenario3.list <- list()
  plt1.scenario3.list <- list()
  plt2.scenario3.list <- list()
  plt3.scenario3.list <- list()
  
  set.seed(1)
  er2 <- erdos.renyi.game(n=40, p=0.1, type="gnp")

  set.seed(100*which(sig_list==sig))
  for(k in 1:20){
    print(paste("iteration:", k))
    scenario3.net <- igraphtoGNAR(er2)
    
    scenario3.nts <- GNARsim.custom(n=220, net=scenario3.net, initsignal = rnorm(40,0,1), 
                                    alphaParams = list(rep(0.6,40)), betaParams = list(c(0.3)), sigma=sig, 
                                    drift = rep(2,40))
    
    scenario3.nts[101:110,6] <- scenario3.nts[101:110,6] + 3*c(c(1:5), c(4:0))
    scenario3.nts[161:170,6] <- scenario3.nts[161:170,6] + 3*c(c(1:5), c(4:0))
    colnames(scenario3.nts) <- 1:40
    
    h3 <- 20
    N3 <- 40
    
    scenario3 <- list()
    A <- matrix(0, nrow=N3, ncol=N3)
    
    for(i in 1:length(scenario3.net$edges)){
      A[i,scenario3.net$edges[[i]]] <- scenario3.net$dist[[i]]
    }
    scenario3$A <- A
    
    train_scenario3 <- t(scenario3.nts[1:200,])
    true_scenario3 <- t(scenario3.nts[201:220,])
    
    pred.gnar.scenario3 <- forecast_narima0(vts = ts(t(train_scenario3)), h = h3, N = N3, 
                                            net = scenario3.net, max.alpha = 5, 
                                            max.beta = 3, globalalpha = TRUE, centering=FALSE)
    
    
    rownames(pred.gnar.scenario3) <- paste("Node",1:N3, sep="")
    colnames(pred.gnar.scenario3) <- paste("t", 1:h3, sep="")
    
    
    L3 <- diag(rowSums(scenario3$A)) - scenario3$A 
    val3 <- eigensort(L3)
    evalues.scenario3 <- val3$evalues
    evectors.scenario3 <- val3$evectors
    #- largest eigenvalue
    lmax.scenario3 <- max(evalues.scenario3)
    #- parameter that controls the scale number b <- 2
    b.scenario3 <- 3.5
    tf.scenario3 <- tight_frame(evalues.scenario3, evectors.scenario3, b=b.scenario3)
    
    # plot_filter(lmax.scenario3, b=b.scenario3)
    
    m3 <- nrow(tf.scenario3)/N3
    
    wcf.scenario3 <- c() 
    
    # b: hyperparameter
    for(i in 1:ncol(train_scenario3)){
      f <- train_scenario3[,i]
      wcf <- forward_sgwt(f, evalues.scenario3, evectors.scenario3, b=b.scenario3)
      wcf.scenario3 <- cbind(wcf.scenario3, wcf)
    }
    colnames(wcf.scenario3) <- paste("t", 1:ncol(train_scenario3))
    rownames(wcf.scenario3) <- paste(paste("Node",1:N3, sep=""), rep(1:m3, each=N3), sep="")
    
    
    pred.wc.arima.scenario3 <- forecast_sgwt.arima(wcf.scenario3, h=h3, period=1)
    
    pred.sgwt.arima.scenario3 <- matrix(0, nrow=N3, ncol=h3) 
    rownames(pred.sgwt.arima.scenario3) <- paste("Node",1:N3, sep="")
    colnames(pred.sgwt.arima.scenario3) <- paste("t", 1:h3, sep="")
    
    for(i in 1:h3){
      pred.sgwt.arima.scenario3[,i] <- inverse_sgwt(pred.wc.arima.scenario3[,i], evalues.scenario3, evectors.scenario3,  b=b.scenario3)
    }
    
    
    fc.data.scenario3 <- compute_gft(train_scenario3, evectors.scenario3)
    
    true_fc.scenario3 <- compute_gft(true_scenario3, evectors.scenario3)
    
    # fc.data: N  x  train date   (for train data)
    # predicted fourier coefficients
    pred.fc.scenario3 <- forecast_gft(fc.data.scenario3, h=h3, period=1) # N by h (predicted data for h time points)
    
    # predict
    pred.gft.scenario3 <- inverse_gft(pred.fc.scenario3, evectors.scenario3) 
    rownames(pred.gft.scenario3) <- paste("Node",1:N3, sep="")
    colnames(pred.gft.scenario3) <- paste("t", 1:h3, sep="")
    
    # nodewise ARIMA
    res <- matrix(0, nrow=N3, ncol=h3)
    for(i in 1:N3){
      model.fit <- auto.arima(ts(train_scenario3[i,], frequency=1))
      res[i,] <- as.numeric(forecast(model.fit, h=h3)$mean)
    }
    pred.nodewise.arima.scenario3 <- res
    
    
    plt0.scenario3 <- vector(length=20)
    for(i in 1:20){
      plt0.scenario3[i] <- rmse(c(true_scenario3[,i]), c(pred.nodewise.arima.scenario3[,i]))
    }
    
    plt1.scenario3 <- vector(length=20)
    for(i in 1:20){
      plt1.scenario3[i] <- rmse(c(true_scenario3[,i]), c(pred.gnar.scenario3[,i]))
    }
    
    plt2.scenario3 <- vector(length=20)
    for(i in 1:20){
      plt2.scenario3[i] <- rmse(c(true_scenario3[,i]), c(pred.sgwt.arima.scenario3[,i]))
    }
    
    plt3.scenario3 <- vector(length=20)
    for(i in 1:20){
      plt3.scenario3[i] <- rmse(c(true_scenario3[,i]), c(pred.gft.scenario3[,i]))
    }
    
    plt0.scenario3.list[[k]] <- pred.nodewise.arima.scenario3
    plt1.scenario3.list[[k]] <- pred.gnar.scenario3
    plt2.scenario3.list[[k]] <- pred.sgwt.arima.scenario3
    plt3.scenario3.list[[k]] <- pred.gft.scenario3
    
    plt0.scenario3.mat <- rbind(plt0.scenario3.mat, plt0.scenario3)
    plt1.scenario3.mat <- rbind(plt1.scenario3.mat, plt1.scenario3)
    plt2.scenario3.mat <- rbind(plt2.scenario3.mat, plt2.scenario3)
    plt3.scenario3.mat <- rbind(plt3.scenario3.mat, plt3.scenario3)
  }
  
  if(sig == sig_list[1]){
    plt0.scenario3.mat0 <- plt0.scenario3.mat
    plt1.scenario3.mat0 <- plt1.scenario3.mat
    plt2.scenario3.mat0 <- plt2.scenario3.mat
    plt3.scenario3.mat0 <- plt3.scenario3.mat
    
    plt0.scenario3.list0 <- plt0.scenario3.list
    plt1.scenario3.list0 <- plt1.scenario3.list
    plt2.scenario3.list0 <- plt2.scenario3.list
    plt3.scenario3.list0 <- plt3.scenario3.list
  }
  if(sig == sig_list[2]){
    plt0.scenario3.mat1 <- plt0.scenario3.mat
    plt1.scenario3.mat1 <- plt1.scenario3.mat
    plt2.scenario3.mat1 <- plt2.scenario3.mat
    plt3.scenario3.mat1 <- plt3.scenario3.mat
    
    plt0.scenario3.list1 <- plt0.scenario3.list
    plt1.scenario3.list1 <- plt1.scenario3.list
    plt2.scenario3.list1 <- plt2.scenario3.list
    plt3.scenario3.list1 <- plt3.scenario3.list
  }
  if(sig == sig_list[3]){
    plt0.scenario3.mat2 <- plt0.scenario3.mat
    plt1.scenario3.mat2 <- plt1.scenario3.mat
    plt2.scenario3.mat2 <- plt2.scenario3.mat
    plt3.scenario3.mat2 <- plt3.scenario3.mat
    
    plt0.scenario3.list2 <- plt0.scenario3.list
    plt1.scenario3.list2 <- plt1.scenario3.list
    plt2.scenario3.list2 <- plt2.scenario3.list
    plt3.scenario3.list2 <- plt3.scenario3.list
  }
}



par(mfrow=c(1,3))
k <- 0
sig_list <- c(0.5, 1, 1.5)
train_scenario3.list <- list()
true_scenario3.list <- list()
for(sig in sig_list){
  k <- k +1
  set.seed(1)
  scenario3.net <- igraphtoGNAR(er2)
  
  set.seed(100*which(sig_list==sig))
  scenario3.nts <- GNARsim.custom(n=220, net=scenario3.net, initsignal = rnorm(40,0,1), 
                                  alphaParams = list(rep(0.6,40)), betaParams = list(c(0.3)), sigma=sig, 
                                  drift = rep(2,40))
  
  scenario3.nts[101:110,6] <- scenario3.nts[101:110,6] + 3*c(c(1:5), c(4:0))
  scenario3.nts[161:170,6] <- scenario3.nts[161:170,6] + 3*c(c(1:5), c(4:0))
  colnames(scenario3.nts) <- 1:40
  
  h3 <- 20
  N3 <- 40
  
  train_scenario3.list[[k]] <- t(scenario3.nts[1:200,])
  true_scenario3.list[[k]] <- t(scenario3.nts[201:220,])
  plot(scenario3.nts[165,], type="l", main=bquote(sigma ==.(sig)), cex.main =2, ylab="", xlab="Node index", cex.lab=1.5, ylim=c(14, 40))
}

par(mfrow=c(3,1))
k <- 0
sig_list <- c(0.5, 1, 1.5)
train_scenario3.list <- list()
true_scenario3.list <- list()
for(sig in sig_list){
  k <- k +1
  set.seed(1)
  scenario3.net <- igraphtoGNAR(er2)
  
  set.seed(100*which(sig_list==sig))
  scenario3.nts <- GNARsim.custom(n=220, net=scenario3.net, initsignal = rnorm(40,0,1), 
                                  alphaParams = list(rep(0.6,40)), betaParams = list(c(0.3)), sigma=sig, 
                                  drift = rep(2,40))
  
  scenario3.nts[101:110,6] <- scenario3.nts[101:110,6] + 3*c(c(1:5), c(4:0))
  scenario3.nts[161:170,6] <- scenario3.nts[161:170,6] + 3*c(c(1:5), c(4:0))
  colnames(scenario3.nts) <- 1:40
  
  h3 <- 20
  N3 <- 40
  
  train_scenario3.list[[k]] <- t(scenario3.nts[1:200,])
  true_scenario3.list[[k]] <- t(scenario3.nts[201:220,])
  plot(scenario3.nts[,6], type="l", main=bquote(sigma ==.(sig)), cex.main =2, ylab="", xlab="time", cex.lab=1.8, ylim=c(14, 40))
}


############ sigma = 0.05 #################
MM <-list()
plotdf <- data.frame(time = as.factor(rep(1:h3, each=20)), 
                     val = c(plt2.scenario3.mat0 - plt0.scenario3.mat0))

MM[[1]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="black", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 4: Proposed method RMSE - Nodewise ARIMA RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h3, each=20)), 
                     val = c(plt2.scenario3.mat0 - plt1.scenario3.mat0))

MM[[2]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="green", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 4: Proposed method RMSE - GNARI RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h3, each=20)), 
                     val = c(plt2.scenario3.mat0 - plt3.scenario3.mat0))

MM[[3]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="deepskyblue3", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 4: Proposed method RMSE - GFT method RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


grid.arrange(MM[[1]], MM[[2]], MM[[3]], nrow=3)




############ sigma = 0.5 #################
MM <-list()
plotdf <- data.frame(time = as.factor(rep(1:h3, each=20)), 
                     val = c(plt2.scenario3.mat1 - plt0.scenario3.mat1))

MM[[1]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="black", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 4: Proposed method RMSE - Nodewise ARIMA RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h3, each=20)), 
                     val = c(plt2.scenario3.mat1 - plt1.scenario3.mat1))

MM[[2]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="green", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 4: Proposed method RMSE - GNARI RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h3, each=20)), 
                     val = c(plt2.scenario3.mat1 - plt3.scenario3.mat1))

MM[[3]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="deepskyblue3", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 4: Proposed method RMSE - GFT method RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


grid.arrange(MM[[1]], MM[[2]], MM[[3]], nrow=3)


############ sigma = 2 #################
MM <-list()
plotdf <- data.frame(time = as.factor(rep(1:h3, each=20)), 
                     val = c(plt2.scenario3.mat2 - plt0.scenario3.mat2))

MM[[1]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="black", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 4: Proposed method RMSE - Nodewise ARIMA RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h3, each=20)), 
                     val = c(plt2.scenario3.mat2 - plt1.scenario3.mat2))

MM[[2]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="green", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 4: Proposed method RMSE - GNARI RMSE") + scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


plotdf <- data.frame(time = as.factor(rep(1:h3, each=20)), 
                     val = c(plt2.scenario3.mat2 - plt3.scenario3.mat2))

MM[[3]] <- ggplot(plotdf, aes(x = time, y = val)) + geom_boxplot(fill="deepskyblue3", alpha=0.3)+ geom_hline(yintercept=0, col="red", linetype="dashed", size=1.5) +
  labs(title="Scenario 4: Proposed method RMSE - GFT method RMSE")+ scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(plot.title = element_text(size=20,face="bold"), axis.text = element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_blank())


grid.arrange(MM[[1]], MM[[2]], MM[[3]], nrow=3)



