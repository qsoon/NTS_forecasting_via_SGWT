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
source("method.R", chdir=TRUE)


#### Seoul covid19 data #####


###################
#### load data ####
###################

covid.seoul <- read.csv("../Data/covidseoul/covid_seoul.csv", 
                        header=TRUE, fileEncoding = "euc-kr")

covid.seoul <- covid.seoul[,c(1,2*c(1:25)+1)]

district.name_kor <- as.character(sapply(colnames(covid.seoul)[2:26],
                                         function(x) {strsplit(x, split="[.]")[[1]][1]}))
district.name_kor[22] <- "서초구"

district.name <- c("Jongno-gu", "Jung-gu", "Yongsan-gu", "Seongdong-gu",
                   "Gwangjin-gu", "Dongdaemun-gu", "Jungnang-gu", "Seongbuk-gu",
                   "Gangbuk-gu", "Dobong-gu", "Nowon-gu", "Eunpyeong-gu",
                   "Seodaemun-gu", "Mapo-gu", "Yangcheon-gu", "Gangseo-gu",
                   "Guro-gu", "Geumcheon-gu", "Yeongdeungpo-gu", "Dongjak-gu",
                   "Gwanak-gu", "Seocho-gu", "Gangnam-gu", "Songpa-gu", "Gangdong-gu")

colnames(covid.seoul) <- c("date", district.name)

covid.seoul <- covid.seoul[nrow(covid.seoul):1,]

covid.seoul$date <- as.character(covid.seoul$date)

covid.seoul <- covid.seoul[!duplicated(covid.seoul),]

covid.seoul$date[1:23] <- paste("20", covid.seoul$date[1:23], sep="")

covid.seoul$date <- as.character(sapply(covid.seoul$date, function(x) {paste(
  strsplit(x, "[.]")[[1]][1], strsplit(x, "[.]")[[1]][2],
  strsplit(x, "[.]")[[1]][3], sep="-")}))

rownames(covid.seoul) <- covid.seoul$date
covid.seoul <- t(covid.seoul)
covid.seoul <- covid.seoul[-1,]
covid.seoul <- apply(covid.seoul, 2, as.numeric)
rownames(covid.seoul) <- district.name

covid.seoul <- covid.seoul[,636:965]


# include gyeonggi data
covid.gyeonggi <- read.xlsx("../Data/covidseoul/covid_gyeonggi.xlsx", 1,
                            header=TRUE)

district.name_kor2 <- paste(colnames(covid.gyeonggi)[-1], "시", sep="")
district.name_kor2[27] <- "양평군"
district.name_kor2[30] <- "가평군"
district.name_kor2[31] <- "연천군"

district.name2 <- c("Suwon-si", "Goyang-si", "Yongin-si", "Seongnam-si", "Bucheon-si",
                    "Ansan-si", "Hwaseong-si", "Namyangju-si", "Anyang-si",
                    "Pyeongtaek-si", "Uijeongbu-si", "Paju-si", "Siheung-si",
                    "Gimpo-si", "Gwangmyeong-si", "Gwangju-si", "Gunpo-si",
                    "Icheon-si", "Osan-si", "Hanam-si", "Yangju-si",
                    "Guri-si", "Anseong-si", "Pocheon-si", "Uiwang-si",
                    "Yeoju-si", "Yangpyeong-gun", "Dongducheon-si", "Gwacheon-si",
                    "Gapyeong-gun", "Yeoncheon-gun")

colnames(covid.gyeonggi) <- c("date", district.name2)

covid.gyeonggi <- covid.gyeonggi[nrow(covid.gyeonggi):1,]

covid.gyeonggi$date <- as.character(covid.gyeonggi$date)

covid.gyeonggi <- covid.gyeonggi[!duplicated(covid.gyeonggi),]

covid.gyeonggi$date <- as.character(sapply(covid.gyeonggi$date, function(x) {paste(
  strsplit(x, "[.]")[[1]][1], strsplit(x, "[.]")[[1]][2],
  strsplit(x, "[.]")[[1]][3], sep="-")}))

rownames(covid.gyeonggi) <- covid.gyeonggi$date
covid.gyeonggi <- t(covid.gyeonggi)
covid.gyeonggi <- covid.gyeonggi[-1,]
covid.gyeonggi[is.na(covid.gyeonggi)] <- 0

whitespace <- strsplit(covid.gyeonggi[1,241], "")[[1]][4]
covid.gyeonggi <- apply(covid.gyeonggi, 2, function(x) gsub(whitespace,"",x))
covid.gyeonggi <- apply(covid.gyeonggi, 2, function(x) as.numeric(gsub(",","",x)))
rownames(covid.gyeonggi) <- district.name2


# impute missing
covid.gyeonggi[,132] <- (covid.gyeonggi[,125] + covid.gyeonggi[,139]) / 2
covid.gyeonggi[,273] <- (covid.gyeonggi[,266] + covid.gyeonggi[,280]) / 2
covid.gyeonggi[,286] <- (covid.gyeonggi[,279] + covid.gyeonggi[,293]) / 2
covid.gyeonggi[,312] <- (covid.gyeonggi[,305] + covid.gyeonggi[,319]) / 2


train_date.covid <- colnames(covid.seoul[,1:242]) 
pred_date.covid <- colnames(covid.seoul[,243:256])


# load location info
seoul.district.loc <- read.csv("../Data/covidseoul/seoul_district_loc.csv", 
                               header=TRUE, fileEncoding = "euc-kr")

seoul.district.loc <- seoul.district.loc[,c(2,6,7)]

colnames(seoul.district.loc) <- c("district", "lon", "lat")

seoul.district.loc <- seoul.district.loc[
  unlist(lapply(district.name_kor, 
                function(x) which(seoul.district.loc$district %in% x))),]

gyeonggi.district.loc <- read.xlsx("../Data/covidseoul/gyeonggi_district_loc.xlsx", 1, 
                                   header=TRUE)
colnames(gyeonggi.district.loc) <- colnames(seoul.district.loc)


# aggregate seoul and gyeonggi
total.district.loc <- rbind(seoul.district.loc, gyeonggi.district.loc)
covid.total <- rbind(covid.seoul, covid.gyeonggi)


district.name.total <- c(district.name, district.name2)
N.covid <- length(district.name.total)

distmat.covid <- distm(total.district.loc[,2:3], fun = distHaversine) / 1000

A.covid <- c()
for(i in 1:(nrow(distmat.covid)-1)){
  for(j in (i+1):ncol(distmat.covid)){
    val <- distmat.covid[i,j]
    A.covid <- rbind(A.covid, c(i,j,val))
  }
}

G <- graph.data.frame(A.covid[,1:2], directed=FALSE)
E(G)$weight <- A.covid[,3]
mst <- minimum.spanning.tree(G)


# as_data_frame(mst, what="vertices")
edge.wt <- igraph::as_data_frame(mst, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)

wmat <- matrix(0, nrow=N.covid, ncol=N.covid)

colnames(wmat) <- district.name.total  
rownames(wmat) <- district.name.total

for(i in 1:nrow(edge.wt)){
  wmat[edge.wt[i,1], edge.wt[i,2]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
  wmat[edge.wt[i,2], edge.wt[i,1]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
}  

L.covid <- diag(rowSums(wmat)) - wmat 
val2 <- eigensort(L.covid)
evalues.covid <- val2$evalues
evectors.covid <- val2$evectors
#- largest eigenvalue
lmax.covid <- max(evalues.covid)
#- parameter that controls the scale number b <- 2
b.covid <- 2
tf.covid <- tight_frame(evalues.covid, evectors.covid, b=b.covid)

plot_filter(lmax.covid, b=b.covid)


N.covid <- nrow(L.covid)
m.covid <- nrow(tf.covid)/N.covid

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat[edge.wt[i,1], edge.wt[i,2]]))
}

h.covid <- length(pred_date.covid)

train_data.covid <- covid.total[,1:length(train_date.covid)]
true_data.covid <- covid.total[,pred_date.covid]



##############################
#### make covid data ####
##############################
covidseoul <- list()

# location
covidseoul$xy <- as.matrix(total.district.loc[,2:3])
rownames(covidseoul$xy) <- district.name.total

# weight matrix
covidseoul$A <- wmat

# sparse weight matrix
covidseoul$sA <- sp.wmat

covidseoul$dist <- distmat.covid
covidseoul$sdist <- A.covid

covid.GNARnet <- matrixtoGNAR(covidseoul$A)

#######################
#### visualization ####
#######################
plot_graph_custom2 <- function (z, e.size=1, v.size=3, vertex_color=NULL) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                  xend = y1, yend = y2, color=w), size=e.size, 
                                              data = df2) +
    scale_color_gradient(low="grey", high="black", name="Edge Weight")+
    geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
    scale_fill_gradient(low="yellow", high="red", na.value = "yellow", name = "Confirmed cases") +
    geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
    theme_void() +
    theme(legend.margin = margin(10,10,10,10), plot.margin = margin(10,10,10,10))
  print(p1)
}


plot_graph(covidseoul)

plot_graph_custom2(covidseoul, e.size=1.3, v.size=6, vertex_color = train_data.covid[,1])




#####################################
###  USE log & DIFFERENCED DATA   ###
#####################################
train_data.covid.log <- log(1+train_data.covid)
true_data.covid.log <- log(1+true_data.covid)

train_data.covid.log.diff <- t(diff(t(train_data.covid.log), lag = 7, differences=1))


#############################
## Forecast via SGWT
#############################
wcf.data.covid.log.diff <- c() 

# b: hyperparameter
for(i in train_date.covid[8:length(train_date.covid)]){
  f <- train_data.covid.log.diff[,i]
  wcf <- forward_sgwt(f, evalues.covid, evectors.covid, b=b.covid)
  wcf.data.covid.log.diff <- cbind(wcf.data.covid.log.diff, wcf)
}
colnames(wcf.data.covid.log.diff) <- train_date.covid[8:length(train_date.covid)]
rownames(wcf.data.covid.log.diff) <- paste(district.name.total, rep(1:m.covid, each=N.covid), sep="")


pred.wc.arima.covid.log.diff <- forecast_sgwt.arima(wcf.data.covid.log.diff, h=h.covid, period=1) # m*N by h (predicted data for h time points)

pred.sgwt.arima.covid.log.diff <- matrix(0, nrow=N.covid, ncol=h.covid)
rownames(pred.sgwt.arima.covid.log.diff) <- district.name.total
colnames(pred.sgwt.arima.covid.log.diff) <- pred_date.covid

for(i in 1:h.covid){
  pred.sgwt.arima.covid.log.diff[,i] <- inverse_sgwt(pred.wc.arima.covid.log.diff[,i], evalues.covid, evectors.covid, b=b.covid)
}

pred.sgwt.arima.covid.log.diff <- t(diffinv(t(pred.sgwt.arima.covid.log.diff), lag=7, differences=1,
                                            xi=t(train_data.covid.log)[(ncol(train_data.covid.log)-6):ncol(train_data.covid.log),])[8:(7+h.covid),])


#############################
## Forecast via GFT
#############################
fc.data.covid.log.diff <- compute_gft(train_data.covid.log.diff, evectors.covid)

# fc.data: N  x  train date   (for train data)
# predicted fourier coefficients
pred.fc.covid.log.diff <- forecast_gft(fc.data.covid.log.diff, h=h.covid, period=1) # N by h (predicted data for h time points)

# predict
pred.gft.covid.log.diff <- inverse_gft(pred.fc.covid.log.diff, evectors.covid) 
rownames(pred.gft.covid.log.diff) <- district.name.total
colnames(pred.gft.covid.log.diff) <- pred_date.covid

pred.gft.covid.log.diff <- t(diffinv(t(pred.gft.covid.log.diff), lag=7, differences=1,
                                     xi=t(train_data.covid.log)[(ncol(train_data.covid.log)-6):ncol(train_data.covid.log),])[8:(7+h.covid),])



#############################
## GNARI
#############################
pred.gnar.covid.log.diff0 <- forecast_narima0(vts = ts(t(train_data.covid.log.diff)), h = h.covid, N = N.covid, 
                                              net = covid.GNARnet, max.alpha = 5, 
                                              max.beta = 3, globalalpha = TRUE, centering=FALSE)
rownames(pred.gnar.covid.log.diff0) <- district.name.total
colnames(pred.gnar.covid.log.diff0) <- pred_date.covid


pred.gnar.covid.log.diff0 <- t(diffinv(t(pred.gnar.covid.log.diff0), lag=7, differences=1,
                                       xi=t(train_data.covid.log)[(ncol(train_data.covid.log)-6):ncol(train_data.covid.log),])[8:(7+h.covid),])


#############################
## GNARI + LOCAAT
#############################
train_data.covid.log.diff.spdiff <- c()
for(i in train_date.covid[8:length(train_date.covid)]){
  tmp2 <- rep(0, N.covid)
  res2 <- LOCAAT(train_data.covid.log.diff[,i], covidseoul, stop = 3, given=FALSE)
  tmp2[res2$D] <- res2$d
  tmp2[res2$S] <- res2$c
  train_data.covid.log.diff.spdiff <- cbind(train_data.covid.log.diff.spdiff, tmp2) 
} 

colnames(train_data.covid.log.diff.spdiff) <- train_date.covid[8:length(train_date.covid)]


pred.gnar.covid.log.diff.spdiff <- forecast_narima0(vts = t(train_data.covid.log.diff.spdiff), h = h.covid, N = N.covid, 
                                                    net = covid.GNARnet, max.alpha = 5, 
                                                    max.beta = 3, globalalpha = TRUE, centering=FALSE)


info.covid.spdiff <- LOCAAT(train_data.covid.log.diff.spdiff[,1], covidseoul, stop = 3, given=FALSE)

for(i in 1:ncol(pred.gnar.covid.log.diff.spdiff)){
  pred.gnar.covid.log.diff.spdiff[,i] <- inverse_LOCAAT(pred.gnar.covid.log.diff.spdiff[,i], info.covid.spdiff, stop=3)
}


pred.gnar.covid.log.diff.spdiff <- t(diffinv(t(pred.gnar.covid.log.diff.spdiff), lag=7, differences=1,
                                             xi=t(train_data.covid.log)[(ncol(train_data.covid.log)-6):ncol(train_data.covid.log),])[8:(7+h.covid),])

rownames(pred.gnar.covid.log.diff.spdiff) <- district.name.total
colnames(pred.gnar.covid.log.diff.spdiff) <- pred_date.covid


#############################
## LOCAAT-based
#############################
# Nunes et al. (2015)
pred.covid.log.diff.spdiff <- matrix(0, nrow=N.covid, ncol=h.covid)
for(i in info.covid.spdiff$D){
  model.fit <- auto.arima(ts(train_data.covid.log.diff.spdiff[i,], frequency=1))
  pred.covid.log.diff.spdiff[i,] <- as.numeric(forecast(model.fit, h=h.covid)$mean)
}

for(i in info.covid.spdiff$S){
  df = data.frame(val = train_data.covid.log.diff.spdiff[i,], day=1:length(train_data.covid.log.diff.spdiff[i,]))
  mod.smsp <- smooth.spline(df$day, df$val, nknots = 5)
  pred.covid.log.diff.spdiff[i,] <- predict(mod.smsp, data.frame(day = (length(df$day)+1):(length(df$day)+h.covid)))$y$day
}

for(i in 1:ncol(pred.covid.log.diff.spdiff)){
  pred.covid.log.diff.spdiff[,i] <- inverse_LOCAAT(pred.covid.log.diff.spdiff[,i], info.covid.spdiff, stop=3)
}


pred.covid.log.diff.spdiff <- t(diffinv(t(pred.covid.log.diff.spdiff), lag=7, differences=1,
                                        xi=t(train_data.covid.log)[(ncol(train_data.covid.log)-6):ncol(train_data.covid.log),])[8:(7+h.covid),])

rownames(pred.covid.log.diff.spdiff) <- district.name.total
colnames(pred.covid.log.diff.spdiff) <- pred_date.covid


#############################
## Nodewise-ARIMA
#############################
res <- matrix(0, nrow=N.covid, ncol=h.covid)
for(i in 1:N.covid){
  model.fit <- auto.arima(ts(train_data.covid.log.diff[i,], frequency=1))
  res[i,] <- as.numeric(forecast(model.fit, h=h.covid)$mean)
}

pred.nodewise.arima.covid.log.diff <- t(diffinv(t(res), lag=7, differences=1,
                                                xi=t(train_data.covid.log)[(ncol(train_data.covid.log)-6):ncol(train_data.covid.log),])[8:(7+h.covid),])




#####################
## results
#####################
plt1.covid.rmse <- vector(length=14)
for(i in 1:14){
  plt1.covid.rmse[i] <- rmse(c(true_data.covid.log[,i]), c(pred.nodewise.arima.covid.log.diff[,i]))
}

plt2.covid.rmse <- vector(length=14)
for(i in 1:14){
  plt2.covid.rmse[i] <- rmse(c(true_data.covid.log[,i]), c(pred.covid.log.diff.spdiff[,i]))
}

plt3.covid.rmse <- vector(length=14)
for(i in 1:14){
  plt3.covid.rmse[i] <- rmse(c(true_data.covid.log[,i]), c(pred.gnar.covid.log.diff0[,i]))
}

plt4.covid.rmse <- vector(length=14)
for(i in 1:14){
  plt4.covid.rmse[i] <- rmse(c(true_data.covid.log[,i]), c(pred.gnar.covid.log.diff.spdiff[,i]))
}

plt5.covid.rmse <- vector(length=14)
for(i in 1:14){
  plt5.covid.rmse[i] <- rmse(c(true_data.covid.log[,i]), c(pred.gft.covid.log.diff[,i]))
}

plt6.covid.rmse <- vector(length=14)
for(i in 1:14){
  plt6.covid.rmse[i] <- rmse(c(true_data.covid.log[,i]), c(pred.sgwt.arima.covid.log.diff[,i]))
}



plt1.covid.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.covid.log[,j]) - c(pred.nodewise.arima.covid.log.diff[,j]))^2)
  }
  plt1.covid.crmse[i] <- sqrt(tmp/N.covid)
}


plt2.covid.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.covid.log[,j]) - c(pred.covid.log.diff.spdiff[,j]))^2)
  }
  plt2.covid.crmse[i] <- sqrt(tmp/N.covid)
}

plt3.covid.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.covid.log[,j]) - c(pred.gnar.covid.log.diff0[,j]))^2)
  }
  plt3.covid.crmse[i] <- sqrt(tmp/N.covid)
}


plt4.covid.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.covid.log[,j]) - c(pred.gnar.covid.log.diff.spdiff[,j]))^2)
  }
  plt4.covid.crmse[i] <- sqrt(tmp/N.covid)
}


plt5.covid.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.covid.log[,j]) - c(pred.gft.covid.log.diff[,j]))^2)
  }
  plt5.covid.crmse[i] <- sqrt(tmp/N.covid)
}

plt6.covid.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.covid.log[,j]) - c(pred.sgwt.arima.covid.log.diff[,j]))^2)
  }
  plt6.covid.crmse[i] <- sqrt(tmp/N.covid)
}



par(mfrow=c(1,2),
    mar=c(5, 4, 4, 2) + 0.1,
    oma=c(4, 0, 0, 0))

plot(plt1.covid.rmse, col="black", pch=16,
     ylim = c(min(plt1.covid.rmse, plt2.covid.rmse, plt3.covid.rmse, plt4.covid.rmse, plt5.covid.rmse, plt6.covid.rmse), 
              max(plt1.covid.rmse, plt2.covid.rmse, plt3.covid.rmse, plt4.covid.rmse, plt5.covid.rmse, plt6.covid.rmse)), 
     main = "", ylab="RMSE", xlab = "day", cex.lab=1.3)
lines(1:14, plt1.covid.rmse, col="black")
lines(1:14, plt2.covid.rmse, col="magenta")
lines(1:14, plt3.covid.rmse, col="green")
lines(1:14, plt4.covid.rmse, col="darkorange1")
lines(1:14, plt5.covid.rmse, col="blue")
lines(1:14, plt6.covid.rmse, col="red")

points(1:14, plt2.covid.rmse, col="magenta", pch=16)
points(1:14, plt3.covid.rmse, col="green", pch=16)
points(1:14, plt4.covid.rmse, col="darkorange1", pch=16)
points(1:14, plt5.covid.rmse, col="blue", pch=16)
points(1:14, plt6.covid.rmse, col="red", pch=16)


plot(plt1.covid.crmse, col="black", pch=16,
     ylim = c(min(plt1.covid.crmse, plt2.covid.crmse, plt3.covid.crmse, plt4.covid.crmse, plt5.covid.crmse, plt6.covid.crmse), 
              max(plt1.covid.crmse, plt2.covid.crmse, plt3.covid.crmse, plt4.covid.crmse, plt5.covid.crmse, plt6.covid.crmse)), 
     main = "", ylab="cRMSE", xlab = "day", cex.lab=1.3)
lines(1:14, plt1.covid.crmse, col="black")
lines(1:14, plt2.covid.crmse, col="magenta")
lines(1:14, plt3.covid.crmse, col="green")
lines(1:14, plt4.covid.crmse, col="darkorange1")
lines(1:14, plt5.covid.crmse, col="blue")
lines(1:14, plt6.covid.crmse, col="red")

points(1:14, plt2.covid.crmse, col="magenta", pch=16)
points(1:14, plt3.covid.crmse, col="green", pch=16)
points(1:14, plt4.covid.crmse, col="darkorange1", pch=16)
points(1:14, plt5.covid.crmse, col="blue", pch=16)
points(1:14, plt6.covid.crmse, col="red", pch=16)

par(fig=c(0,1,0,1), oma = c(0, 0, 0, 0), mar = c(0, 5, 0, 0), new = TRUE)
# plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
plot.new()
legend("bottom", legend=c("Nodewise ARIMA", "LOCAAT", "GNARI", "GNARI + LOCAAT", "GFT", "SGWT"), 
       col=c("black", "magenta", "green", "darkorange1", "blue", "red"), pch=16, horiz=TRUE, lty=1, 
       lwd=2, pt.cex=1.2, cex=1, bty="n", x.intersp=1, text.width=0.11)