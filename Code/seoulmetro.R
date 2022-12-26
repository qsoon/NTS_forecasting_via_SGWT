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


#### Seoul metro data #####


###################
#### load data ####
###################

# station Korean, English name
station.name <- read.csv("/home/kyu9510/NTS_forecasting_via_SGWT/Data/station_name.csv", 
                         header=TRUE, fileEncoding = "euc-kr")
colnames(station.name) <- c("code", "name_kor", "name", "line", "external_code")


# station latitude, longitude info
station.loc <- read.csv("/home/kyu9510/NTS_forecasting_via_SGWT/Data/station_location.csv", 
                        header=TRUE, fileEncoding = "euc-kr")
colnames(station.loc) <- c("ID", "name_kor", "line", "lon", "lat")

station.loc2 <- read.csv("/home/kyu9510/NTS_forecasting_via_SGWT/Data/station_location2.csv", 
                         header=TRUE, fileEncoding = "euc-kr")
station.loc2 <- station.loc2[,c(2:6)]
colnames(station.loc2) <- c("line", "station_num", "name_kor", "lon", "lat")


# distance between stations
station.distance <- read.csv("/home/kyu9510/NTS_forecasting_via_SGWT/Data/station_distance.csv", 
                             header=TRUE, fileEncoding = "euc-kr")

station.distance <- station.distance[, c(2:5)]
colnames(station.distance) <- c("line", "name_kor", "btwdist", "cumdist")

# 2021 hourly getting on/off info for each station
hourly.pplnum.2021 <- read.csv("/home/kyu9510/NTS_forecasting_via_SGWT/Data/hourly_pplnum_2021.csv", 
                               header=TRUE, fileEncoding = "euc-kr")
hourly.pplnum.2021 <- hourly.pplnum.2021[, -1]
colnames(hourly.pplnum.2021) <- c("date", "line", "station_num", "name_kor", "type",
                                  paste(rep("t",19), c(1:19), sep=""))


# We have getting on/off info only for line 1~8, target line = line 1~8
target_line <- c("01호선", "02호선", "03호선", "04호선",
                 "05호선", "06호선", "07호선", "08호선",
                 "1호선", "2호선", "3호선", "4호선",
                 "5호선", "6호선", "7호선", "8호선")


############################
#### Data preprocessing ####
############################

station.name$name_kor <- as.character(station.name$name_kor)
station.name <- station.name[station.name$name_kor != "이수",]
station.name$name_kor[which(station.name$name_kor == "4?19민주묘지")] <- "4.19민주묘지"
station.name <- station.name[station.name$line %in% target_line,]
station.name$line <- as.character(station.name$line)
station.name$line <- as.character(sapply(station.name$line,
                                         function(x) {strsplit(x, split="")[[1]][2]}))

station.loc$name_kor <- as.character(sapply(as.character(station.loc$name_kor), 
                                            function(x) {strsplit(x, split="\\(")[[1]][1]}))
station.loc <- station.loc[station.loc$line %in% target_line,]
station.loc <- as.data.frame(dplyr::select(station.loc, name_kor, lon, lat))

station.loc$name_kor[which(station.loc$name_kor=="이수")] <- "총신대입구"
# averaging location if there are several line passing the same station
station.loc <- as.data.frame(station.loc %>% group_by(name_kor) %>% 
                               summarise(lat = mean(lat), lon = mean(lon)))

# included in station.loc. lack of info 
station.loc2$name_kor <- as.character(station.loc2$name_kor)
station.loc2$name_kor[which(station.loc2$name_kor=="서울")] <- "서울역"

station.distance$name_kor <- as.character(station.distance$name_kor)
station.distance$name_kor[which(station.distance$name_kor=="이수")] <- "총신대입구"

# remove "(", ")"  in the station name
hourly.pplnum.2021$name_kor <- as.character(sapply(as.character(hourly.pplnum.2021$name_kor), 
                                                   function(x) {strsplit(x, split="\\(")[[1]][1]}))
hourly.pplnum.2021$name_kor[which(hourly.pplnum.2021$name_kor=="이수")] <- "총신대입구"
hourly.pplnum.2021 <- dplyr::select(hourly.pplnum.2021,-c("line", "station_num", "type"))
hourly.pplnum.2021 <- as.data.frame(hourly.pplnum.2021)
hourly.pplnum.2021[,3:ncol(hourly.pplnum.2021)] <- 
  sapply(hourly.pplnum.2021[,3:ncol(hourly.pplnum.2021)], as.numeric)

# summing the numbers of people across type, line
hourly.pplnum.2021 <- as.data.frame(hourly.pplnum.2021 %>% group_by(date, name_kor) %>% 
                                      summarise(across(everything(),sum)))


hourly.pplnum.2021 %>% group_by(date) %>% summarise(count = n())

# stations with no info
station_removed <- c("신내", "강일", "하남검단산", "하남시청")

hourly.pplnum.2021 <- hourly.pplnum.2021[!(hourly.pplnum.2021$name_kor %in% station_removed), ]

hourly.pplnum.2021 %>% group_by(date) %>% summarise(count = n())

########################
#### aggregate info ####
########################

# add location to getting on/off info
station.info <- inner_join(hourly.pplnum.2021, station.loc, by='name_kor')

# add station's English name
station.info <- inner_join(station.info, station.name[, c("name_kor", "name")][
  !duplicated(station.name[, c("name_kor", "name")]),], by='name_kor')

station.info$name <- as.character(station.info$name)

# column reordering
station.info <- dplyr::select(station.info, 1,2, ncol(station.info):(ncol(station.info)-2),
                              3:(ncol(station.info)-3))

# our target stations
target_station <- unique(station.info$name) # total 244 stations
target_station_kor <- unique(station.info$name_kor)

# fill NA values
tmp <- target_station
tmp.kor <- target_station_kor
tmp.loc <- sapply(station.info[!duplicated(station.info[,c("lon","lat")]),
                               c("lon","lat")], as.numeric)

datenum <- length(unique(station.info$date))

tmp.mat <- as.data.frame(cbind(rep(tmp.kor, datenum), rep(tmp, datenum),
                               rep(tmp.loc[,1], datenum), rep(tmp.loc[,2], datenum)))

colnames(tmp.mat) <- c("name_kor", "name", "lon", "lat")
tmp.mat$name_kor <- as.character(tmp.mat$name_kor)
tmp.mat$name <- as.character(tmp.mat$name)
tmp.mat$lon <- as.numeric(as.character(tmp.mat$lon))
tmp.mat$lat <- as.numeric(as.character(tmp.mat$lat))
tmp.mat$date <- as.factor(as.character(rep(seq(as.Date("2021-01-01"), by = "day", length.out = datenum),
                                           each = length(target_station))))

tmp2 <- left_join(tmp.mat, as.data.frame(station.info), 
                  by=c("date"="date","name_kor"="name_kor", 
                       "name"="name"))


station.info <- dplyr::select(tmp2, 5,1,2,3,4,8:ncol(tmp2))
colnames(station.info)[4:5] <- c("lon", "lat")

station.info$total <- rowSums(station.info[,c(6:24)])

#########################################
#### edge weight matrix construction ####
#########################################

station.distance <- as.data.frame(station.distance) 
# remove "(", ")" in the name
station.distance$name_kor <- as.character(sapply(as.character(station.distance$name_kor), 
                                                 function(x) {strsplit(x, split="\\(")[[1]][1]}))
# remove blank at the end of name
station.distance$name_kor <- as.character(sapply(as.character(station.distance$name_kor), 
                                                 function(x) {strsplit(x, split=" ")[[1]][1]}))
station.distance$name_kor[which(station.distance$name_kor=="신내역")] <- "신내"
# add station's English name
station.distance <- station.distance[station.distance$name_kor %in% target_station_kor,]

station.distance <- inner_join(station.distance, 
                               station.info[!duplicated(station.info[,c("name_kor","name")]),
                                            c("name_kor", "name")], by='name_kor')

station.distance[station.distance$name_kor=="산성",]$btwdist <-
  station.distance[station.distance$name_kor=="산성",]$btwdist + 1.5 # 남위례역이 hourly data 정보에 없어서 빠지기 때문에 간격 더해줌


station.distance[station.distance$name_kor=="미사",]$btwdist <-
  station.distance[station.distance$name_kor=="미사",]$btwdist + 0.8 # 강일이 hourly data 정보에 없어서 빠지기 때문에 간격 더해줌

station.distance$name <- as.character(station.distance$name)

e.weight <- matrix(0,nrow=length(target_station), ncol=length(target_station))
colnames(e.weight) <- target_station
rownames(e.weight) <- target_station

for(i in 1:8){
  tmp <-station.distance[station.distance$line==i,]
  if(i==2){
    n <- 44 # circular line. 54th line : City Hall again
    e.weight["Seongsu", "Yongdap"] <- tmp$btwdist[45]
    e.weight["Yongdap", "Seongsu"] <- tmp$btwdist[45]
    for(j in (n+1):47){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
    e.weight["Sindorim", "Dorimcheon"] <- tmp$btwdist[49]
    e.weight["Dorimcheon", "Sindorim"] <- tmp$btwdist[49]
    for(j in 49:51){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
  } else if(i==5){
    n <- 46 # Hanam Pungsan is the end station of one line in Line5
    e.weight["Gangdong", "Dunchon-dong"] <- tmp$btwdist[47]
    e.weight["Dunchon-dong", "Gangdong"] <- tmp$btwdist[47]
    for(j in (n+1):52){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
  } else{
    n <- nrow(tmp)
  }
  for(j in 1:(n-1)){
    e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
    e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
  }
}

e.weight.old <- e.weight
e.weight.old[e.weight.old!=0] <- exp(-e.weight.old[e.weight.old!=0])
e.weight[e.weight!=0] <- exp(-(e.weight[e.weight!=0])^2 / mean(e.weight[e.weight!=0])^2)

e.sp.weight <- NULL
e.color <- c() # for line color
color.cand <- c("blue", "yellowgreen", "orangered", "cyan",
                "darkorchid", "chocolate3", "darkolivegreen", "hotpink")
for(i in 1:8){
  tmp <-station.distance[station.distance$line==i,]
  if(i==2){
    n <- 44 # circular line. 54th line : City Hall again
    e.sp.weight <- rbind(e.sp.weight, c("Seongsu", "Yongdap", tmp$btwdist[45]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, c("Yongdap", "Seongsu", tmp$btwdist[45]))
    e.color <- c(e.color, color.cand[i])
    for(j in (n+1):47){
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
    e.sp.weight <- rbind(e.sp.weight, c("Sindorim", "Dorimcheon", tmp$btwdist[49]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, c("Dorimcheon", "Sindorim", tmp$btwdist[49]))
    e.color <- c(e.color, color.cand[i])
    
    for(j in 49:51){
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
  } else if(i==5){
    n <- 46 # Hanam Pungsan is the end station of one line in Line5
    e.sp.weight<- rbind(e.sp.weight, c("Gangdong", "Dunchon-dong", tmp$btwdist[47]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight<- rbind(e.sp.weight, c("Dunchon-dong", "Gangdong", tmp$btwdist[47]))
    e.color <- c(e.color, color.cand[i])
    for(j in (n+1):52){
      e.sp.weight <- rbind(e.sp.weight, c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
  } else{
    n <- nrow(tmp)
  }
  for(j in 1:(n-1)){
    e.sp.weight <- rbind(e.sp.weight, 
                         c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, 
                         c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
    e.color <- c(e.color, color.cand[i])
  }
}

tmp <- as.data.frame(e.sp.weight)
tmp$V1 <- as.character(tmp$V1)
tmp$V2 <- as.character(tmp$V2)
tmp2 <- exp(-as.numeric(as.character(tmp$V3)))
tmp$V3 <- exp(-(as.numeric(as.character(tmp$V3)))^2 / (mean(as.numeric(as.character(tmp$V3))))^2)

e.sp.weight <- tmp
e.sp.weight.old <- e.sp.weight
e.sp.weight.old[,3] <- tmp2

write.csv(station.info,file="/home/kyu9510/NTS_forecasting_via_SGWT/Data/master_station_info.csv", row.names=FALSE)
write.csv(station.distance,file="/home/kyu9510/NTS_forecasting_via_SGWT/Data/master_station_distance.csv", row.names=FALSE)


############################################################################################################
############################################################################################################
#### start from this line!!!!
############################################################################################################
############################################################################################################

station.info = read.csv("/home/kyu9510/NTS_forecasting_via_SGWT/Data/master_station_info.csv")
station.distance = read.csv("/home/kyu9510/NTS_forecasting_via_SGWT/Data/master_station_distance.csv")

##############################
#### make seoulmetro data ####
##############################
seoulmetro <- list()

# location
seoulmetro$xy <- sapply(station.info[!duplicated(station.info[,c("lon","lat")]),
                                     c("lon","lat")], as.numeric)
rownames(seoulmetro$xy) <- target_station

# weight matrix
seoulmetro$A <- e.weight

# sparse weight matrix
seoulmetro$sA <- e.sp.weight

# dist matrix
distmat <- e.weight.old
distmat[distmat!=0] <- -log(e.weight.old[e.weight.old!=0])

tmp <- e.sp.weight.old
tmp[,3] <- -log(tmp[,3])

seoulmetro$dist <- distmat
seoulmetro$sdist <- tmp

metro.GNARnet <- matrixtoGNAR(seoulmetro$A)

#######################
#### visualization ####
#######################
plot_graph_custom <- function (z, size = 0.75, edge_color, vertex_color=NULL) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$color <- factor(edge_color, levels = unique(edge_color))
  p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                  xend = y1, yend = y2, colour = color), 
                                              size = 2, data = df2) +
    scale_color_manual(values=color.cand, labels = paste("Line", 1:8, sep=""),
                       name = "Line number") + 
    geom_point(aes(fill=vertex_color), size = size, shape=21) + 
    scale_fill_gradient(low="white", high="black", na.value = "yellow", name = "People") +
    theme_void() +
    theme(legend.margin = margin(10,10,10,10))
  print(p1)
}

plot_graph_custom(seoulmetro, size = 3, edge_color = e.color)


###########################
#### spectral analysis ####
###########################

L.metro <- laplacian_mat(seoulmetro$A) # laplacian matrix
val1 <- eigensort(L.metro)
evalues.metro <- val1$evalues
evectors.metro <- val1$evectors
#- largest eigenvalue
lmax.metro <- max(evalues.metro)
#- parameter that controls the scale number b <- 2
b.metro <- 2.2 # m = 4 when b = 2, m = 6 when b = 1.4
tf.metro <- tight_frame(evalues.metro, evectors.metro, b=b.metro)

plot_filter(lmax.metro,b=b.metro)

N.metro <- nrow(L.metro)
m.metro <- nrow(tf.metro)/N.metro


#############################
#### daily data analysis ####
#############################
train_date.metro <- as.character(unique(station.info$date)) # train date 
train_date.metro <- train_date.metro[-1] # remove 2021-01-01 since it is holiday and outlier for initial value leads to serious error
train_date.metro <- train_date.metro[-c((length(train_date.metro)-13):length(train_date.metro))] 

train_data.metro <- c() # col: date, row: station

# b: hyperparameter
for(i in train_date.metro){
  f <- station.info[station.info$date == i,]$total
  train_data.metro <- cbind(train_data.metro, f)
}
colnames(train_data.metro) <- train_date.metro
rownames(train_data.metro) <- target_station

# graph signal snapshot example when day = 2021-01-02
plot_graph_custom(seoulmetro, size = 3, edge_color = e.color, 
                  vertex_color = train_data.metro[,1])


##########################################
## time series cleansing for each station
##########################################
data <- train_data.metro[209,] # time series of node 1
autoplot(ts(data)) + ylab("Time series of node 1") + xlab("day")  # has a periodicity of 7 and several outliers

# let's detect outliers by differencing with period 7 
autoplot(ts(diff(data, lag = 7, differences=1))) # lag: integer indicating which lag to use / difference: how many times will you perform differencing?
# can see where the outliers exist

which(abs(diff(data, lag = 7, differences=1)) > 5000) # negative outliers on holidays, positive outliers 7 days after holidays

# holiday leads to outliers 
# We need to treat holidays

# 1/1: New Year's Day
# 2/11,12,13: the Lunar New year
# 3/1: Independence Movement Day (March 1)
# 5/5: Children's day
# 5/19: Buddha's Birthday
# 8/16: replaced holiday of National Liberation Day
# 9/20,21,22: Chuseok
# 10/4: replaced holiday of National foundation day of Korea
# 10/11: replaced holiday of Hangul proclamation day
holidays <- c("2021-01-01", "2021-02-11", "2021-02-12", "2021-02-13",
              "2021-03-01", "2021-05-05", "2021-05-19", "2021-08-16",
              "2021-09-20", "2021-09-21", "2021-09-22", "2021-10-04",
              "2021-10-11")
holidays.ind <- which(train_date.metro %in% holidays)

# highlight holidays with red dots, which are outliers
autoplot(ts(data)) + ylab("peoplenum") + xlab("day") +
  geom_point(aes(y=ts(data))) +
  geom_point(data=data.frame(time=t[holidays.ind], val=ts(data)[holidays.ind]), aes(x=time, y=val), color="red")

# cleansing the differenced data
z1 <- diff(data, lag = 7, differences=1) # lag: integer indicating which lag to use / difference: how many times will you perform differencing?
autoplot(ts(z1)) + ylab("diffrenced data") + xlab("day") +
  geom_point(aes(y=ts(z1))) +
  geom_point(data=data.frame(time=t[holidays.ind-7], val=ts(z1)[holidays.ind-7]), aes(x=time, y=val), color="red")

tsoutliers(z1)
# 4   6  33  34(holiday)  35(holiday)  36(holiday)  41  42  43  52(holiday)  
# 59 117(holiday) 124 131(holiday) 138 185 186 187 220(holiday) 227 253
# 255(holiday) 256(holiday) 257(holiday) 258 262 263 264 265 269(holiday) (276, holiday) 283

autoplot(ts(z1)) + ylab("diffrenced data") + xlab("day") +
  geom_point(aes(y=ts(z1))) +
  geom_point(data=data.frame(time=t[tsoutliers(z1)$index], 
                             val=ts(z1)[tsoutliers(z1)$index]), aes(x=time, y=val), color="blue") +
  geom_point(data=data.frame(time=t[holidays.ind-7], val=ts(z1)[holidays.ind-7]), aes(x=time, y=val), color="red")

# cleaning z1 values on holidays and holidays+7 
autoplot(ts(z1)) + ylab("diffrenced data") + xlab("day") +
  geom_point(aes(y=ts(z1))) +
  geom_point(data=data.frame(time=t[tsoutliers(z1)$index], 
                             val=ts(z1)[tsoutliers(z1)$index]), aes(x=time, y=val), color="blue") +
  geom_point(data=data.frame(time=t[c(holidays.ind-7, holidays.ind)], 
                             val=ts(z1)[c(holidays.ind-7, holidays.ind)]), aes(x=time, y=val), color="red")


replace.ind <- sort(c(holidays.ind, holidays.ind-7))
replace.ind2 <- intersect(replace.ind[!duplicated(replace.ind)], tsoutliers(z1)$index) # 276 does not classified as outliers in tsoutliers
replace.ind3 <- which(tsoutliers(z1)$index %in% replace.ind2)
z2 <- z1
z2[replace.ind2] <- tsoutliers(z1)$replacements[replace.ind3] # cleansed z1

autoplot(ts(z2)) + ylab("cleansed diffrenced data") + xlab("day") +
  geom_point(aes(y=z2)) +
  geom_point(data=data.frame(time=t[replace.ind], val=z2[replace.ind]), aes(x=time, y=val), color="red")


# preprocess data considering holidays effect
data.cleansed <- data
data.cleansed[holidays.ind] <- (ts(diffinv(z2, lag=7, differences=1, xi=data[1:7])))[holidays.ind]
autoplot(ts(data.cleansed)) + ylab("cleansed data of node 1") + xlab("day") +
  geom_point(aes(y=ts(data.cleansed))) +
  geom_point(data=data.frame(time=t[holidays.ind], val=ts(data.cleansed)[holidays.ind]), aes(x=time, y=val), color="red")


# holiday_cleansing
holiday_cleansing <- function(signal, holidays.ind){
  signal.cleansed <- matrix(0, nrow=nrow(signal), ncol=ncol(signal))
  for(i in 1:nrow(signal)){
    data <- signal[i,]
    z1 <- diff(data, lag = 7, differences=1)
    replace.ind <- sort(c(holidays.ind, holidays.ind-7))
    replace.ind2 <- intersect(replace.ind[!duplicated(replace.ind)], tsoutliers(z1)$index) # 276 does not classified as outliers in tsoutliers
    replace.ind3 <- which(tsoutliers(z1)$index %in% replace.ind2)
    z2 <- z1
    z2[replace.ind2] <- tsoutliers(z1)$replacements[replace.ind3] # cleansed z1
    data.cleansed <- data
    data.cleansed[holidays.ind] <- (ts(diffinv(z2, lag=7, differences=1, xi=data[1:7])))[holidays.ind]
    signal.cleansed[i,] <- data.cleansed  
  }
  rownames(signal.cleansed) <- rownames(signal)
  colnames(signal.cleansed) <- colnames(signal)
  return(signal.cleansed)
}

train_data.metro.cleansed <- holiday_cleansing(train_data.metro, holidays.ind)



################################
### USE DIFFERENCED DATA     ###
################################
train_data.metro.cleansed.diff <- t(diff(t(train_data.metro.cleansed), lag = 7, differences=1))
autoplot(ts(train_data.metro.cleansed.diff[1,])) + geom_point(aes(y=train_data.metro.cleansed.diff[1,]))
apply(train_data.metro.cleansed.diff, 1, mean)

#############################
## Forecast via SGWT
#############################
wcf.data.metro.diff <- c() # col: date, row: station, scale 

# b: hyperparameter
for(i in train_date.metro[8:length(train_date.metro)]){
  f <- train_data.metro.cleansed.diff[,i]
  wcf <- forward_sgwt(f, evalues.metro, evectors.metro, b=b.metro)
  wcf.data.metro.diff <- cbind(wcf.data.metro.diff, wcf)
}
colnames(wcf.data.metro.diff) <- train_date.metro[8:length(train_date.metro)]
rownames(wcf.data.metro.diff) <- paste(target_station, rep(1:m.metro, each=N.metro), sep="")


pred.wc.arima.metro.diff <- forecast_sgwt.arima(wcf.data.metro.diff, h=h.metro, period=1) # m*N by h (predicted data for h time points)


# predict
# sgwt.arima
pred.sgwt.arima.metro.diff <- matrix(0, nrow=N.metro, ncol=h.metro) # predicted people num, row: stations, col: pred_date
rownames(pred.sgwt.arima.metro.diff) <- target_station
colnames(pred.sgwt.arima.metro.diff) <- pred_date.metro

for(i in 1:h.metro){
  pred.sgwt.arima.metro.diff[,i] <- inverse_sgwt(pred.wc.arima.metro.diff[,i], evalues.metro, evectors.metro, b=b.metro)
}

pred.sgwt.arima.metro.diff <- t(diffinv(t(pred.sgwt.arima.metro.diff), lag=7, differences=1,
                                        xi=t(train_data.metro.cleansed)[(ncol(train_data.metro.cleansed)-6):ncol(train_data.metro.cleansed),])[8:(7+h.metro),])



#############################
## Forecast via GFT
#############################
fc.data.metro.diff <- compute_gft(train_data.metro.cleansed.diff, evectors.metro)

# fc.data: N  x  train date   (for train data)
# predicted fourier coefficients
pred.fc.metro.diff <- forecast_gft(fc.data.metro.diff, h=h.metro, period=1) # N by h (predicted data for h time points)
# takes about 10 minutes

# predict
pred.gft.metro.diff <- inverse_gft(pred.fc.metro.diff, evectors.metro) # predicted people num, row: stations, col: pred_date
rownames(pred.gft.metro.diff) <- target_station
colnames(pred.gft.metro.diff) <- pred_date.metro

pred.gft.metro.diff <- t(diffinv(t(pred.gft.metro.diff), lag=7, differences=1,
                                 xi=t(train_data.metro.cleansed)[(ncol(train_data.metro.cleansed)-6):ncol(train_data.metro.cleansed),])[8:(7+h.metro),])


#############################
## GNARI
#############################

metro.GNARnet <- matrixtoGNAR(seoulmetro$A)

metroTS.diff <- ts(t(train_data.metro.cleansed.diff))

## without spatial differencing 

# we set globalalpha = TRUE b/c it shows better performance than globalalpha = FALSE 


pred.gnar.metro.diff0 <- forecast_narima0(vts = ts(t(train_data.metro.cleansed.diff)), h = h.metro, N = N.metro, 
                                          net = metro.GNARnet, max.alpha = 5, 
                                          max.beta = 3, globalalpha = TRUE, centering=FALSE)
rownames(pred.gnar.metro.diff0) <- target_station
colnames(pred.gnar.metro.diff0) <- pred_date.metro

pred.gnar.metro.diff0 <- t(diffinv(t(pred.gnar.metro.diff0), lag=7, differences=1,
                                   xi=t(train_data.metro.cleansed)[(ncol(train_data.metro.cleansed)-6):ncol(train_data.metro.cleansed),])[8:(7+h.metro),])





#############################
## GNARI + LOCAAT
#############################

train_data.metro.cleansed.diff.spdiff <- c()
for(i in train_date.metro[8:length(train_date.metro)]){
  tmp2 <- rep(0, N.metro)
  res2 <- LOCAAT(train_data.metro.cleansed.diff[,i], seoulmetro, stop = 3, given=TRUE)
  tmp2[res2$D] <- res2$d
  tmp2[res2$S] <- res2$c
  train_data.metro.cleansed.diff.spdiff <- cbind(train_data.metro.cleansed.diff.spdiff, tmp2) 
} # takes about 15 min

colnames(train_data.metro.cleansed.diff.spdiff) <- train_date.metro[8:length(train_date.metro)]

pred.gnar.metro.cleansed.diff.spdiff0 <- forecast_narima0(vts = t(train_data.metro.cleansed.diff.spdiff), h = h.metro, N = N.metro, 
                                                          net = metro.GNARnet, max.alpha = 5, 
                                                          max.beta = 3, globalalpha = TRUE, centering=FALSE)


info.metro.cleansed.diff.spdiff <- LOCAAT(train_data.metro.cleansed.diff.spdiff[,1], seoulmetro, stop = 3, given=TRUE)

for(i in 1:ncol(pred.gnar.metro.cleansed.diff.spdiff0)){
  pred.gnar.metro.cleansed.diff.spdiff0[,i] <- inverse_LOCAAT(pred.gnar.metro.cleansed.diff.spdiff0[,i], 
                                                              info.metro.cleansed.diff.spdiff, stop=3)
}


rownames(pred.gnar.metro.cleansed.diff.spdiff0) <- target_station
colnames(pred.gnar.metro.cleansed.diff.spdiff0) <- pred_date.metro

pred.gnar.metro.cleansed.diff.spdiff0 <- t(diffinv(t(pred.gnar.metro.cleansed.diff.spdiff0), lag=7, differences=1,
                                                   xi=t(train_data.metro.cleansed)[(ncol(train_data.metro.cleansed)-6):ncol(train_data.metro.cleansed),])[8:(7+h.metro),])

rownames(pred.gnar.metro.cleansed.diff.spdiff0) <- target_station
colnames(pred.gnar.metro.cleansed.diff.spdiff0) <- pred_date.metro



############################
## LOCAAT-based
############################

# Nunes et al. (2015)
pred.metro.cleansed.diff.spdiff <- matrix(0, nrow=N.metro, ncol=h.metro)
for(i in info.metro.cleansed.diff.spdiff$D){
  model.fit <- auto.arima(ts(train_data.metro.cleansed.diff.spdiff[i,], frequency=1))
  pred.metro.cleansed.diff.spdiff[i,] <- as.numeric(forecast(model.fit, h=h.metro)$mean)
}

for(i in info.metro.cleansed.diff.spdiff$S){
  df <- data.frame(val = train_data.metro.cleansed.diff.spdiff[i,], day=1:length(train_data.metro.cleansed.diff.spdiff[i,]))
  mod.smsp <- smooth.spline(df$day, df$val, nknots = 5)
  pred.metro.cleansed.diff.spdiff[i,] <- predict(mod.smsp, data.frame(day = (length(df$day)+1):(length(df$day)+h.metro)))$y$day
}

for(i in 1:ncol(pred.metro.cleansed.diff.spdiff)){
  pred.metro.cleansed.diff.spdiff[,i] <- inverse_LOCAAT(pred.metro.cleansed.diff.spdiff[,i], info.metro.cleansed.diff.spdiff, stop=3)
}


pred.metro.cleansed.diff.spdiff <- t(diffinv(t(pred.metro.cleansed.diff.spdiff), lag=7, differences=1,
                                             xi=t(train_data.metro.cleansed)[(ncol(train_data.metro.cleansed)-6):ncol(train_data.metro.cleansed),])[8:(7+h.metro),])

rownames(pred.metro.cleansed.diff.spdiff) <- target_station
colnames(pred.metro.cleansed.diff.spdiff) <- pred_date.metro


##########################
## Node-wise ARIMA
##########################

# nodewise ARIMA
res <- matrix(0, nrow=N.metro, ncol=h.metro)
for(i in 1:N.metro){
  model.fit <- auto.arima(ts(train_data.metro.cleansed.diff[i,], frequency=1))
  res[i,] <- as.numeric(forecast(model.fit, h=h.metro)$mean)
}

pred.nodewise.arima.metro.cleansed.diff <- t(diffinv(t(res), lag=7, differences=1,
                                                     xi=t(train_data.metro.cleansed)[(ncol(train_data.metro.cleansed)-6):ncol(train_data.metro.cleansed),])[8:(7+h.metro),])



#####################
## results
#####################

plt1.metro.rmse <- vector(length=14)
for(i in 1:14){
  plt1.metro.rmse[i] <- rmse(c(true_data.metro[,i]), c(pred.nodewise.arima.metro.cleansed.diff[,i]))
}

plt2.metro.rmse <- vector(length=14)
for(i in 1:14){
  plt2.metro.rmse[i] <- rmse(c(true_data.metro[,i]), c(pred.metro.cleansed.diff.spdiff[,i]))
}

plt3.metro.rmse <- vector(length=14)
for(i in 1:14){
  plt3.metro.rmse[i] <- rmse(c(true_data.metro[,i]), c(pred.gnar.metro.diff0[,i]))
}

plt4.metro.rmse <- vector(length=14)
for(i in 1:14){
  plt4.metro.rmse[i] <- rmse(c(true_data.metro[,i]), c(pred.gnar.metro.cleansed.diff.spdiff[,i]))
}

plt5.metro.rmse <- vector(length=14)
for(i in 1:14){
  plt5.metro.rmse[i] <- rmse(c(true_data.metro[,i]), c(pred.gft.metro.diff[,i]))
}

plt6.metro.rmse <- vector(length=14)
for(i in 1:14){
  plt6.metro.rmse[i] <- rmse(c(true_data.metro[,i]), c(pred.sgwt.arima.metro.diff[,i]))
}



plt1.metro.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.metro[,j]) - c(pred.nodewise.arima.metro.cleansed.diff[,j]))^2)
  }
  plt1.metro.crmse[i] <- sqrt(tmp/N.metro)
}


plt2.metro.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.metro[,j]) - c(pred.metro.cleansed.diff.spdiff[,j]))^2)
  }
  plt2.metro.crmse[i] <- sqrt(tmp/N.metro)
}

plt3.metro.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.metro[,j]) - c(pred.gnar.metro.diff0[,j]))^2)
  }
  plt3.metro.crmse[i] <- sqrt(tmp/N.metro)
}


plt4.metro.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.metro[,j]) - c(pred.gnar.metro.cleansed.diff.spdiff[,j]))^2)
  }
  plt4.metro.crmse[i] <- sqrt(tmp/N.metro)
}


plt5.metro.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.metro[,j]) - c(pred.gft.metro.diff[,j]))^2)
  }
  plt5.metro.crmse[i] <- sqrt(tmp/N.metro)
}

plt6.metro.crmse <- vector(length=14)
for(i in 1:14){
  tmp <- 0
  for(j in 1:i){
    tmp <- tmp + sum((c(true_data.metro[,j]) - c(pred.sgwt.arima.metro.diff[,j]))^2)
  }
  plt6.metro.crmse[i] <- sqrt(tmp/N.metro)
}



par(mfrow=c(1,2),
    mar=c(5, 4, 4, 2) + 0.1,
    oma=c(4, 0, 0, 0))

plot(plt1.metro.rmse, col="black", pch=16,
     ylim = c(min(plt1.metro.rmse, plt2.metro.rmse, plt3.metro.rmse, plt4.metro.rmse, plt5.metro.rmse, plt6.metro.rmse), 
              max(plt1.metro.rmse, plt2.metro.rmse, plt3.metro.rmse, plt4.metro.rmse, plt5.metro.rmse, plt6.metro.rmse)), 
     main = "", ylab="RMSE", xlab = "day", cex.lab=1.3)
lines(1:14, plt1.metro.rmse, col="black")
lines(1:14, plt2.metro.rmse, col="magenta")
lines(1:14, plt3.metro.rmse, col="green")
lines(1:14, plt4.metro.rmse, col="darkorange1")
lines(1:14, plt5.metro.rmse, col="blue")
lines(1:14, plt6.metro.rmse, col="red")

points(1:14, plt2.metro.rmse, col="magenta", pch=16)
points(1:14, plt3.metro.rmse, col="green", pch=16)
points(1:14, plt4.metro.rmse, col="darkorange1", pch=16)
points(1:14, plt5.metro.rmse, col="blue", pch=16)
points(1:14, plt6.metro.rmse, col="red", pch=16)


plot(plt1.metro.crmse, col="black", pch=16,
     ylim = c(min(plt1.metro.crmse, plt2.metro.crmse, plt3.metro.crmse, plt4.metro.crmse, plt5.metro.crmse, plt6.metro.crmse), 
              max(plt1.metro.crmse, plt2.metro.crmse, plt3.metro.crmse, plt4.metro.crmse, plt5.metro.crmse, plt6.metro.crmse)), 
     main = "", ylab="RMSE", xlab = "day", cex.lab=1.3)
lines(1:14, plt1.metro.crmse, col="black")
lines(1:14, plt2.metro.crmse, col="magenta")
lines(1:14, plt3.metro.crmse, col="green")
lines(1:14, plt4.metro.crmse, col="darkorange1")
lines(1:14, plt5.metro.crmse, col="blue")
lines(1:14, plt6.metro.crmse, col="red")

points(1:14, plt2.metro.crmse, col="magenta", pch=16)
points(1:14, plt3.metro.crmse, col="green", pch=16)
points(1:14, plt4.metro.crmse, col="darkorange1", pch=16)
points(1:14, plt5.metro.crmse, col="blue", pch=16)
points(1:14, plt6.metro.crmse, col="red", pch=16)

par(fig=c(0,1,0,1), oma = c(0, 0, 0, 0), mar = c(0, 5, 0, 0), new = TRUE)
# plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
plot.new()
legend("bottom", legend=c("Nodewise ARIMA", "LOCAAT", "GNARI", "GNARI + LOCAAT", "GFT", "SGWT"), 
       col=c("black", "magenta", "green", "darkorange1", "blue", "red"), pch=16, horiz=TRUE, lty=1, 
       lwd=2, pt.cex=1.2, cex=1, bty="n", x.intersp=1, text.width=0.11)
