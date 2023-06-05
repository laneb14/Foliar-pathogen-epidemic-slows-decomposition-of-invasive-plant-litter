
# Packages
library(ggplot2)
library(grid)
library(dplyr)
library(tidyr)
library(expss)

# Set working directory as needed
setwd("C:/Users/brett/Dropbox/Decomposition Experiment/FinalCodeAndData/UsingDataFromEDI")

MonthNumbers = c('0'='December 2019','1'='January 2020','2'="February 2020",'3'="March 2020",'4'="April 2020",'5'='May 2020','6'="June 2020",'7'="July 2020",
                 '8'='August 2020','9'='September 2020','10'='October 2020','11'='November 2020','12'='December 2020')

DataFolder <- "edi.1434.1"
InputLogger <- read.table(file = paste(sep='/',DataFolder, 'LoggerData.tsv'), sep = '\t', header = TRUE) # Load metadata
colnames(InputLogger) <- c("LoggerNumber","Date.Time","Temp","RH")

head(InputLogger)
# InputLogger$Date.Time
InputLogger$Date <- as.Date(sub(" .*",'',InputLogger$Date.Time))
InputLogger$Hour24 <- sub(":.*",'', sub(".* ",'',InputLogger$Date.Time) )
tail(InputLogger)
LoggerDF <- InputLogger[InputLogger$Date!='2020-12-16',] # This day was not part of the experiment, the loggers were in the lab already
LoggerDF <- LoggerDF[(LoggerDF$Date!='2020-12-15')|(LoggerDF$Hour24<=12),] # Keep up until noon on this day
LoggerDF <- LoggerDF[(LoggerDF$Date!='2019-12-10')|(LoggerDF$Hour24>=12),] # Keep after noon on this day
head(LoggerDF)

InputRain <- read.csv(paste(sep='/',DataFolder, "RainDataForSubmission.csv"))
InputRain$Year <- sub("-.*",'',InputRain$Date)
InputRain$Month <- sub("....-","", sub("-..$",'',InputRain$Date))
InputRain$MonthNumber <- as.numeric(ifelse(InputRain$Year==2019,0,InputRain$Month))
InputRain$DayOfMonth <- as.numeric(sub(".*-",'',InputRain$Date))
RainDF <- InputRain[InputRain$MonthNumber>0|InputRain$DayOfMonth>=10,]
RainDF <- RainDF[RainDF$MonthNumber<12|RainDF$DayOfMonth<=15,]


## This loop will provide a unique observation number to each hour as well as a number to each day
for (x in 1:nrow(LoggerDF)){
  if(x==1){
    LoggerDF$Obs.Number <- NA # Give them a order for easier manipulation
    LoggerDF$DayNumber <- NA
    DayCount <- 1
  }
  if(LoggerDF$LoggerNumber[x]==1){
    if(LoggerDF$Hour24[x]=="00"){
      DayCount <- DayCount+1
    }
    LoggerDF$Obs.Number[x] <- x
    LoggerDF$DayNumber[x] <- DayCount
  }else{
    LoggerDF$Obs.Number[x] <-LoggerDF$Obs.Number[LoggerDF$Date.Time==LoggerDF$Date.Time[x]&LoggerDF$LoggerNumber==1]
    LoggerDF$DayNumber[x] <-LoggerDF$DayNumber[LoggerDF$Date.Time==LoggerDF$Date.Time[x]&LoggerDF$LoggerNumber==1]
  }
  if(x==nrow(LoggerDF)){
    rm(x,DayCount)
  }
}
head(LoggerDF)

## The RainDF has a few missing days, so loop through it to copy day numbers over from the LoggerDF
x <- 1
for (x in 1:nrow(RainDF) ){
  if(x==1){RainDF$DayNumber<-NA}
  RainDF$DayNumber[x] <- LoggerDF$DayNumber[LoggerDF$Date==as.Date(RainDF$Date[x])][1]
  rm(x)
}
head(RainDF)

## This loop will calculate the daily high, daily low, and average RH of each day.
## Daily high and low will be found for each data logger and the median will be used.
## Average RH will be computed for each data logger and the median selected
for (x in 1:max(LoggerDF$DayNumber)){
  if (x == 1){
    DailyWeather <- data.frame(Date=vector(),DayNumber=numeric(),DailyHigh=numeric(),DailyLow=numeric(),DailyRH=numeric())
  }
  TempDF1 <- LoggerDF[LoggerDF$DayNumber==x&LoggerDF$LoggerNumber==1,]
  TempDF2 <- LoggerDF[LoggerDF$DayNumber==x&LoggerDF$LoggerNumber==2,]
  TempDF3 <- LoggerDF[LoggerDF$DayNumber==x&LoggerDF$LoggerNumber==3,]
  TempDF4 <- LoggerDF[LoggerDF$DayNumber==x&LoggerDF$LoggerNumber==4,]
  TempDF5 <- LoggerDF[LoggerDF$DayNumber==x&LoggerDF$LoggerNumber==5,]
  
  High <-mean(max(TempDF1$Temp),max(TempDF2$Temp),max(TempDF3$Temp),max(TempDF4$Temp),max(TempDF5$Temp))
  Low <- median(min(TempDF1$Temp),min(TempDF2$Temp),min(TempDF3$Temp),min(TempDF4$Temp),min(TempDF5$Temp))
  RH <- median(mean(TempDF1$RH),mean(TempDF2$RH),mean(TempDF3$RH),mean(TempDF4$RH),mean(TempDF5$RH))
  
  DailyWeather <- rbind(DailyWeather, data.frame(Date=TempDF1$Date[1],DayNumber=x,DailyHigh=High,DailyLow=Low,DailyRH=RH)  )
  rm(TempDF1,TempDF2,TempDF3,TempDF4,TempDF5,High,Low,RH,x)
}

head(DailyWeather)
DailyWeather$Year <- sub("-.*",'',DailyWeather$Date)
DailyWeather$Month <- sub("....-",'', sub("-..$",'',DailyWeather$Date))
DailyWeather$MonthNumber <- as.numeric(ifelse(DailyWeather$Year==2019,0,DailyWeather$Month))

##### Table S2 #####
x <- 3
for (x in 0:12){
  if(x==0){
    TableS2 <- data.frame(Month=numeric(),High=numeric(),Low=numeric(),RH=numeric(),Percip=vector())
  }
  TempDF <- DailyWeather[DailyWeather$MonthNumber==x,]
  H <- round(((mean(TempDF$DailyHigh)-32)*5/9),1)
  L <- round(((mean(TempDF$DailyLow)-32)*5/9),1)
  R <- round(mean(TempDF$DailyRH),0)
  P <- round((sum(RainDF$PRCP[RainDF$MonthNumber==x])*2.54),1)
  TableS2 <- rbind(TableS2, data.frame(Month=MonthNumbers[names(MonthNumbers)==x],High=H,Low=L,RH=R,Percip=P) )
  rm(H,L,R,P,x,TempDF)
  
}

TableS2


##### Copy Rain amounts into the DailyWeather data frame, translate missing days into zeros #####

for (x in 1:nrow(DailyWeather)){
  if(x==1){DailyWeather$PRCP <- NA}
  TempDF <- RainDF[RainDF$DayNumber==DailyWeather$DayNumber[x],]
  if(nrow(TempDF)>0){
    if(nrow(TempDF)>1){ print("ERROR")  } # This shouldn't happen, but double check
    DailyWeather$PRCP[x] <- TempDF$PRCP[1]
  }else{
    DailyWeather$PRCP[x] <- 0
  }
}

##### compute 7 day sliding window for daily highs, lows, RH #####

for (x in 1:nrow(DailyWeather)){
  if(x==1){
    DailyWeather$SlideHigh <- NA
    DailyWeather$SlideLow <- NA
    DailyWeather$SlideRH <- NA
  }
  a <- max((x-13),0) # Adjusts sliding window for the first 14 days of data
  TempDF <- DailyWeather[DailyWeather$DayNumber>=a & DailyWeather$DayNumber<=x,]
  DailyWeather$SlideHigh[x] <-  mean(TempDF$DailyHigh)
  DailyWeather$SlideLow[x] <-  mean(TempDF$DailyLow)
  DailyWeather$SlideRH[x] <- mean(TempDF$DailyRH)
  if(x==nrow(DailyWeather)){
    DailyWeather$SlideHigh.C <- (DailyWeather$SlideHigh-32)*5/9
    DailyWeather$SlideLow.C <- (DailyWeather$SlideLow-32)*5/9
    rm(x,a,TempDF)
  }
}

##### Figure S1 ######

TempGraph <-
ggplot(DailyWeather, aes(x=DayNumber,y=SlideHigh.C))+geom_line(color = "red", cex = 0.75) + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  # ggtitle(paste(WindowLength,"Day Window")) + 
  geom_line(aes(x=DayNumber,y=SlideLow.C),cex=0.75, color = "blue") + ylab("Temperature (°C)")+ 
  scale_y_continuous(sec.axis = sec_axis(~.*9/5+32, name = "Temperature (°F)")   )+
  xlab("") +
  # theme(axis.text.x = element_blank() )+
  theme(axis.title.y = element_text(vjust = 3) ) +
  scale_x_continuous(breaks=c(23,114,205,297,389),labels=c("Jan 2020","Apr 2020","Jul 2020","Oct 2020","Jan 2021"))+
  geom_vline(xintercept =  1   , linetype=2 ) +
  geom_vline(xintercept =  77  , linetype=2 ) +
  geom_vline(xintercept =  137 , linetype=2 ) +
  geom_vline(xintercept =  199 , linetype=2 ) +
  geom_vline(xintercept =  260 , linetype=2 ) +
  geom_vline(xintercept =  324 , linetype=2 ) +
  geom_vline(xintercept =  372 , linetype=2 )

TempGraph


RHadjMultiple <- 25

RainRHGraph <- 
ggplot(DailyWeather,aes(x=DayNumber,y=SlideRH/RHadjMultiple))+geom_line(color="black", cex = 0.75)+ theme_bw() + ylab("Precipitation (inches)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  geom_col(aes(y=as.numeric(PRCP)))+
  scale_y_continuous(breaks = c(0,1,2,3,4), labels = c(0,25,50,75,100) )+
  scale_y_continuous(sec.axis = sec_axis(~.*RHadjMultiple, name = "Relative Humidity (%)")   )+
  theme(axis.title.y = element_text(vjust = 3) )+
  theme(axis.title.y.right = element_text(vjust=1))+
  xlab("")+
  scale_x_continuous(breaks=c(23,114,205,297,389),labels=c("Jan 2020","Apr 2020","Jul 2020","Oct 2020","Jan 2021"))+
  geom_vline(xintercept =  1   , linetype=2 ) +
  geom_vline(xintercept =  77  , linetype=2 ) +
  geom_vline(xintercept =  137 , linetype=2 ) +
  geom_vline(xintercept =  199 , linetype=2 ) +
  geom_vline(xintercept =  260 , linetype=2 ) +
  geom_vline(xintercept =  324 , linetype=2 ) +
  geom_vline(xintercept =  372 , linetype=2 )
  
grid.draw(rbind(ggplotGrob(TempGraph), ggplotGrob(RainRHGraph), size = "last"))


