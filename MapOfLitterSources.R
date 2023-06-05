
##### Load packages #####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(ggmap)

# Set working directory as needed
setwd("C:/Users/brett/Dropbox/Decomposition Experiment/FinalCodeAndData/UsingDataFromEDI")

##### Input files #####
DataFolder <- "edi.1434.1"
OriginalCoor <- read.csv(paste(sep='/',DataFolder, "Big Oaks Coordinates.csv"))

# Offset the experimental site ever so slightly to allow to be visible on the map
CoorList <- OriginalCoor
CoorList$Lat <- ifelse(OriginalCoor$Site=='Common Garden',OriginalCoor$Lat+0.00175,OriginalCoor$Lat)

##### Make Map of Big Oaks National Wild Refuge Sites #####
myLocation <- c(-85.47, 38.9, -85.375, 39.06) #west most, south most, east most, north most
myMap <- get_map(location=myLocation, source="stamen", maptype= 'terrain', crop=FALSE) # get the map using locations input
WestPoint <- (-85.375) # Point used for calculating scale
EastPoint <- (-85.35191428) # Point used for calculating scale
Mid.EW.Point <- mean(EastPoint,WestPoint) # Get the mean east/west point of these points
VertPoint <- 38.9 # Set a location for the scale

BigOaksMap <- 
ggmap(myMap)+
  xlab("Longitude")+ylab('Latitude')+ theme(panel.border = element_rect(colour = "black", fill=NA), panel.background = element_blank() )+
  scale_x_continuous(breaks=c(-85.475,-85.425,-85.375))+ # This will produce a warning 'scale for x is already present', but sets the boundaries where I want
  geom_point(data=CoorList,inherit.aes = F,aes(x=Long,y=Lat,color=as.factor(Infection),shape=as.factor(Infection)),size=2)+
  scale_color_manual(values=c('Non-infected'='blue','Infected'='red','Experimental Site'='black'),
                     labels = c("Non-infected"='Non-infected site',"Infected"='Infected site',"Experimental Site"='Experimental site'), name='')+
  scale_shape_manual(values=c('Non-infected'=16,'Infected'=16,'Experimental Site'=17),
                     labels = c("Non-infected"='Non-infected site',"Infected"='Infected site',"Experimental Site"='Experimental site'), name='')+
  geom_point(data=CoorList[CoorList$Infection=='Experimental Site',],inherit.aes = F,
             aes(x=Long,y=Lat,color=as.factor(Infection),shape=as.factor(Infection)),size=2,shape=17, show.legend = FALSE)+
  theme(legend.key = element_rect(colour = "transparent", fill = "transparent"))+
  geom_segment(aes(x=WestPoint,xend=EastPoint,y=VertPoint,yend=VertPoint),size=1.3)+ geom_text(aes(x=WestPoint+0.01,y=VertPoint-0.0032),label='2 km')+
  theme(legend.position = 'bottom',legend.box = 'vertical' )

BigOaksMap
# ggsave(paste(sep='',format(Sys.time(),"%d%b%Y"),'-',"BigOaksMap.jpeg"),Plot,dp=1000)

