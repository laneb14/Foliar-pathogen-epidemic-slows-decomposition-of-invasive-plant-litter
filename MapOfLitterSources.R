
##### Load packages #####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(ggmap)
library(cowplot)

# Set working directory as needed
setwd("C:/Users/brett/Dropbox/Decomposition Experiment/FinalCodeAndData/UsingDataFromEDI")

##### Input files #####
DataFolder <- "edi.1434.1"
OriginalCoor <- read.csv(paste(sep='/',DataFolder, "Big Oaks Coordinates.csv"))
My.API.Key <- # Insert your API key here
register_google(My.API.Key)

# Offset the experimental site ever so slightly to allow to be visible on the map
CoorList <- OriginalCoor
CoorList$Lat <- ifelse(OriginalCoor$Site=='Common Garden',OriginalCoor$Lat+0.00175,OriginalCoor$Lat)

##### Make Map of Big Oaks National Wild Refuge Sites #####
myLocation <- c(-85.47, 38.9, -85.375, 39.06) #west most, south most, east most, north most
LongCenter <- -85.40
LatCenter <- 38.975

myMap1 <- get_googlemap(center=c(lon=LongCenter,lat=LatCenter),zoom=12,scale=2,maptype='terrain',style = c(feature = "all", element = "labels", visibility = "off"),)

ggmap(myMap1)+  geom_point(data=CoorList,inherit.aes = F,aes(x=Long,y=Lat,color=as.factor(Infection),shape=as.factor(Infection)),size=2)

WestPoint <- (-85.375) # Point used for calculating scale
EastPoint <- (-85.35191428) # Point used for calculating scale
Mid.EW.Point <- mean(EastPoint,WestPoint) # Get the mean east/west point of these points
VertPoint <- 38.9 # Set a location for the scale

BigOaksMap <-
ggmap(myMap1)+
  theme(panel.border = element_rect(colour = "black", fill=NA), panel.background = element_blank() )+
  ylim(38.90,39.05)+
  scale_x_continuous(breaks=c(-85.475,-85.425,-85.375),limits=c(-85.475,-85.35))+
  geom_point(data=CoorList,inherit.aes = F,aes(x=Long,y=Lat,color=as.factor(Infection),shape=as.factor(Infection)),size=2)+
  scale_color_manual(values=c('Non-infected'='blue','Infected'='red','Experimental Site'='black'),
                     labels = c("Non-infected"='Healthy site',"Infected"='Infected site',"Experimental Site"='Experimental site'), name='')+
  scale_shape_manual(values=c('Non-infected'=16,'Infected'=16,'Experimental Site'=17),
                     labels = c("Non-infected"='Healthy site',"Infected"='Infected site',"Experimental Site"='Experimental site'), name='')+  
  theme(legend.key = element_rect(colour = "transparent", fill = "transparent"))+  
  theme(legend.position = 'bottom' )+
  geom_segment(aes(x=WestPoint,xend=EastPoint,y=VertPoint,yend=VertPoint),size=1.3)+
  geom_text(aes(x=WestPoint+0.01,y=VertPoint+0.0032),label='2 km')+
  xlab("Longitude")+ylab('Latitude')


BigOaksMap
# ggsave(paste(sep='',format(Sys.time(),"%d%b%Y"),'-',"BigOaksMap.jpeg"),BigOaksMap,dp=1000,height = (3.5)*2 ,width = (2.25)*2)


States <- map_data("state")
RegionMap <-
ggplot(States,aes(long,lat))+theme_bw()+theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill='grey90'))+
  geom_polygon(aes(group=group),fill='white',color='black')+
  coord_map(xlim=c(-95,-75), ylim=c(34,44))+
  geom_point(x=LongCenter,y=LatCenter,color='red',size=2)+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())

# ggsave(paste(sep='',format(Sys.time(),"%d%b%Y"),'-',"RegionMap.jpeg"),RegionMap,dp=1000,height = (1.75)*2 ,width = (2.25)*2)


##### Combine all of the images #####
library(magick)
# Width X Height
RegionImage <- image_read("30Oct2024-RegionMap.jpeg") # 4500x3500
Crop.RegionImage <- image_crop(RegionImage,"4500x3000+0+250") # 4500x3000
BigOaksImage <- image_read("30Oct2024-BigOaksMap.jpeg") # 4500x7000
# Combined Height has been cropped to 10,000

PlotImage <- image_read('PlotPhoto.png') # 1038x779
Resized.PlotImage <- image_resize(PlotImage,"x1000") # 1332x1000
Resized.PlotImage.Padded <- image_blank(1332+50, 1000, color = 'white') %>% image_composite(Resized.PlotImage,gravity='east')
Resized.PlotImage.Annotated <- image_annotate(Resized.PlotImage.Padded,'C',gravity='northwest',size=48)

Resize.Crop.RegionImage <- image_resize(Crop.RegionImage,"x300") # 450x300
Resize.Crop.RegionImage.Padded <- image_blank(500, 300, color = 'white') %>% image_composite(Resize.Crop.RegionImage,gravity='east')
Resize.Crop.RegionImage.Annotated <- image_annotate(Resize.Crop.RegionImage.Padded,'A',gravity='northwest',size=48)

Resize.BigOaksImage <- image_resize(BigOaksImage,"x700") # 450x700
Resize.BigOaksImage.Padded <- image_blank(500, 700, color = 'white') %>% image_composite(Resize.BigOaksImage,gravity='east')
Resize.BigOaksImage.Annotated <- image_annotate(Resize.BigOaksImage.Padded,'B',gravity='northwest',size=48)

CombinedImage <- 
image_blank(1832+50, 1000, color = 'white') %>%
  image_composite(Resize.Crop.RegionImage.Annotated,gravity='northwest') %>%
  image_composite(Resize.BigOaksImage.Annotated,gravity='southwest') %>%
  image_composite(Resized.PlotImage.Annotated,gravity='east')
  
# image_write(CombinedImage, path = "Combined.png", format = "png")
