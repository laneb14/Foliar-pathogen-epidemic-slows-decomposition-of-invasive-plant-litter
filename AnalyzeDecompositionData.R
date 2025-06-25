
# Note, throughout this code "Non-infected" and "Healthy" are used as synonyms. These are often corrected to "Healthy" in post processing

##### Load packages #####
library(dplyr)
library(lmerTest)
library(expss)
library(cowplot)
library(ggplot2)
library(ggResidpanel)
library(lsr)

#### Set working directory as needed #####
# Using rstudioapi::getSourceEditorContext()$path will set the working directory to where this R file is saved
# Taking this approach allows the user to be able to access code saved on a cloud directory (i.e. dropbox) from multiple computers without having to change setwd() each time (e.g. C:\ vs D:\)
setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()

##### Adding a function to calculate SE #####
# Note, this was added later and may not be used for all SE calculations
FunctionMeanSE <- function(measurements){ measurements <- measurements[!is.na(measurements)] ; if(length(measurements)==0){stop('length was zero after removing NA values')};  m <- mean(measurements) ; se <- sd(measurements)/sqrt(length(measurements)) ; return(paste(m,"±",se))  }

##### Persistent variables #####
Months <- c("February","April","June","August","October","December")
AllMonths <- c("Start",Months)

LevelOrder <- c('Infected','Non-infected','Common Garden')
LevelNames <- c('Non-infected'="Healthy\nSites", 'Infected'="Infected\nSites",'Common Garden'="Experimental\nSite")
LevelNames.NoNewLine <- c('Non-infected'="Healthy Sites", 'Infected'="Infected Sites",'Common Garden'="Common Garden Sites")
LevelColors <- c('Non-infected'="#619CFF", 'Infected'="#F8766D",'Common Garden'="#00BA38")
LevelOrder.NoCG <- c('Infected','Non-infected')
LevelColors.NoCG <- LevelColors[names(LevelColors)!='Common Garden']
LevelNames.NoCG <- LevelNames[names(LevelNames)!='Common Garden']
RemChartMax <- 125
PointSize <- 1
TitleSize <- 10
CG.RemChartMax <- 125
Overall.CG.RemChartMax <- 125
CG.PointSize <- .75
CGAdj <- 5
LevelOrder.CG <- LevelOrder[LevelOrder!='Infected']
LevelColors.CG <- LevelColors[names(LevelColors)!='Infected']
LevelNames.CG <- c('Non-infected'="Non-infected\nSites", 'Common Garden'="Experimental Site")
AxisFontSize <- theme(axis.text.x = element_text(size=8,angle=0), legend.text=element_text(size=7))

##### Load and format data #####
DataFolder <- "edi.1434.1"
Chart <- read.table(file = paste(sep='/',DataFolder,'MicrostegiumDecompositionDataCollected2020.tsv'), sep = '\t', header = TRUE) # Load data
head(Chart) # Confirm the data read in correctly

# Descriptions of each column included as comment
Chart$Month <- as.factor(Chart$Month) # Month of sampling
Chart$Month.Number <- as.numeric(Chart$Month.Number) # Months since beginning of experiment, December 2019 is 0
Chart$Litter.Source <- as.factor(Chart$Litter.Source) # Source of litter
Chart$Infection <- as.factor(Chart$Infection) # Infection status of litter
Chart$Starting.Dry.Mass <- as.numeric(Chart$Starting.Dry.Mass) # Starting dry mass of each litter bag in grams
Chart$Final.Dry.Mass <- as.numeric(Chart$Final.Dry.Mass) # Final dry mass of each litter bag in grams
Chart$InitialStemProp <- as.numeric(Chart$InitialStemProp) # Starting proportion of each litter bag that is stem material, measured by Litter.Source
Chart$StemProp <- as.numeric(Chart$StemProp) # Proportion of sample (as dry mass) which was comprised of stem tissue
Chart$Replicate <- as.factor(Chart$Replicate) # Experimental replicate plot
Chart$ADLom <- as.numeric(Chart$ADLom) # Acid detergent lignin organic matter basis
Chart$pC <- as.numeric(Chart$pC) # percent carbon
Chart$pN <- as.numeric(Chart$pN) # percent nitrogen
Chart$CN <- as.numeric(Chart$CN) # carbon:nitrogen ratio

# head(Chart)

##### Calculate values #####
Chart$Mass.Lost <- Chart$Starting.Dry.Mass-Chart$Final.Dry.Mass # Calculated in grams, negative values would represent an increase in mass over time
Chart$Prop.Mass.Rem <- Chart$Final.Dry.Mass/Chart$Starting.Dry.Mass # Proportion of the original sample remaining
Chart$StemMass <- Chart$StemProp*Chart$Final.Dry.Mass
Chart$InitialStemMass <- Chart$InitialStemProp*Chart$Starting.Dry.Mass
Chart$Prop.StemMass.Rem <- Chart$StemMass/Chart$InitialStemMass
Chart$LeafMass <- (1-Chart$StemProp)*Chart$Final.Dry.Mass
Chart$InitialLeafMass <- (1-Chart$InitialStemProp)*Chart$Starting.Dry.Mass
Chart$Prop.LeafMass.Rem <- Chart$LeafMass/Chart$InitialLeafMass

##### Create Data Frames #####
Chart.NoStart <- Chart[Chart$Month.Number>0,]
Chart.NoStart.NoCG <- Chart.NoStart[Chart.NoStart$Litter.Source!='Common Garden',]
Chart.NoStart.NoCG.LeafStem <- Chart.NoStart.NoCG[Chart.NoStart.NoCG$Month.Number<=8,] # We didn't measure Leaf and Stem past the August sampling
Chart.AllTimes.NoCG <- Chart[Chart$Litter.Source!='Common Garden',]

##### Table 1: Anova of Infected vs Noninfected, without common garden #####

head(Chart.NoStart.NoCG)

### Overall ###
hist(Chart.NoStart.NoCG$Prop.Mass.Rem) ; shapiro.test(Chart.NoStart.NoCG$Prop.Mass.Rem) # 1.545e-6
# While this data is not normal, a test of the below transformations did not significant improve the fit of the data. We opted to retain the original data
# Transformations tried: sqrt, cubic root, squared, cubed, log, exp. All had no improvement

OverallAllMonthsModel <- lmer(Prop.Mass.Rem~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.NoStart.NoCG )
resid_panel(OverallAllMonthsModel)
anova(OverallAllMonthsModel)

### Leaf ###

hist(Chart.NoStart.NoCG$Prop.LeafMass.Rem) ; shapiro.test(Chart.NoStart.NoCG$Prop.LeafMass.Rem) # 0.003558
hist(sqrt(Chart.NoStart.NoCG$Prop.LeafMass.Rem)) ; shapiro.test(sqrt(Chart.NoStart.NoCG$Prop.LeafMass.Rem)) # normal
Chart.NoStart.NoCG$sqrt.Prop.LeafMass.Rem <- sqrt(Chart.NoStart.NoCG$Prop.LeafMass.Rem)

LeafAllMonthsModel <-  lmer(sqrt.Prop.LeafMass.Rem~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.NoStart.NoCG )
resid_panel(LeafAllMonthsModel)
anova(LeafAllMonthsModel)

# Effect size, leaf mass remaining in June
lsr::cohensD(sqrt.Prop.LeafMass.Rem~Infection,data=droplevels(Chart.NoStart.NoCG[Chart.NoStart.NoCG$Month.Number==6,]))

# post hoc of the month:infeciton interaction
Means.LeafAllMonths <- emmeans::emmeans(LeafAllMonthsModel, ~ Infection*Month)
pairs(Means.LeafAllMonths,simple='Infection',adjust='tukey')


### Stem ###

hist(Chart.NoStart.NoCG$Prop.StemMass.Rem) ; shapiro.test(Chart.NoStart.NoCG$Prop.StemMass.Rem) # 9.6859e-7
hist((Chart.NoStart.NoCG$Prop.StemMass.Rem)^2) ; shapiro.test((Chart.NoStart.NoCG$Prop.StemMass.Rem)^2) # 0.01316, as close to normal as we can achieve, residuals below are okay
Chart.NoStart.NoCG$sq.Prop.StemMass.Rem <- (Chart.NoStart.NoCG$Prop.StemMass.Rem)^2

StemAllMonthsModel <- lmer(sq.Prop.StemMass.Rem~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.NoStart.NoCG )
resid_panel(StemAllMonthsModel)
anova(StemAllMonthsModel)

### Lignin (Acid detergent lignin organic basis) ###
hist(Chart.AllTimes.NoCG$ADLom) ; shapiro.test(Chart.AllTimes.NoCG$ADLom) # 1.47e-6
hist((Chart.AllTimes.NoCG$ADLom)^(1/3)) ; shapiro.test((Chart.AllTimes.NoCG$ADLom)^(1/3)) # normal
Chart.AllTimes.NoCG$CubRoot.ADLom <- (Chart.AllTimes.NoCG$ADLom)^(1/3)

ADLom.Model <- lmer(CubRoot.ADLom~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.AllTimes.NoCG )
resid_panel(ADLom.Model)
anova(ADLom.Model)

# what percent higher is lignin in infected litter vs non-infected by month, includes min, max, and mean
for (x in c(0,2,4,6,8,10,12)){
  if(x==0){temp<-vector()}
  TempDF <- Chart.AllTimes.NoCG[Chart.AllTimes.NoCG$Month.Number==x,]
  TempDF <- TempDF[which(!is.na(TempDF$ADLom)),]
  m <- droplevels(TempDF$Month[1])
  y<-(mean(TempDF$ADLom[TempDF$Infection=='Infected'])/mean(TempDF$ADLom[TempDF$Infection=='Non-infected']))
  print(paste(m,round((y-1)*100,1))  ) ; temp <- c(temp,y)
  if(x==12){print(''); print(paste('Min',round((min(temp)-1)*100,1)) )}
  if(x==12){ print(paste('Max',round((max(temp)-1)*100,1)) )}
  if(x==12){ print(paste('Average',round((mean(temp)-1)*100,1)) ); rm(temp)  }
  rm(x,y,m,TempDF)
}

# Effect sizes of lignin in infected litter vs non-infected by month, includes min, max, and mean
for (x in c(0,2,4,6,8,10,12)){
  if(x==0){temp<-vector()}
  TempDF <- Chart.AllTimes.NoCG[Chart.AllTimes.NoCG$Month.Number==x,]
  TempDF <- TempDF[which(!is.na(TempDF$ADLom)),]
  m <- droplevels(TempDF$Month[1])
  y <- lsr::cohensD(TempDF$ADLom[TempDF$Infection=='Infected'],TempDF$ADLom[TempDF$Infection=='Non-infected'])
  print(paste(m,round(y,3))  ) ; temp <- c(temp,y)
  if(x==12){print(''); print(paste('Min',round(min(temp),3)) )}
  if(x==12){ print(paste('Max',round(max(temp),3)) )}
  if(x==12){ print(paste('Average',round(mean(temp),3)) ); rm(temp)  }
  rm(x,y,m,TempDF)
}

### CN Ratio ###
hist(Chart.AllTimes.NoCG$CN) ; shapiro.test(Chart.AllTimes.NoCG$CN) # 8.831e-15
hist(log10(Chart.AllTimes.NoCG$CN)) ; shapiro.test(log10(Chart.AllTimes.NoCG$CN)) # 3.465e-5, not perfect, but it is by far the best fit of all transformations tried
Chart.AllTimes.NoCG$Log.CN <- log10(Chart.AllTimes.NoCG$CN)

CN.Model <- lmer(Log.CN~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.AllTimes.NoCG )
resid_panel(CN.Model)
anova(CN.Model)

##### Table S3 #####
for (z in c('Overall','Leaf','Stem')){
  # if(z=='Overall'){MassRemDF <- data.frame()}
  for(x in Months){
    if(x=='October'&z!='Overall'){next}
    if(x=='December'&z!='Overall'){next}
    TempDF <- Chart[Chart$Month==x,]
    if(z=='Overall'){TempDF$REM <- TempDF$Prop.Mass.Rem}
    if(z=='Leaf'){TempDF$REM <- TempDF$Prop.LeafMass.Rem}
    if(z=='Stem'){TempDF$REM <- TempDF$Prop.StemMass.Rem}
    i<-TempDF$REM[TempDF$Infection=='Infected']
    n<-TempDF$REM[TempDF$Infection=='Non-infected']
    c<-TempDF$REM[TempDF$Infection=='Common Garden']
    print(paste(sep=',',z,x,paste(sep='',(mean(i,na.rm=T)*100)," (",(sd(i,na.rm=T)*100)/sqrt(length(i)),')'),paste(sep='',(mean(n,na.rm=T)*100)," (",(sd(n,na.rm=T)*100)/sqrt(length(n)),')'),paste(sep='',(mean(c,na.rm=T)*100)," (",(sd(c,na.rm=T)*100)/sqrt(length(c)),')')))
    # This can be copy and pasted from the output as a tab delimited, then cleaned up in excel
    rm(TempDF,i,n,c)
  }
  rm(x,z)
}

##### Table S4, first two columns #####

for (y in c("Infected",'Non-infected','Common Garden')){
  for(x in c('Start',Months)){
    TempDF <- Chart[Chart$Month==x & Chart$Infection==y,]
    print(paste(sep=',',y,x,paste(sep="",mean(TempDF$ADLom,na.rm=T)," (",(sd(TempDF$ADLom,na.rm=T)/sqrt(length(TempDF$ADLom))),")"),paste(sep="",mean(TempDF$CN,na.rm=T)," (",(sd(TempDF$CN,na.rm=T)/sqrt(length(TempDF$CN))),")")))
    # This can be copy and pasted from the output as a tab delimited, then cleaned up in excel
  }
  rm(x,TempDF,y)
}

##### Table S8, ANOVA results comparing common garden with non-infected litter #####

Chart.NoInf <- Chart[Chart$Infection!='Infected',]
Chart.NoInf.NoStart <- Chart.NoInf[Chart.NoInf$Month.Number!=0,]

### Overall ###
hist(Chart.NoInf.NoStart$Prop.Mass.Rem) ; shapiro.test(Chart.NoInf.NoStart$Prop.Mass.Rem) # 0.001543
# While this data is not normal, a test of the below transformations did not significant improve the fit of the data. We opted to retain the original data
# Transformations tried: sqrt, cubic root, squared, cubed, log, exp. All had no improvement

CG.Overall.Model <- lmer(Prop.Mass.Rem~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.NoInf.NoStart )
resid_panel(CG.Overall.Model)
anova(CG.Overall.Model)

# percent 
for (x in c(2,4,6,8,10,12)){
  if(x==2){temp <-vector()}
  TempDF <- Chart.NoInf.NoStart[Chart.NoInf.NoStart$Month.Number==x,] ; TempDF <- TempDF[which(!is.na(TempDF$Prop.Mass.Rem)),]
  (m <- TempDF$Month[1]) ; cg <- TempDF[TempDF$Infection=='Common Garden',] ; ni <- TempDF[TempDF$Infection=='Non-infected',]
  (y<-(((1-mean(cg$Prop.Mass.Rem))/(1-mean(ni$Prop.Mass.Rem))-1)*100))
  print(paste(m,round(y,2))) ; temp <- c(temp,y)
  if(x==12){
    print('') ; print(paste('Min',min(round(temp,2))))
    print(paste('Max',max(round(temp,2))));print(paste('Average',round(mean(temp),2)));rm(temp)
  }
  rm(x,y,TempDF,m,cg,ni)
}

# effect size 
for (x in c(2,4,6,8,10,12)){
  if(x==2){temp <-vector()}
  TempDF <- Chart.NoInf.NoStart[Chart.NoInf.NoStart$Month.Number==x,] ; TempDF <- TempDF[which(!is.na(TempDF$Prop.Mass.Rem)),]
  (m <- TempDF$Month[1]) ; cg <- TempDF[TempDF$Infection=='Common Garden',] ; ni <- TempDF[TempDF$Infection=='Non-infected',]
  (y <- lsr::cohensD(1-ni$Prop.Mass.Rem,1-cg$Prop.Mass.Rem))
  print(paste(m,round(y,2))) ; temp <- c(temp,y)
  if(x==12){
    print('') ; print(paste('Min',min(round(temp,2))))
    print(paste('Max',max(round(temp,2))));print(paste('Average',round(mean(temp),2)));rm(temp)
  }
  rm(x,y,TempDF,m,cg,ni)
}

### Leaf ###
hist(Chart.NoInf.NoStart$Prop.LeafMass.Rem) ; shapiro.test(Chart.NoInf.NoStart$Prop.LeafMass.Rem) # normal

CG.Leaf.Model <- lmer(Prop.LeafMass.Rem~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.NoInf.NoStart )
resid_panel(CG.Leaf.Model)
anova(CG.Leaf.Model)

### Stem ###
hist(Chart.NoInf.NoStart$Prop.StemMass.Rem) ; shapiro.test(Chart.NoInf.NoStart$Prop.StemMass.Rem) # 0.0001697
hist(Chart.NoInf.NoStart$Prop.StemMass.Rem^2) ; shapiro.test(Chart.NoInf.NoStart$Prop.StemMass.Rem^2) # normal
Chart.NoInf.NoStart$sq.Prop.StemMass.Rem <- Chart.NoInf.NoStart$Prop.StemMass.Rem^2

CG.Stem.Model <- lmer(sq.Prop.StemMass.Rem~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.NoInf.NoStart )
resid_panel(CG.Stem.Model)
anova(CG.Stem.Model)

### Lignin ###

hist(Chart.NoInf$ADLom) ; shapiro.test(Chart.NoInf$ADLom) # 1.697e-5
hist(sqrt(Chart.NoInf$ADLom)) ; shapiro.test(sqrt(Chart.NoInf$ADLom)) # normal
Chart.NoInf$sqrt.ADLom <- sqrt(Chart.NoInf$ADLom)

CG.Lignin.Model <- lmer(sqrt.ADLom~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.NoInf )
resid_panel(CG.Lignin.Model)
anova(CG.Lignin.Model)

### CN ###

hist(Chart.NoInf$CN) ; shapiro.test(Chart.NoInf$CN) # 9.601e-13
hist(log10(Chart.NoInf$CN)) ; shapiro.test(log10(Chart.NoInf$CN)) # 3.418e-6, not perfect, but as close as it can be transformed
Chart.NoInf$log.CN <- log10(Chart.NoInf$CN)

CG.CN.Model <- lmer(log.CN~Infection*Month + (1|Replicate)+(1|Litter.Source),Chart.NoInf )
resid_panel(CG.CN.Model)
anova(CG.CN.Model)

Means.CNAllMonths.CG.Healthy <- emmeans::emmeans(CG.CN.Model, ~ Infection*Month)
pairs(Means.CNAllMonths.CG.Healthy,simple='Infection',adjust='tukey')

lsr::cohensD(log.CN~Infection,data=droplevels(Chart.NoInf[Chart.NoInf$Month.Number==0,]))

for (x in c(0,2,4,6,8,10,12)){
  if(x==0){temp <-vector()}
  TempDF <- Chart.NoInf[Chart.NoInf$Month.Number==x,] ; TempDF <- TempDF[which(!is.na(TempDF$CN)),]
  (m <- TempDF$Month[1]) ; cg <- TempDF[TempDF$Infection=='Common Garden',] ; ni <- TempDF[TempDF$Infection=='Non-infected',]
  (y <- ((mean(cg$CN)/mean(ni$CN))-1)*100)
  print(paste(m,round(y,2))) ; temp <- c(temp,y)
  if(x==12){
    print('') ; print(paste('Min',min(round(temp,2))))
    print(paste('Max',max(round(temp,2))));print(paste('Average',round(mean(temp),2)));rm(temp)
  }
  rm(x,y,TempDF,m,cg,ni)
}

# effect size 
for (x in c(0,2,4,6,8,10,12)){
  if(x==0){temp <-vector()}
  TempDF <- Chart.NoInf[Chart.NoInf$Month.Number==x,] ; TempDF <- TempDF[which(!is.na(TempDF$log.CN)),]
  (m <- TempDF$Month[1]) ; cg <- TempDF[TempDF$Infection=='Common Garden',] ; ni <- TempDF[TempDF$Infection=='Non-infected',]
  (y <- lsr::cohensD(1-ni$log.CN,1-cg$log.CN))
  # stop()
  print(paste(m,round(y,2))) ; temp <- c(temp,y)
  if(x==12){
    print('') ; print(paste('Min',min(round(temp,2))))
    print(paste('Max',max(round(temp,2))));print(paste('Average',round(mean(temp),2)));rm(temp)
  }
  rm(x,y,TempDF,m,cg,ni)
}

##### Will be part of Figure 2 #####

head(Chart)

ADL.Chart <-
ggplot(Chart[Chart$Infection!='Common Garden',],aes(x=factor(Month,levels=AllMonths),y=ADLom,fill=factor(Infection,levels=c('Non-infected','Infected')) ) )+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0.2))+
  scale_fill_manual(values=LevelColors,labels=LevelNames.NoNewLine,name='Litter Source')+theme(legend.position = 'bottom')+
  xlab('')+ylab('Acid detergent lignin\n(as percent of organic mass)')

CN.Chart <- 
ggplot(Chart[Chart$Infection!='Common Garden',],aes(x=factor(Month,levels=AllMonths),y=CN,fill=factor(Infection,levels=c('Non-infected','Infected')) ) )+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(jitter.width = 0.15))+
  scale_fill_manual(values=LevelColors,labels=LevelNames.NoNewLine,name='Litter Source')+theme(legend.position = 'bottom')+
  xlab('')+ylab('Carbon/Nitrogen ratio')

Chemistry.Chart <- plot_grid(ADL.Chart+theme(legend.position='none'),CN.Chart,ncol=1)

# ggsave(paste(sep='',"ChemistryChart",format(Sys.time(),"%e%b%Y"),'.jpg'),Chemistry.Chart,dpi=1000,width=6,height=6)

##### Prep data frame #####

for (x in c(2,4,6,8,10,12)){
  if(x==2){
    DF.For.Figures <- data.frame(Month.Number=vector(),Material=vector(),Infection=vector(),Avg=vector(),Low=vector(),High=vector())
  }
  for (y in LevelOrder){
    TempDF <- Chart[Chart$Month.Number==x&Chart$Infection==y,]
    head(TempDF)
    margin <- (qt(0.975,df=(nrow(TempDF)-1))*sd(TempDF$Prop.Mass.Rem)/sqrt(nrow(TempDF)))
    l <- mean(TempDF$Prop.Mass.Rem,na.rm=T)-margin
    h <- mean(TempDF$Prop.Mass.Rem,na.rm=T)+margin
    DF.For.Figures <- rbind(DF.For.Figures, data.frame(Month.Number=x,Material='Overall',Infection=y,Avg=mean(TempDF$Prop.Mass.Rem,na.rm=T),Low=l,High=h) )
    rm(margin,l,h)
    if(x<9){
      margin <- (qt(0.975,df=(nrow(TempDF)-1))*sd(TempDF$Prop.LeafMass.Rem)/sqrt(nrow(TempDF)))
      l <- mean(TempDF$Prop.LeafMass.Rem)-margin
      h <- mean(TempDF$Prop.LeafMass.Rem)+margin
      DF.For.Figures <- rbind(DF.For.Figures, data.frame(Month.Number=x,Material="Leaf",Infection=y,Avg=mean(TempDF$Prop.LeafMass.Rem ) ,Low=l, High=h )  )
      rm(margin,l,h)
      
      margin <- (qt(0.975,df=(nrow(TempDF)-1))*sd(TempDF$Prop.StemMass.Rem)/sqrt(nrow(TempDF)))
      l <- mean(TempDF$Prop.StemMass.Rem)-margin
      h <- mean(TempDF$Prop.StemMass.Rem)+margin
      DF.For.Figures <- rbind(DF.For.Figures, data.frame(Month.Number=x,Material="Stem",Infection=y,Avg=mean(TempDF$Prop.StemMass.Rem ) ,Low=l, High=h )  )
      rm(margin,l,h)
    }
  }
  if(x==12){
    DF.For.Figures$Avg2 <- DF.For.Figures$Avg*100
    DF.For.Figures$High2 <- DF.For.Figures$High*100
    DF.For.Figures$Low2 <- DF.For.Figures$Low*100
    rm(x,y)
  }
}

##### Figure S4 #####
# I think this is now Fig S2 __
head(DF.For.Figures)
head(Chart)

CG.NI.Overall.Plot <- 
ggplot(DF.For.Figures[DF.For.Figures$Infection!='Infected'&DF.For.Figures$Material=='Overall',],
       aes(x=(Month.Number),y=Avg2,color=factor(Infection,levels=c('Non-infected','Common Garden'))) )+
  scale_color_manual(values=LevelColors ,labels=LevelNames,name='Litter Source')+
  scale_fill_manual(values=LevelColors  ,labels=LevelNames,name='Litter Source')+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2,4,6,8,10,12),labels=c('February','April','June','August','October','December'))+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  # Also need to add the jittering of points
  geom_point(inherit.aes = F, data=Chart[Chart$Infection!='Infected'&Chart$Month.Number!=0,],
             aes(x=Month.Number,y=Prop.Mass.Rem*100,fill=factor(Infection,levels=c('Non-infected','Common Garden'))),
             pch=21,size=1,position=position_jitterdodge(jitter.width = 0.2) )+
  scale_y_continuous(breaks = c(0,25,50,75,100,125),limits=c(-3,128))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab('Mass remaining (%)')+xlab('')+ggtitle("Overall")+theme(plot.title = element_text(size=TitleSize))


CG.NI.Leaf.Plot <-
ggplot(DF.For.Figures[DF.For.Figures$Infection!='Infected'&DF.For.Figures$Material=='Leaf',],
       aes(x=(Month.Number),y=Avg2,color=factor(Infection,levels=c('Non-infected','Common Garden'))) )+
  scale_color_manual(values=LevelColors ,labels=LevelNames,name='Litter Source')+
  scale_fill_manual(values=LevelColors  ,labels=LevelNames,name='Litter Source')+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2,4,6,8,10,12),labels=c('February','April','June','August','October','December'))+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  # Also need to add the jittering of points
  geom_point(inherit.aes = F, data=Chart[Chart$Infection!='Infected'&Chart$Month.Number!=0&Chart$Month.Number<9,],
             aes(x=Month.Number,y=Prop.LeafMass.Rem*100,fill=factor(Infection,levels=c('Non-infected','Common Garden'))),
             pch=21,size=1,position=position_jitterdodge(jitter.width = 0.2) )+
  scale_y_continuous(breaks = c(0,25,50,75,100,125),limits=c(-3,128))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab('Mass remaining (%)')+xlab('')+ggtitle('Leaf')+theme(plot.title = element_text(size=TitleSize))

CG.NI.Stem.Plot <- 
ggplot(DF.For.Figures[DF.For.Figures$Infection!='Infected'&DF.For.Figures$Material=='Stem',],
       aes(x=(Month.Number),y=Avg2,color=factor(Infection,levels=c('Non-infected','Common Garden'))) )+
  scale_color_manual(values=LevelColors ,labels=LevelNames,name='Litter Source')+
  scale_fill_manual(values=LevelColors  ,labels=LevelNames,name='Litter Source')+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2,4,6,8,10,12),labels=c('February','April','June','August','October','December'))+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  # Also need to add the jittering of points
  geom_point(inherit.aes = F, data=Chart[Chart$Infection!='Infected'&Chart$Month.Number!=0&Chart$Month.Number<9,],
             aes(x=Month.Number,y=Prop.StemMass.Rem*100,fill=factor(Infection,levels=c('Non-infected','Common Garden'))),
             pch=21,size=1,position=position_jitterdodge(jitter.width = 0.2) )+
  scale_y_continuous(breaks = c(0,25,50,75,100,125),limits=c(-3,128))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab('Mass remaining (%)')+xlab('')+ggtitle('Stem')+theme(plot.title = element_text(size=TitleSize))

CG.NI.Chart <-
  plot_grid(CG.NI.Overall.Plot,
          plot_grid(CG.NI.Leaf.Plot+theme(legend.position='none'),
                    CG.NI.Stem.Plot+theme(legend.position='none')+ylab(''),
                    nrow=1,labels = c('B','C')),
          ncol=1,labels = ('A') )

# ggsave(paste(sep='',"Common Garden vs Healthy",format(Sys.time(),"%e%b%Y"),'.jpg'),CG.NI.Chart,dpi=1000,width=6,height=5)

##### Figure 2 #####

Fig2.Overall.Plot <-
  ggplot(DF.For.Figures[DF.For.Figures$Infection!='Common Garden'&DF.For.Figures$Material=='Overall',],
         aes(x=(Month.Number),y=Avg2,color=factor(Infection,levels=c('Infected','Non-infected'))) )+
  scale_color_manual(values=LevelColors ,labels=LevelNames,name='Litter Source')+
  scale_fill_manual(values=LevelColors  ,labels=LevelNames,name='Litter Source')+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2,4,6,8,10,12),labels=c('February','April','June','August','October','December'))+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  geom_point(inherit.aes = F, data=Chart[Chart$Infection!='Common Garden'&Chart$Month.Number!=0,],
             aes(x=Month.Number,y=Prop.Mass.Rem*100,fill=factor(Infection,levels=c('Infected','Non-infected'))),
             pch=21,size=1,position=position_jitterdodge(jitter.width = 0.2) )+
  scale_y_continuous(breaks = c(0,25,50,75,100,125),limits=c(-3,128))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab('Mass remaining (%)')+xlab('')+ggtitle("Overall")+theme(plot.title = element_text(size=TitleSize))


Fig2.Leaf.Plot <-
  ggplot(DF.For.Figures[DF.For.Figures$Infection!='Common Garden'&DF.For.Figures$Material=='Leaf',],
         aes(x=(Month.Number),y=Avg2,color=factor(Infection,levels=c('Infected','Non-infected'))) )+
  scale_color_manual(values=LevelColors ,labels=LevelNames,name='Litter Source')+
  scale_fill_manual(values=LevelColors  ,labels=LevelNames,name='Litter Source')+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2,4,6,8,10,12),labels=c('February','April','June','August','October','December'))+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  # Also need to add the jittering of points
  geom_point(inherit.aes = F, data=Chart[Chart$Infection!='Common Garden'&Chart$Month.Number!=0&Chart$Month.Number<9,],
             aes(x=Month.Number,y=Prop.LeafMass.Rem*100,fill=factor(Infection,levels=c('Infected','Non-infected'))),
             pch=21,size=1,position=position_jitterdodge(jitter.width = 0.2) )+
  scale_y_continuous(breaks = c(0,25,50,75,100,125),limits=c(-3,128))+
  geom_segment(aes(x=5.75,xend=6.25,y=122,yend=122),color='black',size=1)+
  geom_text(aes(x=6,y=124,label='*'),size=8,color='black')+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab('Mass remaining (%)')+xlab('')+ggtitle('Leaf')+theme(plot.title = element_text(size=TitleSize))

Fig2.Stem.Plot <-
  ggplot(DF.For.Figures[DF.For.Figures$Infection!='Common Garden'&DF.For.Figures$Material=='Stem',],
         aes(x=(Month.Number),y=Avg2,color=factor(Infection,levels=c('Infected','Non-infected'))) )+
  scale_color_manual(values=LevelColors ,labels=LevelNames,name='Litter Source')+
  scale_fill_manual(values=LevelColors  ,labels=LevelNames,name='Litter Source')+
  geom_line(size=1)+
  scale_x_continuous(breaks=c(2,4,6,8,10,12),labels=c('February','April','June','August','October','December'))+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  # Also need to add the jittering of points
  geom_point(inherit.aes = F, data=Chart[Chart$Infection!='Common Garden'&Chart$Month.Number!=0&Chart$Month.Number<9,],
             aes(x=Month.Number,y=Prop.StemMass.Rem*100,fill=factor(Infection,levels=c('Infected','Non-infected'))),
             pch=21,size=1,position=position_jitterdodge(jitter.width = 0.2) )+
  scale_y_continuous(breaks = c(0,25,50,75,100,125),limits=c(-3,128))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab('Mass remaining (%)')+xlab('')+ggtitle('Stem')+theme(plot.title = element_text(size=TitleSize))

Figure2 <-
  plot_grid(Fig2.Overall.Plot,
            plot_grid(Fig2.Leaf.Plot+theme(legend.position='none'),
                      Fig2.Stem.Plot+theme(legend.position='none')+ylab(''),
                      nrow=1,labels = c('B','C')),
            ncol=1,labels = ('A') )

# Figure2
# Not using this one in publication, but save because I may use it in seminars
### ggsave(paste(sep='',"Infected vs Healthy",format(Sys.time(),"%e%b%Y"),'.jpg'),Figure2,dpi=1000,width=6,height=5)

##### One large combined litter chart #####

Figure2
plot_grid(Fig2.Overall.Plot,
          plot_grid(Fig2.Leaf.Plot+theme(legend.position='none'),
                    Fig2.Stem.Plot+theme(legend.position='none')+ylab(''),
                    nrow=1,labels = c('B','C')),
          ncol=1,labels = ('A') )

CombinedChartLitter <-
plot_grid(Figure2,
          ADL.Chart+theme(legend.position='none'),
          CN.Chart+theme(legend.position='none'),
          ncol=1,rel_heights = c(2,1,1),labels=c('','D','E'),vjust=(0))

CombinedChartLitter
# ggsave(paste(sep='',"Combined Litter Chart",format(Sys.time(),"%e%b%Y"),'.jpg'),CombinedChartLitter,dpi=1000,width=6,height=10)

