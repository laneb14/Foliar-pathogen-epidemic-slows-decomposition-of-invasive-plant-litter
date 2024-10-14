
##### Load packages #####
library(dplyr)
library(lmerTest)
library(expss)
library(cowplot)
library(ggplot2)
library(ggResidpanel)

# Set working directory as needed
setwd("C:/Users/brett/Dropbox/Decomposition Experiment/FinalCodeAndData/7Oct2024-UsingEDIData")

##### Persistent variables #####
Months <- c("February","April","June","August","October","December")

##### Load and format data #####
DataFolder <- "edi.1434.1"
Chart <- read.table(file = paste(sep='/',DataFolder,'MicrostegiumDecompositionDataCollected2020.tsv'), sep = '\t', header = TRUE) # Load data

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
Chart.NoStart.NoCG

### Overall ###
hist(Chart.NoStart.NoCG$Prop.Mass.Rem)
shapiro.test(Chart.NoStart.NoCG$Prop.Mass.Rem) # 1.545e-6
# While this data is not normal, a test of the below transformations did not significant improve the fit of the data. We opted to retain the original data
# Transformations tried: sqrt, cubic root, squared, cubed, log, exp. All had no improvement

OverallAllMonthsModel <- lmer(Prop.Mass.Rem~Infection*Month + (1|Replicate),Chart.NoStart.NoCG )
resid_panel(OverallAllMonthsModel)
anova(OverallAllMonthsModel)

### Leaf ###

hist(Chart.NoStart.NoCG$Prop.LeafMass.Rem)
shapiro.test(Chart.NoStart.NoCG$Prop.LeafMass.Rem) # 0.003558
hist(sqrt(Chart.NoStart.NoCG$Prop.LeafMass.Rem))
shapiro.test(sqrt(Chart.NoStart.NoCG$Prop.LeafMass.Rem)) # normal
Chart.NoStart.NoCG$sqrt.Prop.LeafMass.Rem <- sqrt(Chart.NoStart.NoCG$Prop.LeafMass.Rem)

LeafAllMonthsModel <- lmer(sqrt.Prop.LeafMass.Rem~Infection*Month + (1|Replicate),Chart.NoStart.NoCG )
resid_panel(LeafAllMonthsModel)
anova(LeafAllMonthsModel)

Means.LeafAllMonths <- emmeans::emmeans(LeafAllMonthsModel, ~ Infection*Month)
pairs(Means.LeafAllMonths,simple='Infection')


### Stem ###


hist(Chart.NoStart.NoCG$Prop.StemMass.Rem)
shapiro.test(Chart.NoStart.NoCG$Prop.StemMass.Rem) # 9.6859e-7
hist((Chart.NoStart.NoCG$Prop.StemMass.Rem)^2)
shapiro.test((Chart.NoStart.NoCG$Prop.StemMass.Rem)^2) # 0.01316, as close to normal as we can achieve, residuals below are okay
Chart.NoStart.NoCG$sq.Prop.StemMass.Rem <- (Chart.NoStart.NoCG$Prop.StemMass.Rem)^2

StemAllMonthsModel <- lmer(sq.Prop.StemMass.Rem~Infection*Month + (1|Replicate),Chart.NoStart.NoCG )
resid_panel(StemAllMonthsModel)
anova(StemAllMonthsModel)

### Lignin (Acid detergent lignin organic basis) ###
hist(Chart.AllTimes.NoCG$ADLom)
shapiro.test(Chart.AllTimes.NoCG$ADLom) # 1.47e-6
hist((Chart.AllTimes.NoCG$ADLom)^(1/3))
shapiro.test((Chart.AllTimes.NoCG$ADLom)^(1/3)) # normal
Chart.AllTimes.NoCG$CubRoot.ADLom <- (Chart.AllTimes.NoCG$ADLom)^(1/3)

ADLom.Model <- lmer(CubRoot.ADLom~Infection*Month + (1|Replicate),Chart.AllTimes.NoCG )
resid_panel(ADLom.Model)
anova(ADLom.Model)

### CN Ratio ###
hist(Chart.AllTimes.NoCG$CN)
shapiro.test(Chart.AllTimes.NoCG$CN) # 8.831e-15
hist(log10(Chart.AllTimes.NoCG$CN))
shapiro.test(log10(Chart.AllTimes.NoCG$CN)) # 3.465e-5, not perfect, but it is by far the best fit of all transformations tried
Chart.AllTimes.NoCG$Log.CN <- log10(Chart.AllTimes.NoCG$CN)

CN.Model <- lmer(Log.CN~Infection*Month + (1|Replicate),Chart.AllTimes.NoCG )
# resid_panel(CN.Model)
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
hist(Chart.NoInf.NoStart$Prop.Mass.Rem)
shapiro.test(Chart.NoInf.NoStart$Prop.Mass.Rem) # 0.001543
# While this data is not normal, a test of the below transformations did not significant improve the fit of the data. We opted to retain the original data
# Transformations tried: sqrt, cubic root, squared, cubed, log, exp. All had no improvement

CG.Overall.Model <- lmer(Prop.Mass.Rem~Infection*Month + (1|Replicate),Chart.NoInf.NoStart )
resid_panel(CG.Overall.Model)
anova(CG.Overall.Model)

### Leaf ###
hist(Chart.NoInf.NoStart$Prop.LeafMass.Rem)
shapiro.test(Chart.NoInf.NoStart$Prop.LeafMass.Rem) # normal

CG.Leaf.Model <- lmer(Prop.LeafMass.Rem~Infection*Month + (1|Replicate),Chart.NoInf.NoStart )
resid_panel(CG.Leaf.Model)
anova(CG.Leaf.Model)

### Stem ###
hist(Chart.NoInf.NoStart$Prop.StemMass.Rem)
shapiro.test(Chart.NoInf.NoStart$Prop.StemMass.Rem) # 0.0001697
hist(Chart.NoInf.NoStart$Prop.StemMass.Rem^2)
shapiro.test(Chart.NoInf.NoStart$Prop.StemMass.Rem^2) # normal
Chart.NoInf.NoStart$sq.Prop.StemMass.Rem <- Chart.NoInf.NoStart$Prop.StemMass.Rem^2

CG.Stem.Model <- lmer(sq.Prop.StemMass.Rem~Infection*Month + (1|Replicate),Chart.NoInf.NoStart )
resid_panel(CG.Stem.Model)
anova(CG.Stem.Model)

### Lignin ###

hist(Chart.NoInf$ADLom)
shapiro.test(Chart.NoInf$ADLom) # 1.697e-5
hist(sqrt(Chart.NoInf$ADLom))
shapiro.test(sqrt(Chart.NoInf$ADLom)) # normal
Chart.NoInf$sqrt.ADLom <- sqrt(Chart.NoInf$ADLom)

CG.Lignin.Model <- lmer(sqrt.ADLom~Infection*Month + (1|Replicate),Chart.NoInf )
resid_panel(CG.Lignin.Model)
anova(CG.Lignin.Model)

### CN ###

hist(Chart.NoInf$CN)
shapiro.test(Chart.NoInf$CN) # 9.601e-13
hist(log10(Chart.NoInf$CN))
shapiro.test(log10(Chart.NoInf$CN)) # 3.418e-6, not perfect, but as close as it can be transformed
Chart.NoInf$log.CN <- log10(Chart.NoInf$CN)

CG.CN.Model <- lmer(log.CN~Infection*Month + (1|Replicate),Chart.NoInf )
resid_panel(CG.CN.Model)
anova(CG.CN.Model)




##### Figure 2 #####

# This code creates a dataframe which is in a more manageable form to create Figure 2
for (x in c(2,4,6,8,10,12) ){
  if(x==2){
    Fig2DF <- data.frame(Month.Number=vector(),Material=vector(),Infection=vector(),Avg=vector(),Low=vector(),High=vector())
  }
  for (y in LevelOrder.NoCG){
    TempDF <- Chart[Chart$Month.Number==x&Chart$Infection==y,]
    head(TempDF)
    margin <- (qt(0.975,df=(nrow(TempDF)-1))*sd(TempDF$Prop.Mass.Rem)/sqrt(nrow(TempDF)))
    l <- mean(TempDF$Prop.Mass.Rem)-margin
    h <- mean(TempDF$Prop.Mass.Rem)+margin
    Fig2DF <- rbind(Fig2DF, data.frame(Month.Number=x,Material="Overall",Infection=y,Avg=mean(TempDF$Prop.Mass.Rem),Low=l,High=h) )
    # next
    if(x<9){
      margin <- (qt(0.975,df=(nrow(TempDF)-1))*sd(TempDF$Prop.LeafMass.Rem)/sqrt(nrow(TempDF)))
      l <- mean(TempDF$Prop.LeafMass.Rem)-margin
      h <- mean(TempDF$Prop.LeafMass.Rem)+margin
      Fig2DF <- rbind(Fig2DF, data.frame(Month.Number=x,Material="Leaf",Infection=y,Avg=mean(TempDF$Prop.LeafMass.Rem ) ,Low=l, High=h )  )
      
      margin <- (qt(0.975,df=(nrow(TempDF)-1))*sd(TempDF$Prop.StemMass.Rem)/sqrt(nrow(TempDF)))
      l <- mean(TempDF$Prop.StemMass.Rem)-margin
      h <- mean(TempDF$Prop.StemMass.Rem)+margin
      Fig2DF <- rbind(Fig2DF, data.frame(Month.Number=x,Material="Stem",Infection=y,Avg=mean(TempDF$Prop.StemMass.Rem ) ,Low=l, High=h )  )
    }
  }
  if(x==12){
    Fig2DF$Avg2 <- Fig2DF$Avg*100
    Fig2DF$Low2 <- Fig2DF$Low*100
    Fig2DF$High2 <- Fig2DF$High*100 
  }
  rm(h,l,x,y,margin,TempDF)
}
Fig2DF


OverallRemChart <-
  ggplot(data=Fig2DF[Fig2DF$Material=='Overall',],aes(x=Month.Number,y=Avg2,color=Infection))+
  scale_color_manual(values = LevelColors.NoCG,labels=LevelNames.NoCG,name='Litter Source')+
  scale_fill_manual(values = LevelColors.NoCG,labels=LevelNames.NoCG,name='Litter Source')+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  geom_line(size=1)+xlab('')+ylab('Mass Remaining (%)')+coord_cartesian(ylim=c(0,RemChartMax))+scale_y_continuous(breaks=c(0,25,50,75,100,125))+
  scale_x_continuous(breaks=c(2,4,6,8,10,12),labels=c('February','April','June','August','October','December'))+
  theme_bw()+theme(panel.grid = element_blank() )+
  geom_point(data=Chart[Chart$Month.Number<=12&Chart$Month.Number>0&(Chart$Infection!='Common Garden'),], inherit.aes = F,
             aes(x=Month.Number,y=(Prop.Mass.Rem*100),fill=Infection ),shape=21,color='black', size=PointSize,
             position = position_jitterdodge(jitter.width = 0.05 ,dodge.width = 0.4)  )+
  ggtitle("Overall")+theme(plot.title = element_text(size=TitleSize))


LeafRemChart <-
  ggplot(data=Fig2DF[Fig2DF$Material=='Leaf',],aes(x=Month.Number,y=Avg2,color=Infection))+
  scale_color_manual(values = LevelColors.NoCG,labels=LevelNames.NoCG,name='Litter Source')+
  scale_fill_manual(values = LevelColors.NoCG,labels=LevelNames.NoCG,name='Litter Source')+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  geom_line(size=1)+xlab('')+ylab('Mass Remaining (%)')+coord_cartesian(ylim=c(0,RemChartMax))+scale_y_continuous(breaks=c(0,25,50,75,100,125))+
  scale_x_continuous(breaks=c(2,4,6,8),labels=c('February','April','June','August'))+
  geom_segment(aes(x=5.80,xend=6.20,y=123,yend=123),size=1,color='black')+ geom_text(size=7,color='black',aes(x=6,y=125),label='*')+
  theme_bw()+theme(panel.grid = element_blank() )+
  theme(legend.position = "none")+
  geom_point(data=Chart[Chart$Month.Number<9&Chart$Month.Number>0&(Chart$Infection!='Common Garden')&Chart$Prop.LeafMass.Rem<(125/100),], inherit.aes = F,
             aes(x=Month.Number,y=(Prop.LeafMass.Rem*100),fill=Infection ),shape=21,color='black',  size=PointSize,
             position = position_jitterdodge(jitter.width = 0.05 ,dodge.width = 0.4)  )+
  ggtitle("Leaf")+theme(plot.title = element_text(size=TitleSize))


StemRemChart <-
  ggplot(data=Fig2DF[Fig2DF$Material=='Stem',],aes(x=Month.Number,y=Avg2,color=Infection))+
  scale_color_manual(values = LevelColors.NoCG,labels=LevelNames.NoCG,name='Litter Source')+
  scale_fill_manual(values = LevelColors.NoCG,labels=LevelNames.NoCG,name='Litter Source')+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  geom_line(size=1)+xlab('')+ylab('')+coord_cartesian(ylim=c(0,RemChartMax))+scale_y_continuous(breaks=c(0,25,50,75,100,125))+
  scale_x_continuous(breaks=c(2,4,6,8),labels=c('February','April','June','August'))+
  theme_bw()+theme(panel.grid = element_blank() )+
  theme(legend.position = "none")+
  geom_point(data=Chart[Chart$Month.Number<9&Chart$Month.Number>0&(Chart$Infection!='Common Garden'),], inherit.aes = F,
             aes(x=Month.Number,y=(Prop.StemMass.Rem*100),fill=Infection ),shape=21,color='black',  size=PointSize,
             position = position_jitterdodge(jitter.width = 0.05 ,dodge.width = 0.4)  )+
  ggtitle("Stem")+theme(plot.title = element_text(size=TitleSize))


CombinedRemChart <- 
  plot_grid(OverallRemChart + AxisFontSize,
            plot_grid(LeafRemChart,StemRemChart,ncol=2,labels=c("B","C")
                      ),
            labels="A",ncol = 1
  )

CombinedRemChart

# ggsave(paste(sep='',"CombinedRemChart",format(Sys.time(),"%e%b%Y"),'.jpg'),CombinedRemChart,dpi=1000,width=6.5,height=6)
