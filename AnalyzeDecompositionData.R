
##### Load packages #####
library(dplyr)
library(lmerTest)
library(expss)
library(ggpubr)
library(cowplot)
library(ggplot2)

# Set working directory as needed
setwd("C:/Users/brett/Dropbox/Decomposition Experiment/FinalCodeAndData/UsingDataFromEDI")

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

for (x in 1:nrow(Chart)){
  if (x==1){ 
    Chart$Rate.Overall <- NA  # In grams. Rate determined as proportion of original mass lost since previous sampling. Rates are compared to the average of the previous sampling to account for random variance (except when comparing to initial)
    Chart$Rate.Leaf <- NA # In grams. Rate determined as proportion of original mass lost since previous sampling. Rates are compared to the average of the previous sampling to account for random variance (except when comparing to initial)
    Chart$Rate.Stem <- NA # In grams. Rate determined as proportion of original mass lost since previous sampling. Rates are compared to the average of the previous sampling to account for random variance (except when comparing to initial)
  }
  Chart$Rate.Overall[x] <- mean(Chart$Prop.Mass.Rem[Chart$Month.Number==(Chart$Month.Number[x]-2)&Chart$Litter.Source==Chart$Litter.Source[x]])-Chart$Prop.Mass.Rem[x]
}

for (x in 1:nrow(Chart)){
  if (x==1){Chart$Rate.Stem <- NA}
  if(Chart$Month.Number[x]<=2 | Chart$Month.Number[x]>8){
    Chart$Rate.Stem[x] <- NA
  }else{
    PreviousDF <- Chart[Chart$Litter.Source==Chart$Litter.Source[x]&Chart$Month.Number==(Chart$Month.Number[x]-2) ,]
    Chart$Rate.Stem[x] <- mean(PreviousDF$Prop.StemMass.Rem)-Chart$Prop.StemMass.Rem[x]
  }
  rm(x)
}

for (x in 1:nrow(Chart)){
  if (x==1){Chart$Rate.Leaf <- NA}
  if(Chart$Month.Number[x]<=2 | Chart$Month.Number[x]>8){
    Chart$Rate.Leaf[x] <- NA
  }else{
    PreviousDF <- Chart[Chart$Litter.Source==Chart$Litter.Source[x]&Chart$Month.Number==(Chart$Month.Number[x]-2) ,]
    Chart$Rate.Leaf[x] <- mean(PreviousDF$Prop.LeafMass.Rem)-Chart$Prop.LeafMass.Rem[x]
  }
  rm(x)
}

##### Table S3 #####

# Inspection of residuals for overall mass remaining revealed no need for transformation

# For table S3 none of the data needs to be transformed
for (x in Months){
  if (x==Months[1]){ #Initiate data frame
    TableS3 <- data.frame(Tissue=factor(),Month=factor(), # Tissue (Overall, Leaf, Stem), and Month
                          CG.PerRem.mean=numeric(), CG.PerRem.sd=numeric(), # Stats on litter from Common Garden
                          I.PerRem.mean=numeric(), I.PerRem.sd=numeric(), # Stats on litter from Infected sites
                          NI.PerRem.mean=numeric(), NI.PerRem.sd=numeric(), # Stats on litter from Non-infected sites
                          pval.I.NI=numeric(), qval.I.NI=numeric(), # Comparison of litter from Infected and Non-infected sites
                          pval.NI.CG=numeric(), qval.NI.CG=numeric(),
                          Fval.I.NI=numeric(),Fval.NI.CG=numeric() ) # Comparison of litter from Non-infected and Common Garden sites
  }
  ### Start with Overall tissue ###
  MonthChart <- Chart[Chart$Month==x,]
  I.NI.res <- anova(lmer(Prop.Mass.Rem~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Common Garden',])) # Comparison of infected and non-infected
  NI.CG.res <- anova(lmer(Prop.Mass.Rem~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Infected',])) # Comparison of non-infected and common garden
  
  TableS3 <- rbind( TableS3,
    data.frame(Tissue='Overall',Month=x,CG.PerRem.mean=(mean(MonthChart$Prop.Mass.Rem[MonthChart$Infection=='Common Garden'])*100),CG.PerRem.sd=(sd(MonthChart$Prop.Mass.Rem[MonthChart$Infection=='Common Garden'])*100),
             I.PerRem.mean=(mean(MonthChart$Prop.Mass.Rem[MonthChart$Infection=='Infected'])*100), I.PerRem.sd=(sd(MonthChart$Prop.Mass.Rem[MonthChart$Infection=='Infected'])*100),
             NI.PerRem.mean=(mean(MonthChart$Prop.Mass.Rem[MonthChart$Infection=='Non-infected'])*100), NI.PerRem.sd=(sd(MonthChart$Prop.Mass.Rem[MonthChart$Infection=='Non-infected'])*100),
             pval.I.NI=I.NI.res$`Pr(>F)`[1]%>%round(3), qval.I.NI=p.adjust(I.NI.res$`Pr(>F)`[1],'bonferroni',n=6)%>%round(3),
             pval.NI.CG=NI.CG.res$`Pr(>F)`[1]%>%round(3), qval.NI.CG=p.adjust(NI.CG.res$`Pr(>F)`[1],'bonferroni',n=6)%>%round(3),
             Fval.I.NI=round(I.NI.res$`F value`[1],3),Fval.NI.CG=round(NI.CG.res$`F value`[1],3) )
  )
  
  rm(I.NI.res,NI.CG.res)
  
  if(x!="October" & x!="December"){ # Skip Oct and Dec for Leaf and Stem tissue, this data was not collected
    ### Leaf tissue ###
    I.NI.res <- anova(lmer(Prop.LeafMass.Rem~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Common Garden',])) # Comparison of infected and non-infected
    NI.CG.res <- anova(lmer(Prop.LeafMass.Rem~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Infected',])) # Comparison of non-infected and common garden
    TableS3 <- rbind( TableS3,
                      data.frame(Tissue='Leaf',Month=x,CG.PerRem.mean=(mean(MonthChart$Prop.LeafMass.Rem[MonthChart$Infection=='Common Garden'])*100),CG.PerRem.sd=(sd(MonthChart$Prop.LeafMass.Rem[MonthChart$Infection=='Common Garden'])*100),
                                 I.PerRem.mean=(mean(MonthChart$Prop.LeafMass.Rem[MonthChart$Infection=='Infected'])*100), I.PerRem.sd=(sd(MonthChart$Prop.LeafMass.Rem[MonthChart$Infection=='Infected'])*100),
                                 NI.PerRem.mean=(mean(MonthChart$Prop.LeafMass.Rem[MonthChart$Infection=='Non-infected'])*100), NI.PerRem.sd=(sd(MonthChart$Prop.LeafMass.Rem[MonthChart$Infection=='Non-infected'])*100),
                                 pval.I.NI=I.NI.res$`Pr(>F)`[1]%>%round(3), qval.I.NI=p.adjust(I.NI.res$`Pr(>F)`[1],'bonferroni',n=4)%>%round(3),
                                 pval.NI.CG=NI.CG.res$`Pr(>F)`[1]%>%round(3), qval.NI.CG=p.adjust(NI.CG.res$`Pr(>F)`[1],'bonferroni',n=4)%>%round(3),
                                 Fval.I.NI=round(I.NI.res$`F value`[1],3),Fval.NI.CG=round(NI.CG.res$`F value`[1],3) )
    )
    rm(I.NI.res,NI.CG.res)

    ### Stem tissue ###
    I.NI.res <- anova(lmer(Prop.StemMass.Rem~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Common Garden',])) # Comparison of infected and non-infected
    NI.CG.res <- anova(lmer(Prop.StemMass.Rem~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Infected',])) # Comparison of non-infected and common garden
    TableS3 <- rbind( TableS3,
                      data.frame(Tissue='Stem',Month=x,CG.PerRem.mean=(mean(MonthChart$Prop.StemMass.Rem[MonthChart$Infection=='Common Garden'])*100),CG.PerRem.sd=(sd(MonthChart$Prop.StemMass.Rem[MonthChart$Infection=='Common Garden'])*100),
                                 I.PerRem.mean=(mean(MonthChart$Prop.StemMass.Rem[MonthChart$Infection=='Infected'])*100), I.PerRem.sd=(sd(MonthChart$Prop.StemMass.Rem[MonthChart$Infection=='Infected'])*100),
                                 NI.PerRem.mean=(mean(MonthChart$Prop.StemMass.Rem[MonthChart$Infection=='Non-infected'])*100), NI.PerRem.sd=(sd(MonthChart$Prop.StemMass.Rem[MonthChart$Infection=='Non-infected'])*100),
                                 pval.I.NI=I.NI.res$`Pr(>F)`[1]%>%round(3), qval.I.NI=p.adjust(I.NI.res$`Pr(>F)`[1],'bonferroni',n=4)%>%round(3),
                                 pval.NI.CG=NI.CG.res$`Pr(>F)`[1]%>%round(3), qval.NI.CG=p.adjust(NI.CG.res$`Pr(>F)`[1],'bonferroni',n=4)%>%round(3),
                                 Fval.I.NI=round(I.NI.res$`F value`[1],3),Fval.NI.CG=round(NI.CG.res$`F value`[1],3) )
    )
    rm(I.NI.res,NI.CG.res)
  } 
  rm(x, MonthChart)
}


##### Table S4 #####

#### Will need to transform some rates here #####
TempDF <- Chart[Chart$Infection!='Common Garden',]
TempDF$Rate.Stem[TempDF$Month.Number== 4]  %>% shapiro.test() # 0.001
TempDF$Rate.Stem[TempDF$Month.Number== 6]  %>% shapiro.test() # 0.003
TempDF$Rate.Stem[TempDF$Month.Number== 8]  %>% shapiro.test() # normal

for (x in c(4,6,8,10,12)){
  if (x==4){ #Initiate data frame
    TableS4 <- data.frame(Tissue=factor(),Month=factor(), # Tissue (Overall, Leaf, Stem), and Month interval
                          I.Rate.mean=numeric(), I.Rate.sd=numeric(), # Stats on litter rates from Infected sites
                          NI.Rate.mean=numeric(), NI.Rate.sd=numeric(), # Stats on litter rates from Non-infected sites
                          pval.I.NI=numeric(), qval.I.NI=numeric(),Fval.I.NI=numeric() ) # Comparison of litter rates from Infected and Non-infected sites
  }
  if(x==6 | x==8 | x==12){
    if(x==6){
      result <- anova(lmer(log(Rate.Overall)~Infection+(1|Replicate/Infection),Chart[Chart$Month.Number==x&Chart$Infection!='Common Garden',]))
    }else{
      result <- anova(lmer(exp(Rate.Overall)~Infection+(1|Replicate/Infection),Chart[Chart$Month.Number==x&Chart$Infection!='Common Garden',]))
    }
  }else{
    result <- anova(lmer(Rate.Overall~Infection+(1|Replicate/Infection),Chart[Chart$Month.Number==x&Chart$Infection!='Common Garden',]))
  }
  qval <- p.adjust(result$`Pr(>F)`[1],'bonferroni',n=5)%>%round(3)
  TableS4 <- rbind(TableS4,
  data.frame(Tissue='Overall',Month=paste(Chart$Month[Chart$Month.Number==(x-2)][1],Chart$Month[Chart$Month.Number==x][1],sep='-'),
             I.Rate.mean=round((mean(Chart$Rate.Overall[Chart$Month.Number==x&Chart$Infection=='Infected'])*100),digits=1), I.Rate.sd=round((sd(Chart$Rate.Overall[Chart$Month.Number==x&Chart$Infection=='Infected'])*100),digits=2),
             NI.Rate.mean=round((mean(Chart$Rate.Overall[Chart$Month.Number==x&Chart$Infection=='Non-infected'])*100),digits=1), NI.Rate.sd=round((sd(Chart$Rate.Overall[Chart$Month.Number==x&Chart$Infection=='Non-infected'])*100),digits=2), 
             pval.I.NI=result$`Pr(>F)`[1]%>%round(3), qval.I.NI=qval,Fval.I.NI=result$`F value`[1]%>%round(3) )
  )
  rm(result,qval)
  if(x<=8){
    ## Leaf tissue ##
    # No transformations necessary here
    result <- anova(lmer(Rate.Leaf~Infection+(1|Replicate/Infection),Chart[Chart$Month.Number==x&Chart$Infection!='Common Garden',]))
    qval <- p.adjust(result$`Pr(>F)`[1],'bonferroni',n=3)%>%round(3)
    TableS4 <- rbind(TableS4,
    data.frame(Tissue='Leaf',Month=paste(Chart$Month[Chart$Month.Number==(x-2)][1],Chart$Month[Chart$Month.Number==x][1],sep='-'),
               I.Rate.mean=round((mean(Chart$Rate.Leaf[Chart$Month.Number==x&Chart$Infection=='Infected'])*100),digits=1), I.Rate.sd=round((sd(Chart$Rate.Leaf[Chart$Month.Number==x&Chart$Infection=='Infected'])*100),digits=2),
               NI.Rate.mean=round((mean(Chart$Rate.Leaf[Chart$Month.Number==x&Chart$Infection=='Non-infected'])*100),digits=1), NI.Rate.sd=round((sd(Chart$Rate.Leaf[Chart$Month.Number==x&Chart$Infection=='Non-infected'])*100),digits=2), 
               pval.I.NI=result$`Pr(>F)`[1]%>%round(3), qval.I.NI=qval,Fval.I.NI=result$`F value`[1]%>%round(3) )
    )
    rm(result,qval)
    
    ## Stem tissue ##
    result <- anova(lmer(Rate.Stem~Infection+(1|Replicate/Infection),Chart[Chart$Month.Number==x&Chart$Infection!='Common Garden',]))
    qval <- p.adjust(result$`Pr(>F)`[1],'bonferroni',n=3)%>%round(3)
    TableS4 <- rbind(TableS4,
    data.frame(Tissue='Stem',Month=paste(Chart$Month[Chart$Month.Number==(x-2)][1],Chart$Month[Chart$Month.Number==x][1],sep='-'),
               I.Rate.mean=round((mean(Chart$Rate.Stem[Chart$Month.Number==x&Chart$Infection=='Infected'])*100),digits=1), I.Rate.sd=round((sd(Chart$Rate.Stem[Chart$Month.Number==x&Chart$Infection=='Infected'])*100),digits=2),
               NI.Rate.mean=round((mean(Chart$Rate.Stem[Chart$Month.Number==x&Chart$Infection=='Non-infected'])*100),digits=1), NI.Rate.sd=round((sd(Chart$Rate.Stem[Chart$Month.Number==x&Chart$Infection=='Non-infected'])*100),digits=2), 
               pval.I.NI=result$`Pr(>F)`[1]%>%round(3), qval.I.NI=qval,Fval.I.NI=result$`F value`[1]%>%round(3) )
    )
    rm(result,qval)
    }
  rm(x)
}

##### Table S5 - Litter Chemistry #####

for (x in c(0,2,4,6,8,10,12)){
  if (x == 0 ){
    TableS5 <- data.frame(Level=vector(),Month=vector(),CG.mean=numeric(),CG.sd=numeric(),NI.mean=numeric(),NI.sd=numeric(),I.mean=numeric(),I.sd=numeric(),CG.NI.p=numeric(),CG.NI.q=numeric(),I.NI.p=numeric(),I.NI.q=numeric(),CG.NI.f=numeric(),I.NI.f=numeric() )
  }
  MonthChart <- Chart[Chart$Month.Number==x & (!is.na(Chart$pC)),]
  Month <- MonthChart$Month[1]
  CG.NI <- anova(lmer(ADLom~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Infected',]))
  I.NI <- anova(lmer(ADLom~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Common Garden',]))
  CG.NI.q <- p.adjust(CG.NI$`Pr(>F)`[1],method = 'bonferroni',n=7) %>% round(3)
  I.NI.q <- p.adjust(I.NI$`Pr(>F)`[1],method = 'bonferroni',n=7) %>% round(3)
  TableS5 <- rbind(TableS5,
    data.frame(Level='Lignin',Month=Month,CG.mean=round(mean(MonthChart$ADLom[MonthChart$Infection=='Common Garden']),digits=2),CG.sd=round(sd(MonthChart$ADLom[MonthChart$Infection=='Common Garden']),digits=2),
             NI.mean=round(mean(MonthChart$ADLom[MonthChart$Infection=='Non-infected']),digits=2),NI.sd=round(sd(MonthChart$ADLom[MonthChart$Infection=='Non-infected']),digits=2),
             I.mean=round(mean(MonthChart$ADLom[MonthChart$Infection=='Infected']),digits=2),I.sd=round(sd(MonthChart$ADLom[MonthChart$Infection=='Infected']),digits=2),
             CG.NI.p=CG.NI$`Pr(>F)`[1]%>% round(3),CG.NI.q=CG.NI.q,
             I.NI.p=I.NI$`Pr(>F)`[1]%>% round(3),I.NI.q=I.NI.q,
             CG.NI.f=CG.NI$`F value`[1]%>%round(3),I.NI.f=I.NI$`F value`[1]%>%round(3) )
  )
  
  CG.NI <- anova(lmer(CN~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Infected',]))
  I.NI <- anova(lmer(CN~Infection+(1|Replicate/Infection),MonthChart[MonthChart$Infection!='Common Garden',]))
  CG.NI.q <- p.adjust(CG.NI$`Pr(>F)`[1],method = 'bonferroni',n=7) %>% round(3)
  I.NI.q <- p.adjust(I.NI$`Pr(>F)`[1],method = 'bonferroni',n=7) %>% round(3)
  TableS5 <- rbind(TableS5,
                   data.frame(Level='CN',Month=Month,CG.mean=round(mean(MonthChart$CN[MonthChart$Infection=='Common Garden']),digits=2),CG.sd=round(sd(MonthChart$CN[MonthChart$Infection=='Common Garden']),digits=2),
                              NI.mean=round(mean(MonthChart$CN[MonthChart$Infection=='Non-infected']),digits=2),NI.sd=round(sd(MonthChart$CN[MonthChart$Infection=='Non-infected']),digits=2),
                              I.mean=round(mean(MonthChart$CN[MonthChart$Infection=='Infected']),digits=2),I.sd=round(sd(MonthChart$CN[MonthChart$Infection=='Infected']),digits=2),
                              CG.NI.p=CG.NI$`Pr(>F)`[1]%>% round(3),CG.NI.q=CG.NI.q,
                              I.NI.p=I.NI$`Pr(>F)`[1]%>% round(3),I.NI.q=I.NI.q,
                              CG.NI.f=CG.NI$`F value`[1]%>%round(3),I.NI.f=I.NI$`F value`[1]%>%round(3) )
  )
  rm(I.NI,I.NI.q,CG.NI,CG.NI.q,MonthChart,Month,x)
}

TableS5
TableS5[TableS5$Level=='Lignin',] %>% select(Level,Month,I.NI.f,I.NI.p,I.NI.q)
TableS5[TableS5$Level=='CN',]     %>% select(Level,Month,I.NI.f,I.NI.p,I.NI.q)


TableS5[TableS5$Level=='Lignin',] %>% select(Level,Month,I.NI.f)
TableS5[TableS5$Level=='CN',]     %>% select(Level,Month,I.NI.f)




TableS3[TableS3$Tissue=='Overall',]
TableS3[TableS3$Tissue=='Leaf',]
TableS3[TableS3$Tissue=='Stem',]

TableS4[TableS4$Tissue=='Overall',]
TableS4[TableS4$Tissue=='Leaf',]
TableS4[TableS4$Tissue=='Stem',]

##### Persistent variables for figures #####
LevelOrder <- c('Infected','Non-infected','Common Garden')
LevelNames <- c('Non-infected'="Non-infected\nSites", 'Infected'="Infected\nSites",'Common Garden'="Common Garden\nSites")
LevelColors <- c('Non-infected'="#619CFF", 'Infected'="#F8766D",'Common Garden'="#00BA38")
LevelOrder.NoCG <- c('Infected','Non-infected')
LevelColors.NoCG <- LevelColors[names(LevelColors)!='Common Garden']
LevelNames.NoCG <- LevelNames[names(LevelNames)!='Common Garden']
RemChartMax <- 125
PointSize <- 0.7
TitleSize <- 10
CG.RemChartMax <- 125
Overall.CG.RemChartMax <- 125
CG.PointSize <- .75
CGAdj <- 5
LevelOrder.CG <- LevelOrder[LevelOrder!='Infected']
LevelColors.CG <- LevelColors[names(LevelColors)!='Infected']
LevelNames.CG <- c('Non-infected'="Non-infected\nSites", 'Common Garden'="Experimental Site")
AxisFontSize <- theme(axis.text.x = element_text(size=8,angle=0), legend.text=element_text(size=7))


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
  ggarrange(OverallRemChart + AxisFontSize,
            ggarrange(LeafRemChart,StemRemChart,ncol=2,labels=c("B","C")),
            nrow=2, heights = c(3,5),
            labels = "A"
  )
CombinedRemChart
# ggsave(paste(sep='',"CombinedRemChart",format(Sys.time(),"%e%b%Y"),'.jpg'),CombinedRemChart,dpi=1000,width=6.5,height=6)

##### Figure 3 #####

OverallDotRateChart.NoCG <-
ggplot(Chart[Chart$Month.Number>2&Chart$Infection!='Common Garden',],
                                   aes(x=factor(Month.Number), y= Rate.Overall*100, fill = factor(Infection, level=LevelOrder)))+
  coord_cartesian(ylim=c(-60,60))+ geom_hline(yintercept = 0)+ geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0.3),aes(group=factor(Infection, level=LevelOrder) ),color="#444444",size=.75)+
  theme_bw()+theme(panel.grid=element_blank())+
  scale_fill_manual("Litter Source",values=LevelColors.NoCG,  labels = LevelNames.NoCG )+
  scale_x_discrete(breaks = c(4,6,8,10,12), labels = c("February-April","April-June","June-August","August-October","October-December"))+
  geom_boxplot(outlier.shape = NA,alpha=0)+
  # geom_segment(aes(x=2.875,xend=3.125,y=42.5,yend=42.5),size=1,color='red3')+ geom_text(aes(x=3,y=43.5),label='*',color='red3',size=8)+
  xlab("") + ylab("Mass lost bimonthly (%)") + ggtitle("Overall")+theme(plot.title = element_text(size=TitleSize))

LeafDotRateChart.NoCG <-
ggplot(Chart[Chart$Month.Number>2&Chart$Month.Number<10&Chart$Infection!='Common Garden',],
                                aes(x=factor(Month.Number), y= Rate.Leaf*100, fill = factor(Infection, level=LevelOrder)))+
  geom_hline(yintercept = 0)+
  geom_boxplot(outlier.shape = NA)+ 
  geom_point(position=position_jitterdodge(jitter.width = 0.3),aes(group=factor(Infection, level=LevelOrder)),color="#444444",size=.75 )+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  theme(legend.position = "none")+
  scale_fill_manual("Litter Source",values=LevelColors,  labels = LevelNames )+
  scale_x_discrete(breaks = c(4,6,8,10,12), labels = c("February-April","April-June","June-August","\nAugust-October","October-December"))+
  xlab("") + ggtitle("Leaf")+ ylab("Mass Lost Bimonthly (%)")+
  scale_y_continuous(breaks = c(75,50,25,0,-25,-50,-75))+ coord_cartesian(ylim=c(-60,85))+
  geom_segment(aes(x=2.875,xend=3.125,y=70.5,yend=70.5),size=1,color='red3')+ geom_text(aes(x=3,y=71.5),label='*',color='red3',size=8)+
  geom_boxplot(outlier.shape = NA,alpha=0)+
  theme(plot.title = element_text(size=TitleSize))


StemDotRateChart.NoCG <-
ggplot(Chart[Chart$Month.Number>2&Chart$Month.Number<10&Chart$Infection!='Common Garden',],
                                aes(x=factor(Month.Number), y= Rate.Stem*100, fill = factor(Infection, level=LevelOrder)))+
  geom_hline(yintercept = 0)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0.3),aes(group=factor(Infection, level=LevelOrder)),color="#444444",size=.75 )+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  theme(legend.position = "none")+
  scale_fill_manual("Litter Source",values=LevelColors,  labels = LevelNames )+
  scale_x_discrete(breaks = c(4,6,8,10,12), labels = c("February-April","April-June","June-August","\nAugust-October","October-December"))+
  scale_y_continuous(breaks = c(75,50,25,0,-25,-50,-75))+ coord_cartesian(ylim=c(-60,85))+
  geom_boxplot(outlier.shape = NA,alpha=0)+
  xlab("") + ggtitle("Stem")+ylab('')+theme(plot.title = element_text(size=TitleSize))

Fig3RateDots.NoCG <-
  ggarrange(OverallDotRateChart.NoCG + ylab("Mass lost bimonthly (%)") + AxisFontSize,
            ggarrange(LeafDotRateChart.NoCG+ ylab("Mass lost bimonthly (%)")
                      ,StemDotRateChart.NoCG,ncol=2,labels=c("B","C")),
            nrow=2, heights = c(4,5),
            labels = "A"
  )
Fig3RateDots.NoCG

# ggsave(paste(sep='',"Fig3RateWithDotsAndNoCG",format(Sys.time(),"%e%b%Y"),'.jpg'),Fig3RateDots.NoCG,dpi=1000, width=7,height=8)

##### Figure 4 ####

BoxCN <-
ggplot(Chart[Chart$Infection!='Common Garden'&!is.na(Chart$CN),],aes(x=as.factor(Month.Number),y=CN,fill=Infection)  )+
  theme_bw()+ theme(panel.grid=element_blank() )+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width=0.15),color="#444444",size=0.75)+
  scale_x_discrete(breaks=c(0,2,4,6,8,10,12),labels=c('Start','February','April','June','August','October','December'))+
  scale_fill_manual( values=c('Infected'='#F8766D',"Non-infected"='#619CFF'), labels=c('Infected'='Infected sites','Non-infected'='Non-infected sites'), name='Litter source'  )+
  coord_cartesian(ylim=c(15,130))+ scale_y_continuous(breaks=c(25,50,75,100,125))+
  geom_boxplot(outlier.shape = NA,aes(x=as.factor(Month.Number),y=CN,fill=Infection),alpha=0)+
  theme(legend.position="none")+
  xlab("Month")+ xlab('')+ ylab('Carbon/Nitrogen ratio')


BoxADL <-
ggplot(Chart[Chart$Infection!='Common Garden'&!is.na(Chart$ADLom),],aes(x=as.factor(Month.Number),y=ADLom,fill=Infection)  )+
  theme_bw()+ theme(panel.grid=element_blank() )+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width=0.15),color="#444444",size=0.75)+
  scale_x_discrete(breaks=c(0,2,4,6,8,10,12),labels=c('Start','February','April','June','August','October','December'))+
  scale_fill_manual( values=c('Infected'='#F8766D',"Non-infected"='#619CFF'), labels=c('Infected'='Infected Sites','Healthy'='Non-infected Sites'), name='Litter Source'  )+
  coord_cartesian(ylim=c(0,26))+scale_y_continuous(breaks=c(0,5,10,15,20,25))+
  geom_boxplot(outlier.shape = NA,aes(x=as.factor(Month.Number),y=ADLom,fill=Infection),alpha=0)+
  theme(legend.position="none")+
  ylab('Acid detergent lignin\n(as percent of organic mass)')+  xlab("")

(ADL.CN.Chart <- ggarrange(BoxADL, BoxCN,ncol = 1, nrow = 2,common.legend = T,legend='bottom',labels=c("A",'B')))

# ggsave(paste(sep='',"ADL.CN.Chart.WithGreyPoints",format(Sys.time(),"%e%b%Y"),'.jpg'),ADL.CN.Chart,dpi=1000,height=7,width=7)

##### Figure S2 #####

# Down around 970 in my code

for (x in c(2,4,6,8,10,12) ){
  if(x==2){
    F2S <- data.frame(Month.Number=vector(),Material=vector(),Infection=vector(),Avg=vector(),Low=vector(),High=vector())
  }
  for (y in LevelOrder.CG){
    # break
    TempDF <- Chart[Chart$Month.Number==x&Chart$Infection==y,]
    # head(TempDF)
    margin <- (qt(0.975,df=(nrow(TempDF)-1))*sd(TempDF$Prop.Mass.Rem)/sqrt(nrow(TempDF)))
    l <- mean(TempDF$Prop.Mass.Rem)-margin
    h <- mean(TempDF$Prop.Mass.Rem)+margin
    F2S <- rbind(F2S, data.frame(Month.Number=x,Material="Overall",Infection=y,Avg=mean(TempDF$Prop.Mass.Rem),Low=l,High=h) )
    # next
    if(x<9){
      margin <- (qt(0.975,df=(nrow(TempDF)-1))*sd(TempDF$Prop.LeafMass.Rem)/sqrt(nrow(TempDF)))
      l <- mean(TempDF$Prop.LeafMass.Rem)-margin
      h <- mean(TempDF$Prop.LeafMass.Rem)+margin
      F2S <- rbind(F2S, data.frame(Month.Number=x,Material="Leaf",Infection=y,Avg=mean(TempDF$Prop.LeafMass.Rem ) ,Low=l, High=h )  )
      
      margin <- (qt(0.975,df=(nrow(TempDF)-1))*sd(TempDF$Prop.StemMass.Rem)/sqrt(nrow(TempDF)))
      l <- mean(TempDF$Prop.StemMass.Rem)-margin
      h <- mean(TempDF$Prop.StemMass.Rem)+margin
      F2S <- rbind(F2S, data.frame(Month.Number=x,Material="Stem",Infection=y,Avg=mean(TempDF$Prop.StemMass.Rem ) ,Low=l, High=h )  )
    }
  }
  if(x==12){
    F2S$Avg2 <- F2S$Avg*100
    F2S$Low2 <- F2S$Low*100
    F2S$High2 <- F2S$High*100 
  }
  rm(x,y,margin,l,h,TempDF)
}

head(F2S)

CG.OverallRemChart <-
ggplot(data=F2S[F2S$Material=='Overall',],aes(x=Month.Number,y=Avg2,color=Infection))+
  scale_color_manual(values = LevelColors.CG,labels=LevelNames.CG,name='Litter Source')+
  scale_fill_manual(values = LevelColors.CG,labels=LevelNames.CG,name='Litter Source')+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  geom_line(size=1)+xlab('')+ylab('Mass Remaining (%)')+coord_cartesian(ylim=c(0,Overall.CG.RemChartMax))+
  scale_y_continuous(breaks=c(0,25,50,75,100,125))+
  scale_x_continuous(breaks=c(2,4,6,8,10,12),labels=c('February','April','June','August','October','December'))+
  theme_bw()+theme(panel.grid = element_blank() )+
  geom_segment(aes(x=5.80,xend=6.20,y=123-CGAdj,yend=123-CGAdj),size=1,color='black')+ geom_text(size=7,color='black',aes(x=6,y=125-CGAdj),label='*')+
  geom_segment(aes(x=3.80,xend=4.20,y=123-CGAdj,yend=123-CGAdj),size=1,color='black')+ geom_text(size=7,color='black',aes(x=4,y=125-CGAdj),label='*')+
  geom_point(data=Chart[Chart$Month.Number<=12&Chart$Month.Number>0&(Chart$Infection!='Infected'),], inherit.aes = F,
             aes(x=Month.Number,y=(Prop.Mass.Rem*100),fill=Infection ),shape=21,color='black', size=CG.PointSize,
             position = position_jitterdodge(jitter.width = 0.05 ,dodge.width = -0.4)  )+
  ggtitle("Overall")+theme(plot.title = element_text(size=TitleSize))



CG.LeafRemChart <-
ggplot(data=F2S[F2S$Material=='Leaf',],aes(x=Month.Number,y=Avg2,color=Infection))+
  scale_color_manual(values = LevelColors.CG,labels=LevelNames.CG,name='Litter Source')+
  scale_fill_manual(values = LevelColors.CG,labels=LevelNames.CG,name='Litter Source')+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  geom_line(size=1)+xlab('')+ylab('Mass Remaining (%)')+coord_cartesian(ylim=c(0,CG.RemChartMax))+scale_y_continuous(breaks=c(0,25,50,75,100,125))+
  scale_x_continuous(breaks=c(2,4,6,8),labels=c('February','April','June','August'))+
  geom_segment(aes(x=5.80,xend=6.20,y=123-CGAdj,yend=123-CGAdj),size=1,color='black')+ geom_text(size=7,color='black',aes(x=6,y=125-CGAdj),label='*')+
  geom_segment(aes(x=7.80,xend=8.20,y=123-CGAdj,yend=123-CGAdj),size=1,color='black')+ geom_text(size=7,color='black',aes(x=8,y=125-CGAdj),label='*')+
  theme_bw()+theme(panel.grid = element_blank() )+
  theme(legend.position = "none")+
  geom_point(data=Chart[Chart$Month.Number<9&Chart$Month.Number>0&(Chart$Infection!='Infected')&Chart$Prop.LeafMass.Rem,], inherit.aes = F,
             aes(x=Month.Number,y=(Prop.LeafMass.Rem*100),fill=Infection ),shape=21,color='black',  size=CG.PointSize,
             position = position_jitterdodge(jitter.width = 0.05 ,dodge.width = -0.4)  )+
  ggtitle("Leaf")+theme(plot.title = element_text(size=TitleSize))



CG.StemRemChart <-
ggplot(data=F2S[F2S$Material=='Stem',],aes(x=Month.Number,y=Avg2,color=Infection))+
  scale_color_manual(values = LevelColors.CG,labels=LevelNames.CG,name='Litter Source')+
  scale_fill_manual(values = LevelColors.CG,labels=LevelNames.CG,name='Litter Source')+
  geom_ribbon(aes(ymin=Low2,ymax=High2, fill=Infection ),alpha=0.25   )+
  geom_line(size=1)+xlab('')+ylab('')+coord_cartesian(ylim=c(0,CG.RemChartMax))+scale_y_continuous(breaks=c(0,25,50,75,100,125))+
  scale_x_continuous(breaks=c(2,4,6,8),labels=c('February','April','June','August'))+
  theme_bw()+theme(panel.grid = element_blank() )+
  theme(legend.position = "none")+
  geom_point(data=Chart[Chart$Month.Number<9&Chart$Month.Number>0&(Chart$Infection!='Infected'),], inherit.aes = F,
             aes(x=Month.Number,y=(Prop.StemMass.Rem*100),fill=Infection ),shape=21,color='black',  size=CG.PointSize,
             position = position_jitterdodge(jitter.width = 0.05 ,dodge.width = -0.4)  )+
  ggtitle("Stem")+theme(plot.title = element_text(size=TitleSize))


CG.CombinedRemChart <- 
  ggarrange(CG.OverallRemChart + AxisFontSize,
            ggarrange(CG.LeafRemChart,CG.StemRemChart,ncol=2,labels=c("B","C")),
            nrow=2, heights = c(3,5),
            labels = "A"  )
CG.CombinedRemChart

# ggsave(paste(sep='',"CG.CombinedRemChart.ForSupp",format(Sys.time(),"%e%b%Y"),'.jpg'),CG.CombinedRemChart,dpi=1000,width=6.5,height=6)

