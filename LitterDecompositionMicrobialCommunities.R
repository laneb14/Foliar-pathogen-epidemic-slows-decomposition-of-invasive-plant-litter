
##### Load Packages #####
library(tidyr)
library(dplyr)
library(phyloseq)
library(decontam)
library(vegan)
library(expss)
library(pairwiseAdonis)
library(ape)
library(ggplot2)
library(ggResidpanel)

# Note: Throughout I use Shannon's Evenness, or "SEven" in place of Pielou's. These metrics are the same, but named differently.

# Set working directory as needed
setwd("C:/Users/brett/Dropbox/Decomposition Experiment/FinalCodeAndData/UsingDataFromEDI")
DataFolder <- "edi.1434.1"

##### Persistent Variables #####
MonthOrder <- c('Start','February','April','June','August','October','December')
Lvalues <- c("Non-infected"="Non-infected Sites","Infected"='Infected Sites')

##### Load data #####
InputMeta <- read.table(file = paste(sep='/',DataFolder, 'MetadataMicrostegiumDecomposition2020Sequencing.tsv'), sep = '\t', header = TRUE) # Load metadata
InputOTUTable.LongFormat <- read.table(file = paste(sep='/',DataFolder, 'OTU_Table_InLongFormat.tsv'), sep = '\t', header = TRUE) # Load otu table
# OTU Table had to be uploaded to data base in long format, but phyloseq will want it in wide format
InputOTUTable <-spread(InputOTUTable.LongFormat,key='Sample',value='Reads')
rownames(InputOTUTable) <- InputOTUTable$ASV
InputOTUTable$ASV <- NULL
InputOTUTable <- as.matrix(InputOTUTable) # phyloseq wants it as a matrix
head(InputOTUTable)
colnames(InputOTUTable) <- gsub("\\.","-",colnames(InputOTUTable)) # May be needed but not always. Sometimes the dashes get replaced by dots when bringing the data into R. This puts it back if needed
rownames(InputMeta) <- InputMeta$SampleName # Copy sample names into the row names
OriginalPhyloObject <- phyloseq(otu_table(InputOTUTable,taxa_are_rows = T),sample_data(InputMeta)) # Create the phyloseq object

##### Run decontam package #####
# Note, we don't have the data for frequency
contamdf.prev <- isContaminant(OriginalPhyloObject, method="prevalence", neg="Negative") # run decontam package; note, will have a warning due to a negative control with no sequences
table(contamdf.prev$contaminant) # No contaminants found
contamdf.prev[contamdf.prev$contaminant==T,] # Confirm no contaminating taxa found


##### Format data #####
# Since no contaminates are found we can just simply remove the negative controls from further analyses
PhyloObject <- subset_samples(OriginalPhyloObject,Negative==0)
# Also remove taxa which were only present in the negative controls
PhyloObject <- prune_taxa(taxa_sums(PhyloObject)>0,PhyloObject)
# Finally, create a phyloseq object based on relative proportion of reads, this will be used later
EvenObject <- transform_sample_counts(PhyloObject,function(x) (x/sum(x))*100 )

##### Confirm adequate sequencing depth #####
# Adequate sequencing depth confirmed through manual inspection of rarefaction plots formed from original data (not transformed)
# Viewing them in sets of 25 to facilitate inspection
# All curves should approach they asysmptote
if(1==2){ # These take a while to run and only need to be viewed once. Add this step so they have to be viewed manually
  rarecurve(t(otu_table(PhyloObject))[1:25,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[26:50,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[51:75,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[76:100,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[101:125,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[126:150,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[151:175,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[176:200,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[201:225,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[226:250,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[251:275,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[276:300,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[301:325,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[326:350,], step=1000, cex=0.5)
  rarecurve(t(otu_table(PhyloObject))[351:368,], step=1000, cex=0.5)
}
##### Calculate alpha diversity #####
# Alpha diversity metrics are calculated using data rarified to an equal sampling depth
set.seed(55735)
RarePhylo <- rarefy_even_depth(PhyloObject)
Indicies <- estimate_richness(RarePhylo,measures=c('Observed','Chao1','Shannon')) # Calculate Chao1 and Shannon's Diversity
rownames(Indicies) <- gsub("\\.","-",rownames(Indicies)) # I'm not sure why it replaces dashes with periods, but this fixes it
Indicies$SEven <- Indicies$Shannon/log(Indicies$Observed) # Calculate Shannon's Evenness
Indicies$Infection <- vlookup(rownames(Indicies),InputMeta,'Infection','SampleName')
Indicies$Month <- vlookup(rownames(Indicies),InputMeta,'Month','SampleName')
Indicies$Replicate <- vlookup(rownames(Indicies),InputMeta,'Replicate','SampleName')

##### Table S4, right three columns ####

for (z in 1:3){
    if(z==1){print('Status,Month,Chao,Shannon,Pielou')}
  for (m in c('Start','February','April','June','August','October','December')){
    if(z==1)(TempDF <- Indicies[Indicies$Month==m & Indicies$Infection=='Infected',])
    if(z==2)(TempDF <- Indicies[Indicies$Month==m & Indicies$Infection=='Non-infected',])
    if(z==3)(TempDF <- Indicies[Indicies$Month==m & Indicies$Infection=='Common Garden',])
    if(nrow(TempDF)==0){stop('no rows found')}
    a <- round(mean(TempDF$Chao1),2)
    b <- round(sd(TempDF$Chao1)/sqrt(length(TempDF$Chao1)),2)
    chao1 <- paste(sep='',a,' (',b,')')
    a <- round(mean(TempDF$Shannon),2)
    b <- round(sd(TempDF$Shannon)/sqrt(length(TempDF$Shannon)),3)
    Shan <- paste(sep='',a,' (',b,')')
    a <- round(mean(TempDF$SEven),2)
    b <- round(sd(TempDF$SEven)/sqrt(length(TempDF$SEven)),3)
    Piel <- paste(sep='',a,' (',b,')')
    print(paste(TempDF$Infection[1],TempDF$Month[1],chao1,Shan,Piel,sep=','))
    rm(a,b,chao1,Shan,Piel,m)
  }
  rm(z)
}

##### Table S5 #####

## Compare Infected and Non Infected

hist(Indicies$Chao1[Indicies$Infection!='Common Garden'])
shapiro.test(Indicies$Chao1[Indicies$Infection!='Common Garden'])
hist(sqrt(Indicies$Chao1[Indicies$Infection!='Common Garden']))
shapiro.test(sqrt(Indicies$Chao1[Indicies$Infection!='Common Garden']))
Chao1.NI.I <- lmer(sqrt(Chao1)~Infection*Month+(1|Replicate),Indicies[Indicies$Infection!='Common Garden',])
resid_panel(Chao1.NI.I)
anova(Chao1.NI.I)


hist(Indicies$Shannon[Indicies$Infection!='Common Garden'])
shapiro.test(Indicies$Shannon[Indicies$Infection!='Common Garden'])
hist((Indicies$Shannon[Indicies$Infection!='Common Garden'])^(3))
shapiro.test((Indicies$Shannon[Indicies$Infection!='Common Garden'])^(3))
Shannon.NI.I <- lmer((Shannon)^3~Infection*Month+(1|Replicate),Indicies[Indicies$Infection!='Common Garden',])
resid_panel(Shannon.NI.I)
anova(Shannon.NI.I)

hist(Indicies$SEven[Indicies$Infection!='Common Garden'])
shapiro.test(Indicies$SEven[Indicies$Infection!='Common Garden'])
hist((Indicies$SEven[Indicies$Infection!='Common Garden'])^3)
shapiro.test((Indicies$SEven[Indicies$Infection!='Common Garden'])^3)
Pielou.NI.I<- lmer((SEven)^3~Infection*Month+(1|Replicate),Indicies[Indicies$Infection!='Common Garden',])
resid_panel(Pielou.NI.I)
anova(Pielou.NI.I)

## Same as above but compare Common Garden with Non-Infected

CG.NI.DF <- Indicies[Indicies$Infection!='Infected',]

hist(CG.NI.DF$Chao1)
shapiro.test(CG.NI.DF$Chao1)
Chao1.CG <- lm(Chao1~Infection*Month, CG.NI.DF)
resid_panel(Chao1.CG)
anova(Chao1.CG)


hist(CG.NI.DF$Shannon)
shapiro.test(CG.NI.DF$Shannon)
hist(CG.NI.DF$Shannon^2)
shapiro.test(CG.NI.DF$Shannon^2)
CG.NI.DF$Shannon.sq <- CG.NI.DF$Shannon^2
Shannon.CG <- lm(Shannon.sq~Infection*Month, CG.NI.DF)
resid_panel(Shannon.CG)
anova(Shannon.CG)


hist(CG.NI.DF$SEven)
shapiro.test(CG.NI.DF$SEven)
hist(CG.NI.DF$SEven^2)
shapiro.test(CG.NI.DF$SEven^2)
CG.NI.DF$SEven.sq <- CG.NI.DF$SEven^2
Pielous.CG <- lm(SEven.sq~Infection*Month, CG.NI.DF)
resid_panel(Pielous.CG)
anova(Pielous.CG)

##### Table S6, Pairwise permanova #####

for (x in MonthOrder){
  head(sample_data(EvenObject))
  TempObject <- subset_samples(EvenObject, Infection!='Common Garden' & Month==x)
  TempDist <- as(vegdist(t(otu_table(TempObject)),method='bray'),'matrix')
  TempDF <- as(sample_data(TempObject),'data.frame')
  Res <- adonis2(TempDist~Infection,data=TempDF,permutations = 10000)
  y <- Res$`Pr(>F)`[1]
  y2 <- ifelse(y<=0.001,"p<0.001",y)
  z <- ifelse(y<=0.001,"p<0.001",p.adjust(y,'bonferroni',n=7) )
  print(paste(sep=',','Non vs Infected','Bray-Curtis',x,Res$F[1],y2,z))
  rm(TempDist,y,y2,z,Res)
  TempDist <- as(vegdist(t(otu_table(TempObject)),method='jaccard'),'matrix')
  Res <- adonis2(TempDist~Infection,data=TempDF,permutations = 10000)
  y <- Res$`Pr(>F)`[1]
  y2 <- ifelse(y<=0.001,"p<0.001",y)
  z <- ifelse(y<=0.001,"p<0.001",p.adjust(y,'bonferroni',n=7) )
  print(paste(sep=',','Non vs Infected','Jaccard',x,Res$F[1],y2,z))
  rm(x,y,y2,z,Res,TempDist,TempDF,TempObject)
  # Can copy and paste the results into excel for manipulation, comma delimited
}

for (x in MonthOrder){
  head(sample_data(EvenObject))
  TempObject <- subset_samples(EvenObject, Infection!='Infected' & Month==x)
  TempDist <- as(vegdist(t(otu_table(TempObject)),method='bray'),'matrix')
  TempDF <- as(sample_data(TempObject),'data.frame')
  Res <- adonis2(TempDist~Infection,data=TempDF,permutations = 10000)
  y <- Res$`Pr(>F)`[1]
  y2 <- ifelse(y<=0.001,"p<0.001",y)
  z <- ifelse(y<=0.001,"p<0.001",p.adjust(y,'bonferroni',n=7) )
  print(paste(sep=',','Common Garden vs Non','Bray-Curtis',x,Res$F[1],y2,z))
  rm(TempDist,y,y2,z,Res)
  TempDist <- as(vegdist(t(otu_table(TempObject)),method='jaccard'),'matrix')
  Res <- adonis2(TempDist~Infection,data=TempDF,permutations = 10000)
  y <- Res$`Pr(>F)`[1]
  y2 <- ifelse(y<=0.001,"p<0.001",y)
  z <- ifelse(y<=0.001,"p<0.001",p.adjust(y,'bonferroni',n=7) )
  print(paste(sep=',','Common Garden vs Non','Jaccard',x,Res$F[1],y2,z))
  rm(x,y,y2,z,Res,TempDist,TempDF,TempObject)
  # Can copy and paste the results into excel for manipulation, comma delimited
}

##### Table S7, Pairwise permanova #####

TempObject <- subset_samples(EvenObject, Infection=='Infected')
TempDist <- as(vegdist(t(otu_table(TempObject)),method='bray'),'matrix')
TempDF <- as(sample_data(TempObject),'data.frame')
pairwise.adonis2(TempDist~Month,TempDF,permutations=10000)
rm(TempObject,TempDist,TempDF)

TempObject <- subset_samples(EvenObject, Infection=='Non-infected')
TempDist <- as(vegdist(t(otu_table(TempObject)),method='bray'),'matrix')
TempDF <- as(sample_data(TempObject),'data.frame')
pairwise.adonis2(TempDist~Month,TempDF,permutations=10000)
rm(TempObject,TempDist,TempDF)

TempObject <- subset_samples(EvenObject, Infection=='Infected')
TempDist <- as(vegdist(t(otu_table(TempObject)),method='jaccard'),'matrix')
TempDF <- as(sample_data(TempObject),'data.frame')
pairwise.adonis2(TempDist~Month,TempDF,permutations=10000)
rm(TempObject,TempDist,TempDF)

TempObject <- subset_samples(EvenObject, Infection=='Non-infected')
TempDist <- as(vegdist(t(otu_table(TempObject)),method='jaccard'),'matrix')
TempDF <- as(sample_data(TempObject),'data.frame')
pairwise.adonis2(TempDist~Month,TempDF,permutations=10000)
rm(TempObject,TempDist,TempDF)

##### Table S9, Beta Disper #####

for (I in c('Non-infected','Infected')){
  if (I == 'Non-infected'){ # on first iteration only
    TableS9 <- data.frame(Infection=vector(),Month=vector(),BC.dist=numeric(),J.dist=numeric() )
  }
  for (x in MonthOrder){
    ps <- prune_samples(sample_data(EvenObject)$Month==x,EvenObject)
    ps2 <- prune_samples(sample_data(ps)$Infection==I,ps)
    ps3 <- prune_taxa((taxa_sums(ps2)>0),ps2)
    rm(ps,ps2) #only need phyloseq object 3
    df <- sample_data(ps3) %>% as('data.frame')
    dist.BC <-t(otu_table(ps3)) %>% vegdist(method='bray')
    dist.J <-t(otu_table(ps3)) %>% vegdist(method='jaccard')
    result.BC<-betadisper(dist.BC, df$Infection, type = "median")
    result.J<-betadisper(dist.J, df$Infection, type = "median")
    TableS9 <- rbind(TableS9,data.frame(Infection=I,Month=x,BC.dist=result.BC$group.distances,J.dist=result.J$group.distances) )
    rm(ps3,df,dist.BC,dist.J,result.BC,result.J)
  }
  if(I=="Infected"){
    TableS9$BC.dist <- round(TableS9$BC.dist,digits=4)
    TableS9$J.dist <- round(TableS9$J.dist,digits=4)
    rownames(TableS9) <- 1:nrow(TableS9)
  }
  rm(I,x)
}

TableS9
