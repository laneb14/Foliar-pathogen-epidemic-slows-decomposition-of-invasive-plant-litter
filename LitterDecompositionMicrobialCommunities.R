
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

##### Calculate alpha diversity #####
# Alpha diversity metrics are calculated using data rarified to an equal sampling depth
RarePhylo <- rarefy_even_depth(PhyloObject)
Indicies <- estimate_richness(RarePhylo,measures=c('Observed','Chao1','Shannon')) # Calculate Chao1 and Shannon's Diversity
rownames(Indicies) <- gsub("\\.","-",rownames(Indicies)) # I'm not sure why it replaces dashes with periods, but this fixes it
Indicies$SEven <- Indicies$Shannon/log(Indicies$Observed) # Calculate Shannon's Evenness
Indicies$Infection <- vlookup(rownames(Indicies),InputMeta,'Infection','SampleName')
Indicies$Month <- vlookup(rownames(Indicies),InputMeta,'Month','SampleName')
Indicies$Replicate <- vlookup(rownames(Indicies),InputMeta,'Replicate','SampleName')

##### Table S6 #####

head(Indicies)

for (x in MonthOrder){
  MonthDF <- Indicies[Indicies$Month==x,]
  # break
  if (x==MonthOrder[1]){
    #Initiate dataframe
    TableS6 <- data.frame(Metric=vector(),Month=vector(),CG=numeric(),NI=numeric(),I=numeric(),CG.NI.p=numeric(),CG.NI.q=numeric(),NI.I.p=numeric(),NI.I.q=numeric(),CG.NI.f=numeric(),NI.I.f=numeric() )
  }
  ## Chao1 ##
  CG.NI <- anova(lm(Chao1~Infection,MonthDF[MonthDF$Infection!='Infected',]))
  NI.I <- anova(lm(Chao1~Infection,MonthDF[MonthDF$Infection!='Common Garden',]))
  CG.NI.q <- p.adjust(CG.NI$`Pr(>F)`[1],'bonferroni',n=7)
  NI.I.q <- p.adjust(NI.I$`Pr(>F)`[1],'bonferroni',n=7)
  TableS6 <- rbind(TableS6,
  data.frame(Metric='Chao1',Month=x,CG=round(mean(MonthDF$Chao1[MonthDF$Infection=='Common Garden']),digits=2),
             NI=round(mean(MonthDF$Chao1[MonthDF$Infection=='Non-infected']),digits=2),I=round(mean(MonthDF$Chao1[MonthDF$Infection=='Infected']),digits=2),
             CG.NI.p=round(CG.NI$`Pr(>F)`[1],3),CG.NI.q=round(CG.NI.q,3),NI.I.p=round(NI.I$`Pr(>F)`[1],3), NI.I.q=round(NI.I.q,3),
             CG.NI.f=CG.NI$`F value`[1]%>%round(3),NI.I.f=NI.I$`F value`[1]%>%round(3) )
  )
  ## Shannon Diversity ##
  CG.NI <- anova(lm(Shannon~Infection,MonthDF[MonthDF$Infection!='Infected',]))
  NI.I <- anova(lm(Shannon~Infection,MonthDF[MonthDF$Infection!='Common Garden',]))
  CG.NI.q <- p.adjust(CG.NI$`Pr(>F)`[1],'bonferroni',n=7)
  NI.I.q <- p.adjust(NI.I$`Pr(>F)`[1],'bonferroni',n=7)
  TableS6 <- rbind(TableS6,
  data.frame(Metric='Shannon Diversity',Month=x,CG=round(mean(MonthDF$Shannon[MonthDF$Infection=='Common Garden']),digits=2),
             NI=round(mean(MonthDF$Shannon[MonthDF$Infection=='Non-infected']),digits=2),I=round(mean(MonthDF$Shannon[MonthDF$Infection=='Infected']),digits=2),
             CG.NI.p=round(CG.NI$`Pr(>F)`[1],3),CG.NI.q=round(CG.NI.q,3),NI.I.p=round(NI.I$`Pr(>F)`[1],3), NI.I.q=round(NI.I.q,3),
             CG.NI.f=CG.NI$`F value`[1]%>%round(3),NI.I.f=NI.I$`F value`[1]%>%round(3) )
  )
  ## Shannon Evenness ##
  CG.NI <- anova(lm(SEven~Infection,MonthDF[MonthDF$Infection!='Infected',]))
  NI.I <- anova(lm(SEven~Infection,MonthDF[MonthDF$Infection!='Common Garden',]))
  CG.NI.q <- p.adjust(CG.NI$`Pr(>F)`[1],'bonferroni',n=7)
  NI.I.q <- p.adjust(NI.I$`Pr(>F)`[1],'bonferroni',n=7)
  TableS6 <- rbind(TableS6,
  data.frame(Metric='Shannon Evenness',Month=x,CG=round(mean(MonthDF$SEven[MonthDF$Infection=='Common Garden']),digits=2),
             NI=round(mean(MonthDF$SEven[MonthDF$Infection=='Non-infected']),digits=2),I=round(mean(MonthDF$SEven[MonthDF$Infection=='Infected']),digits=2),
             CG.NI.p=round(CG.NI$`Pr(>F)`[1],3),CG.NI.q=round(CG.NI.q,3),NI.I.p=round(NI.I$`Pr(>F)`[1],3), NI.I.q=round(NI.I.q,3),
             CG.NI.f=CG.NI$`F value`[1]%>%round(3),NI.I.f=NI.I$`F value`[1]%>%round(3) )
  )
  
  rm(x,CG.NI.q,NI.I.q,CG.NI,NI.I,MonthDF)
}

TableS6[TableS6$Metric=='Chao1',]             %>% select(Metric,Month,CG.NI.f,CG.NI.p)
TableS6[TableS6$Metric=='Shannon Diversity',] %>% select(Metric,Month,CG.NI.f,CG.NI.p)
TableS6[TableS6$Metric=='Shannon Evenness',]  %>% select(Metric,Month,CG.NI.f,CG.NI.p)

TableS6[TableS6$Metric=='Chao1',]             %>% select(Metric,Month,NI.I.f,NI.I.p,NI.I.q)
TableS6[TableS6$Metric=='Shannon Diversity',] %>% select(Metric,Month,NI.I.f,NI.I.p,NI.I.q)
TableS6[TableS6$Metric=='Shannon Evenness',]  %>% select(Metric,Month,NI.I.f,NI.I.p,NI.I.q)

##### Table S7 #####
# Need to do pairwise adonis for both Bray+Jaccard by Healthy+Infected

for (I in c('Non-infected','Infected')){
  if (I == 'Non-infected'){ # on first iteration only
    TableS7 <- data.frame(Infection=vector(),Distance=vector(),Months=vector(),qval=numeric())
  }
  for ( D in c('bray','jaccard')){
    # break
    ps <- prune_samples(sample_data(EvenObject)$Infection==I,EvenObject) # subset samples by infection, iteration of phyloseq used here has a problem using variables in subset_samples
    ps <- prune_taxa(taxa_sums(ps)>0,ps)
    dist <- t(otu_table(ps)) %>% vegdist(method=D)
    # month <- sample_data(ps)$Month
    df <- sample_data(ps) %>% as('data.frame')
    result <- pairwise.adonis(x=dist,factors=df$Month,sim.method=D,perm=10000)
    result2 <- result %>% select(Months=pairs,qval=p.adjusted)
    result2$Infection <- I
    result2$Distance <- D
    TableS7<-rbind(TableS7,result2)
    rm(ps,dist,df,result,result2)
  }
  rm(D,I)
}

head(TableS7)

##### Table S8 #####

  # if (I == 'Non-infected'){ # on first iteration only
    # TableS7 <- data.frame(Infection=vector(),Distance=vector(),Months=vector(),qval=numeric())
  # }
for ( D in c('bray','jaccard')){
  if (D == "bray"){ # on first iteration only
    TableS8 <- data.frame(Distance=vector(),Month=vector(),CG.NI.p=numeric(),CG.NI.q=numeric(),NI.I.p=numeric(),NI.I.q=numeric() )
  }
  for ( x in MonthOrder){
    ## First subset the selected month
    ps <- prune_samples(sample_data(EvenObject)$Month==x,EvenObject) # subset samples by Month, iteration of phyloseq used here has a problem using variables in subset_samples
    ## Common garden vs non-infected
    NoI <- subset_samples(ps,Infection!='Infected')
    NoI <- prune_taxa(taxa_sums(NoI)>0,NoI)
    dist.NoI <-t(otu_table(NoI)) %>% vegdist(method=D)
    df.NoI <- sample_data(NoI) %>% as("data.frame")
    NoI.result <- adonis2(dist.NoI~Infection,data=df.NoI,permutations=10000)
    NoI.q <- NoI.result$`Pr(>F)`[1] %>% p.adjust(method='bonferroni',n=7)
    
    ## Non-infected vs infected
    NoCG <- subset_samples(ps,Infection!='Common Garden')
    NoCG <- prune_taxa(taxa_sums(NoCG)>0,NoCG)
    dist.NoCG <-t(otu_table(NoCG)) %>% vegdist(method=D)
    df.NoCG <- sample_data(NoCG) %>% as("data.frame")
    NoCG.result <- adonis2(dist.NoCG~Infection,data=df.NoCG,permutations=10000)
    NoCG.q <- NoCG.result$`Pr(>F)`[1] %>% p.adjust(method='bonferroni',n=7)
    
    TableS8<- rbind(TableS8,data.frame(Distance=D,Month=x, CG.NI.p=NoI.result$`Pr(>F)`[1],CG.NI.q=NoI.q,NI.I.p=NoCG.result$`Pr(>F)`[1],NI.I.q=NoCG.q ) )
    rm(ps,NoI,NoCG,NoCG.result,NoI.result,NoI.q,NoCG.q,dist.NoCG,dist.NoI,df.NoCG,df.NoI,D,x)
  }
}

##### Table S9 #####

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
TableS9[TableS9$Infection=='Non-infected',]
TableS9[TableS9$Infection=='Infected',]

##### Figure 5 Main ordination #####

EvenObject.NoCG <- subset_samples(EvenObject, Infection!='Common Garden')
data <- sample_data(EvenObject.NoCG) %>% as("data.frame")
BrayDist <- vegdist(t(otu_table(EvenObject.NoCG)), method='bray')
res <- pcoa(BrayDist)
vectors <- res$vectors %>% as.data.frame()
vectors$Month <- vlookup(row.names(vectors),data,'Month','SampleName')
vectors$Infection <- vlookup(row.names(vectors),data,'Infection','SampleName')
# vectors$LitterSource <- vlookup(row.names(vectors),data,'Litter.Source','SampleName')
vectors$Replicate <- vlookup(row.names(vectors),data,'Replicate','SampleName')
OrdNoCG <- vectors

res$values$Relative_eig[1]*100 # Copy these values into the axis labels below
res$values$Relative_eig[2]*100 # Copy these values into the axis labels below

ggplot(OrdNoCG,aes(x=Axis.1*(-1),y=Axis.2,color=factor(Month,levels=MonthOrder), shape=Infection))+
  scale_color_discrete(name='Month')+ theme_bw()+  geom_point()+  stat_ellipse(aes(group=Month))+
  theme(panel.grid = element_blank())+
  scale_shape_manual(values = c('Non-infected'='circle',"Infected"='triangle'),labels=Lvalues,name="Litter Source")+
  ylab("Axis 1  (15.32% Variation)")+xlab("Axis 2  (6.01% Variation)")

##### Figure 6 ######

# These values originated from above
Facet67SigDF <- data.frame(Month=c('Start','February','April','June','August','October','December'),
                           Axis.6=rep(0.23,n=7),Axis.7=rep(-0.3,n=7),
                           # pval=c('p<0.001','p<0.001','p<0.001','p=0.126','p=0.049','p=0.351','p=0.819')  )
                           pval=c('q<0.001','q<0.001','q<0.001','q=0.126','q=0.049','q=0.351','q=0.819')  )

ggplot(OrdNoCG,aes(x=Axis.6,y=Axis.7,color=Infection))+theme_bw()+geom_point()+
  stat_ellipse()+xlab("Axis 6  (2.25% Variation)")+ylab("Axis 7  (2.11% Variation)")+
  theme(panel.grid = element_blank() )+
  facet_wrap(~factor(Month,levels = c('Start','February','April','June','August','October','December')), ncol=2)+
  geom_text(data=Facet67SigDF,inherit.aes=F,
            aes(x=Axis.6,y=Axis.7,label=pval),color='black',size=4 )+
  # theme(legend.position = "none")+
  theme(panel.spacing = unit(1, "lines"))+
  theme(legend.position = c(.65, .15))+
  scale_color_manual(values=c('Non-infected'="#619CFF",'Infected'="#F8766D") ,labels=Lvalues,name='Litter Source')



