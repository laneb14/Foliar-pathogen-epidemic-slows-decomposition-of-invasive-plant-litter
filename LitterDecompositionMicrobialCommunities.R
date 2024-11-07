
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
library(lmerTest)
library(cowplot)

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
# Also remove two samples which manual review showed to be problmatic
PhyloObject <- prune_samples(setdiff(sample_names(PhyloObject),c('Sample-2-67','Sample-3-35')),PhyloObject)
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
pairs(emmeans::emmeans(Shannon.NI.I,~Infection*Month),simple='Infection')
# Sig difference in Shannon's in last timepoint


hist(Indicies$SEven[Indicies$Infection!='Common Garden'])
shapiro.test(Indicies$SEven[Indicies$Infection!='Common Garden'])
hist((Indicies$SEven[Indicies$Infection!='Common Garden'])^3)
shapiro.test((Indicies$SEven[Indicies$Infection!='Common Garden'])^3)
Pielou.NI.I<- lmer((SEven)^3~Infection*Month+(1|Replicate),Indicies[Indicies$Infection!='Common Garden',])
resid_panel(Pielou.NI.I)
anova(Pielou.NI.I)
pairs(emmeans::emmeans(Pielou.NI.I,~Infection*Month),simple='Infection')
# Sig difference in Pielou's in last timepoint

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
res <- pairwise.adonis2(TempDist~Month,TempDF,permutations=10000)
for (x in names(res) ){
  if(x=='parent_call'){next}; temp <- unlist(res[x]);   temp2 <- temp[13]
  temp2 <- round(p.adjust(p=temp2,method='bonferroni',n=21),3)
  if(!grepl(x=names(temp2),pattern='Pr\\(>F\\)')){stop("ERROR IN CODE")}else{ if(temp2!=0.002){ print(paste(x,temp2)) }  }
  # if(!grepl(x=names(temp2),pattern='Pr\\(>F\\)')){stop("ERROR IN CODE")}else{ if(temp2>=0.001){ print(paste(x,temp2)) }  }
}
rm(TempObject,TempDist,TempDF,res)

TempObject <- subset_samples(EvenObject, Infection=='Non-infected')
TempDist <- as(vegdist(t(otu_table(TempObject)),method='bray'),'matrix')
TempDF <- as(sample_data(TempObject),'data.frame')
res <- pairwise.adonis2(TempDist~Month,TempDF,permutations=10000)
for (x in names(res) ){
  if(x=='parent_call'){next}; temp <- unlist(res[x]);   temp2 <- temp[13]
  temp2 <- round(p.adjust(p=temp2,method='bonferroni',n=21),3)
  if(!grepl(x=names(temp2),pattern='Pr\\(>F\\)')){stop("ERROR IN CODE")}else{ if(temp2!=0.002){ print(paste(x,temp2)) }  }
  # if(!grepl(x=names(temp2),pattern='Pr\\(>F\\)')){stop("ERROR IN CODE")}else{ if(temp2>=0.001){ print(paste(x,temp2)) }  }
}
rm(TempObject,TempDist,TempDF,res)

TempObject <- subset_samples(EvenObject, Infection=='Infected')
TempDist <- as(vegdist(t(otu_table(TempObject)),method='jaccard'),'matrix')
TempDF <- as(sample_data(TempObject),'data.frame')
res <- pairwise.adonis2(TempDist~Month,TempDF,permutations=10000)
for (x in names(res) ){
  if(x=='parent_call'){next}; temp <- unlist(res[x]);   temp2 <- temp[13]
  temp2 <- round(p.adjust(p=temp2,method='bonferroni',n=21),3)
  if(!grepl(x=names(temp2),pattern='Pr\\(>F\\)')){stop("ERROR IN CODE")}else{ if(temp2!=0.002){ print(paste(x,temp2)) }  }
  # if(!grepl(x=names(temp2),pattern='Pr\\(>F\\)')){stop("ERROR IN CODE")}else{ if(temp2>=0.001){ print(paste(x,temp2)) }  }
}
rm(TempObject,TempDist,TempDF,res)

TempObject <- subset_samples(EvenObject, Infection=='Non-infected')
TempDist <- as(vegdist(t(otu_table(TempObject)),method='jaccard'),'matrix')
TempDF <- as(sample_data(TempObject),'data.frame')
res <- pairwise.adonis2(TempDist~Month,TempDF,permutations=10000)
for (x in names(res) ){
  if(x=='parent_call'){next}; temp <- unlist(res[x]);   temp2 <- temp[13]
  temp2 <- round(p.adjust(p=temp2,method='bonferroni',n=21),3)
  if(!grepl(x=names(temp2),pattern='Pr\\(>F\\)')){stop("ERROR IN CODE")}else{ if(temp2!=0.002){ print(paste(x,temp2)) }  }
  # if(!grepl(x=names(temp2),pattern='Pr\\(>F\\)')){stop("ERROR IN CODE")}else{ if(temp2>=0.001){ print(paste(x,temp2)) }  }
}
rm(TempObject,TempDist,TempDF,res)

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

##### Figure 3 #####

EvenObject.NoCG <- subset_samples(EvenObject,Infection!='Common Garden')
BrayDist <- as.matrix(vegdist(t(otu_table(EvenObject.NoCG)),method='bray'))

if(all(rownames(BrayDist)==rownames(sample_data(EvenObject.NoCG)))){
  df <- as(sample_data(EvenObject.NoCG),'data.frame')
}else{
  stop("Row names didn't match up, this needs corrected")
}

PCOA <- pcoa(BrayDist)
if(all(rownames(df)==rownames(PCOA$vectors))){
  df2 <- cbind(df,PCOA$vectors)
}else{
  stop("Row names didn't match")  
}

# Axis percent values
PCOA$values$Relative_eig[1]*100
PCOA$values$Relative_eig[2]*100

OrderMonths <- c("Start",'February','April','June','August','October','December')
LevelNames <- c('Non-infected'="Healthy Sites", 'Infected'="Infected Sites",'Common Garden'="Common Garden Sites")
LevelColors <- c('Non-infected'="#619CFF", 'Infected'="#F8766D",'Common Garden'="#00BA38")
head(df2)


PlotPCOA <-
ggplot(df2,aes(x=(Axis.1*(-1)),y=Axis.2,color=factor(Month,levels=OrderMonths),shape=Infection ) )+
  geom_point()+ stat_ellipse(inherit.aes=F, aes(x=(Axis.1*(-1)),y=Axis.2,color=factor(Month,levels=OrderMonths)))+
  scale_color_discrete(name='Month')+
  scale_shape_discrete(name='Litter Source',labels=LevelNames)+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Axis 1 (15.31% Variation)")+
  ylab("Axis 2 (6.09% Variation)")

PlotPCOA

# ggsave(paste(sep='',"PCOA",format(Sys.time(),"%e%b%Y"),'.jpg'),PlotPCOA,dpi=1000,width=8,height=6)

##### Figure S3 #####

PCOA$values$Relative_eig[6]*100
PCOA$values$Relative_eig[7]*100

Fig.S3.Text.DF <- data.frame(Text=c('q<0.001','q<0.001','q<0.001','q=0.175','q=0.027','q=0.449','q=0.832'),Month=c("Start",'February','April','June','August','October','December'))

df2$Infection

FigureS3 <-
ggplot(df2,aes(x=(Axis.6),y=Axis.7,color=Infection ) )+
  facet_wrap(~factor(Month,levels =OrderMonths ),ncol=2)+
  geom_point()+  stat_ellipse()+
  scale_color_manual(values=LevelColors,name='Litter Source',labels=LevelNames)+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(legend.position = c(0.75,0.1))+
  geom_text(inherit.aes=F,data=Fig.S3.Text.DF,aes(y=(-0.275),x=0.265,label=Text))+
  xlab("Axis 6 (2.27% Variation)")+
  ylab("Axis 7 (2.13% Variation)")

# ggsave(paste(sep='',"Facet PCOA ",format(Sys.time(),"%e%b%Y"),'.jpg'),FigureS3,dpi=1000,width=5,height=6)

##### Combine plots #####

FigureS3v3 <-
  ggplot(df2,aes(x=(Axis.6),y=Axis.7,color=Infection ) )+
  facet_wrap(~factor(Month,levels =OrderMonths ),nrow=2)+
  geom_point()+  stat_ellipse()+
  scale_color_manual(values=LevelColors,name='Litter Source',labels=LevelNames)+
  theme_bw()+theme(panel.grid = element_blank())+
  # theme(legend.position = c(0.75,0.1))+
  theme(legend.position = c(0.875,0.225))+
  geom_text(inherit.aes=F,data=Fig.S3.Text.DF,aes(y=(-0.275),x=0.2,label=Text))+
  xlab("Axis 6 (2.27% Variation)")+
  ylab("Axis 7 (2.13% Variation)")

CombinedPCOA.Vert <-
plot_grid(PlotPCOA+theme(legend.position='bottom')+
            guides(color = guide_legend(nrow = 1),legend.position=c(0,0) ) ,
          FigureS3v3,
          ncol=1,labels = c("A","B"),label_size = 18)

# ggsave(paste(sep='',"CombinedPCOA-Vertical",format(Sys.time(),"%e%b%Y"),'.jpg'),CombinedPCOA.Vert,dpi=1000,width=10,height=12)

