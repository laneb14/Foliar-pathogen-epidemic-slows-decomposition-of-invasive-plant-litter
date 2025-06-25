
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
library(ANCOMBC)
library(microbiome)

# Note: Throughout I use Shannon's Evenness, or "SEven" in place of Pielou's. These metrics are the same, but named differently.

# Using rstudioapi::getSourceEditorContext()$path will set the working directory to where this R file is saved
# Taking this approach allows the user to be able to access code saved on a cloud directory (i.e. dropbox) from multiple computers without having to change setwd() each time (e.g. C:\ vs D:\)
setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()

DataFolder <- "edi.1434.1"

##### Persistent Variables #####
MonthOrder <- c('Start','February','April','June','August','October','December')
Lvalues <- c("Non-infected"="Non-infected Sites","Infected"='Infected Sites')
MonthConversion <- c('0'='Start','2'='February','4'='April','6'='June','8'='August','10'='October',"12"='December')

##### Load data #####
InputMeta <- read.table(file = paste(sep='/',DataFolder, 'MetadataMicrostegiumDecomposition2020Sequencing.tsv'), sep = '\t', header = TRUE) # Load metadata
InputOTUTable.LongFormat <- read.table(file = paste(sep='/',DataFolder, 'OTU_Table_InLongFormat.tsv'), sep = '\t', header = TRUE) # Load otu table
length(unique(InputOTUTable.LongFormat$ASV))
# OTU Table had to be uploaded to data base in long format, but phyloseq will want it in wide format
InputOTUTable <-spread(InputOTUTable.LongFormat,key='Sample',value='Reads')
rownames(InputOTUTable) <- InputOTUTable$ASV
InputOTUTable$ASV <- NULL
InputOTUTable <- as.matrix(InputOTUTable) # phyloseq wants it as a matrix
head(InputOTUTable)
colnames(InputOTUTable) <- gsub("\\.","-",colnames(InputOTUTable)) # May be needed but not always. Sometimes the dashes get replaced by dots when bringing the data into R. This puts it back if needed
rownames(InputMeta) <- InputMeta$SampleName # Copy sample names into the row names

##### Clean taxa table #####
InputTaxonomy <- read.table(file = paste('taxonomy2.tsv.txt'), sep = '\t', header = TRUE) # Load metadata

head(InputTaxonomy)
for (x in 1:nrow(InputTaxonomy)){
  if(x==1){ OrganizeTaxMat <- matrix(nrow=nrow(InputOTUTable),ncol=7) ; colnames(OrganizeTaxMat) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species') ; rownames(OrganizeTaxMat) <- rownames(InputOTUTable) }
  if(x%%500==0 | x==1 | x==nrow(InputTaxonomy)){ print(paste("Starting",x,'of',nrow(InputTaxonomy)))  }
  temp <- InputTaxonomy$Taxon[x]
  OrganizeTaxMat[InputTaxonomy$ASV[x],'Kingdom'] <- ifelse(grepl(pattern='k__',x=temp),sub(pattern=';.*','',x=sub(pattern='.*k__','',x=temp)),NA)
  OrganizeTaxMat[InputTaxonomy$ASV[x],'Phylum'] <- ifelse(grepl(pattern='p__',x=temp),sub(pattern=';.*','',x=sub(pattern='.*p__','',x=temp)),NA)
  OrganizeTaxMat[InputTaxonomy$ASV[x],'Class'] <- ifelse(grepl(pattern='c__',x=temp),sub(pattern=';.*','',x=sub(pattern='.*c__','',x=temp)),NA)
  OrganizeTaxMat[InputTaxonomy$ASV[x],'Order'] <- ifelse(grepl(pattern='o__',x=temp),sub(pattern=';.*','',x=sub(pattern='.*o__','',x=temp)),NA)
  OrganizeTaxMat[InputTaxonomy$ASV[x],'Family'] <- ifelse(grepl(pattern='f__',x=temp),sub(pattern=';.*','',x=sub(pattern='.*f__','',x=temp)),NA)
  OrganizeTaxMat[InputTaxonomy$ASV[x],'Genus'] <- ifelse(grepl(pattern='g__',x=temp),sub(pattern=';.*','',x=sub(pattern='.*g__','',x=temp)),NA)
  OrganizeTaxMat[InputTaxonomy$ASV[x],'Species'] <- ifelse(grepl(pattern='s__',x=temp),sub(pattern=';.*','',x=sub(pattern='.*s__','',x=temp)),NA)
  # stop('Manual stop called')
  rm(x,temp)
}

# The genus Xenasmatella is an "incertae sedis", but that has been resolved in the literature.
# Update the taxonomy chart for Xenasmatella because it has significant downstream impacts
for (x in 1:nrow(OrganizeTaxMat)){
  if(is.na(OrganizeTaxMat[x,"Genus"])){rm(x);next}
  if(OrganizeTaxMat[x,"Genus"]=='Xenasmatella'){ OrganizeTaxMat[x,"Family"] <- 'Xenasmataceae' ; OrganizeTaxMat[x,"Order"] <- 'Polyporales' ; rm(x) }
}

OriginalPhyloObject <- phyloseq(otu_table(InputOTUTable,taxa_are_rows = T),sample_data(InputMeta),tax_table(OrganizeTaxMat)) # Create the phyloseq object

##### Remove non-fungal reads #####

TempDF <- data.frame(ASV=as.vector(rownames(tax_table(OriginalPhyloObject)[,'Kingdom'])), K=as.vector(tax_table(OriginalPhyloObject)[,'Kingdom']) )
OriginalPhyloObject <- prune_taxa(TempDF$ASV[which(TempDF$K=='Fungi')],x=OriginalPhyloObject)
rm(TempDF)

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
Indicies$Litter.Source <- vlookup(rownames(Indicies),InputMeta,'Litter.Source','SampleName')

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

hist(Indicies$Chao1[Indicies$Infection!='Common Garden']) ; shapiro.test(Indicies$Chao1[Indicies$Infection!='Common Garden'])
hist(sqrt(Indicies$Chao1[Indicies$Infection!='Common Garden'])) ; shapiro.test(sqrt(Indicies$Chao1[Indicies$Infection!='Common Garden']))
Chao1.NI.I <- lmer(sqrt(Chao1)~Infection*Month+(1|Replicate)+(1|Litter.Source),Indicies[Indicies$Infection!='Common Garden',])
resid_panel(Chao1.NI.I)
anova(Chao1.NI.I)


hist(Indicies$Shannon[Indicies$Infection!='Common Garden'])
shapiro.test(Indicies$Shannon[Indicies$Infection!='Common Garden'])
hist((Indicies$Shannon[Indicies$Infection!='Common Garden'])^(3))
shapiro.test((Indicies$Shannon[Indicies$Infection!='Common Garden'])^(3))
Shannon.NI.I <- lmer((Shannon)^3~Infection*Month+(1|Replicate)+(1|Litter.Source),Indicies[Indicies$Infection!='Common Garden',])
resid_panel(Shannon.NI.I)
anova(Shannon.NI.I)
pairs(emmeans::emmeans(Shannon.NI.I,~Infection*Month),simple='Infection')
# Sig difference in Shannon's in last timepoint


hist(Indicies$SEven[Indicies$Infection!='Common Garden'])
shapiro.test(Indicies$SEven[Indicies$Infection!='Common Garden'])
hist((Indicies$SEven[Indicies$Infection!='Common Garden'])^3)
shapiro.test((Indicies$SEven[Indicies$Infection!='Common Garden'])^3)
Pielou.NI.I<- lmer((SEven)^3~Infection*Month+(1|Replicate)+(1|Litter.Source),Indicies[Indicies$Infection!='Common Garden',])
resid_panel(Pielou.NI.I)
anova(Pielou.NI.I)
pairs(emmeans::emmeans(Pielou.NI.I,~Infection*Month),simple='Infection')
# Sig difference in Pielou's in last timepoint

##### Same as above but compare Common Garden with Non-Infected #####

CG.NI.DF <- Indicies[Indicies$Infection!='Infected',]

hist(CG.NI.DF$Chao1) ; shapiro.test(CG.NI.DF$Chao1)
Chao1.CG <- lmer(Chao1~Infection*Month+(1|Replicate)+(1|Litter.Source), CG.NI.DF)
resid_panel(Chao1.CG)
anova(Chao1.CG)


hist(CG.NI.DF$Shannon) ; shapiro.test(CG.NI.DF$Shannon)
hist(CG.NI.DF$Shannon^2) ; shapiro.test(CG.NI.DF$Shannon^2)
CG.NI.DF$Shannon.sq <- CG.NI.DF$Shannon^2
Shannon.CG <- lmer(Shannon.sq~Infection*Month+(1|Replicate)+(1|Litter.Source), CG.NI.DF)
resid_panel(Shannon.CG)
anova(Shannon.CG)

hist(CG.NI.DF$SEven) ; shapiro.test(CG.NI.DF$SEven)
hist(CG.NI.DF$SEven^2) ; shapiro.test(CG.NI.DF$SEven^2)
CG.NI.DF$SEven.sq <- CG.NI.DF$SEven^2
Pielous.CG <- lmer(SEven.sq~Infection*Month+(1|Replicate)+(1|Litter.Source), CG.NI.DF)
resid_panel(Pielous.CG)
anova(Pielous.CG)



##### Table S6, Pairwise permanova #####

for (x in MonthOrder){
  head(sample_data(EvenObject))
  if(x=='Start'){print(paste(sep=',','Comparison','Metric','Month','DF','Residual DF','R2','F value','p value','q value'))}
  TempObject <- subset_samples(EvenObject, Infection!='Common Garden' & Month==x)
  TempDist <- as(vegdist(t(otu_table(TempObject)),method='bray'),'matrix')
  TempDF <- as(sample_data(TempObject),'data.frame')
  Res <- adonis2(TempDist~Infection,data=TempDF,permutations = 10000)
  y <- Res$`Pr(>F)`[1]
  y2 <- ifelse(y<=0.001,"<0.001",y)
  z <- ifelse(y<=0.001,"<0.001",p.adjust(y,'bonferroni',n=7) )
  extra <- paste(Res$Df[1],Res$Df[2],Res$R2[1],sep=',')
  print(paste(sep=',','Non vs Infected','Bray-Curtis',x,extra,Res$F[1],y2,z))
  rm(TempDist,y,y2,z,Res,extra)
  TempDist <- as(vegdist(t(otu_table(TempObject)),method='jaccard'),'matrix')
  Res <- adonis2(TempDist~Infection,data=TempDF,permutations = 10000)
  y <- Res$`Pr(>F)`[1]
  y2 <- ifelse(y<=0.001,"<0.001",y)
  z <- ifelse(y<=0.001,"<0.001",p.adjust(y,'bonferroni',n=7) )
  extra <- paste(Res$Df[1],Res$Df[2],Res$R2[1],sep=',')
  print(paste(sep=',','Non vs Infected','Jaccard',x,extra,Res$F[1],y2,z))
  rm(x,y,y2,z,Res,TempDist,TempDF,TempObject,extra)
  # Can copy and paste the results into excel for manipulation, comma delimited
}

for (x in MonthOrder){
  head(sample_data(EvenObject))
  TempObject <- subset_samples(EvenObject, Infection!='Infected' & Month==x)
  TempDist <- as(vegdist(t(otu_table(TempObject)),method='bray'),'matrix')
  TempDF <- as(sample_data(TempObject),'data.frame')
  Res <- adonis2(TempDist~Infection,data=TempDF,permutations = 10000)
  y <- Res$`Pr(>F)`[1]
  y2 <- ifelse(y<=0.001,"<0.001",y)
  z <- ifelse(y<=0.001,"<0.001",p.adjust(y,'bonferroni',n=7) )
  extra <- paste(Res$Df[1],Res$Df[2],Res$R2[1],sep=',')
  print(paste(sep=',','Common Garden vs Non','Bray-Curtis',x,extra,Res$F[1],y2,z))
  rm(TempDist,y,y2,z,Res,extra)
  TempDist <- as(vegdist(t(otu_table(TempObject)),method='jaccard'),'matrix')
  Res <- adonis2(TempDist~Infection,data=TempDF,permutations = 10000)
  y <- Res$`Pr(>F)`[1]
  y2 <- ifelse(y<=0.001,"<0.001",y)
  z <- ifelse(y<=0.001,"<0.001",p.adjust(y,'bonferroni',n=7) )
  extra <- paste(Res$Df[1],Res$Df[2],Res$R2[1],sep=',')
  print(paste(sep=',','Common Garden vs Non','Jaccard',x,extra,Res$F[1],y2,z))
  rm(x,y,y2,z,Res,TempDist,TempDF,TempObject,extra)
  # Can copy and paste the results into excel for manipulation, comma delimited
}

##### Table S7, Pairwise permanova #####

for (INFECT in unique(sample_data(EvenObject)$Infection )){
  if(INFECT=='Common Garden'){rm(INFECT);next}
  for(METRIC in c('bray','jaccard')){
    TempObject <- subset_samples(EvenObject, Infection==INFECT)
    TempDist <- as(vegdist(t(otu_table(TempObject)),method=METRIC),'matrix')
    TempDF <- as(sample_data(TempObject),'data.frame')
    res <- pairwise.adonis2(TempDist~Month,TempDF,permutations=10000)
    for (x in 1:length(res) ){
      if(names(res)[x]=='parent_call'){rm(x);next}
      tempname <- names(res)[x] ; m1 <- sub(pattern='_vs.*',replacement = '',x=tempname) ; m2 <- sub(pattern='.*vs_',replacement = '',x=tempname)
      m1num <- as.numeric(TempDF$MonthNumber[TempDF$Month==m1][1]) ; m2num <- as.numeric(TempDF$MonthNumber[TempDF$Month==m2][1])
      # m1num and m2num are for sorting purposes and will be deleted in post processing
      (temp <- res[x][[1]]) ; p.val <- ifelse(temp$`Pr(>F)`[1] < 0.001,'<0.001',temp$`Pr(>F)`[1])
      q.val <- p.adjust(temp$`Pr(>F)`[1],n=21,method='bonferroni') ; q.val2 <- ifelse(temp$`Pr(>F)`[1] < 1e-4,"<0.002",q.val)
      print(paste(INFECT,METRIC,m1,m2,m1num,m2num,temp$Df[1],temp$Df[2],temp$R2[1],temp$F[1],p.val,q.val2,sep=',' ))
      rm(x,tempname,m1,m2,m1num,m2num,temp,p.val,q.val,q.val2)
    }
    rm(METRIC)
  }
  rm(INFECT)
}

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
# Note that this is no longer Figure S3, this is now Figure 3 part B

PCOA$values$Relative_eig[6]*100
PCOA$values$Relative_eig[7]*100

# Fig.S3.Text.DF <- data.frame(Text=c('q<0.001','q<0.001','q<0.001','q=0.189','q=0.032','q=0.420','q=0.883'),Month=c("Start",'February','April','June','August','October','December'))


FigureS3 <-
ggplot(df2,aes(x=(Axis.6),y=Axis.7,color=Infection ) )+
  facet_wrap(~factor(Month,levels =OrderMonths ),ncol=2)+
  geom_point()+  stat_ellipse()+
  scale_color_manual(values=LevelColors,name='Litter Source',labels=LevelNames)+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(legend.position = c(0.75,0.1))+
  # geom_text(inherit.aes=F,data=Fig.S3.Text.DF,aes(y=(-0.275),x=0.265,label=Text))+
  xlab("Axis 6 (2.27% Variation)")+
  ylab("Axis 7 (2.13% Variation)")

# ggsave(paste(sep='',"Facet PCOA ",format(Sys.time(),"%e%b%Y"),'.jpg'),FigureS3,dpi=1000,width=5,height=6)

##### Combine plots #####

(FigureS3v3 <-
  ggplot(df2,aes(x=(Axis.6),y=Axis.7,color=Infection ) )+
  facet_wrap(~factor(Month,levels =OrderMonths ),ncol=3)+
  geom_point()+  stat_ellipse()+
  scale_color_manual(values=LevelColors,name='Litter Source',labels=LevelNames)+
  theme_bw()+theme(panel.grid = element_blank())+
  # theme(legend.position = c(0.875,0.225))+
  # geom_text(inherit.aes=F,data=Fig.S3.Text.DF,aes(y=(-0.275),x=0.2,label=Text))+
  xlab("Axis 6 (2.27% Variation)")+
  ylab("Axis 7 (2.13% Variation)")
)

# PlotPCOA
# PlotPCOA+coord_fixed()
?plot_grid

# p1 <- PlotPCOA+coord_fixed()+#theme(legend.position='bottom')+
#   guides(color = guide_legend(ncol = 1),legend.position=c(0,0), )
# 
# p2 <- FigureS3v3+coord_fixed()+theme(legend.position = c(0.69,0.2))+
#   guides(color=guide_legend(nrow=1))
# rm(p1,p2)

(CombinedPCOA.Vert <-
plot_grid(PlotPCOA+coord_fixed()+#theme(legend.position='bottom')+
            guides(color = guide_legend(ncol = 1),legend.position=c(0,0), ),
          FigureS3v3+coord_fixed()+theme(legend.position = c(0.46,0.225))+guides(color=guide_legend(ncol=1)),
          # FigureS3v3+coord_fixed()+theme(legend.position = c(0.525,0.15))+
          # FigureS3v3+coord_fixed()+theme(legend.position = c(0.575,0.25))+guides(color=guide_legend(nrow=1)),
          ncol=1,labels = c("A","B"),label_size = 18,align = "v",axis='lr')
)


# ggsave('CombinedPCOA-Vertical.10June2025v1.jpg',CombinedPCOA.Vert,dpi=1000,width=8,height=12)
getwd()

# ggsave(paste(sep='',"CombinedPCOA-Vertical",format(Sys.time(),"%e%b%Y"),'.jpg'),CombinedPCOA.Vert,dpi=1000,width=10,height=12)

##### Changes in the microbiome #####

PhyloObject.NoCG <- subset_samples(PhyloObject,CommonGarden==0) ; PhyloObject.NoCG <- prune_taxa(taxa_sums(PhyloObject.NoCG)>0,PhyloObject.NoCG)
EvenObject.NoCG <- subset_samples(EvenObject,CommonGarden==0) ; EvenObject.NoCG <- prune_taxa(taxa_sums(EvenObject.NoCG)>0,EvenObject.NoCG)

for (x in c(0,2,4,6,8,10,12)){
  if(x==0 & exists('LargeRes')){stop("Don't overwrite")} # this takes a while to run, prevents accidental overwrite
  print(paste("Starting month",x))
  TempObject <- subset_samples(PhyloObject.NoCG,MonthNumber==x) ; TempObject <- prune_taxa(taxa_sums(TempObject)>0,TempObject)  
  res <- ancombc2(TempObject,fix_formula='Infection',pseudo_sens =T) %>% suppressMessages() # suppressMessages becuase this script is very vocal
  # res <- ancombc2(TempObject,fix_formula='Infection',pseudo_sens =F) %>% suppressMessages() # use pseudo_sens false for tests, but use true for actual results
  temp <- res$res ; temp$MonthNumber <- x
  temp2 <- (taxa_sums(TempObject)[temp$taxon]/sum(taxa_sums(TempObject)))*100
  if(!all(temp$taxon==names(temp2))){stop('Names dont match')}
  temp$PercentReads <- temp2
  if(x==0){LargeRes <- temp}else{LargeRes <-rbind(LargeRes,temp) }
  rm(x,TempObject,res,temp,temp2)
}

if(length(unique(LargeRes$MonthNumber))!=7){stop("Missing a month")}

# next seek out the best ID
for (x in 1:nrow(LargeRes)){
  if(x==1 | x==nrow(LargeRes)|x%%250==0){  print(paste('Starting',x,'of',nrow(LargeRes))) }
  temp <- as.vector(tax_table(PhyloObject.NoCG)[LargeRes$taxon[x],'Genus'])
  temp <- ifelse((is.na(temp) | temp=='unidentified'),as.vector(tax_table(PhyloObject.NoCG)[LargeRes$taxon[x],'Family']),temp)
  temp <- ifelse((is.na(temp) | temp=='unidentified'),as.vector(tax_table(PhyloObject.NoCG)[LargeRes$taxon[x],'Order']),temp)
  temp <- ifelse((is.na(temp) | temp=='unidentified'),as.vector(tax_table(PhyloObject.NoCG)[LargeRes$taxon[x],'Class']),temp)
  temp <- ifelse((is.na(temp) | temp=='unidentified'),as.vector(tax_table(PhyloObject.NoCG)[LargeRes$taxon[x],'Phylum']),temp)
  temp <- ifelse((is.na(temp) | temp=='unidentified'),as.vector(tax_table(PhyloObject.NoCG)[LargeRes$taxon[x],'Kingdom']),temp)
  LargeRes$ID[x] <- temp ; rm(temp,x)
}

# First I want to identify anything that is lfc>2 and a qvalue that passes per ANCOMBC2
LargeRes$ColorPoint <- LargeRes$`diff_InfectionNon-infected` & (abs(LargeRes$`lfc_InfectionNon-infected`)>2)
LargeRes$ID[LargeRes$ColorPoint==T] %>% unique ; cat('\n',LargeRes$ID[LargeRes$ColorPoint==T] %>% unique %>% length, ' total taxa\n',LargeRes$ColorPoint %>% sum,' total ASVs',sep='')
LargeRes[LargeRes$ColorPoint==T,] %>% dplyr::select(MonthNumber,taxon,`lfc_InfectionNon-infected`,PercentReads,ID)
LargeRes$ColorPoint %>% sum
LargeRes$MonthNumber[LargeRes$ColorPoint==T] %>%table

# Function List
LongFunctionList <- data.frame(ID=c("Periconia","Pleomassariaceae","Didymellaceae","Magnaporthaceae","Mycosphaerellaceae","Pleosporales","Tetraplosphaeria","Apodus","Articulospora","Ceramothyrium","Chloridium","Clitocybe","Deconica","Flagellospora","Helotiaceae","Helotiales","Hyaloscyphaceae","Melanommataceae","Phallus","Phialocephala","Pseudodactylaria","Spirosphaera","Alfaria","Ascomycota","Aureobasidium","Bulleribasidium","Ceratocystis","Cercospora","Chaetothyriales","Cistella","Cylindrium","Didymocyrtis","Keissleriella","Leotiomycetes","Metacordyceps","Neoascochyta","Neodevriesia","Phaeosphaeria","Phialophora","Plectosphaerella","Sordariomycetes","Sporobolomyces","Tetracladium","Tremellales","Trichomeriaceae","Tripospermum","Tubaria","Venturia"),
                               Fun=c("Both","Saprophytic","Pathogenic","Pathogenic","Pathogenic","Unknown","Pathogenic","Saprophytic","Saprophytic","Both","Saprophytic","Saprophytic","Saprophytic","Saprophytic","Saprophytic","Unknown","Saprophytic","Saprophytic","Saprophytic","Saprophytic","Saprophytic","Saprophytic","Saprophytic","Unknown","Both","Saprophytic","Both","Pathogenic","Unknown","Saprophytic","Saprophytic","Other","Saprophytic","Unknown","Other","Pathogenic","Saprophytic","Pathogenic","Saprophytic","Saprophytic","Unknown","Saprophytic","Saprophytic","Unknown","Pathogenic","Saprophytic","Saprophytic","Pathogenic") )
LongFunctionList$Fun <- ifelse(LongFunctionList$Fun=='Pathogenic',"Plant-associated",LongFunctionList$Fun)
# Pathogenic was my shorthand, but it should really say Plant-associated

head(LongFunctionList)

# Add function to LargeRes
LargeRes$Function <- vlookup(LargeRes$ID,LongFunctionList,'Fun','ID')
# There will still be some missing because I only looked up for ID's which are colored
head(LargeRes)

# Next label the ones to make a bidirectional of
LargeRes$BiDirect <- (LargeRes$PercentReads>0.05)&LargeRes$ColorPoint
LargeRes$BiDirect %>% sum

# Estimate where the holm line is
estimate.holm <- mean(max(LargeRes$`p_InfectionNon-infected`[LargeRes$`diff_InfectionNon-infected`]),min(LargeRes$`p_InfectionNon-infected`[!LargeRes$`diff_InfectionNon-infected`]))

##### Supplementary Table #####
head(LargeRes,1)
ColToPrint <- c('MonthNumber','ID','Function','PercentReads','lfc_(Intercept)','lfc_InfectionNon-infected','se_(Intercept)','se_InfectionNon-infected','W_(Intercept)','W_InfectionNon-infected','p_(Intercept)','p_InfectionNon-infected','q_(Intercept)','q_InfectionNon-infected')
length(ColToPrint)
for (x in 1:nrow(LargeRes)){
  if(x==1){ print(paste(ColToPrint[1],ColToPrint[2],ColToPrint[3],ColToPrint[4],ColToPrint[5],ColToPrint[6],ColToPrint[7],ColToPrint[8],ColToPrint[9],ColToPrint[10],ColToPrint[11],ColToPrint[12],ColToPrint[13],ColToPrint[14],sep=',')) }
  if(!LargeRes$`diff_InfectionNon-infected`[x]){rm(x);next}
  if(abs(LargeRes$`lfc_InfectionNon-infected`[x])<2){rm(x);next}
  temp <- LargeRes[x,] %>% dplyr::select(all_of(ColToPrint))
  print(paste(temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],temp[7],temp[8],temp[9],temp[10],temp[11],temp[12],temp[13],temp[14],sep=','))
  rm(x,temp)
  # stop()
}


##### Plot volcano plots #####

FunctionColors <- c("#00ff00","#0000ff","#008b8b","#ffa500","#ffd700")
names(FunctionColors) <- c("Plant-associated",'Both','Unknown','Saprophytic','Other')


VolcanoPlot <-
ggplot(LargeRes,aes(x=(`lfc_InfectionNon-infected`),y=(-log10(`p_InfectionNon-infected`))))+theme_bw()+
  geom_vline(xintercept = 2    ,linetype='dashed',color='grey40')+
  geom_vline(xintercept = (-2) ,linetype='dashed',color='grey40')+
  geom_hline(yintercept=(-log10(estimate.holm)),linetype='dashed',color='grey40')+
  geom_point(data=LargeRes[LargeRes$ColorPoint==F,],color='grey',size=2)+
  geom_point(data=LargeRes[LargeRes$ColorPoint==T & LargeRes$BiDirect==F,],aes(color=Function),size=2)+
  geom_point(data=LargeRes[LargeRes$BiDirect==T,],aes(color=Function),size=2,pch=17)+
  scale_color_manual(values=FunctionColors,breaks=names(FunctionColors))+
  scale_fill_manual(values=FunctionColors,breaks=names(FunctionColors))+
  xlab('Natural log fold change')+
  ylab(expression(-log[10]*"(p-value)"))+
  facet_wrap(~MonthNumber,labeller=labeller(MonthNumber=MonthConversion),ncol=2)+
  scale_x_continuous(limits=c(-5.4,5.4),breaks=c(-5,-2.5,0,2.5,5),labels=c("5",'2.5\nProminent in\ninfected litter','0','2.5\nProminent in\nhealthy litter','5'))+
  theme(panel.grid = element_blank())

##### Plot bidrectional #####

Label.Y.BiDirect <- LargeRes$ID[LargeRes$BiDirect] ; names(Label.Y.BiDirect) <- LargeRes$taxon[LargeRes$BiDirect]

head(LargeRes)

BiDirectPlot <-
ggplot(data=LargeRes[LargeRes$BiDirect,],aes(y=taxon,x=`lfc_InfectionNon-infected`))+theme_bw()+
  geom_bar(stat='identity',aes(fill=Function),color='black')+
  geom_errorbar(aes(xmin=`lfc_InfectionNon-infected`-`se_InfectionNon-infected`,xmax=`lfc_InfectionNon-infected`+`se_InfectionNon-infected`),width=0.4)+
  facet_grid(MonthNumber~.,scales='free_y',space='free_y',labeller=labeller(MonthNumber=MonthConversion))+
  theme(strip.text.y = element_text(angle = 0))+
  scale_fill_manual(values=FunctionColors)+
  scale_y_discrete(labels=Label.Y.BiDirect)+
  theme(panel.grid = element_blank())+
  scale_x_continuous(limits=c(-5.4,5.4),breaks=c(-5,-2.5,0,2.5,5),labels=c("5",'2.5\nProminent in\ninfected litter','0','2.5\nProminent in\nhealthy litter','5'))+
  xlab('Natural log fold change')+ylab('')

VolLegend <- get_legend(VolcanoPlot+guides(color = guide_legend(title="Function of taxa",nrow=3,override.aes = list(shape = 22,size = 5,fill = FunctionColors, colour = "black"))) )


CombinedVolPlot <- 
plot_grid(VolcanoPlot+theme(legend.position = 'none'),
          BiDirectPlot+theme(legend.position = 'none'),
          labels = c("A","B"))


CombinedVolAndBidirPlot <- 
ggdraw()+ 
  draw_plot(CombinedVolPlot, 0, 0, 1, 1) +
  draw_plot(VolLegend, -.1, .1, 1, 0.15)
  # draw_plot(VolLegend, -.1, .1, 1, 0.12)

ggsave(filename='VolcanoAndBidirectional.jpg',plot=CombinedVolAndBidirPlot,dpi=1000,height=6.5,width=10)
getwd()


##### Make a chart of the top 10 orders #####

for (x in c(0,2,4,6,8,10,12)){
  LEVEL <- 'Order'
  if(x==0){ GraphTaxaChart <- data.frame(Level=vector(),Month=vector(),MonthNumber=vector(),Infection=vector(),ASV=vector(),Taxa=vector(),PercentReads=vector()) }
  for (i in unique(sample_data(EvenObject.NoCG)$Infection) ){
    i2 <- ifelse(i=='Non-infected','Healthy',i)
    TempObject <- subset_samples(EvenObject.NoCG,MonthNumber==x & Infection==i) ; TempObject <- prune_taxa(taxa_sums(TempObject)>0,TempObject)  
    GlomObject <- tax_glom(TempObject,taxrank = LEVEL); SortSums <- sort(taxa_sums(GlomObject),decreasing = T); SortSums <- ((SortSums/sum(SortSums))*100)
    TempDF <- data.frame(Level=LEVEL,Month=sample_data(TempObject)$Month[1],MonthNumber=x,Infection=i2,ASV=names(SortSums[1:10]),Taxa=NA,PercentReads=SortSums[1:10])
    for (y in 1:10){ TempDF$Taxa[y] <- as.vector(tax_table(GlomObject)[names(SortSums)[y],LEVEL]) }
    GraphTaxaChart <- rbind(GraphTaxaChart,TempDF)
    rm(i,i2,TempObject,GlomObject,TempDF,SortSums)
  }
  rm(x,LEVEL)
}

GraphTaxaChart$Taxa %>% length # should be 2 infection status x 7 timepoints x 10 taxa = 140 in length
GraphTaxaChart$Taxa %>% unique %>% length
GraphTaxaChart$Taxa %>% unique %>% sort
OrderColorList <- c("Pleosporales"="#8b4513","Capnodiales"="#ee82ee","Hypocreales"="#006400","Magnaporthales"="#808000","Chaetothyriales"="#483d8b","Tremellales"="#000080","Glomerellales"="#9acd32","Helotiales"="#20b2aa","Trichosphaeriales"="#8b008b","Xylariales"="#ff4500","Diaporthales"="#ffa500","Agaricales"="#ffff00","Auriculariales"="#7cfc00","Cantharellales"="#8a2be2","Sordariales"="#00ff7f","Tremellodendropsidales"="#dc143c","Tubeufiales"="#00bfff","Phallales"="#ff00ff","Chaetosphaeriales"="#1e90ff","Polyporales"="#db7093","Pseudodactylariales"="#eee8aa","Trechisporales"="#ff1493","Geastrales"="#ffa07a")

head(GraphTaxaChart)

OrderPlot <-
ggplot(GraphTaxaChart,aes(y=PercentReads,x=as.factor(Infection),fill=Taxa))+theme_bw()+
  geom_hline(yintercept = 0)+
  theme(panel.grid = element_blank())+
  geom_bar(stat = "identity",color='black') +
  ylab('Percent of reads')+xlab('')+
  # theme(legend.position = 'bottom')+  guides(fill = guide_legend(ncol = 7))+
  theme(legend.position = 'right')+  guides(fill = guide_legend(ncol = 1))+
  scale_fill_manual(values=OrderColorList,name='Order')+
  # scale_fill_manual(values=OrderColorList,name='')+
  facet_wrap(~MonthNumber,nrow=1,labeller=labeller(MonthNumber=MonthConversion)) #+ theme(legend.position = 'none')


# ggsave(filename='OrderPlot.jpg',plot=OrderPlot,dpi=1000,height=7,width=10)


##########
TaxaThreeChart <- 
plot_grid(CombinedVolPlot,
          OrderPlot+ggtitle('Fungal orders'),
          ncol=1, labels=c("","C"),rel_heights = c(1, 1)
          )

ThreeChartWithLegend <-
  ggdraw()+ 
  draw_plot(TaxaThreeChart, 0, 0, 1, 1) +
  draw_plot(VolLegend, -.1, .52, 1, height=0.125)

# ggsave(filename='temp2.jpg',plot=ThreeChartWithLegend,dpi=1000,height=12,width=10)
getwd()
