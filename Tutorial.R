##########################################
### Bison-Fly: Plant Breeding Pipeline ###
##########################################

### Necessary packages ###
library(FIELDimageR)
library(raster)
library(rgdal)
library(ggplot2)
library(DescTools)
library(lme4)
library(emmeans)
library(reshape2)
library(car)
library(plyr)
library(factoextra)
library(ggrepel)
library(agricolae)
library(corrplot)
library(RStoolbox)
library(gridExtra)

#####################################
### NDSU - Spring Wheat UAV Data  ###
#####################################

### List of orthomosaics ###
Field <- list.files("./5band/") # 15 5band orthomosaics
Field_RGB <- list.files("./RGB/") # 15 RGB orthomosaics
Field_DSM <- list.files("./DSM/") # 15 DSM orthomosaics

######################################################
### Using one orthomosaic as base to draw the grid ###
######################################################

### Suggested strategy of using lids to highlight where to click ###
Test <- stack(paste("./5band/",Field[4],sep = ""))
plotRGB(FIELDimageR:::RGB.rescale(Test,3))

##########################
### Plot Polygons Grid ###
##########################

Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Map <- read.csv("EX_MAP.csv",header = F,fileEncoding="UTF-8-BOM")
# x11()
Shapefile <- fieldShape(mosaic = Test,
                        ncols = 11,
                        nrows = 20,
                        fieldData = Data,
                        ID = "PLOT",
                        fieldMap = Map)
Shapefile<-Shapefile$fieldShape
dev.off()
plotRGB(FIELDimageR:::RGB.rescale(Test,3))
plot(Shapefile, border="green",add=T)

### Saving Shapefile ###
# library(rgdal)
# writeOGR(Shapefile, ".", "Shapefile", driver="ESRI Shapefile")
# Shapefile <- readOGR("Shapefile.shp") # Reading the saved shapefile.

### If you downloaded the example in the descriptions above use this code to read the shapefile:
# unzip("Shapefile.zip")
# Shapefile <- readOGR("./Shapefile/Shapefile.shp") # Reading the saved shapefile.

############
### Mask ###
############

### Removing soil if necessary ###
Mask <- fieldMask(Test,
                  Blue = 1, Green = 2, Red = 3, RedEdge = 4, NIR = 5, # Layers position according to each sensor
                  index = "NDVI",
                  cropValue = 0.7,
                  cropAbove = FALSE)

################################################
### Extracting data for 14 flights in a loop ###
################################################

# DataTotal<-NULL
# for(i in 2:length(Field)){
#   print(paste("===",Field[i],"==="))
#   EX <- stack(paste("./5band/",Field[i],sep = ""))
#   EX.1 <- fieldMask(EX,plot = F)
#   EX.I <- fieldIndex(mosaic = EX,Red = 3,Green = 2,Blue = 1,RedEdge = 4,NIR = 5,
#                      index = c("NGRDI","BGI","GLI","NDVI","NDRE","CIG","CIRE"),plot = F)
#   crs(Shapefile)<-crs(EX.I)
#   EX.I<- fieldInfo(mosaic = EX.I,
#                    fieldShape = Shapefile,
#                    buffer = -0.05,
#                    n.core = 3)
#   # Canopy Cover
#   EX.I<-fieldArea(mosaic = EX.1$mask, 
#                   fieldShape = EX.I$fieldShape,plot = F)
#   # EPH
#   DSM0 <- stack(paste("./DSM/",Field_DSM[1],sep = ""))
#   DSM1 <- stack(paste("./DSM/",Field_DSM[i],sep = ""))
#   
#   # Canopy Height Model (CHM):
#   DSM0 <- resample(DSM0, DSM1)
#   CHM <- DSM1-DSM0
#   CHM <- fieldMask(CHM, mask = Mask$mask, plot=F)
#   CHM <- CHM$newMosaic
#   
#   # Extracting the estimate plant height average (EPH):
#   EPH <- fieldInfo(CHM, fieldShape = EX.I$fieldShape, fun = "quantile",n.core = 3,plot = F)
#   DataTotal1<-data.frame(Date=Field[i],EPH$fieldShape@data)
#   EPH.1090.A<-extract(x = CHM, y = EPH$fieldShape)
#   EPH.1090<-do.call(rbind,lapply(EPH.1090.A, quantile, probs = c(0.1,0.9), na.rm=TRUE))
#   DataTotal1$'Height_10'<-EPH.1090[,1]
#   DataTotal1$'Height_90'<-EPH.1090[,2]
#   
#   # Data Table:
#   DataTotal<-rbind(DataTotal,DataTotal1)
#   
#   # Making plots:
#   fieldPlot(EX.I$fieldShape,
#             mosaic = EX,
#             color = c("red","green"),
#             #min.lim = 0.45,max.lim = 0.75,
#             fieldAttribute = "NDVI")
# }
# 
# ### Correcting Names ###
# Data.names<-gsub("layer",'Height_0',colnames(DataTotal))
# Data.names<-gsub("NA..1",'Height_50',Data.names)
# Data.names<-gsub("NA..2",'Height_75',Data.names)
# Data.names<-gsub("NA..3",'Height_100',Data.names)
# Data.names<-gsub("\\NA.",'Height_25',Data.names)
# Data.names<-gsub("objArea",'Canopy',Data.names)
# colnames(DataTotal)<-Data.names
# DataTotal<-subset(DataTotal, select = -ID.1)
# DataTotal$DAP<-as.numeric(do.call(rbind,strsplit(DataTotal$Date,split = "_"))[,4])
# head(DataTotal)

### Saving extracted data in a .CSV ###
# write.csv(DataTotal,"DataTotal.csv",row.names = F,col.names = T)
# DataTotal<-read.csv("DataTotal.csv",header = T)

##########################
### Agronomical Traits ###
##########################

### Field data collected manually ###
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)

### Mixed model: getting adjusted means and heritability (H2) ###
Trait<-c("MAT_DAY","HT","DH","LODG","YLD")
H2.AG<-NULL
for(t in 1:length(Trait)){
  Data1<-droplevels(Data[!is.na(Data[,colnames(Data)==Trait[t]]),])
  # mod<-lmer(eval(parse(text = paste(Trait.UAV[t2]," ~ RANGE+ROW+(1|NAME)",sep=""))),data = Data)
  mod<-lmer(eval(parse(text = paste(Trait[t],"~(1|NAME)",sep=""))),data = Data1)
  Var1<-as.data.frame(VarCorr(mod))$vcov
  names(Var1)<-as.data.frame(VarCorr(mod))$grp
  H2.AG<-rbind(H2.AG,data.frame(Trait=Trait[t],
                                H2=round(c(Var1[1]/sum(Var1[1],
                                                       Var1[2]/2 #2 replicates
                                )),3)
  ))
  mod<-lm(eval(parse(text = paste(Trait[t],"~NAME",sep=""))),data = Data1)
  Adj.Mean<-emmeans(mod, ~ NAME)
  if(t==1){
    Pheno.AG<-as.data.frame(Adj.Mean)[,c(1,2)]
  }
  if(t!=1){
    Pheno.AG<-merge(Pheno.AG,as.data.frame(Adj.Mean)[,c(1,2)],by="NAME")
  }
}
colnames(Pheno.AG)<-c("NAME",Trait)
head(Pheno.AG)

### Agronomical traits heritability ###
ggplot(data = H2.AG, 
       aes(x = Trait,
           y = H2*100,
           fill=as.factor(Trait))) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_grey(start=0.2, end=0.8)+
  ylim(c(0,100))+
  labs(y="H2 (%)",
       x="", 
       fill="Agronomical Traits") +
  geom_text(aes(label=paste(H2*100,"%")),size=5, position=position_dodge(width=0.9), vjust=-0.25)+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 

##################
### UAV Traits ###
##################

### UAV data extracted above (DataTotal) ###
DataTotal<-read.csv("DataTotal.csv",header = T)
DataTotal$RANGE<-as.factor(DataTotal$RANGE)
DataTotal$ROW<-as.factor(DataTotal$ROW)
DataTotal$NAME<-as.factor(DataTotal$NAME)

### Preparing information for running the statistic models below in a loop ###
DAP<-unique(DataTotal$DAP)
Trait.UAV<-c("Blue","Green","Red","RedEdge","NIR", # Single Bands
             "NGRDI","BGI","GLI","NDVI","NDRE","CIG","CIRE", # Vegetation indices
             "Canopy", # Canopy cover
             "Height_0","Height_10","Height_25","Height_50","Height_75","Height_90","Height_100" #Estimated Plant Height
)

### Mixed model: getting adjusted means and heritability (H2) ###
H2.UAV<-NULL
Pheno.UAV<-list()
for(t1 in 1:length(DAP)){
  Data<-droplevels(DataTotal[DataTotal$DAP==DAP[t1],])
  for(t2 in 1:length(Trait.UAV)){
    # mod<-lmer(eval(parse(text = paste(Trait.UAV[t2]," ~ RANGE+ROW+(1|NAME)",sep=""))),data = Data)
    mod<-lmer(eval(parse(text = paste(Trait.UAV[t2]," ~ (1|NAME)",sep=""))),data = Data)
    H2.UAV<-rbind(H2.UAV,cbind(DAP=DAP[t1],
                               Trait=Trait.UAV[t2],
                               H2=round(as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov[1],
                                                                                as.data.frame(VarCorr(mod))$vcov[2]/2),3)))
    # mod<-lm(eval(parse(text = paste(Trait.UAV[t2],"~RANGE+ROW+NAME",sep=""))),data = Data)
    mod<-lm(eval(parse(text = paste(Trait.UAV[t2],"~NAME",sep=""))),data = Data)
    Adj.Mean<-emmeans(mod, ~ NAME)
    if(t2==1){
      Pheno.UAV.1<-as.data.frame(Adj.Mean)[,c(1,2)]
      colnames(Pheno.UAV.1)<-c("NAME",Trait.UAV[1:t2])
    }
    if(t2!=1){
      Pheno.UAV.1<-merge(Pheno.UAV.1,as.data.frame(Adj.Mean)[,c(1,2)],by="NAME")
      colnames(Pheno.UAV.1)<-c("NAME",Trait.UAV[1:t2])
    }
  }
  Pheno.UAV[[t1]]<-Pheno.UAV.1
}
names(Pheno.UAV)<-DAP

### UAV traits heritability ###
H2.UAV<-as.data.frame(H2.UAV)
H2.UAV$H2<-as.numeric(as.character(H2.UAV$H2))
H2.UAV$DAP<-as.numeric(as.character(H2.UAV$DAP))
H2.UAV$Trait<-factor(H2.UAV$Trait,levels = Trait.UAV)
ggplot(data = H2.UAV, 
       aes(x = DAP,
           y = H2*100)) +
  facet_wrap(~Trait)+
  geom_line(size=1) +
  geom_point(size=2)+
  ylim(c(0,100))+
  labs(y="H2 (%)",
       x="Days After Planting (DAP)") +
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=10),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white"))

##################################
### Area Under the Curve (AUC) ###
##################################

### Choosing some days after planting (DAP) to investigate ###
unique(DataTotal$DAP)[order(unique(DataTotal$DAP))] # Options: c(16,19,25,31,34,41,46,49,52,55,67,73,74,80) 
DAP<-c(19,34,46,55,73,80)
Data<-DataTotal[DataTotal$DAP%in%DAP,]

### Extracted data visualization per DAP - NDVI ###
ggplot(Data, 
       aes(x = NDVI,
           fill=as.factor(DAP))) +
  geom_density(alpha=.5,position = 'identity') +
  facet_wrap(~DAP, ncol = 1)+
  scale_fill_grey(start=1, end=0)+
  labs(y="#genotypes",
       x="NDVI", 
       fill="DAP") +
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=14),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=14),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 

#########################
### AUC Visualization ###
#########################

### Choosing genotypes to highlight ###
unique(Data$NAME) # "FALLER","SYINGMAR","NDVITPRO"
Data1<-Data[as.character(Data$NAME)%in%c("FALLER","SYINGMAR","NDVITPRO"),]
Data2<-ddply(Data1,NAME~DAP,summarise,NDVI=mean(NDVI))
ggplot(data=Data2, aes(x=as.numeric(DAP), y= NDVI, col= NAME, group=NAME)) +
  geom_point(size=6)+
  geom_line(size=1.2) +
  scale_color_grey(start=0.8, end=0.2)+
  labs(x="Days After Planting (DAP)", fill="", col="")+
  theme_linedraw()+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=18),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 

#######################
### Calculating AUC ###
#######################

### Choosing some UAV traits ###
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90") 
DataAUC<-fieldAUC(data = Data,
                  trait = Trait,
                  keep.columns = colnames(Data)[2:27],
                  frame = "long")
DataAUC$AUC<-as.numeric(as.character(DataAUC$AUC))
DataAUC$TRAIT<-factor(DataAUC$TRAIT,levels = Trait)
DataAUC$NAME<-as.factor(DataAUC$NAME)
DataAUC$RANGE<-as.factor(DataAUC$RANGE)
DataAUC$ROW<-as.factor(DataAUC$ROW)
write.csv(DataAUC,"DataAUC.csv",row.names = F,col.names = T)
ggplot(DataAUC, aes(x = AUC)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.4,position = 'identity', fill="gray") +
  facet_wrap(~TRAIT,scales = "free")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

####################
### AUC analysis ###
####################

DataAUC<-read.csv("DataAUC.csv",header = T)
DataAUC$RANGE<-as.factor(DataAUC$RANGE)
DataAUC$ROW<-as.factor(DataAUC$ROW)
DataAUC$NAME<-as.factor(DataAUC$NAME)

### Mixed model: getting adjusted means and heritability (H2) ###
H2.AUC<-NULL
for(t in 1:length(Trait)){
  Data1<-droplevels(DataAUC[as.character(DataAUC$TRAIT)==Trait[t],])
  # mod<-lmer(AUC~RANGE+ROW+(1|NAME),data = Data1)
  mod<-lmer(AUC~(1|NAME),data = Data1)
  Var1<-as.data.frame(VarCorr(mod))$vcov
  names(Var1)<-as.data.frame(VarCorr(mod))$grp
  H2.AUC<-rbind(H2.AUC,data.frame(Trait=Trait[t],
                                  H2=round(c(Var1[1]/sum(Var1[1],
                                                         Var1[2]/2 #2 replicates
                                  )),3)
  ))
  # mod<-lm(AUC~RANGE+ROW+NAME,data = Data1)
  mod<-lm(AUC~NAME,data = Data1)
  Adj.Mean<-emmeans(mod, ~ NAME)
  if(t==1){
    Pheno.AUC<-as.data.frame(Adj.Mean)[,c(1,2)]
    colnames(Pheno.AUC)<-c("NAME",Trait[1:t])
  }
  if(t!=1){
    Pheno.AUC<-merge(Pheno.AUC,as.data.frame(Adj.Mean)[,c(1,2)],by="NAME")
    colnames(Pheno.AUC)<-c("NAME",Trait[1:t])
  }
}
head(Pheno.AUC)

### AUC traits heritability ###
H2.AUC<-as.data.frame(H2.AUC)
H2.AUC$Trait<-factor(H2.AUC$Trait,Trait)
H2.AUC$H2<-as.numeric(H2.AUC$H2)
ggplot(data = H2.AUC, 
       aes(x = Trait,
           y = H2*100,
           fill=as.factor(Trait))) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_grey(start=0.2, end=0.8)+
  ylim(c(0,100))+
  labs(y="H2 (%)",
       x="", 
       fill="UAV_AUC") +
  geom_text(aes(label=paste(H2*100,"%")), position=position_dodge(width=0.9), vjust=-0.25)+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white"))

##########################################
### Principal component analysis (PCA) ###
##########################################

### Merging Agro and UAV adjusted means ###
Pheno.PCA<-merge(Pheno.AG,Pheno.AUC,by="NAME")
Pheno.PCA.1<-Pheno.PCA[,c("MAT_DAY","HT","DH","LODG","YLD",
                          "NDRE","Canopy","Height_90")]
rownames(Pheno.PCA.1)<-Pheno.PCA$NAME

### PCA ###
Pheno.PCA.2 <- prcomp(Pheno.PCA.1, center = TRUE, scale = TRUE)

### Highlighting checks ###
checks<-c("NDSW0932","NDSW14098","NDVITPRO","SYINGMAR","ALPINE","BARLOW","ELGIN-ND","FALLER","GLENN","MAX")
groups <- as.character(Pheno.PCA$NAME)
groups[groups%in%checks]<-"Checks"
groups[groups!="Checks"]<-"Lines"
groups.text <- as.character(Pheno.PCA$NAME)
groups.text[!groups.text%in%checks]<-""

### PCA visualization ###
fviz_pca_biplot(Pheno.PCA.2,
                col.ind = groups, # color by groups
                legend.title = "",
                palette = c("black", "gray55"),
                repel = TRUE,
                geom.ind = "point",mean.point=F,
                pointshape = 21,
                pointsize = 3,
                fill.ind = groups
)+
  geom_text_repel(aes(label = groups.text),size = 3.5)+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=14),
        # axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=14),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white"))

###################
### Correlation ###
###################

### Specific flight ###
# 31 DAP (MAT_DAY)
# 46 DAP (YLD)
# 49 DAP (DH)
# 74 DAP (LODG) and 73 DAP (1 day before heavy storm)

### 73 DAP ###
Pheno.UAV.2<-Pheno.UAV$`73`[,c("NAME","NGRDI", "NDRE", "CIRE","Canopy","Height_50","Height_90")]
Pheno.COR<-merge(Pheno.AG,Pheno.UAV.2,by="NAME")
Pheno.COR.1<-scale(Pheno.COR[,-1],scale = T)
rownames(Pheno.COR.1)<-Pheno.COR[,1]

### r (73 DAP) ###
r<-correlation(Pheno.COR.1)
r$correlation
round(r$pvalue,2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r$correlation, 
         p.mat = r$pvalue,
         sig.level = 0.05, # Level of significance 5%
         method="color", col=col(200),  
         type="upper", order="alphabet",addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         insig = "blank", 
         diag=FALSE)
         
### 74 DAP ###
Pheno.UAV.2<-Pheno.UAV$`74`[,c("NAME","NGRDI", "NDRE", "CIRE","Canopy","Height_50","Height_90")]
Pheno.COR<-merge(Pheno.AG,Pheno.UAV.2,by="NAME")
Pheno.COR.1<-scale(Pheno.COR[,-1],scale = T)
rownames(Pheno.COR.1)<-Pheno.COR[,1]

### r (74 DAP) ###
r<-correlation(Pheno.COR.1)
r$correlation
round(r$pvalue,2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r$correlation, 
         p.mat = r$pvalue,
         sig.level = 0.05, # Level of significance 5%
         method="color", col=col(200),  
         type="upper", order="alphabet",addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         insig = "blank", 
         diag=FALSE)

### AUC ###
Pheno.COR<-merge(Pheno.AG,Pheno.AUC,by="NAME")
Pheno.COR.1<-scale(Pheno.COR[,-1],scale = T)
rownames(Pheno.COR.1)<-Pheno.COR[,1]

### r (AUC) ###
r<-correlation(Pheno.COR.1)
r$correlation
round(r$pvalue,2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r$correlation, 
         p.mat = r$pvalue,
         sig.level = 0.1,
         method="color", col=col(200),  
         type="upper", order="alphabet",addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         insig = "blank", 
         diag=FALSE)

###############
### Lodging ###
###############

### Flight before: 07/19/2021 (73 DAP) ###
### Flight after_1: 07/20/2021 (74 DAP) ###
### Flight after_2: 07/26/2021 (80 DAP) ###

### Choosing flights around the Lodging event ###
DAP<-c(73,74,80)
Data<-DataTotal[as.numeric(DataTotal$DAP)%in%DAP,]

### Choosing UAV traits to compare with LODG ###
Data<-Data[,c("DAP","LODG","Height_90")] # Other options: c("CIG","CIRE","Canopy","NGRDI","BGI","GLI","NDVI","NDRE","Height_50","Height_75","Height_90") 
Data$DAP<-as.factor(Data$DAP)

### Preparing data to make plots ###
Data.1<-melt(Data,
             value.name = "Index",
             measure.vars = c("Height_90"))
Data.2<-melt(Data.1,
             value.name = "Trait",
             measure.vars = c("LODG"))
colnames(Data.2)<-c("DAP","Index","Index.var","Trait","Trait.var")
Data.2$DAP<-as.factor(Data.2$DAP)
Data.2$Index<-as.factor(Data.2$Index)
Data.2$Trait<-as.factor(Data.2$Trait)
Data.2$Index.var<-as.numeric(as.character(Data.2$Index.var))
Data.2$Trait.var<-as.factor(as.character(Data.2$Trait.var))

### Simple boxplot visualization ###
ggplot(data = Data.2, 
       aes(y = Index.var,
           x = Trait.var,
           fill=Index)) + 
  facet_grid(DAP~Trait, scales = "free")+
  geom_boxplot() +
  labs(y="Estimate Plant Height",
       x="Lodging (Score)") +
  scale_fill_grey(start=0.8, end=0.2)+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=10),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 

##########################
### Maturity (MAT_DAY) ###
##########################

### Choosing flights around the MAT_DAY ###
DAP<-c(49,52,67,74,80)
Data<-DataTotal[as.numeric(DataTotal$DAP)%in%DAP,]

### Choosing UAV traits to compare with MAT_DAY ###
Data<-Data[,c("DAP","MAT_DAY","CIG","CIRE")] # Other options: c("CIG","CIRE","Canopy","NGRDI","BGI","GLI","NDVI","NDRE","Height_50","Height_75","Height_90")
Data$DAP<-as.factor(Data$DAP)

### Preparing data to make plots ###
Data.1<-melt(Data,
             value.name = "Index",
             measure.vars = c("CIG","CIRE"))
Data.2<-melt(Data.1,
             value.name = "Trait",
             measure.vars = c("MAT_DAY"))
colnames(Data.2)<-c("DAP","Index","Index.var","Trait","Trait.var")
Data.2$DAP<-as.factor(Data.2$DAP)
Data.2$Index<-as.factor(Data.2$Index)
Data.2$Trait<-as.factor(Data.2$Trait)
Data.2$Index.var<-as.numeric(as.character(Data.2$Index.var))
Data.2$Trait.var<-as.numeric(as.character(Data.2$Trait.var))

### Simple regression visualization ###
ggplot(data = Data.2, 
       aes(x = Index.var,
           y = Trait.var,
           colour=Index)) + 
  facet_grid(DAP~Trait, scales = "free")+
  geom_smooth(method=lm) + 
  geom_point(size = 2) +
  scale_color_grey(start=0.4, end=0.7)+
  labs(y="Maturity (day of the year)",
       x="Index") +
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=10),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 

##################################
### A. Covariates in the model ###
##################################

dev.off()

### 1) No covariate ###  
# Preparing the data:
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
# Mixed model:
# mod<-lmer(YLD~RANGE+ROW+(1|NAME),data = Data)
mod<-lmer(YLD~(1|NAME),data = Data)
# AIC (Comparing models):
(Yield.AIC<-AIC(mod))
# Residuals visualization:
qqPlot(residuals(mod))

# Making a loop:
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90")
Data.AIC<-NULL
for(i in 1:length(Trait)){
  
  ### 2) Single flight co-variate ###
  # Only one flight (e.g., 55 DAP)
  DataTotal<-read.csv("DataTotal.csv",header = T)
  DAP<-c(55) 
  # Preparing the data:
  Data<-DataTotal[DataTotal$DAP%in%DAP,]
  Data$RANGE<-as.factor(Data$RANGE)
  Data$ROW<-as.factor(Data$ROW)
  Data$NAME<-as.factor(Data$NAME)
  # Mixed model:
  # mod<-lmer(eval(parse(text = paste("YLD~",Trait[i],"+RANGE+ROW+(1|NAME)",sep=""))),data = Data)
  mod<-lmer(eval(parse(text = paste("YLD~",Trait[i],"+(1|NAME)",sep=""))),data = Data)
  # AIC (Comparing models):
  Data.AIC<-rbind(Data.AIC,cbind(Trait=Trait[i],AIC=AIC(mod), Model="55DAP")) 
  # Residuals visualization:
  qqPlot(residuals(mod)) 
  
  ### 3) AUC co-variate ### 
  # Preparing the data:
  DataAUC<-read.csv("DataAUC.csv",header = T)
  Data<-DataAUC[as.character(DataAUC$TRAIT)%in%Trait[i],]
  Data$RANGE<-as.factor(Data$RANGE)
  Data$ROW<-as.factor(Data$ROW)
  Data$NAME<-as.factor(Data$NAME)
  # Mixed model:
  # mod<-lmer(YLD~AUC+RANGE+ROW+(1|NAME),data = Data)
  mod<-lmer(YLD~AUC+(1|NAME),data = Data)
  # AIC (Comparing models):
  Data.AIC<-rbind(Data.AIC,cbind(Trait=Trait[i],AIC=AIC(mod), Model="AUC")) 
  # Residuals visualization:
  qqPlot(residuals(mod))
}
Data.AIC<-as.data.frame(Data.AIC)
Data.AIC$AIC<-as.numeric(Data.AIC$AIC)

### AIC visualization ###
ggplot(data = Data.AIC, 
       aes(x = Trait,
           y = AIC,
           fill=as.factor(Trait))) +
  geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~Model,scales = "free")+
  geom_hline(yintercept = Yield.AIC, col="red",linetype = "dashed", size=1)+
  scale_fill_grey(start=0.2, end=0.8)+
  labs(y="AIC",
       x="", 
       fill="UAV Traits") +
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 

################################################################
### B. UAV as a main trait in the model (Indirect Selection) ###
################################################################

### 1) Yield selection based on BLUP (observed data): ###
# Preparing the data:
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
# Mixed model:
# mod<-lmer(YLD~RANGE+ROW+(1|NAME),data = Data)
mod<-lmer(YLD~(1|NAME),data = Data)
# BLUPs:
BLUP.Pheno<-as.matrix(ranef(mod)$NAME)
# Ranking:
Sel.Pheno<-rownames(BLUP.Pheno)[order(BLUP.Pheno,decreasing = T)]

### 2) Yield selection based on BLUP (UAV AUC data): ###
DataAUC<-read.csv("DataAUC.csv",header = T)
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90")

# Making a loop:
Data.SC<-NULL
for(i in 1:length(Trait)){
  # Preparing the data:
  Data<-droplevels(DataAUC[as.character(DataAUC$TRAIT)%in%Trait[i],])
  Data$RANGE<-as.factor(Data$RANGE)
  Data$ROW<-as.factor(Data$ROW)
  Data$NAME<-as.factor(Data$NAME)
  # Mixed model:
  # mod<-lmer(AUC~RANGE+ROW+(1|NAME),data = Data)
  mod<-lmer(AUC~(1|NAME),data = Data)
  # BLUPs:
  BLUP.AUC<-as.matrix(ranef(mod)$NAME)
  # Ranking:
  Sel.AUC<-rownames(BLUP.AUC)[order(BLUP.AUC,decreasing = T)]
  # Selection coincidence (25%):
  n.sel<-round(length(Sel.Pheno)*0.25,0)
  Data.SC<-rbind(Data.SC,cbind(Trait=Trait[i],SC=sum(Sel.Pheno[1:n.sel]%in%Sel.AUC[1:n.sel])/n.sel))
}
Data.SC<-as.data.frame(Data.SC)
Data.SC$SC<-as.numeric(Data.SC$SC)

# Indirect selection coincidence:
ggplot(data = Data.SC, 
       aes(x = Trait,
           y = SC*100,
           fill=as.factor(Trait))) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_grey(start=0.2, end=0.8)+
  ylim(c(0,100))+
  labs(y="Indirect selection coincidence (%)",
       x="", 
       fill="UAV Traits") +
  geom_text(aes(label=paste(round(SC*100,2),"%")),size=5, position=position_dodge(width=0.9), vjust=-0.25)+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 

###############################################################################
### C. Indirect Selection using UAV traits and DAP as cofactor in the model ###
###############################################################################

### 1) Yield selection based on BLUP (observed data): ###
# Preparing the data:
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
# Mixed model:
# mod<-lmer(YLD~RANGE+ROW+(1|NAME),data = Data)
mod<-lmer(YLD~(1|NAME),data = Data)
# BLUPs:
BLUP.Pheno<-as.matrix(ranef(mod)$NAME)
# Ranking:
Sel.Pheno<-rownames(BLUP.Pheno)[order(BLUP.Pheno,decreasing = T)]

### 2) Yield selection based on BLUP (UAV data): ###
DataTotal<-read.csv("DataTotal.csv",header = T)
DAP<-c(31, 41, 49, 55, 67, 73)
DataTotal<-droplevels(DataTotal[DataTotal$DAP%in%DAP,])
# Preparing the data:
DataTotal$DAP<-as.factor(DataTotal$DAP)
DataTotal$RANGE<-as.factor(DataTotal$RANGE)
DataTotal$ROW<-as.factor(DataTotal$ROW)
DataTotal$NAME<-as.factor(DataTotal$NAME)
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90")

# Making a loop:
Data.SC<-NULL
for(i in 1:length(Trait)){
  # Mixed model:
  # mod<-lmer(eval(parse(text = paste(Trait[i]," ~ DAP+RANGE+ROW+(1|NAME)",sep=""))),data = DataTotal)
  mod<-lmer(eval(parse(text = paste(Trait[i]," ~ DAP+(1|NAME)",sep=""))),data = DataTotal)
  # BLUPs:
  BLUP.UAV<-as.matrix(ranef(mod)$NAME)
  # Ranking:
  Sel.UAV<-rownames(BLUP.UAV)[order(BLUP.UAV,decreasing = T)]
  # Selection coincidence (25%):
  n.sel<-round(length(Sel.Pheno)*0.25,0)
  Data.SC<-rbind(Data.SC,cbind(Trait=Trait[i],SC=sum(Sel.Pheno[1:n.sel]%in%Sel.UAV[1:n.sel])/n.sel))
}
Data.SC<-as.data.frame(Data.SC)
Data.SC$SC<-as.numeric(Data.SC$SC)

# Indirect selection coincidence:
ggplot(data = Data.SC, 
       aes(x = Trait,
           y = SC*100,
           fill=as.factor(Trait))) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_grey(start=0.2, end=0.8)+
  ylim(c(0,100))+
  labs(y="Indirect selection coincidence (%)",
       x="", 
       fill="UAV Traits") +
  geom_text(aes(label=paste(round(SC*100,2),"%")),size=5, position=position_dodge(width=0.9), vjust=-0.25)+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 
