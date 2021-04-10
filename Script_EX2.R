############################
### FIELDimageR pipeline ###
############################

################
### Packages ### 
################

library(FIELDimageR)
library(raster)
library(agricolae)
library(reshape2)
library(ggplot2)
library(lme4)
library(plyr)
library(DescTools)
library(ggrepel)

## Uploading hyperspectral file with 474 bands (EX_HYP.tif):
EX.HYP<- stack("EX_HYP.tif")

## Wavelengths (namesHYP.csv):
NamesHYP<- as.character(read.csv("namesHYP.csv")$NameHYP)
names(EX.HYP)<- NamesHYP
head(NamesHYP)

# RGB from hyperspectral:
R<- EX.HYP$X654.788879 # 654nm (Red)
G<- EX.HYP$X552.598572 # 552nm (Green)
B<- EX.HYP$X450.408295 # 450nm (Blue)
RGB<- stack(c(R,G,B))
plotRGB(RGB, stretch="lin")

## Removing soil using RGB (index NGRDI):
RGB.S<- fieldMask(RGB, index="NGRDI", cropValue = 0.0, cropAbove = F)

## Data frame with field information to make the Map:
Data<- as.data.frame(read.csv("DataHYP.csv"))
head(Data)

Map<- fieldMap(fieldPlot = as.character(Data$Plot),
              fieldRow = as.character(Data$Range),
              fieldColumn = as.character(Data$Row), decreasing = T)
Map

## Building plot shapefile using RGB as base (14 columns and 14 rows):
plotFile<- fieldShape(RGB.S, ncols = 14, nrows = 14, 
                     fieldMap = Map, fieldData = Data, ID = "Plot")
# plotFile<- fieldShape(RGB, ncols = 14, nrows = 14, 
#                       fieldMap = Map, fieldData = Data, ID = "Plot")

## Removing soil from the hyperspectral:
# EX.HYP.S<- fieldMask(EX.HYP,mask = RGB.S$mask, plot = F)
# writeRaster(EX.HYP.S$newMosaic, filename="EX_HYP_S.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
EX.HYP.S<- stack("EX_HYP_S.tif")
names(EX.HYP.S)<- NamesHYP
plot(EX.HYP.S$X750.592224, col=grey(100:1/100), axes=FALSE, box=FALSE)
plot(plotFile$fieldShape, add=T)

## Reducing resolution to accelerate data extraction por plot (e.g., average values):
EX.HYP.S<- aggregate(EX.HYP.S, 2)
plot(EX.HYP.S$X750.592224, col=grey(100:1/100), axes=FALSE, box=FALSE)
plot(plotFile$fieldShape, add=T)

## Extracting data (474 bands):
EX.HYP.I<- fieldInfo(EX.HYP.S, # EX.HYP.S$newMosaic, 
                    fieldShape = plotFile$fieldShape,
                    n.core = 1)

## Saving the new csv with hyperspectral information per plot:
DataHYP<- EX.HYP.I$fieldShape@data
colnames(DataHYP)<- c(colnames(DataHYP)[1:9], NamesHYP)
DataHYP[1:5, 1:12]
# write.csv(DataHYP,"DataHypNew.csv",col.names = T,row.names = F)

## Visualizing the extracted data from FIELDimageR:
dev.off()
DataHYP1<- EX.HYP.I$plotValue[,-1]

## Plotting the data:
plot(x=as.numeric(NamesHYP), y=as.numeric(DataHYP1[1,]), type = "l",
     xlab = "Wavelength (nm)", ylab = "Reflectance", 
     col="black", 
     ylim=c(0,0.8),
     lwd=2,cex.lab=1.2)
for(i in 2:dim(DataHYP1)[1]){
  lines(x=as.numeric(NamesHYP), y=as.numeric(DataHYP1[i,]), type = "l",
        col=i, lwd=2)
}
legend("topright",
       c("Blue (400-480nm)",
         "Green (480-600nm)",
         "Red (600-680nm)",
         "RedEdge (680-750nm)",
         "NIR (750-1300nm)",
         "SWIR-1 (1500-1800nm)",
         "SWIR-2 (2000-2400nm)"),
       col=c("black"), box.lty=0, cex = 0.8)

## Removing atmospheric absorption regions: 
DataHYP2<- DataHYP[,-c(1:9)]
Wavelenght<- as.numeric(colnames(DataHYP2))
DataHYP2<- DataHYP2[,!(Wavelenght>1290&Wavelenght<1510|
                         Wavelenght>1790&Wavelenght<2040|
                         Wavelenght>2350&Wavelenght<2550)]

## Vector normalization to reduce the effects of illumination conditions:
DataHYP2<- as.matrix(t(apply(DataHYP2, 1, function(x){
  x/sqrt(sum(as.numeric(x)^2))
})))

## Plotting the new data:
plot(x=as.numeric(colnames(DataHYP2)), y=as.numeric(DataHYP2[1,]), type ="p",
     xlab = "Wavelength (nm)", ylab = "Reflectance", 
     col="black", pch=19, cex.lab=1.2,
     #ylim=c(0,0.1),
     cex=0.5)
for(i1 in 2:dim(DataHYP2)[1]){
  lines(x=as.numeric(colnames(DataHYP2)), y=as.numeric(DataHYP2[i1,]), type ="p",
        col=i1, pch=19, cex=0.5)
}
legend("topright",
       c("Blue (400-480nm)",
         "Green (480-600nm)",
         "Red (600-680nm)",
         "RedEdge (680-750nm)",
         "NIR (750-1300nm)",
         "SWIR-1 (1500-1800nm)",
         "SWIR-2 (2000-2400nm)"),
       col =c("black"), box.lty=0, cex = 0.8)

## Saving the new csv table:
DataHYP3<- cbind(DataHYP[,c(1:9)], DataHYP2)
# write.csv(DataHYP3, "DataHypNew2.csv", col.names = T, row.names = F)

## Scale the wavelengths values to correct for any unit effects:
DataHYP4<- scale(DataHYP2,scale = T)

## Correlation between wavelengths and Trait_1:
r.Hyp<- NULL
for(r in 1:dim(DataHYP4)[2]){
  r.Hyp<- c(r.Hyp, cor(DataHYP4[,r], scale(DataHYP$Trait_1,scale = T),
                       use = "pairwise.complete.obs"))
}

## Plotting correlation values: 
plot(x=as.numeric(colnames(DataHYP4)), y=r.Hyp, type = "o",
     xlab = "Wavelength (nm)",
     ylab = "r", 
     col="black", pch=19, cex.lab=1.2
)
abline(h=0, col="red", lty=2, lwd=3)
legend("topright",
       c("Blue (400-480nm)",
         "Green (480-600nm)",
         "Red (600-680nm)",
         "RedEdge (680-750nm)",
         "NIR (750-1300nm)",
         "SWIR-1 (1500-1800nm)",
         "SWIR-2 (2000-2400nm)"),
       col=c("black"), box.lty=0, cex = 0.8)

## Preparing the data:
str(DataHYP3)
DataHYP3$Name<- as.factor(DataHYP3$Name)
DataHYP3$Block<- as.factor(DataHYP3$Block)
Trait<- c("Trait_1", paste0("`", colnames(DataHYP3)[-(1:9)],"`", sep="")) 

## Using the package "lme4" for mixed modeling:
H2<- NULL
for(t in 1:length(Trait)){
  mod<- lmer(eval(parse(text = paste(Trait[t],"~Block+(1|Name)", sep=""))), data = DataHYP3)
  H2<-c(H2, c(as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov)))
}

Trait_1.H2<- H2[1]
Trait_1.H2 # H2=81%

## Plotting heritability values:
plot(x=as.numeric(colnames(DataHYP2)), y=as.numeric(H2[-1]), type = "o",
     xlab = "Wavelength (nm)", ylab = "H2", 
     col="black", pch=19, cex.lab=1.2,
     ylim=c(0,0.9),
     cex=0.5)
abline(h=Trait_1.H2, col="blue", lty=2, lwd=3)
text(x = 1300, y = (Trait_1.H2+0.03), paste("Trait_1 (H2=", round(Trait_1.H2,2),")", sep=""), col="blue")
legend(list(x = 1600, y = (Trait_1.H2-0.03)),
       c("Blue (400-480nm)",
         "Green (480-600nm)",
         "Red (600-680nm)",
         "RedEdge (680-750nm)",
         "NIR (750-1300nm)",
         "SWIR-1 (1500-1800nm)",
         "SWIR-2 (2000-2400nm)"),
       col =c("black"), box.lty=0, cex = 0.8)

