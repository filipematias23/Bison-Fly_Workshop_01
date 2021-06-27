############################
### FIELDimageR pipeline ###
############################

################
### Packages ### 
################

library(FIELDimageR)
library(raster)
library(ggplot2)
library(agricolae)
library(reshape2)
library(lme4)
library(readxl)

# Uploading one image as example and decreasing the resolution
EX3<-stack("soybean/11.jpg")
# EX3<-aggregate(EX3, fact= 4)
plotRGB(EX3)

# Shapefile with extent=T (The whole image area will be the shapefile)
EX.shapeFile<-fieldPolygon(EX3,extent = T)

# Select one index to identify leaves and remove the background
EX3.I1<- fieldIndex(mosaic = EX3,index = c("SI","BGI","BI"))

# Thresholding
dev.off()
par(mfrow=c(1,2))
hist(EX3.I1$BGI)
plot(EX3.I1$BGI)

# Removing the background
EX3.R<- fieldMask(mosaic = EX3, index = "BGI",
                   cropValue = 0.7,
                   cropAbove = T)

# Counting the total number of seeds
EX.P.Total<-fieldCount(mosaic = EX3.R$mask, 
                       fieldShape = EX.shapeFile$fieldShape, 
                       minSize = 0.1,
                       cex = 1.5,
                       na.rm = T)

# Select one index to identify green seeds
EX3.I2<- fieldIndex(mosaic = EX3.R$newMosaic,index = c("SI","BGI","BI"))

#BI index
plot(EX3.I2$BI)

# Selecting green seeds
EX3.R2<- fieldMask(mosaic = EX3.R$newMosaic, 
                   index = "BI",
                   #myIndex = "Blue",
                   cropValue = 130,
                   cropAbove = T)

# Counting the number of green seeds
EX.Green<-fieldCount(mosaic = EX3.R2$mask, 
                     fieldShape = EX.shapeFile$fieldShape, 
                     minSize = 0.07,
                     cex = 1.5,
                     na.rm = T)

# Joying information
data.frame(Total=EX.P.Total$fieldCount,
           Green=EX.Green$fieldCount,
           Percentage=round(EX.Green$fieldCount/EX.P.Total$fieldCount,2))

################
### Parallel ###
################

# Required packages
library(parallel)
library(foreach)
library(doParallel)

# Images names (folder directory: "./soybean/")
pics<-list.files("./soybean/")

# Number of cores
n.core<-2

# Starting parallel
cl <- makeCluster(n.core, output = "")
registerDoParallel(cl)
EX.Table.Parallel <- foreach(i = 1:length(pics), .packages = c("raster","FIELDimageR"), 
                             .combine = rbind) %dopar% {
                               EX3<-stack(paste("./soybean/",pics[i],sep = ""))
                               # EX3<-aggregate(EX3, fact= 4)
                               EX.shapeFile<-fieldPolygon(EX3,extent = T,
                                                          plot = F)
                               EX3.R1<- fieldMask(mosaic = EX3, index = "BGI",
                                                  cropValue = 0.7,
                                                  cropAbove = T,
                                                  plot = F)
                               EX.P.Total<-fieldCount(mosaic = EX3.R1$mask, 
                                                      fieldShape = EX.shapeFile$fieldShape, 
                                                      minSize = 0.1,
                                                      cex = 1.5,
                                                      na.rm = T)
                               EX3.R2<- fieldMask(mosaic = EX3.R1$newMosaic, 
                                                  index = "BI",
                                                  cropValue = 130,
                                                  cropAbove = T,
                                                  plot = F)
                               EX.Green<-fieldCount(mosaic = EX3.R2$mask, 
                                                    fieldShape = EX.shapeFile$fieldShape, 
                                                    minSize = 0.07,
                                                    cex = 1.5,
                                                    na.rm = T)
                               data.frame(Total=EX.P.Total$fieldCount,
                                          Green=EX.Green$fieldCount,
                                          Percentage=round(EX.Green$fieldCount/EX.P.Total$fieldCount,2))
                             }
rownames(EX.Table.Parallel)<-pics
EX.Table.Parallel 

# Taking individual seed measurements (Remove artifacts by changing the parameter *minArea* and observing the values on EX3.D$Dimension$area)
dev.off()
EX3.D<-fieldObject2(mosaic = EX3.R$mask, 
                      watershed = T, 
                      minArea = 100)

# Measurement Output:
EX3.D$numObjects
EX3.D$Dimension

# Individual seed visualization:
plotRGB(EX3)
plot(EX.shapeFile$fieldShape, add=T)
plot(EX3.D$Objects, add=T, border="red")
plot(EX3.D$Polygons, add=T, border="blue")
plot(EX3.D$single.obj[[1]], add=T, col="yellow")
lines(EX3.D$x.position[[1]], col="red", lty=2)
lines(EX3.D$y.position[[1]], col="red", lty=2)

# Calculating indices per seed:
EX3.I<- fieldIndex(mosaic = EX3,index = c("SI","BGI","BI"))
EX3.Data<-fieldInfo(mosaic = EX3.I[[c("SI","BGI","BI")]], fieldShape = EX3.D$Objects, projection = F)
EX3.Data$fieldShape@data

# Perimeter:
# install.packages("spatialEco")
library(spatialEco)
perimeter<-polyPerimeter(EX3.D$Objects)
box<-polyPerimeter(EX3.D$Polygons)
Data.Obj<-cbind(EX3.Data$fieldShape@data,EX3.D$Dimension,perimeter=perimeter,box=box)
Data.Obj

# Data visualization: 
library(reshape2)
Data.Obj1<-melt(Data.Obj[,c("SI","BGI","BI","area","x.dist","y.dist","perimeter","box")])

ggplot(Data.Obj1, aes(x=value, fill=variable)) +
  geom_histogram(aes(y=..density..), colour="black")+
  geom_density(alpha=.2)+
  facet_wrap(~variable, scales = "free")

###########
### END ###
###########


