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

############
### Data ###
############

MOSAIC<-list.files("./MOSAIC/")
DSM<-list.files("./DSM/")

###################
### Basic steps ###
###################

# Uploading an example mosaic
Test <- stack(paste("./MOSAIC/",MOSAIC[3],sep = ""))

# Cropping and reducing siza
Test.Crop <- fieldCrop(mosaic = Test)

# Rotating field
Test.Rotate<-fieldRotate(Test.Crop, theta = 2.3)
#Test.Rotate<-fieldRotate(Test.Crop, clockwise = F)

# Removing soil and making a mask
Test.RemSoil<-fieldMask(Test.Rotate) 

##################
### Shape file ###
##################

# Reading FieldData.csv
Data<-read.csv(file = "EX1_Data.csv",header = T)

# Making the field Map
Map<-fieldMap(fieldPlot = Data$Plot,
              fieldColumn = Data$Column,
              fieldRow = Data$Row,
              decreasing = T)
Map

rotate <- function(x) t(apply(x, 2, rev))

Map<-rotate(rotate(Map))
Map

# Building the plot shapefile (ncols = 14 and nrows = 10)
plotShape<-fieldShape(mosaic = Test.RemSoil$newMosaic, 
                      ncols = 14, 
                      nrows = 10, 
                      fieldData = Data, 
                      ID = "Plot", 
                      fieldMap = Map)

##########################
### Vegetation indices ###
##########################

Test.Indices<- fieldIndex(mosaic = Test.RemSoil$newMosaic, 
                         Red = 1, Green = 2, Blue = 3, 
                         index = c("NGRDI","BGI", "GLI","VARI"), 
                         myIndex = c("(Red-Blue)/Green","2*Green/Blue"))

############################
### Extracting plot data ###
############################

Test.Info<- fieldInfo(mosaic = Test.Indices[[c("NGRDI","BGI", "GLI","VARI","myIndex.1","myIndex.2")]],
                      fieldShape = plotShape$fieldShape, 
                      n.core = 3)
Test.Info$fieldShape@data

plot(Test.Indices$myIndex.2)
plot(plotShape$fieldShape,add=T)

###############################
### Estimating plant height ###
###############################

# Uploading files from soil base and vegetative growth:
DSM0 <- stack(paste("./DSM/",DSM[1],sep = ""))
DSM1 <- stack(paste("./DSM/",DSM[3],sep = ""))

# Cropping the image using the previous shape from step 2:
DSM0.Crop <- fieldCrop(mosaic = DSM0,fieldShape = Test.Crop)
DSM1.Crop <- fieldCrop(mosaic = DSM1,fieldShape = Test.Crop)

# Canopy Height Model (CHM):
DSM0.R <- resample(DSM0.Crop, DSM1.Crop)
CHM <- DSM1.Crop-DSM0.R
plot(CHM)

# Rotating the image using the same theta from step 3:
CHM.Rotate<-fieldRotate(CHM, theta = 2.3)

# Removing the soil using mask from step 4:
CHM.RemSoil <- fieldMask(CHM.Rotate, mask = Test.RemSoil$mask)

# Observing EPH profile for 2 ranges:
CHM.Draw <- fieldDraw(mosaic = CHM.Rotate,
                            ndraw = 2)
dev.off()
par(mfrow=c(1,3))
plot(x = CHM.Draw$Draw1$drawData$x-100230, 
     y = CHM.Draw$Draw1$drawData$layer, 
     type="l", col="red",lwd=1,ylim=c(0,1),
     xlab="Distance (m)", ylab="EPH (m)")
lines(x = CHM.Draw$Draw2$drawData$x-100230, 
      y = CHM.Draw$Draw2$drawData$layer,
      type="l", col="blue",lwd=1,add=T)

plot(CHM.Rotate, col = grey(1:100/100), axes = FALSE, box = FALSE,legend=F)
lines(CHM.Draw$Draw1$drawData$x,CHM.Draw$Draw1$drawData$y, type="l", col="red",lwd=2)
lines(CHM.Draw$Draw2$drawData$x,CHM.Draw$Draw2$drawData$y, type="l", col="blue",lwd=2)

plotRGB(Test.Rotate)
lines(CHM.Draw$Draw1$drawData$x,CHM.Draw$Draw1$drawData$y, type="l", col="red",lwd=2)
lines(CHM.Draw$Draw2$drawData$x,CHM.Draw$Draw2$drawData$y, type="l", col="blue",lwd=2)

# Extracting the estimate plant height average (EPH):
EPH <- fieldInfo(CHM.RemSoil$newMosaic, 
                 fieldShape = Test.Info$fieldShape, 
                 fun = "mean") #fun="quantile"
EPH$plotValue

########################################
### Evaluating all mosaics in a loop ###
########################################

DataTotal<-NULL
for(i in 2:length(MOSAIC)){
  EX1 <- stack(paste("./MOSAIC/",MOSAIC[i],sep = ""))
  EX1.Crop <- fieldCrop(mosaic = EX1,
                    fieldShape = Test.Crop,
                    plot = F)
  EX1.Rotate<-fieldRotate(EX1.Crop,
                    theta = 2.3,
                    plot = F)
  EX1.RemSoil<-fieldMask(EX1.Rotate, plot = F)
  EX1.Indices<- fieldIndex(mosaic = EX1.RemSoil$newMosaic, 
                            Red = 1, Green = 2, Blue = 3, 
                            index = c("NGRDI"), # c("NGRDI","BGI", "GLI","VARI")
                            myIndex = c("(Red-Blue)/Green"), # c("(Red-Blue)/Green","2*Green/Blue")
                            plot = F)
  EX1.Info<- fieldInfo(mosaic = EX1.Indices[[c("NGRDI","myIndex")]], # c("NGRDI","BGI", "GLI","VARI","myIndex.1","myIndex.2")
                        fieldShape = plotShape$fieldShape, 
                        n.core = 3)
  DSM0 <- stack(paste("./DSM/",DSM[1],sep = ""))
  DSM1 <- stack(paste("./DSM/",DSM[i],sep = ""))
  DSM0.Crop <- fieldCrop(mosaic = DSM0,fieldShape = Test.Crop,plot = F)
  DSM1.Crop <- fieldCrop(mosaic = DSM1,fieldShape = Test.Crop,plot = F)
  DSM0.R <- resample(DSM0.Crop, DSM1.Crop)
  CHM <- DSM1.Crop-DSM0.R
  CHM.Rotate<-fieldRotate(CHM, theta = 2.3,plot = F)
  CHM.RemSoil <- fieldMask(CHM.Rotate, mask = EX1.RemSoil$mask,plot = F)
  EPH <- fieldInfo(CHM.RemSoil$newMosaic, 
                   fieldShape = EX1.Info$fieldShape, 
                   fun = "mean")
  DataTotal<-rbind(DataTotal,
                   data.frame(DAP=as.character(do.call(c,strsplit(MOSAIC[i],split = "_"))[2]),
                              EPH$fieldShape@data))
  # Making map plots
  fieldPlot(fieldShape=EPH$fieldShape,
            fieldAttribute="NGRDI", 
            mosaic=EX1.Rotate, 
            color=c("red","green"),
            # min.lim = 0,
            # max.lim = 0.35,
            alpha = 0.5,
            round = 2)
  print(paste("### Completed: ", "Mosaic_",i," ###",sep=""))
  }

colnames(DataTotal)<-c(colnames(DataTotal)[-c(dim(DataTotal)[2])],"EPH") # layer=EPH
DataTotal<-DataTotal[,!colnames(DataTotal)%in%c("ID","ID.1","PlotName")] # Removing column 12 ("ID.1")
DataTotal

#write.csv(DataTotal,"DataTotal.csv",row.names = F,col.names = T)

################
### Graphics ###
################

DataTotal$Name<-as.factor(as.character(DataTotal$Name))
DataTotal$Row<-as.factor(as.character(DataTotal$Row))
DataTotal$Column<-as.factor(as.character(DataTotal$Column))
DataTotal$DAP<-as.numeric(as.character(DataTotal$DAP))
DataTotal$NGRDI<-as.numeric(as.character(DataTotal$NGRDI))
DataTotal$EPH<-as.numeric(as.character(DataTotal$EPH))

ggplot(DataTotal, aes(x = NGRDI,fill=as.factor(DAP))) +
  geom_density(alpha=.5,position = 'identity') +
  facet_wrap(~DAP,ncol = 1)+
  scale_fill_grey(start=1, end=0)+
  labs(y="#genotypes",x="NGRDI", fill="DAP") +
  theme_bw() 

####################
### Heritability ###
####################

DAP<-unique(DataTotal$DAP)

H2<-NULL
for(h in 1:length(DAP)){
  mod<-lmer(NGRDI~Row+Column+(1|Name),DataTotal[as.character(DataTotal$DAP)==DAP[h],])
  H2.a<-c("NGRDI",DAP[h],as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov))
  
  mod<-lmer(EPH~Row+Column+(1|Name),DataTotal[as.character(DataTotal$DAP)==DAP[h],])
  H2.b<-c("EPH",DAP[h],as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov))
  
  H2<-rbind.data.frame(H2,rbind(H2.a,H2.b))
}
colnames(H2)<-c("Trait","DAP","H2")
H2$H2<-as.numeric(as.character(H2$H2))
H2$DAP<-as.numeric(as.character(H2$DAP))

ggplot(H2,aes(x=as.factor(DAP),y=H2,fill=as.factor(DAP)))+
  geom_bar(stat="identity")+
  facet_wrap(~Trait)+
  scale_fill_grey(start=0.8, end=0.2)+
  labs(x="Days After Planting (DAP)", fill="")+
  geom_text(aes(label=round(H2,2)), vjust=1.6, color="white", size=6)+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))


##################################
### Area under the curve (AUC) ###
##################################

DataTotal1<-DataTotal[as.character(DataTotal$Name)%in%c("G43","G44","G45"),]

ggplot(data=DataTotal1, aes(x=as.numeric(DAP), y= NGRDI, col= Name, group=Name)) +
  geom_point(size=6)+
  geom_line(size=1.2) +
  scale_color_grey(start=0.8, end=0.2)+
  labs(x="Days After Planting (DAP)", fill="", col="")+
  theme_linedraw()

Trait<-c("NGRDI") # c("GLI","EPH")
Plot<-as.character(unique(DataTotal$Plot))

DataAUC<-NULL
for(a1 in 1:length(Plot)){
  D1<-DataTotal[as.character(DataTotal$Plot)==Plot[a1],]
  x1<-c(0,as.numeric(D1$DAP))
  y1<-c(0,as.numeric(D1[,Trait]))
  DataAUC <- rbind(DataAUC,
                   c(NGRDI_AUC=AUC(x = x1[!is.na(y1)], y = y1[!is.na(y1)]),
                     Name=unique(as.character(D1$Name)),
                     Trait=unique(D1$Trait),
                     Row=unique(D1$Row),
                     Column=unique(D1$Column)))}

DataAUC<-as.data.frame(DataAUC)
DataAUC$NGRDI_AUC<-as.numeric(as.character(DataAUC$NGRDI_AUC))
DataAUC$Name<-as.factor(DataAUC$Name)
DataAUC$Row<-as.factor(DataAUC$Row)
DataAUC$Column<-as.factor(DataAUC$Column)
DataAUC

### AUC Heritability ###

mod<-lmer(NGRDI_AUC~Row+Column+(1|Name),DataAUC)
H2<-as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov)
H2

ggplot(DataAUC, aes(x = NGRDI_AUC)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.5,position = 'identity', fill="cadetblue") +
  labs(y="#genotypes",x=paste("Area under the curve (AUC_NGRDI: H2=",round(H2,2),")",sep="")) +
  theme_bw() 

### Linear Regression ###

DataTotal.Reg<-subset(DataTotal,DAP=="40")
DataTotal.Reg$Check<-as.character(DataTotal.Reg$Name)
DataTotal.Reg$Check[!DataTotal.Reg$Check%in%c("G43","G44","G45")]<-""

ggplot(DataTotal.Reg,aes(y=NGRDI, x=EPH)) + 
  geom_point() +
  geom_smooth(method=lm)+
  labs(y="EPH",x="NGRDI",fill="",alpha="")+
  geom_vline(aes(xintercept=0.02),col="red", linetype = 2, size=0.7) +
  theme_bw()+
  geom_text_repel(aes(label = Check),
                  size = 3.5, col="red",
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

















