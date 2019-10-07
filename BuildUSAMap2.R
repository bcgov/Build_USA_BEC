## Script to create RF model from training points and predict + map US zones/ subzones (predicted within each zone)
### Kiri Daust, based on Will's original script Build_BGC_RndFor_Model_ver7b but has the following differences
## Loops to predict subzones within zones
### Translates predictions to a grid point map
.libPaths("E:/R packages351")

require(rattle)
require(rpart)
require(rpart.plot)
require(plyr)
require(reshape)
require(reshape2)
#require(VSURF)
#require(reports)
#require(rpart.utils)
#require(rfUtilities)
require(parallel)
require(foreach)
require(doParallel)
require(ggplot2)
#require(functional)
require(plot3D)
require(dplyr)
require(tcltk)
require(caret)
require(randomForest)
#require(ranger)
#require (OpenMP)
#require (randomForestSRC)
require (tools)
require(data.table)
require(spatstat)
require(spatialEco)
require(survey)
require(scales)
require(UBL)

rm(list=ls())
options(stringsAsFactors = FALSE)

wd=tk_choose.dir(); setwd(wd)
##########Read USA ClimateBC training dataset to build model
#fname <- "July12_USATrainingPts_AllDat.csv"
##New training points from initial
#fname <- "USA_GridTraining_500pts_Normal_1961_1990MSY.csv"
fplot=(file.choose()) 
Columns = c("ID1", "ID2", "Latitude", "Longitude", "Elevation", "AHM", "bFFP",
            "CMD07","DD5_sp","EMT","Eref_sm","EXT","FFP","MCMT","MSP",
            "PPT07","PPT08", "PPT05","PPT06","PPT09", "SHM","TD","Tmax_sp","Tmin_at",
            "Tmin_sm","Tmin_wt", "PPT_at","PPT_wt", "PAS","eFFP",
            "Eref09","MAT","Tmin_sp","CMD", "Tmin02")#"GCM", 

X1 <- fread(fplot,  stringsAsFactors = FALSE,data.table = FALSE)#, select = Columns
X1$PPT_MJ <- X1$PPT05 + X1$PPT06 # MaY/June precip
X1$PPT_JAS <- X1$PPT07 + X1$PPT08 + X1$PPT09 # July/Aug/Sept precip
X1$PPT.dormant <- X1$PPT_at + X1$PPT_wt # for calculating spring deficit
X1$CMD.def <- 500 - (X1$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
X1$CMD.def [X1$CMD.def < 0] <- 0 #negative values set to zero = no deficit
X1$CMDMax <- X1$CMD07
X1$CMD.total <- X1$CMD.def + X1$CMD
X1save = X1
X1 = X1save
#VarList = c("CMD.total", "bFFP","Eref_sm","MCMT","PPT_JAS",
 #           "SHM","TD","PAS","DD5_sp","PPT06")
############New 16 Variable Set to Use 2019
########test value of Tmin02 versus MCMT and replace if better.
VarList = c("TD",	"Eref",	"AHM",	"Tmax07",	"CMD.total",	"Tmax_sp","Tmin_sm","CMD_sp","PAS",	
            "PPT_sp","PPT06","CMD.def","DD5_sp","Tmin_at","bFFP","MCMT")
#List = c("rn", "Region", "Latitude", "Longitude", "Elevation")
List = c("ID1", "ID2", "Latitude", "Longitude", "Elevation")
X1.sub = X1
X1.sub=X1[,names(X1) %in% c(List,VarList)]
X1 <- X1.sub
###########################
#X1 <- X1[!is.na(X1$Latitude),]
####modify variable names
Y1 <- X1
colnames(X1)[1]=c("PlotNo")
colnames(X1)[2]=c("BGC")
records <- nrow(X1)
X1$PlotNo <- as.character (X1$PlotNo)
attr(X1, "row.names") <- (X1$PlotNo)
X2 <- X1 [, c("PlotNo", "BGC", "Latitude", "Longitude", "Elevation")]
X1 <- X1[,-c(1,3:5)] #c("PlotNo","Latitude", "Longitude", "Elevation")]
####for the lHC training data from below.
X2 <- X1 [, c("ID1", "Zone", "Subzone", "Latitude", "Longitude", "Elevation")]
drop <- c("ID1", "Zone", "Latitude", "Longitude", "Elevation")
X1= X1 [, !(names(X1) %in% drop)]
colnames (X1) [1] = "BGC" 
X1save2 = X1
#X1<- as.data.frame(X1)
X1$BGC <- gsub("[[:space:]]","",X1$BGC) ## for subzone
############## Drop subzone for zone model
#X1$BGC <- gsub("[[:lower:]]|[[:space:]]|[[:digit:]]", "", X1$BGC) ###for zone model
X1 <- X1[X1$BGC != "",]
X1$BGC <- as.factor(X1$BGC)
#colnames(X1)[1]<- "Zone"
##########to be tested - training point balancing
X1.sub2 <- SmoteClassif(BGC ~ ., X1, C.perc = "balance", k= 5 , repl = TRUE, dist = "Euclidean")
prop.table(table(X1.sub2$BGC))
prop.table(table(X1$BGC))
###############Build RandomForest Model t#################
require(doParallel)
set.seed(123321)
coreNo <- makeCluster(detectCores() - 1)
registerDoParallel(coreNo, cores = detectCores() - 1)
Cores <- as.numeric(detectCores()-1)
clusterEvalQ(coreNo, .libPaths("E:/R packages351"))
##may need to change BGC for Zone where running for entire model of for subzone by zone foreach loop
ptm <- proc.time()
BGCmodel <- randomForest(BGC ~ ., data=X1, nodesize = 5, do.trace = 10,
                         ntree=71, na.action=na.fail, importance=TRUE, proximity=FALSE)
proc.time() - ptm

##############Parallel version
ptm <- proc.time()

Zonemodel <- foreach(ntree = rep(15,7), .combine = randomForest::combine, .multicombine = TRUE, .packages = "randomForest") %dopar% {
  randomForest(X1[-1], X1[1],
               ntree=ntree, na.action=na.fail, importance=TRUE, proximity=FALSE)
}
proc.time() - ptm
   #### Save random forest model
file=paste("USA_Subzones_RFModel_16Var_SMOTE",".Rdata",sep="")
save(BGCmodel,file=file)
load("USA_Subzones_RFModel_10VarLHC_v11_4_SMOTE.Rdata")
#load("USA_Subzones_RFModel_27Var.Rdata")
#load(file.choose())
#####Predict BGC units##########################
load("BGCv11_AB_USA_LHC_10VAR_SubZone_RFmodel.Rdata") ##load random forest model
### Read in US Hex Grid file with normal data. SOMETIMES THE FIRST 2 COLUMN NAMES NEED TO BE CHANGED
Columns = c("rn", "Region", "Latitude", "Longitude", "Elevation", "AHM", "bFFP",
            "CMD07","DD5_sp","EMT","Eref_sm","EXT","FFP","MCMT","MSP",
            "PPT07","PPT08", "PPT05","PPT06","PPT09", "SHM","TD","Tmax_sp","Tmin_at",
            "Tmin_sm","Tmin_wt", "PPT_at","PPT_wt", "PAS","eFFP",
            "Eref09","MAT","Tmin_sp","CMD")#"GCM", 

pointDat <- fread(file.choose(), stringsAsFactors = FALSE,  data.table = FALSE)#, select = Columns
#pointDat <- pointDat[pointDat$DD_0_01 != -9999,]
pointDat <- pointDat[pointDat$Latitude < 49,]
X1 <-pointDat
### add additional varaibles to match rF model
X1$PPT_MJ <- X1$PPT05 + X1$PPT06 # MaY/June precip
X1$PPT_JAS <- X1$PPT07 + X1$PPT08 + X1$PPT09 # July/Aug/Sept precip
X1$PPT.dormant <- X1$PPT_at + X1$PPT_wt # for calculating spring deficit
X1$CMD.def <- 500 - (X1$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
X1$CMD.def [X1$CMD.def < 0] <- 0 #negative values set to zero = no deficit
X1$CMDMax <- X1$CMD07
X1$CMD.total <- X1$CMD.def + X1$CMD
X1save = X1
X1 = X1save
colnames(X1)[1]=c("ID1")
colnames(X1)[2]=c("ID2")
#####New 16 Variable list to Use.
VarList = c("TD",	"Eref",	"AHM",	"Tmax07",	"CMD.total",	"Tmax_sp","Tmin_sm","CMD_sp","PAS",	
            "PPT_sp","PPT06","CMD.def","DD5_sp","Tmin_at","bFFP","MCMT")#16Var
#VarList = c("AHM", "bFFP","CMD.total","DD5_sp","EMT","Eref_sm","EXT","FFP","MCMT","MSP",
#            "PPT_JAS","PPT_MJ","PPT06","SHM","TD","Tmax_sp","Tmin_at","Tmin_sm","Tmin_wt",
#            "PAS","CMD.def","CMDMax","eFFP","Eref09","MAT","PPT07","Tmin_sp")
List = c("ID1","ID2", "Latitude", "Longitude", "Elevation")
X1.sub = X1
X1.sub=X1[,names(X1) %in% c(List,VarList)]
X1 <- X1.sub
pointDat <-X1.sub
#pointDat <- Y1
###Predict BGC units
pointDat$Zone <- predict(BGCmodel, newdata = pointDat[,-c(1:5)])
 pointDat.saved <- pointDat
 pointDat <-pointDat.saved
 #points <- pointDat[,c("Latitude","Longitude","Zone")]
 pointDat$BGC.pred <- as.character(pointDat$Zone)
 pointDat$Subzone <- gsub("[[:space:]]","",pointDat$BGC.pred)##for subzones
 pointDat$Zone <- gsub("[[:lower:]]|[[:space:]]|[[:digit:]]","",pointDat$BGC.pred)##for just zone
points <- pointDat[,c("ID1", "Zone", "Subzone", "Latitude","Longitude","Elevation")]
List2 = c("ID1", "Zone", "Subzone", "Latitude","Longitude","Elevation")
List2 = c("ID1", "Subzone")
pointLHC <- pointDat[,names(pointDat) %in% c(List2,VarList)]###includes all the climate variables for LHC procedure below
colnames (pointLHC) [18] <- "ID2"
col_order <-c("ID1", "Zone", "Subzone", "Latitude","Longitude","Elevation","CMD.total", "bFFP","Eref_sm","MCMT","PPT_JAS",
              "SHM","TD","PAS","DD5_sp","PPT06")
pointLHC <- pointLHC[,col_order]
#colnames(points)[2]<- "ZoneSubzone"
#colnames(points)[1]<- "ID1"
###Write to CSV to import into QGIS and join to HexGrid
write.csv (points, "USA_FAIPlots_BGCPredicted.csv")
write.csv (pointLHC, "USA_1k_HexGridPts_BGC_predicted_16Var_LHC_May3_w_BGCdata.csv")
countZone <- points %>% count(Zone)
countSubzone <- points %>% count(Subzone)
countSubzone$logn <- log(countSubzone$n, 10)
countSubzone$rs <- as.integer (rescale(countSubzone$logn, to = c(300, 1000), from = range(countSubzone$logn, na.rm = TRUE, finite = TRUE)))
countSubzone$sample <- ifelse(countSubzone$rs > countSubzone$n, countSubzone$n, countSubzone$rs )
allUnits <- unique(points$Subzone)
write.csv (countSubzone, "USA_ptPerBECLHC2.csv")
###########From here add CSV to QGIS and Join by ID1 to HexPolygon Shape file
  
##############select a subset of training points   
library(clhs)
BGC = "ICHxwz"
allUnits <- unique(pointLHC$Subzone)
pointLHC$ID1 <- row.names(pointLHC)
LHCtraining <- foreach(BGC = allUnits, .combine = rbind, .packages = c("clhs")) %dopar% {
temp <- pointLHC[(pointLHC$Subzone %in% BGC),]
Num <- countSubzone$sample[(countSubzone$Subzone %in% BGC)]
  samples <- clhs(temp[,-c(1:6)], 
                size = Num,           # Test a range of sample sizes
                iter = 100,        # Arbitrarily large number of iterations of optimization procedure. Default=10,000 but a larger number may be used
                progress = TRUE, 
                simple = FALSE)

cLHS <- samples$sampled_data
cLHS$ID1 <- row.names (cLHS)
cLHS_Points <- merge (pointLHC[,c(1:6)], cLHS, by = "ID1")
cLHS_Points
}

write.csv(LHCtraining, "USA_Training_LHC.csv")
X1 <- LHCtraining #feed these back into building a new rF model above

#Older Code
#############################################################################
###Colour schema
#USNames <- read.csv("US_Only_TrainingPoints_12July2018.csv")
USNames<- points[,3:4]
USNames <- unique(USNames)
#USNames$ZoneSubzone <- gsub("[[:space:]]","",USNames$ZoneSubzone)
zoneCols <- read.csv("USZone_Colours.csv", stringsAsFactors = FALSE)
zoneCols <- zoneCols[,-2]
US_BGC_Colours <- merge(USNames, zoneCols,  by = "Zone", all.x = TRUE) ##for zone
write.csv (US_BGC_Colours, "US_BGC_Colours.csv")
#points <- points[points$BGC.pred %in% USNames$ZoneSubzone,]

####Create maps
require(sp)
require(raster)
require(rgdal)
require(rgeos)
###Create Shape File of Hex grid points
CRS.albers <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
states <- readOGR(dsn = "BC_AB_US_Shp", layer = "USA_States")
states <- spTransform(states, CRS.albers)
US_BGC_Colours <- read.csv("US_BGC_Colours.csv", stringsAsFactors = FALSE)
US_BGC_Colours <- US_BGC_Colours[,-c(1:2)]
points2 <- merge(points, US_BGC_Colours, by.x = "Subzone", by.y = "Subzone", all.x = TRUE)
write.csv(points, "US_1KHexGrid_BGCpred.CSV")
points3 <- points2
points2 <- points3
points2$Latitude <- as.numeric(as.character(points2$Latitude))
points2$Longitude <- as.numeric(as.character(points2$Longitude))
#coordinates(points2) <- ~ Latitude + Longitude
coordinates(points2) <-  c("Longitude","Latitude")
proj4string(points2) <- CRS("+init=epsg:4617") ##NAD83
#proj4string(points2) <- CRS("+init=epsg:3005")
points2 <- spTransform(points2, CRS.albers)  # standard albers projection for BC gov't data
writeOGR(points2, dsn = "USAShapefiles", layer = "1kHexGridAllSubzonesPredictednew3", driver = "ESRI Shapefile") ###write as shapefile



#####Now create or load hexagon polygons to overlay with points
library(dplyr)
library(tidyr)
library(sp)
library(raster)
library(rgeos)
library(rgbif)
library(viridis)
library(gridExtra)
library(rasterVis)
# Read in Hexagon shape file and overlay
USA1kHex <- readOGR(dsn = "BC_AB_US_Shp", layer = "USA_1kmHexPoly_predictedBGC_30Jan2018")
USA1kHex <- spTransform(USA1kHex, CRS.albers)
USA1kHexBGC <-over(USA1kHex, points2)
size <- 1
#points4 <- as.SpatialPolygons.GridTopology(points2 , proj4string = CRS.albers)
#points4 <- points2[,c(2:3)]
#points4 <- as.matrix(points3)
#coordinates(points4) <- ~ Latitude + Longitude
points3 <- points
points3$Latitude <- as.numeric(points3$Latitude)
points3$Longitude <- as.numeric(points3$Longitude)
coordinates(points3) <- ~ Latitude + Longitude
size <- 100
pointtest <-spsample(US, type = "hexagonal", cellsize = size, offset = c(0,0))
US_hex_grid <- HexPoints2SpatialPolygons(pointtest)
#plot(study_area, col = "grey50", bg = "light blue", axes = TRUE)
#plot(hex_points, col = "black", pch = 20, cex = 0.5, add = T)
#plot(hex_grid, border = "orange", add = T)
US_hex_grid3 <- spTransform(US_hex_grid2, CRS.albers)  # standard albers projection for BC gov't data
writeOGR(US_hex_grid3, dsn = "USAShapefiles", layer = "1kUSAHexagons", driver = "ESRI Shapefile") ###write as shapefile

############################################################################################
#BC <- readOGR(dsn = "BC_AB_US_Shp", layer = "ProvincialOutline")###read in shapefiles
#AB <- readOGR(dsn = "BC_AB_US_Shp", layer = "AlbertaSubRegions")
US <- readOGR(dsn = "BC_AB_US_Shp", layer = "USA_States")
#BC <- spTransform(BC,CRS.albers)##project to albers
#AB <- spTransform(AB,CRS.albers)
US <- spTransform(US,CRS.albers)


#####Predict subzones within each zone########################
Zones <- as.character(unique(gsub("[[:lower:]]|[[:space:]]|[[:digit:]]", "", X1$BGC)))#X1$BGC <- gsub("[[:lower:]]|[[:space:]]|[[:digit:]]", "", X1$BGC) ###for zone model
sort(Zones)
X1 <- X1save
X1 <- X1[X1$BGC != "",]
Y1 <- X1
pointDat <- pointDat[,c(length(pointDat), 1:(length(pointDat)-1))]
colnames(Y1)[1] <- "BGC"
#Y1 <- Y1[,-1]
Z = list("ICH")
####this loop creates a RF model within each predicted zone for subzones, and output a map and dataset
SZPred <- foreach(Z = Zones, .combine = rbind) %do% {
  trainSub <- Y1[grep(paste(Z,"",sep = ""),Y1$BGC),]###subset training points to only include selected zone
  if(length(unique(trainSub$BGC)) >= 2){ ###somtimes there aren't any subzones
    trainSub$BGC <- as.factor(trainSub$BGC)
    
    SZmodel <- randomForest(BGC ~ ., data=trainSub, nodesize = 5, do.trace = 10, 
                            ntree=101, na.action=na.fail, importance=TRUE, proximity=FALSE)
    
    pointDat <- pointDat.saved
    pointSub <- pointDat[pointDat$Zone == Z,] ###subset grid based on previously predicted zones
    pointSub$Subzone <- predict(SZmodel, newdata = pointSub[,-c(1:5)]) ### predict subzones
    
    points <- pointSub[,c("Latitude","Longitude","Subzone")]
    colnames(points)[3]<- "BGC.pred"
    points <- points[points$Latitude < 49,]
    coordinates(points) <- c("Longitude","Latitude")
    proj4string(points) <- CRS("+init=epsg:4326")
    points <- spTransform(points, CRS.albers)  # standard albers projection for BC gov't data
    
    ###create map
    plot_pal <- rainbow(length(unique(points$BGC.pred)))
    plot_col <- plot_pal[points$BGC.pred]
    
    pdf(paste("USA_",Z,".pdf",sep = ""),paper = "letter",width = 8)
    plot(states, main = Z)
    plot(points, col = plot_col, pch = 15, cex = 0.08, add = TRUE)
    
    #pdf(paste("USA_",Z,".pdf",sep = ""),paper = "letter",width = 8)
    #plot(states, main = unique(trainSub$BGC))
    #plot(points, col = "purple", pch = 15, cex = 0.08, add = TRUE)
    #dev.off()
    
    plotExtent <- extent(states)
    xEx <- plotExtent@xmax - 100000
    yEx <- plotExtent@ymax + 100000
    legend(x = xEx, y = yEx,
           legend = levels(points$BGC.pred),
           fill = plot_pal,
           xpd = TRUE,
           bty = "n",
           cex = 0.8)
    dev.off()
    
    out <- pointSub[,c("Latitude","Longitude","Zone", "Subzone")]
    out
  }else{ ##if only one subzone, plot 
    pointSub <- pointDat[pointDat$Zone == Z,]
    pointSub$Subzone <- pointSub$Zone
    points <- pointSub[,c("Latitude","Longitude","Zone","Subzone")]
    colnames(points)[3]<- "BGC.pred"
    points <- points[points$Latitude < 49,]
    out <- data.frame(Latitude = points$Latitude, Longitude = points$Longitude, Zone = Z, Subzone = Z)
    coordinates(points) <- c("Longitude","Latitude")
    proj4string(points) <- CRS("+init=epsg:4326")
    points <- spTransform(points, CRS.albers)  # standard albers projection for BC gov't data
    
    pdf(paste("USA_",Z,".pdf",sep = ""),paper = "letter",width = 8)
    plot(states, main = unique(trainSub$BGC))
    plot(points, col = "purple", pch = 15, cex = 0.08, add = TRUE)
    dev.off()
    out <- pointSub[,c("Latitude","Longitude","Zone", "Subzone")]
    out
  }
  
}
##########################
####OUtput a complete shape file by running all
USNames <- as.data.table (USNames)
sbzNo <- USNames[, (Number = uniqueN(ZoneSubzone)), by = USNames$Zone]
colnames(sbzNo) [1:2] <- c("Zone", "Number")
withSBZ <- sbzNo[sbzNo$Number >1,]
noSBZ <- sbzNo[sbzNo$Number == 1,]
X1 <- X1save
X1 <- X1[X1$BGC != "",]
Y2 <- X1
pointDat <- pointDat[,c(length(pointDat), 1:(length(pointDat)-1))]
colnames(Y2)[1] <- "BGC"
#### need line to add Zone to Y2
Y2$withSBZ <- (gsub("[[:lower:]]|[[:space:]]|[[:digit:]]", "", Y2$BGC))###limit dataset to only those zones with multiple subzones
SZPred <- foreach(Z = withSBZ$Zone, .combine = rbind) %do% {
  trainSub <- Y2[Y2$withSBZ == Z,]###subset training points to only include selected zone
  trainSub$BGC <- as.factor(trainSub$BGC)
    SZmodel <- randomForest(BGC ~ ., data=trainSub[-249], nodesize = 5, do.trace = 10, 
                            ntree=101, na.action=na.fail, importance=TRUE, proximity=FALSE)
    
    pointDat <- pointDat.saved
    pointSub <- pointDat[pointDat$Zone == Z,] ###subset grid based on previously predicted zones
    pointSub$Subzone <- predict(SZmodel, newdata = pointSub[,-c(1:5)]) ### predict subzones
    
    points <- pointSub[,c("Latitude","Longitude","Subzone")]
    colnames(points)[3]<- "BGC.pred"
   out <- pointSub[,c("Latitude","Longitude","Zone", "Subzone")]
    out
}
    
  ###need to merge back in the single subzone Zone Points  
    pointDat.nosub <- pointDat.saved[pointDat.saved$Zone %in% noSBZ$Zone,]
    pointDat.nosub$Subzone <- pointDat.nosub$Zone
    points.nosub <- pointDat.nosub[,c("Latitude","Longitude","Zone","Subzone")]
    all.grid <- rbind(SZPred, points.nosub)

# Then convert entire new point set to shape file
coordinates(all.grid) <- c("Longitude","Latitude")
proj4string(all.grid) <- CRS("+init=epsg:4326")
all.grid2 <- spTransform(all.grid, CRS.albers)
writeOGR(all.grid2, dsn = "USAShapefiles", layer = "AllSubzones", driver = "ESRI Shapefile") ###write as shapefile

rnd.grid <- stratified.random(all.grid2, strata = "Subzone", n = 500, reps = 1, replace = FALSE)
 writeOGR(rnd.grid, dsn = "USAShapefiles", layer = "Training_f_Grid", driver = "ESRI Shapefile") ###write as shapefile
 
#################NEW TEST SCRIPT FOR GENERATING HEXAGONAL 
#FUNCTIONS
  make_grid <- function(x, cell_diameter, cell_area, clip = FALSE) {
   if (missing(cell_diameter)) {
     if (missing(cell_area)) {
       stop("Must provide cell_diameter or cell_area")
     } else {
       cell_diameter <- sqrt(2 * cell_area / sqrt(3))
     }
   }
   ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
   projection(ext) <- projection(x)
   # generate array of hexagon centers
   g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter, 
                 offset = c(0.5, 0.5))
   # convert center points to hexagons
   g <- HexPoints2SpatialPolygons(g, dx = cell_diameter)
   # clip to boundary of study area
   if (clip) {
     g <- gIntersection(g, x, byid = TRUE)
   } else {
     g <- g[x, ]
   }
   # clean up feature IDs
   row.names(g) <- as.character(1:length(g))
   return(g)
 }
 
#sCRIPT 
 study_area_utm <- CRS("+proj=utm +zone=44 +datum=WGS84 +units=km +no_defs") %>% 
   spTransform(study_area, .)
 # without clipping
 hex_grid <- make_grid(study_area_utm, cell_area = 625, clip = FALSE)
 plot(study_area_utm, col = "grey50", bg = "light blue", axes = FALSE)
 plot(hex_grid, border = "orange", add = TRUE)
 box()
 # with clipping
 hex_grid <- make_grid(study_area_utm, cell_area = 625, clip = TRUE)
 plot(study_area_utm, col = "grey50", bg = "light blue", axes = FALSE)
 plot(hex_grid, border = "orange", add = TRUE)
 box()
