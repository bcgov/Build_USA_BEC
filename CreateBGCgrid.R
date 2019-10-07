###This script imports BEC shapefiles and a DEM to create training point data. It also creates a 4km grid
###for producing maps
###Kiri Daust, August 2018

.libPaths("E:/R packages")

library(dplyr)
library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(maptools)
library(magrittr)
library(tibble)
library(tidyr)
library(sf)
library(tcltk)
library(foreach)
library(httr)
library(jsonlite)
library(randomForest)
library(data.table)
library(mapview)

wd <- tk_choose.dir(); setwd(wd)

dem <- raster("E:/Kiri Clean/CCISS Stuff/Data/bc25fill") ###Read DEM
bec11 <- st_read(dsn="bgc.v11.gdb",layer="bgcv11_bc") ##read bec file
bec11 <- bec11[c("ZONE","MAP_LABEL")]
becZone <- bec11 %>%
  group_by(ZONE) %>%
  summarise(geometry = sf::st_union(SHAPE)) %>%
  ungroup()


CRS.albers <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
allUnits <- unique(as.character(bec11$MAP_LABEL))###What units are in BEC?

temp <- bec11[bec11$MAP_LABEL == "IDFdk3",]

t1 <- st_cast(temp, "MULTIPOLYGON")
st_write(t1, dsn = "ForColin",layer = "IDFdk3_BEC11", driver = "ESRI Shapefile")
t1 <- t1 %>%
  group_by(ZONE) %>%
  summarise(geometry = sf::st_union(SHAPE)) %>%
  ungroup()

##set up for parallel processing
require(doParallel)
set.seed(123321)
coreNum <- as.numeric(detectCores()-1)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)
##clusterEvalQ(coreNo, .libPaths("E:/R packages"))

fireShp <- st_read(dsn = "C:/Users/Kiri Daust/Desktop/PortfolioKiri/FiresIDFdk3", layer = "BCWildfires")
fireShp <- fireShp[fireShp$name %in% c("K20637","C50647","C20729"),c("name","geometry")]
fireShp <- st_transform(fireShp, crs = 3005)
firePnts <- st_sample(fireShp, size = c(80,80,80), type = "random", exact = TRUE)
temp <- data.frame(Fire = rep(c("Ele","Plat","WillLk"), each = 80))
temp <- cbind(temp, st_coordinates(firePnts))
temp <- st_as_sf(temp, coords = c("X","Y")) %>% st_set_crs(3005)
bec11 <- st_transform(bec11, 3005)
firePnts <- st_join(temp, bec11, join = st_intersects)
firePnts <- st_transform(firePnts, st_crs(dem))
firePnts$Elevation <- raster::extract(dem, firePnts)
firePnts <- st_transform(firePnts, 4326) %>% cbind(.,st_coordinates(.)) %>% st_drop_geometry()
temp <- rep(1:80,3)
firePnts$Fire <- paste(firePnts$Fire,temp, sep = "_")
firePnts <- firePnts[,c("ID1","ID2","lat","long","el")]
colnames(firePnts) <- c("ID1","ID2","long","lat","el")
write.csv(firePnts, "C:/Users/Kiri Daust/Desktop/PortfolioKiri/FirePoints.csv", row.names = F)

####good loop
out <- foreach(Shp = fireShp$name, .combine = rbind, .packages = c("sf","sp","raster")) %dopar% {
  temp <-  ###Extract polygons for each subzones
  p <- st_sample(temp, size = 250, type = "regular")
  p <- st_transform(p, "+init=epsg:4326") #to lat long
  coords <- st_coordinates(p)
  coords <- as.data.frame(coords)
  
  p2 <- st_transform(p, st_crs(dem))
  p2 <- st_sf(p2)
  coords$Elevation <- raster::extract(dem,p2) ##get elevation for each point
  coords$BGC <- BGC
  coords
}

###randomly select 2000 points within each BEC unit and get elevation data
out <- foreach(BGC = allUnits, .combine = rbind, .packages = c("sf","sp","raster")) %dopar% {
  temp <- bec11$SHAPE[bec11$MAP_LABEL == BGC] ###Extract polygons for each subzones
  p <- st_sample(temp, size = 250, type = "regular")
  p <- st_transform(p, "+init=epsg:4326") #to lat long
  coords <- st_coordinates(p)
  coords <- as.data.frame(coords)
  
  p2 <- st_transform(p, st_crs(dem))
  p2 <- st_sf(p2)
  coords$Elevation <- raster::extract(dem,p2) ##get elevation for each point
  coords$BGC <- BGC
  coords
}

out$ID1 <- seq_along(out$Elevation)
out <- out[,c("ID1", "BGC","Y","X","Elevation")]
colnames(out) <- c("ID1","ID2","lat","long", "el")

write.csv(out,"BECv11_250Pt.csv", row.names = FALSE)


####4 km grid BC, AB, US########################
CRS.albers <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")

BC <- readOGR(dsn = "BC_AB_US_Shp", layer = "ProvincialOutline")
AB <- readOGR(dsn = "BC_AB_US_Shp", layer = "AlbertaSubRegions")
US <- readOGR(dsn = "BC_AB_US_Shp", layer = "USA_States")

BC <- spTransform(BC, CRS.albers)
AB <- spTransform(AB, CRS.albers)
US <- spTransform(US, CRS.albers)

BCdem <- raster("bc25fill")

###combine alberta and US DEMs
####Don't need to do this
ABdem <- raster("AlbertaDEM.tif")
USdem1 <- raster("USA_WestDEM.tif")
USdem2 <- raster("USA_WestCoastDEM.tif")
USdem <- raster::merge(USdem1,USdem2)
allDEM <- raster::merge(ABdem,USdem)
writeRaster(allDEM, filename = "US_AB_DEM.tif")
####################################################33

shp <- st_read(dsn = "BC_AB_US_Shp", layer = "ProvincialOutline")
p <- st_make_grid(shp, )

allDEM <- raster("US_AB_DEM.tif")

#####Loop through BC, Alberta, and US to create grid
names <- c("BC","AB","USA")
i <-0
grid4k <- foreach(X = c("ProvincialOutline", "AlbertaSubRegions", "USA_States"), .combine = rbind) %do% {
  i <- i+1
  shp <- readOGR(dsn = "BC_AB_US_Shp", layer = X)
  shp <- spTransform(shp, CRS.albers)
  p <- spsample(shp, cellsize = c(2000), type = "hexagonal")
  p <- spTransform(p, CRS("+init=epsg:4326")) #to lat long
  coords <- p@coords
  coords <- as.data.frame(coords)
  if(i == 1){
    dem <- BCdem
  }else{
    dem <- allDEM
  }
  p <- spTransform(p, CRS(proj4string(dem)))
  coords$elev <- raster::extract(dem,p)
  colnames(coords) <- c("long","lat","el")
  coords$Region <- names[i]
  coords
}

rownames(grid4k) <- NULL
gridBC <- grid4k[grid4k$Region == "BC",]##subset just for BC
#############################################################################3
####Assign BGCs and subregions
##################################################################################3
###for BC (input 4 km grid, output is just BC points with BGC assigned)

BCwBGC <- foreach(BGC = allUnits, .combine = rbind, .packages = c("sp","sf","raster")) %dopar%{
  dat <- gridBC
  pointsOrig <- dat
  coordinates(dat) <- c("lat","long")
  proj4string(dat) <- CRS("+init=epsg:4326")
  dat <- spTransform(dat, CRS.albers)  # standard albers projection for BC gov't data
  
  tempPoly <- bec11[bec11$MAP_LABEL == BGC,]
  tempPoly <- as(tempPoly, "Spatial") ##conver to sp
  tempPoly <- spTransform(tempPoly, CRS.albers) 
  dat <- over(dat, tempPoly) ###which ones are inside the BGC
  pointsOrig <- pointsOrig[!is.na(dat$BGC_LABEL),] ###Remove points not inside BGC
  if(nrow(pointsOrig) > 0){ ###check that some points fall inside BGC
    pointsOrig$BGC <- BGC
    pointsOrig
  }
}

###for AB (input 4 km grid, output is just Alberta points with BGC assigned)
abUnits <- as.character(AB$NSRNAME)

ABwBGC <- foreach(BGC = abUnits, .combine = rbind) %do%{
  dat <- grid4k
  pointsOrig <- dat
  coordinates(dat) <- c("lat","long")
  proj4string(dat) <- CRS("+init=epsg:4326")
  dat <- spTransform(dat, CRS.albers)  # standard albers projection for BC gov't data
  
  tempPoly <- AB[AB$NSRNAME == BGC,]
  tempPoly <- spTransform(tempPoly, CRS.albers) 
  dat <- over(dat, tempPoly) ###which ones are inside the BGC
  pointsOrig <- pointsOrig[!is.na(dat),] ###Remove points not inside BGC
  pointsOrig$BGC <- BGC
  pointsOrig
}

##########################
A1 <- readShapePoly("A1_bndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
A1 <- spTransform(A1, CRS.albers)
A2 <- readShapePoly("A2_BlkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
A2 <- spTransform(A2, CRS.albers)
A3 <- readShapePoly("A3_blkbndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
A3 <- spTransform(A3, CRS.albers)
A4 <- readShapePoly("A4_disply_bndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
A4 <- spTransform(A4, CRS.albers)
B1 <- readShapePoly("B1_BlkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
B1 <- spTransform(B1, CRS.albers)
B2 <- readShapePoly("B2_BlkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
B2 <- spTransform(B2, CRS.albers)
B3 <- readShapePoly("B3_BlkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
B3 <- spTransform(B3, CRS.albers)
B4 <- readShapePoly("B4_blkBndy.shp",  proj4string = CRS("+proj=longlat +ellps=WGS84"))
B4 <- spTransform(B4, CRS.albers)
B5 <- readShapePoly("B5_blkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
B5 <- spTransform(B5, CRS.albers)
C1 <- readShapePoly("C1_Blkbndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
C1 <- spTransform(C1, CRS.albers)
C2 <- readShapePoly("C2_blkbndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
C2 <- spTransform(C2, CRS.albers)
C3 <- readShapePoly("C3_blkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
C3 <- spTransform(C3, CRS.albers)
D2 <- readShapePoly("D2_blkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
D2 <- spTransform(D2, CRS.albers)
D3 <- readShapePoly("D3_BlkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
D3 <- spTransform(D3, CRS.albers)
D4 <- readShapePoly("D4_BlkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
D4 <- spTransform(D4, CRS.albers)
D5 <- readShapePoly("D5_blkBndry.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
D5 <- spTransform(D5, CRS.albers)

plotList <- c("D5","D4","D3","D2","C3","C2","C1","B5","B4","B3","B2","B1","A4","A3","A2","A1")

out <- foreach(Plot = plotList, .combine = rbind) %do% {
    elev <- raster::extract(dem, get(Plot))
    max <- max(elev[[1]])
    min <- min(elev[[1]])
    ave <- mean(elev[[1]])
    df <- data.frame(PlotID = Plot, Max = max, Min = min, Mean = ave)
    df
}
elev <- raster::extract(dem,A1)
