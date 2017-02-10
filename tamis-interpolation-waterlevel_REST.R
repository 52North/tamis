# wps.des: tamis-rest-interpolation-wasserstand, title = Interpolation of Wasserstand at Bevertalsperre;

# wps.in: timeseries, string, set of TS URIs, whitespace " " seperated, 
# abstract = timeseries URI as data source,
# value = "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/514 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/515 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/470 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/473 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/474 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/476 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/479 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/482";

# wps.in: timespan, string, timespan of input data, 
# abstract = timeseries URI for the interpolation variable,
# value = "2016-01-01T/2016-09-30TZ";

# wps.in: target, type = geotiff, 
# abstract = geotiff defining the interpolation grid (only non NAs will be interpolated),
# value = "https://github.com/BenGraeler/tamis/raw/master/geotiff.tiff";

## tamis

# updateStatus("Loading libraries")

library(httr)
library(rjson)
library(lattice)
library(sp)
library(spacetime)
library(gstat)
library(rgdal)
library(RCurl)

# "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/444 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/491 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/492 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/435 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/436 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/437 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/438 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/439 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/440 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/441 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/442";

# updateStatus("helper functions")

readTSdata <- function(ts_URI, timespan, .opts, ...) {
  if(!missing(.opts))
    meta <- GET(ts_URI, do.call(config, .opts), ...)
  else
    meta <- GET(ts_URI)
  
  meta <- memDecompress(meta$content, "none", asChar = T)
  meta <- substr(meta, gregexpr("\"id", meta)[[1]][1]-1, nchar(meta))
  meta <- fromJSON(meta)
  
  if(!missing(.opts))
    ts <- GET(paste(ts_URI, "/getData?timespan=", timespan, sep=""), do.call(config, .opts), ...)
  else 
    ts <- GET(paste(ts_URI, "/getData?timespan=", timespan, sep=""))
  
  ts <- memDecompress(ts$content, "none", asChar = T)
  ts <- substr(ts, gregexpr("\"values", ts)[[1]][1]-1, nchar(ts))
  ts <- fromJSON(ts)
  
  ts <- do.call(rbind, ts$values)
  ts <- data.frame(time = as.POSIXct(as.numeric(ts[,1])/1e3, origin="1970-01-01"),
                   val = unlist(ts[,2]))
  
  if(!is.null(meta$parameters$phenomenon$label))
    colnames(ts)[2] <- meta$parameters$phenomenon$label  
  
  coords <- matrix(as.numeric(unlist(meta$station$geometry$coordinates)), nrow=1)
  coords <- coords[,!is.nan(coords) & !is.na(coords), drop=F]
  
  if (length(ts[,1]) == 1)
    return(STFDF(SpatialPoints(coords), ts[,1], ts[,-1,drop=F], ts[,1]))
  STFDF(SpatialPoints(coords), ts[,1], ts[,-1,drop=F])
}

readTSmeta <- function(ts_URI, .opts, ...) {
  if(!missing(.opts))
    meta <- GET(ts_URI, do.call(config, .opts), ...)
  else
    meta <- GET(ts_URI, ...)
  meta <- memDecompress(meta$content, "none", asChar = T)
  meta <- substr(meta, gregexpr("\"id", meta)[[1]][1]-1, nchar(meta))
  fromJSON(meta)
}

###

# updateStatus("Loading secOpts")

source("~/52North/secOpts.R")

# wps.off;

timeseries <- "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/514 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/515 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/470 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/473 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/474 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/476 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/479 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/482"
timespan <-  "2016-08-06T13:44:33.471Z/2017-02-06T14:44:33.471Z"
target <- "geotiff.tiff" 

# wps.on;

isGrid <- FALSE

fileType <- tail(strsplit(target,split =  ".", fixed = T)[[1]],1)

if (fileType == "tiff" | fileType == "tif") {
  isGrid <-TRUE
  target <- readGDAL(target)
  target <- as(target,"SpatialPointsDataFrame")
} else {
  library(rgeos)
  target <- readWKT(target)
}

timeseries <- strsplit(timeseries, split = " ", fixed = T)[[1]]
timeseries <- timeseries[nchar(timeseries) > 0]

# updateStatus("loading TS data")

heightCorrection <- TRUE

dataObs_STFDF <- NULL
for (ts in timeseries) { # ts <- timeseries[3]
  source <- as(readTSdata(ts, timespan, .opts), "data.frame")
  tsId <- tail(strsplit(ts, "/")[[1]],1)
  if (heightCorrection)
    source$Wasserstand.im.Damm <- switch(tsId,
                                         "513"=297.46, # MQA1: 
                                         "514"=298.32, # MQA3: 
                                         "515"=287.07, # MQA4: 
                                         "516"=286.75, # MQA5: 
                                         "470"=275.67, # MQA7: 
                                         "472"=296.71, # MQB1: 
                                         "473"=298.43, # MQB2: 
                                         "474"=287.11, # MQB4: 
                                         "475"=287.21, # MQB5: 
                                         "476"=275.96, # MQB6:
                                         "477"=276.01, # MQB7: 
                                         "478"=296.26, # MQC1: 
                                         "479"=298.36, # MQC2: 
                                         "480"=287.07, # MQC4: 
                                         "481"=287.06, # MQC5: 
                                         "482"=275.96, # MQC6
                                         "483"=275.99) - source$Wasserstand.im.Damm # MQC7
  source$distanceToDam <- switch(tsId,
                                 "513"=-2.66, # MQA1: 
                                 "514"=7.17,  # MQA3: 
                                 "515"=28.47, # MQA4: 
                                 "516"=28.47, # MQA5: 
                                 "470"=52.74, # MQA7: 
                                 "472"=-5.04, # MQB1: 
                                 "473"=7.10,  # MQB2: 
                                 "474"=28.41, # MQB4: 
                                 "475"=28.10, # MQB5: 
                                 "476"=51.90, # MQB6:
                                 "477"=51.96, # MQB7: 
                                 "478"=-5.17, # MQC1: 
                                 "479"=6.98,  # MQC2: 
                                 "480"=28.39, # MQC4: 
                                 "481"=28.48, # MQC5: 
                                 "482"=52.36, # MQC6
                                 "483"=52.19) # MQC7
                                       
  
  dataObs_STFDF <- rbind(dataObs_STFDF, source)
}

dataObs_STFDF <- stConstruct(dataObs_STFDF[,c(1:2,4,7:8)], space = c("coords.x1", "coords.x2"), time = "time")
dataObs_STFDF <- as(dataObs_STFDF, "STFDF")
colnames(dataObs_STFDF@data) <- c("Wasserstand","distanceToDam")

tsIds <- sapply(timeseries, function(x) tail(strsplit(x, "/")[[1]],1))
tsNames <- sapply(tsIds, function(x) switch(x,
                  "513"="MQA1",
                  "514"="MQA3",
                  "515"="MQA4",
                  "516"="MQA5",
                  "470"="MQA7",
                  "472"="MQB1",
                  "473"="MQB2",
                  "474"="MQB4",
                  "475"="MQB5",
                  "476"="MQB6",
                  "477"="MQB7",
                  "478"="MQC1",
                  "479"="MQC2",
                  "480"="MQC4",
                  "481"="MQC5",
                  "482"="MQC6",
                  "483"="MQC7"))
rownames(dataObs_STFDF@sp@coords) <- tsNames

# updateStatus("Setting CRS")

dataObs_STFDF@sp@proj4string <- CRS("+init=epsg:4326")

# updateStatus("Requested Data successfully")

n.time <- length(dataObs_STFDF@time)

dataPos <- "dataPos.png"
png(file = dataPos)
tmpDataPos <- stplot(dataObs_STFDF[,1:min(12, n.time)])
print(tmpDataPos)
graphics.off()

# wps.out: dataPos, png;

dataTS <- "dataTS.png"
png(file = dataTS)
tmpDataTS <- stplot(dataObs_STFDF[c("MQB4", "MQB5", "MQB6", "MQB7"),, "Wasserstand", drop=F], mode="ts")
print(tmpDataTS)
graphics.off()

# wps.out: dataTS, png;

# updateStatus("Calculate empirical variogram")

dataObs_STFDF@sp <- spTransform(dataObs_STFDF@sp, target@proj4string)
colnames(dataObs_STFDF@sp@coords) <- c("x","y")

dataObs_STSDF <- as(dataObs_STFDF, "STSDF")

linMod <- lm(Wasserstand ~ x+y, as(dataObs_STSDF[,,drop=F], "data.frame"))
summary(linMod)

dataObs_STSDF@data$resid <- linMod$residuals

empVgm <- variogram(resid ~ 1, dataObs_STSDF, tlags=0, na.omit = TRUE, cutoff=50)
empVgm <- cbind(empVgm, data.frame(dir.hor=rep(0,nrow(empVgm)), dir.ver=rep(0,nrow(empVgm))))
class(empVgm) <- c("gstatVariogram","data.frame")

fitVgm <- fit.variogram(empVgm, vgm(median(empVgm$gamma)*0.5, 
                                    "Lin", 
                                    mean(empVgm$dist, na.rm = T),
                                    median(empVgm$gamma)*0.5))


# updateStatus("Theoretical variogram has been fitted.")

tmpPlot <- plot(empVgm, fitVgm)

# wps output
vgmFit <- "vgmFit.png"
png(file = vgmFit)
print(tmpPlot)
graphics.off()
# wps.out: vgmFit, png;

# updateStatus("Theoretical variogram has been plotted.")

###
targetData <- NULL
targetVar <- NULL

for (day in 1:n.time) { # day <- 1
  pred <- krige(resid ~ 1, dataObs_STSDF[,day], target, model=fitVgm)@data
  targetData <- cbind(targetData, pred$var1.pred)
  targetVar <- cbind(targetVar, pred$var1.var)
}

target <- addAttrToGeom(geometry(target), as.data.frame(cbind(targetData, targetVar)))
gridded(target) <- TRUE

linModPred <- predict(linMod, as.data.frame(target@coords))
target@data <- cbind(target@data, linModPred)

for (i in 1:((ncol(target@data)-1)/2)) { # i <- 1
  target@data[,i] <- linModPred + target@data[,i]
}

# updateStatus("Interpolation of residuls has been performed.")

# n.plots <- min(12,n.time)
# 
# p1 <- spplot(target[,1:n.plots],
#              sp.layout=list("sp.points",
#                             dataObs_STSDF@sp,
#                             col="black"),
#              strip=strip.custom(factor.levels=as.character(as.Date(index(dataObs_STSDF@time)))[1:n.plots]),
#              as.table=T)

# png("IDW_interpolation_series_waterlevel.png", width = 2400, height = 1600, res=200)
# print(p1)
# dev.off()

if(isGrid) {
  predictions <- "predictions.tiff"
  predVar <- "predVar.tiff"
} else { 
  predictions <- "predictions.csv"
}

if (isGrid) {
  gridded(target) <- TRUE
  writeGDAL(target[,1:n.time], predictions, drivername="GTiff")
  writeGDAL(target[,1:n.time + n.time], predVar, drivername="GTiff")
} else {
  write.csv(as.data.frame(target[,1:n.time]), predictions)
  write.csv(as.data.frame(target[,1:n.time + n.time]), predVar)
}

# wps.out: predictions, type = geotiff;
# wps.out: predVar, type = geotiff;

# if a grid is presented, a NetCDF file can be returned

if(isGrid) {
  gridded(target) <- TRUE
  target_fullGrid <- target
  fullgrid(target_fullGrid) <- TRUE
  
  target_STFDF <- STFDF(geometry(target_fullGrid),
                        dataObs_STFDF@time, 
                        data.frame(var1.pred = unlist(target_fullGrid@data[, 1:n.time]),
                                   var1.var = unlist(target_fullGrid@data[, 1:n.time + n.time])),
                        if(length(index(dataObs_STFDF@time)) > 1) {
                          delta(dataObs_STFDF@time)
                        } else {
                          index(dataObs_STFDF@time)
                        })

  # wps output
  predMap <- "predMap.png"
  png(file = predMap)
  tmpPlot <- stplot(target_STFDF[,sort(sample(n.time, min(n.time,12)))],
                    sp.layout=list("sp.points", dataObs_STFDF@sp))
  print(tmpPlot)
  graphics.off()
  # wps.out: predMap, png;
  
  ## netCDF
  library(RNetCDF)
  ncFile <- "targetnetcdf.nc"
  
  nc <- create.nc(ncFile, prefill = FALSE)
  
  # x
  dim.def.nc(nc, "x", target_STFDF@sp@grid@cells.dim[1])
  var.def.nc(nc, "x", "NC_DOUBLE", "x")
  var.put.nc(nc, "x", unique(target_STFDF@sp@coords[,"x"]))
  
  att.put.nc(nc, "x", "units", "NC_CHAR", "meter")
  att.put.nc(nc, "x", "axis", "NC_CHAR", "x")
  att.put.nc(nc, "x", "long_name", "NC_CHAR", "x")
  att.put.nc(nc, "x", "standard_name", "NC_CHAR", "projection_x_coordinate")
  
  # y
  dim.def.nc(nc, "y", target_STFDF@sp@grid@cells.dim[2])
  var.def.nc(nc, "y", "NC_DOUBLE", "y")
  var.put.nc(nc, "y", unique(target_STFDF@sp@coords[,"y"]))
  
  att.put.nc(nc, "y", "units", "NC_CHAR", "meter")
  att.put.nc(nc, "y", "axis", "NC_CHAR", "y")
  att.put.nc(nc, "y", "long_name", "NC_CHAR", "y")
  att.put.nc(nc, "y", "standard_name", "NC_CHAR", "projection_y_coordinate")
  
  # time
  dim.def.nc(nc, "time", length(target_STFDF@time))
  var.def.nc(nc, "time", "NC_INT", "time")
  var.put.nc(nc, "time", round(as.numeric(index(target_STFDF@time))/60,0))
  
  att.put.nc(nc, "time", "units", "NC_CHAR", "minutes since 1970-01-01")
  att.put.nc(nc, "time", "axis", "NC_CHAR", "t")
  att.put.nc(nc, "time", "calendar", "NC_CHAR", "gregorian")
  att.put.nc(nc, "time", "long_name", "NC_CHAR", "time")
  att.put.nc(nc, "time", "standard_name", "NC_CHAR", "time")
  
  # CRS
  var.def.nc(nc, "crs", "NC_INT", NA)
  var.put.nc(nc, "crs", 31466)
  att.put.nc(nc, "crs", "EPSG_code", "NC_CHAR", "EPSG:31466")
  att.put.nc(nc, "crs", "proj4_params", "NC_CHAR", "+proj=tmerc +lat_0=0 +lon_0=6 +k=1 +x_0=2500000 +y_0=0 +ellps=bessel +units=m +no_defs")
  
  obsProp <- colnames(dataObs_STFDF@data)[1]
  
  var.def.nc(nc, obsProp, "NC_DOUBLE", NA)
  att.put.nc(nc, obsProp, "ancillary_variables", "NC_CHAR", "var1pred var1var")
  att.put.nc(nc, obsProp, "ref", "NC_CHAR", "http://www.uncertml.org/distributions/normal")
  att.put.nc(nc, obsProp, "shape", "NC_CHAR", "x y t")
  
  UncertML <- list(var1.pred = "http://www.uncertml.org/distributions/normal#mean",
                   var1.var = "http://www.uncertml.org/distributions/normal#variance")

  for (var.name in colnames(target_STFDF@data)) {
    varname <- gsub(".","", var.name, fixed=TRUE)
    var.def.nc(nc, varname, "NC_DOUBLE", dimensions = c(2,1,0))
    att.put.nc(nc, varname, "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc, varname, "ref", "NC_CHAR", UncertML[[var.name]])
    att.put.nc(nc, varname, "grid_mapping", "NC_CHAR", "crs")
    dArray <- array(target_STFDF@data[[var.name]], c(target_STFDF@sp@grid@cells.dim[1],
                                                     target_STFDF@sp@grid@cells.dim[2],
                                                     length(target_STFDF@time)))
    dArray <- aperm(dArray, c(3,2,1))
    var.put.nc(nc, varname, dArray, na.mode = 0)
  }
  
  att.put.nc(nc, "NC_GLOBAL", "Conventions", "NC_CHAR", "CF-1.5 UW-1.0")
  att.put.nc(nc, "NC_GLOBAL", "primary_variables", "NC_CHAR", obsProp)
  
  close.nc(nc)
}
# wps.out: ncFile, netcdf;
