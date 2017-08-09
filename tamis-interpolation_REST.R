# wps.des: tamis-rest-interpolation-schuettmenge, title = Interpolation of Schuettmenge at Bevertalsperre;

# wps.in: timeseries, string, set of TS URIs, whitespace " " seperated, 
# abstract = timeseries URI as data source,
# value = "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/450 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/451 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/452 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/453 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/454 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/455";

# wps.in: timespan, string, timespan of input data, 
# abstract = timeseries URI for the interpolation variable,
# value = "2016-01-01T/2016-01-14TZ";

# wps.in: target, type = geotiff, 
# abstract = geotiff defining the interpolation grid (only non NAs will be interpolated),
# value = "https://github.com/BenGraeler/tamis/raw/master/geotiff.tiff";

## tamis

# updateStatus("Loading libraries")

library(httr)
library(rjson)
library(sp)
library(spacetime)
library(gstat)
library(rgdal)
library(RCurl)

# "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/464 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/465 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/466 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/467 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/468 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/469";

# updateStatus("helper functions")

readTSdata <- function(ts_URI, timespan, .opts, ...) {
  if(!missing(.opts) & !is.null(.opts))
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
  if(!missing(.opts) & !is.null(.opts))
    meta <- GET(ts_URI, do.call(config, .opts), ...)
  else
    meta <- GET(ts_URI, ...)
  meta <- memDecompress(meta$content, "none", asChar = T)
  meta <- substr(meta, gregexpr("\"id", meta)[[1]][1]-1, nchar(meta))
  fromJSON(meta)
}

checkCredentials <- function(ts_URI, key="sos2-tamis", secOpts=.opts) {
  if(any(strsplit(ts_URI, "/")[[1]] == key))
    return(.opts)
  else
    return(NULL)
}

###

# updateStatus("Loading secOpts")

source("~/52North/secOpts.R")

# wps.off;

timeseries <- "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/592 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/593 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/594 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/595 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/584 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/585 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/586 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/587 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/588 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/589 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/590 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/596 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/597 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/578 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/579 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/580 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/581 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/582 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/583 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/577 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/591 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/575 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/576"
  # "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/514 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/515 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/470 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/473 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/474 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/476 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/479 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/482" # "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/450 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/451 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/452 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/453 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/454 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/455"
timespan <-  "2017-01-31T09:26:42.339Z/2017-07-31T08:26:42.339Z"
# "2016-12-28T12:10:15.944Z/2017-06-28T11:10:15.944Z"
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

dataObs_STFDF <- NULL
for (ts in timeseries) { # ts <- timeseries[1]
  source <- readTSdata(ts, timespan, checkCredentials(ts))
  if(is.null(dataObs_STFDF)) {
    dataObs_STFDF <- source
    next;
  }
  dataObs_STFDF <- rbind(dataObs_STFDF, source)
}

# updateStatus("Setting CRS")

dataObs_STFDF@sp@proj4string <- CRS("+init=epsg:4326")

# updateStatus("Requested Data successfully")

n.time <- length(dataObs_STFDF@time)

dataPos <- "dataPos.png"
png(file = dataPos)
tmpDataPos <- stplot(dataObs_STFDF[,1:min(6, n.time)])
print(tmpDataPos)
graphics.off()

# wps.out: dataPos, png;

# updateStatus("Calculate empirical variogram")

dataObs_STFDF@sp <- spTransform(dataObs_STFDF@sp, target@proj4string)
colnames(dataObs_STFDF@sp@coords) <- c("x","y")

colnamesOld <- colnames(dataObs_STFDF@data)
colnames(dataObs_STFDF@data) <- "targetVar"

empVgm <- variogram(targetVar ~ 1, dataObs_STFDF, tlags=0)
empVgm <- empVgm[-1,]
empVgm <- cbind(empVgm, data.frame(dir.hor=rep(0,nrow(empVgm)), dir.ver=rep(0,nrow(empVgm))))
class(empVgm) <- c("gstatVariogram","data.frame")

empVgm <- empVgm[empVgm$np>0,]

if(n.time >= 10) {
  fitVgm <- fit.variogram(empVgm, vgm(median(empVgm$gamma), "Lin", 50))
  tmpPlot <- plot(empVgm, fitVgm)
} else {
  tmpPlot <- plot(empVgm)
}
  
# wps output
vgmFit <- "vgmFit.png"
png(file = vgmFit)
print(tmpPlot)
graphics.off()
# wps.out: vgmFit, png;

# updateStatus("Plot fitted or selected variogram")

targetData <- NULL
targetVar <- NULL

if (n.time >= 10 & !attributes(fitVgm)$singular) {
  for (day in 1:n.time) {
    pred <- krige0(targetVar ~ 1, dataObs_STFDF[,day], target, model=fitVgm, computeVar = T)
    targetData <- cbind(targetData, pred$pred)
    targetVar <- cbind(targetVar, pred$var)
  }
} else {
  for (day in 1:n.time) {
    pred <- idw0(targetVar ~ 1, dataObs_STFDF[,day], target)
    targetData <- cbind(targetData, pred)
    targetVar <- cbind(targetVar, rep(NA, length(pred)))
  }
}

# updateStatus("Interpolated data")

target <- addAttrToGeom(geometry(target), as.data.frame(cbind(targetData, targetVar)))

# updateStatus("Interpolated grid")

if(isGrid) {
  predictions <- "predictions.tiff"
  predVar <- "predVar.tiff"
} else { 
  predictions <- "predictions.csv"
}

if (isGrid) {
  gridded(target) <- TRUE
  writeGDAL(target[,1:(ncol(target@data)/2)], predictions, drivername="GTiff")
  writeGDAL(target[,-(1:(ncol(target@data)/2))], predVar, drivername="GTiff")
} else {
  write.csv(as.data.frame(target[,1:(ncol(target@data)/2)]), predictions)
  write.csv(as.data.frame(target[,-(1:(ncol(target@data)/2))]), predVar)
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
  
  if (is.projected(target_STFDF@sp)) {
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
  } else {
    # lon
    dim.def.nc(nc, "lon", target_STFDF@sp@grid@cells.dim[1])
    var.def.nc(nc, "lon", "NC_DOUBLE", "lon")
    var.put.nc(nc, "lon", unique(target_STFDF@sp@coords[,1]))
    
    att.put.nc(nc, "lon", "units", "NC_CHAR", "degree_east")
    att.put.nc(nc, "lon", "axis", "NC_CHAR", "lon")
    att.put.nc(nc, "lon", "long_name", "NC_CHAR", "longitude")
    
    # lat
    dim.def.nc(nc, "lat", target_STFDF@sp@grid@cells.dim[2])
    var.def.nc(nc, "lat", "NC_DOUBLE", "lat")
    var.put.nc(nc, "lat", unique(target_STFDF@sp@coords[,2]))
    
    att.put.nc(nc, "lat", "units", "NC_CHAR", "degree_north")
    att.put.nc(nc, "lat", "axis", "NC_CHAR", "lat")
    att.put.nc(nc, "lat", "long_name", "NC_CHAR", "latitude")
  }
  
  # time
  dim.def.nc(nc, "time", length(target_STFDF@time))
  var.def.nc(nc, "time", "NC_INT", "time")
  var.put.nc(nc, "time", as.numeric(index(target_STFDF@time)) - as.numeric(index(target_STFDF@time)[1]))
  
  att.put.nc(nc, "time", "units", "NC_CHAR", paste("secs since",  index(target_STFDF@time)[1]))
  att.put.nc(nc, "time", "axis", "NC_CHAR", "t")
  att.put.nc(nc, "time", "calendar", "NC_CHAR", "gregorian")
  att.put.nc(nc, "time", "long_name", "NC_CHAR", "time")
  att.put.nc(nc, "time", "standard_name", "NC_CHAR", "time")
  
  # CRS
  epsgNum <- as.numeric(showEPSG(target_STFDF@sp@proj4string@projargs))
  
  var.def.nc(nc, "crs", "NC_INT", NA)
  att.put.nc(nc, "crs", "missing_value", "NC_INT", 0)
  if(is.na(epsgNum))
    epsgNum <- 0
  var.put.nc(nc, "crs", epsgNum)
  att.put.nc(nc, "crs", "EPSG_code", "NC_CHAR", paste("EPSG", epsgNum, sep=":"))
  att.put.nc(nc, "crs", "proj4_params", "NC_CHAR", target_STFDF@sp@proj4string@projargs)
  
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