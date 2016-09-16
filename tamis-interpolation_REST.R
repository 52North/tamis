# wps.des: tamis-rest-interpolation, title = Interpolation of Schüttmenge at Bevertalsperre;

# POST oder GETrequests, wie werden die dAten geliefert?

# wps.in: timeseries, string, set of TS URIs, whitespace " " seperated, 
# abstract = timeseries URI as datasource,
# value = "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/464 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/465 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/466 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/467 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/468 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/469"

# wps.in: timespan, string, timespan of input data, 
# abstract = timeseries URI for the interpolation variable,
# value = "2016-01-01T/2016-01-07TZ";

# wps.in: target, type = geotiff, 
# abstract = geotiff defining the interpolation grid (only non NAs will be interpolated);

# wps.in: targetSOS, string, SOS, 
# abstract = target SOS-URL for the output;

# updateStatus("Requesting SOS")

## tamis
library(httr)
library(rjson)
library(sp)
library(spacetime)
library(gstat)
library(rgdal)
library(RCurl)

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

source("~/52North/secOpts.R")

# wps.off;
timeseries <- "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/464 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/465 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/466 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/467 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/468 http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/469"
timespan <-  "2016-01-01T/2016-01-14TZ"
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

#

timeseries <- strsplit(timeseries, split = " ", fixed = T)[[1]]
timeseries <- timeseries[nchar(timeseries) > 0]

dataObs_STFDF <- NULL
for (ts in timeseries) { # ts <- timeseries[1]
  source <- readTSdata(ts, timespan, .opts)
  if(is.null(dataObs_STFDF)) {
    dataObs_STFDF <- source
    next;
  }
  dataObs_STFDF <- rbind(dataObs_STFDF, source)
}

dataObs_STFDF@sp@proj4string <- CRS("+init=epsg:4326")

# targetVarMeta <- readTSmeta(timeseries_Zielvariable, .opts)  

# updateStatus("Requested Data successfully")

n.time <- length(dataObs_STFDF@time)

dataPos <- "dataPos.png"
png(file = dataPos)
tmpDataPos <- stplot(dataObs_STFDF[,1:min(6, n.time)])
print(tmpDataPos)
graphics.off()

# wps.out: dataPos, png;

empVgm <- variogram(Schüttmenge ~ 1, dataObs_STFDF, tlags=0)
empVgm <- empVgm[-1,]
empVgm <- cbind(empVgm, data.frame(dir.hor=rep(0,nrow(empVgm)), dir.ver=rep(0,nrow(empVgm))))
class(empVgm) <- c("gstatVariogram","data.frame")

empVgm <- empVgm[empVgm$np>0,]

if(n.time >= 10) {
  fitVgm <- fit.variogram(empVgm, vgm(median(empVgm$gamma),"Lin",60))
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

# updateStatus("Fitted/selected variogram")

dataObs_STFDF@sp <- spTransform(dataObs_STFDF@sp, target@proj4string)
colnames(dataObs_STFDF@sp@coords) <- c("x","y")

targetData <- NULL
targetVar <- NULL

if (n.time >= 10) {
  for (day in 1:n.time) {
    pred <- krige(Schüttmenge ~ 1, dataObs_STFDF[,day], target, model=fitVgm)@data
    targetData <- cbind(targetData, pred$var1.pred)
    targetVar <- cbind(targetVar, pred$var1.var)
  }
} else {
  for (day in 1:n.time) {
    pred <- krige(Schüttmenge ~ 1, dataObs_STFDF[,day], target)@data # , model=fitVgm
    targetData <- cbind(targetData, pred$var1.pred)
    targetVar <- cbind(targetVar, pred$var1.var)
  }
}

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
  target_fullGrid <- target
  fullgrid(target_fullGrid) <- TRUE
  
  target_STFDF <- STFDF(geometry(target_fullGrid),
                        dataObs_STFDF@time, 
                        data.frame(var1.pred = unlist(target_fullGrid@data[,1:(ncol(target_fullGrid)/2)]),
                                   var1.var = unlist(target_fullGrid@data[,-(1:(ncol(target_fullGrid)/2))])),
                        index(dataObs_STFDF@time)+24*3600)

  # wps output
  predMap <- "predMap.png"
  png(file = predMap)
  tmpPlot <- stplot(target_STFDF,
                    sp.layout=list("sp.points", dataObs_STFDF@sp))
  print(tmpPlot)
  graphics.off()
  # wps.out: predMap, png;
  
  ## netCDF
  library(RNetCDF)
  ncFile <- "targetnetcdf.nc"
  
  # ncBioTemp <- open.nc("biotemperature_normalDistr.nc")
  # att.inq.nc(ncBioTemp,
  #            ,)
  #            
  # RNetCDF::file.inq.nc(ncBioTemp)
  # 
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
  var.put.nc(nc, "time", 1:length(target_STFDF@time))
  
  att.put.nc(nc, "time", "units", "NC_CHAR", paste("days since", index(target_STFDF@time)[1]) )
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

# wps.off;
targetSOS <- NA # "https://tamis.dev.52north.org/sos/service"
# wps.on;
