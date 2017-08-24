# wps.des: interpolation-demo, title = Demo interpolation service;

# wps.in: endpoint, type = string,
# abstract = endpoint of the SOS REST-API,
# value = "http://fluggs.wupperverband.de/sos2/api/v1/";

# wps.in: phenomenon, type = string,
# abstract = phenomenon of interest,
# value = "Lufttemperatur";

# wps.in: timespan, type = string,
# abstract = timespan for the interpolation variable,
# value = "2013-01-01T00:00:01Z/2013-01-31T23:59:59Z";

# wps.in: gridColumns, type = integer,
# abstract = number of columns in the grid,
# value = 20;

# wps.in: gridRows, type = integer,
# abstract = number of rows in the grid,
# value = 10;

# updateStatus("Loading libraries")

library(sp)
library(spacetime)
library(rgdal)
library(sensorweb4R)

# updateStatus("Set-up REST-API")


## pre-process 

# wps.off;

# define defaults for local testing
endpoint <- "http://fluggs.wupperverband.de/sos2/api/v1"
phenomenon <- "Lufttemperatur"
timespan <- "2013-01-01T00:00:01Z/2013-01-31T23:59:59Z"
gridColumns <- 10
gridRows <- 5
# wps.on;

# updateStatus("define endpoint")
endpoint <- as.Endpoint(endpoint)

# fetch meta-data
tsMetaData <- fetch(timeseries(endpoint))

tsIds <- which(tsMetaData@phenomenon@label == phenomenon)

# stop if no data is available
if(length(tsIds) == 0)
  stop("No timeseries Ids could be idefied for: \"", phenomenon, "\"")

tsSelected <- tsMetaData[tsIds]

# updateStatus("get data")
tsData <- getData(tsSelected, timespan = timespan)

df <- do.call(rbind, lapply(tsData, as.data.frame))
if(nrow(df) == 0)
  stop("No data avaialable for phenomenon \"", phenomenon, "\" during the timespan \"", timespan, "\"")

# build STSDF
sp <- addAttrToGeom(geometry(tsSelected@station), 
                    data.frame(id = tsSelected@station@id,
                               label = tsSelected@station@label))
row.names(sp) <- tsSelected@station@label
nonEmptyLoc <- which(sapply(tsData, length) >= 0)
sp <- sp[nonEmptyLoc,]
time <- floor(as.numeric(df[,1])/60/15)
ordTime <- order(time)

time <- time*15*60+3600
time <- as.factor(as.POSIXct(time, origin="1970-01-01", tz="GMT"))
index <- cbind(rep(1:length(sp), sapply(tsData, length)),
               as.numeric(time))[ordTime,]
time <- as.POSIXct(levels(time), tz="GMT")

stsdf <- STSDF(sp, time, df[ordTime, 2, drop=F], index)
stfdf <- as(stsdf,"STFDF")

# stplot(stfdf, mode="ts")

# updateStatus("aggregate data")

aggStfdf <- aggregate(stsdf, by="day", FUN=function(x) mean(x, na.rm=T))

# stplot(aggStfdf, mode="ts")

dfFull <- as.data.frame(aggStfdf)
dfFull$timeIndex <- as.factor(dfFull$timeIndex)

summary(lm(value~x+timeIndex, dfFull))

library(gstat)
# updateStatus("calculate empirical pooled spatio-temporal variogram")
empVgm <- variogram(value~x+y, aggStfdf, tlags=0)

# make spatio-temporal variogram pure spatial = "pooled" variogram over all time slices
empVgm <- cbind(empVgm, data.frame(dir.hor=rep(0,nrow(empVgm)), dir.ver=rep(0,nrow(empVgm))))
class(empVgm) <- c("gstatVariogram","data.frame")
empVgm <- empVgm[empVgm$np>39,]

# plot(empVgm)
# stplot(aggStfdf)

qPar <- quantile(empVgm$gamma, probs = c(0.9,0.05))

# updateStatus("fit spatial variogram")

fitVgm <- fit.variogram(empVgm, vgm(qPar[1], "Lin", quantile(empVgm$dist, probs = 0.7), qPar[2]))
# plot(empVgm, fitVgm)

cSize <- apply(aggStfdf@sp@bbox,1,diff)/c(gridColumns-1, gridRows-1)

tarGeom <- GridTopology(aggStfdf@sp@bbox[,1], 
                        cellsize = cSize,
                        cells.dim = c(gridColumns, gridRows))
tarGeom <- SpatialGrid(tarGeom, proj4string = aggStfdf@sp@proj4string)

tarGeom <- as(tarGeom, "SpatialPixels")

fullgrid(tarGeom) <- TRUE

res <- NULL 

aggStsdf <- as(aggStfdf, "STSDF")

# updateStatus("interpolate grid")

if (any(fitVgm$psill < 0) | attributes(fitVgm)$singular) {
  for (t in aggStfdf@time) {
    res <- rbind(res, 
                 idw(value~x+y, aggStsdf[,t], tarGeom)@data)
  }
} else {
  for (t in aggStfdf@time) {
    res <- rbind(res, 
                 krige(value~x+y, aggStsdf[,t], tarGeom, model=fitVgm)@data)
  }
}

resStfdf <- STFDF(tarGeom, aggStfdf@time, res)

# stplot(resStfdf)

# updateStatus("produce output")

# wps output
predMap <- "predMap.png"
png(file = predMap)
tmpPlot <- stplot(resStfdf,# xlim=bbox(tarGeom)[1,], ylim=bbox(tarGeom)[2,],
                  sp.layout=list("sp.points", aggStfdf@sp))
print(tmpPlot)
graphics.off()
# wps.out: predMap, png;

resSpPix <- addAttrToGeom(resStfdf@sp, as.data.frame(matrix(resStfdf@data$var1.pred, ncol = length(aggStfdf@time))))

predictions <- "predictions.tiff"

writeGDAL(resSpPix, predictions, drivername="GTiff")

# wps.out: predictions, type = geotiff;