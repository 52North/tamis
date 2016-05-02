# wps.des: tamis-interpolation, title = Interpolation of Wasserstand at Bevertalsperre;

# wps.in: sosInputData, string, SOS-request, 
# abstract = SOS-request for data,
# value = http://fluggs.wupperverband.de/sos2-tamis/service?service%3DSOS&version%3D2.0.0&request%3DGetObservation&responseformat%3Dhttp://www.opengis.net/om/2.0&observedProperty%3DSchuettmenge&procedure%3DTageswert_Prozessleitsystem&namespaces%3Dxmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter%3Dom%3AphenomenonTime%2C2016-02-01T10:00:00.00Z%2F2016-02-28T10:00:00.00Z;

# wps.in: target, type = geotiff, 
# abstract = geotiff defining the interpolation grid (only non NAs will be interpolated);

#   <ows:Title>Plot of the target observations</ows:Title>
#   <ows:Identifier>targetObs_plot</ows:Identifier>
#   <ows:Title>Diagrams with model parameters</ows:Title>
#   <ows:Identifier>model_diagnostics</ows:Identifier>
#   <ows:Title>Relations between observed properties</ows:Title>
#   <ows:Identifier>relations</ows:Identifier>
#   <ows:Title>Interpolated Values</ows:Title>
#   <ows:Identifier>interpolated-values</ows:Identifier>                           


as.Spatial.MonitoringPoint <- function(obj, ...) {
  
  .extractCRS <- function(obj) {
    chars <- strsplit(obj@shape@point@pos@srsName, "/", fixed=T)[[1]]
    stopifnot(any(c("EPSG","epsg") %in% chars))
    CRS(paste("+init=epsg:", tail(chars,1), sep = ""))
  }
  
  if("point" %in% slotNames(obj@shape))
    return(SpatialPoints(matrix(rev(as.numeric(strsplit(obj@shape@point@pos@pos, " ", fixed=T)[[1]])), ncol = 2),
                         proj4string = .extractCRS(obj)))
}

as.STFDF.list.Om_OMObservation <- function (obs) {
  sp <- do.call(rbind, lapply(obs, function(x) as.Spatial.MonitoringPoint(x@featureOfInterest@feature)))
  
  res <- lapply(obs, function(x) x@result)
  
  ids <- lapply(obs, function(x) x@featureOfInterest@feature@id)
  
  if(any(sapply(res, is.null))) {
    dropIds <- which(sapply(res, is.null))
    warning("The following ids have been dropped as they did not contain any data:", paste(ids[dropIds], "\n"))
    res <- res[-dropIds]
    sp <- sp[-dropIds]
    ids <- ids[-dropIds]
  }
  
  data <- res[[1]]
  colnames(data)[-1] <- ids[[1]]
  if(length(res)>1) {
    for (df in 2:length(res)) {
      colnames(res[[df]])[-1] <- ids[[df]]
      
      data <- merge(data, res[[df]])
    }
  }
  
  time <- as.POSIXct(data[,1])
  
  data <- data.frame(as.numeric(t(as.matrix(data[,-1]))))
  colnames(data) <- tail(names(obs[[1]]@result),1)
  
  STFDF(sp, time, data)
}

as.SpatialPointsDataFrame.list.OmOM_Observation <- function (obs) {
  sp <- do.call(rbind, lapply(obs, function(x) as.Spatial.MonitoringPoint(x@featureOfInterest@feature)))
  
  res <- as.numeric(sapply(obs, 
                           function(x) {
                             res <- x@result
                             if(is.null(res))
                               res <- NA_character_
                             res
                           }))
  sp <- addAttrToGeom(sp, data.frame(obs=res))
  colnames(sp@data) <- obs[[1]]@observedProperty@href
  return(sp)
}

# updateStatus("Requesting SOS")

## tamis
library(sos4R)
library(spacetime)
library(gstat)
library(rgdal)

source("~/52North/secOpts.R")

# wps.off;
# 2016-02
sosInputData <- "http://fluggs.wupperverband.de/sos2-tamis/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Schuettmenge&procedure=Tageswert_Prozessleitsystem&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-02-01T10:00:00.00Z%2F2016-02-28T10:00:00.00Z"
# 2016-01-01
# sosInputData <- "http://fluggs.wupperverband.de/sos2-tamis/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Schuettmenge&procedure=Tageswert_Prozessleitsystem&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T23:59:00.00Z"
# wps.on;

SOSreqBreakup <- function(sosReq) {
  parList <- NULL
  parList$observedProperty <- list(sosReq[[match("observedProperty", sapply(sosReq, function(x) x[1]))]][2])
  parList$proc <- sosReq[[match("procedure", sapply(sosReq, function(x) x[1]))]][2]
  FOIid <- match("featureOfInterest", sapply(sosReq, function(x) x[1]))
  if(!is.na(FOIid))
    parList$featureOfInterest <- sosReq[[match("featureOfInterest", sapply(sosReq, function(x) x[1]))]][2]
  parList$responseFormat  <- sosReq[[match("responseformat", sapply(sosReq, function(x) x[1]))]][2]
  parList$eventTime <- sosReq[[match("temporalFilter", sapply(sosReq, function(x) x[1]))]][2]
  
  return(parList)
}
sosInputData <- gsub("%3D","=", sosInputData)#=-signs must always be escaped using WPS4R 
sosInputData <- gsub("&amp;","&", sosInputData)#ampersands might be encoded
sosInputData <- gsub("req_quest","request", sosInputData)#request will be replaced by req_est after sent to Rserve
sosInputData <- gsub("s_system","system", sosInputData)#request will be replaced by req_est after sent to Rserve TODO: check other filtered strings

dataBreakUp <- strsplit(sosInputData,split = "?", fixed = T)[[1]]
dataURL <- dataBreakUp[1] 

dataBreakUp <- lapply(strsplit(dataBreakUp[2], "&", fixed=T)[[1]], function(x) strsplit(x, "=", fixed=T)[[1]])
dataVersion <- dataBreakUp[[match("version", sapply(dataBreakUp, function(x) x[1]))]][2]

parList <- SOSreqBreakup(dataBreakUp)

source("~/52North/secOpts.R")
TaMIS_SOS <- SOS(url = dataURL,
                 version = dataVersion, binding = "KVP", curlOptions = .opts)

parList$sos <- TaMIS_SOS

if (length(strsplit(parList$eventTime, split = "%2F")[[1]]) == 1) {
# add temporal buffer +/- 1 day and 5 minutes
  parList$eventTime <- paste("om%3AphenomenonTime%2C",
                             paste(as.character(as.POSIXct(strsplit(parList$eventTime, "%2C")[[1]][2], 
                                                           format="%Y-%m-%dT%H:%M:%S") + c(-24*3600-300,24*3600+300),
                                                format="%Y-%m-%dT%H:%M:%S"), collapse = "%2F"),
                             sep="")
}

if(is.null(parList$offering)) {
  TaMIS_offs <- sosOfferings(TaMIS_SOS)
  fstOccur <- min(which(sapply(TaMIS_offs, function(x) any(x@observableProperty == parList$observedProperty[[1]]))))
  parList$offering <- names(TaMIS_offs)[4]
}

# updateStatus("Requested SOS successfully")

dataObs <- do.call(getObservation, parList)

# updateStatus("Requested observedproperty observations successfully")

dataObs_STFDF <- as.STFDF.list.Om_OMObservation(dataObs)

dataPos <- "dataPos.png"
png(file = dataPos)
tmpDataPos <- plot(dataObs_STFDF@sp)
print(tmpDataPos)
graphics.off()

# wps.out: dataPos, png;

library(gstat)

n.time <- length(dataObs_STFDF@time)

empVgm <- variogram(Schuettmenge ~ 1, dataObs_STFDF, tlags=0)
empVgm <- empVgm[-1,]
empVgm <- cbind(empVgm, data.frame(dir.hor=rep(0,nrow(empVgm)), dir.ver=rep(0,nrow(empVgm))))
class(empVgm) <- c("gstatVariogram","data.frame")

empVgm <- empVgm[empVgm$np>0,]

if(n.time >= 20) {
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

# # updateStatus("Fitted variogram")

# wps.off;
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

dataObs_STFDF@sp <- spTransform(dataObs_STFDF@sp, target@proj4string)
colnames(dataObs_STFDF@sp@coords) <- c("x","y")

targetData <- NULL
targetVar <- NULL

if (n.time >= 28) {
  for (day in 1:n.time) {
    pred <- krige(Schuettmenge ~ 1, dataObs_STFDF[,day], target, model=fitVgm)@data
    targetData <- cbind(targetData, pred$var1.pred)
    targetVar <- cbind(targetVar, pred$var1.var)
  }
} else {
  for (day in 1:n.time) {
    pred <- krige(Schuettmenge ~ 1, dataObs_STFDF[,day], target)@data # , model=fitVgm
    targetData <- cbind(targetData, pred$var1.pred)
    targetVar <- cbind(targetVar, pred$var1.var)
  }
}

target <- addAttrToGeom(geometry(target), data.frame(var1.pred=apply(targetData,1,mean, na.rm=T)))

# updateStatus("Interpolated grid")

if(isGrid) {
  predictions <- "predictions.tiff"
} else { 
  predictions <- "predictions.csv"
}

if (isGrid) {
  gridded(target) <- TRUE
  writeGDAL(target, predictions, drivername="GTiff")
} else {
  write.csv(as.data.frame(target), predictions)
}

# wps.out: predictions, type = geotiff;

# wps output
predMap <- "predMap.png"
png(file = predMap)
tmpPlot <- spplot(target,"var1.pred",
                  sp.layout=list("sp.points",dataObs_STFDF@sp))
print(tmpPlot)
graphics.off()
# wps.out: predMap, png;
