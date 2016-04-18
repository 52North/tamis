# wps.des: tamis.interpolation, title = Interpolation of Wasserstand at Bevertalsperre;

# wps.in: SOSreqData, string, SOS-request, 
# abstract = SOS-request for data,
# value = http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Niederschlagshoehe&procedure=Tagessumme&featureOfInterest=Bever-Talsperre&&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2015-03-10T13:58:07.519Z%2F2016-03-10T13:58:07.519Z;

# wps.in: target, type = geotiff, 
# abstract = SOS-request for Fuellstand,
# value = http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Speicherfuellstand&procedure=Einzelwert&featureOfInterest=Bever-Talsperre_Windenhaus&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-03-01T10:00:00.00Z%2F2016-03-10T13:00:00.000Z;

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



## 
# updateStatus("Requesting SOS")



# updateStatus(print(class(target)))

## tamis
library(sos4R)
library(spacetime)
library(gstat)
library(rgdal)

source("~/52North/secOpts.R")


SOSreqData <- "http://fluggs.wupperverband.de/sos2-tamis/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Schuettmenge&procedure=Tageswert_Prozessleitsystem&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T10:00:00.00Z"

# updateStatus(print(SOSreqData))

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

dataBreakUp <- strsplit(SOSreqData,split = "?", fixed = T)[[1]]
dataURL <- dataBreakUp[1] 

dataBreakUp <- lapply(strsplit(dataBreakUp[2], "&", fixed=T)[[1]], function(x) strsplit(x, "=", fixed=T)[[1]])
dataVersion <- dataBreakUp[[match("version", sapply(dataBreakUp, function(x) x[1]))]][2]

parList <- SOSreqBreakup(dataBreakUp)

source("~/52North/secOpts.R")
TaMIS_SOS <- SOS(url = dataURL,
                 version = dataVersion, binding = "KVP", curlOptions = .opts)

parList$sos <- TaMIS_SOS

# add temporal buffer +/- 12 hours
parList$eventTime <- paste("om%3AphenomenonTime%2C",
                           paste(as.character(as.POSIXct(strsplit(parList$eventTime, "%2C")[[1]][2], 
                                                         format="%Y-%m-%dT%H:%M:%S") + c(-24*3600-5,24*3600+5),
                                              format="%Y-%m-%dT%H:%M:%S"), collapse = "%2F"),
                           sep="")

if(is.null(parList$offering)) {
  TaMIS_offs <- sosOfferings(TaMIS_SOS)
  fstOccur <- min(which(sapply(TaMIS_offs, function(x) any(x@observableProperty == parList$observedProperty[[1]]))))
  parList$offering <- names(TaMIS_offs)[4]
}


# lapply(dataObs, function(x) x@result)

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

empVgm <- variogram(Schuettmenge~1, dataObs_STFDF[,sample(n.time,min(30,n.time)),drop=F],
                    tlags=0, boundaries=c(1:16*5-2.5))
empVgm <- empVgm[-1,]
empVgm <- cbind(empVgm, data.frame(dir.hor=rep(0,nrow(empVgm)), dir.ver=rep(0,nrow(empVgm))))
class(empVgm) <- c("gstatVariogram","data.frame")

empVgm <- empVgm[empVgm$np>0,]
fitVgm <- fit.variogram(empVgm, vgm(median(empVgm$gamma),"Lin",60))

#### wps output
vgmFit <- "vgmFit.png"
png(file = vgmFit)
tmpPlot <- plot(empVgm, fitVgm)
print(tmpPlot)
graphics.off()
# wps.out: vgmFit, png;

# updateStatus("Fitted variogram")

target <- "geotiff.tiff"

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

targetData <- NULL
targetVar <- NULL
for (day in 1:length(dataObs_STFDF@time)) {
  pred <- krige(Schuettmenge ~ 1, dataObs_STFDF[,day], target)@data # , model=fitVgm
  targetData <- cbind(targetData, pred$var1.pred)
  targetVar <- cbind(targetVar, pred$var1.var)
}

target_STFDF <- STFDF(target, dataObs_STFDF@time, 
                      data.frame(var1.pred=as.numeric((targetData)),
                                 var1.var =as.numeric((targetVar))))

# updateStatus("Interpolated grid")

if(isGrid) {
  predictions <- "predictions.tiff"
} else { 
  predictions <- "predictions.csv"
}

if (isGrid) {
  gridded(target_STFDF@sp) <- TRUE
  writeGDAL(target_STFDF[,1], predictions,
            drivername="GTiff")#, type="Byte", options=NULL)
} else {
  write.csv(as.data.frame(target_STFDF), predictions)
}

# wps.out: predictions, type = geotiff;

#   <ows:Title>Plot of the target observations</ows:Title>
#   <ows:Identifier>targetObs_plot</ows:Identifier>

#### wps output
predMap <- "predMap.png"
png(file = predMap)
tmpPlot <- stplot(target_STFDF[,order(sample(n.time, min(n.time, 9))),"var1.pred"])
print(tmpPlot)
graphics.off()
# wps.out: predMap, png;