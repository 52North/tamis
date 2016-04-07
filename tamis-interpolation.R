# wps.des: tamis.interpolation, title = Interpolation of Wasserstand at Bevertalsperre;


# wps.in: observedproperty, string, observed property, 
# abstract = the observed property default: Wasserstand_im_Damm,
# value = Wasserstand_im_Damm;




# wps.in: phenomenontime, string, phenomenon time, 
# abstract = the phenomenon time;

##########
## Eingabe Punktgeometrien WKT
##         Gitter: geoTiff
## Eingabe Daten: valid sos-request


## tamis
library(sos4R)
library(spacetime)

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
    warning("The following ids have been dropped as they did not contain any data:", paste(ids[dropIds], "\n"))
    dropIds <- which(sapply(res, is.null))
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

source("~/52North/secOpts.R")
TaMIS_SOS <- SOS(url = "http://fluggs.wupperverband.de/sos2-tamis/service",
                 version = "2.0.0", binding = "KVP", curlOptions = .opts)

TaMIS_offs <- sosOfferings(TaMIS_SOS)
lapply(TaMIS_offs, function(x) x@procedure)
lapply(TaMIS_offs, function(x) x@observableProperty)
TaMIS_offs[[1]]@
# updateStatus("Requested SOS successfully")

# updateStatus("Requesting observedproperty observations")

# wps.off
observedproperty <- "Sohlenwasserdruck"
# foi <- list("Bever-Talsperre_MQA7_Piezometer_Kalkzone",
#             "Bever-Talsperre_MQA1_Piezometer_Wasserseite_Schuettkoerper",
#             "Bever-Talsperre_MQA3_Piezometer_Luftseite",
#             "Bever-Talsperre_MQA4_Piezometer_Luftseite",
#             "Bever-Talsperre_MQA5_Piezometer_Berme")
phenomenontime <- "2015-10-01T00:00/2016-03-10T23:59"
# observedproperty <- "Schuettmenge"
# foi <- "Bever-Talsperre_Sickerwassermessstelle_S2A"
# phenomenontime <- "2016-01-01T00:00/2016-03-10T23:59"
# wps.on

offs <- "Zeitreihen_Tageswert_Prozessleitsystem"
proc <- "Tageswert_Prozessleitsystem"

targetObs <- getObservation(TaMIS_SOS,
                            offering = offs, 
                            procedure = proc,
                            # featureOfInterest=SosFeatureOfInterest(foi),
                            observedProperty = list(observedproperty),
                            eventTime = paste("om:phenomenonTime", phenomenontime, sep=","), 
                            responseFormat = "http://www.opengis.net/om/2.0")
str(targetObs)

targetObs_STFDF <- as.STFDF.list.Om_OMObservation(targetObs)

plot(targetObs_STFDF@sp)
hist(targetObs_STFDF[4,,drop=F]@data$Sohlenwasserdruck)

library(gstat)

empVgm <- variogram(Sohlenwasserdruck~1, targetObs_STFDF[,sample(162,30),drop=F], tlags=0, boundaries=c(0.001,1:5*2))
empVgm <- empVgm[-1,]
empVgm <- cbind(empVgm, data.frame(dir.hor=rep(0,nrow(empVgm)), dir.ver=rep(0,nrow(empVgm))))
class(empVgm) <- c("gstatVariogram","data.frame")

plot(empVgm)

fitVgm <- fit.variogram(empVgm, vgm(150,"Sph",5), fit.method=6)
gstat:::plot.gstatVariogram(empVgm, fitVgm)


fit.variogram(variogram(Sohlenwasserdruck~1, targetObs_STFDF[,sample(162,1)], cutoff=20), #, boundaries=c(0.001,1:5*2)),
              vgm(150, "Sph", 5))

plot(variogram(Sohlenwasserdruck~1, targetObs_STFDF[,sample(162,1)], cutoff=20))

plot(v, m)

demo(krige)

v
empVgm

length(targetObs_STFDF[1,,drop=F]@data$Sohlenwasserdruck)

targetObs_STFDF@time
