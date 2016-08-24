# wps.des: tamis-regression, title = TaMIS Regression Model for Wasserstand_im_Damm or Schuettmenge at Bever-Talsperre;

# wps.in: sosInputNiederschlag, string, SOS-request, 
# abstract = SOS-request for Niederschlag, minOccurs = 0, maxOccurs = 1,
# value = http://www.fluggs.de/sos2/sos?service%3DSOS&version%3D2.0.0&request%3DGetObservation&responseformat%3Dhttp://www.opengis.net/om/2.0&observedProperty%3DNiederschlagshoehe&procedure%3DTagessumme&featureOfInterest%3DBever-Talsperre&&namespaces%3Dxmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter%3Dom%3AphenomenonTime%2C2016-01-01T10:00:00.00Z%2F2016-04-30T23:59:00.000Z;

# wps.in: sosInputFuellstand, string, SOS-request, 
# abstract = SOS-request for Fuellstand, minOccurs = 0, maxOccurs = 1,
# value = http://www.fluggs.de/sos2/sos?service%3DSOS&version%3D2.0.0&request%3DGetObservation&responseformat%3Dhttp://www.opengis.net/om/2.0&observedProperty%3DSpeicherfuellstand&procedure%3DEinzelwert&featureOfInterest%3DBever-Talsperre_Windenhaus&namespaces%3Dxmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter%3Dom%3AphenomenonTime%2C2016-01-01T10:00:00.00Z%2F2016-04-30T23:59:00.000Z;

# wps.in: sosInputTarget, string, SOS-request, 
# abstract = SOS-request for the target variable,
# value = http://fluggs.wupperverband.de/sos2-tamis/service?service%3DSOS&version%3D2.0.0&request%3DGetObservation&responseformat%3Dhttp://www.opengis.net/om/2.0&observedProperty%3DWasserstand_im_Damm&procedure%3DHandeingabe&featureOfInterest%3DBever-Talsperre_MQA7_Piezometer_Kalkzone&namespaces%3Dxmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpati-al%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter%3Dom%3AphenomenonTime%2C2016-01-01T00:01:00.00Z%2F2016-04-30T23:59:00.000Z;

# wps.in: sosInputNiederschlagPred, string, SOS-request, 
# abstract = SOS-request for prediction values: Niederschlag, minOccurs = 0, maxOccurs = 1;

# wps.in: sosInputFuellstandPred, string, SOS-request, 
# abstract = SOS-request for prediction values: Fuellstand, minOccurs = 0, maxOccurs = 1;

# wps.in: singleInputNiederschlagPred, integer, single value, 
# abstract = single value for prediction values: Niedreschlag, minOccurs = 0, maxOccurs = 1;

# wps.in: singleInputFuellstandPred, integer, single value, 
# abstract = single value for prediction values: Fuellstand, minOccurs = 0, maxOccurs = 1;


## tamis
library(sos4R)
library(spacetime)
library(rjson)

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
    
  data <- as.data.frame(res[[1]])
  if (ncol(data) > 1) {
    colnames(data)[-1] <- ids[[1]]
    if(length(res)>1) {
      for (df in 2:length(res)) {
        colnames(res[[df]])[-1] <- ids[[df]]
        
        data <- merge(data, res[[df]])
      }
    }
    
    # check for NAs in time column (daylight saving time artefacts)
    data <- data[!is.na(data[,1]),]
    
    time <- as.POSIXct(data[,1])
    
    data <- data.frame(as.numeric(t(as.matrix(data[,-1]))))
    colnames(data) <- tail(names(obs[[1]]@result),1)
    
    ret <- STFDF(sp, time, data)
  } else {
    time <- do.call(c, lapply(obs, function(x) (x@phenomenonTime@timePosition@time)))
    data <- data.frame(do.call(c, res))
    ret <- STIDF(sp, time, data)
    ret <- as(ret,"STFDF")
  }
  
  ret
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

SOSreqBreakup <- function(sosReq) {
  observedproperty <- sosReq[[match("observedProperty", sapply(sosReq, function(x) x[1]))]][2]
  proc <- sosReq[[match("procedure", sapply(sosReq, function(x) x[1]))]][2]
  featureOfInterest <- sosReq[[match("featureOfInterest", sapply(sosReq, function(x) x[1]))]][2]
  responseFormat  <- sosReq[[match("responseformat", sapply(sosReq, function(x) x[1]))]][2]
  eventTime <- sosReq[[match("temporalFilter", sapply(sosReq, function(x) x[1]))]][2]
  
  return(list(procedure = proc,
              featureOfInterest = SosFeatureOfInterest(list(featureOfInterest)),
              observedProperty = list(observedproperty),
              eventTime = eventTime, 
              responseFormat = responseFormat))
}

########### Script
## 
# updateStatus("Requesting SOS")

# TerraTransfer:
# http://tamis-sos.de:8080/52n-sos/service

# wps.off;
sosInputNiederschlag <- "http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Niederschlagshoehe&procedure=Tagessumme&featureOfInterest=Bever-Talsperre&&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T10:00:00.00Z%2F2016-04-30T23:59:00.000Z"

sosInputFuellstand <- "http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Speicherfuellstand&procedure=Einzelwert&featureOfInterest=Bever-Talsperre_Windenhaus&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T10:00:00.00Z%2F2016-04-30T23:59:00.000Z"

# TT:
# sosInputTarget <- "http://tamis-sos.de:8080/52n-sos/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=http://www.52north.org/test/observableProperty/waterLevel&procedure=WaterLevelE10006&featureOfInterest=water%20level%20sensor&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T00:00:01.00Z%2F2016-06-30T23:59:59.000Z&mergeObservationsIntoDataArray=true"
# sosInputTarget <- "http://tamis-sos.de:8080/52n-sos/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=http://www.52north.org/test/observableProperty/waterLevel&procedure=WaterLevelE10013&featureOfInterest=water%20level%20sensor&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T00:00:01.00Z%2F2016-06-30T23:59:59.000Z&mergeObservationsIntoDataArray=true"
# sosInputTarget <- "http://tamis-sos.de:8080/52n-sos/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=http://www.52north.org/test/observableProperty/waterLevel&procedure=WaterLevelE10014&featureOfInterest=water%20level%20sensor&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T00:00:01.00Z%2F2016-06-30T23:59:59.000Z&mergeObservationsIntoDataArray=true"
# sosInputTarget <- "http://tamis-sos.de:8080/52n-sos/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=http://www.52north.org/test/observableProperty/waterLevel&procedure=WaterLevelE10015&featureOfInterest=water%20level%20sensor&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T00:00:01.00Z%2F2016-06-30T23:59:59.000Z&mergeObservationsIntoDataArray=true"
sosInputTarget <- "http://tamis-sos.de:8080/52n-sos/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=http://www.52north.org/test/observableProperty/waterLevel&procedure=WaterLevelE10019&featureOfInterest=water%20level%20sensor&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T00:00:01.00Z%2F2016-06-30T23:59:59.000Z&mergeObservationsIntoDataArray=true"

# WV: 
# sosInputTarget <- "http://fluggs.wupperverband.de/sos2-tamis/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Wasserstand_im_Damm&procedure=Handeingabe&featureOfInterest=Bever-Talsperre_MQA7_Piezometer_Kalkzone&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T00:01:00.00Z%2F2016-04-30T23:59:00.000Z"
# Bever-Talsperre_MQA1_Piezometer_Wasserseite_Schuettkoerper
# Bever-Talsperre_MQA3_Piezometer_Luftseite
# Bever-Talsperre_MQA4_Piezometer_Luftseite
# Bever-Talsperre_MQA5_Piezometer_Berme
# Bever-Talsperre_MQA7_Piezometer_Kalkzone

# sosInputTarget <- "http://fluggs.wupperverband.de/sos2-tamis/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Schuettmenge&procedure=Tageswert_Prozessleitsystem&featureOfInterest=Bever-Talsperre_Sickerwassermessstelle_S2B&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2015-01-01T10:00:00.00Z%2F2016-04-30T23:59:00.000Z"
# Bever-Talsperre_Sickerwassermessstelle_S2A
# Bever-Talsperre_Sickerwassermessstelle_S2B

sosInputNiederschlagPred <- NA # "http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Niederschlagshoehe&procedure=Tagessumme&featureOfInterest=Bever-Talsperre&&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T10:00:00.00Z%2F2016-04-30T23:59:00.000Z"

sosInputFuellstandPred <- NA # "http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Speicherfuellstand&procedure=Einzelwert&featureOfInterest=Bever-Talsperre_Windenhaus&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T10:00:00.00Z%2F2016-04-30T23:59:00.000Z"

singleInputNiederschlagPred <- NA # 10

singleInputFuellstandPred <- NA # 293
# wps.on;

TT <- FALSE
if (substr(sosInputTarget, 8, 19) == "tamis-sos.de")
  TT <- TRUE

WV <- FALSE
if (substr(sosInputTarget, 8, 30) == "fluggs.wupperverband.de")
  WV <- TRUE


if(is.na(sosInputNiederschlag) & is.na(sosInputNiederschlagPred) & is.na(singleInputNiederschlagPred))
  stop("At least one input type for Niederschlag must be present.")
# TODO: stop WPS4R execution on first error??

if(is.na(sosInputFuellstand) & is.na(sosInputFuellstandPred) & is.na(singleInputFuellstandPred))
  stop("At least one input type for Fuellstand must be present.")

if ((is.na(sosInputFuellstand) | is.na(sosInputNiederschlag)) & (is.na(sosInputFuellstandPred) | is.na(sosInputNiederschlagPred)) & (is.na(singleInputFuellstandPred) | is.na(singleInputNiederschlagPred)))
  stop("At least one pair of input data must be present.")

sosInputTarget <- gsub("%3D","=", sosInputTarget)
sosInputTarget <- gsub("&amp;","&", sosInputTarget)
sosInputTarget <- gsub("req_quest","request", sosInputTarget)#request will be replaced by req_est after sent to Rserve
sosInputTarget <- gsub("s_system","system", sosInputTarget)#request will be replaced by req_est after sent to Rserve TODO: check other filtered strings

targetBreakUp <- strsplit(sosInputTarget,split = "?", fixed = T)[[1]]
targetURL <- targetBreakUp[1] 

targetBreakUp <- lapply(strsplit(targetBreakUp[2], "&", fixed=T)[[1]], function(x) strsplit(x, "=", fixed=T)[[1]])
targetVersion <- targetBreakUp[[match("version", sapply(targetBreakUp, function(x) x[1]))]][2]

if(WV) {
  source("~/52North/secOpts.R")
  TaMIS_SOS <- SOS(url = targetURL,
                   version = targetVersion, binding = "KVP", 
                   dataFieldConverters = SosDataFieldConvertingFunctions("l/s" = sosConvertDouble),
                   curlOptions = .opts)
}

if (TT) {
  TaMIS_SOS <- SOS(url = targetURL,
                   version = targetVersion, binding = "KVP",
                   dataFieldConverters = SosDataFieldConvertingFunctions("mNHN" = sosConvertDouble))
}

# updateStatus("Requested SOS successfully")

# updateStatus("Requesting observedproperty observations")

parList <- SOSreqBreakup(targetBreakUp)

# both WV scenarios
if(parList$observedProperty == "Wasserstand_im_Damm") {
  parList$offering <- "Zeitreihen_Handeingabe"
}

if(parList$observedProperty == "Schuettmenge") {
  TaMIS_SOS@dataFieldConverters$'l/s' <- function (x, sos) {
    return(as.double(x = x))}
  parList$offering <- "Zeitreihen_Tageswert_Prozessleitsystem"
}

# TT scenario
if(substr(parList$procedure, 1, 10) == "WaterLevel") {
  parList$offering <- paste("WaterLevelOffering", 
                            substr(parList$procedure, 11, nchar(parList$procedure)), sep="")
}

# merge obs for TT, remove after default!
if (TT)
  TaMIS_SOS@additionalKVPs <- list(mergeObservationsIntoDataArray = TRUE)

# add SOS to parameters
parList$sos <- TaMIS_SOS

# get obs for the list of the parameters/arguments in parList
targetObs <- do.call(getObservation, parList)

# updateStatus("Requested observedproperty observations successfully")

targetObs_STFDF <- as.STFDF.list.Om_OMObservation(targetObs)

# updateStatus("targetObs_STFDF succesfully created.")

# updateStatus("Plotting observedproperty:")

if(substr(parList$procedure, 1, 10) == "WaterLevel" & TT) {
  colnames(targetObs_STFDF@data) <- "WaterLevel"

## aggregate values to daily averages
targetObs_STFDF <- STFDF(targetObs_STFDF@sp,
                         unique(as.Date(index(targetObs_STFDF@time), format = "%D")),
                         as.data.frame(aggregate(targetObs_STFDF, STF(targetObs_STFDF@sp, 
                                                  unique(as.Date(index(targetObs_STFDF@time), format = "%D")),
                                                  endTime = as.POSIXct(unique(as.Date(index(targetObs_STFDF@time), format = "%D"))+1)),
                             mean)[,1]))
}

targetObs_plot <- "targetObs_plot.png"
png(file = targetObs_plot)

tmp <- stplot(targetObs_STFDF, mode="ts", typ="b")
print(tmp)

graphics.off()

# wps.out: targetObs_plot, png;

if(parList$observedProperty == "Wasserstand_im_Damm" | TT) {
  times <- as.character(index(targetObs_STFDF@time))
} else {
  times <- as.character(index(targetObs_STFDF@time))
  times <- sapply(strsplit(times, " "), function(x) x[1])
} 

## Fuellstand
if(!is.na(sosInputFuellstand)) {
  sosInputFuellstand <- gsub("%3D","=", sosInputFuellstand)
  sosInputFuellstand <- gsub("&amp;","&", sosInputFuellstand)
  sosInputFuellstand <- gsub("req_quest","request", sosInputFuellstand)#request will be replaced by req_est after sent to Rserve
  sosInputFuellstand <- gsub("s_system","system", sosInputFuellstand)#request will be replaced by req_est after sent to Rserve TODO: check other filtered strings
  
  fuellstandBreakUp <- strsplit(sosInputFuellstand,split = "?", fixed = T)[[1]]
  
  fuellstandURL <- fuellstandBreakUp[1] 
  
  fuellstandBreakUp <- lapply(strsplit(fuellstandBreakUp[2], "&", fixed=T)[[1]], function(x) strsplit(x, "=", fixed=T)[[1]])
  
  fuellstandVersion <- fuellstandBreakUp[[match("version", sapply(fuellstandBreakUp, function(x) x[1]))]][2]
  
  FLUGGS_SOS <- SOS(url = fuellstandURL,
                    version = fuellstandVersion, binding = "KVP")
  
  parList <- SOSreqBreakup(fuellstandBreakUp)
  
  parList$sos <- FLUGGS_SOS
  parList$offering <- "Zeitreihen_Einzelwert"
  
  if(is.na(sosInputFuellstandPred)) {
    fuellstandPred <- do.call(getObservation, parList)
  }
  
  fuellstand <- vector("list", length(times))
  for(i in 1:length(times)) {
    parList$eventTime <- paste("om:phenomenonTime,",times[i],"T10:00",sep="")
    fuellstand[[i]] <- do.call(getObservation, parList)
  }
  
  fuellstandVec <- as.numeric(sapply(fuellstand, function(x) {
    res <- x[[1]]@result
    if(is.null(res))
      return(NA)
    else 
      return(res)
  }))
} else {
  fuellstandVec <- NULL
}

# check for prediction, prefer single value
if (!is.na(singleInputFuellstandPred)) {
  fuellstandPredVec <- singleInputFuellstandPred
} else {
  if (!is.na(sosInputFuellstandPred)) {
    sosInputFuellstandPred <- gsub("%3D","=", sosInputFuellstandPred)
    sosInputFuellstandPred <- gsub("&amp;","&", sosInputFuellstandPred)
    sosInputFuellstandPred <- gsub("req_quest","request", sosInputFuellstandPred)#request will be replaced by req_est after sent to Rserve
    sosInputFuellstandPred <- gsub("s_system","system", sosInputFuellstandPred)#request will be replaced by req_est after sent to Rserve TODO: check other filtered strings
    fuellstandBreakUp <- strsplit(sosInputFuellstandPred,split = "?", fixed = T)[[1]]
    fuellstandURL <- fuellstandBreakUp[1] 
    fuellstandBreakUp <- lapply(strsplit(fuellstandBreakUp[2], "&", fixed=T)[[1]], function(x) strsplit(x, "=", fixed=T)[[1]])
    fuellstandVersion <- fuellstandBreakUp[[match("version", sapply(fuellstandBreakUp, function(x) x[1]))]][2]
    FLUGGS_SOS <- SOS(url = fuellstandURL,
                      version = fuellstandVersion, binding = "KVP")
    
    parList <- SOSreqBreakup(fuellstandBreakUp)
    parList$sos <- FLUGGS_SOS
    parList$offering <- "Zeitreihen_Einzelwert"
    
    fuellstandPred <- do.call(getObservation, parList)
  }
}

### Niederschlag
if (!is.na(sosInputNiederschlag)) {
  sosInputNiederschlag <- gsub("%3D","=", sosInputNiederschlag)
  sosInputNiederschlag <- gsub("&amp;","&", sosInputNiederschlag)
  sosInputNiederschlag <- gsub("req_quest","request", sosInputNiederschlag)#request will be replaced by req_est after sent to Rserve
  sosInputNiederschlag <- gsub("s_system","system", sosInputNiederschlag)#request will be replaced by req_est after sent to Rserve TODO: check other filtered strings
  
  niederschlagBreakUp <- strsplit(sosInputNiederschlag,split = "?", fixed = T)[[1]]
  
  niederschlagURL <- niederschlagBreakUp[1] 
  
  niederschlagBreakUp <- lapply(strsplit(niederschlagBreakUp[2], "&", fixed=T)[[1]], function(x) strsplit(x, "=", fixed=T)[[1]])
  
  niederschlagVersion <- niederschlagBreakUp[[match("version", sapply(niederschlagBreakUp, function(x) x[1]))]][2]
  
  FLUGGS_SOS <- SOS(url = niederschlagURL,
                    version = niederschlagVersion, binding = "KVP")
  
  parList <- SOSreqBreakup(niederschlagBreakUp)
  
  parList$sos <- FLUGGS_SOS
  parList$offering <- "Zeitreihen_Tagessumme"
  
  if (is.na(sosInputNiederschlagPred) & is.na(singleInputNiederschlagPred)) {
    niederschlagPred <- do.call(getObservation, parList)
  }
  
  # get data
  niederschlag <- vector("list",length = length(times))
  for ( i in 1:length(times)) {
    parList$eventTime <- paste("om:phenomenonTime,", times[i], "T00:00", sep="")
    niederschlag[[i]] <- do.call(getObservation, parList)
  }
  
  niederschlagVec <-as.numeric(sapply(niederschlag, 
                                      function(x) {
    res <- x[[1]]@result
    if(is.null(res))
      return(NA)
    else 
      return(res)
    }))
  
  niederschlagVec[niederschlagVec > 200] <- NA
} else {
  niederschlagVec <- NULL
}

# check for prediction, prefer single value
if (!is.na(singleInputNiederschlagPred)) {
  niederschlagPredVec <- singleInputNiederschlagPred
} else {
  if (!is.na(sosInputNiederschlagPred)) {
    sosInputNiederschlagPred <- gsub("%3D","=", sosInputNiederschlagPred)
    sosInputNiederschlagPred <- gsub("&amp;","&", sosInputNiederschlagPred)
    sosInputNiederschlagPred <- gsub("req_quest","request", sosInputNiederschlagPred)#request will be replaced by req_est after sent to Rserve
    sosInputNiederschlagPred <- gsub("s_system","system", sosInputNiederschlagPred)#request will be replaced by req_est after sent to Rserve TODO: check other filtered strings
    niederschlagBreakUp <- strsplit(sosInputNiederschlagPred,split = "?", fixed = T)[[1]]
    niederschlagURL <- niederschlagBreakUp[1] 
    niederschlagBreakUp <- lapply(strsplit(niederschlagBreakUp[2], "&", fixed=T)[[1]], function(x) strsplit(x, "=", fixed=T)[[1]])
    niederschlagVersion <- niederschlagBreakUp[[match("version", sapply(niederschlagBreakUp, function(x) x[1]))]][2]
    FLUGGS_SOS <- SOS(url = niederschlagURL,
                      version = niederschlagVersion, binding = "KVP")
    
    parList <- SOSreqBreakup(niederschlagBreakUp)
    parList$sos <- FLUGGS_SOS
    parList$offering <- "Zeitreihen_Tagessumme"
    
    niederschlagPred <- do.call(getObservation, parList)
  }
}

# match covariates
if(is.na(singleInputNiederschlagPred) & is.na(singleInputFuellstandPred)) {
  mIds <- match(strptime(niederschlagPred$observationData@result$phenomenonTime, 
                         format="%Y-%m-%d %H:%M"), 
                strptime(fuellstandPred$observationData@result$phenomenonTime, 
                format="%Y-%m-%d %H:%M"))
  
  precipIds <- which(niederschlagPred$observationData@result$Niederschlagshoehe > 200)
  mIds[precipIds] <- NA
  
  fuellstandPredVec <- as.numeric(fuellstandPred$observationData@result$Speicherfuellstand)
  fuellstandPredVec <- fuellstandPredVec[mIds[!is.na(mIds)]]
  
  niederschlagPredVec <- as.numeric(niederschlagPred$observationData@result$Niederschlagshoehe)
  niederschlagPredVec <- niederschlagPredVec[!is.na(mIds)]

  predTimes <- strptime(fuellstandPred$observationData@result$phenomenonTime, 
                        format="%Y-%m-%d %H:%M")[mIds]
  predTimes <- as.numeric(predTimes[!is.na(predTimes)])
} else {
  predTimes <- as.numeric(Sys.time())
}

niederschlagPredVec[niederschlagPredVec > 200] <- NA
# modelling

targetVec <- targetObs_STFDF@data[[1]]

if (WV & (length(times) < 10 | any(c(is.na(sosInputFuellstand),
                              is.na(sosInputNiederschlag))))) {
  load("preDefModel.RData")
  lmMod <- switch(targetBreakUp[[which(lapply(targetBreakUp, function(x) x[[1]]) == "featureOfInterest")]][2],
                  'Bever-Talsperre_MQA1_Piezometer_Wasserseite_Schuettkoerper' = MQA1mod,
                  'Bever-Talsperre_MQA3_Piezometer_Luftseite' = MQA3mod,
                  'Bever-Talsperre_MQA4_Piezometer_Luftseite' = MQA4mod,
                  'Bever-Talsperre_MQA5_Piezometer_Berme' = MQA5mod,
                  'Bever-Talsperre_MQA7_Piezometer_Kalkzone' = MQA7mod,
                  'Bever-Talsperre_Sickerwassermessstelle_S2A' = S2Amod,
                  'Bever-Talsperre_Sickerwassermessstelle_S2B' = S2Bmod)
}
if (TT & (length(times) < 10 | any(c(is.na(sosInputFuellstand),
                                     is.na(sosInputNiederschlag))))) {
  load("preDefModel.RData")
  lmMod <- switch(targetBreakUp[[which(lapply(targetBreakUp, function(x) x[[1]]) == "procedure")]][2],
                  'WaterLevelE10006' = E10006,
                  'WaterLevelE10013' = E10013,
                  'WaterLevelE10014' = E10014,
                  'WaterLevelE10015' = E10015,
                  'WaterLevelE10019' = E10019)
} else {
  lmMod <- lm(targetVec ~ fuellstandVec + niederschlagVec)
  # summary(lmMod)
  
  # MQA1mod <- lmMod
  # MQA3mod <- lmMod
  # MQA4mod <- lmMod
  # MQA5mod <- lmMod
  # MQA7mod <- lmMod
  
  # E10006 <- lmMod
  # E10013 <- lmMod
  # E10014 <- lmMod
  # E10015 <- lmMod
  # E10019 <- lmMod
  
  # S2Amod <- lmMod
  # S2Bmod <- lmMod
  # save(MQA1mod, MQA3mod, MQA4mod, MQA5mod, MQA7mod,
  #      S2Amod, S2Bmod,
  #      E10006, E10013, E10014, E10015, E10019, 
  #      file="preDefModel.RData")
}

model_diagnostics <- "model_diagnostics.png"
png(file = model_diagnostics)#TODO vielleicht mehr parameter wie groesse etc

par(mfrow=c(2,2))
plot(lmMod)

graphics.off()

df <- data.frame(niederschlagVec=niederschlagPredVec, fuellstandVec=fuellstandPredVec)

df$targetVec <- predict(lmMod, df)

df <- df[!apply(df, 1, function(x) any(is.na(x))),]
  
# wps.out: model_diagnostics, png;

library(lattice)
p1 <- xyplot(df$targetVec ~ df$niederschlagVec,
       xlab="Niederschlag",
       ylab=colnames(targetObs_STFDF@data))
p2 <- xyplot(df$targetVec ~  df$fuellstandVec, 
       xlab="Fuellstand",
       ylab=colnames(targetObs_STFDF@data))

p3 <- xyplot(predict(lmMod, data.frame(niederschlagVec=niederschlagVec, fuellstandVec=fuellstandVec)) ~ targetVec,
     xlab=paste("Model", colnames(targetObs_STFDF@data)),
     ylab=paste("gemessen", colnames(targetObs_STFDF@data)),
     panel = function(x, y) {
       panel.xyplot(x, y)
       panel.abline(lm(x~y, data.frame(x=20:30, y=20:30)))
     })

if (is.null(niederschlagVec))
  niederschlagVec <- niederschlagPredVec

if (is.null(fuellstandVec))
  fuellstandVec <- fuellstandPredVec

grid <- data.frame(niederschlagVec = rep(seq(min(niederschlagVec, na.rm = T), max(niederschlagVec, na.rm = T),length.out = 20), each=20))
grid$fuellstandVec <- rep(seq(min(fuellstandVec, na.rm = T), max(fuellstandVec, na.rm = T),length.out = 20), 20)
grid$targetVec <- predict(lmMod, grid)

p4 <- wireframe(targetVec ~ niederschlagVec + fuellstandVec, grid, 
                xlab="Niederschlag",
                ylab="Fuellstand",
                zlab=list(colnames(targetObs_STFDF@data), rot=93),
                scales = list(arrows=F),
                # xlim=c(0,15), ylim=c(294.6, 295.1), zlim=c(290, 302),
                panel =  function(x, y, z, ...) {
                  panel.wireframe(x, y, z, ...)
                  panel.cloud(x=df$niederschlagVec, y=df$fuellstandVec, z=df$targetVec, ...)
                })

relations <- "relations.png"
png(file = relations)

print(p1, position=c(0,0.5, 0.5,1), more=T)
print(p2, position=c(0.5,0.5, 1,1), more=T)
print(p4, position=c(0,0,0.5,0.5), more=T)
print(p3, position=c(0.5,0,1,0.5))

graphics.off()

# wps.out: relations, png;

model_prediction <- "model_prediction.csv"
if(is.na(singleInputNiederschlagPred) & is.na(singleInputFuellstandPred)) {
  write.csv(cbind(niederschlagPred$observationData@result$phenomenonTime[!is.na(mIds)],
                  df), file = model_prediction)
} else {
  write.csv(df, file = model_prediction)
}

# wps.out: model_prediction, csv;

##### json

# metaJson <- fromJSON(file="TimeSeriesMetadataSimple.json")
# 
# readLines("TimeSeriesMetadataSimple.json")
# str(readLines("TimeSeriesMetadataSimple.json"))

statLabel <- targetBreakUp[[which(lapply(targetBreakUp, function(x) x[[1]]) == "featureOfInterest")]][2]
rndId <- function() paste("ts", paste(sample(c(0:9,letters[1:6]), 32, replace = T),collapse=""), sep="_")

rndIdInst <- rndId()
targetJsonMeta <- list(id = rndIdInst,
                       label = colnames(targetObs_STFDF@data)[1],
                       station = list(properties = list(id = rndIdInst,
                                                        label = statLabel),
                                      geometry = list(coordinates = coordinates(targetObs_STFDF@sp),
                                                      type="point"),
                                      type = "Feature"))

metaJson <- "metaJson.json"
writeLines(toJSON(targetJsonMeta), metaJson)

# wps.out: metaJson, json; 

bindTimeData <- cbind(predTimes, df$targetVec)
colnames(bindTimeData) <- NULL
bindTimeData <- apply(bindTimeData, 1, function(x) list(timestamp=x[1], value=x[2]))

targetDataJson <- list(id=list(values=bindTimeData))
names(targetDataJson) <- rndIdInst

dataJson <- "dataJson.json"
writeLines(toJSON(targetDataJson), dataJson)

# wps.out: dataJson, json;

# wps.off;
targetSOS <- "https://tamis.dev.52north.org/sos/service"
# wps.on;

## push values into target SOS

if(!is.na(targetSOS)) {
  # insert new sensor
  
  sensor <- list(list("key:uniqueID", paste("http://www.52north.org/test/procedure/linearRegression", 
                                             as.numeric(Sys.time()), sep="")),
                 list("key:longName", "52°North Initiative for Geospatial Open Source Software GmbH (http://52north.org)"),
                 list("key:shortName", "52°North GmbH"),
                 list("key:fieldName", "\"Offering for a linear regression\""),
                 list("key:offeringID:name", "Offering for a linear regression"),
                 list("key:offeringID:value", "http://www.52north.org/test/offering/linearRegression"),
                 list("key:featureOfInterestID", "http://www.52north.org/test/featureOfInterest/linearRegression"),
                 list("key:easting", targetObs_STFDF@sp@coords[1]),
                 list("key:northing", targetObs_STFDF@sp@coords[2]),
#                 list("key:altitude", 60),
                 inputList=list(list(name="linearRegression",
                                    definition="http://www.52north.org/test/observableProperty/linearRegression")),
                outputList=list(list(name="test_observable_property_44_1",
                                     scale="Category",
                                     definition="http://www.52north.org/test/observableProperty/44_1",
                                     codeSpace="<swe:codeSpace xlink:href=\"NOT_DEFINED\"/>"),
                                list(name="test_observable_property_44_2",
                                     scale="Count",
                                     definition="http://www.52north.org/test/observableProperty/44_2"),
                                list(name="test_observable_property_44_3",
                                     scale="Quantity",
                                     definition="http://www.52north.org/test/observableProperty/44_3",
                                     uom="<swe:uom code=\"NOT_DEFINED\"/>"),
                                list(name="test_observable_property_44_4",
                                     scale="Text",
                                     definition="http://www.52north.org/test/observableProperty/44_4"),
                                list(name="test_observable_property_44_5",
                                     scale="Boolean",
                                     definition="http://www.52north.org/test/observableProperty/44_5")),
                observableProperties=list("http://www.52north.org/test/observableProperty/44_1",
                                          "http://www.52north.org/test/observableProperty/44_2",
                                          "http://www.52north.org/test/observableProperty/44_3",
                                          "http://www.52north.org/test/observableProperty/44_4",
                                          "http://www.52north.org/test/observableProperty/44_5"))

                 insSenRet <- insertSensor("https://tamis.dev.52north.org/sos/service", sensor,
                                           template = "inst/templates/InsertSensor.xml", update=FALSE,
                                           add_headers(Authorization=tamis.dev.auth))

                 cat(memDecompress(insSenRet$content, type="none", asChar = T))
  
  
  inputDf <- data.frame(phenomenonTime = as.POSIXct(predTimes, origin = "1970-01-01"),)
  
                        
  
  # testDf <- data.frame(phenomenonTime=Sys.time()+1:10*60,
  #                      test_observable_property_46_3=runif(10))
  # testDf$phenomenonTime <- format(testDf$phenomenonTime, format = "%Y-%m-%dT%H:%M%:%S")
  # 
  # res <- insertMeasurements("https://tamis.dev.52north.org/sos/service",
  #                    ts=testDf,
  #                    coords=c(52,7),
  #                    meta = metaMeasure,
  #                    fieldDefs=list(phenomenonTime=c("<swe:Time definition=\"http://www.opengis.net/def/property/OGC/0/PhenomenonTime\">",
  #                                        "<swe:uom xlink:href=\"http://www.opengis.net/def/uom/ISO-8601/0/Gregorian\"/>",
  #                                        "</swe:Time>"),
  #                                   test_observable_property_46_3=c("<swe:Quantity definition=\"http://www.52north.org/test/observableProperty/46_3\">",
  #                                        "<swe:uom code=\"NOT_DEFINED\"/>",
  #                                        "</swe:Quantity>")),
  #                    srsName = "http://www.opengis.net/def/crs/EPSG/0/4326",
  #                    template = "inst/templates/InsertMeasurement.xml",
  #                    header=add_headers(Authorization=tamis.dev.auth))
  # 
  # cat(memDecompress(res$content, type = "none", asChar = T))
  # 
  # 
  # list(phenomenonTime=c("<ns:Time definition=\"http://www.opengis.net/def/property/OGC/0/PhenomenonTime\">",
  #                       "<ns:uom xlink:href=\"http://www.opengis.net/def/uom/ISO-8601/0/Gregorian\"/>",
  #                       "</ns:Time>"),
  #      Wasserstand_im_Damm=c("<ns:Quantity definition=\"Wasserstand_im_Damm\">",
  #                            "<ns:uom code=\"m\"/>",
  #                            "</ns:Quantity>"))
}
