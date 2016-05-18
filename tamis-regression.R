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

# wps.off;
sosInputNiederschlag <- "http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Niederschlagshoehe&procedure=Tagessumme&featureOfInterest=Bever-Talsperre&&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T10:00:00.00Z%2F2016-04-30T23:59:00.000Z"

sosInputFuellstand <- "http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Speicherfuellstand&procedure=Einzelwert&featureOfInterest=Bever-Talsperre_Windenhaus&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T10:00:00.00Z%2F2016-04-30T23:59:00.000Z"

sosInputTarget <- "http://fluggs.wupperverband.de/sos2-tamis/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Wasserstand_im_Damm&procedure=Handeingabe&featureOfInterest=Bever-Talsperre_MQA7_Piezometer_Kalkzone&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpati-al%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T00:01:00.00Z%2F2016-04-30T23:59:00.000Z"
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

source("~/52North/secOpts.R")
TaMIS_SOS <- SOS(url = targetURL,
                 version = targetVersion, binding = "KVP", curlOptions = .opts)

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

# updateStatus("Requested SOS successfully")

# updateStatus("Requesting observedproperty observations")

parList <- SOSreqBreakup(targetBreakUp)

if(parList$observedProperty == "Wasserstand_im_Damm") {
  parList$offering <- "Zeitreihen_Handeingabe"
}
if(parList$observedProperty == "Schuettmenge") {
  parList$offering <- "Zeitreihen_Tageswert_Prozessleitsystem"
}

parList$sos <- TaMIS_SOS

targetObs <- do.call(getObservation, parList)

# updateStatus("Requested observedproperty observations successfully")

# updateStatus("Plotting observedproperty")

targetObs_STFDF <- as.STFDF.list.Om_OMObservation(targetObs)

targetObs_plot <- "targetObs_plot.png"
png(file = targetObs_plot)

tmp <- stplot(targetObs_STFDF, mode="ts")
print(tmp)

graphics.off()

# wps.out: targetObs_plot, png;

if(parList$observedProperty == "Wasserstand_im_Damm") {
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

if (length(times) < 10 | any(c(is.na(sosInputFuellstand),
                              is.na(sosInputNiederschlag)))) {
  load("preDefModel.RData")
  lmMod <- switch(targetBreakUp[[which(lapply(targetBreakUp, function(x) x[[1]]) == "featureOfInterest")]][2],
                  'Bever-Talsperre_MQA1_Piezometer_Wasserseite_Schuettkoerper' = MQA1mod,
                  'Bever-Talsperre_MQA3_Piezometer_Luftseite' = MQA3mod,
                  'Bever-Talsperre_MQA4_Piezometer_Luftseite' = MQA4mod,
                  'Bever-Talsperre_MQA5_Piezometer_Berme' = MQA5mod,
                  'Bever-Talsperre_MQA7_Piezometer_Kalkzone' = MQA7mod,
                  'Bever-Talsperre_Sickerwassermessstelle_S2A' = S2Amod,
                  'Bever-Talsperre_Sickerwassermessstelle_S2B' = S2Bmod)
} else {
  lmMod <- lm(targetVec ~ fuellstandVec + niederschlagVec)
  # summary(lmMod)
  
  # MQA1mod <- lmMod
  # MQA3mod <- lmMod
  # MQA4mod <- lmMod
  # MQA5mod <- lmMod
  # MQA7mod <- lmMod
  
  # S2Amod <- lmMod
  # S2Bmod <- lmMod
  # save(MQA1mod, MQA3mod, MQA4mod, MQA5mod, MQA7mod,
  #      S2Amod, S2Bmod,
  #      file="preDefModel.RData")
}

model_diagnostics <- "model_diagnostics.png"
png(file = model_diagnostics)#TODO vielleicht mehr parameter wie groesse etc

par(mfrow=c(2,2))
plot(lmMod)

graphics.off()

df <- data.frame(niederschlagVec=niederschlagPredVec, fuellstandVec=fuellstandPredVec)

df$targetVec <- predict(lmMod, df)
  
# wps.out: model_diagnostics, png;

library(lattice)
p1 <- xyplot(df$targetVec ~ df$niederschlagVec,
       xlab="Niederschlag",
       ylab=colnames(targetObs_STFDF@data))
p2 <- xyplot(df$targetVec ~  df$fuellstandVec, 
       xlab="Fuellstand",
       ylab=colnames(targetObs_STFDF@data))

p3 <- xyplot(predict(lmMod, niederschlagVec=niederschlagVec, fuellstandVec=fuellstandVec) ~ targetVec,
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
                panel =  function(x, y, z, ...) {
                  panel.wireframe(x, y, z, ...)
                  panel.cloud(df$niederschlagVec, df$fuellstandVec, df$targetVec, ...)
                },
                scales = list(arrows=F))

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

# fromJSON(file=metaJson)

# dataJson <- fromJSON(file="TimeSeriesRawData.json")
# str(dataJson)
# 
# readLines("TimeSeriesRawData.json")
# str(readLines("TimeSeriesRawData.json"))
# 
# dataJson[[1]]
# 
# str(df)

bindTimeData <- cbind(predTimes, df$targetVec)
colnames(bindTimeData) <- NULL
bindTimeData <- apply(bindTimeData, 1, function(x) list(timestamp=x[1], value=x[2]))

targetDataJson <- list(id=list(values=bindTimeData))
names(targetDataJson) <- rndIdInst

dataJson <- "dataJson.json"
writeLines(toJSON(targetDataJson), dataJson)

# wps.out: model_prediction, csv;

# fromJSON(file=dataJson)
