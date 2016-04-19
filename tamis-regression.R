# wps.des: tamis-regression, title = TaMIS Regression Model for Wasserstand_im_Damm or Schuettmenge at Bever-Talsperre;

# wps.in: SOSreqNiederschlag, string, SOS-request, 
# abstract = SOS-request for Niederschlag,
# value = http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Niederschlagshoehe&procedure=Tagessumme&featureOfInterest=Bever-Talsperre&&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2015-03-10T13:58:07.519Z%2F2016-03-10T13:58:07.519Z;

# wps.in: SOSreqFuellstand, string, SOS-request, 
# abstract = SOS-request for Fuellstand,
# value = http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Speicherfuellstand&procedure=Einzelwert&featureOfInterest=Bever-Talsperre_Windenhaus&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-03-01T10:00:00.00Z%2F2016-03-10T13:00:00.000Z;

# wps.in: SOSreqTarget, string, SOS-request, 
# abstract = SOS-request for the target variable,
# value = http://fluggs.wupperverband.de/sos2-tamis/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Wasserstand_im_Damm&procedure=Handeingabe&featureOfInterest=Bever-Talsperre_MQA1_Piezometer_Wasserseite_Schuettkoerper&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpati-al%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2016-01-01T10:00:00.00Z%2F2016-03-10T13:00:00.000Z;

# wps.in: SOSreqNiederschlagPred, string, SOS-request, 
# abstract = SOS-request for prediction values: Niederschlag, minOccurs = 0, maxOccurs = 1;

# wps.in: SOSreqFuellstandPred, string, SOS-request, 
# abstract = SOS-request for prediction values: Fuellstand, minOccurs = 0, maxOccurs = 1;

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

# wps.off;
SOSreqNiederschlag <- "http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Niederschlagshoehe&procedure=Tagessumme&featureOfInterest=Bever-Talsperre&&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2015-10-01T13:58:07.519Z%2F2016-03-10T13:58:07.519Z"

SOSreqFuellstand <- "http://www.fluggs.de/sos2/sos?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Speicherfuellstand&procedure=Einzelwert&featureOfInterest=Bever-Talsperre_Windenhaus&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpatial%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2015-10-01T10:00:00.00Z%2F2016-03-10T13:00:00.000Z"

SOSreqTarget <- "http://fluggs.wupperverband.de/sos2-tamis/service?service=SOS&version=2.0.0&request=GetObservation&responseformat=http://www.opengis.net/om/2.0&observedProperty=Wasserstand_im_Damm&procedure=Handeingabe&featureOfInterest=Bever-Talsperre_MQA1_Piezometer_Wasserseite_Schuettkoerper&namespaces=xmlns%28sams%2Chttp%3A%2F%2Fwww.opengis.net%2FsamplingSpati-al%2F2.0%29%2Cxmlns%28om%2Chttp%3A%2F%2Fwww.opengis.net%2Fom%2F2.0%29&temporalFilter=om%3AphenomenonTime%2C2015-10-01T10:00:00.00Z%2F2016-03-10T13:00:00.000Z"

SOSreqNiederschlagPred <- NA

SOSreqFuellstandPred <- NA
# wps.on;

targetBreakUp <- strsplit(SOSreqTarget,split = "?", fixed = T)[[1]]
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

fuellstandBreakUp <- strsplit(SOSreqFuellstand,split = "?", fixed = T)[[1]]

fuellstandURL <- fuellstandBreakUp[1] 

fuellstandBreakUp <- lapply(strsplit(fuellstandBreakUp[2], "&", fixed=T)[[1]], function(x) strsplit(x, "=", fixed=T)[[1]])

fuellstandVersion <- fuellstandBreakUp[[match("version", sapply(fuellstandBreakUp, function(x) x[1]))]][2]

FLUGGS_SOS <- SOS(url = fuellstandURL,
                  version = fuellstandVersion, binding = "KVP")

parList <- SOSreqBreakup(fuellstandBreakUp)

parList$sos <- FLUGGS_SOS
parList$offering <- "Zeitreihen_Einzelwert"

if(is.na(SOSreqFuellstandPred)) {
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

fuellstandVec[fuellstandVec < 291] <- NA

# check for prediction
if (!is.na(SOSreqFuellstandPred)) {
  fuellstandBreakUp <- strsplit(SOSreqFuellstandPred,split = "?", fixed = T)[[1]]
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

### Niederschlag
niederschlagBreakUp <- strsplit(SOSreqNiederschlag,split = "?", fixed = T)[[1]]

niederschlagURL <- niederschlagBreakUp[1] 

niederschlagBreakUp <- lapply(strsplit(niederschlagBreakUp[2], "&", fixed=T)[[1]], function(x) strsplit(x, "=", fixed=T)[[1]])

niederschlagVersion <- niederschlagBreakUp[[match("version", sapply(niederschlagBreakUp, function(x) x[1]))]][2]

FLUGGS_SOS <- SOS(url = niederschlagURL,
                  version = niederschlagVersion, binding = "KVP")

parList <- SOSreqBreakup(niederschlagBreakUp)

parList$sos <- FLUGGS_SOS
parList$offering <- "Zeitreihen_Tagessumme"

if (is.na(SOSreqNiederschlagPred)) {
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

niederschlagVec[niederschlagVec > 100] <- NA

if (!is.na(SOSreqNiederschlagPred)) {
  niederschlagBreakUp <- strsplit(SOSreqNiederschlagPred,split = "?", fixed = T)[[1]]
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

# match covariates
mIds <- match(niederschlagPred$observationData@result$phenomenonTime, 
              fuellstandPred$observationData@result$phenomenonTime)

fuellstandPredVec <- as.numeric(fuellstandPred$observationData@result$Speicherfuellstand)
fuellstandPredVec <- fuellstandPredVec[mIds[!is.na(mIds)]]
fuellstandPredVec[fuellstandPredVec < 291] <- NA

niederschlagPredVec <- as.numeric(niederschlagPred$observationData@result$Niederschlagshoehe)
niederschlagPredVec <- niederschlagPredVec[!is.na(mIds)]
niederschlagPredVec[niederschlagPredVec > 100] <- NA

# modelling

targetVec <- targetObs_STFDF@data[[1]]

lmMod <- lm(targetVec ~ niederschlagVec + fuellstandVec)
summary(lmMod)

model_diagnostics <- "model_diagnostics.png"
png(file = model_diagnostics)#TODO vielleicht mehr parameter wie groesse etc

par(mfrow=c(2,2))
plot(lmMod)

graphics.off()
# 
# if (!SOSreqFuellstandPred == "NULL" & !SOSreqNiederschlagPred == "NULL") {
#   df <- data.frame(niederschlagVec=niederschlagPredVec, fuellstandVec=fuellstandPredVec)
# } else {
  df <- data.frame(niederschlagVec=niederschlagPredVec, fuellstandVec=fuellstandPredVec)
# }

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
write.csv(cbind(niederschlagPred$observationData@result$phenomenonTime[!is.na(mIds)],
                df), file = model_prediction)

# wps.out: model_prediction, csv;