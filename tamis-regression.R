# wps.des: tamis-regression, title = TaMIS Regression Model for observedproperty_im_Damm at Bevertalsperre;

# wps.in: observedproperty, string, observed property, 
# abstract = the observed property default: Wasserstand_im_Damm,
# value = Wasserstand_im_Damm;

# wps.in: foi, string, feature of interest, 
# abstract = the feature of interest default: Bever-Talsperre_MQA7_Piezometer_Kalkzone,
# value = Bever-Talsperre_MQA7_Piezometer_Kalkzone;

# "Bever-Talsperre_MQA7_Piezometer_Kalkzone",
# "Bever-Talsperre_MQA1_Piezometer_Wasserseite_Schuettkoerper",
# "Bever-Talsperre_MQA3_Piezometer_Luftseite",
# "Bever-Talsperre_MQA4_Piezometer_Luftseite",
# "Bever-Talsperre_MQA5_Piezometer_Berme"

# "Bever-Talsperre_Sickerwassermessstelle_S2A"
# "Bever-Talsperre_Sickerwassermessstelle_S2B"

# wps.in: phenomenontime, string, phenomenon time, 
# abstract = the phenomenon time;

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

# updateStatus("Requested SOS successfully")

# updateStatus("Requesting observedproperty observations")

# wps.off
# observedproperty <- "Wasserstand_im_Damm"
# foi <- "Bever-Talsperre_MQA7_Piezometer_Kalkzone"
observedproperty <- "Schuettmenge"
foi <- "Bever-Talsperre_Sickerwassermessstelle_S2A"
phenomenontime <- "2016-01-01T00:00/2016-03-10T23:59"
# wps.on

if(observedproperty == "Wasserstand_im_Damm") {
  offs <- "Zeitreihen_Handeingabe"
  proc <- "Handeingabe"
}
if(observedproperty == "Schuettmenge") {
  offs <- "Zeitreihen_Tageswert_Prozessleitsystem"
  proc <- "Tageswert_Prozessleitsystem"
}

targetObs <- getObservation(TaMIS_SOS,
                            offering = offs, 
                            procedure = proc,
                            featureOfInterest=SosFeatureOfInterest(list(foi)),
                            observedProperty = list(observedproperty),
                            eventTime = paste("om:phenomenonTime", phenomenontime, sep=","), 
                            responseFormat = "http://www.opengis.net/om/2.0")

# updateStatus("Requested observedproperty observations successfully")

# updateStatus("Plotting observedproperty")

targetObs_STFDF <- as.STFDF.list.Om_OMObservation(targetObs)

targetObs_plot <- "targetObs_plot.png"
png(file = targetObs_plot)

tmp <- stplot(targetObs_STFDF, mode="ts")
print(tmp)

graphics.off()

# wps.out: targetObs_plot, png;


if(observedproperty == "Wasserstand_im_Damm") {
  times <- as.character(index(targetObs_STFDF@time))
} else {
  times <- as.character(index(targetObs_STFDF@time))
  times <- sapply(strsplit(times, " "), function(x) x[1])
} 

FLUGGS_SOS <- SOS(url = "http://www.fluggs.de/sos2/sos",
                  version = "2.0.0", binding = "KVP")

niederschlag <- vector("list",length = length(times))
for ( i in 1:length(times)) {
niederschlag[[i]] <- getObservation(FLUGGS_SOS,
                                   offering = "Zeitreihen_Tagessumme",
                                   featureOfInterest=SosFeatureOfInterest(list("Bever-Talsperre")),
                                   observedProperty = list("Niederschlagshoehe"),
                                   procedure = "Tagessumme",
                                   eventTime = paste("om:phenomenonTime,", times[i], "T00:00", sep=""),
                                   responseFormat = "http://www.opengis.net/om/2.0")
}

niederschlagVec <-as.numeric(sapply(niederschlag, 
                                    function(x) {
  res <- x[[1]]@result
  if(is.null(res))
    return(NA)
  else 
    return(res)
  }))

if(observedproperty == "Schuettmenge")
  niederschlagVec[niederschlagVec > 100] <- NA

fuellstand <- vector("list", length(times))
for(i in 1:length(times)) {
  fuellstand[[i]] <- getObservation(FLUGGS_SOS,
                               offering = "Zeitreihen_Einzelwert",
                               featureOfInterest=SosFeatureOfInterest(list("Bever-Talsperre_Windenhaus")),
                               observedProperty = list("Speicherfuellstand"),
                               # eventTime = paste("om:phenomenonTime,",times[i],"T00:00/",times[i],"T23:59", sep=""),
                               eventTime = paste("om:phenomenonTime,",times[i],"T10:00",sep=""),
                               responseFormat = "http://www.opengis.net/om/2.0")
}

fuellstandVec <- as.numeric(sapply(fuellstand, function(x) {
  res <- x[[1]]@result
  if(is.null(res))
    return(NA)
  else 
    return(res)
}))

if(observedproperty == "Schuettmenge")
  fuellstandVec[fuellstandVec < 291] <- NA

targetVec <- targetObs_STFDF@data[[1]]

lmMod <- lm(targetVec ~ niederschlagVec + fuellstandVec)
summary(lmMod)

model_diagnostics <- "model_diagnostics.png"
png(file = model_diagnostics)#TODO vielleicht mehr parameter wie groesse etc

par(mfrow=c(2,2))
plot(lmMod)

graphics.off()

# wps.out: model_diagnostics, png;

library(lattice)
p1 <- xyplot(targetVec ~ niederschlagVec,
       xlab="Niederschalg",
       ylab=observedproperty)
p2 <- xyplot(targetVec ~  fuellstandVec,
       xlab="Fuellstand",
       ylab=observedproperty)

p3 <- xyplot(predict(lmMod, data.frame(niederschlagVec=niederschlagVec, fuellstandVec=fuellstandVec)) ~ targetVec,
     xlab=paste("model", observedproperty),
     ylab=paste("gemessener", observedproperty),
     panel = function(x, y) {
       panel.xyplot(x, y)
       panel.abline(lm(x~y, data.frame(x=20:30, y=20:30)))
     })

p4 <- cloud(targetVec ~ niederschlagVec + fuellstandVec,
            xlab="Niederschlag",
            ylab="Fuellstand",
            zlab=list(observedproperty, rot=93))

relations <- "relations.png"
png(file = relations)

print(p1, position=c(0,0.5, 0.5,1), more=T)
print(p2, position=c(0.5,0.5, 1,1), more=T)
print(p4, position=c(0,0,0.5,0.5), more=T)
print(p3, position=c(0.5,0,1,0.5))

graphics.off()

# wps.out: relations, png;