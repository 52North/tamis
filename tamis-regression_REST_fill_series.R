# wps.des: tamis-rest-fill-series, title = TaMIS Prediction Regression Model for Wasserstand_im_Damm or Schuettmenge at Bever-Talsperre;

# wps.in: timeseriesNiederschlag, string, TS URI, 
# abstract = timeseries Id for Niederschlag,
# value = "http://www.fluggs.de/sos2/api/v1/timeseries/427";

# wps.in: timeseriesFuellstand, string, TS URI, 
# abstract = timeseries URI for Fuellstand,
# value = "http://www.fluggs.de/sos2/api/v1/timeseries/26";

# wps.in: timeseriesZielvariable, string, TS URI, 
# abstract = timeseries URI for the target variable,
# value = "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/451";

# wps.in: timespan, string, timespan of reference period, 
# abstract = timeseries URI for the target variable,
# value = "2016-11-11T/2016-11-16TZ";

#################################

# A)
# Ableiten des Wasserstand im Damm über den Sohlenwasserdruck (SWD) aus Niederschlag und Füllhöhe
# evtl zweiter Schritt um den Wasserstand im Damm über Piezometer aus dem Sohlenwasserdruck abzuleiten.

# B)
# Ableiten des Sickerwassers aus Niederschlag und Füllhöhe
# evtl zweiter Schritt um den Wasserstand im Damm über Piezometer aus dem Sohlenwasserdruck abzuleiten.

# Wasserstand im Damm aus den Piezometern (abgeleitet; alle 14 Tage und bei "Bedarf")
# http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries?phenomenon=32&offering=277
#
# Wasserstand im Damm aus den Piezometern (Handeingabe der "tiefe"; alle 14 Tage und bei "Bedarf")
# http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries?phenomenon=32&offering=269
# 
# Wasserstand im Damm aus dem Sohlenwasserdruck (täglich)
# http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries?phenomenon=35&offering=272
# 
# Das Sickerwasser gemssen als Schüttmenge in l/s
# http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries?phenomenon=34&offering=272

# Hier: Verdichtung der Zeitreihen über Niederschlag und Füllhöhe 

# wps.off;
timeseriesNiederschlag <- "http://www.fluggs.de/sos2/api/v1/timeseries/427"
timeseriesFuellstand <- "http://www.fluggs.de/sos2/api/v1/timeseries/26"

timeseriesZielvariable <- "http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/451"

timespan <- "2016-11-11T/2016-11-16TZ"
# wps.on;

library(RCurl)
library(httr)
library(rjson)
library(sp)
library(xts)
library(spacetime)

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

source("~/52North/secOpts.R")

precip <- readTSdata(timeseriesNiederschlag, timespan)
fillLevel <- readTSdata(timeseriesFuellstand, timespan)
targetVar <- readTSdata(ts_URI = timeseriesZielvariable, timespan, .opts)

as.POSIXct(1477263600, origin="1970-01-01 00:00:00")

as.numeric(Sys.Date())

debugonce(readTSdata)

precipMeta <- readTSmeta(timeseriesNiederschlag)
fillLevelMeta <- readTSmeta(timeseriesFuellstand)
targetVarMeta <- readTSmeta(timeseriesZielvariable, .opts)

# synchronise data sets
# 
# #precipitation
# precipAgg <- aggregate(precip[1,,drop=F]@data[[1]],
#                        by = list(sapply(strsplit(as.character(time(precip[1,])),
#                                                  " ", fixed=T), function(x) x[1])),
#                        function(x) mean(x, na.rm = T))
# precipAgg[precipAgg[,2] > 200,2] <- NA
# colnames(precipAgg) <- colnames(precip@data)
# precip <- STFDF(precip@sp, as.POSIXct(precipAgg[,1]), precipAgg[,2,drop=F])
# 
# # fill level
# fillLevelAgg <- aggregate(fillLevel[1,,drop=F]@data[[1]],
#                           by = list(sapply(strsplit(as.character(time(fillLevel[1,])),
#                                                     " ", fixed=T), function(x) x[1])),
#                           function(x) mean(x, na.rm = T))
# colnames(fillLevelAgg) <- colnames(fillLevel@data)
# fillLevel <- STFDF(fillLevel@sp, as.POSIXct(fillLevelAgg[,1]), fillLevelAgg[,2,drop=F])
# 
# target
targetVarAgg <- aggregate(targetVar[1,,drop=F]@data[[1]],
                          by = list(sapply(strsplit(as.character(time(targetVar[1,])),
                                                    " ", fixed=T), function(x) x[1])),
                          function(x) mean(x, na.rm = T))
colnames(targetVarAgg) <- colnames(targetVar@data)
targetVar <- STFDF(targetVar@sp, as.POSIXct(targetVarAgg[,1], tz = "CET"), targetVarAgg[,2,drop=F])

joinTimeStamps <- unique(time(precip), time(fillLevel))

df <- data.frame(time = joinTimeStamps)
df$precip <- precip@data[match(joinTimeStamps, time(precip)),1]
df$fillLevel <- fillLevel@data[match(joinTimeStamps, time(fillLevel)),1]
df$targetVar <- targetVar@data[match(joinTimeStamps, time(targetVar)),1]


# drop missing data
df <- df[apply(df,1, function(x) !any(is.na(x[c("precip", "fillLevel")]))),]

# modelling

# wps.res: /tmp/gitrepositories/preDefTSModel.RData;

if (sum(!is.na(df$targetVar)) < 10) {
  # load("preDefTSModel.RData")
  lmMod <- modList[[tail(strsplit(timeseriesZielvariable,"/", fixed=T)[[1]],1)]]
} else {
  lmMod <- lm(targetVar ~ fillLevel + precip, df)
}

modelDiagnostics <- "modelDiagnostics.png"
png(file = modelDiagnostics)#TODO vielleicht mehr parameter wie groesse etc

par(mfrow=c(2,2))
plot(lmMod)

graphics.off()

df$predVar <- predict(lmMod, df)

library(lattice)
p1 <- xyplot(df$targetVar ~ df$precip,
             xlab = precipMeta$parameters$phenomenon$label,
             ylab = targetVarMeta$parameters$phenomenon$label)
p2 <- xyplot(df$targetVar ~  df$fillLevel, 
             xlab = fillLevelMeta$parameters$phenomenon$label,
             ylab = targetVarMeta$parameters$phenomenon$label)

p3 <- xyplot(predVar ~ targetVar, df,
             xlab="Modell",
             ylab="Beobachtung",
             main=targetVarMeta$parameters$phenomenon$label,
             panel = function(x, y) {
               panel.xyplot(x, y)
               panel.abline(lm(x~y, data.frame(x=20:30, y=20:30)))
             })

grid <- data.frame(precip = rep(seq(min(df$precip, na.rm = T), 
                                    max(df$precip, na.rm = T), 
                                    length.out = 20), each=20))
grid$fillLevel <- rep(seq(min(df$fillLevel, na.rm = T), 
                          max(df$fillLevel, na.rm = T),
                          length.out = 20), 20)
grid$predVar <- predict(lmMod, grid)

p4 <- wireframe(predVar ~ precip + fillLevel, grid, 
                xlab=precipMeta$parameters$phenomenon$label,
                ylab=fillLevelMeta$parameters$phenomenon$label,
                zlab=list(targetVarMeta$parameters$phenomenon$label, rot=93),
                scales = list(arrows=F),
                panel =  function(x, y, z, ...) {
                  panel.wireframe(x, y, z, ...)
                  panel.cloud(x=df$precip, y=df$fillLevel, z=df$targetVar, ...)
                })

relations <- "relations.png"
png(file = relations)

print(p1, position=c(0,0.5, 0.5,1), more=T)
print(p2, position=c(0.5,0.5, 1,1), more=T)
print(p4, position=c(0,0,0.5,0.5), more=T)
print(p3, position=c(0.5,0,1,0.5))

graphics.off()

df$res <- df$targetVar
df$res[is.na(df$res)] <- df$predVar[is.na(df$res)]


# wps.out: modelDiagnostics, png;

# wps.out: relations, png;

modelPrediction <- "modelPrediction.csv"
write.csv(df[,c(1,2,3,6)], file = modelPrediction)

# wps.out: modelPrediction, csv;

##### json

# metaJson <- fromJSON(file="TimeSeriesMetadataSimple.json")
# 
# readLines("TimeSeriesMetadataSimple.json")
# str(readLines("TimeSeriesMetadataSimple.json"))

statLabel <- targetVarMeta$station$properties$label
rndId <- function() paste("ts", paste(sample(c(0:9,letters[1:6]), 32, replace = T), collapse=""), sep="_")

rndIdInst <- rndId()
targetJsonMeta <- list(id = rndIdInst,
                       label = targetVarMeta$parameters$phenomenon$label,
                       station = list(properties = list(id = rndIdInst,
                                                        label = statLabel),
                                      geometry = list(coordinates = coordinates(targetVar@sp),
                                                      type="point"),
                                      type = "Feature"))

metaJson <- "metaJson.json"
writeLines(toJSON(targetJsonMeta), metaJson)

# wps.out: metaJson, json; 

targetDataJson <- list(id=list(values=df[,c(1,6)]))
names(targetDataJson) <- rndIdInst

dataJson <- "dataJson.json"
writeLines(toJSON(targetDataJson), dataJson)

# wps.out: dataJson, json;

# # wps.off;
# targetSOS <- "https://tamis.dev.52north.org/sos/service"
# # wps.on;
# 
## push values into target SOS
# 
# if(!is.na(targetSOS)) {
#   nowSecs <- round(as.numeric(Sys.time()), 0)
#   
#   targetFoi <- targetVarMeta$parameters$phenomenon$label
#   targetObsPropURI <- targetVarMeta$parameters$phenomenon$id
#   targetObsProp <- tail(strsplit(targetObsPropURI, "/")[[1]], 1)
# 
#   foiRet <- getFeatureOfInterest(TaMIS_SOS, featureOfInterest = targetFoi)
#   
#   targetCoords <- as.numeric(strsplit(foiRet$featureMember@feature@shape@point@pos@pos, " ")[[1]])
#   
#   # check whether appropriate sensor exists in sos, if not, insert new sensor "linearRegression_obsProp"
#   # then push data
#   
#   push_SOS <- SOS(url = targetSOS,
#                    version = "2.0.0", binding = "KVP")
#   
#   push_Cap <- getCapabilities(push_SOS)
#   
#   push_off <- paste("linearRegression", targetObsProp, sep="_")
#   
#   if (length(grep(push_off, names(push_Cap@contents@observationOfferings))) == 0) {
#     sensor <- list(list("key:uniqueID", paste("http://www.52north.org/test/procedure/linearRegression_b", 
#                                               targetObsProp, sep="_")),
#                    list("key:longName", "52°North Initiative for Geospatial Open Source Software GmbH (http://52north.org)"),
#                    list("key:shortName", "52°North GmbH"),
#                    list("key:fieldName", "\"Offering of linear regression WPS\""),
#                    list("key:offeringID:name", "Offering of linear regression WPS"),
#                    list("key:offeringID:value", paste("http://www.52north.org/test/offering/linearRegression_b",
#                                                       targetObsProp, sep="_")),
#                    list("key:featureOfInterestID", paste("http://www.52north.org/test/featureOfInterest/",
#                                                          targetFoi, sep="")),
#                    list("key:northing", targetCoords[1]),
#                    list("key:easting",  targetCoords[2]),
#                    inputList=list(list(name=targetObsProp,
#                                        definition=targetObsPropURI)),
#                    outputList=list(list(name=targetObsProp,
#                                         scale="Quantity",
#                                         definition=targetObsPropURI,
#                                         uom="<swe:uom code=\"NOT_DEFINED\"/>")),
#                    observableProperties=list(targetObsPropURI))
#     
#     insSenRet <- insertSensor("https://tamis.dev.52north.org/sos/service", sensor,
#                               update=FALSE,
#                               add_headers(Authorization=tamis.dev.auth),
#                               template="InsertSensorLinearRegression.xml")
#   }
#   
#   inputDf <- data.frame(phenomenonTime = format(as.POSIXct(predTimes, 
#                                                            origin = "1970-01-01"), 
#                                                 format="%Y-%m-%dT%H:%M:%S"),
#                         targetVar = df[,3])
#   colnames(inputDf) <- c("phenomenonTime", targetObsProp)
#   
#   fieldDefs <- list(phenomenonTime=c("<swe:Time definition=\"http://www.opengis.net/def/property/OGC/0/PhenomenonTime\">",
#                                      "<swe:uom xlink:href=\"http://www.opengis.net/def/uom/ISO-8601/0/Gregorian\"/>",
#                                      "</swe:Time>"))
#   fieldDefs[[targetObsProp]] <- c(paste("<swe:Quantity definition=\"", targetObsPropURI, "\">", sep=""),
#                                   "<swe:uom code=\"NOT_DEFINED\"/>",
#                                   "</swe:Quantity>")
#   
#   proj4split <- sapply(strsplit(targetObs_STFDF@sp@proj4string@projargs, c(" "))[[1]],
#                        function(x) strsplit(x, "="))
#   epsgID <- which(sapply(proj4split, function(x) "+init" %in% x))
#   epsgCode <- strsplit(proj4split[[epsgID]][2], ":")[[1]][2]
#   
#   metaMeasure <- list(list("key:offering", paste("http://www.52north.org/test/offering/linearRegression_b",
#                                                  targetObsProp, sep="_")),
#                       list("key:description", "predictions based on a linear regression"),
#                       list("key:identifier", paste("http://www.52north.org/test/observation/",
#                                                    nowSecs, sep="")),
#                       list("key:procedure", paste("\"http://www.52north.org/test/procedure/linearRegression_b_",
#                                                   targetObsProp,"\"", sep="")),
#                       list("key:obsProp", paste("\"",targetObsProp,"\"", sep="")),
#                       list("key:wml:gml:id", paste("\"", foiRet$featureMember@feature@id, "\"", sep="")),
#                       list("key:resultTime", as.character(format(as.POSIXct(nowSecs, origin="1970-01-01"),
#                                                                      format="%Y-%m-%dT%H:%M:%SZ"))),
#                       list("key:gml:identifier", paste("http://www.52north.org/test/featureOfInterest/",
#                                                       targetFoi, sep="")),
#                       list("key:sf", paste("\"http://www.52north.org/test/featureOfInterest/sf_",
#                                            targetFoi, "\"", sep="")),
#                       list("key:point:gml:id",paste("\"", targetFoi, "\"", sep="")))
#   
#   insMeasRet <- insertMeasurements(targetSOS, coords = targetCoords,
#                                  ts=inputDf, 
#                                    meta = metaMeasure,
#                                    fieldDefs = fieldDefs,
#                                    srsName = "http://www.opengis.net/def/crs/EPSG/0/31466",
#                                    template = "InsertMeasurementLinearRegression.xml",
#                                    header=add_headers(Authorization=tamis.dev.auth))
# 
#   cat(memDecompress(insMeasRet$content, type = "none", asChar = T))
# }
