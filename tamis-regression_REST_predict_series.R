# wps.des: tamis-prediction-series, title = TaMIS Prediction Regression Model for Wasserstand_im_Damm or Schuettmenge at Bever-Talsperre;

# wps.in: timeseries_Niederschlag, string, TS URI, 
# abstract = timeseries Id for Niederschlag,
# value = "http://www.fluggs.de/sos2/api/v1/timeseries/427";

# wps.in: timeseries_Fuellstand, string, TS URI, 
# abstract = timeseries URI for Fuellstand,
# value = "http://www.fluggs.de/sos2/api/v1/timeseries/26";

# wps.in: timeseries_Zielvariable, string, TS URI, 
# abstract = timeseries URI for the target variable,
# value = "http://www.fluggs.de/sos2/api/v1/timeseries/194";

# wps.in: timespan, string, timespan of reference period, 
# abstract = timeseries URI for the target variable,
# value = "2016-01-01T/2016-02-29TZ";

# wps.in: timeseries_predict_Niederschlag, string, TS URI, 
# abstract = timeseries URI for prediction values: Niederschlag,
# value = "http://www.fluggs.de/sos2/api/v1/timeseries/427";

# wps.in: timeseries_predict_Fuellstand, string, TS URI, 
# abstract = timeseries URI for prediction values: Fuellstand,
# value = "http://www.fluggs.de/sos2/api/v1/timeseries/26";

# wps.in: timespan_predict, string, prediction time stamp, 
# abstract = timestamp for which the prediction shall be made,
# value = "2016-03-01T/2016-03-31TZ";

#################################

# wps.off;
timeseries_Niederschlag <- "http://www.fluggs.de/sos2/api/v1/timeseries/427"
timeseries_Fuellstand <- "http://www.fluggs.de/sos2/api/v1/timeseries/26"

timeseries_Zielvariable <- "http://www.fluggs.de/sos2/api/v1/timeseries/194"

timespan <- "2016-01-01T/2016-02-29TZ"
#wps.on;

source("tamis-regression-common.R")

# wps.out: model_diagnostics, png;

# wps.out: relations, png;

#############################
## predict target variable ##
#############################


# wps.off;
timeseries_predict_Niederschlag <- "http://www.fluggs.de/sos2/api/v1/timeseries/427"
timeseries_predict_Fuellstand <- "http://www.fluggs.de/sos2/api/v1/timeseries/26"

timespan_predict <- "2016-03-01T/2016-03-31TZ"
#wps.on;


precip <- readTSdata(timeseries_Niederschlag, timespan_predict)
fillLevel <- readTSdata(timeseries_Fuellstand, timespan_predict)
targetVar <- readTSdata(timeseries_Zielvariable, timespan_predict)

precipMeta <- readTSmeta(timeseries_Niederschlag)
fillLevelMeta <- readTSmeta(timeseries_Fuellstand)
targetVarMeta <- readTSmeta(timeseries_Zielvariable)
# synchronise data sets

#precipitation
precipAgg <- aggregate(precip[1,,drop=F]@data[[1]],
                       by = list(sapply(strsplit(as.character(time(precip[1,])),
                                                 " ", fixed=T), function(x) x[1])),
                       function(x) mean(x, na.rm = T))
precipAgg[precipAgg[,2] > 200,2] <- NA
colnames(precipAgg) <- colnames(precip@data)
precip <- STFDF(precip@sp, as.POSIXct(precipAgg[,1]), precipAgg[,2,drop=F])

# fill level
fillLevelAgg <- aggregate(fillLevel[1,,drop=F]@data[[1]],
                          by = list(sapply(strsplit(as.character(time(fillLevel[1,])),
                                                    " ", fixed=T), function(x) x[1])),
                          function(x) mean(x, na.rm = T))
colnames(fillLevelAgg) <- colnames(fillLevel@data)
fillLevel <- STFDF(fillLevel@sp, as.POSIXct(fillLevelAgg[,1]), fillLevelAgg[,2,drop=F])

joinTimeStamps <- unique(time(precip), time(fillLevel))

df <- data.frame(time = joinTimeStamps)
df$precip <- precip@data[match(joinTimeStamps, time(precip)),1]
df$fillLevel <- fillLevel@data[match(joinTimeStamps, time(fillLevel)),1]

# drop missing data
predDf <- df[apply(df,1, function(x) !any(is.na(x))),]

predDf$predVar <- predict(lmMod, predDf)

plot(predVar ~ time, predDf)
plot(precip ~ time, predDf)
plot(fillLevel ~ time, predDf)

model_prediction <- "model_prediction.csv"
write.csv(predDf, file = model_prediction)

# wps.out: model_prediction, csv;

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

targetDataJson <- list(id=list(values=predDf[,c(1,4)]))
names(targetDataJson) <- rndIdInst

dataJson <- "dataJson.json"
writeLines(toJSON(targetDataJson), dataJson)

# wps.out: dataJson, json;

# wps.off;
targetSOS <- "https://tamis.dev.52north.org/sos/service"
# wps.on;

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