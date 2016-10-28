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
# value = "2016-01-01T/2016-02-29TZ";

#################################

# A)
# Ableiten des Wasserstand im Damm über den Sohlenwasserdruck (SWD) aus Niederschlag und Füllhöhe
# evtl zweiter Schritt um den Wasserstand im Damm über Piezometer aus dem Sohlenwasserdruck abzuleiten.

# B)
# Ableiten des Sickerwassers aus Niederschlag und Füllhöhe
# evtl zweiter Schritt um den Wasserstand im Damm über Piezometer aus dem Sohlenwasserdruck abzuleiten.

# Wasserstand im Damm aus den Piezometern (alle 14 Tage und bei "Bedarf")
# http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries?phenomenon=32&offering=277
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

timespan <- "2016-01-01T/2016-09-30TZ"
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

tsIDs <- c(444, # SWD1.1L
           491, # SWD1.2L
           492, # SWD1.3L
           435, # SWD1.4L
           436, # SWD1.5L
           437, # SWD1.6L
           438, # SWD1.7L
           439, # SWD1.8L
           440, # SWD1.9L
           441, # SWD1.10L
           442, # SWD1.11L
           456, # SWD1.12L
           457, # SWD1.13L
           510, # SWD1W
           511, # SWD1.1W
           512, # SWD1L
           517, # SWD2A
           518) # SWD3A

modList <- NULL

pb <- txtProgressBar(0, length(tsIDs), style = 3)
for (ts in tsIDs) {
  timeseriesZielvariable <- paste("http://fluggs.wupperverband.de/sos2-tamis/api/v1/timeseries/",
                                  ts, sep="/")
  
  targetVar <- readTSdata(ts_URI = timeseriesZielvariable, timespan, .opts)
  
  precipMeta <- readTSmeta(timeseriesNiederschlag)
  fillLevelMeta <- readTSmeta(timeseriesFuellstand)
  targetVarMeta <- readTSmeta(timeseriesZielvariable, .opts)

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

  lmMod <- lm(targetVar ~ fillLevel + precip, df)
  
  modList[[as.character(ts)]] <- lmMod
  cat("\n This is timeseries ID:", ts, "\n")
  print(summary(lmMod))
  setTxtProgressBar(pb, which(tsIDs == ts))
}

close(pb)
save(modList, file="preDefTSModel.RData")

## diagnostics
hist(sapply(modList, function(x) summary(x)$adj.r.squared), n=20)

which(sapply(modList, function(x) summary(x)$adj.r.squared) < 0.1)

summary(modList[["518"]])
