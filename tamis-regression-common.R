## tamins regression common part

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

precip <- readTSdata(timeseries_Niederschlag, timespan)
fillLevel <- readTSdata(timeseries_Fuellstand, timespan)
targetVar <- readTSdata(ts_URI = timeseries_Zielvariable, timespan, .opts)

precipMeta <- readTSmeta(timeseries_Niederschlag)
fillLevelMeta <- readTSmeta(timeseries_Fuellstand)
targetVarMeta <- readTSmeta(timeseries_Zielvariable, .opts)

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

# target
targetVarAgg <- aggregate(targetVar[1,,drop=F]@data[[1]],
                          by = list(sapply(strsplit(as.character(time(targetVar[1,])),
                                                    " ", fixed=T), function(x) x[1])),
                          function(x) mean(x, na.rm = T))
colnames(targetVarAgg) <- colnames(targetVar@data)
targetVar <- STFDF(targetVar@sp, as.POSIXct(targetVarAgg[,1]), targetVarAgg[,2,drop=F])

joinTimeStamps <- unique(time(precip), time(fillLevel), time(targetVar))

df <- data.frame(time = joinTimeStamps)
df$precip <- precip@data[match(joinTimeStamps, time(precip)),1]
df$fillLevel <- fillLevel@data[match(joinTimeStamps, time(fillLevel)),1]
df$targetVar <- targetVar@data[match(joinTimeStamps, time(targetVar)),1]

# drop missing data
df <- df[apply(df,1, function(x) !any(is.na(x))),]

# modelling

if (nrow(df) < 10) {
  load("preDefTSModel.RData")
  lmMod <- switch(tail(strsplit(timeseries_Zielvariable,"/", fixed=T)[[1]],1),
                  '194' = mod194,
                  stop("No model present, use longer time series to fit one."))
} else {
  lmMod <- lm(targetVar ~ fillLevel + precip, df)
  # summary(lmMod)
  
  # mod194 <- lmMod
  # save(mod194,
  #      file="preDefTSModel.RData")
}

model_diagnostics <- "model_diagnostics.png"
png(file = model_diagnostics)#TODO vielleicht mehr parameter wie groesse etc

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
