# define new grid

library(rgdal)
target <- readGDAL("geotiff.tiff")

str(target)

target <- spTransform(target, CRS("+init=epsg:4326"))

plot(target)

round(apply(target@bbox, 1, diff)*1e4, 0)

apply(target@bbox, 1, mean)

newGRid <- GridTopology(c(7.37-70*0.5e-4+0.25e-4, 51.14-0.5e-4), cellsize = c(1e-4,1e-4)/2, cells.dim = c(140,66))

newSpGrid <- SpatialGrid(newGRid, "+init=epsg:4326")

plot(newSpGrid)
points(target)

str(target)

hist(target$band1)

newTarget <- SpatialGridDataFrame(newSpGrid, data.frame(data=runif(140*66)))

spplot(newTarget)

writeGDAL(newTarget, "geotiff.tiff", drivername="GTiff")