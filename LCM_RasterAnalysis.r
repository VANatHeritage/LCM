## This code applies "Landscape Condition" methods developed 
## by Colorado NHP/Co State (See SHRP 2 C21A April 2012). This 2019 revision
## (by David Bucklin, Virginia NHP) adds methods to help with large extent/high 
## resolution analyses at the state level.
## Aissa Feldmann, modifying Tim Howard's code from June 2012. Date 21 May 2013, NYNHP.

# Make sure input rasters, weights, etc, are assigned in rules.xlsx
library(here)
library(raster)
library(snow)
library(arcgisbinding)
library(sf)
library(readxl)
library(fasterize)
arc.check_product()
delete.temp <- TRUE # delete intermediate files after completion?

dt <- gsub("-","",Sys.Date())
rasterOptions(tmpdir = "L:/scratch/raster") # set this to a folder with enough disk space for temp rasters
removeTmpFiles(1)
rules1 <- read_excel(here("rules.xlsx"))

#set a feature extent/mask
ext <- st_read(here("inputs","LCM_domain.shp"))

# location of distance rasters. If geodatabase, make sure to add extension
ras.dir <- "L:/David/projects/LCM_temp/LCM_data_full_10m_20190829.gdb"
proj.name <- gsub(".gdb$", "", basename(ras.dir))
if (grepl(".gdb$", ras.dir)) ras <- arc.open(ras.dir)@children$RasterDataset else ras <- list.files(ras.dir, pattern = ".tif$", full.names = F)

# create output directory
out.dir <- here("outputs", paste0(proj.name))#, "_run", dt))
dir.create(out.dir, showWarnings = F)
setwd(out.dir)

# make rules table for this analysis
rules <- data.frame(filename_ext = ras, filename = gsub(".tif","",ras))
rules <- merge(rules, rules1, by = "filename")
rules$full_path <- paste0(ras.dir, "/", rules$filename_ext)
rules$filename_ext <- NULL
print(rules$filename)
write.csv(rules, file = paste0("rules.csv"), row.names = F)

blksz <- 9000 # side-length of raster processing blocks, in pixels. Optimized for using 9 cores with 32 GB memory

obj <- arc.raster(arc.open(rules$full_path[1]))
rext <- obj$extent
sf <- st_sf(id = 1, geom = st_sfc(st_polygon(list(rbind(c(rext[1],rext[2]), c(rext[1],rext[4]), c(rext[3],rext[4]),
                                    c(rext[3],rext[2]), c(rext[1], rext[2]))))))
st_crs(sf) <- st_crs(wkt = obj$sr$WKT)
blocks <- st_sf(geom = st_make_grid(sf, blksz*obj$cellsize[1]))
blocks <- suppressWarnings(st_intersection(blocks, sf))
ext <- st_transform(ext, st_crs(blocks))

blocks <- blocks[unlist(lapply(st_intersects(blocks, ext), any)),]
blocks$start_row <- NA
blocks$start_col <- NA
blocks$end_row <- NA
blocks$end_col <- NA
blocks$id <- 1:length(blocks$geom)

for (i in 1:nrow(blocks)) {
  bb <- st_bbox(blocks[i,])
  blocks$start_row[i] <- (abs(bb[4] - obj$extent[4])) / obj$cellsize[2]
  blocks$start_col[i] <- (bb[1] - obj$extent[1]) / obj$cellsize[1]
  blocks$end_row[i] <- (abs(bb[2] - obj$extent[4])) / obj$cellsize[2]
  blocks$end_col[i] <- (bb[3] - obj$extent[1]) / obj$cellsize[1]
}
plot(blocks["id"], reset = FALSE, main = "processing blocks")
plot(ext$geometry, add = T, col = NA, reset = T)
st_write(blocks, "blocks.shp", delete_layer = T)
suppressWarnings(st_write(ext, "ext.shp", delete_layer = T))
dir.create("rawImpact")

# process blocks
print(system.time({
for (blk in 1:nrow(blocks)) {

  block <- blocks[blk,]
  sub <- c(block$start_row, block$end_row, block$start_col, block$end_col)
  blkrast <- paste0(tempfile(), ".tif")
  tempr = arc.raster(NULL, path=blkrast, extent=st_bbox(block), 
                 nrow=block$end_row-block$start_row, ncol=block$end_col-block$start_col, 
                 nband=1, pixel_type="F32", sr=obj$sr, overwrite=TRUE)
  blkid <- block$id

  # arc method w/ chunk processing
  if (blk == 1) {
    tdir <- paste0(dirname(ras.dir), "/", proj.name, "_temp")
    dir.create(tdir, showWarnings = F)
    
    # cluster apply, should be faster when many cores available
    # using too many cores + large blocks can bog down system due to high memory usage
    cl <- makeCluster(min(9, nrow(rules)), type = "SOCK")  # max arc clusters calls appears to be 9.
    
    clusterEvalQ(cl, library(arcgisbinding))
    clusterCall(cl, arc.check_product) # takes a while, can error if too many cores (maybe due to an Arc license issue)
    clusterExport(cl, c("ras.dir","tdir"))
  }
  clusterExport(cl, c("tempr","sub", "blkid"))
  
  # process raster values
  message("Creating impact rasters for block ", blkid, "...")
  out <- parApply(cl, rules, 1, FUN = function(x) {
    flnm <- x["filename"]
    shift <- as.numeric(x["a"])
    spread <- as.numeric(x["b"])
    dist <- as.numeric(x["decdist"])
    scal <- as.numeric(x["scalar"])
    wt <- as.numeric(x["weight"])
    r <- arc.raster(arc.open(x["full_path"]))
    file <- paste0(tdir, "/", flnm , "_calc", blkid, ".tif")
    r2 = arc.raster(NULL, path = file, dim=dim(tempr), pixel_type="F32", 
                    nodata=tempr$nodata, extent=tempr$extent, sr=tempr$sr, overwrite = T)
      
    v <- r$pixel_block(ul_y = sub[1], nrow = sub[2]-sub[1], ul_x = sub[3], ncol = sub[4] - sub[3])
      
    # v1 decdist function  
    v[v > dist] <- NA
    v[!is.na(v)] <- (1/(1+(exp(((v[!is.na(v)]/scal)-shift)*spread))))*wt
    
    r2$write_pixel_block(v, nrow = sub[2]-sub[1])
    r2$commit(opt = c("build-pyramid"))
    return(file)
  })
  arc.delete(blkrast)
  
  stk <- stack(out)
  stk <- stackSave(stk, paste0("ImpactRasters", blkid, ".stk"))

  # sum impact rasters
  # this is kinda slow
  message("Summing impact rasters for block ", blkid, "...")
  fun <- function(x) {
    if (all(is.na(x))) return(0) else return(sum(x, na.rm=T))
  }
  # uses cluster from above
  rasSum <- clusterR(stk, calc, args = list(fun = fun), cl = cl, 
                     filename = paste0("rawImpact/sumImpact", blkid, ".tif"), overwrite = TRUE)
  if (blk == nrow(blocks)) stopCluster(cl)
}
}))

if (delete.temp) {
  unlink(tdir, recursive = T, force = T)
  unlink(list.files(pattern=".stk$"))
}

#############
# restart from here if already finished previous steps, setting project folder below
#############

library(here)
delete.temp <- TRUE # delete intermediate files after completion?
# proj.name <- "LCM_data_1mlrg_1m_20190729"
setwd(here("outputs", proj.name))
blocks <- st_read("blocks.shp")
ext <- st_read("ext.shp")

#  process sum rasters
print(system.time({
  
lr <- list.files(path = "rawImpact", pattern = "^sumImpact.*tif$", full.names = T)
rlist <- lapply(lr, raster)
cl <- makeCluster(min(length(rlist), parallel::detectCores()-1), type = "SOCK") 
clusterEvalQ(cl, library(raster))
clusterExport(cl, c("rlist"))
rmx <- parLapply(cl, rlist, fun = function(x) cellStats(x, stat = "max", na.rm = T))
mx <- max(unlist(rmx), na.rm=T)
mx2 <- max(rules$weight) 

clusterEvalQ(cl, library(fasterize))
clusterEvalQ(cl, library(sf))
clusterExport(cl, c("mx","mx2","blocks","ext"))

# rescale; divide by overall maximum value (absmax), and the maximum of the stressor weights (rulesmax). 
# Rounds to 3 decimals and multiplies by 1000 to save filesize
# also masks (only if necessary), using 'ext' polygon
dir.create("final")
rasSumRescl <- parLapply(cl, rlist, fun = function(x) {
  nm <- names(x)
  blkid <- as.integer(gsub("\\D", "", nm))
  if (st_area(st_intersection(ext, blocks[blocks$id == blkid,])) != st_area(blocks[blocks$id == blkid,])) {
    # mask
    m1 <- fasterize(st_intersection(ext, blocks[blocks$id == blkid,]), x) 
    out <- mask(calc(x, fun = function(x) round((1 - (x / mx))*10000)), m1, filename = paste0("final/", names(x), "_resc_absmax.tif"), datatype = "INT2U")
    out2 <- mask(calc(x, fun = function(x) {
      x[x>mx2] <- mx2
      round((1 - (x / mx2))*10000)
      }), m1, filename = paste0("final/", names(x), "_resc_rulesmax.tif"), datatype = "INT2U")
  } else {
    # no mask
    out <- calc(x, fun = function(x) round((1 - (x / mx))*10000), filename = paste0("final/", names(x), "_resc_absmax.tif"), datatype = "INT2U")
    out2 <- calc(x, fun = function(x) {
      x[x>mx2] <- mx2
      round((1 - (x / mx2))*10000)
    }, filename = paste0("final/", names(x), "_resc_rulesmax.tif"), datatype = "INT2U")
  }
  return(out)
  }
)
stopCluster(cl)

}))

# create pyramids/stats
message("Building pyramids/statistics...")
writeLines(text = gsub("%SUB%", paste0(getwd(), "/final"), readLines(here("inputs","pyr.py"))), con = here("inputs","pyr2.py"))
system2("C:/Program Files/ArcGIS/Pro/bin/Python/Scripts/propy.bat", here("inputs","pyr2.py"), wait = T)
file.remove(here("inputs","pyr2.py"))


# delete temp dirs
if (delete.temp) {
  unlink("rawImpact", recursive = T,  force = T)
  removeTmpFiles(0)
}

#########################

# mosaic (optional)
system.time({
names(rasSumRescl)[1:2] <- c('x', 'y')
rasSumRescl$fun <- max
rasSumRescl$na.rm <- TRUE
rasSumRescl$filename <- "LandscapeCondition.tif"
rasSum <- do.call(mosaic, rasSumRescl)
})
plot(rasSum)

#########################


# ras <- raster(file)
# plot(ras)


# sub <- c(block$start_row, block$end_row, block$start_col, block$end_col)

###### STARS
# library(stars)
# fun <- function(x) { 
#   if (x > dist) {
#     x <- NA
#   } else {
#     if (dist == 0) x <- wt else x <- (1/(1+(exp(((x/scal)-shift)*spread))))*wt
#   }
#   return(x)
# }   
# 
# i <- 1
# flnm <- rules$filename[i]
# shift <- rules$a[i]
# spread <- rules$b[i]
# dist <- rules$decdist[i]
# scal <- rules$scalar[i]
# wt <- rules$weight[i]
# flnmNoX <- gsub(".tif","", flnm)
# 
# star <- read_stars(rules$filename[i])
# cl <- makeCluster(11, type = "SOCK") 
# clusterExport(cl, c("dist","wt","scal","shift","spread")) 
# system.time(
# tst <- st_apply(star, 1:2, fun, CLUSTER = cl)
# )
# stopCluster(cl)

#######


#write out stack layers so we don't have to calculate again
# for(i in 1:nlayers(stk)){
#     writeRaster(stk[[i]], filename = paste(stk[[i]]@data@names, "_impct2",sep="") , format = "GTiff", overwrite=TRUE)
# }                

#### save/reload the stack
stk <- stackSave(stk, "ImpactRasters.stk")
stk <- stackOpen("ImpactRasters.stk")
#####

# rasSum <- sum(stk[[1]], stk[[2]], stk[[3]], stk[[4]], stk[[5]], 
# 				stk[[6]], stk[[7]], stk[[8]], stk[[9]], stk[[10]], 
# 				stk[[11]], stk[[12]], stk[[13]], 
#                  na.rm = TRUE)

# this is kinda slow
# fun <- function(x) {
#   if (all(is.na(x))) return(0) else return(sum(x, na.rm=T))
# }
# beginCluster()
# system.time(
# rasSum <- clusterR(stk, calc, args = list(fun = fun))
# )

# T


# rescale from 0 (low integrity) to 1 (highest integerity)
mx <- cellStats(rasSum, max)
system.time(
rasSum2 <- clusterR(rasSum, calc, args = list(fun = function(x) 1 - (x / mx)), export = "mx")
)
endCluster()

# save summed raster with datestamp
# writeRaster(rasSum2, filename = paste0("LandscapeCondition_", dt, ".tif"), format = "GTiff", overwrite=TRUE)

# ouptut final raster (masked)
rasSumMsk <- mask(crop(rasSum2, msk), msk, filename = paste0("LandscapeCondition.tif"), format = "GTiff", overwrite=TRUE)

# clean up (delete raster temp dir)
# unlink(tdir)

