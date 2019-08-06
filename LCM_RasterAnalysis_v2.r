## This code is for applying - here in NY - the "Landscape Condition" methods developed 
## by Colorado NHP/Co State. See SHRP 2 C21A April 2012. 
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

#set a feature extent/mask for final product
# ext <- st_read(here("inputs","LCM_domain.shp"))
ext <- st_read("C:/David/scratch/jurisbnd_lam_clipbound.shp")

# location of distance rasters. If geodatabase, make sure to add extension
ras.dir <- "L:/David/projects/LCM_temp/LCM_data_full_5m_20190802.gdb"
proj.name <- gsub(".gdb$", "", basename(ras.dir))
proj.name <- paste0(proj.name, "_v2")
if (grepl(".gdb$", ras.dir)) ras <- arc.open(ras.dir)@children$RasterDataset else ras <- list.files(ras.dir, pattern = ".tif$", full.names = F)

# create output directory/reload if exists
out.dir <- here("outputs", paste0(proj.name))

if (dir.exists(out.dir)) {
  message("Using existing project files...")
  setwd(out.dir)
  ext <- st_read("ext.shp")
  blocks <- st_read("blocks.shp")
  rules <- read.csv("rules.csv")
  obj <- arc.raster(arc.open(rules$full_path[1]))
  rext <- obj$extent
} else {
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
  st_write(blocks, "blocks.shp", delete_layer = T)
  suppressWarnings(st_write(ext, "ext.shp", delete_layer = T))
}

plot(blocks["id"])
dir.create("rawImpact")

# process blocks
print(system.time({
  for (blk in 1:nrow(blocks)) {
    # arc method w/ chunk processing
    if (blk == 1) {
      message("Starting cluster...")
      tdir <- paste0(dirname(ras.dir), "/", proj.name, "_temp")
      dir.create(tdir, showWarnings = F)
      
      # cluster apply, should be faster when many cores available
      # using too many cores + large blocks can bog down system due to high memory usage
      cl <- makeCluster(min(9, nrow(rules)), type = "SOCK")  # max arc clusters appears to be 9.
      
      clusterEvalQ(cl, library(arcgisbinding))
      clusterCall(cl, arc.check_product) # takes a while, can error if too many cores (>9; maybe due to an Arc license issue)
      clusterExport(cl, c("ras.dir","tdir"))
    }
    
    # get block info
    block <- blocks[blk,]
    blkid <- block$id
    if (file.exists(paste0("rawImpact/prodImpact", blkid, ".tif"))) next
    sub <- c(block$start_row, block$end_row, block$start_col, block$end_col)
    blkrast <- paste0(tempfile(), ".tif")
    tempr <- arc.raster(NULL, path=blkrast, extent=st_bbox(block), 
                       nrow=block$end_row-block$start_row, ncol=block$end_col-block$start_col, 
                       nband=1, pixel_type="F32", sr=obj$sr, overwrite=TRUE)
    clusterExport(cl, c("tempr","sub", "blkid"))
    
    # process raster values
    message("Creating impact rasters for block ", blkid, "...")
    out <- parApply(cl, rules, 1, FUN = function(x) {
      flnm <- x["filename"]
      r <- arc.raster(arc.open(x["full_path"]))
      file <- paste0(tdir, "/", flnm , "_calc", blkid, ".tif")
      r2 = arc.raster(NULL, path = file, dim=dim(tempr), pixel_type="F32", 
                      nodata=tempr$nodata, extent=tempr$extent, sr=tempr$sr, overwrite = T)
      v <- r$pixel_block(ul_y = sub[1], nrow = sub[2]-sub[1], ul_x = sub[3], ncol = sub[4] - sub[3])
      
      # v2 decdist function (see LCM_impact_fn.R)
      si <- as.numeric(x["siteindex"])
      ds <- as.numeric(x["distthresh"])
      # if (!all(is.na(v))) {
        v[!is.na(v)] <- (1/(1+exp(-2*((v[!is.na(v)]-(ds*0.5))/(ds*0.25))))) ^ (-0.25*log(si)) #-8*10-16)
        r2$write_pixel_block(v, nrow = sub[2]-sub[1])
      # }
      
      r2$commit(opt = c("build-pyramid"))
      return(file)
    })
    
    arc.delete(blkrast)
    
    stk <- stack(out)
    stk <- stackSave(stk, paste0("ImpactRasters", blkid, ".stk"))
    # stk <- stackOpen(paste0("ImpactRasters", blkid, ".stk"))
    
    # v2: product of impact rasters
    message("Combining impact rasters...")
    fun <- function(x) {
      if (all(is.na(x))) return(1) else return(prod(x, na.rm=T))
    }
    # uses cluster from above
    ras2 <- clusterR(stk, calc, args = list(fun = fun), cl = cl, 
                     filename = paste0("rawImpact/prodImpact", blkid, ".tif"), overwrite = TRUE)
    
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
delete.temp <- T # delete intermediate files after completion?
proj.name <- "LCM_data_full_30m_20190802_v2"
setwd(here("outputs", proj.name))
blocks <- st_read("blocks.shp")
ext <- st_read("ext.shp")

#  process product rasters
lr <- list.files(path = "rawImpact", pattern = "Impact.*tif$", full.names = T)
rlist <- lapply(lr, raster)
cl <- makeCluster(min(length(rlist), parallel::detectCores()-1), type = "SOCK") 

print(system.time({
  clusterEvalQ(cl, library(raster))
  clusterExport(cl, c("rlist"))
  clusterEvalQ(cl, library(fasterize))
  clusterEvalQ(cl, library(sf))
  clusterExport(cl, c("blocks","ext","delete.temp"))
  
  # rescale to integer (multiply values)
  # also masks (only if necessary), using 'ext' polygon
  dir.create("final")
  rasSumRescl <- parLapply(cl, rlist, fun = function(x) {
    nm <- names(x)
    blkid <- as.integer(gsub("\\D", "", nm))
    if (st_area(st_intersection(ext, blocks[blocks$id == blkid,])) != st_area(blocks[blocks$id == blkid,])) {
      # mask
      m1 <- fasterize(st_intersection(ext, blocks[blocks$id == blkid,]), x) 
      out <- mask(calc(x, fun = function(x) round(x*10000)), m1, filename = paste0("final/", names(x), "_resc.tif"), datatype = "INT2U")
    } else {
      # no mask
      out <- calc(x, fun = function(x) round(x*10000), filename = paste0("final/", names(x), "_resc.tif"), datatype = "INT2U")
    }
    if (delete.temp) unlink(slot(x@file, "name"), force = T)
    return(out)
  }
  )
  stopCluster(cl)
}))


#########################

# mosaic (optional; would take a long time)
# system.time({
#   names(rasSumRescl)[1:2] <- c('x', 'y')
#   rasSumRescl$fun <- max
#   rasSumRescl$na.rm <- TRUE
#   rasSumRescl$filename <- "LandscapeCondition.tif"
#   rasSum <- do.call(mosaic, rasSumRescl)
# })
# plot(rasSum)

#########################

# create pyramids/stats
message("Building pyramids/statistics...")
writeLines(text = gsub("%SUB%", getwd(), readLines(here("inputs","pyr.py"))), con = here("inputs","pyr2.py"))
system2("C:/Program Files/ArcGIS/Pro/bin/Python/Scripts/propy.bat", here("inputs","pyr2.py"), wait = T)
file.remove(here("inputs","pyr2.py"))
if (delete.temp) unlink("rawImpact", recursive = T)