
library(raster)
library(tidyverse)


### summarize climate variables by US county ###

# climate rasters
f <- list.files("F:/PRISM_AN81m", pattern="\\.bil", full.names=T)
f <- f[!grepl("aux.xml", f)]
f <- f[!grepl("_2018", f)] # drop partial year

# counties
counties <- getData("GADM", country="USA", level=2) %>%
      spTransform(crs(raster(f[1])))

# county climate means and variances
z <- rasterize(counties, raster(f[1]), "OBJECTID")
for(v in c("tmin", "tmax", "tmean")){
      r <- stack(f[grepl(v, f)])
      avg <- zonal(r, z, "mean")
      stdev <- zonal(r, z, "sd")
      saveRDS(avg, paste0("e:/vaccinium/data/derived/", v, "_county_means.rds"))
      saveRDS(stdev, paste0("e:/vaccinium/data/derived/", v, "_county_stdevs.rds"))
}





