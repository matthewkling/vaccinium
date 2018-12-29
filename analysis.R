
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

# county climate means and variances -- SLOW
z <- rasterize(counties, raster(f[1]), "OBJECTID")
for(v in c("tmin", "tmax", "tmean")){
      r <- stack(f[grepl(v, f)])
      avg <- zonal(r, z, "mean")
      stdev <- zonal(r, z, "sd")
      saveRDS(avg, paste0("e:/vaccinium/data/derived/", v, "_county_means.rds"))
      saveRDS(stdev, paste0("e:/vaccinium/data/derived/", v, "_county_stdevs.rds"))
}



### calculate anomalies for each occurrence ###


window_anomaly <- function(x, # a vector of the focal year, focal julian day, and then followed by climate time series 
                           lag=0, # vector of month offsets. 0 for the focal month, 1:3 for the three mos preceeding, etc
                           fun=mean,
                           start_year=1950, # these must match the climate data date range
                           end_year=2010
){
      # sort components
      if(any(is.na(x[1:3]))) return(NA)
      y <- x[1]
      j <- x[2]
      x <- x[3:length(x)]
      if(y < start_year+1 | y > end_year-1) return(NA) # buffer needed so offsets can spill into adjoining years
      
      # monthly means
      x <- matrix(x, ncol=length(x)/12)
      m <- apply(x, 1, mean)
      
      # anomalies
      x <- x[,y-start_year+1 + c(-1,0,1)] # isolate data from focal year and neighbors
      x <- x - cbind(m, m, m) # compute anomaly
      x <- as.vector(x) # 36-month window, needed so offsets can spill into adjoining years
      
      # isolate values for focal months
      j <- floor(j/365*12)
      x <- x[j-lag+12]
      fun(x)
}



window_mean <- function(x, # climate time series 
                        lag=0 # vector of month offsets. 0 for the focal month, 1:3 for the three mos preceeding, etc
){
      if(is.na(x[1])) return(NA)
      x <- matrix(x, ncol=length(x)/12)
      lag <- lag + 1
      lag[lag <= 0] <- lag[lag <= 0] + 12
      lag[lag > 12] <- lag[lag > 12] - 12
      x <- x[lag,]
      mean(x)
}



# define seasonal temperature variables
vars <- data.frame(name=c("annual_mean", "winter_min", "spring_mean", "summer_max"),
                   var=c("tmean", "tmin", "tmean", "tmax"))
vars$months <- list(0:11, -1:1, 0:2, c(-7:-5, 5:7)) # lags, relative to january


# county ID for each species occurrence
s <- read.csv("e:/vaccinium/data/raw/Blueberries_Meineke_raw.csv", stringsAsFactors=F)
coordinates(s) <- c("decimalLongitude", "decimalLatitude")
s$OBJECTID <- extract(z, s)
s <- as.data.frame(s)


for(i in 1:nrow(vars)){
      
      # load county climate data and array rows to match species occurrence order
      clim <- readRDS(paste0("e:/vaccinium/data/derived/", vars$var[i], "_county_means.rds"))
      clim <- clim[match(s$OBJECTID, clim[,"zone"]), 2:ncol(clim)]
      
      # temporal anomaly
      s[,paste0(vars$name[i], "_temporal")] <- apply(cbind(s$year, 0, clim), 1, 
                                                          FUN=window_anomaly, start_year=1895, end_year=2017, lag=vars$months[[i]])
      
      # long-term mean climate at each location
      means <- apply(clim, 1, FUN=window_mean, lag=vars$months[[i]])
      
      # spatial anomaly (species-specific)
      for(sp in unique(s$plant_genus_species)){
            spi <- s$plant_genus_species==sp
            s[spi, paste0(vars$name[i], "_spatial")] <- means[spi] - mean(means[spi], na.rm=T)
      }
      
      # spatial variance within county during the average year -- a measure of uncertainty
      clim <- readRDS(paste0("e:/vaccinium/data/derived/", vars$var[i], "_county_stdevs.rds"))
      clim <- clim[match(s$OBJECTID, clim[,"zone"]), 2:ncol(clim)] ^ 2
      s[,paste0(vars$name[i], "_spatial_variance")] <- apply(clim, 1, FUN=window_mean, lag=vars$months[[i]])
      
}

# export results
s <- select(s, -OBJECTID)
write.csv(s, "e:/vaccinium/output/Blueberries_Meineke_anomalies.csv", row.names=F)

