# PROJECT: Closed multi-session SCR
# SCRIPT: 02 - Trap locations and state space
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 14 Jan 2026
# LAST MODIFIED: 16 Jan 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# Here I'll read in the trap locations (xj),
# keeping a shapefile and a .csv of x-y coordinates.
# Then, I'll prescribe a rectangular state space
# so that I can use the unif(x, y) prior on ACs si
# I'll save the xlim and ylim for each unit's S
# assume constant across years (reasonable)

# ______________________________________________________________________________
# 2. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)

# ______________________________________________________________________________
# 3. Read in data ----
# ______________________________________________________________________________

# trap locations
trap.kml <- st_read(dsn = paste0(getwd(), "/Spatial data/", lyr = "final_live.kml"))

# unit boundaries (for visualization)
units <- st_read(dsn = paste0(getwd(), "/Spatial data/units_fixed_utm/", lyr = "units_fixed_utm.shp"))

# ______________________________________________________________________________
# 4. Clean trap location data ----
# ______________________________________________________________________________

traps.sf <- trap.kml %>%
  
  # add site name
  mutate(site = rep(c("1A", "1C", "1B",      # 1B and 1C were flipped
                      "2A", "2B", "2C",
                      "3A", "3B", "3C",
                      "4A", "4B", "4C"),
                    each = 36)) %>%
  
  # keep only site and geometry
  dplyr::select(site, geometry) %>%
  
  # arrange by site
  arrange(site) %>%
  
  # transform to UTM
  st_transform(crs = "epsg:32611")

# ______________________________________________________________________________
# 5. Clean unit boundaries ----
# ______________________________________________________________________________

units.sf <- units %>%
  
  # keep only necessary columns
  dplyr::select(name, geometry) %>%
  
  # rename name to site
  rename(site = name) %>%
  
  # arrange by site
  arrange(site)

# ______________________________________________________________________________
# 6. Plot both ----

library(cowplot)

# ______________________________________________________________________________

all.units <- list()

for (i in 1:12) {
  
  focal.unit <- units.sf %>% slice(i)
  focal.traps <- traps.sf %>% filter(site == focal.unit$site)
  
  all.units[[i]] <- ggplot() +
    
    theme_bw() +
    
    geom_sf(data = focal.unit) +
    
    geom_sf(data = focal.traps,
            shape = 21) +
    
    coord_sf(datum = st_crs(32611)) +
    
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  
}

plot_grid(plotlist = all.units)

# ______________________________________________________________________________
# 7. Rectangular state space ----

# buffering the bounding box of the trap grid seems the most reasonable here
# we'll pack these into a list with lapply
traps.split <- split(traps.sf, traps.sf$site)

# ______________________________________________________________________________

make_state_space <- function (
  
  x,  
  buffer = 150
  
) {
  
  # bounding box of traps
  focal.bbox <- spatialEco::bbox_poly(x)
  
  # buffer and extract THAT bbox to make sure it's actually a rectangle
  focal.S <- spatialEco::bbox_poly(st_buffer(focal.bbox, dist = buffer))
  
  return(focal.S)
  
}

# use function
all.S <- lapply(traps.split, make_state_space)

# sanity check
all.S.plots <- list()

for (i in 1:12) {
  
  focal.unit <- units.sf %>% slice(i)
  focal.traps <- traps.sf %>% filter(site == focal.unit$site)
  
  all.S.plots[[i]] <- ggplot() +
    
    theme_bw() +
    
    geom_sf(data = all.S[[i]]) +
    
    geom_sf(data = focal.unit) +
    
    geom_sf(data = focal.traps,
            shape = 21) +
    
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  
}

plot_grid(plotlist = all.S.plots, nrow = 4)

# ______________________________________________________________________________
# 8. Center S and traps on (0, 0) ----
# ______________________________________________________________________________

all.S.center <- list()
traps.center <- list()

for (i in 1:12) {
  
  focal.S <- all.S[[i]]
  focal.traps <- traps.split[[i]]
  
  # center coordinates (what will end up being (0, 0))
  cent.x = (max(st_coordinates(focal.S)[ , 1]) + min(st_coordinates(focal.S)[ , 1])) / 2
  cent.y = (max(st_coordinates(focal.S)[ , 2]) + min(st_coordinates(focal.S)[ , 2])) / 2
  
  # S
  # center S coords
  S.coords.uncenter <- st_coordinates(focal.S)[ , c(1:2)]
    
  S.coords.center <- matrix(c(S.coords.uncenter[ , 1] - cent.x,
                              S.coords.uncenter[ , 2] - cent.y),
                            ncol = 2)
  
  # promote to sf
  all.S.center[[i]] <- st_as_sf(as.data.frame(S.coords.center),
                                coords = c("V1", "V2")) %>%
    
    summarize(geometry = st_combine(geometry)) %>%
    
    st_cast("POLYGON")
  
  # TRAPS
  # center trap coords
  traps.coords.uncenter <- st_coordinates(focal.traps)[ , c(1:2)]
    
  traps.coords.center <- matrix(c(traps.coords.uncenter[ , 1] - cent.x,
                                  traps.coords.uncenter[ , 2] - cent.y),
                                ncol = 2)
    
  # promote to sf
  traps.center[[i]] <- st_as_sf(as.data.frame(traps.coords.center),
                                coords = c("V1", "V2"))
  
}

# ______________________________________________________________________________
# 9. x and y limits for S ----

# each will be a matrix [site, 2]
# we can then append to an array with depth 2
# we'll then save as an .rds

# ______________________________________________________________________________

# xlim
S.xlim <- matrix(
  
  data = c(
    
    st_bbox(all.S.center[[1]])[c(1, 3)],
    st_bbox(all.S.center[[2]])[c(1, 3)],
    st_bbox(all.S.center[[3]])[c(1, 3)],
    st_bbox(all.S.center[[4]])[c(1, 3)],
    st_bbox(all.S.center[[5]])[c(1, 3)],
    st_bbox(all.S.center[[6]])[c(1, 3)],
    st_bbox(all.S.center[[7]])[c(1, 3)],
    st_bbox(all.S.center[[8]])[c(1, 3)],
    st_bbox(all.S.center[[9]])[c(1, 3)],
    st_bbox(all.S.center[[10]])[c(1, 3)],
    st_bbox(all.S.center[[11]])[c(1, 3)],
    st_bbox(all.S.center[[12]])[c(1, 3)]
    
  ),
  
  nrow = 12,
  ncol = 2,
  byrow = T
  
)

# ylim
S.ylim <- matrix(
  
  data = c(
    
    st_bbox(all.S.center[[1]])[c(2, 4)],
    st_bbox(all.S.center[[2]])[c(2, 4)],
    st_bbox(all.S.center[[3]])[c(2, 4)],
    st_bbox(all.S.center[[4]])[c(2, 4)],
    st_bbox(all.S.center[[5]])[c(2, 4)],
    st_bbox(all.S.center[[6]])[c(2, 4)],
    st_bbox(all.S.center[[7]])[c(2, 4)],
    st_bbox(all.S.center[[8]])[c(2, 4)],
    st_bbox(all.S.center[[9]])[c(2, 4)],
    st_bbox(all.S.center[[10]])[c(2, 4)],
    st_bbox(all.S.center[[11]])[c(2, 4)],
    st_bbox(all.S.center[[12]])[c(2, 4)]
    
  ),
  
  nrow = 12,
  ncol = 2,
  byrow = T
  
)

# make array
S.lim <- array(data = c(S.xlim, S.ylim),
               dim = c(12, 2, 2))

# write file
saveRDS(S.lim, file = "Data for model/S_lim.rds")

# ______________________________________________________________________________
# 12. S areas ----

# each will be in ha

# ______________________________________________________________________________

S.area <- unlist(lapply(all.S.center, function (x) {st_area(x) * 0.0001}))

# how much bigger are S than unit areas?
unit.area <- unlist(lapply(split(units.sf, units.sf$site), function (x) {st_area(x) * 0.0001}))

S.area / unit.area

# round this to give data augmentation factors
S.area.df <- data.frame(
  
  area = S.area,
  aug.factor = round(S.area / unit.area)
  
)

# write file
write.csv(S.area.df, "Data for model/S_area.csv", row.names = F)

# ______________________________________________________________________________
# 11. Save trap coordinates ----

# each will be a matrix [36, 2]
# we can then append to an array with depth 12
# we'll then save as an .rds

# ______________________________________________________________________________

trap.coords <- array(data = NA, dim = c(36, 2, 12))

for (i in 1:12) {
  
  trap.coords[ , , i] <- st_coordinates(traps.center[[i]])
  
}

# write file
saveRDS(trap.coords, file = "Data for model/trap_coords.rds")

# ______________________________________________________________________________
# 12. Final sanity check ----
# ______________________________________________________________________________

all.final.plots <- list()

for (i in 1:12) {
  
  all.final.plots[[i]] <- ggplot() +
    
    theme_bw() +
    
    geom_sf(data = all.S.center[[i]]) +
    
    geom_sf(data = traps.center[[i]]) +
    
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
}

plot_grid(plotlist = all.final.plots)

