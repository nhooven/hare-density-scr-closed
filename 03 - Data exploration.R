# PROJECT: Closed multi-session SCR
# SCRIPT: 03 - Data exploration
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 
# LAST MODIFIED: 14 Jan 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# It will be helpful to examine spatial capture histories, as well as total
# spatial captures, as well as summarizing various information on recaptures, etc.

# ______________________________________________________________________________
# 2. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(mefa4)

# ______________________________________________________________________________
# 3. Read in data ----
# ______________________________________________________________________________

# trap locations
trap.kml <- st_read(dsn = paste0(getwd(), "/Spatial data/", lyr = "final_live.kml"))

# capture histories
ch <- read.csv("Data for model/ch.csv")

# covariates
indiv.covs <- read.csv("Data for model/indivs_covs.csv")

# site lookup table
site.lookup <- read.csv("Data for model/site_lookup.csv")

# ______________________________________________________________________________
# 4. Cleaning ----
# ______________________________________________________________________________
# 4a. Trap location data ----
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
# 4b. Add year name to lookup ----
# ______________________________________________________________________________

site.lookup <- site.lookup %>%
  
  mutate(yearName = case_when(
    
    year == 1 ~ "PRE 1",
    year == 2 ~ "PRE 2",
    year == 3 ~ "POST 1",
    year == 4 ~ "POST 2"
    
  )
  
  )

# ______________________________________________________________________________
# 5. Total captures by trap ----

# We'll allow this to vary by occasion, or look at the total

# ______________________________________________________________________________
# 5a. Function ----
# ______________________________________________________________________________

total_caps <- function (
    
  session,
  occ            # total, all, or specific occasions
  
) {
  
  # FILTER BY SESSION
    # indivs
    focal.indivs <- indiv.covs %>% filter(sessionID == session)
    
    # CH
    focal.ch <- ch %>% slice(min(focal.indivs$indivID.1):max(focal.indivs$indivID.1))
    
    # traps
    focal.traps <- traps.sf %>% filter(site == site.lookup$SiteName[site.lookup$sessionID == session])
  
  # compute totals by occasion (as a list)
  total.list <- apply(focal.ch, 2, table)
  
  # loop through occasions
  all.occ <- data.frame()
  
  for (i in 1:length(total.list)) {
    
    focal.occ <- as.data.frame(total.list[[i]])
    
    # coerce to integer
    focal.occ$Var1 <- as.integer(as.character(focal.occ$Var1))
    
    # add traps that did not catch anything
    focal.occ.1 <- rbind(
      
      focal.occ,
      
      data.frame(Var1 = which(1:nrow(focal.traps) %notin% focal.occ$Var1),
                 Freq = 0)
      
    ) %>%
      
      # rename
      rename(trap.id = Var1) %>%
      
      # arrange by trap
      arrange(trap.id) %>%
      
      # remove "no captures"
      filter(trap.id != nrow(focal.traps) + 1) %>%
      
      # add occasion variable
      mutate(occ = i)
    
    # bind in
    all.occ <- rbind(all.occ, focal.occ.1)
  
  }
  
  # pivot - by occasion
  trap.cap.freq.occ <- all.occ %>%
    
    pivot_wider(names_from = occ,
                values_from = Freq)
  
  # summarize in total
  trap.cap.freq.all <- data.frame(trap.id = trap.cap.freq.occ$trap.id,
                                  freq = apply(trap.cap.freq.occ[ , c(2:ncol(trap.cap.freq.occ))], 1, sum))
  
  # plot
  # for now, no need to spatialize any of this
  # maybe if I want to include the unit boundary later or something I can
  
  # bounding box for traps
  traps.bbox <- st_bbox(focal.traps)
  
  # total - all in one map
  if (occ == "total") {
    
    # add in coords
    total.forPlot <- trap.cap.freq.all %>%
      
      mutate(x = st_coordinates(focal.traps)[ , 1],
             y = st_coordinates(focal.traps)[ , 2])
    
    # plot
    out.plot <- ggplot() +
      
      theme_bw() +
      
      geom_point(data = total.forPlot,
                 aes(x = x,
                     y = y,
                     fill = freq),
                 size = 4,
                 shape = 21,
                 stroke = 1) +
      
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            legend.title = element_blank()) +
      
      scale_fill_gradient(low = "white", high = "darkgreen") +
      
      # expand the plotting window a bit
      coord_cartesian(xlim = c(traps.bbox[1] - 40,
                               traps.bbox[3] + 40),
                      ylim = c(traps.bbox[2] - 40,
                               traps.bbox[4] + 40)) +
      
      # title with session info
      ggtitle(paste0(
        
        # unit
        site.lookup$SiteName[site.lookup$sessionID == session],
        
        " - ",
        
        # year
        site.lookup$yearName[site.lookup$sessionID == session]
        
      )
      
      )
    
  }
  
  # add in coords
  all.forPlot <- trap.cap.freq.occ %>%
    
    mutate(x = st_coordinates(focal.traps)[ , 1],
           y = st_coordinates(focal.traps)[ , 2]) %>%
    
    # pivot for plotting
    pivot_longer(cols = 2:(ncol(focal.ch) + 1)) %>%
    
    # factors
    mutate(name = factor(name,
                         labels = paste0("k == ", unique(name))))
  
  # "all" - facetted by occasion
  if (occ == "all") {
    
    # facetted plot
    # plot
    out.plot <- ggplot(all.forPlot) +
      
      theme_bw() +
      
      facet_wrap(~ name) +
      
      geom_point(aes(x = x,
                     y = y,
                     fill = as.factor(value)),
                 size = 2.5,
                 shape = 21,
                 stroke = 0.25) +
      
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            strip.background = element_rect(fill = "white"),
            strip.text = element_text(hjust = 0)) +
      
      scale_fill_manual(values = c("white", "darkgreen")) +
      
      # expand the plotting window a bit
      coord_cartesian(xlim = c(traps.bbox[1] - 40,
                               traps.bbox[3] + 40),
                      ylim = c(traps.bbox[2] - 40,
                               traps.bbox[4] + 40)) +
      
      # title with session info
      ggtitle(paste0(
        
        # unit
        site.lookup$SiteName[site.lookup$sessionID == session],
        
        " - ",
        
        # year
        site.lookup$yearName[site.lookup$sessionID == session]
        
      )
      
      )
    
  }
  
  if (occ %in% 1:ncol(focal.ch)) {
    
    # plot
    out.plot <- ggplot(all.forPlot %>% filter(name == paste0("k == ", occ))) +
      
      theme_bw() +
      
      geom_point(aes(x = x,
                     y = y,
                     fill = as.factor(value)),
                 size = 4,
                 shape = 21,
                 stroke = 1) +
      
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none") +
      
      scale_fill_manual(values = c("white", "darkgreen")) +
      
      # expand the plotting window a bit
      coord_cartesian(xlim = c(traps.bbox[1] - 40,
                               traps.bbox[3] + 40),
                      ylim = c(traps.bbox[2] - 40,
                               traps.bbox[4] + 40)) +
      
      ggtitle(paste0(
        
        # unit
        site.lookup$SiteName[site.lookup$sessionID == session],
        
        " - ",
        
        # year
        site.lookup$yearName[site.lookup$sessionID == session],
        
        " - ",
        
        "k == ", occ
        
        )
    
    )
    
  }
  
  return(out.plot)
  
}

# ______________________________________________________________________________
# 5b. Use function ----
# ______________________________________________________________________________

total_caps(41, "total")


# NEXT:
# Individual capture histories, all spatial capture histories
# Might be better to spend my time on running the model
# and visualize/summarize for the open model