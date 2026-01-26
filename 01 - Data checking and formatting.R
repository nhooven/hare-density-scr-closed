# PROJECT: Closed multi-session SCR
# SCRIPT: 01 - Data checking and formatting
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 12 Jan 2026
# COMPLETED: 12 Jan 2026
# LAST MODIFIED: 16 Jan 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(mefa4)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________
# 2a. Open capture histories ----

# this includes all individuals and their mark-recapture IDs (MRID)

# ______________________________________________________________________________

open.ch <- read.csv("Trapping data/hare_open_CapHis.csv")

# keep needed columns and rename to make everything easier
open.ch <- open.ch %>%
  
  dplyr::select(AnimalID.1,
                MRID,
                Sex,
                MR22f,
                MR23,
                MR24,
                MR25) %>%
  
  rename(y1 = MR22f,
         y2 = MR23,
         y3 = MR24,
         y4 = MR25)

# ______________________________________________________________________________
# 2b. Mark-recapture data ----

# these are J x K frames of MRIDs, that I will need to manipulate for SCR
# we'll extract both individual capture histories and a trap operation matrix from this
# we can also pull session-specific Ks from here

# ______________________________________________________________________________

mr.data.2022 <- read.csv("Trapping data/trapping_2022.csv")
mr.data.2023 <- read.csv("Trapping data/trapping_2023.csv")
mr.data.2024 <- read.csv("Trapping data/trapping_2024.csv")
mr.data.2025 <- read.csv("Trapping data/trapping_2025.csv")

# ______________________________________________________________________________
# 3. Site / session / treatment lookup table ----
# ______________________________________________________________________________

site.lookup <- data.frame(
  
  # integer session ID
  sessionID = 1:41,
  
  # character site name
  SiteName = c(c("2A", "2B", "2C", "3A", "3C"),
               rep(c("1A", "1B", "1C", "2A", "2B", "2C",
                     "3A", "3B", "3C", "4A", "4B", "4C"),
                   times = 3)),
  
  # integer site ID
  siteID = c(c(4, 5, 6, 7, 9),
             rep(1:12,
                 times = 3)),
  
  # integer year T
  year = c(rep(1, times = 5),
           rep(2, times = 12),
           rep(3, times = 12),
           rep(4, times = 12)),
  
  # treatment variable
  trt = c(rep("untreated", times = 17),
          rep(c("ret", "pil", "untreated",
                "pil", "ret", "untreated",
                "pil", "ret", "untreated",
                "ret", "pil", "untreated"),
              times = 2))

)

# write to .csv for reference
write.csv(site.lookup, "Data for model/site_lookup.csv", row.names = F)

# ______________________________________________________________________________
# 4. Bind MR data together with unique session IDs ----
# ______________________________________________________________________________

# 2022
mr.data.2022.1 <- mr.data.2022 %>%
  
  # subset only sessions we want
  filter(Site %in% c("2A", "2B", "2C", "3A", "3C")) %>%
  
  # remove the last two columns
  dplyr::select(-c("D9", "D10")) %>%
  
  # replace missing data explictly with NAs
  mutate_at(vars(matches("D")), ~ifelse(. == "", NA, .)) %>%
  
  # add year variable
  mutate(year = 1) %>%
  
  # rename Site
  rename(SiteName = Site) %>%
  
  # join in other index variables from site.lookip
  left_join(site.lookup, by = c("SiteName", "year")) %>%
  
  # rearrange sensibly
  dplyr::select(sessionID, siteID, year, SiteName, trt, Trap, D1:D8)

# 2023
mr.data.2023.1 <- mr.data.2023 %>%
  
  # replace missing data explictly with NAs
  mutate_at(vars(matches("D")), ~ifelse(. == "", NA, .)) %>%
  
  # add year variable
  mutate(year = 2) %>%
  
  # rename Site
  rename(SiteName = Site) %>%
  
  # join in other index variables from site.lookip
  left_join(site.lookup, by = c("SiteName", "year")) %>%
  
  # rearrange sensibly
  dplyr::select(sessionID, siteID, year, SiteName, trt, Trap, D1:D6) %>%
  
  # arrange by sessionID
  arrange(sessionID)

# 2024
mr.data.2024.1 <- mr.data.2024 %>%
  
  # replace missing data explictly with NAs
  mutate_at(vars(matches("D")), ~ifelse(. == "", NA, .)) %>%
  
  # add year variable
  mutate(year = 3) %>%
  
  # rename Site
  rename(SiteName = Site) %>%
  
  # join in other index variables from site.lookip
  left_join(site.lookup, by = c("SiteName", "year")) %>%
  
  # rearrange sensibly
  dplyr::select(sessionID, siteID, year, SiteName, trt, Trap, D1:D8) %>%
  
  # arrange by sessionID
  arrange(sessionID)

# 2025
mr.data.2025.1 <- mr.data.2025 %>%
  
  # replace missing data explictly with NAs
  mutate_at(vars(matches("D")), ~ifelse(. == "", NA, .)) %>%
  
  # add year variable
  mutate(year = 4) %>%
  
  # rename Site
  rename(SiteName = Site) %>%
  
  # join in other index variables from site.lookip
  left_join(site.lookup, by = c("SiteName", "year")) %>%
  
  # rearrange sensibly
  dplyr::select(sessionID, siteID, year, SiteName, trt, Trap, D1:D6) %>%
  
  # arrange by sessionID
  arrange(sessionID)

# bind all together (should have 8 total occasion columns)
mr.data <- bind_rows(mr.data.2022.1, mr.data.2023.1, mr.data.2024.1, mr.data.2025.1)

# ______________________________________________________________________________
# 5. Check that MRIDs match exactly ----
# ______________________________________________________________________________
# 5a. Check function ----
# ______________________________________________________________________________

check_mrid <- function (
    
  focalID
  
) {
  
  # DATA SETUP
  
  # extract correct indices
  focal.indices <- site.lookup %>% filter(sessionID == focalID)
  
  # keep correct column (based on primary occasion)
  if (focal.indices$year == 1) {open.ch.1 <- open.ch[ , c(1:3, 4)]}
  if (focal.indices$year == 2) {open.ch.1 <- open.ch[ , c(1:3, 5)]}
  if (focal.indices$year == 3) {open.ch.1 <- open.ch[ , c(1:3, 6)]}
  if (focal.indices$year == 4) {open.ch.1 <- open.ch[ , c(1:3, 7)]}
  
  # change name to "y" for generality
  names(open.ch.1)[4] <- "y"
  
  # focal indivs
  focal.indivs <- open.ch.1[which(open.ch.1$y == focal.indices$SiteName), ]
  
  # subset just the secondary session CHs
  mr <- mr.data %>%
    
    filter(sessionID == focalID) %>%
    
    dplyr::select(D1:D8)
  
  # CHECK IDs
  
  # generate list of MRIDs from entire trapping data
  # bind together into one vector
  all.outcomes <- as.vector(as.matrix(mr))
  
  # remove trap indicators
  all.outcomes.unique <- unique(all.outcomes)
  
  trap.mrid <- all.outcomes.unique[! all.outcomes.unique %in% c("O", "B", "X", "C", "E", NA)]
  
  # compare to focal.indivs
  # which MRIDs are IN MR and NOTIN indivs?
  iMR.niIndivs <- trap.mrid[which(trap.mrid %notin% focal.indivs$MRID)]
  
  # which MRIDs are NOTIN MR and IN indivs?
  niMR.iIndivs <- focal.indivs$MRID[which(focal.indivs$MRID %notin% trap.mrid)]
  
  # dfs to hold information
  # in MR and not in Indivs
  if (length(iMR.niIndivs) > 0) {
    
    iMR.niIndivs.df <- data.frame(
      
      MRID = iMR.niIndivs,
      which.list = "iMR.niIndivs",
      where.else = ifelse(length(open.ch.1$y[which(open.ch.1$MRID %in% iMR.niIndivs)]) > 0,
                          open.ch.1$y[which(open.ch.1$MRID %in% iMR.niIndivs)],
                          NA)
      
    )
    
  } else {
    
    iMR.niIndivs.df <- 
      
      data.frame(
        
        MRID = NA,
        which.list = "iMR.niIndivs",
        where.else = NA
        
      )
    
  }
  
  # not in MR and in Indivs
  if (length(niMR.iIndivs) > 0) {
    
    niMR.iIndivs.df <- data.frame(
      
      MRID = niMR.iIndivs,
      which.list = "niMR.iIndivs",
      where.else = NA
      
    )
    
  } else {
    
    niMR.iIndivs.df <- data.frame(
      
      MRID = NA,
      which.list = "niMR.iIndivs",
      where.else = NA
      
    )
    
  }
  
  check.df <- rbind(iMR.niIndivs.df, niMR.iIndivs.df)
  
  # return
  return(check.df)
  
}

# ______________________________________________________________________________
# 5b. Change conditions ----

# if we have individuals captured originally in (an)other unit(s),
# change their MR entry(ies) to "C"
# this is to account for the trap being partially closed during the night.
# We'll keep those individuals with their original unit and assume any
# their AC is somewhere between the units.

# an alternative is to merge the grids (for 2A/2B, for example)
# but this complicates any inference on treatment effects
# so for simplicity, since grid-switching (within year) is rare,
# we'll ignore it for the most part

# I think it'll be better to do this manually
# we'll keep code for any we change

# ______________________________________________________________________________

# session 1
check_mrid(1)

mr.data[4, 12] <- "C"
mr.data[20, 10] <- "C"

# session 8
check_mrid(8)

mr.data[265, 10] <- "C"

# session 10
mr.data[339, 10] <- "C"

# session 33
mr.data[1159, 10] <- "C"
mr.data[1166, 11] <- "C"

# session 38
mr.data[1333, 10] <- "C"

# ______________________________________________________________________________
# 5c. Sanity check ----
# ______________________________________________________________________________

all.checks <- data.frame()

for (i in 1:41) {
  
  focal.check <- check_mrid(i)
  
  all.checks <- rbind(all.checks, focal.check)
  
}

# ______________________________________________________________________________
# 6. Reformat MR data ----

# each function will output an n-length data.frame

# function to extract individuals as a vector
indiv_vect <- function (x.trap) {
  
  # convert to vector
  x.trap.vector <- as.vector(as.matrix(x.trap))
  
  # extract unique individuals as a vector
  x.mrid <- unique(x.trap.vector)[! unique(x.trap.vector) %in% c("O", "B", "X", "C", "E", NA)]
  
  return(x.mrid)
  
}

# split for lapply
mr.data.split <- split(mr.data, mr.data$sessionID)

# ______________________________________________________________________________
# 6a. Categorical capture histories ----
# ______________________________________________________________________________

# function
catCH <- function (x) {
  
  # subset trapping histories only
  x.trap <- x %>% dplyr::select(D1:D8)
  
  # extract unique individuals as a vector
  x.mrid <- indiv_vect(x.trap)
  
  # how many occasions?
  K = sum(apply(apply(x.trap, 2, is.na), 2, sum) == 0)
  
  # set up matrix
  x.trap.mat <- matrix(data = NA,
                       nrow = length(x.mrid),
                       ncol = K)
  
  # loop through MRIDs (i)
  for (i in 1:length(x.mrid)) {
    
    # loop through occasions (K)
    for (k in 1:K) {
      
      # find which trap (if any) the focal individual was found in [i, k]
      if (length(which(x.trap[, k] == x.mrid[i])) > 0) {
        
        x.trap.mat[i, k] <- which(x.trap[, k] == x.mrid[i])
        
        # if not, add the "not trapped" index (i.e., 37)
      } else {
        
        x.trap.mat[i, k] <- nrow(x.trap) + 1
        
      }
      
    }
    
  }
  
  # coerce to data.frame
  x.trap.df <- as.data.frame(x.trap.mat)
  
  return(x.trap.df)
  
}

# use function
all.catCH <- bind_rows(lapply(mr.data.split, catCH))

# ______________________________________________________________________________
# 6b. Binary trap response covariate ----
# ______________________________________________________________________________

# function
prev_cap <- function (x) {
  
  # subset trapping histories only
  x.trap <- x %>% dplyr::select(D1:D8)
  
  # extract unique individuals as a vector
  x.mrid <- indiv_vect(x.trap)
  
  # how many occasions?
  K = sum(apply(apply(x.trap, 2, is.na), 2, sum) == 0)
  
  # extract mr.mat from previous function
  mr.mat <- catCH(x)
  
  # set up matrix
  prev.cap.mat <- matrix(data = NA,
                         nrow = length(x.mrid),
                         ncol = K)
  
  # first occasion must be zero
  prev.cap.mat[ , 1] <- 0
  
  # loop through MRIDs (i)
  for (i in 1:length(x.mrid)) {
    
    # loop through occasions (K)
    for (k in 2:K) {
      
      # if was captured on ANY previous occasions, add a 1
      if (any(mr.mat[i, 1:(k - 1)] %in% 1:nrow(x.trap))) {
        
        prev.cap.mat[i, k] <- 1
        
        # if not, add a 0
      } else {
        
        prev.cap.mat[i, k] <- 0
        
      }
      
    }
    
  }
  
  # coerce to data.frame
  prev.cap.df <- as.data.frame(prev.cap.mat)
  
  return(prev.cap.df)
  
}

# use function
all.prev.cap <- bind_rows(lapply(mr.data.split, prev_cap))

# ______________________________________________________________________________
# 6c. Individual indices and covariates ----

# these will be sessionID, siteID, clusterID, treatment, year (for the session)
# and sex (for the individual)

# ______________________________________________________________________________

indiv_covs <- function (x) {
  
  # subset trapping histories only
  x.trap <- x %>% dplyr::select(D1:D8)
  
  # correct cluster
  clusterID = case_when(x$siteID[1] %in% c(1:3) ~ 1,
                        x$siteID[1] %in% c(4:6) ~ 2,
                        x$siteID[1] %in% c(7:9) ~ 3,
                        x$siteID[1] %in% c(10:12) ~ 4)
  
  # extract unique individuals as a vector
  x.mrid <- indiv_vect(x.trap)
  
  # extract correct sex
  # extract correct indices
  focal.indices <- site.lookup %>% filter(sessionID == x$sessionID[1])
  
  # keep correct column (based on primary occasion)
  if (focal.indices$year == 1) {open.ch.1 <- open.ch[ , c(1:3, 4)]}
  if (focal.indices$year == 2) {open.ch.1 <- open.ch[ , c(1:3, 5)]}
  if (focal.indices$year == 3) {open.ch.1 <- open.ch[ , c(1:3, 6)]}
  if (focal.indices$year == 4) {open.ch.1 <- open.ch[ , c(1:3, 7)]}
  
  # change name to "y" for generality
  names(open.ch.1)[4] <- "y"
  
  # focal indivs
  focal.indivs <- open.ch.1[which(open.ch.1$y == focal.indices$SiteName), ] %>%
    
    dplyr::select(MRID, Sex)
  
  # new df, join in (must be in the same order as CHs, prev. cap, etc)
  indiv.covs <- data.frame(
    
    "MRID" = x.mrid
    
  ) %>%
    
    left_join(focal.indivs) %>%
    
    # add site indices
    mutate(
      
      sessionID = x$sessionID[1],
      siteID = x$siteID[1],
      clusterID = clusterID,
      ret = ifelse(x$trt[1] == "ret", 1, 0),
      pil = ifelse(x$trt[1] == "pil", 1, 0),
      year = x$year[1]
      
    ) %>%
    
    mutate(
      
      # add a unique indivID
      indivID = paste0(MRID, "_", sessionID),
      
      # change Sex to binary (with NAs)
      Sex = case_when(
        
        Sex == "F" ~ 0,
        Sex == "M" ~ 1,
        Sex == "U" ~ NA
        
      )
      
    )
  
  return(indiv.covs)
  
}

# use function
all.indivs.covs <- bind_rows(lapply(mr.data.split, indiv_covs))

# add integer indivID
all.indivs.covs$indivID.1 <- 1:nrow(all.indivs.covs)

# ______________________________________________________________________________
# 6d. Write to file ----
# ______________________________________________________________________________

write.csv(all.catCH, "Data for model/ch.csv", row.names = F)
write.csv(all.prev.cap, "Data for model/prev_cap.csv", row.names = F)
write.csv(all.indivs.covs, "Data for model/indivs_covs.csv", row.names = F)

# ______________________________________________________________________________
# 7. Trap operation matrix ----

# J x K matrix
# The values in here will be multiplied within the model to allow for or "zero out"
# trap-specific capture probabilities
# Cheekily, I will include 0.5 as a possible value to encode our uncertainty about
# trap availability given the trap was:
# (1) found closed or 
# (2) caught a bycatch species

# we won't know how long the trap was available for a hare capture, so we consider
# its availability to be a Bernoulli trial. Maybe it was open all night and closed
# or caught another critter in the morning after no hare was going to go in. Maybe
# it closed or caught another critter shortly after we set it.

# In effect, the multiplication by 0.5 allows the full capture probability for
# half of all MCMC draws, and zeroes it out for the other half

# any "escaped" hares that did not receive an AnimalID will also induce a 0.5 here

# ______________________________________________________________________________

# function
make_trap_op <- function (x) {
  
  # subset trapping histories only
  x.trap <- x %>% dplyr::select(D1:D8)
  
  # how many occasions?
  K = sum(apply(apply(x.trap, 2, is.na), 2, sum) == 0)
  
  # unlist and assign
  focal.mr.vect <- unlist(x.trap)
  
  focal.mr.vect[which(focal.mr.vect %notin% c("X", "C", "B", "E"))] <- 1
  focal.mr.vect[which(focal.mr.vect %in% c("X"))] <- 0
  focal.mr.vect[which(focal.mr.vect %in% c("C", "B", "E"))] <- 0.5
  
  trap.op <- as.data.frame(matrix(as.numeric(focal.mr.vect),
                                  nrow = nrow(x.trap),
                                  ncol = K))
  
  return(trap.op)
  
}

# use function
all.trap.op <- bind_rows(lapply(mr.data.split, make_trap_op))

# save to file
write.csv(all.trap.op, "Data for model/trap_op.csv", row.names = F)

# ______________________________________________________________________________
# 8. Occasions by session ----
# ______________________________________________________________________________

occ_sess <- function (x) {
  
  # subset trapping histories only
  x.trap <- x %>% dplyr::select(D1:D8)
  
  # how many occasions?
  K = sum(apply(apply(x.trap, 2, is.na), 2, sum) == 0)
  
  return(K)
  
}

all.occ.sess <- data.frame(sessionID = 1:41,
                           K = unlist(lapply(mr.data.split, occ_sess)))

# save to file
write.csv(all.occ.sess, "Data for model/occ_sess.csv", row.names = F)
