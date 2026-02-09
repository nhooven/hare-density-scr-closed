# PROJECT: Closed multi-session SCR
# SCRIPT: 04 - Prepare all data for modeling
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 20 Jan 2026
# LAST MODIFIED: 05 Feb 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 0. Purpose ----
# ______________________________________________________________________________

# I would like all of this to be done in a separate script so all I need to 
# feed to Kamiak is the modeling script and the constant and data lists

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# capture histories
ch <- read.csv("Data for model/ch.csv")

# previous capture
prev.cap <- read.csv("Data for model/prev_cap.csv")

# covariates
indiv.covs <- read.csv("Data for model/indivs_covs.csv")

# site lookup table
site.lookup <- read.csv("Data for model/site_lookup.csv")

# occasions by session
occ.sess <- read.csv("Data for model/occ_sess.csv")

# trap operation matrix
trap.op <- read.csv("Data for model/trap_op.csv")

# trap coords
trap.coords <- readRDS("Data for model/trap_coords.rds")

# S limits
S.lim <- readRDS("Data for model/S_lim.rds")

# S areas
S.area.df <- read.csv("Data for model/S_area.csv")

# ______________________________________________________________________________
# 3. Data cleaning ----
# ______________________________________________________________________________
# 3a. Trap operation ----

# I want this to be a [J, max(K), G] array

# ______________________________________________________________________________
  
trap.op.1 <- array(data = NA, dim = c(36, 8, 41))

for (i in 1:41) {
  
  trap.op.1[ , , i] <- as.matrix(trap.op[(1:36) + (36 * (i - 1)), ])
  
}

# ______________________________________________________________________________
# 4. Define required quantities ----
# ______________________________________________________________________________

# number of sites U
U = 12

# number of sessions G
G = 41

# number of captured individuals
n.g <- table(indiv.covs$sessionID) # by session
n = sum(n.g)                      # total

# number of traps J
J = 36

# number of occasions by session
K.g <- occ.sess$K

# S limits [U, 2, 2]
S.lim

# trap coords [J, 2, U]
trap.coords

# S areas (in ha) [U]
S.areas <- S.area.df$area

# ______________________________________________________________________________
# 5. Data augmentation ----
# ______________________________________________________________________________

# define how many individuals to add by site
aug.factor.bysession <- data.frame(
  
  sessionID = 1:41
  
) %>%
  
  left_join(site.lookup %>% dplyr::select(sessionID, siteID)) %>%
  
  left_join(
    
    S.area.df %>% dplyr::select(aug.factor) %>% mutate(siteID = 1:12)
    
  )

# totals by session (add 4 to the aug.factor)
n.aug.g <- n.g * (aug.factor.bysession$aug.factor + 4)

# add a few more for sessions: 24, 27, 28, 36, 39, 40, 41
n.aug.g[24] <- 175
n.aug.g[27] <- 100
n.aug.g[28] <- 80
n.aug.g[36] <- 100
n.aug.g[39] <- 160
n.aug.g[40] <- 80
n.aug.g[41] <- 50

n.aug = sum(n.aug.g)

# total individuals in the dataset
M = n + n.aug

# data augmentation variable z
z <- c(rep(1, times = n), rep(NA, times = n.aug))

# zeroes variable (Chandler method)
zeroes <- c(rep(NA, times = n), rep(0, times = n.aug))

# add zero rows to prev.cap
prev.cap.1 <- rbind(
  
  prev.cap,
  
  matrix(0,
         nrow = n.aug,
         ncol = ncol(prev.cap))
  
)

# covariate assignment
# sessionID
aug.sessionID <- vector()

for (i in 1:41) {aug.sessionID <- c(aug.sessionID, rep(i, each = n.aug.g[i]))}

# siteID
aug.siteID <- vector()
  
for (i in 1:length(aug.sessionID)) {aug.siteID <- c(aug.siteID,
                                                    site.lookup$siteID[site.lookup$sessionID == aug.sessionID[i]])}

# data.frame
indiv.covs.aug <- data.frame(
  
  MRID = NA,
  Sex = NA,
  sessionID = aug.sessionID
  
) %>%
  
  # join in siteID and year (and trt)
  left_join(site.lookup %>% dplyr::select(sessionID, siteID, year, trt)) %>%
  
  # add cluster ID
  mutate(clusterID = case_when(siteID %in% c(1:3) ~ 1,
                               siteID %in% c(4:6) ~ 2,
                               siteID %in% c(7:9) ~ 3,
                               siteID %in% c(10:12) ~ 4)) %>%
  
  # add columns
  mutate(
    
    # add ret and pil
    ret = ifelse(trt == "ret", 1, 0),
    pil = ifelse(trt == "pil", 1, 0),
    
    # indivID
    indivID = NA,
    indivID.1 = (n + 1):M) %>%
  
  # keep columns in right order
  dplyr::select(MRID, Sex, sessionID, siteID, clusterID, ret, pil, year, indivID, indivID.1)

# bind all indiv covs together
indiv.covs.M <- rbind(indiv.covs, indiv.covs.aug)

# ______________________________________________________________________________
# 6. Correct covariates for detection parameter linear predictors ----

# these will be: 
  # post1 and post2 for years

# ______________________________________________________________________________

indiv.covs.M <- indiv.covs.M %>%
  
  mutate(
    
    # forest type
    ft = ifelse(clusterID == 4, 1, 0),
    
    # post1
    post1 = ifelse(year == 3, 1, 0),
    
    # post2
    post2 = ifelse(year == 4, 1, 0)
    
  )


# ______________________________________________________________________________
# 6. Build lists for model ----
# ______________________________________________________________________________
# 6a. Constants ----
# ______________________________________________________________________________  

constant.list <- list(
  
  # scalar constants
  n = n,
  M = M,
  J = J,
  G = G,
  U = U,
  
  # site-specific constants [U]
  S.areas = S.areas,
  
  # session specific constants [G]
  K = K.g,
  site.g = site.lookup$siteID[site.lookup$sessionID == 1:41],
  
  # indices for each individual [M]
  session = indiv.covs.M$sessionID,
  site = indiv.covs.M$siteID,
  cluster = indiv.covs.M$clusterID,
  ft = indiv.covs.M$ft
  
)

# ______________________________________________________________________________
# 6b. Data ----
# ______________________________________________________________________________  
  
data.list <- list (
  
  
  # individual data [n]
  ch = ch,
  
  # individual data [M]
  prev.cap = prev.cap.1,
  sex = indiv.covs.M$Sex,
  
  # treatment and year indicators
  ret = indiv.covs.M$ret,
  pil = indiv.covs.M$pil,
  post1 = indiv.covs.M$post1,
  post2 = indiv.covs.M$post2,
  
  # trap operation [J, max(K), G]
  trap.op = trap.op.1,
  
  # state space limits [U, 2, 2]
  S.lim = S.lim,
  
  # trap coordinates [J, 2, U]
  trap.coords = trap.coords,
  
  # data augmentation [M]
  z = z,
  zeroes = zeroes
  
)

# ______________________________________________________________________________
# 7. Save to file ----
# ______________________________________________________________________________    

saveRDS(constant.list, "Data for model/constants.rds")
saveRDS(data.list, "Data for model/data.rds")
