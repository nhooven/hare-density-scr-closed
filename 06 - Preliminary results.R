# PROJECT: Closed multi-session SCR
# SCRIPT: 06 - Preliminary results
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 09 Feb 2026
# COMPLETED: 
# LAST MODIFIED: 09 Feb 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(coda)
library(bayestestR)

# ______________________________________________________________________________
# 2. Read in and bind each chain ----
# ______________________________________________________________________________

model.runs <- mcmc.list(list(readRDS("Runs/SCR_1.rds"),
                             readRDS("Runs/SCR_2.rds"),
                             readRDS("Runs/SCR_3.rds")))

model.runs.df <- as.data.frame(do.call(rbind, model.runs))

# ______________________________________________________________________________
# 3. Site-by-site density ----
# ______________________________________________________________________________
# 3a. Prepare for plotting ----
# ______________________________________________________________________________

# read in site_lookup
site.lookup <- read.csv("Data for model/site_lookup.csv")

# change trt column and add clusterID
site.lookup <- site.lookup %>%
  
  mutate(trt = case_when(siteID %in% c(3, 6, 9, 12) ~ "control",
                         siteID %in% c(1, 5, 8, 10) ~ "retention",
                         siteID %in% c(2, 4, 7, 11) ~ "piling"),
         clusterID = case_when(siteID %in% c(1:3) ~ 1,
                               siteID %in% c(4:6) ~ 2,
                               siteID %in% c(7:9) ~ 3,
                               siteID %in% c(10:12) ~ 4))

d.df <- model.runs.df %>%
  
  # keep only density posteriors
  dplyr::select("D[1]":"D[41]") %>%
  
  # pivot
  pivot_longer(cols ="D[1]":"D[41]",
               names_to = "sessionID",
               names_prefix = "D") %>%
  
  # remove brackets
  mutate(sessionID = as.integer(str_extract(sessionID, pattern = "\\d{1,2}"))) %>%
  
  # summarize with medians and CIs
  group_by(sessionID) %>%
  
  summarize(med = median(value),
            l50 = as.numeric(hdi(value, ci = 0.50)[2]),
            u50 = as.numeric(hdi(value, ci = 0.50)[3]),
            l95 = as.numeric(hdi(value, ci = 0.95)[2]),
            u95 = as.numeric(hdi(value, ci = 0.95)[3])) %>%
  
  # bind in site_lookup
  left_join(site.lookup) %>%
  
  # reorder factors and add life zone
  mutate(trt = factor(trt, levels = c("control", "retention", "piling")),
         lz = ifelse(clusterID == 4, "XMC", "SFL"))

# ______________________________________________________________________________
# 3b. Plot 1 ----

# gridded facets, CIs ONLY, no connecting lines

# ______________________________________________________________________________

ggplot(data = d.df) +
  
  theme_bw() +
  
  facet_grid(clusterID ~ trt) +
  
  # 95% CI
  geom_errorbar(aes(x = year,
                    ymin = l95,
                    ymax = u95,
                    color = lz),
                width = 0,
                linewidth = 2,
                alpha = 0.25) +
  
  # 50% CI
  geom_errorbar(aes(x = year,
                    ymin = l50,
                    ymax = u50,
                    color = lz),
                width = 0,
                linewidth = 2) +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 270,
                                   hjust = 0,
                                   vjust = 0.5,
                                   size = 8),
        strip.text.x = element_text(hjust = 0),
        strip.background.x = element_rect(fill = "white"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") +
  
  # axis title
  ylab("Density (hares/ha)") +
  
  # coordinates
  # let's add some padding on either side
  coord_cartesian(xlim = c(0.75, 4.25)) +
  
  # y axis
  scale_y_continuous(breaks = c(1, 3, 5)) +
  
  # x labels
  scale_x_continuous(breaks = seq(1, 4, 1),
                     labels = c("PRE 1", "PRE 2", "POST 1", "POST 2")) +
  
  # color
  scale_color_manual(values = c("#003300", "#669900"))

# 309 x 474

# ______________________________________________________________________________
# 3c. Plot 2 ----

# gridded facets, points and CIs, connecting lines

# ______________________________________________________________________________

ggplot(data = d.df) +
  
  theme_bw() +
  
  facet_grid(clusterID ~ trt) +
  
  # line
  geom_line(aes(x = year,
                y = med,
                color = lz),
            linewidth = 0.8) +
  
  # 95% CI
  geom_errorbar(aes(x = year,
                    ymin = l95,
                    ymax = u95,
                    color = lz),
                width = 0,
                linewidth = 1.25,
                alpha = 0.25) +
  
  # median
  geom_point(aes(x = year,
                 y = med,
                 color = lz),
             shape = 21,
             size = 1.25,
             stroke = 0.8,
             fill = "white") +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 270,
                                   hjust = 0,
                                   vjust = 0.5,
                                   size = 8),
        strip.text.x = element_text(hjust = 0),
        strip.background.x = element_rect(fill = "white"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") +
  
  # axis title
  ylab("Density (hares/ha)") +
  
  # coordinates
  # let's add some padding on either side
  coord_cartesian(xlim = c(0.75, 4.25)) +
  
  # y axis
  scale_y_continuous(breaks = c(1, 3, 5)) +
  
  # x labels
  scale_x_continuous(breaks = seq(1, 4, 1),
                     labels = c("PRE 1", "PRE 2", "POST 1", "POST 2")) +
  
  # color
  scale_color_manual(values = c("#003300", "#669900"))

# 309 x 474

# ______________________________________________________________________________
# 3d. Plot 3 ----

# all clusters together, facetted by trt, medians and 95% CIs, connecting lines

# ______________________________________________________________________________

ggplot(data = d.df) +
  
  theme_bw() +
  
  facet_wrap(~ trt, nrow = 1) +
  
  # line
  geom_line(aes(x = year,
                y = med,
                color = lz,
                group = clusterID),
            linewidth = 0.8,
            position = position_dodge(width = 0.5)) +
  
  # 95% CI
  geom_errorbar(aes(x = year,
                    ymin = l95,
                    ymax = u95,
                    color = lz,
                    group = clusterID),
                width = 0,
                linewidth = 1.25,
                alpha = 0.25,
                position = position_dodge(width = 0.5)) +
  
  # median
  geom_point(aes(x = year,
                 y = med,
                 color = lz,
                 group = clusterID),
             shape = 21,
             size = 1.25,
             stroke = 1,
             fill = "white",
             position = position_dodge(width = 0.5)) +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 270,
                                   hjust = 0,
                                   vjust = 0.5,
                                   size = 8),
        strip.text.x = element_text(hjust = 0),
        strip.background.x = element_rect(fill = "white"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") +
  
  # axis title
  ylab("Density (hares/ha)") +
  
  # coordinates
  # let's add some padding on either side
  coord_cartesian(xlim = c(0.75, 4.25)) +
  
  # y axis
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6)) +
  
  # x labels
  scale_x_continuous(breaks = seq(1, 4, 1),
                     labels = c("PRE 1", "PRE 2", "POST 1", "POST 2")) +
  
  # color
  scale_color_manual(values = c("#003300", "#669900"))

# 732 x 292
