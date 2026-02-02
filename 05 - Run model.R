# PROJECT: Closed multi-session SCR
# SCRIPT: 05 - Run model
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 26 Jan 2026
# LAST MODIFIED: 02 Feb 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(nimble)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

constant.list <- readRDS("Data for model/constants.rds")
data.list <- readRDS("Data for model/data.rds")

# ______________________________________________________________________________
# 3. Code ----
# ______________________________________________________________________________

model.1.code <- nimbleCode({
  
  # priors
  # detection parameters
  # each one gets a random intercept by cluster
  
    # alpha0 - baseline detection (logit scale)
    # intercept
    alpha0_b0 ~ dunif(-4, 4)
    
    # covariate effects
    alpha0_b1 ~ dnorm(0, sd = 1)     # male effect
    alpha0_b2 ~ dnorm(0, sd = 1)     # XMC effect
    alpha0_b3 ~ dnorm(0, sd = 1)     # RET effect
    alpha0_b4 ~ dnorm(0, sd = 1)     # PIL effect

    # alpha2 - previous capture effect (normal scale)
    # intercept
    alpha2_b0 ~ dnorm(0, sd = 1)
    
    # covariate effects
    alpha2_b1 ~ dnorm(0, sd = 1)     # male effect
    alpha2_b2 ~ dnorm(0, sd = 1)     # XMC effect
    alpha2_b3 ~ dnorm(0, sd = 1)     # PIL effect
    alpha2_b4 ~ dnorm(0, sd = 1)     # PIL effect
  
    # sigma - exponential detection fn spatial scale of movement (log scale)
    # intercept
    sigma_b0 ~ dunif(log(10), log(70))    # Jensen et al. 2022 mean = log(45.27) = 3.81
  
    # covariate effects
    sigma_b1 ~ dnorm(0, sd = 1)     # male effect
    sigma_b2 ~ dnorm(0, sd = 1)     # XMC effect
    sigma_b3 ~ dnorm(0, sd = 1)     # PIL effect
    sigma_b4 ~ dnorm(0, sd = 1)     # PIL effect
  
  # data augmentation - indexed by session [G]
  for (g in 1:G) {
    
    psi[g] ~ dunif(0, 1)
    psi.sex[g] ~ dunif(0, 1)
    
  }
  
  # CALCULATE PROBABILITES FOR CATEGORICAL LIKELIHOOD
  # loop through individuals [M]
  for (i in 1:M) {
    
    # inclusion parameter (z)
    z[i] ~ dbern(psi[session[i]])
    
    # sex
    sex[i] ~ dbern(psi.sex[session[i]])
    
    # detection parameters
    # alpha0 - baseline detection (logit scale)
    alpha0[i] <- alpha0_b0 + 
                 alpha0_b1 * sex[i] +
                 alpha0_b2 * ft[i] +
                 alpha0_b3 * ret[i] +
                 alpha0_b4 * pil[i]
      
    # alpha2 - trap response
    alpha2[i] <- alpha2_b0 + 
                 alpha2_b1 * sex[i] +
                 alpha2_b2 * ft[i] +
                 alpha2_b3 * ret[i] +
                 alpha2_b4 * pil[i]
    
    # spatial scale of movement
    log(sigma[i]) <- sigma_b0 + 
                     sigma_b1 * sex[i] +
                     sigma_b2 * ft[i] +
                     sigma_b3 * ret[i] +
                     sigma_b4 * pil[i]
    
    alpha1[i] <- -1 / sigma[i]
    
    # activity centers (s) [M, 2] - indexed by site [U]
    s[i, 1] ~ dunif(S.lim[site[i], 1, 1], S.lim[site[i], 2, 1])
    s[i, 2] ~ dunif(S.lim[site[i], 1, 2], S.lim[site[i], 2, 2])
    
    # distances (d) between latent activity center and each trap
    # loop through traps [J] - indexed by site [U]
    for (j in 1:J) {
      
      d[i, j] <- sqrt(pow(s[i, 1] - trap.coords[j, 1, site[i]], 2) + 
                      pow(s[i, 2] - trap.coords[j, 2, site[i]], 2))
      
    }
    
    # loop through occasions [K] - indexed by session [G]
    for (k in 1:K[session[i]]) {
      
      # p0 - by-occasion detection probability
      logit(p0[i, k]) <- alpha0[i] + alpha2[i] * prev.cap[i, k]
      
      # and loop through traps [J]
      for (j in 1:J) {
        
        # calculate linear predictor (lp)
        lp[i, k, j] <- p0[i, k] * exp(
          
            # distance term
            alpha1[i] * d[i, j]
          
        ) *
          
          # inclusion
          z[i] *
          
          # trap operation 
          trap.op[j, k, session[i]]
        
        # probability
        p[i, k, j] <- lp[i, k, j] / (1 + sum(lp[i, k, 1:J]))  # sum over all traps
        
      }
      
      # probability of not being captured as the complement of all trap-specific probs
      p[i, k, J + 1] <- 1 - sum(p[i, k, 1:J])
      
    }
    
  }
  
  # likelihood for the observed individuals
  for (i in 1:n) {
    
    # loop through occasions K - indexed by session 
    for (k in 1:(K[session[i]])) {
      
      ch[i, k] ~ dcat(p[i, k, 1:(J + 1)])
      
    }
    
  }
  
  # likelihood for the non-detected individuals
  # inside, the probability of at AT LEAST ONE detection
  for (i in (n + 1):M) {
    
    zeroes[i] ~ dbern(1 - prod(1 - p[i, 1:(K[session[i]]), 1:J]))
    
  }
  
  # derived quantities
  # N
  for (g in 1:G) {
    
    for (i in 1:M) {
      
      N[i, g] <- z[i] * equals(session[i], g)
      
    }
    
    sum_N[g] <- sum(N[1:M , g])
    
    # D
    D[g] <- sum_N[g] / S.areas[site.g[g]]
    
  }
  
  # proportion D change
  # Y2 - Y1
  # untreated
  deltaD[1] <- (D[9] - D[1]) / D[1]    # 2A
  deltaD[2] <- (D[10] - D[2]) / D[2]   # 2B
  deltaD[3] <- (D[11] - D[3]) / D[3]   # 2C
  deltaD[4] <- (D[12] - D[4]) / D[4]   # 3A
  deltaD[5] <- (D[14] - D[5]) / D[5]   # 3C
  
  meanDeltaD[1] <- mean(deltaD[1:5])
  
  # Y3 - Y2
  # untreated
  deltaD[6] <- (D[20] - D[8]) / D[8]   # 1C
  deltaD[7] <- (D[23] - D[11]) / D[11] # 2C
  deltaD[8] <- (D[26] - D[14]) / D[14] # 3C
  deltaD[9] <- (D[29] - D[17]) / D[17] # 4C
  
  meanDeltaD[2] <- mean(deltaD[6:9])
  
  # retention 
  deltaD[10] <- (D[18] - D[6]) / D[6]   # 1A
  deltaD[11] <- (D[22] - D[10]) / D[10] # 2B
  deltaD[12] <- (D[25] - D[13]) / D[13] # 3B
  deltaD[13] <- (D[27] - D[15]) / D[15] # 4A
  
  meanDeltaD[3] <- mean(deltaD[10:13])
  
  # piling
  deltaD[14] <- (D[19] - D[7]) / D[7]   # 1B
  deltaD[15] <- (D[21] - D[9]) / D[9]   # 2A
  deltaD[16] <- (D[24] - D[12]) / D[12] # 3A
  deltaD[17] <- (D[28] - D[16]) / D[16] # 4B
  
  meanDeltaD[4] <- mean(deltaD[14:17])
  
  # Y4 - Y3
  # untreated
  deltaD[18] <- (D[32] - D[20]) / D[20] # 1C
  deltaD[19] <- (D[35] - D[23]) / D[23] # 2C
  deltaD[20] <- (D[38] - D[26]) / D[26] # 3C
  deltaD[21] <- (D[41] - D[29]) / D[29] # 4C
  
  meanDeltaD[5] <- mean(deltaD[18:21])
  
  # retention 
  deltaD[22] <- (D[30] - D[18]) / D[18] # 1A
  deltaD[23] <- (D[34] - D[22]) / D[22] # 2B
  deltaD[24] <- (D[37] - D[25]) / D[25] # 3B
  deltaD[25] <- (D[39] - D[27]) / D[27] # 4A
  
  meanDeltaD[6] <- mean(deltaD[22:25])
  
  # piling
  deltaD[26] <- (D[31] - D[19]) / D[19] # 1B
  deltaD[27] <- (D[33] - D[21]) / D[21] # 2A
  deltaD[28] <- (D[36] - D[24]) / D[24] # 3A
  deltaD[29] <- (D[40] - D[28]) / D[28] # 4B
  
  meanDeltaD[7] <- mean(deltaD[26:29])
  
})

# ______________________________________________________________________________
# 4. Initial values ----
# ______________________________________________________________________________

inits <- list(
  
  # initial s - all within the S bounds
  s = cbind(runif(constant.list$M, -200, 200),
            runif(constant.list$M, -200, 200)),
  
  # alpha0 - baseline detection (logit scale)
  # intercept
  alpha0_b0 = rnorm(1, 0, sd = 1), 
  
  # covariate effects
  alpha0_b1 = rnorm(1, 0, sd = 1),     
  alpha0_b2 = rnorm(1, 0, sd = 1),     
  alpha0_b3 = rnorm(1, 0, sd = 1),     
  alpha0_b4 = rnorm(1, 0, sd = 1),     
  
  # alpha2 - previous capture effect (normal scale)
  # intercept
  alpha2_b0 = rnorm(1, 0, sd = 1),
  
  # covariate effects
  alpha2_b1 = rnorm(1, 0, sd = 1),    
  alpha2_b2 = rnorm(1, 0, sd = 1),     
  alpha2_b3 = rnorm(1, 0, sd = 1),     
  alpha2_b4 = rnorm(1, 0, sd = 1),     
  
  # sigma - exponential detection fn spatial scale of movement (log scale)
  # intercept
  sigma_b0 = runif(1, log(10), log(70)), 
  
  # covariate effects
  sigma_b1 = rnorm(1, 0, sd = 1),     
  sigma_b2 = rnorm(1, 0, sd = 1),     
  sigma_b3 = rnorm(1, 0, sd = 1),     
  sigma_b4 = rnorm(1, 0, sd = 1),     
  
  # data augmentation
  psi = runif(41, 0, 1),
  psi.sex = runif(41, 0, 1),
  z = c(rep(NA, times = constant.list$n), rep(0, times = constant.list$M - constant.list$n)),
  
  # latent covariates
  sex = ifelse(is.na(data.list$sex) == F, NA, 0)
  
)

# ______________________________________________________________________________
# 5. Parameters to monitor ----
# ______________________________________________________________________________

monitor <- c(
  
  # detection
  # intercepts
    "alpha0_b0", "alpha2_b0", "sigma_b0",
  
    # slopes
    "alpha0_b1", "alpha0_b2", "alpha0_b3", "alpha0_b4",
    "alpha2_b1", "alpha2_b2", "alpha2_b3", "alpha2_b4",
    "sigma_b1", "sigma_b2", "sigma_b3", "sigma_b4",
  
  # inclusion
  "psi",
  
  # density
  "D",
  "meanDeltaD"
  
)

# ______________________________________________________________________________
# 6. Set up and run model----
# ______________________________________________________________________________

# set up model (this takes forever!)
model.1 <- nimbleModel(
  
  code = model.1.code,
  constants = constant.list,
  data = data.list,
  inits = inits,
  calculate = FALSE
  
)

# block samplers for sigma_b0 and alpha0_b0
model.1.conf <- configureMCMC(model.1)
model.1.conf$removeSamplers(c("alpha0_b0", "sigma_b0"))
model.1.conf$addSampler(c("alpha0_b0", "sigma_b0"), "RW_block")

# build model
model.1.mcmc <- buildMCMC(
  
  conf = model.1.conf,
  monitors = monitor
  
)

# compile model
compileNimble(model.1)
model.1.comp <- compileNimble(model.1.mcmc)

# run MCMC
model.1.run <- runMCMC(
  
  mcmc = model.1.comp,
  niter = 20000,
  nburnin = 10000,
  nchains = 1,
  thin = 10,
  samplesAsCodaMCMC = TRUE
  
)

saveRDS(model.1.run, "SCR_block.rds")

