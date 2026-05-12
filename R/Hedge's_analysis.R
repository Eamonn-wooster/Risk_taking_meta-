################### Hedges g effect size #####################



#################################################################################
######################### Set Up & Reading in Data ##############################
#################################################################################

# Loading packages 

pacman::p_load(devtools, 
               tidyverse, 
               metafor, 
               patchwork, 
               R.rsp, orchaRd,
               emmeans,
               metafor,
               orchaRd,
               stringr,
               dplyr,
               broom,
               tidyverse,
               here,
               broom,
               ggplot2,
               viridis,
               data.table,
               Matrix, # added
               matrixcalc,# added
               ape, # added
               multcomp 
           #added
               # miWQS #added
) 

# Setting WD

here()

# Reading in Data

dat <- read.csv(here("Data/Analysis_ready.csv")) 

#change a few specie names for analysis

dat$Species <- gsub("Equus_quagga","Equus_burchellii_quagga", dat$Species)
dat$Species <- gsub("Macropus_eugenii","Notamacropus_eugenii",  dat$Species)
dat$Species <- gsub("Hyla_versicolor", "Dryophytes_versicolor",  dat$Species)
dat$Species <- gsub("Bufo_americanus", "Anaxyrus_americanus",  dat$Species)
dat$Species <- gsub( "Gambusiaa_hubbsi","Gambusia_hubbsi",  dat$Species)
dat$Species <- gsub( "Gambusia_hubbsi ","Gambusia_hubbsi",  dat$Species)
dat$Species <- gsub( "Brachyraphis_episcopi","Brachyrhaphis_episcopi",  dat$Species)

# Reading in tree

tree1 <- read.tree(here("Tree/Tree_risk.tre")) 

# Getting branch length and correlation matrix

tree1b <- compute.brlen(tree1)
cor1 <-  vcv(tree1b, corr=T)

# checking the match 

setdiff(dat$Species, tree1$tip.label)
setdiff(tree1$tip.label, dat$Species)

#' [EW - This HAS to return character(0) for the model to run]

# creating non-phylo columns - has to be tree column 

dat$Species2 <- dat$Species

#################################################################################
######################### Converting Error ######################################
#################################################################################

### functions for converting error ------

# Function for converting CI's to SD

ci_to_sd <- function(bound, mean, n, ci_level = 0.95) {
  alpha <- 1 - ci_level
  z <- qnorm(1 - alpha / 2)
  me <- abs(bound - mean)
  sd <- me * sqrt(n) / z
  return(sd)
}

# Function for converting SE to SD

se_to_sd <- function(se, n) {
  # se: standard error
  # n: sample size
  # return: standard deviation
  sd_value <- se * sqrt(n)
  return(sd_value)
}

### Converting error types

unique(dat$Error_type)

### >>> Standard deviation --------

SD <- filter(dat, Error_type == "SD")

SD$SD_cont_trans <- SD$SD_control

SD$SD_exp_trans <- SD$SD_exp

### >>> 95% CI's --------

CI <- filter(dat, Error_type == "CI")

CI$SD_cont_trans <- ci_to_sd(CI$SD_control, CI$Mean_control, CI$n_control)

CI$SD_exp_trans <- ci_to_sd(CI$SD_exp, CI$Mean_exp, CI$n_exp)

### >>> Standard error --------

SE <- filter(dat, Error_type == "SE")

SE$SD_cont_trans <- se_to_sd(SE$SD_control, SE$n_control)

SE$SD_exp_trans <- se_to_sd(SE$SD_exp, SE$n_exp)

# Put them back together!!

data <- rbind(SD, CI, SE)

#################################################################################
############################ Cleaning Data ######################################
#################################################################################

#################################################################################
############################ Overall Effect Size ################################
#################################################################################

# Calculate Hedges' g

es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_exp,
  sd1i = SD_exp_trans,
  n1i = n_exp,
  m2i = Mean_control,
  sd2i = SD_cont_trans,
  n2i = n_control,
  data = data,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(es) 

# obs level random effect

es$Obs_ID <- factor(1:nrow(es))

# Study level random effect # using title for now but will extract extra data for first author/year later 

# es$StudyID <- NULL
# 
# es$Study_ID <- factor(as.numeric(as.factor(es$Firstauthor_Year)))


################################## HEDges Models ##################################

# Variance covariance matrix

vcv <- vcalc(vi, cluster = Title, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             data = es) 

#Overall model - we need this model - important to overall result and for pub bias analysis


mod.overall <- rma.mv(yi = yi, V = vcv,
                             random = list(~1 | Title,
                                           ~1 | Species, # phylo effect 
                                           ~1 | Species2, # non-phylo effect 
                                           ~1 | Obs_ID), 
                             data =  es,
                             # control = list(optimizer="BFGS"),
                             test = "t",
                             sparse = TRUE,
                             R = list(Species = cor1))

  summary(mod.overall)


I2 = round(i2_ml(mod.overall), 2)

overall <- orchard_plot(mod.overall, xlab = "Change in risk-taking (Hedge's g)", group = "Title",
                        angle = 0)
overall

ggsave(here("Outputs", "Figures", "Overall.pdf"), width = 6, height = 4)

#This is the hetrogeneity of the effect sizes - we need to report this in the manuscript 

round(i2_ml(mod.overall), 2)

##Behaviour type####

unique(data$Behaviour_type)

mod.behav <- rma.mv(yi = yi, V = vcv,
                      random = list(~1 | Title,
                                    ~1 | Species, # phylo effect 
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), 
                      data =  es,
                      mods = ~ Behaviour_type -1,
                      test = "t",
                      sparse = TRUE,
                      R = list(Species = cor1))

summary(mod.behav)

#Sex##

mod.sex <- rma.mv(yi = yi, V = vcv,
                    random = list(~1 | Title,
                                  ~1 | Species, # phylo effect 
                                  ~1 | Species2, # non-phylo effect 
                                  ~1 | Obs_ID), 
                    data =  es,
                    mods = ~ Sex -1,
                    test = "t",
                    sparse = TRUE,
                    R = list(Species = cor1))

summary(mod.sex)


#Pred_pressure#######

mod.pp <- rma.mv(yi = yi, V = vcv,
                  random = list(~1 | Title,
                                ~1 | Species, # phylo effect 
                                ~1 | Species2, # non-phylo effect 
                                ~1 | Obs_ID), 
                  data =  es,
                  mods = ~ Predation_pressure -1,
                  test = "t",
                  sparse = TRUE,
                  R = list(Species = cor1))

summary(mod.pp)


#real vs simulated pred#######

#convert es$real_predator to a factor
es$Real_predator <- as.factor(es$Real_predator)

str(es)


mod.rp <- rma.mv(yi = yi, V = vcv,
                 random = list(~1 | Title,
                               ~1 | Species, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  es,
                 mods = ~ Real_predator -1,
                 test = "t",
                 sparse = TRUE,
                 R = list(Species = cor1))

summary(mod.rp)

#low vs no predator###

mod.ln <- rma.mv(yi = yi, V = vcv,
                 random = list(~1 | Title,
                               ~1 | Species, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  es,
                 mods = ~ Low_or_no_pred -1,
                 test = "t",
                 sparse = TRUE,
                 R = list(Species = cor1))

summary(mod.ln)



##############lnVR models##################

# Calculate lnVR

cv <- escalc(
  measure = "VR",       # standardized mean difference
  m1i = Mean_exp,
  sd1i = SD_exp_trans,
  n1i = n_exp,
  m2i = Mean_control,
  sd2i = SD_cont_trans,
  n2i = n_control,
  data = data
)

cv$Obs_ID <- factor(1:nrow(cv))

vcv_cv <- vcalc(vi, cluster = Title, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             data = es) 


mod.overall_cv <- rma.mv(yi = yi, V = vcv_cv,
                      random = list(~1 | Title,
                                    ~1 | Species, # phylo effect 
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), 
                      data =  cv,
                      # control = list(optimizer="BFGS"),
                      test = "t",
                      sparse = TRUE,
                      R = list(Species = cor1))

summary(mod.overall_cv)


I2 = round(i2_ml(mod.overall_cv), 2)

overall_cv <- orchard_plot(mod.overall_cv, xlab = "Heterogeneity in risk-taking (log transformed variability ratio)", group = "Title",
                        angle = 0)
overall_cv

#CV_ behaviour type####

mod.behav_cv <- rma.mv(yi = yi, V = vcv_cv,
                    random = list(~1 | Title,
                                  ~1 | Species, # phylo effect 
                                  ~1 | Species2, # non-phylo effect 
                                  ~1 | Obs_ID), 
                    data =  cv,
                    mods = ~ Behaviour_type -1,
                    test = "t",
                    sparse = TRUE,
                    R = list(Species = cor1))

summary(mod.behav_cv)

#CV_ Sex##


mod.sex.cv <- rma.mv(yi = yi, V = vcv_cv,
                  random = list(~1 | Title,
                                ~1 | Species, # phylo effect 
                                ~1 | Species2, # non-phylo effect 
                                ~1 | Obs_ID), 
                  data =  cv,
                  mods = ~ Sex -1,
                  test = "t",
                  sparse = TRUE,
                  R = list(Species = cor1))

summary(mod.sex)

#CV_real vs sim#######

cv$Real_predator <- as.factor(cv$Real_predator)

str(cv)


mod.rp.cv <- rma.mv(yi = yi, V = vcv_cv,
                 random = list(~1 | Title,
                               ~1 | Species, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  cv,
                 mods = ~ Real_predator -1,
                 test = "t",
                 sparse = TRUE,
                 R = list(Species = cor1))

summary(mod.rp.cv)

#CV low vs no predator###

mod.ln.cv <- rma.mv(yi = yi, V = vcv_cv,
                 random = list(~1 | Title,
                               ~1 | Species, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  cv,
                 mods = ~ Low_or_no_pred -1,
                 test = "t",
                 sparse = TRUE,
                 R = list(Species = cor1))

summary(mod.ln.cv)

