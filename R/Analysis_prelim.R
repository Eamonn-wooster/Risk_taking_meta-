################### Hedges g effect size #####################



#################################################################################
######################### Set Up & Reading in Data ##############################
#################################################################################

# Loading packages 

rm(list = ls())

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
               multcomp,
               taxize
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

# obs level random effect

data$Obs_ID <- factor(1:nrow(data))

# Study level random effect # using Study_ID for now but will extract extra data for first author/year later 

data$Firstauthor_Year <- paste(data$First_author, data$Year, sep = "_")

data$Study_ID <- factor(as.numeric(as.factor(data$Firstauthor_Year)))


#################################################################################
############################ Hedges Effect Size ################################
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




################################## HEDges Models ##################################

# Variance covariance matrix

vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             data = es) 

#Overall model - we need this model - important to overall result and for pub bias analysis


mod.overall <- rma.mv(yi = yi, V = vcv,
                             random = list(~1 | Study_ID,
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

overall <- orchard_plot(mod.overall, xlab = "Change in risk-taking (Hedge's g)", group = "Study_ID",
                        angle = 0)
overall

ggsave(here("Outputs", "Figures", "Overall.pdf"), width = 6, height = 4)

#This is the hetrogeneity of the effect sizes - we need to report this in the manuscript 

round(i2_ml(mod.overall), 2)

##Behaviour type####

unique(data$Behaviour_type)

mod.behav <- rma.mv(yi = yi, V = vcv,
                      random = list(~1 | Study_ID,
                                    ~1 | Species, # phylo effect 
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), 
                      data =  es,
                      mods = ~ Behaviour_type -1,
                      test = "t",
                      sparse = TRUE,
                      R = list(Species = cor1))

summary(mod.behav)

#plotting behaviour type with orchard
behav <- orchard_plot(mod.behav, mod = "Behaviour_type", xlab = "Change in risk-taking (Hedge's g)", group = "Study_ID",
                        angle = 0)

behav



#Sex##

ses <- filter(es, Sex != "N/A" & Sex != "Mixed")

vcv_ses <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             data = ses)

mod.sex <- rma.mv(yi = yi, V = vcv_ses,
                    random = list(~1 | Study_ID,
                                  ~1 | Species, # phylo effect 
                                  ~1 | Species2, # non-phylo effect 
                                  ~1 | Obs_ID), 
                    data =  ses,
                    mods = ~ Sex -1,
                    test = "t",
                    sparse = TRUE,
                    R = list(Species = cor1))

summary(mod.sex)


#Pred_pressure#######

mod.pp <- rma.mv(yi = yi, V = vcv,
                  random = list(~1 | Study_ID,
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

es$Real_predator <- as.factor(es$Real_predator)

str(es)


mod.rp <- rma.mv(yi = yi, V = vcv,
                 random = list(~1 | Study_ID,
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
                 random = list(~1 | Study_ID,
                               ~1 | Species, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  es,
                 mods = ~ Low_or_no_pred -1,
                 test = "t",
                 sparse = TRUE,
                 R = list(Species = cor1))

summary(mod.ln)

#adult vs juvenille###

mod.aj <- rma.mv(yi = yi, V = vcv,
                 random = list(~1 | Study_ID,
                               ~1 | Species, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  es,
                 mods = ~ Adult -1,
                 test = "t",
                 sparse = TRUE,
                 R = list(Species = cor1))

summary(mod.aj)

#comp type###

mod.com <- rma.mv(yi = yi, V = vcv,
                 random = list(~1 | Study_ID,
                               ~1 | Species, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  es,
                 mods = ~ Comparison_type -1,
                 test = "t",
                 sparse = TRUE,
                 R = list(Species = cor1))

summary(mod.com)

###Class##########

unique(es$Class)

mod.class <- rma.mv(yi = yi, V = vcv,
                  random = list(~1 | Study_ID,
                                ~1 | Species, # phylo effect 
                                ~1 | Species2, # non-phylo effect 
                                ~1 | Obs_ID), 
                  data =  es,
                  mods = ~ Class -1,
                  test = "t",
                  sparse = TRUE,
                  R = list(Species = cor1))

summary(mod.class)

class <- orchard_plot(mod.class, mod = "Class", xlab = "Change in risk-taking (Hedge's g)", group = "Study_ID",
                      angle = 0)

class


###publication bias####### only doing mean shifts bc thats what we extracted for 

####Eggers regression - significant intercept - evidence of pub bias

es$effectN <- (es$n_control * es$n_exp) / (es$n_control + es$n_exp)
es$sqeffectN <- sqrt(es$effectN)

mod.egg <- rma.mv(yi = yi, V = vcv,
                  random = list(~1 | Study_ID,
                                ~1 | Species, # phylo effect 
                                ~1 | Species2, # non-phylo effect 
                                ~1 | Obs_ID), 
                  data =  es,
                  mods = ~ sqeffectN,
                  test = "t",
                  sparse = TRUE,
                  R = list(Species = cor1))

summary(mod.egg)

#exploring pub bias 

mod_simple <- rma(yi = yi, vi = vi, data = es)

mod_simple
tf <- trimfill(mod_simple)

range(es$yi)
range(sqrt(es$vi)) 

par(mar = c(5, 5, 4, 2))  # bottom, left, top, right margins
funnel(mod.egg)

##no change with trim and fill. Likely means the funnel assymerty is due to heterogeneity rather than pub bias but we should report this in the manuscript.

#Decline effect 

mod.mad <- rma.mv(yi = yi, V = vcv,
                          random = list(~1 | Study_ID,
                                        ~1 | Species, # phylo effect 
                                        ~1 | Species2, # non-phylo effect 
                                        ~1 | Obs_ID), 
                          data =  es,
                          mods = ~ Year,
                          test = "t",
                          sparse = TRUE,
                          R = list(Species = cor1))


summary(mod.mad) # no effect

# DO LEAVE ONE OUT HERE #######



###########################################
##############lnVR models##################
###########################################


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

vcv_cv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             data = es) 


mod.overall_cv <- rma.mv(yi = yi, V = vcv_cv,
                      random = list(~1 | Study_ID,
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

overall_cv <- orchard_plot(mod.overall_cv, xlab = "Heterogeneity in risk-taking (log transformed variability ratio)", group = "Study_ID",
                        angle = 0)
overall_cv

#CV_ behaviour type####

mod.behav_cv <- rma.mv(yi = yi, V = vcv_cv,
                    random = list(~1 | Study_ID,
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

ses_cv <- filter(cv, Sex != "N/A" & Sex != "Mixed")

vcv_ses_cv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                 data = ses_cv)


mod.sex.cv <- rma.mv(yi = yi, V = vcv_ses_cv,
                  random = list(~1 | Study_ID,
                                ~1 | Species, # phylo effect 
                                ~1 | Species2, # non-phylo effect 
                                ~1 | Obs_ID), 
                  data =  ses_cv,
                  mods = ~ Sex -1,
                  test = "t",
                  sparse = TRUE,
                  R = list(Species = cor1))

summary(mod.sex.cv) #mean-variance relationship in males only

#CV_real vs sim#######

cv$Real_predator <- as.factor(cv$Real_predator)

str(cv)


mod.rp.cv <- rma.mv(yi = yi, V = vcv_cv,
                 random = list(~1 | Study_ID,
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
                 random = list(~1 | Study_ID,
                               ~1 | Species, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  cv,
                 mods = ~ Low_or_no_pred -1,
                 test = "t",
                 sparse = TRUE,
                 R = list(Species = cor1))

summary(mod.ln.cv)

#CV adult vs juvenille###

mod.aj.cv <- rma.mv(yi = yi, V = vcv_cv,
                 random = list(~1 | Study_ID,
                               ~1 | Species, # phylo effect 
                               ~1 | Species2, # non-phylo effect 
                               ~1 | Obs_ID), 
                 data =  cv,
                 mods = ~ Adult -1,
                 test = "t",
                 sparse = TRUE,
                 R = list(Species = cor1))

summary(mod.aj.cv)

#comp type###

mod.com.cv <- rma.mv(yi = yi, V = vcv_cv,
                  random = list(~1 | Study_ID,
                                ~1 | Species, # phylo effect 
                                ~1 | Species2, # non-phylo effect 
                                ~1 | Obs_ID), 
                  data =  cv,
                  mods = ~ Comparison_type -1,
                  test = "t",
                  sparse = TRUE,
                  R = list(Species = cor1))

summary(mod.com.cv)




###END##### 

#Scribles below 

#time exposed### - not sure if this metric makes any sense nor how to test it in a meaningful way....

# #filter out non numbers
# unique(es$time_exposed_days)
# es_time <- es %>% filter(time_exposed_days != "N/A" & time_exposed_days != "Max")
# 
# unique(es_time$time_exposed_days)
# #convert time_exposed_days to numeric
# 
# es_time$time_exposed_days <- as.numeric(es_time$time_exposed_days)
# 
# str(es_time$time_exposed_days)
# 
# #log transform time_exposed_days
# es_time$log_time_exposed_days <- log(es_time$time_exposed_days)
# 
# vcv_time <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
#              data = es_time)
# 
# mod.time <- rma.mv(yi = yi, V = vcv_time,
#                   random = list(~1 | Study_ID,
#                                 ~1 | Species, # phylo effect 
#                                 ~1 | Species2, # non-phylo effect 
#                                 ~1 | Obs_ID), 
#                   data =  es_time,
#                   mods = ~ log_time_exposed_days,
#                   test = "t",
#                   sparse = TRUE,
#                   R = list(Species = cor1))
# 
# summary(mod.time)
# 
# #plot it
# 
# newdat_rat_time <- data.frame(time_exp = seq(min(es_time$log_time_exposed_days, na.rm = TRUE),
#                                               max(es_time$log_time_exposed_days, na.rm = TRUE),
#                                               length.out = 100))
# 
# X_rat_time <- model.matrix(~ time_exp, data = newdat_rat_time)[, -1, drop = FALSE]
# 
# time_rat_time <- predict(mod.time, newmods = X_rat_time)
# 
# plot_data_time <- cbind(newdat_rat_time, time_rat_time)
# 
# setDT(es_time)
# es_time[, wi := 1/sqrt(vi)]
# es_time[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
# es_time
# 
# time_text <- tidy(mod.time)
# time_text
# 
# bottom_margin <- margin(5, 5, 30, 5)
# 
# time_exp <- ggplot()+
#   # now add the rest of the points:
#   geom_jitter(data = es_time, 
#               aes(x = log_time_exposed_days, size = pt_size, #ignore point size
#                   y = yi),
#               position = position_jitter(width = 0.01),
#               inherit.aes = F,  alpha = 0.5, color = "#3d405bff")+
#   geom_hline(yintercept = 0, lty = "dashed")+
#   geom_ribbon(data = plot_data_time,
#               aes(x = time_exp,
#                   ymin = pi.lb, ymax = pi.ub),
#               fill = "transparent", color = "#dad7cd",
#               alpha = .3)+
#   geom_ribbon(data = plot_data_time, 
#               aes(x = time_exp,
#                   ymin = ci.lb, ymax = ci.ub),
#               alpha = .3)+
#   geom_line(data = plot_data_time, 
#             aes(x = time_exp,
#                 y = pred),
#             alpha = .8)+
#   
#   # coord_cartesian(ylim = c(-4, 4))+
#   scale_size_identity()+
#   theme_bw()+
#   theme(
#     plot.margin = bottom_margin)+
#   xlab("Time exposed to predators (log)")+
#   ylab("Risk taking behaviour (Hedge's g)")+
#   theme(legend.position = "none",
#         text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
#         panel.grid = element_blank(),
#         panel.border = element_blank())
# 
# 
# time_exp

