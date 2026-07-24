################### Hedges g effect size #####################




######################### Set Up & Reading in Data ##############################


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
               taxize,
               cowplot
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


######################### Converting Error ######################################


### functions for converting error ------

# Function for converting CI's to SD

me_to_sd <- function(me, n, ci_level = 0.95) {
  alpha <- 1 - ci_level
  tcrit <- qt(1 - alpha/2, df = n - 1)
  me * sqrt(n) / tcrit
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

CI$SD_cont_trans <- me_to_sd(CI$SD_control, CI$n_control)

CI$SD_exp_trans <- me_to_sd(CI$SD_exp, CI$n_exp)

### >>> Standard error --------

SE <- filter(dat, Error_type == "SE")

SE$SD_cont_trans <- se_to_sd(SE$SD_control, SE$n_control)

SE$SD_exp_trans <- se_to_sd(SE$SD_exp, SE$n_exp)

# Put them back together!!

data <- rbind(SD, CI, SE)

# obs level random effect

data$Obs_ID <- factor(1:nrow(data))

# Study level random effect 

data$Firstauthor_Year <- paste(data$First_author, data$Year, sep = "_")

data$Study_ID <- factor(as.numeric(as.factor(data$Firstauthor_Year)))

###Sort of the negative issues

inv <- filter(data, Sign_inverted == "Yes")

not <- filter(data, Sign_inverted == "No")

#Ok lets flip the sign back for inv

inv$Mean_control <- -inv$Mean_control

inv$Mean_exp <- -inv$Mean_exp

data2 <- rbind(inv, not)

############################ Hedges Effect Size ################################

# Calculate Hedges' g

es <- escalc(
  measure = "SMDH",     
  m1i = Mean_exp,
  sd1i = SD_exp_trans,
  n1i = n_exp,
  m2i = Mean_control,
  sd2i = SD_cont_trans,
  n2i = n_control,
  data = data2
)



# Inspect the results

print(es) 

##now we need to reinvert the effect sizes the sign on "yes" on sign-inverted

es_inv <- filter(es, Sign_inverted == "Yes")

es_not <- filter(es, Sign_inverted == "No")

#Ok lets flip the sign back for inv

es_inv$yi <- -es_inv$yi


es2 <- rbind(es_inv, es_not)


################################## Hedges Models ##################################

# Variance covariance matrix

vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             data = es) 


mod.overall.vcv <- rma.mv(yi = yi, V = vcv, #vcv, #' [EJL changed]
                      random = list(#~1 | Study_ID / Obs_ID, #' [EJL changed]
                                    ~1 | Species, # phylo effect 
                                    ~1 | Species2#, # non-phylo effect 
                      ), #' [EJL changed]
                      data =  es2,
                      # control = list(optimizer="BFGS"),
                      test = "t",
                      sparse = TRUE,
                      R = list(Species = cor1))

summary(mod.overall.vcv)


mod.overall <- rma.mv(yi = yi, V = vi, #vcv, #' [EJL changed]
                             random = list(~1 | Study_ID / Obs_ID, #' [EJL changed]
                                           ~1 | Species, # phylo effect 
                                           ~1 | Species2#, # non-phylo effect 
                                           ), #' [EJL changed]
                             data =  es,
                             # control = list(optimizer="BFGS"),
                             test = "t",
                             sparse = TRUE,
                             R = list(Species = cor1))

summary(mod.overall)
AIC(mod.overall.vcv, mod.overall)
# Big improvement without VCV

# >>> Choose optimal random effect (EJL) ----------------------------------
#' [EJL:] Sigmas are super low for some of the levels, which could bias your estimates.
#' I would compare to simpler models before reporting.

mod.overall.simple1 <- rma.mv(yi = yi, V = vi, #vcv, #' [EJL changed]
                              random = list(~1 | Study_ID / Obs_ID, #' [EJL changed]
                                            #~1 | Species, # phylo effect 
                                            ~1 | Species2#, # non-phylo effect 
                              ), #' [EJL changed]
                              data =  es2,
                              # control = list(optimizer="BFGS"),
                              test = "t",
                              sparse = TRUE)

summary(mod.overall.simple1)
AIC(mod.overall, mod.overall.simple1)
# Improved


mod.overall.simple2 <- rma.mv(yi = yi, V = vi, #vcv, #' [EJL changed]
                              random = list(~1 | Study_ID / Obs_ID#, #' [EJL changed]
                                            #~1 | Species, # phylo effect 
                                            #~1 | Species2#, # non-phylo effect 
                              ), #' [EJL changed]
                              data =  es,
                              # control = list(optimizer="BFGS"),
                              test = "t",
                              sparse = TRUE)

summary(mod.overall.simple2)
# Identical, but slightly smaller, so i guess either way haha
AIC(mod.overall.simple1, mod.overall.simple2)
# Improved slightly.


# Copy best model:
mod.overall <- mod.overall.simple2

# >>> Final model summaries/plot ---------------------------------------------------------


I2 = round(i2_ml(mod.overall), 2)

I2

overall <- orchard_plot(mod.overall, xlab = "Difference in risk-taking (Hedge's g)", group = "Study_ID",
                        angle = 0, col = "#eea196")
#' [EJL] This produced an error for me:
# Error in if (colour) as.factor(data_trim$stdy) else data_trim$moderator : 
# argument is not interpretable as logical

overall <- overall + scale_fill_manual(values = "#eea196") + scale_color_manual(values = "#eea196") 

overall

ggsave(here("Outputs", "Figures", "Overall.pdf"), width = 6, height = 4)

#This is the hetrogeneity of the effect sizes - we need to report this in the manuscript 

round(i2_ml(mod.overall), 2)

##Behaviour type####

unique(data$Behaviour_type)

es2$Behaviour_type <- factor(
  es2$Behaviour_type,
  levels = sort(unique(es2$Behaviour_type), decreasing = TRUE)
)


mod.behav <- rma.mv(yi = yi, V = vi,
                      random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                                    # ~1 | Species, # phylo effect 
                                    # ~1 | Species2), 
                      data =  es2,
                      mods = ~ Behaviour_type -1,
                      test = "t",
                      sparse = TRUE)
                      # R = list(Species = cor1))


summary(mod.behav)


behav <- orchard_plot(mod.behav, mod = "Behaviour_type", 
                      xlab = "Difference in risk-taking behaviour (Hedge's g)", 
                      group = "Study_ID",
                      angle = 0) +
  scale_fill_manual(values = rep("#eea196", 7)) +
  scale_color_manual(values = rep("#eea196", 7))

behav



#Sex#####

ses <- filter(es2, Sex != "N/A")

mod.sex <- rma.mv(yi = yi, V = vi,
                  random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                  # ~1 | Species, # phylo effect 
                  # ~1 | Species2),
                    data =  ses,
                    mods = ~ Sex -1,
                    test = "t",
                    sparse = TRUE)

summary(mod.sex)



#Pred_pressure#######

mod.pp <- rma.mv(yi = yi, V = vi,
                 random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                 # ~1 | Species, # phylo effect 
                 # ~1 | Species2),
                  data =  es2,
                  mods = ~ Predation_pressure -1,
                  test = "t",
                  sparse = TRUE)

summary(mod.pp)


#real vs simulated pred#######

es2$Real_predator <- as.factor(es2$Real_predator)

str(es2)


mod.rp <- rma.mv(yi = yi, V = vi,
                 random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                 # ~1 | Species, # phylo effect 
                 # ~1 | Species2), 
                 data =  es2, 
                 mods = ~ Real_predator -1,
                 test = "t",
                 sparse = TRUE)

summary(mod.rp)

anova(mod.rp, L = c(1, -1))

###low vs no predator######

mod.ln <- rma.mv(yi = yi, V = vi,
                 random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                 # ~1 | Species, # phylo effect 
                 # ~1 | Species2), 
                 data =  es2, 
                 mods = ~ Low_or_no_pred -1,
                 test = "t",
                 sparse = TRUE)

summary(mod.ln)

###adult vs juvenille#####

es2$Adult <- factor(es2$Adult, 
                              levels = c("Adult", "Juvenille", "Both"))

mod.aj <- rma.mv(yi = yi, V = vi,
                 random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                 # ~1 | Species, # phylo effect 
                 # ~1 | Species2), 
                 data =  es2, 
                 mods = ~ Adult -1,
                 test = "t",
                 sparse = TRUE)

summary(mod.aj)

###comp type######

mod.com <- rma.mv(yi = yi, V = vi,
                  random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                  # ~1 | Species, # phylo effect 
                  # ~1 | Species2), 
                  data =  es2, 
                 mods = ~ Comparison_type -1,
                 test = "t",
                 sparse = TRUE)

summary(mod.com)

###Class##########

unique(es2$Class)


es2 <- es2 %>%
  mutate(Class = recode(Class,
                        "Actinopterygii" = "Ray-finned fishes",
                        "Amphibia" = "Amphibians",
                        "Mammalia" = "Mammals",
                        "Reptilia" = "Reptiles",
                        "Malacostraca" = "Crustaceans",
                        "Chondrichthyes" = "Cartilaginous fishes",
                        "Aves" = "Birds"
  ))

es2$Class <- factor(
  es2$Class,
  levels = sort(unique(es2$Class), decreasing = TRUE)
)


mod.class <- rma.mv(yi = yi, V = vi,
                    random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                    # ~1 | Species, # phylo effect 
                    # ~1 | Species2),  
                  data =  es2,
                  mods = ~ Class -1,
                  test = "t",
                  sparse = TRUE)
                

summary(mod.class)

class <- orchard_plot(mod.class, mod = "Class", xlab = "Difference in risk-taking (Hedge's g)", group = "Study_ID",
                      angle = 0) + 
  scale_fill_manual(values = rep("#eea196", 7)) +
  scale_color_manual(values = rep("#eea196", 7))
class   


###Publication bias####### 

####Eggers regression - significant intercept - evidence of pub bias

es2$effectN <- (es2$n_control * es2$n_exp) / (es2$n_control + es2$n_exp)
es2$sqeffectN <- sqrt(es2$effectN)

#' [EJL Changed:]
mod.egg <- rma.mv(yi = yi, V = vi,
                  random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                  # ~1 | Species, # phylo effect 
                  # ~1 | Species2),  
                  data =  es2,
                  mods = ~ sqeffectN,
                  test = "t",
                  sparse = TRUE)

summary(mod.egg)

#exploring pub bias 

mod_simple <- rma(yi = yi, vi = vi, data = es2)

mod_simple
tf <- trimfill(mod_simple)
tf

range(es$yi)
range(sqrt(es$vi)) 

par(mar = c(5, 5, 4, 2))  # bottom, left, top, right margins
funnel(mod.egg)

##no change with trim and fill. Likely means the funnel assymerty is due to heterogeneity rather than pub bias but we should report this in the manuscript.

###Decline effect #####
#' [EJL Changed:]
mod.mad <- rma.mv(yi = yi, V = vi,
                  random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                  # ~1 | Species, # phylo effect 
                  # ~1 | Species2),  
                  data =  es2,
                          mods = ~ Year,
                          test = "t",
                          sparse = TRUE)


summary(mod.mad) # no effect

###Leave on out####
#' [EJL: I think cooks.distance() on the model object basically does LOO (and code is one line..It'll give you a value that indicates how much estimates changed with each study)]
#' [If you specify strata I think]

cook.out <- cooks.distance(model = mod.overall,
                           cluster = Study_ID)
# There are a few different thresholds but I've seen a lot that exclude studies with cook > 4/N studies
4 / length(unique(cook.out))

cook.out[cook.out > 4 / length(unique(cook.out))]
# Which is only study 42. So you could rerun those models without study 42.



#' [Back to original code:]
dat <- es %>%
  mutate(leave_out = paste(First_author, Year, sep = "_"))
dat$leave_out <- as.factor(dat$leave_out)

rerun <- F
if(rerun){
  LeaveOneOut_effectsize <- list()
  for (i in 1:length(levels(dat$leave_out))) {
    temp_dat <- dat %>%
      dplyr::filter(leave_out != levels(dat$leave_out)[i])
    
    VCV_leaveout <- vcalc(vi = temp_dat$vi, cluster = temp_dat$Study_ID, obs = temp_dat$Obs_ID, rho = 0.5)
    
    LeaveOneOut_effectsize[[i]] <- rma.mv(yi = yi,
                                          V = VCV_leaveout,
                                          random = list(~1 | Study_ID,
                                                        ~1 | Species,   # phylo effect
                                                        ~1 | Species2,  # non-phylo effect
                                                        ~1 | Obs_ID),
                                          R = list(Species = cor1),
                                          test = "t",
                                          method = "REML",
                                          sparse = TRUE,
                                          data = temp_dat)
  }
  
  # function for extracting est, ci.lb, and ci.ub from all models
  est.func <- function(model) {
    df <- data.frame(est = model$b, lower = model$ci.lb, upper = model$ci.ub)
    return(df)
  }
  
  # form data frame
  MA_oneout <- lapply(LeaveOneOut_effectsize, function(x) est.func(x)) %>%
    bind_rows() %>%
    mutate(left_out = levels(dat$leave_out))
  
  # preserve factor order for ggplot
  MA_oneout$left_out <- factor(MA_oneout$left_out, levels = MA_oneout$left_out)
  
  # save the runs
  saveRDS(MA_oneout, here("R", "MA_oneout.RDS"))
  
} else {
  MA_oneout <- readRDS(here("R", "MA_oneout.RDS"))
}

# plotting
leaveoneout <- ggplot(MA_oneout) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1) +
  geom_hline(yintercept = mod.overall$ci.lb, lty = 3, lwd = 0.75, colour = "black") +
  geom_hline(yintercept = mod.overall$b,     lty = 1, lwd = 0.75, colour = "black") +
  geom_hline(yintercept = mod.overall$ci.ub, lty = 3, lwd = 0.75, colour = "black") +
  geom_pointrange(aes(x = left_out, y = est, ymin = lower, ymax = upper)) +
  xlab("Study left out") +
  ylab("Hedge' g") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 6)
  )

leaveoneout

###Sens analysis for lncvr####

sens.es <- filter(es2, Mean_control > 0)

sens.es <- filter(sens.es, Mean_exp > 0)

sens.es <- filter(sens.es, Can_be_negitive != "Yes")

sens.es$Obs_ID <- factor(1:nrow(sens.es))

mod.behav.sens <- rma.mv(yi = yi, V = vi,
                         random = list(~1 | Study_ID / Obs_ID), #' [EJL Change]
                         # ~1 | Species, # phylo effect 
                         # ~1 | Species2),  
                    data =  sens.es,
                    test = "t",
                    sparse = TRUE)

summary(mod.behav.sens)

#removing those studies does NOT change the overall results!!!


##############lncvr models##################

#removing negatives from the dataset
cv <- filter(data2, Mean_control > 0)

cv <- filter(cv, Mean_exp > 0)

#remove methods where a value could be negative 

cv <- filter(cv, Can_be_negitive != "Yes")


# Calculate lncvr

cv <- escalc(
  measure = "CVR",     
  m1i = Mean_exp,
  sd1i = SD_exp_trans,
  n1i = n_exp,
  m2i = Mean_control,
  sd2i = SD_cont_trans,
  n2i = n_control,
  data = cv
)

cv$Obs_ID <- factor(1:nrow(cv))

vcv_cv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                data = cv) 


mod.overall_cv <- rma.mv(yi = yi, V = vi,
                         random = list(~1 | Study_ID / Obs_ID), 
                         data =  cv,
                         # control = list(optimizer="BFGS"),
                         test = "t",
                         sparse = TRUE)

mod.overall_cv2 <- rma.mv(yi = yi, V = vcv_cv,
                       random = list(~1 | Study_ID,
                                     ~1 | Species, # phylo effect 
                                     ~1 | Species2, # non-phylo effect 
                                     ~1 | Obs_ID), 
                       data =  cv,
                       test = "t",
                       sparse = TRUE,
                       R = list(Species = cor1))

AIC(mod.overall_cv, mod.overall_cv2)

mod.overall_cv3 <- rma.mv(yi = yi, V = vi,
                         random = list(~1 | Study_ID / Obs_ID,
                                       ~1 | Species2), 
                         data =  cv,
                         test = "t",
                         sparse = TRUE)

AIC(mod.overall_cv, mod.overall_cv3)

# Copy best model:
mod.overall_cv <- mod.overall_cv3

summary(mod.overall_cv)


overall_cv <- orchard_plot(mod.overall_cv, xlab = "Heterogeneity in risk-taking (lnCVR)", group = "Study_ID",
                           angle = 0)
overall_cv <- overall_cv + scale_fill_manual(values = "#989aae") + scale_color_manual(values = "#989aae") 

overall_cv

#CV_ behaviour type####

cv$Behaviour_type <- factor(
  cv$Behaviour_type,
  levels = sort(unique(cv$Behaviour_type), decreasing = TRUE)
)

mod.behav_cv <- rma.mv(yi = yi, V = vi,
                       random = list(~1 | Study_ID / Obs_ID,
                                     ~1 | Species2), 
                       data =  cv,
                       mods = ~ Behaviour_type -1,
                       test = "t",
                       sparse = TRUE)

summary(mod.behav_cv)

behav.cv <- orchard_plot(mod.behav_cv, mod = "Behaviour_type", 
                      xlab = "Heterogeneity in risk-taking behaviour (lnCVR)", 
                      group = "Study_ID",
                      angle = 0) 

behav.cv <- behav.cv +
  scale_fill_manual(values = rep("#989aae", 7)) +
  scale_color_manual(values = rep("#989aae", 7))
behav.cv

###CV_ Sex######

ses_cv <- filter(cv, Sex != "N/A")


mod.sex.cv <- rma.mv(yi = yi, V = vi,
                   random = list(~1 | Study_ID / Obs_ID,
                                 ~1 | Species2), 
                     data =  ses_cv,
                     mods = ~ Sex -1,
                     test = "t",
                     sparse = TRUE)

summary(mod.sex.cv) 


###CV_real vs sim#######

cv$Real_predator <- as.factor(cv$Real_predator)

str(cv)


mod.rp.cv <- rma.mv(yi = yi, V = vi,
                    random = list(~1 | Study_ID / Obs_ID,
                                  ~1 | Species2), 
                    data =  cv,
                    mods = ~ Real_predator -1,
                    test = "t",
                    sparse = TRUE,
                    R = list(Species = cor1))

summary(mod.rp.cv)

#only real predators create heterogeneity in risk-taking behaviour.

###CV low vs no predator#####

mod.ln.cv <- rma.mv(yi = yi, V = vi,
                    random = list(~1 | Study_ID / Obs_ID,
                                  ~1 | Species2), 
                    data =  cv,
                    mods = ~ Low_or_no_pred -1,
                    test = "t",
                    sparse = TRUE)

summary(mod.ln.cv)

###CV adult vs juvenille######

cv$Adult <- factor(cv$Adult, 
                   levels = c("Adult", "Juvenille", "Both"))

mod.aj.cv <- rma.mv(yi = yi, V = vi,
                    random = list(~1 | Study_ID / Obs_ID,
                                  ~1 | Species2), 
                    data =  cv,
                    mods = ~ Adult -1,
                    test = "t",
                    sparse = TRUE)

summary(mod.aj.cv)

aj.cv <- orchard_plot(mod.aj.cv, xlab = "Difference in risk-taking (Hedge's g)", group = "Study_ID", mod = "Adult",
                      angle = 0) +
  scale_fill_manual(values = cvcolour) +
  scale_colour_manual(values = cvcolour)

aj.cv


###CV_comp type#######

mod.com.cv <- rma.mv(yi = yi, V = vi,
                     random = list(~1 | Study_ID / Obs_ID,
                                   ~1 | Species2), 
                     data =  cv,
                     mods = ~ Comparison_type -1,
                     test = "t",
                     sparse = TRUE)

summary(mod.com.cv)

###CV_class######

cv <- cv %>%
  mutate(Class = recode(Class,
                        "Actinopterygii" = "Ray-finned fishes",
                        "Amphibia" = "Amphibians",
                        "Mammalia" = "Mammals",
                        "Reptilia" = "Reptiles",
                        "Malacostraca" = "Crustaceans",
                        "Chondrichthyes" = "Cartilaginous fishes",
                        "Aves" = "Birds"
  ))



cv$Class <- factor(
  cv$Class,
  levels = sort(unique(cv$Class), decreasing = TRUE)
)

mod.class.cv <- rma.mv(yi = yi, V = vi,
                       random = list(~1 | Study_ID / Obs_ID,
                                     ~1 | Species2), 
                       data =  cv,
                       mods = ~ Class -1,
                       test = "t",
                       sparse = TRUE)

summary(mod.class.cv)

class.cv <- orchard_plot(mod.class.cv, mod = "Class", xlab = "Heterogeneity in risk-taking (lnCVR)", group = "Study_ID",
                      angle = 0) + 
  scale_fill_manual(values = rep("#989aae", 7)) +
  scale_color_manual(values = rep("#989aae", 7))
class.cv  

###predator pressure.cv######
mod.pp.cv <- rma.mv(yi = yi, V = vi,
                    random = list(~1 | Study_ID / Obs_ID,
                                  ~1 | Species2), 
                    data =  cv,
                    mods = ~ Predation_pressure -1,
                    test = "t",
                    sparse = TRUE,
                    R = list(Species = cor1))

summary(mod.pp.cv)



###figures #######

#put overall and overall_cv together using cowplot



Fig1 <- plot_grid(overall, overall_cv, labels = c("A", "B"), label_size = 12)

Fig1

Fig4 <- plot_grid(s, s.cv, aj,  aj.cv, labels = c("A", "B", "C", "D"), label_size = 12)

Fig4


###Tables#######

# >>> Model lists --------------------------------------------------------
#Hedges table 
mod1 <- coef(summary(mod.overall))

mod2 <- coef(summary(mod.behav)) 

mod3 <- coef(summary(mod.sex))

mod4 <- coef(summary(mod.class))

mod5 <- coef(summary(mod.pp))

mod6 <- coef(summary(mod.ln))

mod7 <- coef(summary(mod.aj))

mod8 <- coef(summary(mod.rp))

mod9 <- coef(summary(mod.com))

mod10 <- coef(summary(mod.behav.sens))

table1 <- rbind(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10)


table1$term <- rownames(table1)

#CV table


mod1.cv <- coef(summary(mod.overall_cv))

mod2.cv <- coef(summary(mod.behav_cv)) 

mod3.cv <- coef(summary(mod.sex.cv))

mod4.cv <- coef(summary(mod.class.cv))

mod5.cv <- coef(summary(mod.rp.cv))

mod6.cv <- coef(summary(mod.ln.cv))

mod7.cv <- coef(summary(mod.aj.cv))

mod8.cv <- coef(summary(mod.com.cv))

mod9.cv <- coef(summary(mod.pp.cv))

table2 <- rbind(mod1.cv, mod2.cv, mod3.cv, mod4.cv, mod5.cv, mod6.cv, mod8.cv, mod9.cv)

table2$term <- rownames(table2)

write.csv(table1, "Sup_table1.csv")

write.csv(table2, "Sup_table2.csv")

###END##### 



