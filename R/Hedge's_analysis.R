################### Effect Size & Meta analysis Modelling #####################
## Authors: EW & CM


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

dat <- read.csv(here("Outputs/Dat_for_analysis-ss-mass-biome-clim-fireact-tax-Rev4.csv")) 

# Reading in tree

tree1 <- read.tree(here("Outputs/Tree_fire-Rev4.tre")) 

# Getting branch length and correlation matrix

tree1b <- compute.brlen(tree1)
cor1 <-  vcv(tree1b, corr=T)

# checking the match 

setdiff(dat$Species_tree, tree1$tip.label)

#' [EW - This HAS to return character(0) for the model to run]

# creating non-phylo columns - has to be tree column 

dat$Species2 <- dat$Species_tree

#################################################################################
######################### Converting Error ######################################
#################################################################################

### functions for converting error ------

# Function for converting CI's to SE

ci_to_se <-function(upper, lower, n){
  # upper: upper limit of CI
  # lower: lower limit of CI
  # n: sample size
  # return: SE
  se_value<-(upper - lower) / (2 * qt(0.975, n-1))
  return(se_value)
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

unique(dat$CI_Type)

### >>> Standard error --------

SD <- filter(dat, CI_Type == "SD")

#Taking the mean from the upper gives us the SD values given that SD is a single number (same applies to SE below)

SD$sd_Unburned <- SD$Upper_Unburned - SD$Mean_Unburned

SD$sd_Fire <- SD$Upper_Fire - SD$Mean_Fire

### >>> 95% CI's --------

CI <- filter(dat, CI_Type == "CI")

#Squaring gets it from SE to SD

CI$sd_Unburned <- ci_to_se(CI$Upper_Unburned, CI$Lower_Unburned, CI$Unburned_n)^2

CI$sd_Fire <- ci_to_se(CI$Upper_Fire, CI$Lower_Fire, CI$Burned_n)^2

### >>> Standard deviation --------

SE <- filter(dat, CI_Type == "SE")

# Taking the upper value from the mean gives us the SD values

SE$sd_Unburned <- SE$Upper_Unburned - SE$Mean_Unburned

SE$sd_Fire <- SE$Upper_Fire - SE$Mean_Fire

# now convert SE to SD

SE$sd_Unburned <- se_to_sd(SE$sd_Unburned, SE$Unburned_n)

SE$sd_Fire <- se_to_sd(SE$sd_Fire, SE$Burned_n)

# Put them back together!!

data <- rbind(SD, CI, SE)

#################################################################################
############################ Cleaning Data ######################################
#################################################################################

## Tidying up movement category 

# looking at unique category counts to identify duplicates 

table(data$Movement_category)

# Standardising all to lower case

data$Movement_category <- tolower(trimws(data$Movement_category))

# Back to Title case

data$Movement_category <- str_to_title(data$Movement_category)

# Cleaning Categories

table(data$Movement_category) 

data <- data %>%
  #filter(Movement_category != "Tortuosity - Opposite Direction") %>% # transforming this one
  mutate(Movement_category = as.character(Movement_category),
         Movement_category = case_when(
           Movement_category == "Linear Distance - Speed" ~ "Linear Distance",
           Movement_category == "Other Space Use - Time" ~ "Other Space Use",
           Movement_category == "Other Space Use - Proportion" ~ "Other Space Use",
           Movement_category == "Tortuosity - Opposite Direction" ~ "Tortuosity",
           TRUE ~ Movement_category
         ),
         Movement_category = as.factor(Movement_category))
         
table(data$Movement_category)    

## Cleaning Fire Category

# looking at unique category counts to identify duplicates 

table(data$Impact_fire_type)

data$Impact_fire_type <- factor(data$Impact_fire_type)

data$class_ncbi = factor(data$class_ncbi)

table(data$class_ncbi)

data <- data %>%
  mutate(class_ncbi = if_else(class_ncbi == "Lepidosauria",
                              "Reptilia",
                              class_ncbi))

# Changing NAs from the torts to reptiles

data$class_ncbi[is.na(data$class_ncbi)] <- "Reptilia"

#################################################################################
############################ Overall Effect Size ################################
#################################################################################

# Calculate Hedges' g

es <- escalc(
  measure = "SMD",       # standardized mean difference
  m1i = Mean_Fire,
  sd1i = sd_Fire,
  n1i = Burned_n,
  m2i = Mean_Unburned,
  sd2i = sd_Unburned,
  n2i = Unburned_n,
  data = data,
  vtype = "UB"           # unbiased estimator = Hedges' g
)

# Inspect the results

print(es) 

# obs level random effect

es$Obs_ID <- factor(1:nrow(es))

# Study level random effect

es$StudyID <- NULL

es$Study_ID <- factor(as.numeric(as.factor(es$Firstauthor_Year)))

# Correcting direction of Space Use Effect Sizes for Lovich 2018, Gasaway 1985 (percentage home range),
# and Tortuousity Harris2020 

es %>%
  dplyr::filter(Firstauthor_Year %in% c("Lovich2018", "Gasaway1985", "Harris2020")) %>%
  dplyr::select(Firstauthor_Year, Comments, Movement_category, Response_measured, Obs_ID)

es <- es %>%
  mutate(yi = ifelse(Obs_ID %in% c(75, 47, 110), yi * -1, yi))

################################## MA Models ##################################

# Variance covariance matrix

vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             data = es) 

#Overall model - we need this model - important to overall result and for pub bias analysis

rerun <- FALSE

if(rerun){
  mod.overall <- rma.mv(yi = yi, V = vcv,
                             random = list(~1 | Study_ID,
                                           ~1 | Species_tree, # phylo effect 
                                           ~1 | Species2, # non-phylo effect 
                                           ~1 | Obs_ID), 
                             data =  es,
                             control = list(optimizer="BFGS"),
                             test = "t",
                             sparse = TRUE,
                             R = list(Species_tree = cor1)
  )
  saveRDS(mod.overall, here("Outputs", "mod.overall.Rds"))
  summary(mod.overall)
}else{
  mod.overall <- readRDS(here("Outputs", "mod.overall.Rds"))
}

I2 = round(i2_ml(mod.overall), 2)

overall <- orchard_plot(mod.overall, xlab = "Change in movement (Hedge's g)", group = "Study_ID",
                        angle = 0)
overall

ggsave(here("Outputs", "Figures", "Overall.pdf"), width = 6, height = 4)

#This is the hetrogeneity of the effect sizes - we need to report this in the manuscript 

round(i2_ml(mod.overall), 2)

###############################################################################
###################### Moderator: Movement Category ###########################
###############################################################################

if(rerun){
  mod.overall_move <- rma.mv(yi = yi, V = vcv,
                             random = list(~1 | Study_ID,
                                           ~1 | Species_tree, # phylo effect 
                                           ~1 | Species2, # non-phylo effect 
                                           ~1 | Obs_ID), 
                             data =  es,
                             mods = ~ Movement_category -1,
                             control = list(optimizer="BFGS"),
                             test = "t",
                             sparse = TRUE,
                             R = list(Species_tree = cor1)
  )
  saveRDS(mod.overall_move, here("Outputs", "mod.overall_move.Rds"))
  summary(mod.overall_move)
}else{
  mod.overall_move <- readRDS(here("Outputs", "mod.overall_move.Rds"))
}

# Plotting 

overall_move <- orchard_plot(mod.overall_move, xlab = "Change in movement (Hedge's g)", group = "Study_ID", mod = "Movement_category",
                             angle = 0)
overall_move

# Saving Plot

ggsave(here("Outputs", "Figures", "Overall_move.pdf"), width = 6, height = 4)

#This is the hetrogeneity of the effect sizes - we need to report this in the manuscript 

round(i2_ml(mod.overall_move), 2)

###############################################################################
###################### Moderator: Species #####################################
###############################################################################

if(rerun){
  mod.species <- rma.mv(yi = yi, V = vcv,
                             random = list(~1 | Study_ID,
                                           ~1 | Species_tree, # phylo effect 
                                           ~1 | Species2, # non-phylo effect 
                                           ~1 | Obs_ID), 
                             data =  es,
                             mods = ~ Species2 -1,
                             control = list(optimizer="BFGS"),
                             test = "t",
                             sparse = TRUE,
                             R = list(Species_tree = cor1)
  )
  saveRDS(mod.species, here("Outputs", "mod_species.Rds"))
  summary(mod.species)
}else{
  mod.species <- readRDS(here("Outputs", "mod_species.Rds"))
}

# Plotting 

# Palette 
n_species <- 40  

n_species

random_palette <- sample(viridis(n_species, option = "D"))

#Plot

Species <- suppressMessages(
  orchard_plot(
    mod.species,
    xlab = "Change in movement (Hedge's g)",
    group = "Study_ID",
    mod = "Species2",
    angle = 0
  ) +
    scale_fill_manual(values = random_palette) +
    scale_color_manual(values = random_palette) +
    theme(axis.text.y = element_text(hjust = 0)) +
    theme(legend.position = "bottom", legend.justification = "right")
)

Species

ggsave(here("Outputs", "Figures", "Species.pdf"), width = 14, height = 12)


###############################################################################
############################# Moderator: Fire Type ############################
###############################################################################

# Filtering out not provided, unclear or both fire types

fire.es <- es %>%
  filter(Impact_fire_type %in% c("Planned", "Wildfire"))

fire.es$Impact_fire_type <- factor(fire.es$Impact_fire_type)

# Variance covariance matrix

fire.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                  data = fire.es) 

# Model with Fire Type

if(rerun){
  mod.fire <- rma.mv(yi = yi, V = fire.vcv,
                     random = list(~1 | Study_ID,
                                   ~1 | Species_tree, # phylo effect 
                                   ~1 | Species2, # non-phylo effect 
                                   ~1 | Obs_ID), 
                     data =  fire.es,
                     mods = ~ Impact_fire_type,
                     control=list(optimizer="BFGS"),
                     test = "t",
                     sparse = TRUE,
                     R = list(Species_tree = cor1)
  )
  saveRDS(mod.fire, here("Outputs", "mod_fire.Rds"))
  summary(mod.fire)
}else{
  mod.fire <- readRDS(here("Outputs", "mod_fire.Rds"))
}

# Plotting 

plot.fire <- orchard_plot(mod.fire, xlab = "Change in movement (Hedge's g)", 
                          group = "Study_ID", mod = "Impact_fire_type", angle = 0)
plot.fire

# Saving Plot

ggsave(here("Outputs", "Figures", "Fire_type.pdf"), width = 6, height = 4)


###############################################################################
######################## Moderator: Class #####################################
###############################################################################

table(es$class_ncbi)

# Filtering out Amphibia due to small sample size - agreed! 

class.es <- es %>%
  filter(class_ncbi %in% c("Aves", "Mammalia", "Reptilia"))

class.es$class_ncbi <- factor(class.es$class_ncbi)

# VCV

class.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                   data = class.es) 

# Model
if(rerun){
  mod.class <- rma.mv(yi = yi, V = class.vcv,
                      random = list(~1 | Study_ID,
                                    ~1 | Species_tree, # phylo effect 
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), 
                      data =  class.es,
                      mods = ~ class_ncbi -1, #remove intercept because we want a value for each group, not comparison between each.
                      control = list(optimizer="BFGS"),
                      test = "t",
                      sparse = TRUE,
                      R = list(Species_tree = cor1)
  )
  saveRDS(mod.class, here("Outputs", "mod_class.Rds"))
  summary(mod.class)
}else{
  mod.class <- readRDS(here("Outputs", "mod_class.Rds"))
}

plot.class <- orchard_plot(mod.class, xlab = "Change in movement (Hedge's g)", 
                           group = "Study_ID", mod = "class_ncbi", angle = 0)
plot.class

ggsave(here("Outputs", "Figures", "class.pdf"), width = 6, height = 4)


###############################################################################
######################### Biome as Moderator ##################################
###############################################################################

table(es$biome_name)

# Filtering out biomes with <10 es

biome.es <- es %>%
  group_by(biome_name) %>%
  filter(n() >= 10) %>%
  ungroup()

table(biome.es$biome_name)

# Filtering out Montane grasslands and shrublands and Tropical and subtropical 
#coniferous forests due to small sample size

#biome.es <- es %>%
  #filter(biome_name %in% c("Boreal forests/taiga ", 
                           "Deserts and xeric shrublands",
                           "Flooded grasslands and savannas", 
                           "Mediterranean forests, woodlands, and scrub or sclerophyll forests", 
                           "Temperate broadleaf and mixed forests",
                           "Temperate coniferous forests",
                           "Temperate grasslands, savannas, and shrublands",
                           "Tropical and subtropical grasslands, savannas, and shrublands",
                           "Tropical and subtropical moist broadleaf forests"))

biome.es$biome_name <- factor(biome.es$biome_name)

# VCV

biome.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                   data = biome.es) 

# Model

if(rerun){
  mod.biome <- rma.mv(yi = yi, V = biome.vcv,
                      random = list(~1 | Study_ID,
                                    ~1 | Species_tree, # phylo effect 
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), 
                      data =  biome.es,
                      mods = ~ biome_name -1, #' [We remove the intercept here bc we want estimates for all values and not a comparison to the first one]
                      control=list(optimizer="BFGS"),
                      test = "t",
                      sparse = TRUE,
                      R = list(Species_tree = cor1) #' [EW - this is the same as above]
  )
  saveRDS(mod.biome, here("Outputs", "mod_biome.Rds"))
  summary(mod.biome)
}else{
  mod.biome <- readRDS(here("Outputs", "mod_biome.Rds"))
}

# Plotting

plot.biome <- orchard_plot(mod.biome, xlab = "Change in movement (Hedge's g)", group = "Study_ID",
                           mod = "biome_name", angle = 0)

plot.biome

#'[EW - check this plot labels and work correctly since i removed the intercept from the model!]

#Fixing labels

wrapped_labels <- function(x) str_wrap(x, width = 25)

plot.biome + 
  scale_x_discrete(labels = wrapped_labels) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 8)) +
  theme(axis.title.x = element_text(size = 10))

# Saving Plot

ggsave(here("Outputs", "Figures", "Biome.pdf"), width = 8, height = 6)

###############################################################################
######################### Veg_clim as Moderator ###############################
###############################################################################

table(es$veg_climate)

# Filtering out closed_Arid and open_Temperate due to small sample size (<10)

veg_clim.es <- es %>%
  filter(veg_climate %in% c("closed_Cold", 
                            "closed_Temperate",
                            "open_Arid", 
                            "open_Tropical" ))


veg_clim.es$veg_climate <- factor(veg_clim.es$veg_climate)

# VCV

veg_clim.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                      data = veg_clim.es) 

# Model

if(rerun){
  mod.veg_clim <- rma.mv(yi = yi, V = veg_clim.vcv,
                         random = list(~1 | Study_ID,
                                       ~1 | Species_tree, # phylo effect 
                                       ~1 | Species2, # non-phylo effect 
                                       ~1 | Obs_ID), 
                         data =  veg_clim.es,
                         mods = ~ veg_climate -1, #' [We remove the intercept here bc we want estimates for all values and not a comparison to the first one]
                         control=list(optimizer="BFGS"),
                         test = "t",
                         sparse = TRUE,
                         R = list(Species_tree = cor1) #' [EW - this is the same as above]
  )
  saveRDS(mod.veg_clim, here("Outputs", "mod_veg_clim.Rds"))
  summary(mod.veg_clim)
}else{
  mod.veg_clim <- readRDS(here("Outputs", "mod_veg_clim.Rds"))
}

# Plotting

plot.veg_clim <- orchard_plot(mod.veg_clim, xlab = "Change in movement (Hedge's g)", group = "Study_ID",
                              mod = "veg_climate", angle = 0)

plot.veg_clim


# Saving Plot

ggsave(here("Outputs", "Figures", "Veg_clim.pdf"), width = 8, height = 6)

###############################################################################
######################### Climate as Moderator ###############################
###############################################################################

table(es$Climate)

es$Climate <- factor(es$Climate)

# Model

if(rerun){
  mod.clim <- rma.mv(yi = yi, V = vcv,
                     random = list(~1 | Study_ID,
                                   ~1 | Species_tree, # phylo effect 
                                   ~1 | Species2, # non-phylo effect 
                                   ~1 | Obs_ID), 
                     data =  es,
                     mods = ~ Climate -1, #' [We remove the intercept here bc we want estimates for all values and not a comparison to the first one]
                     control=list(optimizer="BFGS"),
                     test = "t",
                     sparse = TRUE,
                     R = list(Species_tree = cor1) #' [EW - this is the same as above]
  )
  saveRDS(mod.clim, here("Outputs", "mod_clim.Rds"))
  summary(mod.clim)
}else{
  mod.clim <- readRDS(here("Outputs", "mod_veg_clim.Rds"))
}

# Plotting

plot.clim <- orchard_plot(mod.clim, xlab = "Change in movement (Hedge's g)", group = "Study_ID",
                          mod = "Climate", angle = 0)

plot.clim

# Saving Plot

ggsave(here("Outputs", "Figures", "clim.pdf"), width = 8, height = 6)

###############################################################################
######################## Body Mass model : Overall  ###########################
###############################################################################

#We need to log trans form BM before modelling

es$logmass <- log(es$BodyMass.g)

if(rerun){
  mod.bm <- rma.mv(yi = yi, V = vcv,
                   random = list(~1 | Study_ID,
                                 ~1 | Species_tree, # phylo effect 
                                 ~1 | Species2, # non-phylo effect 
                                 ~1 | Obs_ID), 
                   data =  es,
                   mods = ~ logmass,
                   control=list(optimizer="BFGS"),
                   test = "t",
                   sparse = TRUE,
                   R = list(Species_tree = cor1)
  )
  saveRDS(mod.bm, here("Outputs", "mod1_bm.Rds"))
  summary(mod.bm)
}else{
  mod.bm <- readRDS(here("Outputs", "mod1_bm.Rds"))
}

mod.bm

#Extracting and plotting body mass models 

round(r2_ml(mod.bm)*100, 2) #This is r2 of the line


## replace logpredmass with you body mass column and herb_func with you data frame

newdat <- data.frame(logmass = seq(min(es$logmass, na.rm = TRUE),
                                   max(es$logmass, na.rm = TRUE),
                                   length.out = 100)) ## CM: Why is the length out 100?

X <- model.matrix(~ logmass, data = newdat)[, -1, drop = FALSE]

preds <- predict(mod.bm, newmods = X)

plot_data <- cbind(newdat, preds)

setDT(es)

es[, wi := 1/sqrt(vi)]
es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
es

prey_bm_text <- tidy(mod.bm)
prey_bm_text

#' [EW - i fixed this - you should run a seperate model for each class you can model, mammals etc. We might see a pattern there]
p.bm.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = es, 
              aes(x = logmass, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data,
              aes(x = logmass,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data, 
              aes(x = logmass,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data, 
            aes(x = logmass,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Log body mass (g)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.bm.mean

ggsave(here("Outputs/Figures", "body mass.pdf"), height = 4, width = 8)


###############################################################################
######################## Body Mass model : Mammals  ###########################
###############################################################################

bm.mam.es <- es %>%
  filter(class_ncbi %in% c("Mammalia"))

mam.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                data = bm.mam.es) 

if(rerun){
  mod.bm.mam <- rma.mv(yi = yi, V = mam.vcv,
                   random = list(~1 | Study_ID,
                                 ~1 | Species_tree, # phylo effect 
                                 ~1 | Species2, # non-phylo effect 
                                 ~1 | Obs_ID), 
                   data =  bm.mam.es,
                   mods = ~ logmass,
                   control=list(optimizer="BFGS"),
                   test = "t",
                   sparse = TRUE,
                   R = list(Species_tree = cor1)
  )
  saveRDS(mod.bm.mam, here("Outputs", "mod_bm_mam.Rds"))
  summary(mod.bm.mam)
}else{
  mod.bm.mam <- readRDS(here("Outputs", "mod_bm_mam.Rds"))
}

mod.bm.mam

#Extracting and plotting body mass models 

round(r2_ml(mod.bm.mam)*100, 2) #This is r2 of the line

## replace logpredmass with you body mass column and herb_func with you data frame

newdat.mam <- data.frame(logmass = seq(min(bm.mam.es$logmass, na.rm = TRUE),
                                   max(bm.mam.es$logmass, na.rm = TRUE),
                                   length.out = 100)) ## CM: Why is the length out 100?

X <- model.matrix(~ logmass, data = newdat.mam)[, -1, drop = FALSE]

preds.mam <- predict(mod.bm.mam, newmods = X)

plot_data.mam <- cbind(newdat.mam, preds.mam)

setDT(bm.mam.es)

bm.mam.es[, wi := 1/sqrt(vi)]
bm.mam.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
bm.mam.es

prey_bm.mam_text <- tidy(mod.bm.mam)
prey_bm.mam_text

#' [EW - i fixed this - you should run a seperate model for each class you can model, mammals etc. We might see a pattern there]
p.bm.mam.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = bm.mam.es, 
              aes(x = logmass, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.mam,
              aes(x = logmass,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.mam, 
              aes(x = logmass,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.mam, 
            aes(x = logmass,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Mammal Log body mass (g)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.bm.mam.mean

ggsave(here("Outputs/Figures", "body mass mammal.pdf"), height = 4, width = 8)

## ------------------- CM: Should I look at removing the two large outliers?

#bm.mam2.es <- es %>%
  filter(class_ncbi %in% c("Mammalia"))

#summary(bm.mam2.es$yi)

#bm.mam2.es <- bm.mam2.es %>%
  filter(abs(yi - 5.402831) > 1e-6)

#summary(bm.mam2.es$yi)

#bm.mam2.es <- bm.mam2.es %>%
  filter(abs(yi - 3.702135) > 1e-6)

#mam2.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                 data = bm.mam2.es) 

#if(rerun){
  #mod.bm.mam2 <- rma.mv(yi = yi, V = mam2.vcv,
                       random = list(~1 | Study_ID,
                                     ~1 | Species_tree, # phylo effect 
                                     ~1 | Species2, # non-phylo effect 
                                     ~1 | Obs_ID), 
                       data =  bm.mam2.es,
                       mods = ~ logmass,
                       control=list(optimizer="BFGS"),
                       test = "t",
                       sparse = TRUE,
                       R = list(Species_tree = cor1)
  )
  #saveRDS(mod.bm.mam2, here("Outputs", "mod_bm_mam2.Rds"))
  #summary(mod.bm.mam2)
}else{
  #mod.bm.mam2 <- readRDS(here("Outputs", "mod_bm_mam2.Rds"))
}

#mod.bm.mam2

#Extracting and plotting body mass models 

#round(r2_ml(mod.bm.mam2)*100, 2) #This is r2 of the line

## replace logpredmass with you body mass column and herb_func with you data frame

#newdat.mam2 <- data.frame(logmass = seq(min(bm.mam2.es$logmass, na.rm = TRUE),
                                       max(bm.mam2.es$logmass, na.rm = TRUE),
                                       length.out = 100)) ## CM: Why is the length out 100?

#X <- model.matrix(~ logmass, data = newdat.mam2)[, -1, drop = FALSE]

#preds.mam2 <- predict(mod.bm.mam2, newmods = X)

#plot_data.mam2 <- cbind(newdat.mam2, preds.mam2)

#setDT(bm.mam2.es)

#bm.mam2.es[, wi := 1/sqrt(vi)]
#bm.mam2.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
#bm.mam2.es

#prey_bm.mam2_text <- tidy(mod.bm.mam2)
#prey_bm.mam2_text

#' [EW - i fixed this - you should run a seperate model for each class you can model, mammals etc. We might see a pattern there]
#p.bm.mam2.mean <- ggplot()+
  # now add the rest of the points:
  #geom_jitter(data = bm.mam2.es, 
              aes(x = logmass, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  #geom_hline(yintercept = 0, lty = "dashed")+
 # geom_ribbon(data = plot_data.mam2,
              aes(x = logmass,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  #geom_ribbon(data = plot_data.mam2, 
              aes(x = logmass,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  #geom_line(data = plot_data.mam2, 
            aes(x = logmass,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  #scale_size_identity()+
  theme_bw()+
  xlab("Mammal2 Log body mass (g)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
#p.bm.mam2.mean

#ggsave(here("Outputs/Figures", "body mass mammal2.pdf"), height = 4, width = 8)


###############################################################################
######################## Body Mass model : Birds  #############################
###############################################################################

bm.ave.es <- es %>%
  filter(class_ncbi %in% c("Aves"))

aves.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                 data = bm.ave.es) 

if(rerun){
  mod.bm.ave <- rma.mv(yi = yi, V = aves.vcv,
                       random = list(~1 | Study_ID,
                                     ~1 | Species_tree, # phylo effect 
                                     ~1 | Species2, # non-phylo effect 
                                     ~1 | Obs_ID), 
                       data =  bm.ave.es,
                       mods = ~ logmass,
                       control=list(optimizer="BFGS"),
                       test = "t",
                       sparse = TRUE,
                       R = list(Species_tree = cor1)
  )
  saveRDS(mod.bm.ave, here("Outputs", "mod_bm_ave.Rds"))
  summary(mod.bm.ave)
}else{
  mod.bm.ave <- readRDS(here("Outputs", "mod_bm_ave.Rds"))
}

mod.bm.ave

#Extracting and plotting body mass models 

round(r2_ml(mod.bm.ave)*100, 2) #This is r2 of the line

## replace logpredmass with you body mass column and herb_func with you data frame

newdat.ave <- data.frame(logmass = seq(min(bm.ave.es$logmass, na.rm = TRUE),
                                       max(bm.ave.es$logmass, na.rm = TRUE),
                                       length.out = 100)) ## CM: Why is the length out 100?

X <- model.matrix(~ logmass, data = newdat.ave)[, -1, drop = FALSE]

preds.ave <- predict(mod.bm.ave, newmods = X)

plot_data.ave <- cbind(newdat.ave, preds.ave)

setDT(bm.ave.es)

bm.ave.es[, wi := 1/sqrt(vi)]
bm.ave.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
bm.ave.es

prey_bm.ave_text <- tidy(mod.bm.ave)
prey_bm.ave_text

#' [EW - i fixed this - you should run a seperate model for each class you can model, mammals etc. We might see a pattern there]
p.bm.ave.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = bm.ave.es, 
              aes(x = logmass, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.ave,
              aes(x = logmass,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.ave, 
              aes(x = logmass,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.ave, 
            aes(x = logmass,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Aves Log body mass (g)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.bm.ave.mean

ggsave(here("Outputs/Figures", "body mass aves.pdf"), height = 4, width = 8)

###############################################################################
######################## Body Mass model : Reptiles ###########################
###############################################################################

bm.rep.es <- es %>%
  filter(class_ncbi %in% c("Reptilia"))

rep.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                  data = bm.rep.es) 

if(rerun){
  mod.bm.rep <- rma.mv(yi = yi, V = rep.vcv,
                       random = list(~1 | Study_ID,
                                     ~1 | Species_tree, # phylo effect 
                                     ~1 | Species2, # non-phylo effect 
                                     ~1 | Obs_ID), 
                       data =  bm.rep.es,
                       mods = ~ logmass,
                       control=list(optimizer="BFGS"),
                       test = "t",
                       sparse = TRUE,
                       R = list(Species_tree = cor1)
  )
  saveRDS(mod.bm.rep, here("Outputs", "mod_bm_rep.Rds"))
  summary(mod.bm.rep)
}else{
  mod.bm.rep <- readRDS(here("Outputs", "mod_bm_rep.Rds"))
}

mod.bm.rep

#Extracting and plotting body mass models 

round(r2_ml(mod.bm.rep)*100, 2) #This is r2 of the line

## replace logpredmass with you body mass column and herb_func with you data frame

newdat.rep <- data.frame(logmass = seq(min(bm.rep.es$logmass, na.rm = TRUE),
                                       max(bm.rep.es$logmass, na.rm = TRUE),
                                       length.out = 100)) 

X <- model.matrix(~ logmass, data = newdat.rep)[, -1, drop = FALSE]

preds.rep <- predict(mod.bm.rep, newmods = X)

plot_data.rep <- cbind(newdat.rep, preds.rep)

setDT(bm.rep.es)

bm.rep.es[, wi := 1/sqrt(vi)]
bm.rep.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
bm.rep.es

prey_bm.rep_text <- tidy(mod.bm.rep)
prey_bm.rep_text 

#' [EW - i fixed this - you should run a seperate model for each class you can model, mammals etc. We might see a pattern there]
p.bm.rep.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = bm.rep.es, 
              aes(x = logmass, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.rep,
              aes(x = logmass,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.rep, 
              aes(x = logmass,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.rep, 
            aes(x = logmass,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Reptilia Log body mass (g)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.bm.rep.mean

ggsave(here("Outputs/Figures", "body mass rept.pdf"), height = 4, width = 8)

###############################################################################
############################ Burn Treatment Age Model #########################
###############################################################################

# Coercing variable to integer

es <- es %>%
  mutate(Maximum_treatment_burn_age_days = if_else(
    grepl("^[-]?[0-9]+$", Maximum_treatment_burn_age_days),
    as.integer(Maximum_treatment_burn_age_days),
    NA_integer_
  ))

# Looking for NAs

summary(es$Maximum_treatment_burn_age_days)

# Removing NAs

age.es <- es[!is.na(es$Maximum_treatment_burn_age_days), ]

es <- es %>%
  mutate(Maximum_treatment_burn_age_days = if_else(
    Maximum_treatment_burn_age_days == 60989, 
    609, 
    Maximum_treatment_burn_age_days
  ))

# Removing outlier and transforming data [ignore, this outlier was converted wrong]

# age.es <- age.es %>%
  #filter(Maximum_treatment_burn_age_days != "60989")

summary(age.es$Maximum_treatment_burn_age_days)

hist(age.es$Maximum_treatment_burn_age_days)

age.es$log.age <- log(age.es$Maximum_treatment_burn_age_days)

hist(age.es$log.age)

#' [EW - i'd try rerunning this after dropping the crazy outlier and then consider log transforming if you need to make it a but more normal]
 
#VCV 

age.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                   data = age.es) 
if(rerun){
  mod.age <- rma.mv(yi = yi, V = age.vcv,
                    random = list(~1 | Study_ID,
                                  ~1 | Species_tree, # phylo effect 
                                  ~1 | Species2, # non-phylo effect 
                                  ~1 | Obs_ID), 
                    data =  age.es,
                    mods = ~ log.age, #uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                    test = "t",
                    sparse = TRUE,
                    R = list(Species_tree = cor1)
  )
  saveRDS(mod.age, here("Outputs", "mod1_age.Rds"))
  summary(mod.age)
}else{
  mod.age <- readRDS(here("Outputs", "mod1_age.Rds"))
}

mod.age

#Extracting and plotting

round(r2_ml(mod.age)*100, 2) #This is r2 of the line

## ### replace logpredmass with you body mass column and herb_func with you data frame

newdat.age <- data.frame(log.age = seq(min(age.es$log.age, na.rm = TRUE),
                                                               max(age.es$log.age, na.rm = TRUE),
                                                               length.out = 100)) ## CM: Why is the length out 100?

X.age <- model.matrix(~ log.age, data = newdat.age)[, -1, drop = FALSE]

preds.age <- predict(mod.age, newmods = X.age)

plot_data.age <- cbind(newdat.age, preds.age)

setDT(age.es)

age.es[, wi := 1/sqrt(vi)]
age.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
age.es

prey_age_text <- tidy(mod.age)
prey_age_text

#' [EW - this should work better now, effect size should be on the y not se - lemme know if it doesn't plot like the bm one above]

p.age.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = age.es, ## CM Changed this to age.es to due error "object 'yi' not found" - Matches the BM model
              aes(x = log.age, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.age,
              aes(x = log.age,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.age, 
              aes(x = log.age,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.age, 
            aes(x = log.age,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Log Maximum Burn Age (days)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.age.mean

ggsave(here("Outputs/Figures", "Burn Age.pdf"), height = 4, width = 8)

###############################################################################
########################## Model: Fire Size ###################################
###############################################################################

# Coercing variable to integer

es <- es %>%
  mutate(Impact_fire_size_ha = if_else(
    grepl("^[-]?[0-9]+$", Impact_fire_size_ha),
    as.integer(Impact_fire_size_ha),
    NA_integer_
  ))

# Looking for NAs

summary(es$Impact_fire_size_ha)

# Removing NAs

size.es <- es[!is.na(es$Impact_fire_size_ha), ]
summary(size.es$Impact_fire_size_ha)

# Transforming data

hist(size.es$Impact_fire_size_ha)

size.es$log.size <- log(size.es$Impact_fire_size_ha)

hist(size.es$log.size)

#VCV 

size.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                  data = size.es) 
if(rerun){
  mod.size <- rma.mv(yi = yi, V = size.vcv,
                     random = list(~1 | Study_ID,
                                   ~1 | Species_tree, # phylo effect 
                                   ~1 | Species2, # non-phylo effect 
                                   ~1 | Obs_ID), 
                     data =  size.es,
                     mods = ~ log.size, #uncomment this and add whatever you want to test i.e., movement_category, body mass etc
                     test = "t",
                     sparse = TRUE,
                     R = list(Species_tree = cor1)
  )
  saveRDS(mod.size, here("Outputs", "mod1_size.Rds"))
  summary(mod.size)
}else{
  mod.size <- readRDS(here("Outputs", "mod1_size.Rds"))
}

mod.size

round(r2_ml(mod.size)*100, 2) #This is r2 of the line

## ### replace logpredmass with you body mass column and herb_func with you data frame

newdat.size <- data.frame(log.size = seq(min(size.es$log.size, na.rm = TRUE),
                                         max(size.es$log.size, na.rm = TRUE),
                                         length.out = 100)) ## CM: Why is the length out 100?

X.size <- model.matrix(~ log.size, data = newdat.size)[, -1, drop = FALSE]

preds.size <- predict(mod.size, newmods = X.size)

plot_data.size <- cbind(newdat.size, preds.size)

setDT(size.es)

size.es[, wi := 1/sqrt(vi)]
size.es[, pt_size := 2 + 7 * (wi-min(wi, na.rm = T)) / (max(wi, na.rm = T) - min(wi, na.rm = T))]
size.es

prey_age_text <- tidy(mod.size)
prey_age_text

p.size.mean <- ggplot()+
  # now add the rest of the points:
  geom_jitter(data = size.es, ## CM Changed this to age.es to due error "object 'yi' not found" - Matches the BM model
              aes(x = log.size, size = pt_size, 
                  y = yi),
              position = position_jitter(width = 0.01),
              inherit.aes = F,  alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_ribbon(data = plot_data.size,
              aes(x = log.size,
                  ymin = pi.lb, ymax = pi.ub),
              fill = "transparent", color = "#dad7cd",
              alpha = .3)+
  geom_ribbon(data = plot_data.size, 
              aes(x = log.size,
                  ymin = ci.lb, ymax = ci.ub),
              alpha = .3)+
  geom_line(data = plot_data.size, 
            aes(x = log.size,
                y = pred),
            alpha = .8)+
  
  # coord_cartesian(ylim = c(-4, 4))+
  scale_size_identity()+
  theme_bw()+
  xlab("Log Burn Area (ha)")+
  ylab("Hedge's G")+
  theme(legend.position = "none",
        text = element_text(color = "black", family = "Helvetica"), axis.text = element_text(color = "black", family = "Helvetica"),
        panel.grid = element_blank(),
        panel.border = element_blank())
p.size.mean

ggsave(here("Outputs/Figures", "Burn Size.pdf"), height = 4, width = 8)

################################################################################
################### Pausas2017 Fire Activity ###################################
################################################################################

summary(es$fireactivi)

hist(es$fireactivi)

if(rerun){
  mod.fireact <- rma.mv(yi = yi, V = vcv,
                        random = list(~1 | Study_ID,
                                      ~1 | Species_tree, # phylo effect 
                                      ~1 | Species2, # non-phylo effect 
                                      ~1 | Obs_ID), 
                        data =  es,
                        mods = ~ fireactivi,
                        control=list(optimizer="BFGS"),
                        test = "t",
                        sparse = TRUE,
                        R = list(Species_tree = cor1)
  )
  saveRDS(mod.fireact, here("Outputs", "mod_fireact.Rds"))
  summary(mod.fireact)
}else{
  mod.fireact <- readRDS(here("Outputs", "mod_fireact.Rds"))
}

mod.fireact

#Extracting and plotting model

round(r2_ml(mod.fireact)*100, 2) #This is r2 of the line

## Not much going on, haven't plotted

################################################################################
########################### Publication bias analyses #########################
################################################################################

es$effectN <- (es$Unburned_n * es$Burned_n) / (es$Unburned_n + es$Burned_n)
es$sqeffectN <- sqrt(es$effectN)

if(rerun){
  mod_diff_egg <- rma.mv(yi = yi, V = vcv,
                         random = list(~1 | Study_ID,
                                       ~1 | Species_tree, # phylo effect 
                                       ~1 | Species2, # non-phylo effect 
                                       ~1 | Obs_ID), 
                         data =  es,
                         mods = ~ sqeffectN,
                         control=list(optimizer="BFGS"),
                         test = "t",
                         sparse = TRUE,
                         R = list(Species_tree = cor1)
  )
  summary(mod_diff_egg)
  saveRDS(mod_diff_egg, here("Outputs", "eggers.Rds"))
}else{
  mod_diff_egg <- readRDS(here("Outputs", "eggers.Rds"))
}

mod_diff_egg #No significant intercept

##Funnel

funnel(mod_diff_egg, 
       yaxis="seinv",
       xlab = "Standardized residuals",
       ylab = "Precision (inverse of SE)",
       
)

ggsave("Outputs/funnel_main.pdf", width = 5, height = 5)

### Decline effect for difference

str(es$Publication_Year)

mod_mad <- rma.mv(yi = yi, V = vcv,
                                 random = list(~1 | Study_ID,
                                                  ~1 | Species_tree, # phylo effect 
                                                  ~1 | Species2, # non-phylo effect 
                                                  ~1 | Obs_ID), 
                                    data = es,
                                    mods = ~ Publication_Year,
                                    control=list(optimizer="BFGS"),
                                    test = "t",
                                    sparse = TRUE,
                                    R = list(Species_tree = cor1)
  )
  
  summary(mod_mad)
  
decline <- bubble_plot(mod_mad,
                       mod = "Publication_Year",
                       group = "Study_ID",
                       xlab = "Publication year",
                       ylab = "Hedge's g",
                       g = TRUE,
                       weights = 1/sqrt(es$vi))

decline ## if this is giving an error, just restart R and it should work

ggsave("Outputs/decline.pdf", width = 8, height = 4)

#OK you have a a result where the effect size is declining over time - this is known as a decline effect. 

### >>>> leave-one-out analysis 

#####if rerunning load the model below - this takes hours to run####

es$Firstauthor_Year <- as.factor(es$Firstauthor_Year)

es$leave_out <- as.factor(es$Firstauthor_Year)

LeaveOneOut_effectsize <- list()
for (i in 1:length(levels(es$leave_out))) {
  temp_dat <- es %>%
    dplyr::filter(leave_out != levels(es$leave_out)[i])
  
  VCV_leaveout <- vcalc(vi = temp_dat$vi, cluster = temp_dat$Study_ID, obs = temp_dat$Obs_ID, rho = 0.5)
  
  LeaveOneOut_effectsize[[i]] <-  rma.mv(yi = yi, V = VCV_leaveout,
                                         random = list(~1 | Study_ID,
                                                       ~1 | Species_tree, # phylo effect 
                                                       ~1 | Species2, # non-phylo effect 
                                                       ~1 | Obs_ID), 
                                         control=list(optimizer="BFGS"),
                                         test = "t",
                                         method = "REML", 
                                         sparse = TRUE,
                                         R = list(Species_tree = cor1),
                                         data = temp_dat[temp_dat$leave_out != levels(temp_dat$leave_out)[i], ])
}


# writing function for extracting est, ci.lb, and ci.ub from all models
est.func <- function(model) {
  df <- data.frame(est = model$b, lower = model$ci.lb, upper = model$ci.ub)
  return(df)
}

# using dplyr to form data frame
MA_oneout <- lapply(LeaveOneOut_effectsize,function(x) est.func(x)) %>%
  bind_rows %>%
  mutate(left_out = levels(es$leave_out))

# telling ggplot to stop reordering factors
MA_oneout$left_out <- as.factor(MA_oneout$left_out)
MA_oneout$left_out <- factor(MA_oneout$left_out, levels = MA_oneout$left_out)

# saving the runs
saveRDS(MA_oneout, here("data", "MA_oneout.RDS"))

#Load model here if rerunning 

MA_oneout <- readRDS(here("data", "MA_oneout.RDS"))
summary(MA_oneout)

# plotting
leaveoneout <- ggplot(MA_oneout) + geom_hline(yintercept = 0, lty = 2, lwd = 1) +
  geom_hline(yintercept = mod.overall$ci.lb, lty = 3, lwd = 0.75, colour = "black") +
  geom_hline(yintercept = mod.overall$b, lty = 1, lwd = 0.75, colour = "black") + 
  geom_hline(yintercept = mod.overall$ci.ub,
             lty = 3, lwd = 0.75, colour = "black") + 
  geom_pointrange(aes(x = left_out, y = est,
                      ymin = lower, ymax = upper)) + 
  xlab("Study left out") + 
  ylab("Difference in logit overlap (logit)") +
  coord_flip() + 
  theme(panel.grid.minor = element_blank()) + theme_bw() + theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) + theme(axis.text.y = element_text(size = 6))

leaveoneout

ggsave("Outputs/leaveoneout.pdf", width = 8, height = 4)

#' [END ? I think you are done - make sure to report each of these models - i've pasted code below for you to make tables easily ]

#for making tables for the supplement

mod1 <- coef(summary(mod.bm.rep)) #paste model names into the brackets

mod2 <- coef(summary(mod.bm.ave))

mod3 <- coef(summary(mod.bm.mam))

mod4 <- coef(summary(mod.age))

mod5 = coef(summary(mod.bm))

mod6 <- coef(summary(mod.size))

mod7 <- coef(summary(mod.fire))

mod8 <- coef(summary(mod.biome))

mod9 <- coef(summary(mod.class))

mod10 <- coef(summary(mod.overall_move))

mod11 <- coef(summary(mod.overall))

mod12 = coef(summary(mod_mad))

mod13 = coef(summary(mod_diff_egg))

mod14 = coef(summary(mod.veg_clim))

mod15 = coef(summary(mod.clim))

mod16 = coef(summary(mod.fireact))


table1 <- rbind(mod1,mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15, mod16) #bind the models together 

table1$term <- rownames(table1)

rownames(table1) <- NULL

write.csv(table1, "Outputs/Sup_tableRev5.csv") # Write it out to CSV to then paste into word doc

table2 = rbind(mod14, mod15)

table2$term = rownames(table2)

rownames(table2) = NULL

write.csv(table2, "Outputs/Sup_table3.csv")

