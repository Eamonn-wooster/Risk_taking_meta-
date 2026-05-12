#title: Phylo Tree for multi-level meta-analysis

#Author: EIF Wooster

library("pacman")
#install.packages("pacman")

#install.packages("remotes")
#remotes::install_github("ropensci/rotl", force = TRUE)

#install.packages("rotl")# You'll have to install the package below
pacman::p_load(devtools, 
               tidyverse, 
               R.rsp, 
               dplyr,
               broom,
               tidyverse,
               here,
               broom,
               data.table,
               rotl,
               ape)

#Setting Directory

here()

#Reading in data

dat <- read.csv(here("Data/Analysis_ready.csv"))

###############################################################################
######### Species tree ##########
###############################################################################

#Pull out all species - leaving object as Species 
Species <- dplyr::select(dat, Species) #after column is you latin species name column

Species <- unique(Species)

length(Species$Species)

setDT(Species)

str(Species)

#Few names changes 

Species$Species <- gsub("Equus_quagga","Equus_burchellii_quagga", Species$Species)
Species$Species <- gsub("Macropus_eugenii","Notamacropus_eugenii",  Species$Species)
Species$Species <- gsub("Hyla_versicolor", "Dryophytes_versicolor",  Species$Species)
Species$Species <- gsub("Bufo_americanus", "Anaxyrus_americanus",  Species$Species)
Species$Species <- gsub( "Gambusiaa_hubbsi","Gambusia_hubbsi",  Species$Species)


# Now pass the character vector
resolved_spp <- tnrs_match_names(Species$Species, context_name = "Animals")

# Look at the matches
resolved_spp

# Extract ott_ids
ott_ids <- resolved_spp$ott_id

# Remove NA ott_ids
clean_ott_ids <- ott_ids[!is.na(ott_ids)]

# Now check which are in the tree
is_in_tree(clean_ott_ids)

#Checking for failed matches - NONE
failed_matches <- resolved_spp[is.na(resolved_spp$ott_id), ]
print(failed_matches$search_string)


Species_tree <- tol_induced_subtree(ott_ids = resolved_spp$ott_id)

length(Species_tree$tip.label)

# getting rid of the ott id   
Species_tree$tip.label <- gsub("_ott\\d+", "", Species_tree$tip.label)

#' [ EW This tells us the match between the two columns specified - what it's telling us is there is no match, hence why i've been adding]
setdiff(Species_tree$tip.label, Species$Species) 

# Writing new csvs with updated Study Species and Species_tree

out <- as.data.frame(Species_tree$tip.label)

write.csv(out, here("Outputs/Species_tree.csv"))

write.csv(dat, here("Outputs/Dat_for_analysis-ss-tree.csv")) 

## Species_tree appended to dat manually ##

# writing tree

tree1_Rev4 <- Species_tree

plot(tree1_Rev4, no.margin = TRUE)

write.tree(tree1_Rev4, file = here("Outputs/Tree_fire-Rev4.tre"))

