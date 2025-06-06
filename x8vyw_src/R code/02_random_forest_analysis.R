#  ------------------------------------------------------------------------
# file: 02_random_forest_analysis.R
# author: Jason Grafmiller
# date: 2017/07/24

# Description: R code for conducting random forest analysis as reported in
# Grafmiller, J & B. Szmrecsanyi. 2017. 'Mapping out particle placement in
# Englishes around the world'. Language Variation and Change.

#  ------------------------------------------------------------------------

# Source the setup file
source("R/00_project_setup.R")

# Load additional libraries
library(party)
library(edarf)


# data --------------------------------------------------------------------

# read in the dataframe.
source("R/01_load_data.R")


# Model formula -----------------------------------------------------------

f0_rf <- Response ~ Variety +
  Genre +
  DirObjWordLength +
  Semantics +
  DirObjConcreteness +
  DirObjGivenness +
  DirObjDefiniteness +
  DirObjThematicity +
  DirectionalPP +
  CV.binary +
  Surprisal.P +
  Surprisal.V +
  Rhythm +
  PrimeType


# Random forest analysis --------------------------------------------------

# Set random seed.
set.seed(43214)

# Set control hyperparameters:
# - 'ntree': number of trees
# - 'mtry': number of predictors per tree
for_ctrl <- party::cforest_unbiased(ntree = 1500, mtry = 4)

# Fit tree model process:

t <- proc.time() # track duration of the process
# fit the model
cforest1 <- party::cforest(f0_rf, data = pv, controls = for_ctrl)

# Get forest predictions.
preds <- party::treeresponse(cforest1)

# Get predicted probability of the discontinuous variant.
rf_probs <- sapply(preds, function(x) return(x[2])) %>% as.vector
t2 <- proc.time() - t
# Generate pop up message alert when process has finished
msg <- paste("Process ran in", round(t2[3]/60, 1), "minutes...")
notify(msg)

# Save forest to file.
# saveRDS(cforest1, 'cforest1.rds')


# Model statistics --------------------------------------------------------

# Calculate accuracy measures.
# C index and Somer's Dxy
Hmisc::somers2(rf_probs, as.numeric(pv$Response) - 1) %>%
  round(3)

# Percent correctly predicted.
length(which(pv$Response == ifelse(rf_probs > .5, 'Discontinuous', 'Continuous')))/nrow(pv)

# Baseline accuracy
baseline <- max(table(pv$Response))/nrow(pv)
baseline

# Binomial test for significant improvement in accuracy.
ifelse(rf_probs > .5, 0, 1) %>%
  sum %>% # count number of continuous responses
  binom.test(., nrow(pv), p = baseline)


# Variable importance -----------------------------------------------------

# Compute variable importance based on AUC (same as index of concordance C).
# If process takes extremely long time, set 'conditional = F'

system.time( # track duration (can take a loooong time)
  varimp1 <- party::varimpAUC(cforest1, conditional = F))
# saveRDS(varimp1, 'data/forest_varimp1.rds')
notify() # Generate pop up message alert when process has finished


# Partial dependencies ----------------------------------------------------

# Compute partial dependency dataframes of individual predictors
# edarf package vignette for further info:
# https://cran.r-project.org/web/packages/edarf/vignettes/edarf.html

# Get predictor labels from formula (minus Response and Variety)
predictors <- all.vars(f0_rf)[-c(1:2)]

# Loop through predictors and calculate the partial dependencies of each
# predictor in interaction with Variety.
pd_list <- vector("list")
for (i in 1:length(predictors)){
  p <- predictors[i]
  cat(paste(' working on...', p))
  if (is.factor(pv[, p]) & nlevels(pv[, p]) > 2){
    # if predictor is factor with more than 2 levels
    nlevs <- 9 * nlevels(pv[, p])
  } else if (is.numeric(pv[, p])){
    if (length(unique(pv[, p])) < 10){
      # if predictor is numeric but has fewer than 10 unique values
      nlevs <- 9 * length(unique(pv[, p]))
    } else nlevs <- 90
  } else nlevs <- 18 # Variety * binary predictor
  pd <- partial_dependence(cforest1, c('Variety', p),
    interaction = T, n = c(nlevs, 25))
  pd_list[[i]] <- pd
  rm(nlevs, pd)
}
names(pd_list) <- all.vars(f0_rf)[-c(1:2)] # add names to list

# Save list of partial dependencies to file.
# saveRDS(pd_list, file = 'pd_list.rds')

# For graphs, see code in file 'supplementary_materials.Rmd'
