#  ------------------------------------------------------------------------
# file: 03_glmm_analysis.R
# author: Jason Grafmiller
# date: 2017/07/24

# Description: R code for conducting by-variety generalized linear mixed model as reported in
# Grafmiller, J & B. Szmrecsanyi. 'Mapping out particle placement in
# Englishes around the world'. Language Variation and Change.
#  ------------------------------------------------------------------------

# Source the setup file and read in the dataframe.
source("R/00_project_setup.R")

# Load libraries for modelling.
library(lme4)
library(rms)
library(car)
library(merTools)
library(broom)

# Use sum coding for our categorical predictors, so set the global options
# for the regression models
# options(contrasts = c("contr.sum", "contr.poly"))

source("R/00_project_setup.R")

# Data --------------------------------------------------------------------

# Source data file.
source("R/01_load_data.R")
# glmm_list <- readRDS("glmm_list.rds")

# Model formulas ----------------------------------------------------------

# Main model formula
f0 <- Response ~ Register +
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

# Model formula for standardized predictors
f1 <- Response ~ Register +
  z.DirObjWordLength +
  c.Semantics +
  c.DirObjConcreteness +
  c.DirObjGivenness +
  c.DirObjDefiniteness +
  z.DirObjThematicity +
  c.DirectionalPP +
  CV.binary +
  z.Surprisal.P +
  z.Surprisal.V +
  z.Rhythm +
  PrimeType

# Add random effects
f1a <- update(f1, .~. + (1|Verb) + (1|Particle) + (1|VerbPart) + (1|Genre))


# Model controls ----------------------------------------------------------

# Set controls for glmer model
mer_ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e7))


# Fit model to each dataset -----------------------------------------------

# We split the dataframe, standardize numeric predictors, center binary predictors,
# and fit the model. See Gelman & Hill (2007:57) and Gelman (2008) on standardization
glmm_list <- vector("list")
t <- proc.time()
glmm_list <- pv %>%
  # (Nerd note: I find this list method with plyr more transparent and more versatile than
  # 'modern' methods using dplyr)
  plyr::dlply("Variety", .fun = function(x) {
    # Standardize predictors
    x <- mutate(x, z.Surprisal.P = z.(Surprisal.P, factor = 2),
      z.Surprisal.V = z.(Surprisal.V, factor = 2),
      z.DirObjThematicity = z.(DirObjThematicity, factor = 2),
      z.DirObjWordLength = z.(DirObjWordLength, factor = 2),
      z.Rhythm = z.(Rhythm, factor = 2),
      c.Semantics = c.(Semantics),
      c.DirObjConcreteness = c.(DirObjConcreteness),
      c.DirObjGivenness = c.(DirObjGivenness),
      c.DirObjDefiniteness = c.(DirObjDefiniteness),
      c.DirectionalPP = c.(DirectionalPP))
    # use sum coding with online as the intercept
    contrasts(x$Register) <- contr.sum(5)
    # fit model
    glmer(f1a, data = x, family = binomial, control = mer_ctrl)},
      .progress = "time")
t2 <- proc.time() - t
notify()

# save list to disk
# saveRDS(glmm_list, 'glmm_list.rds')


# drop random terms -------------------------------------------------------
# Not done in our analysis

glmm_list2 <- vector("list")
for(i in 1:9){
  cur_m <- glmm_list[[i]]
  cur_ranefs <- tidy(glmm_list[[i]]) %>%
    subset(group != 'fixed')
  if (min(cur_ranefs$estimate) < .001) {
    zeros <- cur_ranefs %>%
      subset(estimate < 0.001)
    zero_terms <- zeros[,1] %>%
      gsub(".*\\.", "", ., perl = T)
    if(length(zero_terms) > 1){
      new_fmla <- update(f1a,
        paste0(". ~ . -(1|", paste(zero_terms, collapse = ") -(1|"), ")"))
    } else {
      new_fmla <- update(f1a, paste0(". ~ . -(1|", zero_terms, ")"))
    }
    new_m <- update(m, formula = new_fmla, data = m@frame)
    glmm_list2[[i]] <- new_m
  } else {
    glmm_list2[[i]] <- cur_m
  }
}
names(glmm_list2) <- names(glmm_list)

# save list to disk
# saveRDS(glmm_list2, 'glmm_list2.rds')

# Evaluate models ---------------------------------------------------------

# Check models stats
glmm_list %>% plyr::ldply(.fun = broom::glance)

# check models predictions
glmm_list %>% plyr::ldply(.fun = scores_mer)

# Check models for overdispersion
glmm_list %>% plyr::ldply(.fun = overdisp_fun)

# Check predictors' variance inflation factors
glmm_list %>% plyr::ldply(.fun = VIF_mer)


# Extract model statistics ------------------------------------------------

# Extract coefficients and std errors
model_coefs <- glmm_list %>%
  plyr::ldply(.fun = tidy) %>%
  subset(group == 'fixed') %>%
  droplevels

# Creates a dataframe with 162 (9 models * 18 terms) rows
head(model_coefs, 10)

# Extract ranef variances
model_ranef_vars <- glmm_list %>%
  plyr::ldply(.fun = tidy) %>%
  subset(group != 'fixed') %>%
  droplevels

# Creates a dataframe with 36 (9 models * 4 random terms) rows
head(model_ranef_vars)


# Visualize coefficients --------------------------------------------------

# Patterns across all varieties by predictor
# (positive = greater probability of split variant)
model_coefs %>%
  subset(term != "(Intercept)") %>%
  droplevels %>%
  ggplot(aes(x = Variety, estimate, fill = Variety)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = estimate-1.96*std.error,
    ymax = estimate+1.96*std.error), width = .2) +
  facet_wrap(~ term)

# Patterns across all predictors by variety
# (positive = greater probability of split variant)
model_coefs %>%
  subset(term != "(Intercept)") %>%
  droplevels %>%
  ggplot(aes(x = term, estimate)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = estimate-1.96*std.error,
    ymax = estimate+1.96*std.error), width = .2) +
  facet_wrap(~Variety, scales = "free") +
  coord_flip()


# Visualize random effects ------------------------------------------------

# Extract simulated random effects.
# Uses simulation estimates based on the arm::sim() function.
model_ranef_sims <- glmm_list %>%
  plyr::ldply(.fun = merTools::REsim)

# Create simulation plots for each random effect in each model.
ranef_plots <- vector("list")
for(i in 1:9){
  p <- plotREsim(model_ranef_sims[[i]])
  p <- p + ggtitle(names(glmm_list)[i])
  ranef_plots[[i]] <- p
}

# Check plots for GB
ranef_plots[[1]]

do.call("grid.arrange", ranef_plots)


# Wald scores for predictor importance ------------------------------------

# Get the Wald test Chi sq values with Anova
Aov <- function(x) {
  as.data.frame(car::Anova(x)) %>%
    mutate(Predictor = rownames(.))
}

aov_df <- glmm_list %>%
  plyr::ldply(.fun = Aov) %>%
  with(., tapply(Chisq, list(Variety, Predictor), sum))

# The higher the Chisq, the more important that predictor is
aov_df %>%
  melt %>%
  ggplot(aes(Var2, value)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + labs(y = "Wald Chisq", x = "") +
  facet_wrap(~ Var1, scales = "free")


# Permutation variable importance -----------------------------------------

# This is a considerably time consuming process. Two methods
#
# Method 1: "Normal" loop through the data ------------------------------
set.seed(43214)
varimp_glmm_list <- vector('list')
t <- proc.time()
# varimp_glmm_list <- llply(glmm_list, .fun = permute.varimp(., verbose = F, ranef = T))
for (i in 1:9){
  print(names(glmm_list)[i]) # track how it's going...
  varimp_glmm_list[[i]] <- permute_varimp_mer(glmm_list[[i]], verbose = F, ranef = F)
}
t2 <- proc.time() - t
names(varimp_glmm_list) <- names(glmm_list2)
notify(paste("Ran in", t2[3]/60, "minutes"))

# Method 2: Parallel processing -----------------------------------------

library(foreach)
library(doParallel)

n_cores <- detectCores() # use all available cores
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

t <- proc.time()
varimp_glmm_list <- foreach(mod = 1:9, .combine = rbind) %dopar% permute_varimp_mer(glmm_list[[mod]], verbose = F, ranef = F)
stopCluster(cluster)

varimp_glmm_list <- varimp_glmm_list %>%
  mutate(Variety = factor(rep(names(glmm_list), each = 13), levels = vars),
    predictor = gsub("\\d$", "", rownames(.))) %>%
  plyr::dlply(., plyr::.(Variety))
t2 <- proc.time() - t
notify(paste("Ran in", t2[3]/60, "minutes"))

# save to file so we don't have to rerun the analysis
# saveRDS(varimp_glmm_list, file = "varimp_glmm_list.rds")

