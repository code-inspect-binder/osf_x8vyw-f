#  ------------------------------------------------------------------------
# file: 04_prob_distance_analysis.R
# author: Jason Grafmiller
# date: 2017/07/24

# Description: R code for conducting multidimensional scaling, neighborNet, and clustering
# analyses of by-variety generalized linear mixed models as reported in
# Grafmiller, J & B. Szmrecsanyi. 2018. 'Mapping out particle placement in
# Englishes around the world'. Language Variation and Change.
#  ------------------------------------------------------------------------

# Source the setup file and read in the dataframe.
source("00_project_setup.R")

# Load libraries for glmm modelling and NeighborNets
library(lme4)
library(phangorn)


# Read in models ----------------------------------------------------------

# This list contains the 9 per-variety glmm models
glmm_list <- readRDS('glmm_list.rds')
varimp_glmm_list <- readRDS('varimp_glmm_list.rds')


# Distance matrices -------------------------------------------------------

# create quick function for extracting fixed effects coefficients
coeffs <- function(x) summary(x)$coefficients[,1]

# Coefficients
coef_df <- glmm_list %>%
  plyr::ldply(.fun = coeffs, .id = "Variety")
rownames(coef_df) <- coef_df$Variety
coef_dist <- dist(coef_df[, -1])

# AICc varimps
AIC_df <- plyr::ldply(varimp_glmm_list, .fun = function(x) x[, 3], .id = "Variety")
names(AIC_df) <- c("Variety", rownames(varimp_glmm_list[[1]]))
rownames(AIC_df) <- AIC_df$Variety

AIC_dist <- 1 - cor(t(AIC_df[, -1]), method = "spearman") %>%
  as.dist

# Correlations
mcor <- AIC_df[, 3] %>%
  t %>%
  cor(method = "spearman") %>%
  round(2)


# Mantel test -------------------------------------------------------------

# test correlations between matrices

# install.packages("ade4")
ade4::mantel.rtest(coef_dist, AIC_dist, nrepet = 999)

# inspect observed value against the distribution of the simulations
ade4::mantel.rtest(coef_dist, AIC_dist, nrepet = 999) %>%
  plot


# Constraint effect neighborNet -------------------------------------------

library(phangorn)

coef_nnet <- neighborNet(coef_dist)
par("mar" = rep(1, 4))
set.seed(100)
plot(coef_nnet, "2D")


# Contraint ranking neighborNet -------------------------------------------

AIC_nnet <- neighborNet(AIC_dist)
par("mar" = rep(1, 4))
set.seed(20)
plot(AIC_nnet, "2D")


# Cluster analysis --------------------------------------------------------

# Constraints
hclust(coef_dist, method = "ward.D2") %>%
  plot(main = "Coefficient-based dendrogram")

# Ranking
hclust(AIC_dist, method = "ward.D2") %>%
  plot(main = "Ranking-based dendrogram")


# Constraint effect MDS ---------------------------------------------------

coef_mds <- cmdscale(coef_dist, eig = TRUE, k = 3)

# Make dataframe for plotting.
coef_mds_df <- as.data.frame(coef_mds[[1]])
names(coef_mds_df) <- c("x", "y", "z")
coef_mds_df$Variety <- rownames(coef_mds_df)

# Plot in 2 dimensions
ggplot(coef_mds_df, aes(x, y)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_label(aes(label = Variety), fill = 'blue',
             color = "white", fontface = 'bold',
             label.padding = unit(0.3, "lines"))

# Plot in 3 dimensions
library(plotly)

plot_ly(coef_mds_df,
        type = "scatter3d", mode = "markers",
        x = ~x, y = ~y, z = ~z, text = ~Variety) %>%
  add_text(textfont = list(size = 16, color = "black"), textposition = "top right") %>%
  layout(showlegend = F)

# Constraint ranking MDS --------------------------------------------------

# Run MDS
AIC_mds <- cmdscale(AIC_dist, eig = TRUE, k = 3, add = T)

# Make dataframe for plotting.
AIC_mds_df <- as.data.frame(AIC_mds[[1]])
names(AIC_mds_df) <- c("x", "y", "z")
AIC_mds_df <- dplyr::mutate(AIC_mds_df,
                            Variety = rownames(AIC_mds_df))

# plot in 2 dimensions
ggplot(AIC_mds_df, aes(x, y)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_label(aes(label = Variety), fill = "blue",
             color = "white", fontface = 'bold',
             label.padding = unit(0.3, "lines"))

# plot in 3 dimensions
plot_ly(AIC_mds_df,
        type = "scatter3d", mode = "markers",
        x = ~x, y = ~y, z = ~z, text = ~Variety) %>%
  add_text(textfont = list(size = 16, color = "black"), textposition = "top right") %>%
  layout(showlegend = F)

