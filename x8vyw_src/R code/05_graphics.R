#  ------------------------------------------------------------------------
# file: 05_figures_and_tables.R
# author: Jason Grafmiller
# date: 2017/07/24

# Description: R code for creating all the figures and tables reported in
# Grafmiller, J & B. Szmrecsanyi. 'Mapping out particle placement in
# Englishes around the world'. Language Variation and Change.
#  ------------------------------------------------------------------------

# Source the setup file and read in the dataframe.
source("00_project_setup.R")
source("01_load_data.R")

library(ggsci)
library(RColorBrewer)
library(scatterplot3d)

# For conditional inference trees we use the partykit package.
tree_ctrl <- partykit::ctree_control(mincriterion = .99, maxdepth = 4)

# set global theme and color scheme
theme_set(theme_minimal())
theme_update(legend.position = 'top',
             axis.text = element_text(color = "black"))

cols <- brewer.pal(3, "BrBG")[-2]
scale_fill_discrete <- function(...) {
  scale_fill_manual(name = "Variant", values = cols)}


# Figures =================================================================

# Figure 1 ----------------------------------------------------------------

# Distribution of PV variants by Variety and Corpus

ycols <- c("white", "black")

# for saving
tiff("figures/figure1.tiff", res = 300, height = 5, width = 8, units = 'in')
pv %>%
  dplyr::count(Variety, Corpus, Response) %>%
  mutate(Variety = factor(Variety, levels = vars),
    Response = fct_recode(Response, Split = "Discontinuous"),
    Corpus = factor(Corpus, levels = c("ICE", "GloWbE"))) %>%
  group_by(Variety, Corpus) %>%
  mutate(
    Prop = n/sum(n),
    pos = c(.95, .05)) %>%
  ggplot(aes(Variety, Prop)) +
  geom_bar(aes(fill = Response), stat = "identity",
    width = .7, col = "black") +
  geom_text(aes(label = n, y = pos),
    col = rep(ycols, 18), size = 3) +
  facet_wrap(~ Corpus, ncol = 2) +
  labs(x = "", y = "Proportion of tokens") +
  scale_fill_jco(name = "") +
  theme_minimal() + theme(legend.position = 'top',
      axis.text = element_text(color = "black", size = rel(1.1)),
      axis.title = element_text(color = "black", size = rel(1.1)),
      legend.text = element_text(color = "black", size = rel(1)),
      strip.text = element_text(size = rel(1.1)))
dev.off()

# Figure 2 ----------------------------------------------------------------

# Heatmap of distribution of PV variants by Variety and Genre
tiff("figures/figure2.tiff", res = 300, height = 7, width = 7, units = 'in')
pv %>%
  mutate(Genre = fct_recode(Genre,
    Academic = "AcademicWrit", `Press editorials` = "PersuasiveWrit",
    `Instructional\nwriting` = "InstructWrit", `Private\nconversation` = "PrivateDia",
    `Public\nconversation` = "PublicDia", `Scripted\nmonologues` = "ScriptedMono",
    `Unscripted\nmonologues` = "UnscriptMono", `Creative\nwriting` = "CreativeWrit",
    `Popular\nnonfiction` = "PopularWrit", `Press reporting` = "Reportage",
    `General Web` = "GeneralWeb", `Online blogs` = "Blog")
    ) %>%
  group_by(Variety, Genre, Response) %>%
  count %>%
  group_by(Variety, Genre) %>%
  mutate(prop = n/sum(n),
    total = sum(n)) %>%
  subset(Response == "Discontinuous") %>%
  ggplot(aes(Variety, fct_rev(Genre))) +
    geom_tile(aes(fill = prop)) +
    geom_text(aes(label = total)) +
    scale_fill_gradient(name = "Proportion split variant",
      low = "white", high = "green4") +
    labs(x = "", y = "") +
    theme(legend.position = "bottom") +
    scale_x_discrete(position = "top")
dev.off()

# Figure 3 ----------------------------------------------------------------

# Plot of by-variety predictor importance ranking

# Load in varimp list if not already loaded
varimp_glmm_list <- readRDS('varimp_glmm_list.rds')
# remove "c." and "z." prefixes
preds <- rownames(varimp_glmm_list[[1]]) %>%
  gsub("^[cz]\\.", "", .)
rownames(varimp_glmm_list[[1]]) <- preds

df <- plyr::ldply(varimp_glmm_list, .id = 'Variety') %>%
  mutate(Predictor = rep(preds, 9) %>%
    factor(levels = rownames(varimp_glmm_list[[1]][order(varimp_glmm_list[[1]]$AICc), ]))
  )

tiff("figures/figure3.tiff", res = 300, height = 7, width = 7, units = 'in')
ggplot(df, aes(Predictor, AICc)) +
  geom_hline(yintercept = 0) +
  geom_segment(aes(xend = Predictor, yend = 0), col = "blue") +
  geom_point(col = "blue") +
  coord_flip() +
  facet_wrap(~Variety, scales = "free_x") +
  labs(x = "", y = "Increase in AICc")
dev.off()


# Figure 4 ----------------------------------------------------------------

# For figures 4 and 5 run '04_prob_distance_analysis.R'
library(phangorn)

# Neighbor net diagram of constraint effects
tiff("figures/figure4.tiff", res = 1200, height = 5, width = 5, units = 'in')
par("mar" = rep(1, 4))
set.seed(100)
plot(coef_nnet, "2D")
dev.off()

# Figure 5 ----------------------------------------------------------------

# Neighbor net diagram of constraint rankings

tiff("figures/figure5.tiff", res = 1200, height = 5, width = 5, units = 'in')
par("mar" = rep(1, 4))
set.seed(20)
plot(AIC_nnet, "2D")
dev.off()
