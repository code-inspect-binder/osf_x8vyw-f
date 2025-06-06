# ------------------------------------------------------------------------
# file: 01_data_setup.R author: Jason Grafmiller date:
# 2017/07/24

# Description: Load in dataframe for analysis as reported in
# Grafmiller, J & B. Szmrecsanyi. 2017. 'Mapping out particle
# placement in Englishes around the world'. Language
# Variation and Change.

# ------------------------------------------------------------------------

# Source the setup file
source("R/00_project_setup.R")

# Read in the dataframe from .rds file (more efficient)
pv <- readRDS("data/particle_verbs_17-04-2018.rds")
# Or read from tab-delim text file
# pv <- read.delim("particle_verbs_df.txt")
# You can create a .txt file from the .rds as follows:
# write.table(pv, "particle_verbs_df.txt", sep = "\t", quote = F, row.names = F)

# Set variety labels.
vars <- c("GB", "CA", "NZ", "IE", "JA", "SG", "HK", "PH", "IN")

# Set register labels
regs <- c("spok.informal", "spok.formal", "writ.informal", "writ.formal",
  "online")

# Set genre category labels and reorder `Genre`.
cats <- c("PrivateDia", "UnscriptMono", "PublicDia", "ScriptedMono",
          "Letters", "CreativeWrit", "Blog", "GeneralWeb", "StudentWrit",
          "Reportage", "PopularWrit", "InstructWrit", "PersuasiveWrit",
          "AcademicWrit")

# Remove some extreme data points and NAs.
pv <- subset(pv, DirObjThematicity < 37 & !is.na(Rhythm)) %>%
  droplevels

# Convert some columns to character (if loading from .txt), and set the levels for other factors
pv <- pv %>%
  dplyr::mutate(
    Genre = factor(Genre, levels = cats),
    Register = factor(Register, levels = regs),
    PrimeType = relevel(PrimeType, ref = "none"),
    Variety = factor(Variety, levels = vars))
