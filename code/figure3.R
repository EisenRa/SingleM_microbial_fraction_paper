################################################################################
## Script to reproduce figure 3. 
## Raphael Eisenhofer Sept 2023
################################################################################

## load libraries/functions
library(tidyverse)
library(patchwork)

# Load data
div1 <- read_delim("data/pipe/simulated/Novel_genomes_no_reference_read_fraction.tsv")
div2 <- read_delim("data/pipe/simulated/Novel_genomes_placement_ani_read_fraction.tsv")




