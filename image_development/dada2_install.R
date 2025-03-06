#!/usr/bin/env Rscript
install.packages("devtools", repos = "https://cloud.r-project.org/")
install.packages("dplyr", repos = "https://cloud.r-project.org/" )
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.20")
devtools::install_local('/master.zip')