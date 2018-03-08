
if(!require(pacman))install.packages("pacman")
pacman::p_load(rlang,
              data.table,
               ggplot2,
               knitr) # add other packages needed to this list

# Sourcing all R files in the modules subdirectory
sourcefiles <- dir("DA_Functions", pattern="[.]R$", recursive = TRUE, full.names = TRUE)
for(z in sourcefiles)source(z)


