# Packages used in the OMARI "Beyond the AUC" companion documents.
# Run once: source("install_packages.R")
pkgs <- c("dplyr", "tidyr", "forcats", "ggplot2", "glmnet", "glmnetUtils",
          "ranger", "rmarkdown", "knitr", "quarto")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing)
message("All packages installed. Also install Quarto from https://quarto.org if you want to render the .qmd files.")
