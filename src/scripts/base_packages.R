#######################################################################################################
# Load 'essential packages'
#######################################################################################################

# List of required packages
required_packages <- c("magrittr", "tidyverse", "pomp", "here", "ggplot2", "dplyr",  "scales")

# Check and install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg)
  } else {
    message(paste("Package already installed:", pkg))
  }
}

# Load all packages
lapply(required_packages, library, character.only = TRUE)

installed_packages <- installed.packages()[, "Package"]
save(installed_packages, file = "installed_packages.RData")


