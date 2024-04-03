# Load required libraries -------------------------------------------------

## MAKE SURE ALL REQUIREMETS ARE MET AND LOAD LIBRARIES

# Check if the packages is already installed
if (!requireNamespace("readr", quietly = TRUE)) {
  # If not installed, install it
  install.packages("readr")
}
# Check if the package is already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  # If not installed, install it
  install.packages("ggplot2")
}

# Check if the package is already installed
if (!requireNamespace("dplyr", quietly = TRUE)) {
  # If not installed, install it
  install.packages("dplyr")
}

# Check if the package is already installed
if (!requireNamespace("corrr", quietly = TRUE)) {
  # If not installed, install it
  install.packages("corrr")
}

# Check if the package is already installed
if (!requireNamespace("ggcorrplot", quietly = TRUE)) {
  # If not installed, install it
  install.packages("ggcorrplot")
}
# Check if the package is already installed
if (!requireNamespace("FactoMineR", quietly = TRUE)) {
  # If not installed, install it
  install.packages("FactoMineR")
}
# Check if the package is already installed
if (!requireNamespace("factoextra", quietly = TRUE)) {
  # If not installed, install it
  install.packages("factoextra")
}
# Check if the package is already installed
if (!requireNamespace("ggbiplot", quietly = TRUE)) {
  # If not installed, install it
  install.packages("ggbiplot")
}
