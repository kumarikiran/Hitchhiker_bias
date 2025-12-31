#######################################################################################################
# This script runs SEIR model simulations across a grid of parameters
# and stores the simulated trajectories in the "output_data" folder. 
#######################################################################################################

# --- Setup ---
rm(list = ls())                       # Clear the R workspace
library(here)
here::here()                          # Set project root path (using {here})
source(here("scripts", "mod_SEIR.R"))       # Load SEIR model definition

# --- Define parameter grids ---

# Generate 100 values for initial exposure (E10) 
# spaced logarithmically between 1e4 and 1e-6
data_points <- exp(seq(log(1e4), log(1e-6), length.out = 100))

# Generate 10 evenly spaced values of rho2 between 0.01 and 1
rho_points <- seq(0.01, 1, length.out = 10)

# --- Ensure output directory exists ---
if (!dir.exists("output_data")) {
  dir.create("output_data")  # Create "output_data" if missing
}

# --- Run simulations ---
# Loop over rho2 values (10 total)
for (k in 1:10) {
  params_case1["rho2"] <- rho_points[k]   # Set rho2
  
  # Loop over E10 values (100 total)
  for (i in 1:100) {
    params_case1["E10"] <- data_points[i]   # Set E10
    
    # Simulate trajectory with pomp
    tj_1 <- trajectory(object = SIR_matrix_mod, params = params_case1, format = "data.frame") %>%
      select(-c(".id")) %>%                  # Drop pomp trajectory ID
      select(time, everything()) %>%         # Ensure time is first column
      mutate(across(where(is.numeric), ~ round(., 2)))  # Round numeric cols
    
    # Save results
    combined_df <- tj_1
    variable_name <- i
    rho_name <- k
    
    # File name encodes rho2 value and E10 index
    output_file <- sprintf("output_data/tj_combined_rho%s_E1_%s.rds", rho_points[k], variable_name)
    saveRDS(combined_df, file = output_file)
  }
}

# Show any warnings generated during the simulations
warnings()
