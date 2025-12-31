#######################################################################################################
# Purpose: Run SEIR simulations across a grid of initial exposures (E10) and rho2 values,
#          generate negative-binomial “observed” hospitalizations with noise, and save results.
# Assumes: 
#   - `SEIR_matrix_mod` (the simulator/ode object) and `params_case1` are defined via sourced scripts.
#   - Output directory is local to the project root (managed with {here}).
#######################################################################################################

# --- Setup ---
rm(list = ls())                       # Clear the R workspace
library(here)
here::here()                          # Set project root path (using {here})
source(here("scripts", "mod_SEIR.R"))       # Load SEIR model definition

# --- Experiment grid --------------------------------------------------------------------------------------------------
# E10 initial exposure values to sweep over 
data_points <- c(10000, 0.0001, 0.000001)

# Explore a few rho2 values (contact/transmission modifier; adjust as needed)
rho_points <- c(0.01, 0.10, 0.25, 0.50)

# Fixed parameter toggles
params_case1["hitchhikers"] <- 1                 

# --- Output configuration ---------------------------------------------------------------------------------------------
out_dir <- here::here("output_data_noise_4percent")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Number of noisy replicate observations per simulated trajectory
n_rep <- 100

# Negative-binomial noise parameters (reporting model)
# rho_mean: expected reporting multiplier (1 = no underreporting on average)
# rho_k: dispersion 
parms <- list(
  rho_mean = 1,
  rho_k    = 0.04
)



for (k in 1:4) {
  params_case1["rho2"] <- rho_points[k]
  for (i in 1:3) {
    params_case1["E10"] <-data_points[i]
    tj_1 <- trajectory(object = SIR_matrix_mod, params = params_case1, format = "data.frame") %>%
      mutate(.id = as.integer(.id)) %>%
      select(time, everything()) %>%
      mutate(across(where(is.numeric), ~ round(., 2)))
    sim_length = dim(tj_1)[1]
    sim_list <- vector("list", n_rep)

    for (n_i in 1:n_rep) {
      sim <- tj_1 %>%
        group_by(.id) %>%
        mutate(
          X_H2o_obs = rnbinom(n = length(X_H2o), 
                              mu = parms[["rho_mean"]] * X_H2o, 
                              size = 1 / parms[["rho_k"]]),
          X_H1o_obs = rnbinom(n = length(X_H1o), 
                              mu = parms[["rho_mean"]] * X_H1o, 
                              size = 1 / parms[["rho_k"]])
        ) %>%
        ungroup()
      sim_list[[n_i]] <- sim
    }
    
    
    sim_obs <- do.call(rbind, sim_list)
    sims_obs <- sim_obs %>%
      mutate(rep = rep(1:n_rep, times = 1, each = sim_length))
    
    
    ### end of generation the observation with noise 
    combined_df <-sims_obs
    variable_name <- i
    output_file <- sprintf("output_data_noise_4percent/tj_combined_rho%s_E1_%s.rds", rho_points[k], variable_name)
    saveRDS(combined_df, file = output_file)
  
  }
}
warnings()




