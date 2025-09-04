# Load required libraries
library(BayesianTools)  # For creating prior distributions compatible with Bayesian MCMC/SMC samplers
library(tibble)         # For creating and manipulating tidy tabular data (tribble)
library(dplyr)          # For data manipulation (filter, mutate, arrange, etc.)

# --------------------------------------------------------------------------------------------------
# Function: set_prior
# Purpose:  Define prior distributions for a set of model parameters and return a BayesianTools
#           prior object, which includes:
#             - A log-density function (dens_fun) for evaluating priors
#             - A sampling function (samp_fun) for generating draws from priors
#             - Lower and upper parameter bounds
#
# Input:    pars_nm - character vector of parameter names to include in the prior
# Output:   A prior object (class "prior") usable in BayesianTools samplers
# --------------------------------------------------------------------------------------------------

set_prior <- function(pars_nm) {
  
  # --------------------------------------------------------------------------------------------
  # Step 1: Define hard parameter bounds (lower = inf, upper = sup) for all possible parameters
  # --------------------------------------------------------------------------------------------
  bounds <- tribble(
    ~par,      ~inf,  ~sup,
    "Ri1p",     0.1,   1,   # Reproduction parameter for virus 1 (exponential prior support)
    "rho1",     0,     1,   # Reporting probability for virus 1 (bounded [0,1])
    "Ri2p",     0.1,   1,   # Reproduction parameter for virus 2 (exponential prior support)
    "rho2",     0,     1,   # Reporting probability for virus 2 (bounded [0,1])
    "thetap",   0,     1,   # Parameter theta (beta prior support)
    "psi",      0,     1    # Misc parameter psi (uniform prior support)
  )
  
  # Keep only the parameters requested in pars_nm, preserve their order
  bounds <- bounds %>%
    filter(par %in% pars_nm) %>%                  # Subset to requested parameters
    mutate(par = factor(par, levels = pars_nm)) %>% # Maintain input order
    arrange(par)
  
  # --------------------------------------------------------------------------------------------
  # Step 2: Define log-density function for the priors
  #         Each parameter has its own prior distribution:
  #           - Ri1p, Ri2p ~ Exponential(rate=0.1)
  #           - rho1, rho2, thetap ~ Beta(2, 2)
  #           - psi ~ Uniform(0, 1)
  # --------------------------------------------------------------------------------------------
  dens_fun <- function(par) {
    names(par) <- pars_nm    # Assign names for clarity
    d <- numeric(length(pars_nm))  # Initialize vector of log densities
    
    # Add contributions from each parameter if present
    if ("Ri1p" %in% pars_nm) d["Ri1p"]   <- dexp(par["Ri1p"],  rate = 0.1, log = TRUE)
    if ("rho1" %in% pars_nm) d["rho1"]   <- dbeta(par["rho1"], shape1 = 2, shape2 = 2, log = TRUE)
    if ("Ri2p" %in% pars_nm) d["Ri2p"]   <- dexp(par["Ri2p"],  rate = 0.1, log = TRUE)
    if ("rho2" %in% pars_nm) d["rho2"]   <- dbeta(par["rho2"], shape1 = 2, shape2 = 2, log = TRUE)
    if ("thetap" %in% pars_nm) d["thetap"] <- dbeta(par["thetap"], shape1 = 2, shape2 = 2, log = TRUE)
    if ("psi" %in% pars_nm)   d["psi"]    <- dunif(par["psi"], min = 0, max = 1, log = TRUE)
    
    return(sum(d))  # Return joint log-prior (sum of independent components)
  }
  
  # --------------------------------------------------------------------------------------------
  # Step 3: Define sampling function for the priors
  #         Generates random samples from the same distributions as above
  # --------------------------------------------------------------------------------------------
  samp_fun <- function(n = 1) {
    samp_mat <- matrix(NA, nrow = n, ncol = length(pars_nm))
    colnames(samp_mat) <- pars_nm
    
    # Sample each parameter according to its prior distribution
    if ("Ri1p" %in% pars_nm)   samp_mat[, "Ri1p"]   <- rexp(n,  rate = 0.1)
    if ("rho1" %in% pars_nm)   samp_mat[, "rho1"]   <- rbeta(n, 2, 2)
    if ("Ri2p" %in% pars_nm)   samp_mat[, "Ri2p"]   <- rexp(n,  rate = 0.1)
    if ("rho2" %in% pars_nm)   samp_mat[, "rho2"]   <- rbeta(n, 2, 2)
    if ("thetap" %in% pars_nm) samp_mat[, "thetap"] <- rbeta(n, 2, 2)
    if ("psi" %in% pars_nm)    samp_mat[, "psi"]    <- runif(n, min = 0, max = 1)
    
    return(samp_mat)
  }
  
  # --------------------------------------------------------------------------------------------
  # Step 4: Construct and return a BayesianTools prior object
  # --------------------------------------------------------------------------------------------
  prior <- createPrior(
    density = dens_fun,   # Log-density function
    sampler = samp_fun,   # Sampling function
    lower   = bounds$inf, # Lower bounds for each parameter
    upper   = bounds$sup  # Upper bounds for each parameter
  )
  
  return(prior)
}
