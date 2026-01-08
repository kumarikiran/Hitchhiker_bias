# ---------------------- Setup ----------------------
rm(list = ls())  # Start with a clean workspace to avoid object/name collisions

# Load required libraries
library(pomp)          # State-space / POMP modeling framework
library(BayesianTools) # Bayesian inference utilities (priors, MCMC)
library(MCMCvis)       # Convenient MCMC diagnostics/summary/plots
library(coda)          # MCMC objects and diagnostics

# Local project sources:
source(here("scripts", "mod_SEIR.R"))            # Defines model states, flows, measurement, and C snippets
source(here("scripts", "fit_set_prior_MCMC.R"))  # Provides set_prior() for 6 params used below

# ---------------------- Unique compile directory (for C snippets) ----------------------
# Use SLURM job id if available; otherwise fall back to the current R process id.
job_id <- Sys.getenv("SLURM_JOB_ID")
if (job_id == "") job_id <- as.character(Sys.getpid())

# Temporary compilation directory for pomp's C snippets (avoid collisions on clusters)
compile_dir <- file.path(tempdir(), paste0("compile_tmp_", job_id))
dir.create(compile_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------- Parse command line arguments ----------------------
# Args:
#   1) i      : index selecting minial(1), partial (2) and complete(3) overlap value (default 3)
#   2) rho_i  : reporting probability for virus 2 used to locate input file (default 0.01)
args <- commandArgs(trailingOnly = TRUE)
i     <- ifelse(length(args) >= 1, as.numeric(args[1]), 3)
rho_i <- ifelse(length(args) >= 2, as.numeric(args[2]), 0.01)

# ---------------------- Experiment configuration ----------------------
# E10 candidates; select via index i
data_points <- c(10000, 0.0001, 0.000001)

# Hitchhiker bias switch/level (model-specific; defined in mod_SEIR.R)
hitch_i <- 1.0

# Suffix used in saved filenames to track which E10 index was fitted
suffix_save <- sprintf("Hon_%d", i)

# ---------------------- Fit loop over replicated datasets ----------------------
# Expect files produced by run_simulations.R with 100 replicated noisy datasets per setting.
for (replica in 1:100) {
  # ---------- Load one replicate of observed hospitalization data ----------
  data_file <- sprintf("output_data_noise_4percent/tj_combined_rho%s_E1_%s.rds", rho_i, i)
  full_data <- readRDS(data_file)
  obs_data  <- subset(full_data, rep == replica)
  
  # Prepare the observed series for pomp: two observation variables over time
  data_for_pomp <- data.frame(
    time     = obs_data$time,
    cases_v1 = obs_data$X_H1o_obs,  # observed hospitalizations stream 1
    cases_v2 = obs_data$X_H2o_obs   # observed hospitalizations stream 2
  )
  
  # ---------- Build POMP object ----------
  SEIR_mod_fit <- pomp(
    data      = data_for_pomp,
    times     = "time",
    t0        = 0,
    obsnames  = c("cases_v1", "cases_v2"),
    statenames = c(
      "X_SS", "X_ES", "X_IS", "X_RS",
      "X_SE", "X_EE", "X_IE", "X_RE",
      "X_SI", "X_EI", "X_II", "X_RI",
      "X_SR", "X_ER", "X_IR", "X_RR",
      "X_H1", "X_H2",
      "X_H1o", "X_H2o", "X_Ht"
    ),
    accumvars = c("X_H1", "X_H2", "X_H1o", "X_H2o", "X_Ht"),
    paramnames = c(
      "Ri1p", "Ri2p",                 # baseline reproduction multipliers
      "b_amp1", "b_amp2",             # seasonal amplitude virus 1/2
      "b_pha1", "b_pha2",             # seasonal phase virus 1/2
      "sigma1", "sigma2",             # 1/incubation period virus 1/2
      "gamma1", "gamma2",             # 1/infectious period virus 1/2
      "thetap",                       # interaction strength during infectiousness
      "rho1", "rho2",                 # reporting probabilities
      "alpha",                        # (model-specific)
      "N",                            # population size
      "E10", "E20",                   # initial exposed fractions
      "R10", "R20",                   # initial recovered fractions
      "psi",                          # Euler/overdispersion or similar (model-specific)
      "hitchhikers"                   # hitchhiker bias toggle/level
    ),
    # Baseline parameter values (many are fixed during this fit; see pars_est_nm below)
    params = c(
      Ri1p = 2.5/10,
      Ri2p = 2.5/10,
      b_amp1 = 0.0, b_amp2 = 0.0,
      b_pha1 = 0.0, b_pha2 = 0.0,
      sigma1 = 1 / 4, sigma2 = 1 / 4,
      gamma1 = 1 / 5, gamma2 = 1 / 5,
      thetap = 1/10,
      alpha = 0,
      N = pop,
      E10 = data_points[i], E20 = 0.000001,
      R10 = 0, R20 = 0,
      rho1 = 1, rho2 = 0.01,
      psi = 1e-10,
      hitchhikers = hitch_i
    ),
    rinit     = Csnippet(init_SIR_matrix),               # state initializer
    skeleton  = vectorfield(Csnippet(skel_SIR_matrix)),  # deterministic skeleton (for traj_objfun)
    dmeasure  = d_meas,                                  # obs density
    rmeasure  = r_meas,                                  # obs simulator
    # Parameter transformations for optimization/MCMC (keep parameters in valid domains)
    partrans  = parameter_trans(
      log  = c("b_amp1", "b_amp2", "b_pha1", "b_pha2",
               "sigma1", "sigma2", "gamma1", "gamma2",
               "N", "E10", "E20", "R10", "R20"),
      logit = c("Ri1p", "Ri2p", "thetap", "rho1", "rho2", "psi", "alpha")
    ),
    cdir  = compile_dir,   # unique compile directory
    cfile = "COVmod"       # base name for generated C code
  )
  
  mod <- SEIR_mod_fit  # alias used below
  
  # ---------------------- Parameters to estimate & initialization ----------------------
  # We estimate 6 parameters; others remain fixed at the values above.
  pars_est_nm <- c("Ri1p", "rho1", "Ri2p", "rho2", "thetap", "psi")
  n_pars_est  <- length(pars_est_nm)
  
  # Initial values for the estimated parameters (natural scale expected by traj_objfun)
  my_init <- c(Ri1p = 0.25, rho1 = 1, Ri2p = 0.25, rho2 = 0.01, thetap = 0.1, psi = 10^-10)
  coef(mod, pars_est_nm) <- unname(my_init)
  
  # ---------------------- Likelihood function ----------------------
  # traj_objfun returns (by default) a *negative* log-likelihood function.
  # We wrap it so LL_natScale returns a *log-likelihood* as required by BayesianTools.
  minusLL_natScale <- traj_objfun(data = mod, est = pars_est_nm, partrans = NULL)
  LL_natScale <- function(x) -minusLL_natScale(x)
  
  # ---------------------- Priors ----------------------
  # Vague/neutral priors defined in fit_set_prior_MCMC.R for the 6 parameters
  #source(here("scripts", "fit_set_prior_MCMC.R"))
  prior_vague <- set_prior(pars_nm = pars_est_nm)
  
  # ---------------------- MCMC settings ----------------------
  n_chains <- 5           # number of concurrent chains for DEzs
  n_iter   <- 2e4         # iterations per chain (burn-in handled when extracting samples)

  
  # Custom parameter ranges to generate proposal archive Z for DEzs (Differential Evolution sampler)
  ranges_custom <- matrix(c(
    0.1,   1,     # Ri1p
    0.001, 1,     # rho1
    0.1,   1,     # Ri2p
    0.001, 1,     # rho2
    0.01, 0.5,  # thetap
    1e-12, 1      # psi
  ), ncol = 2, byrow = TRUE)
  colnames(ranges_custom) <- c("min", "max")
  rownames(ranges_custom) <- pars_est_nm
  
  # Starting values for each chain (replicate my_init)
  #start_chains <- replicate(n = n_chains, expr = my_init) %>% t()
  ### random ±25% perturbation from the start value + clamp to valid ranges
  start_chains <- t(replicate(n_chains, {
    perturbed <- my_init * (1 + runif(length(my_init), min = -0.25, max = 0.25))
    
    # clamp to parameter ranges
    pmin(pmax(perturbed, ranges_custom[, "min"]), ranges_custom[, "max"])
  }))
  
  colnames(start_chains) <- names(my_init)
  
  # Archive Z: random draws within the specified ranges (500 draws per parameter)
  Z_mat <- matrix(
    runif(500 * n_pars_est,
          rep(ranges_custom[, "min"], each = 500),
          rep(ranges_custom[, "max"], each = 500)),
    ncol = n_pars_est, byrow = FALSE
  )
  
  # ---------------------- BayesianTools setup ----------------------
  bayesianSetup <- createBayesianSetup(
    likelihood = LL_natScale,
    prior      = prior_vague,
    names      = pars_est_nm
    # , best = mle  # (optional) if an MLE were available
  )
  
  # Combine settings for DEzs sampler
  settings <- list(
    iterations = n_chains * n_iter,  # total iterations across chains (DEzs manages multiple chains)
    startValue = start_chains,       # per-chain starting vectors
    Z          = Z_mat               # differential evolution archive
  )
  
  # ---------------------- Run MCMC (cached with bake) ----------------------
  # Saves sampler object to disk; avoids recomputation if file already exists.
  fit_mcmc_det <- bake(
    file = sprintf("MCMC_fit/mcmc_fit_%s_rho%s_r%s.rds", suffix_save, rho_i, replica),
    expr = {
      runMCMC(
        bayesianSetup = bayesianSetup,
        sampler       = "DEzs",   # Differential Evolution with snooker update
        settings      = settings
      )
    }
  )
  
  # ---------------------- Diagnostics & results ----------------------
  # Convert to coda object; discard first half as burn-in when summarizing
  est_mcmc_det <- getSample(sampler = fit_mcmc_det, coda = TRUE, start = n_iter / 2)
  
  # Information criteria & MAP estimate
  dic      <- DIC(sampler = fit_mcmc_det, start = n_iter / 2)
  mle_mcmc <- MAP(fit_mcmc_det)
  
  # Optional trace plotting (disabled)
  # MCMCvis::MCMCtrace(object = est_mcmc_det, pdf = TRUE, ind = TRUE,
  #                    file = "MCMC_fit/my_trace_plot.pdf")
  
  # Console summary (means, SDs, quantiles, Rhat, n_eff)
  MCMCvis::MCMCsummary(object = est_mcmc_det, round = 4)
  
  # Extract posterior samples after burn-in and save as RDS for downstream analysis
  samples_df <- as.data.frame(getSample(fit_mcmc_det, start = n_iter / 2, coda = FALSE))
  saveRDS(samples_df, file = sprintf("MCMC_results/mcmc_results_%s_rho%s_r%s.rds", suffix_save, rho_i, replica))
}
