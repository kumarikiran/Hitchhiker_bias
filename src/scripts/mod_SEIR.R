#######################################################################################################
# SEIR matrix model setup in pomp
# - Builds the deterministic skeleton (ODEs) and measurement model
# - Defines state initialization, parameters, and parameter transformations
# - Produces a pomp object ready for simulation/inference
#######################################################################################################

# Housekeeping
rm(list = ls())

# Resolve project-rooted paths; load base packages (installs if missing)
here::here()
source(here("scripts", "base_packages.R"))

# ---------------------------- Deterministic skeleton (ODEs) ----------------------------
# This C snippet defines:
#  - Seasonality terms for each virus (cosine forcing)
#  - Basic reproduction numbers (transformed/parameterized)
#  - Forces of infection (lambda1, lambda2)
#  - Flows between the 16 epidemiological compartments (X_*)
#  - Accumulator variables for hospitalizations (X_H*) and total hospitalizations (X_Ht)
#
# Notation:
#  - Subscripts 1 and 2 refer to virus 1 and virus 2
#  - Columns (S, E, I, R) correspond to status w.r.t. the other virus (matrix layout)
#  - theta controls interaction during infectiousness (e.g., facilitation)
#  - p1, p2 are prevalences used in force of infection
skel_SIR_matrix <- "
  // ---- Seasonality (cosine forcing) ----
  double omega1 = 2.0 * M_PI / 365.0;
  double seas_t1 = 1.0 + b_amp1 * cos(omega1 * (t - b_pha1)); // Virus 1 seasonal multiplier
  
  double omega2 = 2.0 * M_PI / 365.0;
  double seas_t2 = 1.0 + b_amp2 * cos(omega2 * (t - b_pha2)); // Virus 2 seasonal multiplier
      
  // ---- Parameter scaling (convenience reparameterization) ----
  double Ri1 = 10 * Ri1p;     
  double Ri2 = 10 * Ri2p;     
  double theta = 10 * thetap; // scaled interaction strength

  // ---- Effective reproduction numbers given initial immunity fractions (R10, R20) ----
  double Rn1 = Ri1 / (1.0 - R10); // virus 1
  double Rn2 = Ri2 / (1.0 - R20); // virus 2 
  
  // ---- Prevalence proxies (fraction infectious w.r.t. each virus) ----
  double p1 = (X_IS + X_IE + X_II +  X_IR)/N; // prevalence for virus 1
  double p2 = (X_SI + X_EI + X_II +  X_RI)/N; // prevalence for virus 2
  
  // ---- Forces of infection (FOI), with seasonality ----
  double lambda1 = Rn1 * gamma1 * p1 * seas_t1; // FOI for virus 1
  double lambda2 = Rn2 * gamma2 * p2 * seas_t2; // FOI for virus 2
  
  // ===================== Differential equations (matrix layout) =====================

  // Column: S_* (susceptible to virus 2)
  DX_SS = -(lambda1 + lambda2) * X_SS; 
  DX_ES = lambda1 * X_SS  -(sigma1 + lambda2) * X_ES; 
  DX_IS = sigma1 * X_ES - (gamma1 + theta * lambda2) * X_IS; 
  DX_RS = gamma1 * X_IS - lambda2 * X_RS; 
  
  // Column: E_* (exposed to virus 2)
  DX_SE = lambda2 * X_SS - (lambda1 + sigma2) * X_SE; 
  DX_EE = lambda2 * X_ES + lambda1 * X_SE - (sigma1 + sigma2) * X_EE; 
  DX_IE = (theta * lambda2) * X_IS + sigma1 * X_EE - (gamma1 + sigma2) * X_IE ; 
  DX_RE = lambda2 * X_RS + gamma1 * X_IE - sigma2 * X_RE;  
  
  // Column: I_* (infectious with virus 2)
  DX_SI = sigma2 * X_SE - (theta * lambda1 + gamma2) * X_SI; 
  DX_EI = sigma2 * X_EE + (theta * lambda1) * X_SI - (sigma1 + gamma2) * X_EI; 
  DX_II = sigma2 * X_IE + sigma1 * X_EI - (gamma1 + gamma2) * X_II; 
  DX_RI = sigma2 * X_RE + gamma1 * X_II - gamma2 * X_RI; 
  
  // Column: R_* (recovered from virus 2)
  DX_SR = gamma2 * X_SI - lambda1 * X_SR; 
  DX_ER = gamma2 * X_EI + lambda1 * X_SR - sigma1 * X_ER;
  DX_IR = gamma2 * X_II + sigma1 * X_ER - gamma1 * X_IR; 
  DX_RR = gamma2 * X_RI + gamma1 * X_IR;
  
  // ===================== Accumulators for hospitalizations =====================
  // X_H1, X_H2: hospitalizations attributable to virus 1 and 2 (true causes)
  // X_H1o, X_H2o: 'observed' hospitalizations including hitchhiker contributions via co-infection
  // X_Ht: total hospitalizations (background + v1 + v2)
  DX_H1  = (rho1 * gamma1 * p1 * N);  
  DX_H2  = (rho2 * gamma2 * p2 * N);   
  DX_H1o = (rho1 * gamma1 * p1 * N) + (1-rho1)*rho2*gamma2*X_II;  // virus 2 hitchhiking onto v1
  DX_H2o = (rho2 * gamma2 * p2 * N) + (1-rho2)*rho1*gamma1*X_II;  // virus 1 hitchhiking onto v2
  DX_Ht  = (alpha * N) + (rho1 * gamma1 * p1 * N) + (rho2 * gamma2 * p2 * N); // background + both
"

# ---------------------------- Initial values at t = 0 ----------------------------
# Initializes the 16 state variables and accumulators.
# E10/E20: initial exposed counts to virus 1/2; R10/R20: initial recovered fractions to virus 1/2.
init_SIR_matrix <- "
  X_SS = nearbyint(N - E10 - E20 - R10 - R20) ; 
  X_ES = E10;
  X_IS = 0; 
  X_RS = R10; 
  
  X_SE = E20; 
  X_EE = 0; 
  X_IE = 0; 
  X_RE = 0; 
  
  X_SI = 0; 
  X_EI = 0; 
  X_II = 0; 
  X_RI = 0; 
  
  X_SR = R20; 
  X_ER = 0; 
  X_IR = 0;
  X_RR = 0;
  
  X_H1 = 0;
  X_H2 = 0;
  X_Ht = 0;
  X_H1o = 0;
  X_H2o = 0;
"

# ---------------------------- Measurement model (Negative Binomial) ----------------------------
# r_meas: simulates observed hospitalizations (cases_v1, cases_v2) given latent accumulators.
# d_meas: evaluates the log-likelihood of observed hospitalizations under NB(mu, size=1/psi).
# - 'hitchhikers' flag: if > 0.5, virus 2 observations use X_H2o (hitchhiker-inclusive); else X_H2.
# - psi is the overdispersion parameter. 
r_meas <- Csnippet("
  if (hitchhikers > 0.5) {
    cases_v1 = rnbinom_mu(1.0 / psi, X_H1);
    cases_v2 = rnbinom_mu(1.0 / psi, X_H2o);
  } else {
    cases_v1 = rnbinom_mu(1.0 / psi, X_H1);
    cases_v2 = rnbinom_mu(1.0 / psi, X_H2);
  }
")

d_meas <- Csnippet("
  double logL1, logL2;

  // Handle missing values (NA => contribute zero to log-likelihood)
  if (ISNA(cases_v1)) {
    logL1 = 0.0;
  } else {
    logL1 = dnbinom_mu(nearbyint(cases_v1), 1.0 / psi, X_H1, 1);
  }

  if (ISNA(cases_v2)) {
    logL2 = 0.0;
  } else {
    if (hitchhikers > 0.5) {
      logL2 = dnbinom_mu(nearbyint(cases_v2), 1.0 / psi, X_H2o, 1);
    } else {
      logL2 = dnbinom_mu(nearbyint(cases_v2), 1.0 / psi, X_H2, 1);
    }
  }

  // Return total log-likelihood (or likelihood if give_log == 0)
  lik = (give_log) ? (logL1 + logL2) : exp(logL1 + logL2);
")

# ---------------------------- Build pomp object ----------------------------
pop <- 1e6  # total population size

SIR_matrix_mod <- pomp(
  # Stub data frame so we can simulate/fit; NA means no observations yet
  data = data.frame(time = seq(from = 0.0, to = 500, by = 1),
                    cases_v1 = NA, cases_v2 = NA),
  times = "time",
  t0 = 0,
  
  # Observables (measurement variables)
  obsnames = c("cases_v1", "cases_v2"),
  
  # State variables (16 epidemiological states + 5 accumulators)
  statenames =  c("X_SS", "X_ES", "X_IS", "X_RS", 
                  "X_SE", "X_EE", "X_IE", "X_RE", 
                  "X_SI", "X_EI", "X_II", "X_RI", 
                  "X_SR", "X_ER", "X_IR", "X_RR",
                  "X_H1", "X_H2",
                  "X_H1o", "X_H2o", "X_Ht"),
  
  # Accumulators are integrated over time (non-reset between steps)
  accumvars = c("X_H1", "X_H2", "X_H1o", "X_H2o", "X_Ht"),
  
  # Parameters (naming matches usage in skeleton/init/measurement code)
  paramnames = c("Ri1p", "Ri2p",        # reproduction number proxies (scaled inside skeleton)
                 "b_amp1", "b_amp2",    # seasonal amplitude
                 "b_pha1", "b_pha2",    # seasonal phase shifts
                 "sigma1", "sigma2",    # progression rates E->I
                 "gamma1", "gamma2",    # recovery rates I->R
                 "thetap",              # interaction strength proxy (scaled inside skeleton)
                 "rho1", "rho2",        # hospitalization probabilities (given infection)
                 "alpha",               # background hospitalization rate (per capita)
                 "N",                   # population size
                 "E10", "E20",          # initial exposed counts
                 "R10", "R20",          # initial recovered fractions
                 "psi",                 # NB overdispersion
                 "hitchhikers"),        # flag to toggle hitchhiker-inclusive observations
  
  # Default parameter values (can be overwritten at simulate/fit time)
  params = c(Ri1p = 2.5/10,
             Ri2p = 2.5/10, 
             b_amp1 = 0.0, 
             b_amp2 = 0.0, 
             b_pha1 = 0.0,
             b_pha2 = 0.0,
             sigma1 = 1 / 4, 
             sigma2 = 1 / 4, 
             gamma1 = 1 / 5, 
             gamma2 = 1 / 5, 
             thetap = 1/10, 
             alpha = 0, 
             N = pop,
             E10 = 1, 
             E20 = 1, 
             R10 = 0, 
             R20 = 0,
             rho1 = 1, 
             rho2 = 0.01, 
             psi = 0.04, 
             hitchhikers = 1.0),
  
  # Hooks for pomp
  rinit     = Csnippet(init_SIR_matrix),                # state initialization
  skeleton  = vectorfield(Csnippet(skel_SIR_matrix)),   # deterministic ODEs
  dmeasure  = d_meas,                                   # likelihood evaluator
  rmeasure  = r_meas,                                   # simulator for observations
  
  # Parameter transformations for stability/constraints during inference:
  #  - log: parameters constrained to (0, ∞)
  #  - logit: parameters constrained to (0, 1)
  partrans = parameter_trans(
    log   = c("b_amp1", "b_amp2", "b_pha1", "b_pha2", "sigma1", "sigma2",
              "gamma1", "gamma2", "N", "E10","E20", "R10", "R20"),
    logit = c("Ri1p", "Ri2p", "thetap", "rho1", "rho2", "psi", "alpha")
  ),
  
  # C compilation options
  cdir = '.', cfile = 'COVmod'
)

# ---------------------------- Convenience parameter set for simulations ----------------------------
# params_case1 can be passed to trajectory() or simulate() to generate trajectories

coef(SIR_matrix_mod)['E10'] <- 1e-6
coef(SIR_matrix_mod)['E20'] <- 1e-6
params_case1 = SIR_matrix_mod@params

