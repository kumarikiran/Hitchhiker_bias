#######################################################################################################
# This script create the SEIR modal and stores the parameters
#######################################################################################################
rm(list = ls())
here::here()
source(here("scripts", "base_packages.R"))

# Deterministic skeleton (i.e., set of ordinary differential equations)
skel_SIR_matrix <- "
  double omega1 = 2.0 * M_PI / 365.0;
  double seas_t1 = 1.0 + b_amp1 * cos(omega1 * (t - b_pha1)); // Seasonal term
  
  double omega2 = 2.0 * M_PI / 365.0;
  double seas_t2 = 1.0 + b_amp2 * cos(omega2 * (t - b_pha2)); // Seasonal term
      
  double Ri1 = 10 * Ri1p;
  double Ri2 = 10 * Ri2p;
  double theta = 10 * thetap;

  double Rn1 = Ri1 / (1.0 - R10); // Basic reproduction no of virus 1 
  double Rn2 = Ri2 / (1.0 - R20); // Basic reproduction no of virus 2 
  
  double p1 = (X_IS + X_IE + X_II +  X_IR)/N; // Prevalence of infection with virus 1
  double p2 = (X_SI + X_EI + X_II +  X_RI)/N; // Prevalence of infection with virus 2
  
  double lambda1 = Rn1 * gamma1 * p1 * seas_t1; // Force of infection with virus 1
  double lambda2 = Rn2 * gamma2 * p2 * seas_t2; // Force of infection with virus 2
  
  // Differential equations
  
  // Susceptible (S) to virus 2 (1st column in model schematic)
  DX_SS = -(lambda1 + lambda2) * X_SS; 
  DX_ES = lambda1 * X_SS  -(sigma1 + lambda2) * X_ES; 
  DX_IS = sigma1 * X_ES - (gamma1 + theta * lambda2) * X_IS; 
  DX_RS = gamma1 * X_IS - lambda2 * X_RS; 
  
  // Infected (E) with virus 2 (2nd column in model schematic)
  DX_SE = lambda2 * X_SS - (lambda1 + sigma2) * X_SE; 
  DX_EE = lambda2 * X_ES + lambda1 * X_SE - (sigma1 + sigma2) * X_EE; 
  DX_IE = (theta * lambda2) * X_IS + sigma1 * X_EE - (gamma1 + sigma2) * X_IE ; 
  DX_RE = lambda2 * X_RS + gamma1 * X_IE - sigma2 * X_RE;  
  
  // Infected (I) with virus 2 (3nd column in model schematic)
  DX_SI = sigma2 * X_SE - (theta * lambda1 + gamma2) * X_SI; 
  DX_EI = sigma2 * X_EE + (theta * lambda1) * X_SI - (sigma1 + gamma2) * X_EI; 
  DX_II = sigma2 * X_IE + sigma1 * X_EI - (gamma1 + gamma2) * X_II; 
  DX_RI = sigma2 * X_RE + gamma1 * X_II - gamma2 * X_RI; 
  
  // Recovered (R) with virus 2 (4th column in model schematic)
  DX_SR = gamma2 * X_SI - lambda1 * X_SR; 
  DX_ER = gamma2 * X_EI + lambda1 * X_SR - sigma1 * X_ER;
  DX_IR = gamma2 * X_II + sigma1 * X_ER - gamma1 * X_IR; 
  DX_RR = gamma2 * X_RI + gamma1 * X_IR;
  
  // Accumulator variables
  DX_H1 = (rho1 * gamma1 * p1 * N);  
  DX_H2 = (rho2 * gamma2 * p2 * N);   
  DX_H1o = (rho1*gamma1*p1*N) + (1-rho1)*rho2*gamma2*X_II;
  DX_H2o = (rho2*gamma2*p2*N) + (1-rho2)*rho1*gamma1*X_II;
  DX_Ht = (alpha * N) + (rho1 * gamma1 * p1 * N) + (rho2 * gamma2 * p2 * N); 
"

# Initializer (values of state variables at time t = 0)
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
"



## The below code is for the negative binomial for measuement model.
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

  // Handle missing values
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

  // Combine log-likelihoods
  lik = (give_log) ? (logL1 + logL2) : exp(logL1 + logL2);
")



pop <- 1e6 ## Total number of population
SIR_matrix_mod <- pomp(
  data = data.frame(time = seq(from = 0.0, to = 500, by = 1), cases_v1 = NA, cases_v2 =NA),
  times = "time", # Name of time variable in data
  t0 = 0, 
  obsnames = c("cases_v1", "cases_v2"),
  statenames =  c("X_SS", "X_ES", "X_IS", "X_RS", # Names of state variables
                  "X_SE", "X_EE", "X_IE", "X_RE", 
                  "X_SI", "X_EI", "X_II", "X_RI", 
                  "X_SR", "X_ER", "X_IR",  "X_RR",
                  "X_H1", "X_H2",
                  "X_H1o", "X_H2o", "X_Ht"), 
  accumvars = c("X_H1", "X_H2", "X_H1o", "X_H2o", "X_Ht"),
  paramnames = c("Ri1p", "Ri2p", # Initial reproduction numbers
                 "b_amp1", "b_amp2",
                 "b_pha1", "b_pha2",
                 "sigma1", "sigma2",
                 "gamma1", "gamma2", 
                 "thetap",  # Strengths of interaction during infectious period
                 "rho1", "rho2",
                 "alpha",
                 "N",
                 "E10", "E20", # Initial fractions infected 
                 "R10", "R20", # Initial fractions recovered
                 "psi", 
                 "hitchhikers"),
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
  rinit = Csnippet(init_SIR_matrix), # Initializer for state variables
  skeleton = vectorfield(Csnippet(skel_SIR_matrix)), # Deterministic skeleton for ODE-based model
  dmeasure = d_meas,  # Measurement function
  rmeasure = r_meas,  # Simulation function
  partrans = parameter_trans(log = c("b_amp1", "b_amp2", "b_pha1", "b_pha2", "sigma1", "sigma2",
                                      "gamma1", "gamma2", "N", "E10","E20", "R10", "R20"), 
                             logit = c("Ri1p", "Ri2p", "thetap", "rho1", "rho2", "psi", "alpha")), # Parameter transformation
  cdir = '.', cfile = 'COVmod'
)

params_case1 = c(Ri1p = 2.5/10,
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
           alpha = 0.0, 
           N = pop,
           E10 = 1e-6, 
           E20 = 1e-6, 
           R10 = 0, 
           R20 = 0,
           rho1 = 1, 
           rho2 = 0.01, 
           psi = 0.04, 
           hitchhikers = 1.0)


