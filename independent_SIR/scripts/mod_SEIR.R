#######################################################################################################
# This script create the SEIR modal and stores the parameters
#######################################################################################################
setwd("/Users/kiran/Dropbox/Max_Planck/Hitchhiker")
library(pomp)
skel_SEIR <- "
  // for the seasonality term 
  double omega = 2.0 * M_PI / 365.0;
  double seas_t = 1.0 + b_amp * cos(omega * (t - b_pha)); // Seasonal term
  double Rn = Ri / (1.0 - R0); // Basic reproduction number of virus 1 
  double p1 = I/N;
  double lambda = Rn * gamma * p1 * seas_t; // Force of infection
  
  DS = -lambda * S;           // dS/dt
  DE = lambda * S - sigma * E;
  DI = sigma * E - gamma * I;   // dI/dt
  DR = gamma * I;             // dR/dt 

  // Accumulator variables
  DH = rho * gamma * I;       // Accumulation of H
"


init_SEIR_1 <- "
  S = nearbyint(N - E0 - R0); // S(0)
  E = nearbyint(E0); // E(0)
  I = 0; // I(0)
  R = nearbyint(R0); // R(0)
  H = 0; // H(0)
"

init_SEIR_2 <- "
  S = nearbyint(N - E0 - R0); // S(0)
  E = nearbyint(E0); // E(0)
  I = 0; // I(0)
  R = nearbyint(R0); // R(0)
  H = 0; // H(0)
"

# init_SEIR_2 <- "
#   S = N - E0 - R0;  
#   E = E0;
#   I = 0;
#   R = R0;
#   H = 0; "


## The below code is for the negative binomial for measuement model.
r_meas <- Csnippet("
  cases = rnbinom_mu(1.0 / psi, H);
")
d_meas <- Csnippet("
  double logL = ISNA(cases) ? 0.0 : dnbinom_mu(nearbyint(cases), 1.0 / psi, H, 1);
  lik = (give_log) ? logL : exp(logL);
")



#### The below code is for the normal distribution for measurement model.
# r_meas <- Csnippet("
#   cases = rnorm(H, psi);
# ")
# 
# d_meas <- Csnippet("
#   double logL = ISNA(cases) ? 0.0 : dnorm(cases, H, psi, 1);
#   lik = (give_log) ? logL : exp(logL);
# ")



pop <- 1e6
SEIR_mod1 <- pomp(
  #data = data.frame(time = obs_data$time, cases = obs_data$H_2), # Use observed cases
  data = data.frame(time = seq(from = 0.0, to = 400, by = 1), cases = NA),
  times = "time", # Name of time variable in data
  t0 = 0, 
  obsnames = "cases",
  statenames = c("S", "E", "I", "R", "H"), # Names of state variables
  accumvars = c("H"),
  paramnames = c("Ri", "b_amp", "b_pha", "sigma", "gamma", "N", "E0", "R0", "rho", "psi"),
  params =  c(
    Ri =2.5,
    b_amp = 0.0, #0.5, 
    b_pha = 0.0, #182.5,
    sigma = 1 / 4,
    gamma = 1 / 5,
    N = pop,
    E0 = 1,  #1e-3 * pop,
    R0 = 0,
    rho = 1, 
    psi = 1e-10),
  rinit = Csnippet(init_SEIR_1), # Initializer for state variables
  skeleton = vectorfield(Csnippet(skel_SEIR)), # Deterministic skeleton for ODE-based model
  dmeasure = d_meas,  # Measurement function
  rmeasure = r_meas,  # Simulation function
  partrans = parameter_trans(log = c("Ri", "b_amp", "b_pha", "sigma", "gamma", "N", "E0","R0"), 
                             logit = c("rho", "psi")) # Parameter transformation
)

params_v1 <- c(
  Ri = 2.5,
  b_amp = 0.0, #0.5, 
  b_pha = 0.0, #182.5,
  sigma = 1 / 4,
  gamma = 1 / 5,
  N = pop,
  E0 = 1,
  R0 = 0,
  rho = 1, 
  psi = 1e-10
)

params_v2 <- c(
  Ri = 2.5,
  b_amp = 0.0, 
  b_pha = 0.0,
  sigma = 1 / 4,
  gamma = 1 / 5,
  N = pop,
  E0 = 1,
  R0 = 0,
  rho = 0.01, 
  psi =1e-10
)






