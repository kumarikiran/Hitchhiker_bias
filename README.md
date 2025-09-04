# SEIR Matrix Model Scripts

This folder contains R scripts for building and simulating a **two-virus SEIR model** with the `pomp` package.

## Files
- `scripts/base_packages.R` – installs and loads required packages.  
- `scripts/mod_SEIR.R` – defines the model equations, initial conditions, and measurement model.
- `scripts/make_Fig_3.R` – generates the plots for Fig. 3 in the manuscript:  
  - **Fig. 3a**: Daily infection rates (I₁, I₂) under three scenarios (`minimal`, `partial`, `complete`).  
  - **Fig. 3b**: Daily hospitalization rates (H₁, H₂, H₁ᵒ, H₂ᵒ), showing the impact of hitchhiker bias.
- `scripts/simulate.R` – runs simulations across a **grid of parameter values**:  
  - Varies initial exposure (`E10`, 100 log-spaced values between 1e4 and 1e-6).  
  - Varies reporting probability of virus 2 (`rho2`, 10 evenly spaced values between 0.01 and 1).  
  - Produces **1000 simulated trajectories** in total.  
  - Results are saved as `.rds` files in the `output_data/` folder, automatically created if missing.
- `scripts/make_Fig_4.R` – reads simulation outputs from `output_data/` and generates Fig. 4:  
  - **Fig. 4a**: Change in peak hospitalization amplitude vs. `corr(I₁, I₂)`.  
  - **Fig. 4b**: Change in peak hospitalization timing vs. `corr(I₁, I₂)`.  
  - **Fig. 4c–d**: Heatmaps of these changes over `corr(I₁, I₂)` × pathogen pathogenicity (`rho2`).
- `scripts/simulations_with_noise.R` – runs SEIR simulations with noise injection for hospitalizations:  
  - Sweeps over three values of initial exposure (`E10`: 1e4, 1e-4, 1e-6).  
  - Sweeps over four reporting probabilities of virus 2 (`rho2`: 0.01, 0.10, 0.25, 0.50).  
  - For each combination, generates **100 noisy replicates** using a negative binomial observation model (`rho_mean = 1`, `rho_k = 0.04`).  
  - Saves combined outputs as `.rds` files in the `output_data_noise_4percent/` folder (created if missing).  
  - Useful for sensitivity analyses of noisy hospitalization data.
- `scripts/fit_set_prior_MCMC.R` – defines Bayesian prior distributions for model parameters (for use with **BayesianTools**):  
  - Parameters supported: `Ri1p`, `rho1`, `Ri2p`, `rho2`, `thetap`, `psi`.  
  - Priors:  
    - `Ri1p`, `Ri2p` ~ Exponential(rate = 0.1)  
    - `rho1`, `rho2`, `thetap` ~ Beta(2, 2)  
    - `psi` ~ Uniform(0, 1)  
  - Returns a prior object with:  
    - **density function** (log-scale) for evaluating priors,  
    - **sampling function** for drawing random samples,  
    - **parameter bounds** (lower/upper). 
- `scripts/fit_run_MCMC.R` – fits the SEIR POMP model to simulated noisy hospitalization data using **Bayesian MCMC**:  
  - Loops over 100 replicated noisy datasets (from `run_simulations.R`).  
  - Builds a `pomp` object with observed `cases_v1` and `cases_v2`.  
  - Estimates six parameters: `Ri1p`, `rho1`, `Ri2p`, `rho2`, `thetap`, `psi`.  
  - Likelihood function created via `traj_objfun` (trajectory matching).  
  - Priors specified with `set_prior.R`.  
  - Runs **Differential Evolution MCMC with snooker update (DEzs)** via `BayesianTools`.  
  - Saves chain objects in `MCMC_fit/` and posterior samples in `MCMC_results/`.  
  - Provides diagnostic outputs (DIC, MAP, MCMC summary) for model assessment.




