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


## Usage
1. Load packages:
   ```r
   source(here::here("src", "scripts", "base_packages.R"))


