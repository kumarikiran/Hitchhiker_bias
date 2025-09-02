# SEIR Matrix Model Scripts

This folder contains R scripts for building and simulating a **two-virus SEIR model** with the `pomp` package.

## Files
- `base_packages.R` – installs and loads required packages.  
- `mod_SEIR.R` – defines the model equations, initial conditions, and measurement model.

## Usage
1. Load packages:
   ```r
   source(here::here("src", "scripts", "base_packages.R"))


