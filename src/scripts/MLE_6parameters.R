# --- Setup ---
rm(list = ls())                       # Clear the R workspace
library(here)
here::here()                          # Set project root path (using {here})
source(here("scripts", "mod_SEIR.R"))       # Load SEIR model definition

library(nloptr)
theme_set(theme_bw()) # Set theme for ggplot2
par(mfrow = c(1, 1), lty = 1, bty = "l") # Set the graphical parameters 

i = 1
# Observed hospitalization data
data_file <- sprintf("output_data/tj_combined_rho0.01_E1_%s.rds", i)
obs_data <- readRDS(data_file)

data_for_pomp <- data.frame(
  time = obs_data$time,   
  cases_v1 = obs_data$X_H1, 
  cases_v2 = obs_data$X_H2o    #### the original data 
)

library(ggplot2)
# Create the plot
plot_1 <- ggplot(data = data_for_pomp, aes(x = time)) +
  geom_line(aes(y = cases_v1), color = "blue", size = 1, linetype = "solid") +  # Line for cases_2
  geom_line(aes(y = cases_v2), color = "red", size = 1, linetype = "dashed") +  # Line for cases_1
  labs(title = "Cases Over Time",
       x = "Time",
       y = "Cases") +
  theme_minimal()
print(plot_1)

# Define POMP object
SEIR_mod_fit <- pomp(
  data = data.frame(time = data_for_pomp$time, cases_v1 = data_for_pomp$cases_v1, cases_v2 = data_for_pomp$cases_v2),
  #data = data.frame(time = seq(from = 0.0, to = 400, by = 1), cases = NA),
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
                 "thetap", # Strengths of interaction during infectious period
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
             E10 = 10000, 
             E20 = 0.000001 , 
             R10 = 0, 
             R20 = 0,
             rho1 = 1, 
             rho2 = 0.01, 
             psi = 1e-10, 
             hitchhikers = 1.0),
  rinit = Csnippet(init_SIR_matrix), # Initializer for state variables
  skeleton = vectorfield(Csnippet(skel_SIR_matrix)), # Deterministic skeleton for ODE-based model
  dmeasure = d_meas,  # Measurement function
  rmeasure = r_meas,  # Simulation function
  partrans = parameter_trans(log = c( "b_amp1", "b_amp2", "b_pha1", "b_pha2", "sigma1", "sigma2",
                                      "gamma1", "gamma2",  "N", "E10","E20", "R10", "R20"), 
                             logit = c("Ri1p", "Ri2p", "thetap", "rho1", "rho2", "psi", "alpha")), # Parameter transformation
  cdir = '.', cfile = 'COVmod'
)

# Save as .rds file
saveRDS(SEIR_mod_fit, file = here("_saved", "pomp_mod_Hon_01.rds"))


estpars_nm <- c("Ri1p", "rho1", "Ri2p", "rho2", "thetap", "psi") 


start_values <- sobol_design(lower = c("Ri1p" = 1/10, "rho1" = 0, "Ri2p" = 1/10, "rho2" = 0, "thetap"= 0.1/10, "psi" = 1e-10), 
                             upper = c("Ri1p" = 10/10, "rho1" = 1, "Ri2p" = 10/10, "rho2" = 1, "thetap"= 5/10, "psi" = 1e-5), 
                             nseq = 100)



# Define the negative log-likelihood function
# The partrans argument is used to indicate that the estimation is on the log scale
nLogL <- traj_objfun(data = SEIR_mod_fit, 
                     est = estpars_nm, 
                     partrans = SEIR_mod_fit@partrans,
                     verbose = TRUE)
print(nLogL(coef(SEIR_mod_fit, estpars_nm, transform = T)))  


###### below is the modified code for the mac 
library(nloptr)
library(doMC)
library(foreach)
tm <- bake(file = "_saved/parameter_6_Hoff_3.rds", 
           expr = {
             # Set up parallel processing
             numCores <- min(parallel::detectCores() - 1, 4)  # Use up to 4 cores, or one less than available
             registerDoMC(numCores)
             
             print(getDoParWorkers())
             
             results <- foreach(r = iter(start_values, by = "row"), 
                                .combine = list,
                                .multicombine = TRUE,
                                .inorder = FALSE, 
                                .errorhandling = "pass", 
                                .verbose = TRUE) %dopar% {
                                  
                                  tryCatch({
                                    nloptr::nloptr(x0 = log(unname(unlist(r))), # Starting parameter values (on log-scale) 
                                                   eval_f = nLogL, # Function to minimize (negative log-likelihood)
                                                   opts = list(algorithm = "NLOPT_LN_SBPLX", # Name of the optimizer (subplex algorithm)
                                                               maxtime = 500, # Maximum execution time (in seconds)
                                                               maxeval = -1))
                                  }, error = function(e) {
                                    list(error = e$message)
                                  })
                                }
             
             results  # Return the results
           })
print(attr(tm, "system.time")) # Print how long the estimations took

# Extract convergence status and messages
convergence_info <- lapply(tm, function(res) {
  if (!is.null(res$objective)) {
    list(status = res$status, message = res$message)
  } else {
    list(status = NA, message = res$error)  # error caught in tryCatch
  }
})

# Convert to data.frame for easier inspection
convergence_df <- do.call(rbind, lapply(convergence_info, as.data.frame))

# View status codes
table(convergence_df$status)

successful_tm <- tm[sapply(tm, function(x) !is.null(x$objective) && x$status == 4)]

solutions <- sapply(successful_tm, getElement, "solution") %>% t()
# Create a function to apply the model's transformation
transform_params <- function(params) {
  temp_params <- params
  coef(SEIR_mod_fit, estpars_nm, transform = TRUE) <- temp_params
  return(coef(SEIR_mod_fit, estpars_nm))
}
# Apply the transformation to each row of solutions
pars_est <- t(apply(solutions, 1, transform_params)) %>% as.data.frame()
colnames(pars_est) <- estpars_nm


pars_est <- pars_est %>% 
  mutate(ll = -sapply(successful_tm, getElement, "objective")) %>% # Extract log-likelihood
  arrange(desc(ll)) # Arrange by decreasing log-likelihood (best parameters on top of data frame)
print(pars_est)
# Save the data frame to an RDS file
saveRDS(pars_est, file = here("_saved", "fit_MLE_p6_Hoff_3.rds"))

#### Transform back to original 
#pars_est$Ri1 <- 10 * pars_est$Ri1p
#pars_est$Ri2 <- 10 * pars_est$Ri2p
#pars_est$theta <- 10 * pars_est$thetap
#pars_est <- pars_est[, !(names(pars_est) %in% c("Ri1p", "Ri2p", "thetap"))]
#pars_est <- pars_est[, c("Ri1", "Ri2", "rho1", "rho2", "theta", "ll")] 
#print(pars_est)

# Need to ensure that we have reached the maximum likelihood estimate (MLE)
# Pair plot of estimates
# All parameters seem well identified
pairs(x = pars_est %>% filter(max(ll) - ll <= 20) %>% select(all_of(estpars_nm), "ll"))
### plotting all the data
pairs(x = pars_est %>% select(all_of(estpars_nm), "ll"))


# Define a custom panel function to show correlation
# Custom panel function for correlation with larger font
panel.cor <- function(x, y, digits = 2, prefix = "", ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "complete.obs")
  txt <- formatC(r, digits = digits, format = "f")
  # Set large font size here with cex
  text(0.5, 0.5, txt, cex = 2.0, font = 1)  # Bold & bigger
}

# Prepare data
plot_data <- pars_est %>%
  filter(max(ll) - ll <= 20) %>%
  select(all_of(estpars_nm), ll)

# Plot with custom upper.panel
pairs(
  x = plot_data,
  upper.panel = panel.cor,
  lower.panel = panel.smooth,  # Optional
  pch = 19,
  cex = 0.6
)





mle <- setNames(object = as.numeric(pars_est[1, estpars_nm]), 
                nm = estpars_nm) # MLE

#install.packages("SlicedLHD")
#library(SlicedLHD)
slices <- slice_design(center = mle, 
                       Ri1p = seq(from = 0.5 * mle["Ri1p"], to = 1.5 * mle["Ri1p"], length.out = 20), 
                       rho1 = seq(from = 0.5 * mle["rho1"], to = 1.5 * mle["rho1"], length.out = 20),
                       Ri2p = seq(from = 0.5 * mle["Ri2p"], to = 1.5 * mle["Ri2p"], length.out = 20), 
                       rho2 = seq(from = 0.5 * mle["rho2"], to = 1.5 * mle["rho2"], length.out = 20),
                       thetap = seq(from = 0.5 * mle["thetap"], to = 1.5 * mle["thetap"], length.out = 20), 
                       psi = seq(from = 0.5 * mle["psi"], to = 1.5 * mle["psi"], length.out = 20))%>% 
  mutate(ll = NA)



# Now we can compute the log-likelihood
for(r in 1:nrow(slices)) {
  slices$ll[r] <- -nLogL(par = log(as.numeric(slices[r, estpars_nm])))
}

# Plot the sliced log-likelihood and check that we have reached the maximum
par(mfrow = c(3, 2), bty = "l", las = 2, lwd = 2) # Set a single plot

for(par in estpars_nm) {
  slices_cur <- filter(slices, slice == par)
  plot(slices_cur[[par]], slices_cur$ll, 
       type = "p",  # 'p' for points
       pch = 10,    # Circle shape
       bg = "blue", # Fill color of circles
       col = "black", # Border color of circles
       cex = 0.5,   # Size of circles
       xlab = par, 
       ylab = "Log-Likelihood", 
       main = par)
  abline(v = mle[par], lty = 2, col = "red") # Add vertical line at MLE
  points(mle[par], max(slices_cur$ll), pch = 21, bg = "red", cex = 1.5) # Highlight MLE point
}


# EXERCISE:  MODEL INTERPRETATION AND EVALUATION -----------------------------------------------------------
# As a first evaluation, let's look at the R_effective over time
coef(SEIR_mod_fit, names(mle)) <- unname(mle) # Set the parameters to their MLE

sims_E3 <- trajectory(object = SEIR_mod_fit, 
                      format = "data.frame")

# Now let's compare the data simulated at the MLE to the observed data
# Plot (superimpose the data in red)


pl_E3 <- ggplot(data = sims_E3, 
                mapping = aes(x = time)) + 
  
  # Modeled H2
  geom_line(aes(y = X_H2o, color = "Modeled H2", linetype = "Modeled H2"), size = 1) + 
  
  # Observed H2
  geom_line(data = data_for_pomp, 
            mapping = aes(x = time, y = cases_v2, color = "Observed H2", linetype = "Observed H2"), 
            size = 1.5) + 
  
  # Modeled H1
  geom_line(aes(y = X_H1, color = "Modeled H1", linetype = "Modeled H1"), size = 1) +
  
  # Observed H1
  geom_line(data = data_for_pomp, 
            mapping = aes(x = time, y = cases_v1, color = "Observed H1", linetype = "Observed H1"), 
            size = 1.5) +
  
  labs(x = "Time (days)", 
       y = "Number of Hospitalizations", 
       title = "Observed vs. Modeled Hospitalizations",
       color = "Data",
       linetype = "Data") +
  
  scale_color_manual(values = c(
    "Modeled H2" = "blue", 
    "Observed H2" = "red", 
    "Modeled H1" = "darkgreen", 
    "Observed H1" = "orange"
  )) +
  
  scale_linetype_manual(values = c(
    "Modeled H2" = "solid", 
    "Observed H2" = "dashed", 
    "Modeled H1" = "solid", 
    "Observed H1" = "dashed"
  )) +
  
  scale_x_continuous(limits = c(0, 400)) +  # <-- Set X-axis limits
  
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 14) )

print(pl_E3)
pdf("compare.pdf", width = 5, height = 4); 
print(pl_E3)
dev.off()











