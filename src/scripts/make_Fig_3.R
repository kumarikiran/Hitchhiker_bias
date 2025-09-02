#######################################################################################################
# This script generates the plots used for Fig. 3 of the manuscript.
# Fig. 3a: Daily infection rate under different initial exposure levels (E10)
# Fig. 3b: Daily hospitalization rate under the same scenarios
#######################################################################################################

# --- Setup ---
rm(list = ls())                       # Clear the workspace
here::here()                          # Ensure project root is set with {here}
source(here("scripts", "base_packages.R"))  # Load required packages
source(here("scripts", "mod_SEIR.R"))       # Load SEIR model definition

# --- Define scenarios ---
E10_values <- c(1e4, 1e-4, 1e-6)                     # Values of E10 to test
scenario_labels <- c("minimal", "partial", "complete")  # Human-readable labels
plots_list <- list()                                  # Store infection plots
plots_hosp <- list()                                  # Store hospitalization plots

#######################################################################################################
# Fig. 3a: Infection dynamics (I1 and I2) across scenarios
#######################################################################################################
for (i in seq_along(E10_values)) {
  
  val   <- E10_values[i]
  label <- scenario_labels[i]
  
  # --- Update model parameter ---
  params_case1["E10"] <- val
  
  # --- Simulate model trajectory ---
  tj_sim <- trajectory(object = SIR_matrix_mod, params = params_case1, format = "data.frame") %>%
    select(-c(".id")) %>%
    select(time, everything()) %>%
    mutate(across(where(is.numeric), ~ round(., 2)))
  
  # --- Compute total infections for each virus ---
  combined_df <- tj_sim %>%
    mutate(
      I_1 = X_IS + X_IE + X_II + X_IR,  # Virus 1 infections
      I_2 = X_SI + X_EI + X_II + X_RI   # Virus 2 infections
    ) %>%
    pivot_longer(
      cols = -time,
      names_to = "variable",
      values_to = "value"
    ) %>%
    filter(variable %in% c("I_1", "I_2"))
  
  # --- Build infection plot ---
  p <- ggplot(combined_df, aes(x = time, y = value / pop, color = variable, size = variable)) +
    geom_line() +
    scale_color_manual(
      values = c("I_1" = "navyblue", "I_2" = "red"),
      labels = c("I_1" = expression(I[1]), "I_2" = expression(I[2]))
    ) +
    scale_size_manual(
      values = c("I_1" = 3, "I_2" = 1),
      labels = c("I_1" = expression(I[1]), "I_2" = expression(I[2]))
    ) +
    scale_y_log10(
      limits = c(1e-6, NA),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(limits = c(1, 350)) +
    labs(
      x = "Time (days)",
      y = "Daily infection rate",
      color = "State",
      size = "State"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", color = "black", size = 0.3),
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5)
    )
  
  plots_list[[label]] <- p
  
  # --- Save to PDF ---
  pdf(sprintf("figure_3a_%s.pdf", label), width = 5, height = 4)
  print(p)
  dev.off()
}

# Example: Display one scenario
print(plots_list[["partial"]])

#######################################################################################################
# Fig. 3b: Hospitalization dynamics (H1, H2, H1o, H2o) across scenarios
#######################################################################################################

# Define plotting aesthetics
color_map <- c("X_H1" = "dodgerblue", "X_H1o" = "navyblue",
               "X_H2" = "coral",      "X_H2o" = "red")
linetype_map <- c("X_H1" = "solid", "X_H1o" = "dashed",
                  "X_H2" = "solid", "X_H2o" = "dashed")
label_map <- c("X_H1"  = expression(H[1]),
               "X_H1o" = expression(H[1]^(o)),
               "X_H2"  = expression(H[2]),
               "X_H2o" = expression(H[2]^(o)))

for (i in seq_along(E10_values)) {
  
  val   <- E10_values[i]
  label <- scenario_labels[i]
  
  # --- Update model parameter ---
  params_case1["E10"] <- val
  
  # --- Simulate model trajectory ---
  tj_sim <- trajectory(object = SIR_matrix_mod,
                       params = params_case1,
                       format = "data.frame") %>%
    select(-.id) %>%
    select(time, everything())
  
  # --- Reshape hospitalization data ---
  combined_df_long <- tj_sim %>%
    pivot_longer(cols = c(X_H1, X_H2, X_H1o, X_H2o),
                 names_to = "variable", values_to = "value")
  
  # --- Build hospitalization plot ---
  p_hosp <- ggplot() +
    # Solid lines = true hospitalizations
    geom_line(data = combined_df_long %>% filter(variable %in% c("X_H1","X_H2")),
              aes(x = time, y = (value * 1e5) / pop,
                  color = variable, linetype = variable),
              size = 1.2) +
    # Dashed lines = observed hospitalizations (with hitchhiker bias)
    geom_line(data = combined_df_long %>% filter(variable %in% c("X_H1o","X_H2o")),
              aes(x = time, y = (value * 1e5) / pop,
                  color = variable, linetype = variable),
              size = 2) +
    scale_color_manual(values = color_map, labels = label_map) +
    scale_linetype_manual(values = linetype_map, labels = label_map) +
    scale_y_log10(limits = c(1, NA),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(limits = c(1, 350)) +
    labs(
      x = "Time (days)",
      y = "Daily hospitalization rate × 10^5"
    ) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      legend.text  = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 16),
      axis.text.x  = element_text(size = 16),
      axis.text.y  = element_text(size = 16),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      legend.position = c(0.98, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", color = "black", size = 0.5)
    )
  
  plots_hosp[[label]] <- p_hosp
  
  # --- Save to PDF ---
  pdf(sprintf("figure_3b_%s.pdf", label), width = 5, height = 4)
  print(p_hosp)
  dev.off()
}

# Example: Display one scenario
# print(plots_hosp[["minimal"]])
