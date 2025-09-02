#######################################################################################################
# This script reads simulated trajectories from `output_data/` and produces:
#  - Fig 4a: Δ amplitude in hospitalizations vs corr(I1, I2)
#  - Fig 4b: Δ peak time in hospitalizations vs corr(I1, I2)
#  - Fig 4c: Heatmap of Δ amplitude over corr(I1, I2) × pathogenicity (rho2)
#  - Fig 4d: Heatmap of Δ peak time over corr(I1, I2) × pathogenicity (rho2)
#######################################################################################################

# --- Setup ---
rm(list = ls())                                   # Clear the R workspace
here::here()                                      # Project-rooted paths via {here}
source(here("scripts", "base_packages.R"))        # Load required packages
source(here("scripts", "mod_SEIR.R"))             # Load SEIR model definition (for `pop`, etc.)

# Grid used to label E10 (initial exposure) values
data_points <- exp(seq(log(1e4), log(1e-6), length.out = 100))

# ---------------------------------------------------------------------------------------
# Read a single rho2 slice (rho2 = 0.01) across 100 E10 values and compute summary stats
# ---------------------------------------------------------------------------------------

# Pre-allocate results frame for the rho2 = 0.01 slice
df <- data.frame(
  H2o_peak      = numeric(100),   # peak of observed H2 (with hitchhiker effect)
  H2_peak       = numeric(100),   # peak of true H2 (without hitchhiker effect)
  H2o_peakTime  = numeric(100),   # time of observed peak
  H2_peakTime   = numeric(100),   # time of true peak
  corr          = numeric(100),   # corr(X_H2, X_H2o)
  index         = numeric(100)    # overlap index of I1 and I2 (area of min)
)
# add additional columns used below
df$corr_I1I2    <- numeric(100)   # corr(I1, I2) across time
df$diff_max_I1I2 <- numeric(100)  # time difference between max(I2) and max(I1)
df$JSD          <- numeric(100)   # here used as sum(min(P, Q)) similarity proxy

# Loop over 100 E10 values for a fixed rho2 = 0.01
for (i in 1:100) {
  file_path <- sprintf("output_data/tj_combined_rho0.01_E1_%s.rds", i)
  tj_combined_1 <- readRDS(file_path)
  
  # Add total infections for each virus (summing relevant compartments)
  tj_combined <- tj_combined_1 %>%
    mutate(
      I_1 = X_IS + X_IE + X_II + X_IR,
      I_2 = X_SI + X_EI + X_II + X_RI
    )
  
  # Peaks and peak times for true vs observed hospitalizations of virus 2
  df$H2_peak[i]      <- max(tj_combined$X_H2,  na.rm = TRUE)
  df$H2_peakTime[i]  <- tj_combined$time[which.max(tj_combined$X_H2)]
  df$H2o_peak[i]     <- max(tj_combined$X_H2o, na.rm = TRUE)
  df$H2o_peakTime[i] <- tj_combined$time[which.max(tj_combined$X_H2o)]
  
  # Correlation between true and observed hospitalization series for virus 2
  df$corr[i] <- cor(tj_combined$X_H2, tj_combined$X_H2o,
                    method = "pearson", use = "complete.obs")
  
  # Overlap index of infection curves for I1 and I2 (area of minimum of normalized series)
  I_1_normalized <- tj_combined$I_1 / sum(tj_combined$I_1)
  I_2_normalized <- tj_combined$I_2 / sum(tj_combined$I_2)
  df$index[i]    <- sum(pmin(I_1_normalized, I_2_normalized))
  
  # Correlation and timing difference between I1 and I2
  df$corr_I1I2[i]     <- cor(tj_combined$I_1, tj_combined$I_2,
                             method = "pearson", use = "complete.obs")
  df$diff_max_I1I2[i] <- tj_combined$time[which.max(tj_combined$I_2)] -
    tj_combined$time[which.max(tj_combined$I_1)]
  
  # Similarity proxy between normalized H2 and H2o (here using sum(min(P, Q)))
  # Note: The KL-based JSD code is commented out; we keep a bounded overlap metric instead.
  P <- tj_combined$X_H2  / sum(tj_combined$X_H2)
  Q <- tj_combined$X_H2o
  Q <- Q[!is.na(Q)]
  Q <- Q / sum(Q)
  M <- (P + Q) / 2  # (unused when using overlap proxy)
  # KL_P_M <- sum(P * log(P / M), na.rm = TRUE)
  # KL_Q_M <- sum(Q * log(Q / M), na.rm = TRUE)
  # JSD <- 0.5 * KL_P_M + 0.5 * KL_Q_M
  JSD <- sum(pmin(P, Q))           # bounded similarity in [0,1]
  df$JSD[i] <- JSD
}

# --------------------------
# Fig 4a: Δ amplitude (peak)
# --------------------------
plot_1 <- ggplot(df) +
  geom_line(aes(x = corr_I1I2, y = (H2o_peak - H2_peak) / pop),
            linetype = "dashed", color = "blue3", size = 1.2) +
  geom_point(aes(x = corr_I1I2, y = (H2o_peak - H2_peak) / pop),
             color = "blue3", shape = 16, size = 3) +
  labs(
    x = "correlation (I1, I2)",
    y = expression(Delta ~ "amplitude (daily infection rate)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 17),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

print(plot_1)
pdf("figure_4a.pdf", width = 5, height = 4); print(plot_1); dev.off()

# -----------------------------
# Fig 4b: Δ peak time (in days)
# -----------------------------
plot_2 <- ggplot(df) +
  geom_line(aes(x = corr_I1I2, y = H2o_peakTime - H2_peakTime),
            linetype = "dashed", color = "blue3", size = 1.2) +
  geom_point(aes(x = corr_I1I2, y = H2o_peakTime - H2_peakTime),
             color = "blue3", shape = 16, size = 3) +
  labs(
    x = "correlation (I1, I2)",
    y = expression(Delta ~ "peak time (days)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

print(plot_2)
pdf("figure_4b.pdf", width = 5, height = 4); print(plot_2); dev.off()

#######################################################################################################
# Heatmaps over (corr(I1, I2), rho2)
#  - Fig 4c: Δ amplitude = (H2o_peak - H2_peak)/pop
#  - Fig 4d: Δ peak time = H2o_peakTime - H2_peakTime
#######################################################################################################

# Accumulate results over rho2 ∈ [0.01, 1] and E10 grid
df_results <- data.frame(
  rho               = numeric(0),
  corr_I1I2         = numeric(0),
  E_0               = numeric(0),
  H2_peak_diff      = numeric(0),
  H2_peak_time_diff = numeric(0),
  correlation       = numeric(0)
)

rho_points <- seq(0.01, 1, length.out = 10)

for (i in 1:10) {
  for (j in 1:100) {
    file_path <- sprintf("output_data/tj_combined_rho%s_E1_%s.rds", rho_points[i], j)
    tj_combined_1 <- readRDS(file_path)
    
    # Enrich with total infections I1, I2
    tj_combined <- tj_combined_1 %>%
      mutate(
        I_1 = X_IS + X_IE + X_II + X_IR,
        I_2 = X_SI + X_EI + X_II + X_RI
      )
    
    # Peak differences and peak time differences (observed - true)
    H2_peak         <- max(tj_combined$X_H2,  na.rm = TRUE)
    H2o_peak        <- max(tj_combined$X_H2o, na.rm = TRUE)
    H2_peak_diff    <- H2o_peak - H2_peak
    H2_peakTime     <- tj_combined$time[which.max(tj_combined$X_H2)]
    H2o_peakTime    <- tj_combined$time[which.max(tj_combined$X_H2o)]
    H2_peak_time_diff <- H2o_peakTime - H2_peakTime
    
    # Correlations (true vs observed hospitalizations; I1 vs I2)
    correlation <- cor(tj_combined$X_H2, tj_combined$X_H2o,
                       method = "pearson", use = "complete.obs")
    corr_I1I2   <- cor(tj_combined$I_1, tj_combined$I_2,
                       method = "pearson", use = "complete.obs")
    
    # Append row
    df_results <- rbind(
      df_results,
      data.frame(
        rho               = rho_points[i],
        corr_I1I2         = corr_I1I2,
        E_0               = data_points[j],
        H2_peak_diff      = H2_peak_diff,
        H2_peak_time_diff = H2_peak_time_diff,
        correlation       = correlation
      )
    )
  }
}

# --- Heatmap utilities
library(akima)    # for 2D interpolation
library(viridis)  # color scale

# Clean NAs
df_clean <- df_results[complete.cases(df_results), ]

# ----------------------------
# Fig 4c: Heatmap Δ amplitude
# ----------------------------
interp_data <- with(
  df_clean,
  interp(x = corr_I1I2, y = rho, z = H2_peak_diff / pop, duplicate = "mean")
)
interp_df <- expand.grid(x = interp_data$x, y = interp_data$y)
interp_df$z <- as.vector(interp_data$z)

plt_3 <- ggplot(interp_df, aes(x = x, y = y, z = z)) +
  geom_raster(aes(fill = z)) +
  scale_fill_viridis_c(option = "viridis") +
  geom_contour(aes(z = z), color = "black", size = 0.5) +
  labs(
    x = "correlation (I1, I2)",
    y = "pathogenicity of virus 2",
    fill = "", color = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1.0),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y  = element_text(size = 14)
  ) +
  xlim(-0.10, 1) +
  ylim(0.01, 1)

print(plt_3)
pdf("figure_4c.pdf", width = 5.2, height = 4); print(plt_3); dev.off()

# ----------------------------
# Fig 4d: Heatmap Δ peak time
# ----------------------------
df_clean <- df_results[complete.cases(df_results), ]  # re-ensure cleanliness
interp_data <- with(
  df_clean,
  interp(x = corr_I1I2, y = rho, z = H2_peak_time_diff, duplicate = "mean")
)
interp_df <- expand.grid(x = interp_data$x, y = interp_data$y)
interp_df$z <- as.vector(interp_data$z)

plt_4 <- ggplot(interp_df, aes(x = x, y = y, z = z)) +
  geom_raster(aes(fill = z)) +
  scale_fill_viridis_c(option = "viridis") +
  geom_contour(aes(z = z), color = "black", size = 0.5) +
  labs(
    x = "correlation (I1, I2)",
    y = "pathogenicity of virus 2",
    fill = "", color = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1.0),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y  = element_text(size = 14)
  ) +
  xlim(-0.10, 1) +
  ylim(0.01, 1)

print(plt_4)
pdf("figure_4d.pdf", width = 5.2, height = 4); print(plt_4); dev.off()
