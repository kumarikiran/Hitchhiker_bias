# Load libraries
rm(list = ls())
library(ggplot2)
library(dplyr)
library(ggridges)

#mod <- readRDS("/Users/kiran/Dropbox/Max_Planck/Hitchhiker/matrix_model/_saved/pomp_mod_Hon_01.rds") 
data_points <- c(10000, 0.0001, 0.000001)
i <- 3
hitch_i <- 1.0
replica = 1
suffix_save <- sprintf("Hon_%d", i)
rho_i = 0.01

# Initialize list to store all means
all_param_means_list <- list()

# Loop over both groups
for (group in c("Hon", "Hoff")) {
  for (i in 1:3) {
    suffix_save <- sprintf("%s_%d", group, i)
    
    for (replica in 1:100) {
      file <- sprintf("MCMC_results/mcmc_results_%s_rho%s_r%s.rds", suffix_save, rho_i, replica)
      #file <- sprintf("MCMC_fit/results_%s_r%s.rds", suffix_save, replica)
      if (file.exists(file)) {
        samples <- readRDS(file)
        
        param_means <- colMeans(samples)
        
        all_param_means_list[[paste0(suffix_save, "_", replica)]] <- data.frame(
          group = group,
          setting = i,
          replica = replica,
          Ri1p = param_means["Ri1p"] ## change
        )
      } else {
        warning(sprintf("File not found: %s", file))
      }
    }
  }
}

# Combine into a single dataframe
param_df <- do.call(rbind, all_param_means_list)

param_df$Ri1p <- param_df$Ri1p * 10

# Make 'setting' a factor with meaningful labels
param_df$setting <- factor(param_df$setting, levels = 1:3,
                           labels = c("minimal", " partial ", " complete "))

# Optional: Combine group and setting for y-axis
param_df$group_setting <- paste(param_df$group, param_df$setting, sep = "_")

# Remove "Hon_" and "Hoff_" prefixes from y-axis labels
param_df$group_setting_clean <- gsub("^(Hon_|Hoff_)", "", param_df$group_setting)

# Plot
plot_1<-ggplot(param_df, aes(x = Ri1p, y = group_setting_clean, fill = group)) +
  geom_density_ridges(
    aes(point_color = group, point_fill = group, point_shape = group, color = group),
    alpha = 0.1,           # semi-transparent fill
    point_alpha = 0.6,     # semi-transparent points
    jittered_points = TRUE,
    size = 1             # border line thickness
  ) +
  geom_vline(xintercept = 2.5, color = "black", linetype = "dashed", size = 1.0) +   # <- vertical line
  scale_color_hue(l = 40) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22)) +
  #scale_x_continuous(limits = c(0.08, 0.16)) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  # theme(
  #   legend.position = c(1, 1),
  #   legend.justification = c(1, 1),
  #   legend.text = element_text(size = 8),
  #   legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  # )
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "none" )+
  coord_cartesian(clip = "off", ylim = c(1.2,4.2))
print(plot_1)
pdf("aa_plot.pdf", width = 5, height = 4); 
print(plot_1)
dev.off()

