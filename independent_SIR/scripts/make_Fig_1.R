#######################################################################################################
# This script reads the results from the output folder and plots 
#######################################################################################################
#source("/scripts/base_packages")
library(ggplot2)
library(dplyr)
library(tidyr)
theme_set(theme_bw()) # Set theme for ggplot2
par(mfrow = c(1, 1), lty = 1, bty = "l") # Set the graphical parameters 
rm(list = ls())
file_path <- sprintf("output_2/tj_combined_rho0.01_eo100.rds")
tj_combined <- readRDS(file_path)
####################################
#### Plot only the infection state 
###################################
combined_df_2 <- tj_combined %>%
  pivot_longer(
    cols = -time,  # Only 'time' is kept, all others are pivoted
    names_to = "variable", 
    values_to = "value"
  ) %>%
  filter(variable %in% c("I_1", "I_2")) %>%
  mutate(variable = recode(variable, "I_1" = "I1", "I_2" = "I2"))


plot_virus <- ggplot(combined_df_2, aes(x = time, y = value, color = variable, size = variable)) +
  geom_line() +  # Plot lines with variable-specific sizes
  scale_color_manual(values = c("I1" = "dodgerblue", "I2" = "coral")) +
  scale_size_manual(values = c("I1" = 2, "I2" = 1)) +  # Set thickness for each line
  labs(
    x = "Time (days)",
    y = "Infection state"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 200)) +
  theme(
    legend.position = "none",  # Remove legend
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )
print(plot_virus)
# Save the plot to a PDF
pdf("virus_ifection_plot.pdf", width = 5, height = 4); 
print(plot_virus)
dev.off()



###### plot both hospitalization
combined_df_2 <- tj_combined %>%
  pivot_longer(
    cols = -time,  # Only 'time' is kept, all others are pivoted
    names_to = "variable", 
    values_to = "value"
  ) %>%
  filter(variable %in% c("H_1", "H_2")) %>%
  mutate(variable = recode(variable, "H_1" = "H1", "H_2" = "H2"))
plot_virus <-ggplot(combined_df_2, aes(x = time, y = value, color = variable, size = variable)) +
geom_line() +  # Plot lines with variable-specific sizes
  scale_color_manual(values = c("H1" = "dodgerblue", "H2" = "coral")) +
  scale_size_manual(values = c("H1" = 2, "H2" = 1)) +  # Set thickness for each line
  labs(
    x = "Time (days)",
    y = "Hospitalization",
    color = "Infection"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 200)) +
  #scale_y_continuous(limits = c(0.0, 1)) +
  theme(
    legend.position = "none",  # Remove legend
    #legend.position = "right",
    #legend.box = "vertical",
    #legend.title = element_text(size = 12),
    #legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )
print(plot_virus)
pdf("plot_hospitalization.pdf", width = 5, height = 4); 
print(plot_virus)
dev.off()

############################################
library(gridExtra)  # For arranging the plots side by side

# Prepare the data
combined_df_long <- tj_combined %>%
  gather(key = "variable", value = "value", H_1, H_2, H_1_o, H_2_o)

# Plot for variables ending with _1
plot_1 <- ggplot() +
  geom_line(data = combined_df_long %>% dplyr::filter(variable %in% c("H_1_o")), 
            aes(x = time, y = value, color = variable), 
            linetype = "dashed", 
            size = 2) +
  geom_line(data = combined_df_long %>% dplyr::filter(variable %in% c("H_1")), 
            aes(x = time, y = value, color = variable), 
            size = 1) +
  # geom_line(data = combined_df_long %>% dplyr::filter(variable == "H_t"), 
  #           aes(x = time, y = value), 
  #           color = "black", 
  #           size = 0.5, 
  #           linetype = "solid") +
  scale_color_manual(values = c(
    "H_1_o" = "navyblue", 
    "H_1" = "dodgerblue"
  ),
  labels = c(
    "H_1_o" = "n1 Observed", 
    "H_1" = "H1 True"
  )
  # labels = c("n_1_observed", "H_1")
  ) +
  labs(
    x = "Time (days)",
    y = "values",
    color = "Variable"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 200)) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1.0)
  )
plot_2 <- ggplot() +
  plot_2 <- ggplot() +
  geom_line(data = combined_df_long %>% dplyr::filter(variable %in% c("H_2_o")), 
            aes(x = time, y = value, color = variable), 
            linetype = "dashed", 
            size = 2) +
  geom_line(data = combined_df_long %>% dplyr::filter(variable %in% c("H_2")), 
            aes(x = time, y = value, color = variable), 
            size = 1) +
  scale_color_manual(values = c(
    "H_2_o" = "red",
    "H_2" = "coral"
  )) +
  labs(
    x = "Time (days)",
    y = "Values"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 200)) +
  theme(
    legend.position = "none",  # Remove legend
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )





#############################################################
#####
#############################################################

plot_1 <- ggplot() +
  geom_line(data = combined_df_long %>% dplyr::filter(variable %in% c("H_1_o")), 
            aes(x = time, y = value, color = variable), 
            linetype = "dashed", 
            size = 2) +
  geom_line(data = combined_df_long %>% dplyr::filter(variable %in% c("H_1")), 
            aes(x = time, y = value, color = variable), 
            size = 1) +
  scale_color_manual(values = c(
    "H_1_o" = "navyblue", 
    "H_1" = "dodgerblue"
  )) +
  labs(
    x = "Time (days)",
    y = "Hospitalization"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 200)) +
  theme(
    legend.position = "none",  # Remove legend
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

plot_1
pdf("hospit_1.pdf", width = 5, height = 4)  # Adjust width and height as needed
grid.arrange(plot_1)
dev.off()





# Plot for variables ending with _2
#pdf("hospitalization.pdf", width = 10, height = 3.5); 

plot_2 <- ggplot() +
  geom_line(data = combined_df_long %>% dplyr::filter(variable %in% c("H_2_o")), 
            aes(x = time, y = value, color = variable), 
            linetype = "dashed", 
            size = 2) +
  geom_line(data = combined_df_long %>% dplyr::filter(variable %in% c("H_2")), 
            aes(x = time, y = value, color = variable), 
            size = 1) +
  scale_color_manual(values = c(
    "H_2_o" = "red",
    "H_2" = "coral"
  )) +
  labs(
    x = "Time (days)",
    y = "Hospitalization"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 200)) +
  theme(
    legend.position = "none",  # Remove legend
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

plot_2

pdf("hospit_2.pdf", width = 5, height = 4)  # Adjust width and height as needed
grid.arrange(plot_2)
dev.off()





####### analysis 
library(dplyr)
rm(list = ls())


data_points <- 10^seq(-3, -8, length.out = 100)
data_points

df <- data.frame(
  Column1 = rep(0, 100),
  Column2 = rep(0, 100),
  Column3 = rep(0, 100),
  Column4 = rep(0, 100),
  Column5 = rep(0, 100), 
  Column6 = rep(0, 100)
)
colnames(df) <- c("H2o_peak", "H2_peak", "H2o_peakTime", "H2_peakTime", "corr", "index")

for (i in 1:100) {
  file_path <- sprintf("output_2/tj_combined_rho0.01_eo%s.rds", i)  # Or use paste0 if preferred
  tj_combined <- readRDS(file_path)
  
  df$H2_peak[i]  <- max(tj_combined$H_2, na.rm = TRUE)
  df$H2_peakTime[i]  <- tj_combined$time[which.max(tj_combined$H_2)]
  
  df$H2o_peak[i]  <- max(tj_combined$H_2_o, na.rm = TRUE)
  df$H2o_peakTime[i]  <- tj_combined$time[which.max(tj_combined$H_2_o)]
  
  correlation <- cor(tj_combined$H_2, tj_combined$H_2_o, method = "pearson", use = "complete.obs")
  df$corr[i] <- correlation
  
  I_1_normalized <- tj_combined$I_1 / sum(tj_combined$I_1)
  I_2_normalized <- tj_combined$I_2 / sum(tj_combined$I_2)
  
  # Compute the fraction of overlap by summing the minimum of the two distributions
  overlap_area <- sum(pmin(I_1_normalized, I_2_normalized))
  df$index[i] <- overlap_area
  
  # Normalize the columns (make sure they sum to 1)
  P <- tj_combined$H_2 / sum(tj_combined$H_2)
  Q <- tj_combined$H_2_o[!is.na(tj_combined$H_2_o)]  # Remove NAs
  Q <- Q / sum(Q)  # Normalize the column
  
  # Calculate Jensen-Shannon Divergence
  M <- (P + Q) / 2
  KL_P_M <- sum(P * log(P / M), na.rm = TRUE)
  KL_Q_M <- sum(Q * log(Q / M), na.rm = TRUE)
  
  #JSD <- 0.5 * KL_P_M + 0.5 * KL_Q_M
  #df$JSD[i] <- KL_P_M
  
  JSD <- sum(pmin(P, Q))
  df$JSD[i] <- JSD

}


# Load necessary library
library(ggplot2)
library(dplyr)


# Assuming 'df' has columns H2o_peak, H2_peak, and index
plot_1 <- ggplot(df) +
  geom_line(aes(x = index, y = H2o_peak), 
            linetype = "dashed", 
            color = "red", 
            size = 1.2) +
  geom_point(aes(x = index, y = H2o_peak), 
             color = "red", 
             shape = 16, 
             size = 3) +
  geom_line(aes(x = index, y = H2_peak), 
            color = "coral", 
            size = 1.2) +
  geom_point(aes(x = index, y = H2_peak), 
             color = "coral", 
             shape = 16, 
             size = 3)+
  labs(x = "Index", y = "Peak Value", color = "Legend") +
  theme_minimal()+
  labs(
    x = "overlap (I1-I2)",
    y = "values",
    color = "Variable"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1.0)
  )

# Print the plot
print(plot_1)
pdf("v1.pdf", width = 5, height = 4); 
print(plot_1)
dev.off()








plot_1 <- ggplot(df) +
  geom_point(aes(x = index, y = JSD), 
             color = "brown", 
             shape = 16,  # You can change the shape if needed
             size = 3) +
  labs(x = "Index", y = "Peak Value", color = "Legend") +
  theme_minimal()+
  labs(
    x = "overlap (I1-I2)",
    y = "values",
    color = "Variable"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1.0)
  )

# Print the plot
print(plot_1)
pdf("v1.pdf", width = 5, height = 4); 
print(plot_1)
dev.off()




