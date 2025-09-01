#######################################################################################################
# This script run the simulations and store the data in the output folder 
#######################################################################################################

rm(list = ls())
setwd("/Users/u_kumari/Dropbox/Max_Planck/Hitchhiker")
## Load packages
library(tidyverse)
library(pomp)
library(ggplot2)
library(dplyr)
print(packageVersion("pomp")) # Should be at least 3.6
print(packageVersion("tidyverse")) # Should be 1.3.1
theme_set(theme_bw()) # Set theme for ggplot2
par(mfrow = c(1, 1), lty = 1, bty = "l") # Set the graphical parameters 
source("scripts/mod_SEIR.R")

#data_points <- seq(1000, 1, length.out = 100)
## Create a logrithmic space dataframe 
data_points <- exp(seq(log(5000), log(1), length.out = 100))

#data_points <- [1, 100, 999, 1000]
#### check the simulation 
i = 100
params_v1["E0"] <-1
#### virus 1
tj_1 <- trajectory(object = SEIR_mod1, params = params_v1, format = "data.frame") %>%
  select(-c(".id")) %>%
  select(time, everything()) %>%
  rename_with(~ paste0(., "_1"), -time) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))


# Plot state variables from the simulation
ggplot(tj_1, aes(x = time)) + 
  geom_line(aes(y = S_1, color = "S_1")) +
  geom_line(aes(y = E_1, color = "E_1")) +
  geom_line(aes(y = I_1, color = "I_1")) +
  geom_line(aes(y = R_1, color = "R_1")) +
  labs(x = "Time", y = "Proportion", title = "State Variables Over Time") +
  scale_color_manual(name = "State Variable", values = c("S_1" = "blue", "E_1" = "green", "I_1" = "red", "R_1" = "purple")) +
  theme_minimal() +
  xlim(0, 200)  # Set x-axis limits from 0 to 500 (adjust as needed)

#### virus 2
tj_2 <- trajectory(object = SEIR_mod1,params = params_v2, format = "data.frame") %>% 
  select(-c(".id")) %>% 
  select(time, everything())%>%
  rename_with(~ paste0(., "_2"), -time) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))
ggplot(tj_2, aes(x = time)) + 
  geom_line(aes(y = S_2, color = "S_2")) +
  geom_line(aes(y = E_2, color = "E_2")) +
  geom_line(aes(y = I_2, color = "I_2")) +
  geom_line(aes(y = R_2, color = "R_2")) +
  labs(x = "Time", y = "Proportion", title = "State Variables Over Time") +
  scale_color_manual(name = "State Variable", values = c("S_2" = "blue", "E_2" = "green", "I_2" = "red", "R_2" = "purple")) +
  theme_minimal() +
  xlim(0, 200)  # Set x-axis limits from 0 to 500 (adjust as needed)

# Combine data frames
combined_df <- full_join(tj_1, tj_2, by = "time")
combined_df <- combined_df %>%
  mutate(
    H_t = rowSums(select(., H_1, H_2), na.rm = TRUE))

#### false positive
combined_df <- combined_df %>%
  mutate(
    FP_1 = round((I_1/pop) * H_2, 2),
    FP_2 = round((I_2/pop) * H_1, 2),
    Total_1 = round(H_1 + FP_1, 2), 
    Total_2 = round(H_2 + FP_2, 2),
    H_1_o = round(Total_1, 2),
    H_2_o = round(Total_2, 2)
  )

#### The whole range of E0 is simulated and hospitalization quantities are calculated 
#### and the data are saved file 
rho_points <- seq(0.01, 1, length.out = 10)
for (k in 1:10) {
  params_v2["rho"] <- rho_points[k]
  for (i in 1:100) {
    params_v1["E0"] <-data_points[i]
    tj_1 <- trajectory(object = SEIR_mod1, params = params_v1, format = "data.frame") %>%
      select(-c(".id")) %>%
      select(time, everything()) %>%
      rename_with(~ paste0(., "_1"), -time)
    # Save the data as an RDS file
    saveRDS(tj_1, file = "output_2/tj_1.rds")
    # Check total pop size
    plot(tj_1$time, rowSums(select(tj_1, !c(time, contains("H")))), 
         type = "l", ylim  = c(0,1), 
         xlab  = "Time", ylab = "Total proportion")
    
    # SIR model
    tj_2 <- trajectory(object = SEIR_mod1,params = params_v2, format = "data.frame") %>% 
      select(-c(".id")) %>% 
      select(time, everything())%>%
      rename_with(~ paste0(., "_2"), -time)
    saveRDS(tj_2, file = "output_2/tj_2.rds")
    
    
    
    #### concatenating the data 
    # Combine data frames
    combined_df <- full_join(tj_1, tj_2, by = "time")
    combined_df <- combined_df %>%
      mutate(
        H_t = rowSums(select(., H_1, H_2), na.rm = TRUE))
    
    #### false positive
    combined_df <- combined_df %>%
      mutate(
        #FP_1 = round((I_1/pop) * H_2, 2),
        #FP_2 = round((I_2/pop) * H_1, 2),
        FP_1 = round( (1-params_v1['rho'])* (I_1/pop) * H_2, 2),
        FP_2 = round( (1-params_v2['rho'])* (I_2/pop) * H_1, 2),
        Total_1 = round(H_1 + FP_1, 2), 
        Total_2 = round(H_2 + FP_2, 2),
        H_1_o = round(Total_1, 2),
        H_2_o = round(Total_2, 2))
    # combined_df$H_t <- rowSums(combined_df[, c("H_1", "H_2")], na.rm = TRUE)
    # combined_df$TP_1 <- ifelse(combined_df$H_t == 0, NA, combined_df$H_1 / combined_df$H_t)
    # combined_df$TP_2 <- ifelse(combined_df$H_t == 0, NA, combined_df$H_2 / combined_df$H_t)
    # combined_df$TP_t <- rowSums(combined_df[, c("TP_1", "TP_2")], na.rm = TRUE)
    # 
    # #### false positive
    # combined_df$FP_1<- combined_df$I_1 * combined_df$H_2
    # combined_df$FP_2<- combined_df$I_2 * combined_df$H_1
    # combined_df$Total_1<- combined_df$TP_1 + combined_df$FP_1
    # combined_df$Total_2<- combined_df$TP_2 + combined_df$FP_2
    # 
    # combined_df$H_1_o<- combined_df$Total_1 * combined_df$H_t
    # combined_df$H_2_o<- combined_df$Total_2 * combined_df$H_t
    
    #saveRDS(combined_df, file = "output/tj_combined.rds")
    
    
    variable_name <- i
    rho_name <-k
    #output_file <- sprintf("output_2/tj_combined_%s.rds", variable_name)
    output_file <- sprintf("output_2/tj_combined_rho%s_eo%s.rds", rho_points[k], variable_name)
    saveRDS(combined_df, file = output_file)
  
  }
}
warnings()


