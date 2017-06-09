library(dplyr)
library(ggplot2)
library(plotly)
library(metaflu)


# Code to create the graphs for final presentation, from dataframe files in same folder

#To get the plots for varying 1 intervention parameter
#for i_crit
Icrit_df <- readRDS("inst/experiments/kate_presentation/Icrit_df.rds")
means <- Icrit_df %>%
  group_by(I_crit) %>%
  summarize(mean_prop_loss = mean(total_i))

prop_lost <- ggplot(data = means) +
  geom_point(aes(x = c(1:10), y = mean_prop_loss)) +
  labs(title = "Proportion of Loss", x = "Reporting Threshold", y = "Proportion of Chickens Lost") +
  scale_y_continuous(limits = c(0.2,0.8)) +
  theme(axis.text = element_text(size = 60)) +
  theme_minimal()

plot(prop_lost)

#for cull_time
cull_df <- readRDS("inst/experiments/kate_presentation/cull_df.rds")
means <- cull_df %>%
  group_by(cull_time) %>%
  summarize(mean_prop_loss = mean(total_i))

prop_lost <- ggplot(data = means) +
  geom_point(aes(x = c(1:20), y = mean_prop_loss)) +
  labs(title = "Proportion of Loss", x = "Cull Time", y = "Proportion of Chickens Lost") +
  scale_y_continuous(limits = c(0.2,0.8)) +
  theme_minimal()

plot(prop_lost)

#for pi_report
pi_df <- readRDS("inst/experiments/kate_presentation/pi_df.rds")
means <- pi_df %>%
  group_by(pi_report) %>%
  summarize(mean_prop_loss = mean(total_i))

prop_lost <- ggplot(data = means) +
  geom_point(aes(x = seq(0.05,1,0.05), y = mean_prop_loss)) +
  labs(title = "Proportion of Loss", x = "Reporting Probability", y = "Proportion of Chickens Lost") +
  scale_y_continuous(limits = c(0.2,0.8)) +
  theme_minimal()

plot(prop_lost)


# to produce the 3D graph of cull time, farm size, and proportion lost
three_dim = readRDS("inst/experiments/kate_presentation/3Dvariation.rds")
means <- three_dim %>%
  group_by(cull_time, size) %>%
  summarize(mean_prop_loss = mean(prop_loss))

plot_ly(means, x = ~cull_time, y = ~size, z = ~mean_prop_loss)


# to produce 9-panel graph showing how cull time affects proportion lost when both omega and farm size are varied
final <- readRDS("inst/experiments/kate_presentation/vary_omega_size_df.rds")
means <- final %>%
  group_by(omega, size, cull_time) %>%
  summarize(mean_prop_loss = mean(prop_loss))

ggplot(data = means, aes(x = cull_time, y = mean_prop_loss)) +
  geom_point() +
  facet_grid(omega ~ size) +
  theme_minimal()



