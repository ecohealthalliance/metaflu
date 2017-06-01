library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)
library(purrr)
library(gridExtra)
library(abind)
registerDoMC(cores = 20)
library(doRNG)
library(grid)
set.seed(123)
source(rprojroot::find_rstudio_root_file("R","summarizing_functions.R"))


#Setup (1D)
# threshold = c(1:21)
# varied_param = "cull_time"
#
# results_list <- vary_params(param_name = varied_param, param_values = threshold)
#
# saveRDS(results_list, rprojroot::find_rstudio_root_file("inst/vary_parameters/t_detect_under_1_day.rds"))
#
# print("I AM DONE!")

#Setup (2D)
#num_of_farms = 200 #default
threshold_1 = c(1:7)
varied_param_1 = "cull_time"

threshold_2 = seq(0.1, 0.5, 0.1)
varied_param_2 = "pi_detect"

#debugging
param1_values <- threshold_1
param1_name <- varied_param_1
param2_values <- threshold_2
param2_name <- varied_param_2

results_list <- vary_params_2d(param1_name = varied_param_1, param1_values = threshold_1, param2_name = threshold_2,
                               param2_values = varied_param_2)

saveRDS(results_list, rprojroot::find_rstudio_root_file("inst/vary_parameters/cull_time_vs_pi_firsthalf.rds"))

print("I AM DONE!")

##############################################################
#visualization for 1D case

# create_pres_df <- function(res_array, scenario){
#   mean_prop_loss <- mean(get_proportion_loss(res_array)[,2])
#   mean_proportion_farms <- mean(get_number_farms(res_array)[,2]/num_of_farms)
#   mean_duration <- mean(get_duration_array(res_array)[,2])
#   mean_fraction_exposure <- mean(get_exposure_fraction(res_array)[,2])
#   prop_failure <- proportion_failed(get_failure_array(res_array))
#   df <- data.frame(scenario, mean_prop_loss, mean_proportion_farms, mean_duration,
#                    mean_fraction_exposure, prop_failure, stringsAsFactors = FALSE)
#   return(df)
# }
#
# summarized_runs <- lapply(seq_along(results_list), function(x) create_pres_df(results_list[[x]], threshold[x]))
#
# final_df <- bind_rows(summarized_runs)
#
# loss <- ggplot(data = final_df) +
#   geom_point(aes(x = threshold, y = mean_prop_loss)) +
#   labs(title = "Proportion of Loss", x = "Culling Time (Days)", y = "Proportion of Chickens Lost") +
#   theme_minimal()
# #scale_x_log10() + scale_y_log10()
#
# farms <- ggplot(data = final_df) +
#   geom_point(aes(x = threshold, y = mean_proportion_farms)) +
#   labs(title = "Proportion of Infected Farms", x = "Culling Time (Days)", y = "Proportion of Farms Infected") +
#   theme_minimal()
# #scale_x_log10() + scale_y_log10()
#
# duration <- ggplot(data = final_df) +
#   geom_point(aes(x = threshold, y = mean_duration)) +
#   labs(title = "Duration of Epidemic ", x = "Culling Time (Days)", y = "Days") +
#   theme_minimal()
# #scale_x_log10() + scale_y_log10()
#
# exposure <- ggplot(data = final_df) +
#   geom_point(aes(x = threshold, y = mean_fraction_exposure)) +
#   labs(title = "Exposure Index", x = "Culling Time (Days)", y = "Exposure Index") +
#   theme_minimal()
# #scale_x_log10() + scale_y_log10()
#
# lay <- rbind(c(1,2),
#              c(3,4))
#
# grid.arrange(loss, farms, duration, exposure, layout_matrix = lay, top=textGrob("Cull Time Parameter Analysis", gp = gpar(fontsize = 16)))
#
# saveRDS(graphs, file = "graphs.rds")