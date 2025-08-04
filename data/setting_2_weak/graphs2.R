library(purrr)
library(tidyr)
library(ggplot2)
results_control_all <- get(load("outputs/all_beta_and_coverage/2000_control_reri.Rdata")[1])
results_both_all <- get(load("outputs/all_beta_and_coverage/2000_both_reri.Rdata")[1])
source("scripts/summary_table/summary_fitting_once.R")
source("scripts/kable_graph/making_kable.R")
source("scripts/kable_graph/bias_graph.R")

population_truth <- get(load("outputs/population_truth.Rdata"))

model_names <- c("Default LOM", "MSLOM(X7X8)", "MSLOM(X1278)", "MSLOM(X3478)", "MSLOM(X5678)", "MSLOM(X123478)", "MSLOM(all)")

merged_bias <- model_names %>% 
  map(~ compute_b3_bias(.x, results_both_all, population_truth))  %>% 
  reduce(full_join, by = "beta0")

bias_long <- merged_bias %>% 
  pivot_longer(
    -beta0,
    names_to   = "Model",
    values_to  = "bias"
  )

bias_long$Model[bias_long$Model == "Default LOM"] <- "SLOM"

ggplot(bias_long,
       aes(x = beta0, 
           y = bias,
           colour = Model,
           linewidth= Model,
           size = Model,
           alpha = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c("SLOM" = "black",
                                 "MSLOM(X7X8)" = "red",
                                 "MSLOM(X1278)" = "lightgrey",
                                 "MSLOM(X3478)" = "darkgrey",
                                 "MSLOM(X5678)" = "blue",
                                 "MSLOM(X123478)" = "#E69F00",
                                 "MSLOM(all)" = "#56B4E9")) +
  scale_alpha_manual(values = c("SLOM" = 1,
                                "MSLOM(X7X8)" = 1,
                                "MSLOM(X1278)" = 1,
                                "MSLOM(X3478)" = 1,
                                "MSLOM(X5678)" = 1,
                                "MSLOM(X123478)" = 1,
                                "MSLOM(all)" = 1)) +
  labs(x = expression(beta[0]), y = "Bias (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x  = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        legend.text  = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1.1, "cm"))
