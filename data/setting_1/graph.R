library(purrr)
library(tidyr)
library(ggplot2)
results_control_all <- get(load("outputs/all_beta_and_coverage/2000_control_reri.Rdata")[1])
results_both_all <- get(load("outputs/all_beta_and_coverage/2000_both_reri.Rdata")[1])
source("scripts/summary_table/summary_fitting_once.R")
source("scripts/kable_graph/making_kable.R")
source("scripts/kable_graph/bias_graph.R")

population_truth <- get(load("outputs/population_truth.Rdata"))

model_names_common <- c("Default LOM", "CLOM(all)", "DR")
mslom_name         <- "MSLOM(corresp.)"

# ---- helper already defined: compute_b3_bias() ------------------------------

# 1) Common models (same in both runs)
bias_common <- model_names_common %>%
  map(~ compute_b3_bias(.x, results_control_all, population_truth)) %>%
  reduce(full_join, by = "beta0")

# 2) MSLOM from each run, rename
bias_mslom1 <- compute_b3_bias(mslom_name, results_control_all, population_truth) %>%
  rename(MSLOM_C = !!mslom_name)

bias_mslom2 <- compute_b3_bias(mslom_name, results_both_all,    population_truth) %>%
  rename(MSLOM_A = !!mslom_name)

# 3) Combine wide, then long if you want
bias_wide <- bias_common %>%
  full_join(bias_mslom1, by = "beta0") %>%
  full_join(bias_mslom2, by = "beta0")

bias_long <- bias_wide %>%
  pivot_longer(-beta0, names_to = "Model", values_to = "bias")

bias_long <- bias_long[bias_long$Model != "DR",]
bias_long$Model[bias_long$Model == "Default LOM"] <- "SLOM"

# ggplot(bias_long,
#        aes(x = beta0, 
#            y = bias,
#            colour = Model,
#            linewidth = Model,
#            size = Model,
#            alpha = Model)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 1.2) +
#   scale_colour_manual(values = c("SLOM" = "black",
#                                  "CLOM(all)"   = "lightblue",
#                                  "MSLOM_C"  = "blue",
#                                  "MSLOM_A" = "red")) +
#   scale_alpha_manual(values = c("SLOM" = 1,
#                                 "CLOM(all)" = 1,
#                                 "MSLOM_C" = 1,
#                                 "MSLOM_A" = 1)) +
#   labs(x = expression(beta[0]),
#        y = "Bias (%)") +
#   theme_minimal() +
#   theme(legend.position = "bottom",
#         panel.grid.minor = element_blank(),
#         axis.title.x = element_text(size = 18),
#         axis.text.x  = element_text(size = 18),
#         axis.title.y = element_text(size = 18),
#         axis.text.y  = element_text(size = 18),
#         legend.text  = element_text(size = 18),   # entry text
#         legend.title = element_text(size = 18),   # title text
#         legend.key.size = unit(1.2, "cm"))
# 
ggplot(bias_long,
       aes(x = beta0,
           y = bias,
           colour = Model,
           linewidth = Model,
           size = Model,
           alpha = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  
  ## ---- scales -----------------------------------------------------------
scale_colour_manual(
  values = c("SLOM"      = "black",
             "CLOM(all)" = "lightblue",
             "MSLOM_C"   = "blue",
             "MSLOM_A"   = "red"),
  breaks = c("SLOM", "CLOM(all)", "MSLOM_C", "MSLOM_A"),
  labels = c("SLOM",
             "CLOM(all)",
             expression(plain(MSLOM)[C]),   # keeps “MSLOM” in roman
             expression(plain(MSLOM)[A]))   # <- 𝑀𝑆𝐿𝑂𝑀 with subscript A
) +
  scale_alpha_manual(values = c(SLOM = 1, `CLOM(all)` = 1,
                                MSLOM_C = 1, MSLOM_A = 1)) +
  
  ## ---- theme ------------------------------------------------------------
labs(x = expression(beta[0]), y = "Bias (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(size = 18),
        axis.text.x   = element_text(size = 18),
        axis.title.y  = element_text(size = 18),
        axis.text.y   = element_text(size = 18),
        legend.text   = element_text(size = 18),
        legend.title  = element_text(size = 18),
        legend.key.size = unit(1.2, "cm")) +
  guides(linewidth = "none",      # keep only one combined legend
         size      = "none",
         alpha     = "none")

