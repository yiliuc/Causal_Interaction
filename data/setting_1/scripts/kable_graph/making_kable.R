# load("scripts/population_truth/population_truth.Rdata")
source("scripts/data_generation.R")

making_prevalence_table <- function(beta0, seed = 2025){  
  set.seed(seed)
  
  sim_dat <- simulate_data_intercept(3000, beta0)
  
  # sim_dat %>%
  #   dplyr::group_by(G, E) %>%
  #   dplyr::summarise(
  #     Y0 = sum(Y == 0),
  #     Y1 = sum(Y == 1),
  #     .groups = "drop") %>%
  #   mutate(r = Y1/(Y0 + Y1)) %>% 
  #   dplyr::rename(`Y = 0` = Y0,
  #                 `Y = 1` = Y1,
  #                 `Rate` = r) %>% 
  #   arrange(G, E) %>% 
  #   kable(booktabs = TRUE,
  #         caption = "Summary of the simulated data")
  
  caption_txt <- sprintf(
    "Summary of the simulated data for $\\beta_0 = %s$",
    beta0                          # the value passed to your function
  )
  
  sim_dat %>%
  summarise(
    Y0   = sum(Y == 0),                # number of 0’s
    Y1   = sum(Y == 1),                # number of 1’s
    rate = Y1 / (Y0 + Y1)              # proportion of 1’s
  ) %>% 
    dplyr::rename(`Y = 0` = Y0,
                  `Y = 1` = Y1,
                  `Rate` = rate) %>%
    kable(booktabs = TRUE,
          caption = caption_txt)
}




make_reri_table <- function(summary_table, beta0) {
  ## --- pull the two “true” numbers for this β₀ ----------------------------
  truth      <- population_truth[[as.character(beta0)]]
  truth_or   <- round(truth[["true_marginal_RERI_OR"]], 3)   # used for bias
  truth_rr   <- round(truth[["true_marginal_RERI_RR"]], 3)
  marginal_reri_or <- exp(0.3 + 0.4 + 0.8) - exp(0.3) - exp(0.4) + 1
  
  foot_txt <- sprintf(
    "Conditional: $RERI_{OR}=2.64$; Marginal: $RERI_{OR} = %.3f$, $RERI_{RR} = %.3f$",
    truth_or, truth_rr
  )
  
  ## --- add Bias column ----------------------------------------------------
  summary_table <- summary_table %>% 
    mutate(
      ## pull the “a” part out of "a(b)"
      est_b3 = as.numeric(sub("\\(.*", "", b3)),
      bias   = case_when(
        Model %in% c("CLOM(all)", "DR")   ~ abs(est_b3 - marginal_reri_or),
        TRUE                                           ~ abs(est_b3 - truth_or)
      ),
      bias_per   = case_when(
        Model %in% c("CLOM(all)", "DR")   ~ round((abs(est_b3 - marginal_reri_or))/marginal_reri_or * 100, 2),
        TRUE                                           ~ round((abs(est_b3 - truth_or))/truth_or * 100, 2)
      ),
      bias = round(bias, 3),
      bias_per = round(bias_per, 3)
    ) %>% 
    dplyr::select(-est_b3)                                 # drop the helper column
  
  ## --- make the kable ------------------------------------------------------
  summary_table  %>% 
    dplyr::select(-b0, -b1, -b2) %>% 
    dplyr::rename(
      # `$\\beta_0$` = b0,
      # `$\\beta_1$` = b1,
      # `$\\beta_2$` = b2,
      `$\\beta_3$` = b3,
      `$\\widehat{SE}_{\\beta_3}$` = se_b3,
      `Abs. Bias` = bias,
      `Bias(\\%)` = bias_per
    ) %>% 
    dplyr::filter(Model %in% c("Default LOM", "CLOM(all)", "MSLOM(corresp.)", "DR")) %>% 
    kableExtra::kbl(booktabs = TRUE,
                    escape   = FALSE,
                    align    = c("l", rep("c", 7)))  %>% 
    kableExtra::kable_styling(position      = "center",
                              latex_options = c("hold_position", "threeparttable"),
                              full_width    = FALSE)  %>% 
    kableExtra::footnote(
      general           = foot_txt,
      general_title     = "",
      escape            = FALSE,
      footnote_as_chunk = TRUE
    )
}
