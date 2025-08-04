source("scripts/summary_table/summary_fitting_once.R")
source("scripts/kable_graph/making_kable.R")
library(tibble)
library(purrr)

compute_b3_bias <- function(model_name,
                            results_list,
                            population_truth,
                            const = 2.64) {
  
  purrr::map_dfr(names(results_list), function(beta0_chr) {
    tab <- results_list[[beta0_chr]]
    
    idx <- which(tab$Model == model_name)
    if (length(idx) == 0) {
      warning(sprintf("Model '%s' not found for beta0 = %s – skipped.", 
                      model_name, beta0_chr))
      return(NULL)
    }
    
    est <- as.numeric(sub("\\s*\\(.*", "", tab$b3[idx]))
    
    ref <- if (idx %in% c(1, 3)) {
      population_truth[[beta0_chr]][["true_marginal_RERI_OR"]]
    } else {
      const
    }
    # It is the percentage of bias now
    tibble::tibble(beta0 = as.numeric(beta0_chr),
                   !!model_name := 100 * (abs(est - ref))/ref)
  })
}

# dd <- compute_b3_bias("MSLOM(all)", results_both, population_truth)
# plot(dd)
