library(testthat)

test_that("rri_recovery_metrics returns direction-aware recovery columns", {
  sim <- simulate_redox_holobiont(
    n_plot = 2,
    n_depth = 1,
    n_plant = 2,
    n_time = 12,
    p_micro = 20,
    p_gene = 36,
    gene_mode = "both",
    seed = 1
  )

  micro_block <- cbind(sim$micro_traits, sim$gene_abundance, sim$gene_log2fc)

  res <- rri_pipeline_st(
    ROS_flux = sim$ROS_flux,
    Eh_stability = sim$Eh_stability,
    micro_data = micro_block,
    id = sim$id,
    reducer = "per_domain",
    scaling = "pnorm"
  )

  rec <- rri_recovery_metrics(
    res = res,
    id = sim$id,
    time_col = "time",
    group_cols = c("plot", "depth", "plant_id"),
    perturb_start = 5,
    perturb_end = 7
  )

  required_cols <- c(
    "plot", "depth", "plant_id",
    "x0", "xmin_perturb", "xmax_perturb",
    "x_extreme", "perturb_direction",
    "xmin_recovery", "xmax_recovery", "xeq",
    "A", "A_norm", "tau_lag",
    "O", "O_norm", "I", "I_norm",
    "k", "k_r2", "k_n", "k_flag",
    "tau_r", "t_half", "H", "trajectory_class"
  )

  expect_true(all(required_cols %in% names(rec)))
  expect_equal(nrow(rec), 4)
  expect_true(any(is.finite(rec$A)))
  expect_true(all(!is.na(rec$k_flag)))
})
