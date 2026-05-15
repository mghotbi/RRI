#' @title Quantify perturbation-recovery metrics from an RRI trajectory
#'
#' @description
#' Computes direction-aware perturbation-recovery metrics from a Redox
#' Resilience Index (RRI) trajectory. The RRI is treated as the observed system
#' state, \eqn{x(t)}, and the function estimates amplitude, lag, overshoot,
#' incomplete recovery, exponential recovery rate, recovery-fit diagnostics,
#' recovery timescale, half-recovery time, optional hysteresis, and categorical
#' trajectory class.
#'
#' @param res An object returned by \code{rri_pipeline_st()}. Must contain
#'   \code{res$row_scores}.
#' @param id Data frame containing time, grouping, and optional forcing
#'   variables. Must have the same number of rows as \code{res$row_scores}.
#' @param time_col Character scalar. Name of the time column in \code{id}.
#' @param group_cols Optional character vector of grouping columns in \code{id}.
#'   Metrics are computed independently for each group. If \code{NULL}, all rows
#'   are treated as one trajectory.
#' @param state_col Character scalar. Column in \code{res$row_scores} used as
#'   the system state. Defaults to \code{"RRI"}.
#' @param perturb_start Numeric scalar. Perturbation onset time, \eqn{t_0}.
#' @param perturb_end Numeric scalar. Perturbation termination time, \eqn{t_1}.
#' @param detect_threshold Numeric scalar. Absolute deviation from baseline
#'   required to define detectable response onset.
#' @param equilibrium_window Integer scalar. Number of final recovery
#'   observations used to estimate post-recovery equilibrium,
#'   \eqn{x_{\mathrm{eq}}}.
#' @param forcing_col Optional character scalar. Name of forcing variable in
#'   \code{id}. If supplied, hysteresis is estimated at matched forcing values.
#' @param overshoot_threshold Numeric scalar. Relative threshold for classifying
#'   overshoot.
#' @param incomplete_threshold Numeric scalar. Relative threshold for classifying
#'   incomplete recovery.
#' @param hysteresis_threshold Numeric scalar. Relative threshold for classifying
#'   hysteresis.
#' @param fast_k_threshold Numeric scalar. Minimum \eqn{k} used to classify fast
#'   recovery.
#' @param recovery_r2_threshold Numeric scalar in \eqn{[0, 1]}. Minimum
#'   \eqn{R^2} used to flag the exponential recovery-rate fit as reliable.
#'
#' @details
#' The function interprets an RRI trajectory as a perturbation-response state:
#'
#' \deqn{x(t)}
#'
#' where \eqn{x(t)} is the observed redox-resilience state of the holobiont at
#' time \eqn{t}. In RedoxRRI applications, \eqn{x(t)} usually corresponds to the
#' per-sample \code{RRI} score, but any continuous state column in
#' \code{res$row_scores} can be supplied through \code{state_col}.
#'
#' The trajectory is partitioned into three temporal domains:
#'
#' \enumerate{
#'   \item baseline: \eqn{t < t_0};
#'   \item perturbation: \eqn{t_0 \le t \le t_1};
#'   \item recovery: \eqn{t > t_1}.
#' }
#'
#' Baseline state is estimated as:
#'
#' \deqn{
#' x_0 =
#' \frac{1}{n_0}
#' \sum_{t_i < t_0} x(t_i)
#' }
#'
#' Perturbation amplitude is direction-aware and defined as the largest absolute
#' displacement from baseline during perturbation:
#'
#' \deqn{
#' A =
#' \max_{t_0 \le t_i \le t_1}
#' |x(t_i) - x_0|
#' }
#'
#' The perturbation extremum is:
#'
#' \deqn{
#' x_{\mathrm{extreme}} = x(t^*),
#' \quad
#' t^* =
#' \arg\max_{t_0 \le t_i \le t_1}
#' |x(t_i) - x_0|
#' }
#'
#' The perturbation direction is:
#'
#' \deqn{
#' d =
#' \mathrm{sign}(x_{\mathrm{extreme}} - x_0)
#' }
#'
#' Detectable response lag is:
#'
#' \deqn{
#' \tau_{\mathrm{lag}} =
#' t_{\mathrm{detect}} - t_0
#' }
#'
#' where \eqn{t_{\mathrm{detect}}} is the first time satisfying:
#'
#' \deqn{
#' |x(t) - x_0| > \epsilon
#' }
#'
#' and \eqn{\epsilon} is specified by \code{detect_threshold}.
#'
#' Overshoot is direction-aware. If the perturbation decreases the state,
#' overshoot is recovery above baseline. If the perturbation increases the
#' state, overshoot is recovery below baseline:
#'
#' \deqn{
#' O =
#' \left\{
#' \begin{array}{ll}
#' \max(0, x_{\max,\mathrm{rec}} - x_0), & d < 0 \\
#' \max(0, x_0 - x_{\min,\mathrm{rec}}), & d > 0 \\
#' 0, & d = 0
#' \end{array}
#' \right.
#' }
#'
#' Post-recovery equilibrium is estimated from the final recovery observations:
#'
#' \deqn{
#' x_{\mathrm{eq}} =
#' \frac{1}{m}
#' \sum_{j=1}^{m} x_j
#' }
#'
#' where \eqn{m} is determined by \code{equilibrium_window}.
#'
#' Incomplete recovery is:
#'
#' \deqn{
#' I =
#' |x_{\mathrm{eq}} - x_0|
#' }
#'
#' Recovery kinetics are estimated from the log-linearized exponential
#' relaxation model:
#'
#' \deqn{
#' x(t) =
#' x_{\mathrm{eq}} +
#' [x(t_1) - x_{\mathrm{eq}}]
#' \exp[-k(t - t_1)]
#' }
#'
#' equivalently:
#'
#' \deqn{
#' \log |x(t) - x_{\mathrm{eq}}|
#' =
#' a - k(t - t_1)
#' }
#'
#' The characteristic recovery timescale is:
#'
#' \deqn{
#' \tau_r = \frac{1}{k}
#' }
#'
#' and the half-recovery time is:
#'
#' \deqn{
#' t_{1/2} = \frac{\log(2)}{k}
#' }
#'
#' If \code{forcing_col} is supplied, hysteresis is approximated as:
#'
#' \deqn{
#' H \approx
#' \frac{1}{n}
#' \sum_i
#' |x_{\mathrm{rec}}(F_i) - x_{\mathrm{pert}}(F_i)|
#' }
#'
#' where recovery and perturbation states are compared at matched forcing values.
#'
#' Normalized metrics are returned as:
#'
#' \deqn{
#' A_n = A / |x_0|,
#' \quad
#' O_n = O / |x_0|,
#' \quad
#' I_n = I / |x_0|
#' }
#'
#' Trajectories are classified using ordered decision rules. Incomplete recovery
#' has priority over hysteresis, hysteresis over overshoot, and overshoot over
#' simple fast or slow recovery. Default thresholds are operational and should
#' be calibrated for specific experimental systems when sufficient empirical
#' replication is available.
#'
#' @return
#' A data frame with one row per trajectory. Columns include:
#'
#' \describe{
#'   \item{x0}{Estimated pre-perturbation baseline state.}
#'   \item{xmin_perturb}{Minimum state during perturbation.}
#'   \item{xmax_perturb}{Maximum state during perturbation.}
#'   \item{x_extreme}{Perturbation state farthest from baseline.}
#'   \item{perturb_direction}{Perturbation direction: \code{"increase"},
#'   \code{"decrease"}, or \code{"neutral"}.}
#'   \item{xmin_recovery}{Minimum state during recovery.}
#'   \item{xmax_recovery}{Maximum state during recovery.}
#'   \item{xeq}{Estimated post-recovery equilibrium.}
#'   \item{A}{Direction-aware perturbation amplitude.}
#'   \item{A_norm}{Amplitude normalized by \code{abs(x0)}.}
#'   \item{tau_lag}{Detectable response lag.}
#'   \item{O}{Direction-aware overshoot.}
#'   \item{O_norm}{Overshoot normalized by \code{abs(x0)}.}
#'   \item{I}{Incomplete recovery.}
#'   \item{I_norm}{Incomplete recovery normalized by \code{abs(x0)}.}
#'   \item{k}{Estimated exponential recovery rate constant.}
#'   \item{k_r2}{Coefficient of determination for the log-linear recovery fit.}
#'   \item{k_n}{Number of recovery observations used to estimate \code{k}.}
#'   \item{k_flag}{Diagnostic flag for recovery-rate estimation.}
#'   \item{tau_r}{Characteristic recovery time, \eqn{1/k}.}
#'   \item{t_half}{Half-recovery time, \eqn{\log(2)/k}.}
#'   \item{H}{Hysteresis estimate.}
#'   \item{trajectory_class}{Classified recovery trajectory.}
#' }
#'
#' @references
#' Holling CS (1973). Resilience and stability of ecological systems. Annual
#' Review of Ecology and Systematics, 4, 1--23.
#'
#' Scheffer M, Carpenter S, Foley JA, Folke C, Walker B (2001). Catastrophic
#' shifts in ecosystems. Nature, 413, 591--596.
#'
#' Gunderson LH (2000). Ecological resilience in theory and application. Annual
#' Review of Ecology and Systematics, 31, 425--439.
#'
#' @importFrom stats aggregate approx coef fitted lm residuals
#' @importFrom utils tail
#'
#' @examples
#' sim <- simulate_redox_holobiont(
#'   n_plot = 2,
#'   n_depth = 1,
#'   n_plant = 2,
#'   n_time = 12,
#'   p_micro = 20,
#'   seed = 1
#' )
#'
#' res <- rri_pipeline_st(
#'   ROS_flux = sim$ROS_flux,
#'   Eh_stability = sim$Eh_stability,
#'   micro_data = sim$micro_data,
#'   id = sim$id,
#'   reducer = "per_domain",
#'   scaling = "pnorm"
#' )
#'
#' rec <- rri_recovery_metrics(
#'   res = res,
#'   id = sim$id,
#'   time_col = "time",
#'   group_cols = c("plot", "depth", "plant_id"),
#'   perturb_start = 5,
#'   perturb_end = 7
#' )
#'
#' rec
#' rri_metric_table()
#'
#' @export
rri_recovery_metrics <- function(
    res,
    id,
    time_col,
    group_cols = NULL,
    state_col = "RRI",
    perturb_start,
    perturb_end,
    detect_threshold = 0.05,
    equilibrium_window = 3,
    forcing_col = NULL,
    overshoot_threshold = 0.25,
    incomplete_threshold = 0.25,
    hysteresis_threshold = 0.25,
    fast_k_threshold = 0.25,
    recovery_r2_threshold = 0.2
) {
  validate_recovery_inputs(
    res = res,
    id = id,
    time_col = time_col,
    group_cols = group_cols,
    state_col = state_col,
    perturb_start = perturb_start,
    perturb_end = perturb_end,
    detect_threshold = detect_threshold,
    equilibrium_window = equilibrium_window,
    forcing_col = forcing_col,
    overshoot_threshold = overshoot_threshold,
    incomplete_threshold = incomplete_threshold,
    hysteresis_threshold = hysteresis_threshold,
    fast_k_threshold = fast_k_threshold,
    recovery_r2_threshold = recovery_r2_threshold
  )

  id <- base::as.data.frame(id)

  data <- base::cbind(
    id,
    .state = base::as.numeric(res$row_scores[[state_col]])
  )

  data[[time_col]] <- base::suppressWarnings(
    base::as.numeric(data[[time_col]])
  )

  if (base::is.null(group_cols)) {
    data$.rri_group <- "all"
    group_cols <- ".rri_group"
  }

  group_key <- base::interaction(
    data[, group_cols, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )

  out <- base::lapply(base::unique(group_key), function(group_id) {
    idx <- base::which(group_key == group_id)
    group_data <- data[idx, , drop = FALSE]
    group_data <- group_data[
      base::order(group_data[[time_col]]),
      ,
      drop = FALSE
    ]

    compute_recovery_metrics_one_group(
      data = group_data,
      time_col = time_col,
      group_cols = group_cols,
      perturb_start = perturb_start,
      perturb_end = perturb_end,
      detect_threshold = detect_threshold,
      equilibrium_window = equilibrium_window,
      forcing_col = forcing_col,
      overshoot_threshold = overshoot_threshold,
      incomplete_threshold = incomplete_threshold,
      hysteresis_threshold = hysteresis_threshold,
      fast_k_threshold = fast_k_threshold,
      recovery_r2_threshold = recovery_r2_threshold
    )
  })

  out <- base::do.call(base::rbind, out)
  base::rownames(out) <- NULL
  out
}

validate_recovery_inputs <- function(
    res,
    id,
    time_col,
    group_cols,
    state_col,
    perturb_start,
    perturb_end,
    detect_threshold,
    equilibrium_window,
    forcing_col,
    overshoot_threshold,
    incomplete_threshold,
    hysteresis_threshold,
    fast_k_threshold,
    recovery_r2_threshold
) {
  if (base::is.null(res$row_scores)) {
    base::stop("`res` must contain `row_scores`.", call. = FALSE)
  }

  if (!state_col %in% base::names(res$row_scores)) {
    base::stop("`state_col` must be a column in `res$row_scores`.",
               call. = FALSE)
  }

  if (base::missing(id) || base::is.null(id)) {
    base::stop("`id` must be supplied.", call. = FALSE)
  }

  id <- base::as.data.frame(id)

  if (base::nrow(id) != base::nrow(res$row_scores)) {
    base::stop(
      "`id` and `res$row_scores` must have the same number of rows.",
      call. = FALSE
    )
  }

  if (!base::is.character(time_col) ||
      base::length(time_col) != 1 ||
      !time_col %in% base::names(id)) {
    base::stop("`time_col` must be a single column name in `id`.",
               call. = FALSE)
  }

  if (!base::is.null(group_cols) &&
      base::any(!group_cols %in% base::names(id))) {
    base::stop("All `group_cols` must be columns in `id`.", call. = FALSE)
  }

  if (!base::is.null(forcing_col) &&
      !forcing_col %in% base::names(id)) {
    base::stop("`forcing_col` must be a column in `id`.", call. = FALSE)
  }

  assert_numeric_scalar(perturb_start, "`perturb_start`")
  assert_numeric_scalar(perturb_end, "`perturb_end`")
  assert_numeric_scalar(detect_threshold, "`detect_threshold`")
  assert_numeric_scalar(equilibrium_window, "`equilibrium_window`")
  assert_numeric_scalar(overshoot_threshold, "`overshoot_threshold`")
  assert_numeric_scalar(incomplete_threshold, "`incomplete_threshold`")
  assert_numeric_scalar(hysteresis_threshold, "`hysteresis_threshold`")
  assert_numeric_scalar(fast_k_threshold, "`fast_k_threshold`")
  assert_numeric_scalar(recovery_r2_threshold, "`recovery_r2_threshold`")

  if (perturb_end <= perturb_start) {
    base::stop("`perturb_end` must be greater than `perturb_start`.",
               call. = FALSE)
  }

  if (detect_threshold < 0) {
    base::stop("`detect_threshold` must be non-negative.", call. = FALSE)
  }

  if (equilibrium_window < 1) {
    base::stop("`equilibrium_window` must be >= 1.", call. = FALSE)
  }

  if (overshoot_threshold < 0 ||
      incomplete_threshold < 0 ||
      hysteresis_threshold < 0 ||
      fast_k_threshold < 0) {
    base::stop("Classification thresholds must be non-negative.",
               call. = FALSE)
  }

  if (recovery_r2_threshold < 0 || recovery_r2_threshold > 1) {
    base::stop("`recovery_r2_threshold` must be in [0, 1].",
               call. = FALSE)
  }

  invisible(TRUE)
}

assert_numeric_scalar <- function(x, name) {
  if (!base::is.numeric(x) ||
      base::length(x) != 1 ||
      !base::is.finite(x)) {
    base::stop(name, " must be a single finite numeric value.",
               call. = FALSE)
  }

  invisible(TRUE)
}

compute_recovery_metrics_one_group <- function(
    data,
    time_col,
    group_cols,
    perturb_start,
    perturb_end,
    detect_threshold,
    equilibrium_window,
    forcing_col,
    overshoot_threshold,
    incomplete_threshold,
    hysteresis_threshold,
    fast_k_threshold,
    recovery_r2_threshold
) {
  t <- data[[time_col]]
  x <- data$.state

  ok <- base::is.finite(t) & base::is.finite(x)
  data <- data[ok, , drop = FALSE]
  t <- t[ok]
  x <- x[ok]

  if (base::length(x) < 4) {
    return(empty_recovery_row(data, group_cols))
  }

  baseline_idx <- base::which(t < perturb_start)
  perturb_idx <- base::which(t >= perturb_start & t <= perturb_end)
  recovery_idx <- base::which(t > perturb_end)

  if (base::length(baseline_idx) < 1 ||
      base::length(perturb_idx) < 1 ||
      base::length(recovery_idx) < 2) {
    return(empty_recovery_row(data, group_cols))
  }

  x0 <- base::mean(x[baseline_idx], na.rm = TRUE)

  perturb_x <- x[perturb_idx]
  recovery_t <- t[recovery_idx]
  recovery_x <- x[recovery_idx]

  xmin_perturb <- base::min(perturb_x, na.rm = TRUE)
  xmax_perturb <- base::max(perturb_x, na.rm = TRUE)

  perturb_deviation <- base::abs(perturb_x - x0)
  extreme_idx <- base::which.max(perturb_deviation)
  x_extreme <- perturb_x[extreme_idx]

  A <- base::abs(x_extreme - x0)
  A_norm <- safe_normalize(A, x0)

  perturb_direction_numeric <- base::sign(x_extreme - x0)
  perturb_direction <- direction_label(perturb_direction_numeric)

  xmin_recovery <- base::min(recovery_x, na.rm = TRUE)
  xmax_recovery <- base::max(recovery_x, na.rm = TRUE)

  n_eq <- base::min(
    base::as.integer(equilibrium_window),
    base::length(recovery_x)
  )
  xeq <- base::mean(utils::tail(recovery_x, n_eq), na.rm = TRUE)

  detectable_idx <- base::which(
    t >= perturb_start &
      base::abs(x - x0) > detect_threshold
  )

  tau_lag <- if (base::length(detectable_idx) > 0) {
    base::min(t[detectable_idx], na.rm = TRUE) - perturb_start
  } else {
    NA_real_
  }

  O <- compute_directional_overshoot(
    x0 = x0,
    xmin_recovery = xmin_recovery,
    xmax_recovery = xmax_recovery,
    perturb_direction_numeric = perturb_direction_numeric
  )
  O_norm <- safe_normalize(O, x0)

  I <- base::abs(xeq - x0)
  I_norm <- safe_normalize(I, x0)

  k_fit <- estimate_recovery_rate_diagnostics(
    time = recovery_t,
    state = recovery_x,
    xeq = xeq,
    perturb_end = perturb_end,
    recovery_r2_threshold = recovery_r2_threshold
  )

  k <- k_fit$k
  k_r2 <- k_fit$k_r2
  k_n <- k_fit$k_n
  k_flag <- k_fit$k_flag

  tau_r <- if (base::is.finite(k) && k > 0) 1 / k else NA_real_
  t_half <- if (base::is.finite(k) && k > 0) base::log(2) / k else NA_real_

  H <- if (!base::is.null(forcing_col)) {
    estimate_recovery_hysteresis(
      data = data,
      time_col = time_col,
      state_col = ".state",
      forcing_col = forcing_col,
      perturb_start = perturb_start,
      perturb_end = perturb_end
    )
  } else {
    NA_real_
  }

  trajectory_class <- classify_recovery_trajectory(
    A = A,
    O = O,
    I = I,
    H = H,
    k = k,
    k_flag = k_flag,
    overshoot_threshold = overshoot_threshold,
    incomplete_threshold = incomplete_threshold,
    hysteresis_threshold = hysteresis_threshold,
    fast_k_threshold = fast_k_threshold
  )

  group_data <- data[1, group_cols, drop = FALSE]

  base::cbind(
    group_data,
    base::data.frame(
      x0 = x0,
      xmin_perturb = xmin_perturb,
      xmax_perturb = xmax_perturb,
      x_extreme = x_extreme,
      perturb_direction = perturb_direction,
      xmin_recovery = xmin_recovery,
      xmax_recovery = xmax_recovery,
      xeq = xeq,
      A = A,
      A_norm = A_norm,
      tau_lag = tau_lag,
      O = O,
      O_norm = O_norm,
      I = I,
      I_norm = I_norm,
      k = k,
      k_r2 = k_r2,
      k_n = k_n,
      k_flag = k_flag,
      tau_r = tau_r,
      t_half = t_half,
      H = H,
      trajectory_class = trajectory_class,
      stringsAsFactors = FALSE
    )
  )
}

empty_recovery_row <- function(data, group_cols) {
  if (base::nrow(data) < 1) {
    group_data <- base::as.data.frame(
      base::as.list(
        stats::setNames(
          base::rep(NA_character_, base::length(group_cols)),
          group_cols
        )
      ),
      stringsAsFactors = FALSE
    )
  } else {
    group_data <- data[1, group_cols, drop = FALSE]
  }

  base::cbind(
    group_data,
    base::data.frame(
      x0 = NA_real_,
      xmin_perturb = NA_real_,
      xmax_perturb = NA_real_,
      x_extreme = NA_real_,
      perturb_direction = NA_character_,
      xmin_recovery = NA_real_,
      xmax_recovery = NA_real_,
      xeq = NA_real_,
      A = NA_real_,
      A_norm = NA_real_,
      tau_lag = NA_real_,
      O = NA_real_,
      O_norm = NA_real_,
      I = NA_real_,
      I_norm = NA_real_,
      k = NA_real_,
      k_r2 = NA_real_,
      k_n = NA_integer_,
      k_flag = "insufficient_trajectory",
      tau_r = NA_real_,
      t_half = NA_real_,
      H = NA_real_,
      trajectory_class = NA_character_,
      stringsAsFactors = FALSE
    )
  )
}

direction_label <- function(direction) {
  if (!base::is.finite(direction) || direction == 0) {
    return("neutral")
  }

  if (direction > 0) {
    return("increase")
  }

  "decrease"
}

compute_directional_overshoot <- function(
    x0,
    xmin_recovery,
    xmax_recovery,
    perturb_direction_numeric
) {
  if (!base::is.finite(x0) ||
      !base::is.finite(perturb_direction_numeric) ||
      perturb_direction_numeric == 0) {
    return(0)
  }

  if (perturb_direction_numeric < 0) {
    return(base::max(0, xmax_recovery - x0, na.rm = TRUE))
  }

  base::max(0, x0 - xmin_recovery, na.rm = TRUE)
}

safe_normalize <- function(value, reference) {
  if (!base::is.finite(value) ||
      !base::is.finite(reference) ||
      base::abs(reference) == 0) {
    return(NA_real_)
  }

  base::abs(value) / base::abs(reference)
}

estimate_recovery_rate_diagnostics <- function(
    time,
    state,
    xeq,
    perturb_end,
    recovery_r2_threshold
) {
  ok <- base::is.finite(time) & base::is.finite(state)
  time <- time[ok]
  state <- state[ok]

  if (base::length(state) < 3) {
    return(base::list(
      k = NA_real_,
      k_r2 = NA_real_,
      k_n = base::length(state),
      k_flag = "insufficient_recovery_points"
    ))
  }

  if (!base::is.finite(xeq)) {
    return(base::list(
      k = NA_real_,
      k_r2 = NA_real_,
      k_n = base::length(state),
      k_flag = "nonfinite_equilibrium"
    ))
  }

  delta <- base::abs(state - xeq)

  if (base::all(delta <= 0, na.rm = TRUE)) {
    return(base::list(
      k = NA_real_,
      k_r2 = NA_real_,
      k_n = base::length(state),
      k_flag = "no_remaining_deviation"
    ))
  }

  eps <- base::max(delta, na.rm = TRUE) * 1e-6
  delta <- base::pmax(delta, eps)

  fit_data <- base::data.frame(
    log_delta = base::log(delta),
    time_from_end = time - perturb_end
  )

  fit_data <- fit_data[
    base::is.finite(fit_data$log_delta) &
      base::is.finite(fit_data$time_from_end),
    ,
    drop = FALSE
  ]

  if (base::nrow(fit_data) < 3) {
    return(base::list(
      k = NA_real_,
      k_r2 = NA_real_,
      k_n = base::nrow(fit_data),
      k_flag = "insufficient_fit_points"
    ))
  }

  fit <- stats::lm(log_delta ~ time_from_end, data = fit_data)
  slope <- stats::coef(fit)[["time_from_end"]]
  k <- -base::as.numeric(slope)

  residuals <- stats::residuals(fit)
  ss_res <- base::sum(residuals^2, na.rm = TRUE)
  ss_tot <- base::sum(
    (fit_data$log_delta - base::mean(fit_data$log_delta, na.rm = TRUE))^2,
    na.rm = TRUE
  )

  k_r2 <- if (base::is.finite(ss_tot) && ss_tot > 0) {
    1 - ss_res / ss_tot
  } else {
    NA_real_
  }

  if (!base::is.finite(k) || k <= 0) {
    return(base::list(
      k = NA_real_,
      k_r2 = k_r2,
      k_n = base::nrow(fit_data),
      k_flag = "nonpositive_recovery_rate"
    ))
  }

  k_flag <- if (base::is.finite(k_r2) && k_r2 < recovery_r2_threshold) {
    "low_fit_quality"
  } else {
    "ok"
  }

  base::list(
    k = k,
    k_r2 = k_r2,
    k_n = base::nrow(fit_data),
    k_flag = k_flag
  )
}

estimate_recovery_hysteresis <- function(
    data,
    time_col,
    state_col,
    forcing_col,
    perturb_start,
    perturb_end
) {
  perturb <- data[
    data[[time_col]] >= perturb_start &
      data[[time_col]] <= perturb_end,
    ,
    drop = FALSE
  ]

  recovery <- data[
    data[[time_col]] > perturb_end,
    ,
    drop = FALSE
  ]

  perturb <- perturb[
    base::is.finite(perturb[[forcing_col]]) &
      base::is.finite(perturb[[state_col]]),
    ,
    drop = FALSE
  ]

  recovery <- recovery[
    base::is.finite(recovery[[forcing_col]]) &
      base::is.finite(recovery[[state_col]]),
    ,
    drop = FALSE
  ]

  if (base::nrow(perturb) < 2 || base::nrow(recovery) < 2) {
    return(NA_real_)
  }

  perturb <- collapse_duplicate_forcing(
    data = perturb,
    forcing_col = forcing_col,
    state_col = state_col
  )

  recovery <- collapse_duplicate_forcing(
    data = recovery,
    forcing_col = forcing_col,
    state_col = state_col
  )

  if (base::nrow(perturb) < 2 || base::nrow(recovery) < 2) {
    return(NA_real_)
  }

  force_min <- base::max(
    base::min(perturb[[forcing_col]], na.rm = TRUE),
    base::min(recovery[[forcing_col]], na.rm = TRUE)
  )

  force_max <- base::min(
    base::max(perturb[[forcing_col]], na.rm = TRUE),
    base::max(recovery[[forcing_col]], na.rm = TRUE)
  )

  if (!base::is.finite(force_min) ||
      !base::is.finite(force_max) ||
      force_max <= force_min) {
    return(NA_real_)
  }

  forcing_grid <- base::seq(force_min, force_max, length.out = 50)

  x_perturb <- stats::approx(
    x = perturb[[forcing_col]],
    y = perturb[[state_col]],
    xout = forcing_grid,
    rule = 1
  )$y

  x_recovery <- stats::approx(
    x = recovery[[forcing_col]],
    y = recovery[[state_col]],
    xout = forcing_grid,
    rule = 1
  )$y

  base::mean(base::abs(x_recovery - x_perturb), na.rm = TRUE)
}

collapse_duplicate_forcing <- function(data, forcing_col, state_col) {
  collapsed <- stats::aggregate(
    data[[state_col]],
    by = base::list(forcing = data[[forcing_col]]),
    FUN = mean,
    na.rm = TRUE
  )

  base::names(collapsed) <- c(forcing_col, state_col)
  collapsed <- collapsed[base::order(collapsed[[forcing_col]]), , drop = FALSE]
  collapsed
}

classify_recovery_trajectory <- function(
    A,
    O,
    I,
    H,
    k,
    k_flag,
    overshoot_threshold,
    incomplete_threshold,
    hysteresis_threshold,
    fast_k_threshold
) {
  if (!base::is.finite(A) || A == 0) {
    return("unclassified")
  }

  if (base::is.finite(I) && I > incomplete_threshold * A) {
    return("incomplete_recovery")
  }

  if (base::is.finite(H) && H > hysteresis_threshold * A) {
    return("hysteresis")
  }

  if (base::is.finite(O) && O > overshoot_threshold * A) {
    return("overshoot")
  }

  if (base::identical(k_flag, "ok") &&
      base::is.finite(k) &&
      k >= fast_k_threshold) {
    return("fast_recovery")
  }

  "slow_recovery"
}

#' @title Return the perturbation-recovery metric definition table
#'
#' @description
#' Returns a publication-ready table defining the perturbation-recovery metrics
#' produced by \code{rri_recovery_metrics()}. The table is suitable for
#' R Markdown, Quarto, vignettes, and supplementary methods.
#'
#' @return
#' A data frame with metric symbols, names, equations, and interpretations.
#'
#' @examples
#' rri_metric_table()
#'
#' @export
rri_metric_table <- function() {
  base::data.frame(
    symbol = c(
      "x(t)",
      "x0",
      "xmin_perturb",
      "xmax_perturb",
      "x_extreme",
      "perturb_direction",
      "xmin_recovery",
      "xmax_recovery",
      "xeq",
      "A",
      "A_norm",
      "tau_lag",
      "O",
      "O_norm",
      "I",
      "I_norm",
      "k",
      "k_r2",
      "k_n",
      "k_flag",
      "tau_r",
      "t_half",
      "H"
    ),
    metric = c(
      "System state",
      "Pre-perturbation baseline",
      "Minimum perturbation state",
      "Maximum perturbation state",
      "Perturbation extremum",
      "Perturbation direction",
      "Minimum recovery state",
      "Maximum recovery state",
      "Post-recovery equilibrium",
      "Amplitude of change",
      "Normalized amplitude",
      "Lag time",
      "Direction-aware overshoot",
      "Normalized overshoot",
      "Incomplete recovery",
      "Normalized incomplete recovery",
      "Recovery rate constant",
      "Recovery fit R-squared",
      "Recovery fit sample size",
      "Recovery fit diagnostic flag",
      "Characteristic recovery time",
      "Half-recovery time",
      "Hysteresis"
    ),
    equation = c(
      "x(t)",
      "mean{x(t): t < t0}",
      "min{x(t): t0 <= t <= t1}",
      "max{x(t): t0 <= t <= t1}",
      "x(t*) where t* = arg max |x(t) - x0| during perturbation",
      "sign(x_extreme - x0)",
      "min{x(t): t > t1}",
      "max{x(t): t > t1}",
      "mean{x(t): final recovery window}",
      "max |x(t) - x0| during perturbation",
      "A / |x0|",
      "tdetect - t0",
      "Directional exceedance beyond baseline",
      "O / |x0|",
      "|xeq - x0|",
      "I / |x0|",
      "log|x(t)-xeq| = a - k(t-t1)",
      "1 - SSres / SStot",
      "Number of recovery observations used in k fit",
      "Diagnostic classification of k estimation",
      "1 / k",
      "log(2) / k",
      "H ≈ mean(|xrec(Fi) - xpert(Fi)|)"
    ),
    interpretation = c(
      "Observed RRI or resilience-associated system state.",
      "Estimated steady-state baseline before perturbation.",
      "Lowest observed state during perturbation.",
      "Highest observed state during perturbation.",
      "State farthest from baseline during perturbation.",
      "Direction of perturbation relative to baseline.",
      "Lowest observed state during recovery.",
      "Highest observed state during recovery.",
      "Estimated equilibrium after recovery.",
      "Magnitude of perturbation displacement from baseline.",
      "Dimensionless perturbation amplitude.",
      "Delay between perturbation onset and detectable response.",
      "Transient exceedance beyond baseline opposite the perturbation direction.",
      "Dimensionless overshoot metric.",
      "Persistent displacement from baseline after recovery.",
      "Dimensionless incomplete-recovery metric.",
      "Estimated first-order exponential recovery rate.",
      "Goodness-of-fit for exponential recovery approximation.",
      "Effective sample size used to estimate recovery kinetics.",
      "Diagnostic information describing recovery-rate estimation quality.",
      "Characteristic recovery timescale.",
      "Time required to recover half the remaining deviation.",
      "Path dependence between perturbation and recovery trajectories."
    ),
    stringsAsFactors = FALSE
  )
}
