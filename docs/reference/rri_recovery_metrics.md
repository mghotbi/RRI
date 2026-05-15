# Quantify perturbation-recovery metrics from an RRI trajectory

Computes direction-aware perturbation-recovery metrics from a Redox
Resilience Index (RRI) trajectory. The RRI is treated as the observed
system state, \\x(t)\\, and the function estimates amplitude, lag,
overshoot, incomplete recovery, exponential recovery rate, recovery-fit
diagnostics, recovery timescale, half-recovery time, optional
hysteresis, and categorical trajectory class.

## Usage

``` r
rri_recovery_metrics(
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
)
```

## Arguments

- res:

  An object returned by
  [`rri_pipeline_st()`](https://mghotbi.github.io/RRI/reference/rri_pipeline_st.md).
  Must contain `res$row_scores`.

- id:

  Data frame containing time, grouping, and optional forcing variables.
  Must have the same number of rows as `res$row_scores`.

- time_col:

  Character scalar. Name of the time column in `id`.

- group_cols:

  Optional character vector of grouping columns in `id`. Metrics are
  computed independently for each group. If `NULL`, all rows are treated
  as one trajectory.

- state_col:

  Character scalar. Column in `res$row_scores` used as the system state.
  Defaults to `"RRI"`.

- perturb_start:

  Numeric scalar. Perturbation onset time, \\t_0\\.

- perturb_end:

  Numeric scalar. Perturbation termination time, \\t_1\\.

- detect_threshold:

  Numeric scalar. Absolute deviation from baseline required to define
  detectable response onset.

- equilibrium_window:

  Integer scalar. Number of final recovery observations used to estimate
  post-recovery equilibrium, \\x\_{\mathrm{eq}}\\.

- forcing_col:

  Optional character scalar. Name of forcing variable in `id`. If
  supplied, hysteresis is estimated at matched forcing values.

- overshoot_threshold:

  Numeric scalar. Relative threshold for classifying overshoot.

- incomplete_threshold:

  Numeric scalar. Relative threshold for classifying incomplete
  recovery.

- hysteresis_threshold:

  Numeric scalar. Relative threshold for classifying hysteresis.

- fast_k_threshold:

  Numeric scalar. Minimum \\k\\ used to classify fast recovery.

- recovery_r2_threshold:

  Numeric scalar in \\\[0, 1\]\\. Minimum \\R^2\\ used to flag the
  exponential recovery-rate fit as reliable.

## Value

A data frame with one row per trajectory. Columns include:

- x0:

  Estimated pre-perturbation baseline state.

- xmin_perturb:

  Minimum state during perturbation.

- xmax_perturb:

  Maximum state during perturbation.

- x_extreme:

  Perturbation state farthest from baseline.

- perturb_direction:

  Perturbation direction: `"increase"`, `"decrease"`, or `"neutral"`.

- xmin_recovery:

  Minimum state during recovery.

- xmax_recovery:

  Maximum state during recovery.

- xeq:

  Estimated post-recovery equilibrium.

- A:

  Direction-aware perturbation amplitude.

- A_norm:

  Amplitude normalized by `abs(x0)`.

- tau_lag:

  Detectable response lag.

- O:

  Direction-aware overshoot.

- O_norm:

  Overshoot normalized by `abs(x0)`.

- I:

  Incomplete recovery.

- I_norm:

  Incomplete recovery normalized by `abs(x0)`.

- k:

  Estimated exponential recovery rate constant.

- k_r2:

  Coefficient of determination for the log-linear recovery fit.

- k_n:

  Number of recovery observations used to estimate `k`.

- k_flag:

  Diagnostic flag for recovery-rate estimation.

- tau_r:

  Characteristic recovery time, \\1/k\\.

- t_half:

  Half-recovery time, \\\log(2)/k\\.

- H:

  Hysteresis estimate.

- trajectory_class:

  Classified recovery trajectory.

## Details

The function interprets an RRI trajectory as a perturbation-response
state:

\$\$x(t)\$\$

where \\x(t)\\ is the observed redox-resilience state of the holobiont
at time \\t\\. In RedoxRRI applications, \\x(t)\\ usually corresponds to
the per-sample `RRI` score, but any continuous state column in
`res$row_scores` can be supplied through `state_col`.

The trajectory is partitioned into three temporal domains:

1.  baseline: \\t \< t_0\\;

2.  perturbation: \\t_0 \le t \le t_1\\;

3.  recovery: \\t \> t_1\\.

Baseline state is estimated as:

\$\$ x_0 = \frac{1}{n_0} \sum\_{t_i \< t_0} x(t_i) \$\$

Perturbation amplitude is direction-aware and defined as the largest
absolute displacement from baseline during perturbation:

\$\$ A = \max\_{t_0 \le t_i \le t_1} \|x(t_i) - x_0\| \$\$

The perturbation extremum is:

\$\$ x\_{\mathrm{extreme}} = x(t^\*), \quad t^\* = \arg\max\_{t_0 \le
t_i \le t_1} \|x(t_i) - x_0\| \$\$

The perturbation direction is:

\$\$ d = \mathrm{sign}(x\_{\mathrm{extreme}} - x_0) \$\$

Detectable response lag is:

\$\$ \tau\_{\mathrm{lag}} = t\_{\mathrm{detect}} - t_0 \$\$

where \\t\_{\mathrm{detect}}\\ is the first time satisfying:

\$\$ \|x(t) - x_0\| \> \epsilon \$\$

and \\\epsilon\\ is specified by `detect_threshold`.

Overshoot is direction-aware. If the perturbation decreases the state,
overshoot is recovery above baseline. If the perturbation increases the
state, overshoot is recovery below baseline:

\$\$ O = \left\\ \begin{array}{ll} \max(0, x\_{\max,\mathrm{rec}} -
x_0), & d \< 0 \\ \max(0, x_0 - x\_{\min,\mathrm{rec}}), & d \> 0 \\ 0,
& d = 0 \end{array} \right. \$\$

Post-recovery equilibrium is estimated from the final recovery
observations:

\$\$ x\_{\mathrm{eq}} = \frac{1}{m} \sum\_{j=1}^{m} x_j \$\$

where \\m\\ is determined by `equilibrium_window`.

Incomplete recovery is:

\$\$ I = \|x\_{\mathrm{eq}} - x_0\| \$\$

Recovery kinetics are estimated from the log-linearized exponential
relaxation model:

\$\$ x(t) = x\_{\mathrm{eq}} + \[x(t_1) - x\_{\mathrm{eq}}\]
\exp\[-k(t - t_1)\] \$\$

equivalently:

\$\$ \log \|x(t) - x\_{\mathrm{eq}}\| = a - k(t - t_1) \$\$

The characteristic recovery timescale is:

\$\$ \tau_r = \frac{1}{k} \$\$

and the half-recovery time is:

\$\$ t\_{1/2} = \frac{\log(2)}{k} \$\$

If `forcing_col` is supplied, hysteresis is approximated as:

\$\$ H \approx \frac{1}{n} \sum_i \|x\_{\mathrm{rec}}(F_i) -
x\_{\mathrm{pert}}(F_i)\| \$\$

where recovery and perturbation states are compared at matched forcing
values.

Normalized metrics are returned as:

\$\$ A_n = A / \|x_0\|, \quad O_n = O / \|x_0\|, \quad I_n = I / \|x_0\|
\$\$

Trajectories are classified using ordered decision rules. Incomplete
recovery has priority over hysteresis, hysteresis over overshoot, and
overshoot over simple fast or slow recovery. Default thresholds are
operational and should be calibrated for specific experimental systems
when sufficient empirical replication is available.

## References

Holling CS (1973). Resilience and stability of ecological systems.
Annual Review of Ecology and Systematics, 4, 1–23.

Scheffer M, Carpenter S, Foley JA, Folke C, Walker B (2001).
Catastrophic shifts in ecosystems. Nature, 413, 591–596.

Gunderson LH (2000). Ecological resilience in theory and application.
Annual Review of Ecology and Systematics, 31, 425–439.

## Examples

``` r
sim <- simulate_redox_holobiont(
  n_plot = 2,
  n_depth = 1,
  n_plant = 2,
  n_time = 12,
  p_micro = 20,
  seed = 1
)

res <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
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

rec
#>   plot depth plant_id        x0 xmin_perturb xmax_perturb x_extreme
#> 1   P1    D1   Plant1 0.2384152    0.3713365    0.9644632 0.9644632
#> 2   P2    D1   Plant1 0.1388000    0.6030127    0.9766382 0.9766382
#> 3   P1    D1   Plant2 0.2943933    0.6328227    0.9761766 0.9761766
#> 4   P2    D1   Plant2 0.1923865    0.5299607    0.9750439 0.9750439
#>   perturb_direction xmin_recovery xmax_recovery       xeq         A   A_norm
#> 1          increase     0.2744988     0.7873794 0.3169218 0.7260480 3.045309
#> 2          increase     0.2001049     0.8181152 0.2355018 0.8378382 6.036299
#> 3          increase     0.2633617     0.8740760 0.3776235 0.6817833 2.315893
#> 4          increase     0.1383482     0.7863628 0.3414029 0.7826573 4.068150
#>   tau_lag          O    O_norm          I    I_norm         k       k_r2 k_n
#> 1       0 0.00000000 0.0000000 0.07850661 0.3292852 0.9396725 0.94822998   5
#> 2       0 0.00000000 0.0000000 0.09670178 0.6966988 0.7067229 0.51590052   5
#> 3       0 0.03103156 0.1054085 0.08323021 0.2827177 0.3206229 0.34257221   5
#> 4       0 0.05403836 0.2808843 0.14901633 0.7745673 0.1121424 0.07793873   5
#>            k_flag    tau_r    t_half  H trajectory_class
#> 1              ok 1.064201 0.7376476 NA    fast_recovery
#> 2              ok 1.414982 0.9807907 NA    fast_recovery
#> 3              ok 3.118928 2.1618765 NA    fast_recovery
#> 4 low_fit_quality 8.917231 6.1809532 NA    slow_recovery
rri_metric_table()
#>               symbol                         metric
#> 1               x(t)                   System state
#> 2                 x0      Pre-perturbation baseline
#> 3       xmin_perturb     Minimum perturbation state
#> 4       xmax_perturb     Maximum perturbation state
#> 5          x_extreme          Perturbation extremum
#> 6  perturb_direction         Perturbation direction
#> 7      xmin_recovery         Minimum recovery state
#> 8      xmax_recovery         Maximum recovery state
#> 9                xeq      Post-recovery equilibrium
#> 10                 A            Amplitude of change
#> 11            A_norm           Normalized amplitude
#> 12           tau_lag                       Lag time
#> 13                 O      Direction-aware overshoot
#> 14            O_norm           Normalized overshoot
#> 15                 I            Incomplete recovery
#> 16            I_norm Normalized incomplete recovery
#> 17                 k         Recovery rate constant
#> 18              k_r2         Recovery fit R-squared
#> 19               k_n       Recovery fit sample size
#> 20            k_flag   Recovery fit diagnostic flag
#> 21             tau_r   Characteristic recovery time
#> 22            t_half             Half-recovery time
#> 23                 H                     Hysteresis
#>                                                    equation
#> 1                                                      x(t)
#> 2                                        mean{x(t): t < t0}
#> 3                                  min{x(t): t0 <= t <= t1}
#> 4                                  max{x(t): t0 <= t <= t1}
#> 5  x(t*) where t* = arg max |x(t) - x0| during perturbation
#> 6                                      sign(x_extreme - x0)
#> 7                                         min{x(t): t > t1}
#> 8                                         max{x(t): t > t1}
#> 9                         mean{x(t): final recovery window}
#> 10                      max |x(t) - x0| during perturbation
#> 11                                                 A / |x0|
#> 12                                             tdetect - t0
#> 13                   Directional exceedance beyond baseline
#> 14                                                 O / |x0|
#> 15                                               |xeq - x0|
#> 16                                                 I / |x0|
#> 17                              log|x(t)-xeq| = a - k(t-t1)
#> 18                                        1 - SSres / SStot
#> 19            Number of recovery observations used in k fit
#> 20                Diagnostic classification of k estimation
#> 21                                                    1 / k
#> 22                                               log(2) / k
#> 23                         H ≈ mean(|xrec(Fi) - xpert(Fi)|)
#>                                                               interpretation
#> 1                        Observed RRI or resilience-associated system state.
#> 2                       Estimated steady-state baseline before perturbation.
#> 3                                 Lowest observed state during perturbation.
#> 4                                Highest observed state during perturbation.
#> 5                          State farthest from baseline during perturbation.
#> 6                            Direction of perturbation relative to baseline.
#> 7                                     Lowest observed state during recovery.
#> 8                                    Highest observed state during recovery.
#> 9                                      Estimated equilibrium after recovery.
#> 10                     Magnitude of perturbation displacement from baseline.
#> 11                                     Dimensionless perturbation amplitude.
#> 12                 Delay between perturbation onset and detectable response.
#> 13 Transient exceedance beyond baseline opposite the perturbation direction.
#> 14                                           Dimensionless overshoot metric.
#> 15                     Persistent displacement from baseline after recovery.
#> 16                                 Dimensionless incomplete-recovery metric.
#> 17                          Estimated first-order exponential recovery rate.
#> 18                   Goodness-of-fit for exponential recovery approximation.
#> 19                 Effective sample size used to estimate recovery kinetics.
#> 20       Diagnostic information describing recovery-rate estimation quality.
#> 21                                        Characteristic recovery timescale.
#> 22                    Time required to recover half the remaining deviation.
#> 23           Path dependence between perturbation and recovery trajectories.
```
