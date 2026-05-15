# Plot a recovery landscape from RRI perturbation-recovery metrics

Visualises user-level perturbation-recovery metrics from
[`rri_recovery_metrics()`](https://mghotbi.github.io/RRI/reference/rri_recovery_metrics.md).
Each row is one trajectory and each column is a core resilience metric.
Values are scaled within metric columns to make amplitudes, overshoot,
incomplete recovery, and recovery times comparable.

## Usage

``` r
plot_rri_recovery_landscape(
  rec,
  group_cols = c("plot", "depth", "plant_id"),
  metrics = c("A_norm", "O_norm", "I_norm", "k", "tau_r", "t_half"),
  order_by = "I_norm",
  base_size = 12
)
```

## Arguments

- rec:

  A data frame returned by
  [`rri_recovery_metrics()`](https://mghotbi.github.io/RRI/reference/rri_recovery_metrics.md).

- group_cols:

  Character vector identifying trajectory labels.

- metrics:

  Character vector of recovery metric columns to plot.

- order_by:

  Character scalar. Metric used to order trajectories.

- base_size:

  Numeric. Base font size.

## Value

A `ggplot` object.

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

head(rec)
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
```
