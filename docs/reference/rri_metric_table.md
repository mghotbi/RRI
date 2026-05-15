# Return the perturbation-recovery metric definition table

Returns a publication-ready table defining the perturbation-recovery
metrics produced by
[`rri_recovery_metrics()`](https://mghotbi.github.io/RRI/reference/rri_recovery_metrics.md).
The table is suitable for R Markdown, Quarto, vignettes, and
supplementary methods.

## Usage

``` r
rri_metric_table()
```

## Value

A data frame with metric symbols, names, equations, and interpretations.

## Examples

``` r
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
