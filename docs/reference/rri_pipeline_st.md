# Holobiont Redox Resilience Index (RRI) with Spatio-Temporal Dynamics

Computes a holobiont-level Redox Resilience Index (RRI) by integrating
plant physiological traits, soil redox chemistry, and microbial
resilience into a unified, directionally identifiable index. The
framework supports static (snapshot), rolling-window, and event-based
resilience modes, optional compositional geometry (clr), multiblock
reduction via MFA, and covariance-based compensation.

## Usage

``` r
rri_pipeline_st(
  ROS_flux,
  Eh_stability,
  micro_data = NULL,
  graph = NULL,
  id = NULL,
  time_col = NULL,
  group_cols = NULL,
  mode = c("snapshot", "rolling", "event"),
  window = 3,
  align = c("right", "center", "left"),
  event_col = NULL,
  baseline_label = "pre",
  recovery_labels = "recovery",
  alpha_micro = 0.5,
  method_phys = "pca",
  method_soil = "pca",
  method_micro = "pca",
  direction_phys = c("auto", "higher_is_better", "lower_is_better"),
  direction_soil = c("auto", "higher_is_better", "lower_is_better"),
  direction_micro = c("auto", "higher_is_better", "lower_is_better"),
  direction_anchor_phys = NULL,
  direction_anchor_soil = NULL,
  direction_anchor_micro = NULL,
  scale_by = NULL,
  network_agg = c("equation", "mean"),
  w1 = 0.4,
  w2 = 0.35,
  w3 = 0.25,
  add_coupling = FALSE,
  coupling_weight = 0,
  coupling_fun = c("geometric_mean", "agreement"),
  norm_method = NULL,
  reducer = c("per_domain", "mfa"),
  scaling = c("minmax_legacy", "pnorm"),
  comp_space = c("closure_legacy", "clr"),
  ref_stats = NULL,
  add_compensation = FALSE,
  compensation_weight = 0
)
```

## Arguments

- ROS_flux:

  Data frame of plant physiological variables (rows = samples).

- Eh_stability:

  Data frame of soil redox chemistry variables (rows = samples).

- micro_data:

  Optional data frame of microbial abundance or functional features.

- graph:

  Optional `igraph` object or list of `igraph` objects representing
  microbial network structure.

- id:

  Optional data frame describing experimental design (same number of
  rows as inputs).

- time_col:

  Optional character. Name of time column in `id`.

- group_cols:

  Optional character vector of grouping variables in `id`.

- mode:

  Character. One of `"snapshot"`, `"rolling"`, or `"event"`.

- window:

  Integer \>= 2. Rolling window size (for mode = "rolling").

- align:

  Character. Alignment rule for rolling window: `"right"`, `"center"`,
  or `"left"`.

- event_col:

  Optional character. Column in `id` identifying event phases.

- baseline_label:

  Character. Label identifying baseline phase.

- recovery_labels:

  Character vector identifying recovery phases.

- alpha_micro:

  Numeric between 0 and 1 controlling blending of microbial abundance
  and network components.

- method_phys:

  Character. Reduction method for plant block.

- method_soil:

  Character. Reduction method for soil block.

- method_micro:

  Character. Reduction method for microbial block.

- direction_phys:

  Character. Orientation rule for plant latent dimension.

- direction_soil:

  Character. Orientation rule for soil latent dimension.

- direction_micro:

  Character. Orientation rule for microbial latent dimension.

- direction_anchor_phys:

  Optional character. Anchor variable for plant orientation.

- direction_anchor_soil:

  Optional character. Anchor variable for soil orientation.

- direction_anchor_micro:

  Optional character. Anchor variable for microbial orientation.

- scale_by:

  Optional character vector of grouping variables used for scaling.

- network_agg:

  Character. Network aggregation method: `"equation"` or `"mean"`.

- w1:

  Numeric weight for plant domain.

- w2:

  Numeric weight for soil domain.

- w3:

  Numeric weight for microbial domain. Must sum with w1 and w2 to 1.

- add_coupling:

  Logical. If TRUE, adds cross-domain coherence term.

- coupling_weight:

  Numeric between 0 and 1 controlling weight of coupling term.

- coupling_fun:

  Character. Coupling function: `"geometric_mean"` or `"agreement"`.

- norm_method:

  Optional character. If provided, overrides block-specific methods.

- reducer:

  Character. Reduction strategy: `"per_domain"` or `"mfa"`.

- scaling:

  Character. Scaling rule: `"minmax_legacy"` or `"pnorm"`.

- comp_space:

  Character. Compositional projection method: `"closure_legacy"` or
  `"clr"`.

- ref_stats:

  Optional list of reference statistics used for scaling.

- add_compensation:

  Logical. If TRUE, includes covariance-based compensation term.

- compensation_weight:

  Numeric between 0 and 1 controlling compensation weight.

## Value

A list of class `"RRI"` containing:

- `row_scores`: Raw domain and RRI values.

- `row_scores_comp`: Compositional domain scores and RRI.

- `dyn_scores`: Dynamic resilience metrics (if applicable).

- `meta`: Metadata describing model configuration.

## Details

When `reducer = "mfa"`, blocks are integrated using FactoMineR multiple
factor analysis. If partial coordinates are unavailable, the function
safely falls back to per-domain reduction.

When `comp_space = "clr"`, domain scores are projected into Aitchison
geometry using centered log-ratio transformation and returned in simplex
form for ternary visualization.

## Examples

``` r
# ---- Simulate small holobiont dataset ----
sim <- simulate_redox_holobiont(
  n_plot = 10,
  n_depth = 10,
  n_plant = 4,
  n_time = 8,
  p_micro = 20,
  seed = 1234
)

# ---- Snapshot RRI computation ----
res <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm"
)

# Per-sample domain scores and RRI
head(res$row_scores)
#>       Physio       Soil     Micro        RRI Micro_abundance Micro_network
#> 1 0.05574732 0.07688982 0.2291179 0.08913382       0.2291179            NA
#> 2 0.19742231 0.09473518 0.1449078 0.11625515       0.1449078            NA
#> 3 0.19046019 0.14878026 0.2497579 0.14917434       0.2497579            NA
#> 4 0.01720520 0.04333184 0.1404001 0.06356902       0.1404001            NA
#> 5 0.03673169 0.01515582 0.1693542 0.06595317       0.1693542            NA
#> 6 0.15538785 0.21527039 0.1898927 0.14439287       0.1898927            NA
#>   Micro_mfa
#> 1        NA
#> 2        NA
#> 3        NA
#> 4        NA
#> 5        NA
#> 6        NA

# Compositional (ternary-ready) allocation
head(res$row_scores_comp)
#>       Physio       Soil     Micro        RRI
#> 1 0.15410240 0.21254664 0.6333510 0.08913382
#> 2 0.45169977 0.21675290 0.3315473 0.11625515
#> 3 0.32336283 0.25259876 0.4240384 0.14917434
#> 4 0.08562481 0.21564871 0.6987265 0.06356902
#> 5 0.16602519 0.06850346 0.7654713 0.06595317
#> 6 0.27720557 0.38403357 0.3387609 0.14439287

# ---- Rolling dynamic mode example ----
res_roll <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  mode = "rolling",
  time_col = "time",
  group_cols = c("plot", "depth", "plant_id"),
  window = 2
)

# Note: The first (window - 1) rows per group are NA
# due to right-aligned rolling windows.
head(res_roll$dyn_scores)
#>     P_level P_stability    S_level S_stability    M_level M_stability
#> 1        NA          NA         NA          NA         NA          NA
#> 2 0.1541924   0.7819390 0.02308851   0.9673479 0.04483425   0.9762581
#> 3 0.3927778   0.8806503 0.06959809   0.9668776 0.08220687   0.9234052
#> 4 0.7385854   0.6303039 0.46524255   0.4735966 0.56082645   0.3997245
#> 5 0.9197829   0.8865559 0.91873297   0.8850711 0.99264267   0.9895952
#> 6 0.6595614   0.7454353 0.69405961   0.5673350 0.84787026   0.7848561
#>   Physio_dyn  Soil_dyn Micro_dyn   RRI_dyn
#> 1         NA        NA        NA        NA
#> 2  0.1186893 0.1638692 0.1891791 0.1596825
#> 3  0.3981783 0.2020604 0.1763458 0.2876956
#> 4  0.4772790 0.1210674 0.1389897 0.2813489
#> 5  0.8397562 0.8385868 0.9859771 0.9194183
#> 6  0.5071981 0.3886387 0.6962290 0.5384447

# System-level mean RRI
attr(res$row_scores_comp, "RRI_index")
#> [1] 0.4829054
```
