# Simulate a holobiont redox dataset with spatio-temporal structure

Generates a synthetic plant-soil-microbiome dataset with latent redox
dynamics, hydrological disturbance, plant physiological responses,
radial oxygen loss, root exudation, soil redox chemistry, dissolved
organic carbon, water content, redox-buffering proxies, microbial
abundance, microbial redox functional traits, gene coverage, and
optional MetaT gene expression or upregulation features.

## Usage

``` r
simulate_redox_holobiont(
  n_plot = 4,
  n_depth = 2,
  n_plant = 5,
  n_time = 15,
  p_micro = 60,
  p_gene = 36,
  seed = 123,
  include_graph = FALSE,
  depth_labels = NULL,
  gene_mode = c("both", "abundance", "metat", "none"),
  disturbance_strength = 0.4,
  disturbance_center = NULL,
  disturbance_width = 0.1,
  seasonal_amp = 0.08,
  seasonal_phase = 0,
  stochastic_reassembly = TRUE,
  decoupling = 0.3,
  zero_inflation = 0.4,
  MNAR_strength = 0.4,
  Eh_dropout_threshold = 0.25,
  micro_mean = 8,
  micro_slope = 3,
  micro_lambda_min = 1e-08,
  micro_lambda_max = 1e+06,
  metat_noise = 0.25
)
```

## Arguments

- n_plot:

  Integer \>= 1. Number of plots.

- n_depth:

  Integer \>= 1. Number of soil depth strata.

- n_plant:

  Integer \>= 1. Number of plants per plot-depth combination.

- n_time:

  Integer \>= 2. Number of time points per plant.

- p_micro:

  Integer \>= 1. Number of microbial abundance features.

- p_gene:

  Integer \>= 1. Number of microbial functional gene features.

- seed:

  Optional integer seed. Use `NULL` for stochastic runs.

- include_graph:

  Logical. If `TRUE` and igraph is available, returns a simulated
  microbial network.

- depth_labels:

  Optional character vector of length `n_depth`. If `NULL`, depths are
  labelled `"D1"`, `"D2"`, ...

- gene_mode:

  Character. One of `"both"`, `"abundance"`, `"metat"`, or `"none"`.
  Controls whether synthetic functional gene coverage, MetaT expression,
  both, or neither are returned.

- disturbance_strength:

  Numeric in `[0, 1]`. Magnitude of the disturbance pulse.

- disturbance_center:

  Optional numeric. If `NULL`, the disturbance is centred at the median
  time point.

- disturbance_width:

  Numeric \> 0. Pulse width as a fraction of the season length.

- seasonal_amp:

  Numeric \>= 0. Amplitude of seasonal forcing.

- seasonal_phase:

  Numeric. Phase shift for seasonal forcing in radians.

- stochastic_reassembly:

  Logical. If `TRUE`, microbial counts include a stochastic reassembly
  component.

- decoupling:

  Numeric in `[0, 1]`. Cross-domain decoupling. Higher values weaken
  deterministic coupling between latent redox state and observed
  domains.

- zero_inflation:

  Numeric in `[0, 1]`. Probability that a microbial abundance entry is
  set to zero.

- MNAR_strength:

  Numeric in `[0, 1]`. Strength of missing-not-at-random Eh dropout
  under strongly reduced states.

- Eh_dropout_threshold:

  Numeric in `[0, 1]`. Latent redox threshold below which Eh dropout can
  occur.

- micro_mean:

  Numeric \> 0. Baseline mean intensity for microbial counts.

- micro_slope:

  Numeric \>= 0. Coupling strength between reduced redox state and
  microbial abundance.

- micro_lambda_min:

  Numeric \> 0. Lower bound for microbial Poisson intensity.

- micro_lambda_max:

  Numeric \> 0. Upper bound for microbial Poisson intensity.

- metat_noise:

  Numeric \>= 0. Noise level added to simulated metatranscriptomic
  expression.

## Value

A named list with:

- `id`: Experimental design data frame.

- `ROS_flux`: Plant physiological and rhizosphere boundary-flux
  variables.

- `Eh_stability`: Soil redox, hydrological, carbon, and
  electron-buffering variables.

- `micro_data`: Microbial abundance feature matrix.

- `micro_traits`: Microbial redox functional trait proxies.

- `gene_abundance`: Optional synthetic functional gene
  coverage/abundance matrix.

- `gene_expression`: Optional synthetic MetaT expression matrix.

- `gene_log2fc`: Optional synthetic MetaT log2 fold-change matrix
  relative to the pre-disturbance baseline.

- `latent_truth`: Underlying latent redox state in `[0, 1]`.

- `graph`: Optional `igraph` object.

## Details

The simulator is designed for demonstrations, benchmarking, package
examples, and method development. It represents a simplified
soil-plant-microbiome continuum in which hydrological forcing modifies
latent redox state, soil water content, electron-donor supply,
electron-acceptor availability, plant stress, radial oxygen loss,
rhizodeposition, microbial abundance, redox process potential, and
functional gene response.

The microbial functional layer includes synthetic proxies for aerobic
respiration, nitrate reduction, Fe/Mn/humic extracellular electron
transfer, sulfate reduction, methanogenesis, redox mediator production,
ROS detoxification, AMF connectivity, protistan grazing, viral lysis,
microbial memory, and redox flexibility.

The generated data are synthetic and are not intended to represent any
real ecosystem. Parameter values are chosen to generate biologically
interpretable covariance structure rather than calibrated field realism.
\#' Users are not required to measure every variable generated by the
simulator. The simulator returns a rich synthetic benchmark dataset, but
external datasets may contain fewer variables. Downstream RedoxRRI
functions require only row-aligned data frames for the domains used in
the analysis. For example,
[`rri_pipeline_st()`](https://mghotbi.github.io/RRI/reference/rri_pipeline_st.md)
can be run with a small set of plant, soil, and microbial variables, and
[`rri_recovery_metrics()`](https://mghotbi.github.io/RRI/reference/rri_recovery_metrics.md)
can be run directly on any precomputed RRI trajectory stored in
`res$row_scores`.

## Examples

``` r
sim <- simulate_redox_holobiont(
  n_plot = 2,
  n_depth = 2,
  n_plant = 2,
  n_time = 8,
  p_micro = 20,
  p_gene = 36,
  gene_mode = "both",
  seed = 1
)

names(sim$ROS_flux)
#>  [1] "SPAD"                   "FvFm"                   "PhiPSII"               
#>  [4] "NPQ"                    "ROL"                    "root_exudates"         
#>  [7] "organic_acids"          "phenolics"              "exudate_redox_activity"
#> [10] "aerenchyma"             "ROL_barrier"            "root_oxidative_stress" 
#> [13] "root_redox_buffering"   "Fe_plaque_proxy"       
names(sim$Eh_stability)
#>  [1] "Eh"                             "pH"                            
#>  [3] "water_content"                  "air_filled_porosity"           
#>  [5] "pore_connectivity"              "aqueous_connectivity"          
#>  [7] "oxygen_availability"            "DOC"                           
#>  [9] "dissolved_organic_matter_redox" "EAC"                           
#> [11] "EDC"                            "redox_buffer_capacity"         
#> [13] "Fe2.Fe3"                        "Mn2.Mn4"                       
#> [15] "NH4.NO3"                        "sulfide_risk"                  
#> [17] "methane_potential"              "nitrate_reduction_potential"   
names(sim$micro_traits)
#>  [1] "aerobic_respiration"         "denitrification"            
#>  [3] "Fe_Mn_reduction"             "EET_potential"              
#>  [5] "sulfate_reduction"           "methanogenesis"             
#>  [7] "flavin_mediator"             "phenazine_mediator"         
#>  [9] "quinone_humic_shuttle"       "microbial_ROS_detox"        
#> [11] "AMF_connectivity"            "protist_grazing"            
#> [13] "viral_lysis"                 "microbial_memory"           
#> [15] "microbial_redox_flexibility"
names(sim$gene_abundance)
#>  [1] "coxA_cov"                   "coxB_cov"                  
#>  [3] "cyoA_cov"                   "cyoB_cov"                  
#>  [5] "sodA_cov"                   "katG_cov"                  
#>  [7] "narG_cov"                   "narH_cov"                  
#>  [9] "napA_cov"                   "nirK_cov"                  
#> [11] "nirS_cov"                   "norB_cov"                  
#> [13] "nosZ_cov"                   "mtrA_cov"                  
#> [15] "mtrB_cov"                   "mtrC_cov"                  
#> [17] "omcA_cov"                   "omcS_cov"                  
#> [19] "omcZ_cov"                   "cymA_cov"                  
#> [21] "dsrA_cov"                   "dsrB_cov"                  
#> [23] "aprA_cov"                   "aprB_cov"                  
#> [25] "sat_cov"                    "hdrA_cov"                  
#> [27] "mcrA_cov"                   "mcrB_cov"                  
#> [29] "mcrG_cov"                   "NiFe_hydrogenase_cov"      
#> [31] "flavin_biosynthesis_cov"    "phenazine_biosynthesis_cov"
#> [33] "quinone_reductase_cov"      "fungal_sod1_cov"           
#> [35] "fungal_CAT_cov"             "viral_capsid_proxy_cov"    
names(sim$gene_log2fc)
#>  [1] "coxA_log2FC"                   "coxB_log2FC"                  
#>  [3] "cyoA_log2FC"                   "cyoB_log2FC"                  
#>  [5] "sodA_log2FC"                   "katG_log2FC"                  
#>  [7] "narG_log2FC"                   "narH_log2FC"                  
#>  [9] "napA_log2FC"                   "nirK_log2FC"                  
#> [11] "nirS_log2FC"                   "norB_log2FC"                  
#> [13] "nosZ_log2FC"                   "mtrA_log2FC"                  
#> [15] "mtrB_log2FC"                   "mtrC_log2FC"                  
#> [17] "omcA_log2FC"                   "omcS_log2FC"                  
#> [19] "omcZ_log2FC"                   "cymA_log2FC"                  
#> [21] "dsrA_log2FC"                   "dsrB_log2FC"                  
#> [23] "aprA_log2FC"                   "aprB_log2FC"                  
#> [25] "sat_log2FC"                    "hdrA_log2FC"                  
#> [27] "mcrA_log2FC"                   "mcrB_log2FC"                  
#> [29] "mcrG_log2FC"                   "NiFe_hydrogenase_log2FC"      
#> [31] "flavin_biosynthesis_log2FC"    "phenazine_biosynthesis_log2FC"
#> [33] "quinone_reductase_log2FC"      "fungal_sod1_log2FC"           
#> [35] "fungal_CAT_log2FC"             "viral_capsid_proxy_log2FC"    

res <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = cbind(
    sim$micro_data,
    sim$micro_traits,
    sim$gene_abundance,
    sim$gene_log2fc
  ),
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm"
)

head(res$row_scores)
#>       Physio       Soil     Micro        RRI Micro_abundance Micro_network
#> 1 0.15840610 0.08260258 0.8648679 0.10529250       0.8648679            NA
#> 2 0.08326995 0.15168922 0.9391671 0.12161084       0.9391671            NA
#> 3 0.62719819 0.44492599 0.5552689 0.63060605       0.5552689            NA
#> 4 0.56009260 0.63625168 0.5247466 0.70913975       0.5247466            NA
#> 5 0.07202312 0.12724376 0.8873417 0.08995628       0.8873417            NA
#> 6 0.08680528 0.06476088 0.9619669 0.09292494       0.9619669            NA
#>   Micro_mfa
#> 1        NA
#> 2        NA
#> 3        NA
#> 4        NA
#> 5        NA
#> 6        NA
```
