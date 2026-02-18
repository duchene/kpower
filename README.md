# kpower

**kpower** is an R package for assessing the statistical power of a sequence
alignment to discriminate among phylogenetic mixture models with different
numbers of categories (K).

Rather than simply selecting the best K from an alignment, kpower asks a more
fundamental question: *given data generated from the best-fitting K-category
model, would alignments of this size and composition reliably recover that K?*
This is a direct measure of the data's discriminatory power, using information
criteria (AIC, BIC, AICc) as the criterion.

## How it works

1. **Fit** models +R1, +R2, ..., +RK to the empirical alignment using IQ-TREE
   (`--fast` for efficiency)
2. **Identify** K_best as the model minimising the chosen information criterion
3. **Simulate** B replicate alignments from the K_best model using IQ-TREE's
   AliSim, using the ALISIM COMMAND printed in the `.iqtree` report (with gap
   mimicking by default)
4. **Refit** all K values on each simulated alignment
5. **Report** the proportion of simulations that recover K_best (the empirical
   power estimate) and produce an IC profile figure

## Requirements

- R (>= 4.1.0)
- [IQ-TREE2](http://www.iqtree.org) (>= 2.2.0) installed and available on
  PATH, or pointed to via `options(kpower.iqtree_path = "/path/to/iqtree2")`

## Installation

```r
# Install from GitHub (once the repository is public)
# install.packages("remotes")
remotes::install_github("duchenelab/kpower")

# Or install locally from source
install.packages("path/to/kpower", repos = NULL, type = "source")
```

## Quick start

```r
library(kpower)

# Verify IQ-TREE is found
check_iqtree()

# Run the power assessment
# Fits +R1 through +R5, simulates 1000 replicates from the best model,
# refits all K on each replicate, and returns a figure and power estimate.
result <- kpower(
  alignment  = "my_alignment.fasta",
  K_max      = 5,
  ic         = "BIC",
  B          = 1000,
  mimic_gaps = TRUE      # reproduce empirical gap pattern in simulations
)

# Inspect results
print(result)
#> kpower result
#>   IC used  : BIC
#>   K_best   : 3
#>   Power    : 87.4%

# View the IC profile figure
result$plot

# Save it
ggplot2::ggsave("kpower_figure.pdf", result$plot, width = 7, height = 5)

# Access underlying data
result$empirical   # IC scores from empirical fits
result$sim_ic      # IC scores from all bootstrap fits (long format)
result$K_best      # selected K
result$power       # proportion of simulations recovering K_best
```

## Output figure

The IC profile figure shows:
- **Thin blue lines** — IC profile across K for each of the B simulated
  replicates (the expected range under the true model)
- **Thick red line** — IC profile for the empirical alignment
- **Open circle** — K_best selected from empirical data
- **Legend** — empirical power estimate (% of simulations recovering K_best)

If the empirical line sits within the cloud of simulation lines and the power
is high (e.g., > 80%), the data have sufficient signal to reliably identify
K_best. If the empirical line is an outlier or power is low, results should be
interpreted cautiously.

## Supported model families

| Family | Description | `mix_type` argument |
|---|---|---|
| FreeRate | Rate categories with freely estimated rates and weights | `"+R"` (default) |
| GHOST / Heterotachy | Rate + branch-length classes | `"+H"` |
| Empirical profiles | Amino acid frequency profiles (C10, C20, ...) | planned |
| Tree mixtures (MAST) | Mixture over tree topologies | planned |

## Advanced options

```r
# Use AIC instead of BIC
result <- kpower("alignment.fasta", K_max = 4, ic = "AIC")

# Specify IQ-TREE path explicitly
result <- kpower("alignment.fasta", K_max = 4,
                 iqtree_bin = "/usr/local/bin/iqtree2")

# Parallelise bootstrap refits across 4 cores
result <- kpower("alignment.fasta", K_max = 4, B = 500, n_cores = 4)

# Use a protein model
result <- kpower("proteins.fasta", K_max = 4,
                 base_model = "LG", mix_type = "+R")

# Disable gap mimicking
result <- kpower("alignment.fasta", K_max = 4, mimic_gaps = FALSE)
```

## Citation

If you use kpower in published work, please cite:

> Duchene et al. (in prep). kpower: assessing alignment power to discriminate
> phylogenetic mixture model categories via parametric bootstrap.

and the IQ-TREE paper:

> Minh BQ et al. (2020). IQ-TREE 2: New models and methods for phylogenetic
> inference. *Molecular Biology and Evolution*, 37(5): 1530–1534.
