# kpower — Software Brief

## Purpose

`kpower` assesses the statistical power of an alignment to discriminate
among phylogenetic mixture models with different numbers of categories (K),
using a parametric bootstrap framework built around IQ-TREE.

Given the best-fitting K-category model estimated from an empirical
alignment, `kpower` simulates B replicate alignments under that model and
asks: do data of this size and composition reliably identify K as the
optimal number of categories? This is a direct assessment of the data's
statistical power to discriminate among K values, using information criteria
(AIC, BIC, AICc) as the discriminating criterion.

The package does not merely select the best K — it quantifies whether the
empirical alignment carries enough signal for that selection to be
meaningful.

---

## Statistical Framework

### Step 1 — Fit all K values to empirical data

For a user-specified maximum K, fit models +R1, +R2, ..., +RK to the
empirical alignment using IQ-TREE (`--fast`). Record log-likelihood, AIC,
BIC, and AICc for each K. +R1 corresponds to a single-rate model (no
FreeRate specification); +R2 upward are FreeRate models with increasing
numbers of rate categories.

### Step 2 — Identify K_best

Select K_best from the empirical runs using the user-chosen information
criterion (default: BIC). This is the reference model for simulation.

### Step 3 — Parametric bootstrap simulation

Simulate B alignments (default: 1000) from the K_best model using IQ-TREE's
AliSim, with the fitted tree and all parameter estimates from the K_best
run. Use `--site-rate SAMPLING` to correctly propagate the rate category
structure. All simulated alignments have the same length and number of taxa
as the empirical data.

### Step 4 — Refit all K values on each simulated alignment

For each of the B simulated alignments, refit +R1, +R2, ..., +RK using
IQ-TREE (`--fast`). Record AIC/BIC/AICc/lnL for each K on each replicate.

### Step 5 — Power assessment and output

Two complementary summaries:

**A. Proportion selecting K_best** — across the B simulations, what fraction
of replicates identify K_best as optimal under the chosen IC? This is the
empirical power of the data to recover the true K at this sample size.

**B. IC profile plot** — a figure with K on the x-axis and IC value on the
y-axis showing:
- Each of the B simulated replicates as a thin, semi-transparent line
  (the expected IC profile under the true model)
- The empirical alignment's IC profile as a single thicker line
- The power estimate (proportion of simulations selecting K_best) displayed
  in the figure legend

This plot reveals at a glance whether the empirical IC profile falls within
the range expected under the true model, and whether the minimum at K_best
is a sharp, consistent feature of the simulations or a weak, noisy signal.

---

## Supported Model Families

The framework is general across IQ-TREE mixture model families. The initial
implementation targets `+R{K}`, with the others as natural extensions sharing
the same IQ-TREE wrapper and simulation infrastructure:

| Family | K controls | IQ-TREE syntax | Priority |
|---|---|---|---|
| FreeRate | Rate categories | `+R{K}` | Initial release |
| Heterotachy (GHOST) | Rate + branch-length classes | `+H{K}` | Next |
| Empirical profile mixtures | Frequency profiles | `+C10`, `+C20`, ... | Extension |
| Tree mixtures (MAST) | Tree topologies | `-m "MAST+..."` | Extension |

---

## IQ-TREE Integration

### Empirical fits
```
iqtree2 -s alignment.fasta -m GTR+R{K} --fast --prefix out_K{K} -T AUTO
```

### AliSim simulation from K_best
```
iqtree2 --alisim sim_prefix -m GTR+R{K_best} -t K_best.treefile
        --length {n_sites} --num-alignments {B}
        --site-rate SAMPLING --seed {seed} -T AUTO
```

### Refit on simulated alignments
```
iqtree2 -s sim_prefix_{b}.fasta -m GTR+R{K} --fast
        --prefix sim{b}_K{K} -T AUTO
```

Key flags:
- `--fast`: heuristic tree search; essential for the B x K total runs
- `--prefix`: unique per-run prefix to avoid output file collisions
- `--alisim`: activates AliSim simulation mode
- `--site-rate SAMPLING`: samples per-site rate categories from fitted weights
- `--redo`: overwrites existing checkpoints in automated runs
- Parse `.iqtree` report files for `Log-likelihood`, `AIC`, `BIC`,
  `Number of free parameters`

---

## R Package Architecture

```
kpower/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── R/
│   ├── kpower.R             # Main entry point: kpower()
│   ├── iqtree_runner.R      # Fit one model via IQ-TREE; return parsed results
│   ├── alisim_runner.R      # Simulate B alignments via AliSim
│   ├── parse_iqtree.R       # Parse .iqtree report files (lnL, AIC, BIC, df)
│   ├── power_assessment.R   # Compute power; compare empirical vs simulated IC
│   ├── plot_kpower.R        # IC profile figure (ggplot2)
│   └── utils.R              # IQ-TREE binary detection, temp dir management
├── tests/
│   └── testthat/
│       ├── test-parse_iqtree.R
│       ├── test-iqtree_runner.R
│       └── test-power_assessment.R
└── man/
```

### Key design decisions

- `processx::run()` for all IQ-TREE calls — no shell invocation, stderr
  captured separately, timeout support
- Each IQ-TREE run gets its own subdirectory via `--prefix` to prevent
  output file collisions across parallel bootstrap replicates
- Parallelism over bootstrap replicates via `parallel::mclapply()` or
  `future.apply::future_lapply()` (user-configurable)
- Model family ("+R", "+H", "+C") is a string parameter, enabling extension
  without structural changes to the core functions
- `find_iqtree()` searches PATH then a user option `kpower.iqtree_path`
- `ggplot2` for the IC profile figure: simulated replicates as
  `geom_line(alpha = 0.05)`, empirical data as `geom_line(linewidth = 1.2)`,
  power estimate injected into the legend via a custom label

---

## Primary Output

`kpower()` returns a list (class `kpower_result`) containing:

- `$empirical`: data frame of K, lnL, AIC, BIC, AICc for the empirical runs
- `$simulations`: data frame of replicate, K, lnL, AIC, BIC, AICc for all
  B x K simulation fits
- `$K_best`: the K selected from empirical data
- `$power`: proportion of simulations that recover K_best under the chosen IC
- `$plot`: a ggplot2 object (IC profile figure)

---

## Key References

- Le & Gascuel (2008) — C10–C60 empirical profile models
- Lartillot & Philippe (2004) — Bayesian CAT mixture model
- Goldman (1993) — parametric bootstrap LRT in phylogenetics
- Goldman, Anderson & Rodrigo (2000) — SOWH test framework
- Susko et al. (2023) — AIC/BIC unreliable for mixture model selection
- Banos, Susko & Roger (2024) — consistency of mixture model estimation
- Wong et al. (2024) — MixtureFinder (nearest precedent tool)
- Minh et al. (2020) — IQ-TREE 2 and AliSim
