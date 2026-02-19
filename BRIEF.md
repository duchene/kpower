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
empirical alignment using IQ-TREE. Each K gets its own BioNJ tree
(`-t BIONJ --tree-fix`), which is fast to compute and avoids conflating
model fit with tree topology variation. Record log-likelihood, AIC, BIC,
and AICc for each K. +R1 corresponds to a single-rate model (no FreeRate
specification); +R2 upward are FreeRate models with increasing numbers of
rate categories.

### Step 2 — Identify K_best

Select K_best from the empirical runs using the user-chosen information
criterion (default: BIC). This is the reference model for simulation.

### Step 3 — Parametric bootstrap simulation

Simulate B alignments (default: 1000) from the K_best model using IQ-TREE's
AliSim. The plain AliSim command is extracted from the `ALISIM COMMAND`
section of the K_best `.iqtree` report (which carries the fully-specified
model string and fitted tree path), then the original empirical alignment is
appended via `-s` to reproduce the observed gap pattern. All simulated
alignments have the same length and number of taxa as the empirical data.
Falls back to constructing the AliSim command from scratch if the section is
not found.

### Step 4 — Refit all K values on each simulated alignment

For each of the B simulated alignments, refit +R1, +R2, ..., +RK using
IQ-TREE with a per-K BioNJ tree (`-t BIONJ --tree-fix`). Record
AIC/BIC/AICc/lnL for each K on each replicate. Replicates are processed
in parallel via `parallel::mclapply()` when `n_cores > 1`.

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

### Empirical fits (and simulation refits)
```
iqtree3 -s alignment.fasta -m GTR+R{K} -t BIONJ --tree-fix
        --prefix out_K{K} -T {n_cores} --redo
```
Each K gets its own BioNJ tree; branch lengths are optimised under that K's
model. This avoids the ~16-second multi-threading benchmark triggered by
`-T AUTO` and keeps runs fast without a full heuristic tree search.

### AliSim simulation from K_best
The plain AliSim command is parsed from the `ALISIM COMMAND` section of the
K_best `.iqtree` report, then modified:
```
iqtree3 --alisim sim_prefix -m "GTR{...}+R{K_best}{...}" -t K_best.treefile
        -s empirical_alignment   # gap mimicking
        --length {n_sites} --num-alignments {B}
        --seed {seed} -T {n_cores} --redo
```
The quoted model string (with fitted rate/frequency parameters) is tokenised
character-by-character to avoid shell-quoting issues when passed through
`processx`.

Key flags:
- `-t BIONJ --tree-fix`: compute BioNJ tree and fix topology for likelihood
  optimisation; much faster than heuristic search
- `--prefix`: unique per-run prefix to avoid output file collisions
- `--alisim`: activates AliSim simulation mode
- `-s` (in AliSim): reproduce gap pattern from the empirical alignment
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
- Parallelism over bootstrap replicates via `parallel::mclapply()`; the
  single `n_cores` parameter drives both R-level parallelism and IQ-TREE's
  `-T` thread count
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
- `$sim_ic`: long-format data frame of replicate, K, lnL, AIC, BIC, AICc
  for all B x K simulation fits
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
