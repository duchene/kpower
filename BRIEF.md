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

### Standard pathway (+R, +H, *H)

#### Step 1 — Fit all K values to empirical data

For a user-specified maximum K, fit mixture models to the empirical
alignment using IQ-TREE.  The model family is controlled by `mix_type`:

- `+R{K}`: FreeRate — K rate categories with freely estimated rates and
  weights.
- `+H{K}`: GHOST (linked) — K branch-length classes sharing one
  substitution model.
- `*H{K}`: GHOST (unlinked) — K branch-length classes, each with its own
  substitution model parameters.

Each K gets its own BioNJ tree (`-t BIONJ --tree-fix`), which is fast to
compute and avoids conflating model fit with tree topology variation.
Record log-likelihood, AIC, BIC, and AICc for each K.

For +H and *H models, site-class probabilities are saved via `-wspm`.

#### Step 2 — Identify K_best

Select K_best from the empirical runs using the user-chosen information
criterion (default: BIC). This is the reference model for simulation.

#### Step 3 — Parametric bootstrap simulation

Simulate B alignments (default: 1000) from the K_best model using IQ-TREE's
AliSim. The plain AliSim command is extracted from the `ALISIM COMMAND`
section of the K_best `.iqtree` report (which carries the fully-specified
model string and fitted tree path), then the original empirical alignment is
appended via `-s` to reproduce the observed gap pattern. All simulated
alignments have the same length and number of taxa as the empirical data.
Falls back to constructing the AliSim command from scratch if the section is
not found.

#### Step 4 — Refit all K values on each simulated alignment

For each of the B simulated alignments, refit all K values using IQ-TREE
with a per-K BioNJ tree. Record AIC/BIC/AICc/lnL for each K on each
replicate. Replicates are processed in parallel via `parallel::mclapply()`
when `n_cores > 1`.

#### Step 5 — Power assessment and output

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

### MAST pathway (+T)

When `mix_type = "+T"`, `kpower` assesses the power to discriminate the
number of tree-topology classes (MAST model).  Candidate trees are derived
automatically from the alignment itself.

#### Step 1 — Split alignment into windows

Divide the alignment into K_max non-overlapping windows of near-equal size.

#### Step 2 — Estimate per-window trees

Run IQ-TREE with ModelFinder (`-m MFP`) and full heuristic tree search on
each window.  This yields K_max candidate trees and, as a side product,
the best-fit model for each window.

#### Step 3 — Determine within-class rate heterogeneity

Summarise the rate heterogeneity component (e.g. `+R3`, `+G4`) across the
per-window best models.  The most common FreeRate category count is used;
falls back to `+R4` when no clear pattern emerges.

#### Step 4 — Fit MAST with K_max trees

Fit the full MAST model (`base_model+FO+rate_model+T`) to the empirical
alignment using all K_max candidate trees supplied via `-te`.  Extract
tree weights.

#### Step 5 — Rank trees and build tree sets

Sort trees by descending weight from the K_max run.  For each K in
K_min..K_max, build a tree file containing the top-K ranked trees.

#### Step 6 — Fit all K values

- K = 1: standard single-tree BioNJ fit with `base_model+FO+rate_model`
  (no +T).
- K >= 2: MAST fit with top-K trees via `-te`.

Record IC scores across all K values.

#### Step 7 — Select K_best

Choose K_best by minimising the chosen IC across the K range.

#### Step 8 — Simulate B alignments

For K_best > 1, simulate from the fitted MAST model by running AliSim once
per tree class (site count proportional to class weight), then concatenating
the per-class partial alignments into full replicate alignments.  Model
parameters (substitution rates, frequencies, rate heterogeneity) are parsed
from the MAST `.iqtree` report and supplied as a fully parameterised model
string.  For K_best = 1, standard single-tree AliSim is used with gap
mimicking.

#### Step 9 — Refit all K on each replicate

For each simulated alignment, refit K = 1..K_max using the same tree sets
as the empirical analysis (K = 1 via BioNJ, K >= 2 via MAST with the
top-K ranked trees).

#### Step 10 — Power assessment and output

Identical to the standard pathway: proportion recovering K_best, plus the
IC profile figure.

---

## Supported Model Families

| Family | K controls | `mix_type` | IQ-TREE syntax | Status |
|---|---|---|---|---|
| FreeRate (linked) | Rate categories | `"+R"` | `+R{K}` | Implemented |
| FreeRate (unlinked) | Rate categories + per-class substitution | `"*R"` | `*R{K}` | Implemented |
| GHOST (linked) | Rate + branch-length classes | `"+H"` | `+H{K}` | Implemented |
| GHOST (unlinked) | Per-class substitution models + branch lengths | `"*H"` | `*H{K}` | Implemented |
| Tree mixtures (linked) | Tree topologies | `"+T"` | `model+FO+rate+T` with `-te` | Implemented |
| Tree mixtures (unlinked) | Tree topologies + per-tree substitution | `"*T"` | `MIX{model+FO,...}+rate+T` with `-te` | Implemented |
| Empirical profile mixtures | Frequency profiles | — | `+C10`, `+C20`, ... | Planned |

### Linked vs unlinked

- **Linked** (`+R`, `+H`, `+T`): all mixture classes share the same
  substitution model parameters (e.g. GTR exchangeability rates and base
  frequencies). Classes differ only in rates, branch lengths, or tree
  topology.
- **Unlinked** (`*R`, `*H`, `*T`): each mixture class gets its own
  substitution model parameters. For `*R` and `*H`, IQ-TREE implements
  this natively via the `*` prefix. For `*T`, the package constructs
  `MIX{base+FO,...,base+FO}+rate+T` syntax, which gives each tree class
  its own GTR parameters and frequencies.

---

## IQ-TREE Integration

### Standard fits (+R / +H / *H)
```
iqtree3 -s alignment.fasta -m GTR+R{K} -t BIONJ --tree-fix
        --prefix out_K{K} -T {n_cores} --redo
```
For +H and *H, the flag `-wspm` is appended to save site-class
probabilities.

### MAST fits (+T)
```
iqtree3 -s alignment.fasta -m GTR+FO+R3+T -te candidate_trees.newick
        --prefix mast_K{K} -T {n_cores} -wspm --redo
```
Candidate trees are supplied one per line in the Newick file passed to
`-te`.  Tree weights are parsed from the `Tree weights:` line of the
`.iqtree` report.

### AliSim simulation (standard)
```
iqtree3 --alisim sim_prefix -m "GTR{...}+R{K_best}{...}" -t K_best.treefile
        -s empirical_alignment   # gap mimicking
        --length {n_sites} --num-alignments {B}
        --seed {seed} -T {n_cores} --redo
```

### AliSim simulation (MAST)
For each tree class k in the MAST model, a separate AliSim call generates
B partial alignments with site count proportional to class weight:
```
iqtree3 --alisim sim_prefix_k -m "GTR{...}+FU{...}+R{n}{...}" -t tree_k.nwk
        --length {n_k} --num-alignments {B}
        --seed {seed+k} -T {n_cores} --redo
```
Per-class replicates are concatenated site-wise to produce B full-length
alignments.

Key flags:
- `-t BIONJ --tree-fix`: compute BioNJ tree and fix topology for likelihood
  optimisation; much faster than heuristic search
- `-te`: supply fixed tree topologies for MAST (one per line)
- `--prefix`: unique per-run prefix to avoid output file collisions
- `--alisim`: activates AliSim simulation mode
- `-s` (in AliSim): reproduce gap pattern from the empirical alignment
- `-wspm`: save site-class posterior probabilities
- `--redo`: overwrites existing checkpoints in automated runs

---

## R Package Architecture

```
kpower/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── R/
│   ├── kpower.R             # Main entry point: kpower(), kpower_mast()
│   ├── survey.R             # Cross-model survey: kpower_survey()
│   ├── iqtree_runner.R      # Fit one model via IQ-TREE; return parsed results
│   ├── alisim_runner.R      # Simulate via AliSim (standard + MAST)
│   ├── mast_runner.R        # MAST workflow: windows, trees, ranking, fitting
│   ├── alignment_io.R       # Read/write FASTA & PHYLIP, subset, concatenate
│   ├── parse_iqtree.R       # Parse .iqtree reports (IC, tree weights, model params)
│   ├── power_assessment.R   # Compute power (standard + MAST + all-IC)
│   ├── plot_kpower.R        # IC profile figure (ggplot2)
│   └── utils.R              # IQ-TREE binary detection, model string builders
├── tests/
│   └── testthat/
└── man/
```

### Key design decisions

- `processx::run()` for all IQ-TREE calls — no shell invocation, stderr
  captured separately, timeout support
- Each IQ-TREE run gets its own subdirectory via `--prefix` to prevent
  output file collisions across parallel bootstrap replicates
- Parallelism is controlled by two independent parameters: `n_cores`
  drives R-level parallel bootstrap refits via `parallel::mclapply()`,
  while `threads` controls IQ-TREE's `-T` thread count per run. For
  example, `n_cores = 4, threads = 2` runs 4 parallel refits each using
  2 IQ-TREE threads (8 CPUs total)
- Model family (`"+R"`, `"+H"`, `"*H"`, `"+T"`) is a string parameter;
  `+R`/`+H`/`*H` share the standard pipeline while `+T` dispatches to the
  dedicated MAST pathway
- MAST candidate trees are derived from alignment windows (K_max equal
  partitions); trees are ranked by weight from a K_max run, then
  incrementally added for K = 1, 2, ..., K_max
- MAST simulation decomposes into per-class AliSim calls (proportional
  site counts) followed by concatenation, since AliSim does not natively
  simulate from MAST models
- Within-class rate heterogeneity for MAST is auto-determined from
  per-window ModelFinder selections (modal +R or +G category count)
- `find_iqtree()` searches PATH then a user option `kpower.iqtree_path`
- `ggplot2` for the IC profile figure: simulated replicates as
  `geom_line(alpha = 0.15)`, empirical data as `geom_line(linewidth = 1.2)`,
  power estimate injected into the legend via a custom label

---

## Cross-Model Survey

`kpower_survey()` runs a two-phase pipeline across all requested mixture
families to determine which model family and K the alignment best supports.

### Phase 1 — Empirical fitting

For each mixture family in `mix_types`, fit K = K_min..K_max to the
empirical alignment and record IC profiles.  Identify K_best per family.

### Phase 2 — Parametric bootstrap

For each family, simulate B alignments from the K_best model, refit the
full K range on every replicate, and compute power.

### Output

The result is a ranked comparison table of all families at their respective
K_best, annotated with power estimates under all three ICs.  This directly
answers: *which model family and complexity does the alignment best support,
and how reliably?*

```r
surv <- kpower_survey(aln, K_max = 5, mix_types = c("+R", "*R", "+H", "*H", "+T", "*T"))
print(surv)
#   mix_type  K  df      lnL      AIC     AICc      BIC  power
#         +R  3  47  -12345.6  24785.2  24790.1  25012.3  87.0%
#         *R  2  52  -12340.1  24784.2  24790.8  25044.5  72.0%
#         +H  3  77  -12310.2  24774.4  24785.1  25089.7  65.0%
#         ...
```

---

## All-IC Power Reporting

Both `kpower()` and `kpower_survey()` always compute K_best and power
under all three information criteria (AIC, AICc, BIC), regardless of which
one the user selects as primary.  The `$power_all` element in the result
contains a named list:

```r
res$power_all$AIC   # list(K_best = 4, power = 0.82)
res$power_all$AICc  # list(K_best = 4, power = 0.80)
res$power_all$BIC   # list(K_best = 3, power = 0.91)
```

The primary IC (specified by `ic`) still determines `$K_best` and `$power`
in the top-level result, and controls which K the simulation is generated
from.

---

## Primary Output

`kpower()` returns a list (class `kpower_result`) containing:

- `$empirical`: data frame of K, lnL, df, AIC, AICc, BIC for the empirical
  runs
- `$sim_ic`: long-format data frame of replicate, K, lnL, df, AIC, AICc,
  BIC for all B x K simulation fits
- `$K_best`: the K selected from empirical data under the primary IC
- `$power`: proportion of simulations that recover K_best under the primary IC
- `$power_all`: named list with AIC, AICc, BIC entries, each containing
  K_best and power under that criterion
- `$ic`: the primary IC used
- `$mix_type`: mixture family used
- `$plot`: a ggplot2 object (IC profile figure)

For MAST runs, additionally:
- `$rate_model`: auto-determined rate heterogeneity string (e.g. `"+R3"`)
- `$tree_weights`: numeric vector of K_best tree weights
- `$ranked_trees`: integer vector of tree indices in weight order

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
- Crotty et al. (2020) — GHOST heterotachy model
- Woodhams et al. (2024) — MAST tree-mixture model
