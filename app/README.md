# WGCNAmp GUI

This directory contains the local Shiny interface for the modular WGCNAmp
pipeline.

Launch it from the project root:

```bash
Rscript run_gui.R
```

For a headless launch or a fixed port:

```bash
Rscript run_gui.R --no-browser --port=3838
```

Uploads are limited to 512 MB by default. Set `WGCNAMP_MAX_UPLOAD_MB` before
launching to choose a different limit.

Network construction runs in a separate R process. Set
`WGCNAMP_BACKGROUND_WORKERS` before launch to change the number of concurrent
background workers; the default is one to avoid overlapping memory-intensive
network jobs.

The current milestone provides:

- discovery and selection of existing analysis sessions;
- safe creation of self-contained GUI sessions whose configuration is stored
  in `session_state.rds` and stage reports;
- sidebar navigation for the seven pipeline stages;
- stage status inferred from existing artifacts;
- recent-artifact visibility;
- a complete data-preparation and QC screen with session uploads, parsing
  controls, reusable filtering settings, sample exclusions, and QC grouping;
- execution of `run_data_preparation()` with the legacy
  `out/1_prepared_data.RData` artifact retained for downstream stages;
- persistent preparation metrics that distinguish import, missing-value QC,
  low-expression filtering, and sample exclusion, plus PCA and
  sample-clustering previews;
- a soft-threshold screen with configurable candidate powers, network type,
  scale-free fit target, explicit selection policy, diagnostics, and the legacy
  `out/2_pt.Rdata` artifact;
- background network construction with blockwise and step-by-step methods,
  retained method-specific block-size and adjacency settings, persisted task
  state and logs, method-compatible TOM manifests, downloadable block
  dendrogram/eigengene diagnostics, and standardized network artifacts;
- dynamic module–trait controls populated from each session's traits, metadata,
  and available networks, including reference-level selection, automatic
  single-level exclusion, contrast previews, correlation tables, heatmaps, and
  covariate-adjusted categorical linear models with confidence intervals and
  within-network FDR correction;
- WGCNA hub-gene ranking with dynamic network, module, threshold, and optional
  gene-significance controls, using either scalable adjacency connectivity or
  notebook-compatible stored TOM connectivity, plus bounded TOM heatmaps and
  weighted Cytoscape node/edge exports;
- Arabidopsis GO and KEGG enrichment with selectable modules, database-specific
  controls, TAIR-to-Entrez mapping diagnostics, cutoffs, result tables, and
  separate plots and artifacts;
- dynamic normalized-expression views with selectable genes, notebook-faithful
  individual-gene module trajectories, module summaries, configurable sample
  labels, and metadata-annotated eigengene heatmaps;
- versioned persistent session state that restores the last open stage and
  records each run's parameters, selected source files, upstream artifacts,
  and outputs;
- stale-step detection that propagates upstream input or artifact changes to
  dependent stages while preserving their existing results for inspection.

The Shiny application is the supported user interface for all seven stages.
Legacy Quarto notebooks are retained only as implementation and reproducibility
references.
