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

The current milestone provides:

- discovery and selection of existing analysis sessions;
- safe creation of a new session from `set_parameters.R`;
- sidebar navigation for the seven pipeline stages;
- stage status inferred from existing artifacts;
- recent-artifact visibility;
- a complete data-preparation and QC screen with session uploads, parsing
  controls, reusable filtering settings, sample exclusions, and QC grouping;
- execution of `run_data_preparation()` with the legacy
  `out/1_prepared_data.RData` artifact retained for downstream stages;
- persistent preparation metrics plus PCA and sample-clustering previews;
- placeholders where the remaining notebooks' controls and execution actions
  will be migrated.

The GUI currently executes step 1. Existing notebooks remain available and
continue to use the same pipeline artifacts.
