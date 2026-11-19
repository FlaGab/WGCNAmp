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

The current milestone provides:

- discovery and selection of existing analysis sessions;
- safe creation of a new session from `set_parameters.R`;
- sidebar navigation for the seven pipeline stages;
- stage status inferred from existing artifacts;
- recent-artifact visibility;
- placeholders where each notebook's controls and execution action will be
  migrated.

The GUI does not execute pipeline stages yet. Existing notebooks and analysis
artifacts remain unchanged.
