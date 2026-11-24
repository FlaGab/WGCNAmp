# WGCNAmp

WGCNAmp is a local Shiny application for running a modular weighted gene
co-expression network analysis with the R `WGCNA` package. The graphical
interface is the supported user workflow: an analysis can be configured,
executed, repeated, and inspected without editing notebooks or R source files.

## Launch the application

From the project root, run:

```bash
Rscript run_gui.R
```

The application opens in the default browser. To launch it without opening a
browser, or to use a fixed port, run:

```bash
Rscript run_gui.R --no-browser --port=3838
```

## Canonical workflow

### 1. Create or select a session

Use the **Analysis session** selector in the sidebar to open an existing
analysis. Select **New session** to create an isolated workspace for a new
analysis.

Each session keeps its uploaded inputs, saved settings, intermediate objects,
plots, logs, and downloadable results separate from other analyses. Returning
to a session restores its recorded pipeline state. GUI configuration is stored
in `session_state.rds` and the stage reports; creating a session does not
require or generate an editable parameter script. Parameter files retained by
older sessions are ignored by the GUI.

### 2. Upload and prepare data in Stage 1

Open **1 Data preparation** and select or upload:

- a gene-count table;
- a sample-metadata table;
- an optional trait table for downstream module–trait analysis.

Configure the table separators and sample identifiers, choose the quality-
control grouping variables and filtering rule, and optionally exclude samples.
Run the stage, then review the preparation metrics, PCA, and sample-clustering
diagnostics before continuing.

### 3. Configure and run Stages 2–7

Use the sidebar to move through the pipeline. Every stage exposes its settings
in the GUI and can be run again after changing them.

1. **Data preparation** — import, quality-check, filter, and normalize the
   expression data.
2. **Soft threshold** — evaluate candidate powers and select the WGCNA
   soft-thresholding power.
3. **Network construction** — build a blockwise or step-by-step network and
   inspect modules, dendrograms, and eigengene diagnostics.
4. **Module–trait correlation** — relate module eigengenes to continuous and
   categorical traits, with optional covariate adjustment.
5. **Hub genes** — rank intramodular hubs, inspect connectivity, and export
   Cytoscape networks.
6. **Functional enrichment** — run GO or KEGG enrichment for selected modules.
7. **Gene visualization** — explore expression trajectories, module summaries,
   and metadata-annotated eigengene heatmaps.

Stages that require upstream results remain unavailable until their
prerequisites have completed.

### 4. Respond to stale-stage warnings

WGCNAmp records the inputs, parameters, upstream artifacts, and outputs used by
each completed stage. A stage is marked **stale** when an earlier result or one
of its source inputs has changed.

When this happens:

1. rerun the earliest stale stage after reviewing its settings;
2. inspect the new result;
3. rerun any downstream stages that remain marked stale.

Existing downstream results are preserved for inspection until they are
replaced, but they should not be interpreted as current after an upstream
change.

### 5. Download and locate results

Use the download controls beside plots, tables, logs, and exports in each
stage. The same session artifacts are also retained on disk:

```text
<session>/in/       uploaded source tables
<session>/out/      analysis objects, tables, and exports
<session>/plots/    generated figures
<session>/logs/     background-task logs
```

Running a stage again replaces that stage's current outputs. Download or copy
results that must be retained before rerunning with different settings.

## Notebook status

The Shiny application is the canonical and supported user workflow for all
seven stages. The Quarto notebooks in the repository are retained only as
legacy implementation and reproducibility references; users should not edit
them to configure or run an analysis. They are archived under
[`legacy/notebooks/`](legacy/notebooks/).

## Programmatic API

Advanced users should call the extracted pipeline functions from their own R
scripts rather than modifying the archived notebooks. Each function accepts a
session directory and explicit stage settings, and writes artifacts compatible
with the GUI.

| Stage | Function | Implementation |
|---|---|---|
| Data preparation | `run_data_preparation()` | `src/data_preparation_pipeline.R` |
| Soft threshold | `run_soft_threshold_selection()` | `src/soft_threshold_pipeline.R` |
| Network construction | `run_network_construction()` | `src/network_construction_pipeline.R` |
| Module–trait analysis | `run_module_trait_analysis()` | `src/module_trait_pipeline.R` |
| Hub genes | `run_hub_analysis()` | `src/downstream_pipeline.R` |
| Functional enrichment | `run_module_enrichment()` | `src/downstream_pipeline.R` |
| Gene visualization | `run_gene_visualization()` | `src/downstream_pipeline.R` |

Additional GUI implementation notes are available in
[`app/README.md`](app/README.md).
