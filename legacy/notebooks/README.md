# Legacy Quarto notebooks

> Legacy developer reference. These notebooks are not part of the supported
> user workflow and should not be edited to run an analysis. Use the Shiny GUI.

This directory preserves the original stage-oriented Quarto material for
implementation history and reproducibility comparisons. Notebook execution is
disabled in each file, and the examples no longer load `set_parameters.R` or
contain editable session names.

The supported workflow is documented in the project-level
[`README.md`](../../README.md). Advanced users should call the maintained
programmatic entry points directly:

| Stage | Function |
|---|---|
| 1 | `run_data_preparation()` |
| 2 | `run_soft_threshold_selection()` |
| 3 | `run_network_construction()` |
| 4 | `run_module_trait_analysis()` |
| 5 | `run_hub_analysis()` |
| 6 | `run_module_enrichment()` |
| 7 | `run_gene_visualization()` |
