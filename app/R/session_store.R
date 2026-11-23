wgcnamp_project_root <- function() {
  root <- Sys.getenv("WGCNAMP_PROJECT_ROOT", unset = "")
  if (!nzchar(root)) {
    root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
  }
  normalizePath(root, mustWork = TRUE)
}

stage_definitions <- function() {
  data.frame(
    id = paste0("step", 1:7),
    number = 1:7,
    label = c(
      "Data preparation",
      "Soft threshold",
      "Network construction",
      "Module–trait correlation",
      "Hub genes",
      "Functional enrichment",
      "Gene visualization"
    ),
    description = c(
      "Import counts and metadata, inspect quality, filter genes, and normalize expression.",
      "Estimate and select the soft-thresholding power.",
      "Build the WGCNA network and inspect detected modules.",
      "Relate eigengenes to traits with correlations and covariate-adjusted models.",
      "Calculate candidate hub genes, visualize TOM connectivity, and export Cytoscape networks.",
      "Run GO and KEGG enrichment for selected modules.",
      "Explore gene trajectories, module summaries, and metadata-annotated eigengene heatmaps."
    ),
    notebook = c(
      "1_data_preparation.qmd",
      "2_sft_power_estim.qmd",
      "3_WGCN-auto.qmd",
      "4_module-trait_cor.qmd",
      "5_hub_genes.qmd",
      "6_modules_GO_enrichment.qmd",
      "7_visualization_genes.qmd"
    ),
    stringsAsFactors = FALSE
  )
}

stage_dependencies <- function() {
  list(
    step1 = character(),
    step2 = "step1",
    step3 = c("step1", "step2"),
    step4 = c("step1", "step3"),
    step5 = c("step1", "step3"),
    step6 = c("step1", "step3"),
    step7 = c("step1", "step3")
  )
}

session_state_path <- function(session_path) {
  file.path(session_path, "session_state.rds")
}

new_session_state <- function(name, created_at = Sys.time()) {
  list(
    schema_version = 2L,
    name = name,
    created_at = created_at,
    updated_at = Sys.time(),
    app_version = 2L,
    artifacts_baselined = TRUE,
    ui = list(active_stage = NULL),
    stages = list()
  )
}

write_session_state <- function(session_path, state) {
  state$schema_version <- 2L
  state$app_version <- 2L
  state$updated_at <- Sys.time()
  target <- session_state_path(session_path)
  temporary <- tempfile(pattern = ".session-state-", tmpdir = session_path)
  saveRDS(state, temporary)
  if (!file.rename(temporary, target)) {
    copied <- file.copy(temporary, target, overwrite = TRUE)
    unlink(temporary)
    if (!copied) stop("Could not persist session state.", call. = FALSE)
  }
  invisible(state)
}

relative_session_path <- function(path, session_path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  root <- paste0(normalizePath(session_path, mustWork = TRUE), .Platform$file.sep)
  if (startsWith(normalized, root)) substring(normalized, nchar(root) + 1L) else normalized
}

fingerprint_files <- function(paths, session_path) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  if (length(paths) == 0L) {
    return(data.frame(
      path = character(), exists = logical(), size = numeric(),
      modified = numeric(), stringsAsFactors = FALSE
    ))
  }
  information <- file.info(paths)
  data.frame(
    path = vapply(paths, relative_session_path, character(1), session_path = session_path),
    exists = file.exists(paths),
    size = ifelse(file.exists(paths), information$size, NA_real_),
    modified = ifelse(file.exists(paths), as.numeric(information$mtime), NA_real_),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

stage_artifact_paths <- function(stage_id, session_path) {
  out <- file.path(session_path, "out")
  plots <- file.path(session_path, "plots")
  candidates <- switch(
    stage_id,
    step1 = file.path(out, "1_prepared_data.RData"),
    step2 = file.path(out, "2_pt.Rdata"),
    step3 = file.path(out, c("3_networks.rds", "3_networks.Rdata", "networks.Rdata")),
    step4 = c(
      file.path(out, "4_module_trait_report.rds"),
      file.path(out, "module_trait_binary_contrasts.tsv"),
      file.path(plots, c("EMhc-condition_corr.pdf", "EMkm-condition_corr.pdf"))
    ),
    step5 = c(
      file.path(out, "5_hub_report.rds"),
      file.path(out, c("geneMM_hc.csv", "geneMM_km.csv", "intramodular_connectivity_km.csv"))
    ),
    step6 = c(
      file.path(out, "GO", "6_enrichment_report.rds"),
      list.files(file.path(out, "GO"), full.names = TRUE, all.files = FALSE)
    ),
    step7 = c(
      file.path(out, "7_visualization_report.rds"),
      file.path(plots, "7_gene_visualization.pdf")
    ),
    character()
  )
  existing <- unique(candidates[file.exists(candidates)])
  if (length(existing) == 0L) return(character())
  canonical <- existing[grepl(paste0("/", sub("step", "", stage_id), "_.*report[.]rds$"), existing)]
  if (length(canonical) > 0L) canonical[[1L]] else existing[[1L]]
}

stage_output_fingerprint <- function(stage_id, session_path) {
  fingerprint_files(stage_artifact_paths(stage_id, session_path), session_path)
}

step1_source_paths <- function(session_path, parameters = NULL) {
  configured <- parameters$input_files
  if (!is.null(configured)) {
    paths <- unlist(configured, use.names = FALSE)
    paths <- paths[!is.na(paths) & nzchar(paths)]
    resolved <- vapply(paths, function(path) {
      if (file.exists(path)) return(path)
      session_copy <- file.path(session_path, "in", basename(path))
      if (file.exists(session_copy)) session_copy else path
    }, character(1))
    return(unname(resolved))
  }
  input_dir <- file.path(session_path, "in")
  if (!dir.exists(input_dir)) return(character())
  list.files(
    input_dir,
    pattern = "[.](csv|tsv|txt)$",
    ignore.case = TRUE,
    full.names = TRUE
  )
}

stage_result_parameters <- function(stage_id, result) {
  if (is.null(result)) return(NULL)
  if (!is.null(result$parameters)) return(result$parameters)
  if (identical(stage_id, "step2")) {
    return(list(
      network_type = result$network_type,
      powers = result$powers,
      fit_cutoff = result$fit_cutoff,
      selection_mode = result$selection_mode,
      power_threshold = result$power_threshold
    ))
  }
  NULL
}

read_session_state <- function(session_path, migrate = TRUE) {
  path <- session_state_path(session_path)
  state <- if (file.exists(path)) {
    tryCatch(readRDS(path), error = function(error) NULL)
  } else {
    NULL
  }
  state_was_missing <- is.null(state) || !is.list(state)
  legacy_schema <- !state_was_missing &&
    (is.null(state$schema_version) || state$schema_version < 2L)
  needs_artifact_baseline <- state_was_missing || legacy_schema ||
    !isTRUE(state$artifacts_baselined)
  changed <- FALSE
  if (state_was_missing) {
    state <- new_session_state(basename(session_path))
    changed <- TRUE
  }
  if (is.null(state$schema_version) || state$schema_version < 2L) {
    state$schema_version <- 2L
    state$app_version <- 2L
    state$ui <- state$ui %||state% list(active_stage = NULL)
    state$stages <- state$stages %||state% list()
    changed <- TRUE
  }
  if (is.null(state$ui)) {
    state$ui <- list(active_stage = NULL)
    changed <- TRUE
  }
  if (is.null(state$stages)) {
    state$stages <- list()
    changed <- TRUE
  }

  if (isTRUE(migrate) && needs_artifact_baseline) {
    dependencies <- stage_dependencies()
    for (stage_id in names(dependencies)) {
      if (stage_is_complete(stage_id, session_path) &&
          is.null(state$stages[[stage_id]])) {
        outputs <- stage_output_fingerprint(stage_id, session_path)
        upstream <- stats::setNames(
          lapply(dependencies[[stage_id]], stage_output_fingerprint, session_path = session_path),
          dependencies[[stage_id]]
        )
        state$stages[[stage_id]] <- list(
          completed_at = if (nrow(outputs) > 0L) {
            as.POSIXct(max(outputs$modified), origin = "1970-01-01")
          } else {
            Sys.time()
          },
          parameters = NULL,
          source_inputs = if (identical(stage_id, "step1")) {
            fingerprint_files(step1_source_paths(session_path), session_path)
          } else {
            NULL
          },
          upstream = upstream,
          outputs = outputs,
          legacy_baseline = TRUE
        )
        changed <- TRUE
      }
    }
    state$artifacts_baselined <- TRUE
    changed <- TRUE
  }
  if (changed && isTRUE(migrate)) write_session_state(session_path, state)
  state
}

record_stage_completion <- function(session_path, stage_id, result = NULL) {
  if (!stage_id %in% stage_definitions()$id) {
    stop("Unknown pipeline stage: ", stage_id, call. = FALSE)
  }
  state <- read_session_state(session_path, migrate = TRUE)
  parameters <- stage_result_parameters(stage_id, result)
  dependencies <- stage_dependencies()[[stage_id]]
  state$stages[[stage_id]] <- list(
    completed_at = result$completed_at %||state% Sys.time(),
    parameters = parameters,
    source_inputs = if (identical(stage_id, "step1")) {
      fingerprint_files(step1_source_paths(session_path, parameters), session_path)
    } else {
      NULL
    },
    upstream = stats::setNames(
      lapply(dependencies, stage_output_fingerprint, session_path = session_path),
      dependencies
    ),
    outputs = stage_output_fingerprint(stage_id, session_path),
    legacy_baseline = FALSE
  )
  write_session_state(session_path, state)
}

update_session_ui_state <- function(session_path, active_stage = NULL) {
  state <- read_session_state(session_path, migrate = TRUE)
  if (!is.null(active_stage)) state$ui$active_stage <- active_stage
  write_session_state(session_path, state)
}

saved_stage_parameters <- function(session_path, stage_id) {
  read_session_state(session_path, migrate = TRUE)$stages[[stage_id]]$parameters
}

fingerprints_equal <- function(left, right) {
  is.data.frame(left) && is.data.frame(right) && identical(left, right)
}

is_analysis_session <- function(path) {
  dir.exists(file.path(path, "in")) &&
    dir.exists(file.path(path, "out")) &&
    dir.exists(file.path(path, "plots"))
}

discover_analysis_sessions <- function(project_root = wgcnamp_project_root()) {
  candidates <- list.dirs(project_root, recursive = FALSE, full.names = TRUE)
  candidates <- candidates[vapply(candidates, is_analysis_session, logical(1))]

  if (length(candidates) == 0L) {
    return(data.frame(
      name = character(),
      path = character(),
      parameter_file = character(),
      stringsAsFactors = FALSE
    ))
  }

  names <- basename(candidates)
  parameter_files <- file.path(candidates, "in", paste0(names, "_parameters.R"))

  result <- data.frame(
    name = names,
    path = normalizePath(candidates, mustWork = TRUE),
    parameter_file = parameter_files,
    stringsAsFactors = FALSE
  )
  result[order(tolower(result$name)), , drop = FALSE]
}

validate_session_name <- function(name) {
  name <- trimws(name)
  if (!nzchar(name)) {
    stop("Enter a session name.", call. = FALSE)
  }
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", name)) {
    stop(
      "Use only letters, numbers, underscores, and hyphens; start with a letter or number.",
      call. = FALSE
    )
  }
  name
}

create_analysis_session <- function(name,
                                    project_root = wgcnamp_project_root()) {
  name <- validate_session_name(name)
  target <- file.path(project_root, name)
  if (file.exists(target) || dir.exists(target)) {
    stop("A file or session named '", name, "' already exists.", call. = FALSE)
  }

  template <- file.path(project_root, "set_parameters.R")
  if (!file.exists(template)) {
    stop("The parameter template 'set_parameters.R' was not found.", call. = FALSE)
  }

  staging <- tempfile(pattern = ".wgcnamp-session-", tmpdir = project_root)
  dir.create(staging)
  committed <- FALSE
  on.exit({
    if (!committed && dir.exists(staging)) {
      unlink(staging, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  directories <- c(
    "in",
    "out",
    file.path("out", "TOMfiles"),
    file.path("out", "cytoscape"),
    file.path("out", "GO"),
    "plots",
    "logs"
  )
  for (directory in directories) {
    dir.create(file.path(staging, directory), recursive = TRUE, showWarnings = FALSE)
  }

  parameter_file <- file.path(staging, "in", paste0(name, "_parameters.R"))
  parameters <- readLines(template, warn = FALSE)
  session_line <- grep("^session_name\\s*<-", parameters)
  if (length(session_line) != 1L) {
    stop("The parameter template must define exactly one 'session_name'.", call. = FALSE)
  }
  parameters[session_line] <- paste0("session_name <- ", encodeString(name, quote = '"'))
  writeLines(parameters, parameter_file)

  saveRDS(new_session_state(name), file.path(staging, "session_state.rds"))

  if (!file.rename(staging, target)) {
    stop("Could not move the new session into the project directory.", call. = FALSE)
  }
  committed <- TRUE

  normalizePath(target, mustWork = TRUE)
}

stage_is_complete <- function(stage_id, session_path) {
  out <- file.path(session_path, "out")
  plots <- file.path(session_path, "plots")

  switch(
    stage_id,
    step1 = file.exists(file.path(out, "1_prepared_data.RData")),
    step2 = file.exists(file.path(out, "2_pt.Rdata")),
    step3 = any(file.exists(file.path(
      out,
      c("3_networks.rds", "3_networks.Rdata", "networks.Rdata")
    ))),
    step4 = file.exists(file.path(out, "4_module_trait_report.rds")) ||
      file.exists(file.path(out, "module_trait_binary_contrasts.tsv")) ||
      any(file.exists(file.path(
        plots,
        c("EMhc-condition_corr.pdf", "EMkm-condition_corr.pdf")
      ))),
    step5 = file.exists(file.path(out, "5_hub_report.rds")) ||
      any(file.exists(file.path(
      out,
      c("geneMM_hc.csv", "geneMM_km.csv", "intramodular_connectivity_km.csv")
    ))),
    step6 = {
      go_dir <- file.path(out, "GO")
      file.exists(file.path(go_dir, "6_enrichment_report.rds")) ||
        (dir.exists(go_dir) && length(list.files(go_dir, all.files = FALSE)) > 0L)
    },
    step7 = {
      file.exists(file.path(out, "7_visualization_report.rds")) ||
        length(list.files(
          plots,
          pattern = "^7_.*visualization|gene[_-]expression|expression[_-]profile",
          ignore.case = TRUE
        )) > 0L
    },
    FALSE
  )
}

session_stage_status <- function(session_path) {
  definitions <- stage_definitions()
  complete <- setNames(
    vapply(definitions$id, stage_is_complete, logical(1), session_path = session_path),
    definitions$id
  )
  dependencies <- stage_dependencies()
  state <- read_session_state(session_path, migrate = TRUE)
  status <- stats::setNames(rep("blocked", nrow(definitions)), definitions$id)
  stale_reason <- stats::setNames(rep(NA_character_, nrow(definitions)), definitions$id)

  for (stage_id in definitions$id) {
    dependencies_fresh <- all(status[dependencies[[stage_id]]] == "complete")
    if (!complete[[stage_id]]) {
      status[[stage_id]] <- if (dependencies_fresh) "ready" else "blocked"
      next
    }
    record <- state$stages[[stage_id]]
    reason <- NULL
    if (is.null(record)) {
      reason <- "No persistent run record is available."
    } else if (!fingerprints_equal(
      record$outputs,
      stage_output_fingerprint(stage_id, session_path)
    )) {
      reason <- "The saved stage artifact changed after this run was recorded."
    } else if (identical(stage_id, "step1") &&
               !fingerprints_equal(
                 record$source_inputs,
                 fingerprint_files(
                   step1_source_paths(session_path, record$parameters),
                   session_path
                 )
               )) {
      reason <- "One or more selected input files changed after data preparation."
    } else if (!dependencies_fresh) {
      stale_dependencies <- dependencies[[stage_id]][
        status[dependencies[[stage_id]]] != "complete"
      ]
      reason <- paste(
        "Upstream stage(s) require rerunning:",
        paste(sub("step", "Step ", stale_dependencies), collapse = ", ")
      )
    } else {
      for (dependency in dependencies[[stage_id]]) {
        recorded <- record$upstream[[dependency]]
        current <- stage_output_fingerprint(dependency, session_path)
        if (!fingerprints_equal(recorded, current)) {
          reason <- paste0(
            sub("step", "Step ", dependency),
            " was rerun after this result was created."
          )
          break
        }
      }
    }
    if (is.null(reason)) {
      status[[stage_id]] <- "complete"
    } else {
      status[[stage_id]] <- "stale"
      stale_reason[[stage_id]] <- reason
    }
  }

  transform(
    definitions,
    status = unname(status),
    stale_reason = unname(stale_reason)
  )
}

stage_status_value <- function(session_path, stage_id) {
  stages <- session_stage_status(session_path)
  stages$status[match(stage_id, stages$id)]
}

stage_stale_reason <- function(session_path, stage_id) {
  stages <- session_stage_status(session_path)
  stages$stale_reason[match(stage_id, stages$id)]
}

stage_dependencies_ready <- function(session_path, stage_id) {
  dependencies <- stage_dependencies()[[stage_id]]
  if (length(dependencies) == 0L) return(TRUE)
  stages <- session_stage_status(session_path)
  all(stages$status[match(dependencies, stages$id)] == "complete")
}

session_output_files <- function(session_path, limit = 8L) {
  roots <- c(file.path(session_path, "out"), file.path(session_path, "plots"))
  files <- unlist(lapply(roots[dir.exists(roots)], function(path) {
    list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  }), use.names = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]

  if (length(files) == 0L) {
    return(character())
  }

  modified <- file.info(files)$mtime
  files <- files[order(modified, decreasing = TRUE, na.last = TRUE)]
  head(sub(paste0("^", normalizePath(session_path), "/?"), "", files), limit)
}

`%||state%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L) right else left
}
