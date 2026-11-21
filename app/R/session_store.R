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
      "GO enrichment",
      "Gene visualization"
    ),
    description = c(
      "Import counts and metadata, inspect quality, filter genes, and normalize expression.",
      "Estimate and select the soft-thresholding power.",
      "Build the WGCNA network and inspect detected modules.",
      "Generate categorical contrasts and relate module eigengenes to traits.",
      "Calculate module membership, connectivity, and candidate hub genes.",
      "Run GO and KEGG enrichment for selected modules.",
      "Explore normalized expression for selected genes and modules."
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

  saveRDS(
    list(
      name = name,
      created_at = Sys.time(),
      app_version = 1L
    ),
    file.path(staging, "session_state.rds")
  )

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
    step5 = any(file.exists(file.path(
      out,
      c("geneMM_hc.csv", "geneMM_km.csv", "intramodular_connectivity_km.csv")
    ))),
    step6 = {
      go_dir <- file.path(out, "GO")
      dir.exists(go_dir) && length(list.files(go_dir, all.files = FALSE)) > 0L
    },
    step7 = {
      length(list.files(plots, pattern = "gene|expression", ignore.case = TRUE)) > 0L
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
  dependencies <- list(
    step1 = character(),
    step2 = "step1",
    step3 = c("step1", "step2"),
    step4 = c("step1", "step3"),
    step5 = c("step1", "step3"),
    step6 = "step3",
    step7 = c("step1", "step3")
  )

  status <- vapply(definitions$id, function(stage_id) {
    if (complete[[stage_id]]) {
      "complete"
    } else if (all(complete[dependencies[[stage_id]]])) {
      "ready"
    } else {
      "blocked"
    }
  }, character(1))

  transform(definitions, status = unname(status))
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
