dataPreparationUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::div(
      class = "stage-heading",
      shiny::div(
        shiny::span(class = "stage-number", "Step 1"),
        shiny::h2("Data preparation and quality control"),
        shiny::p(
          class = "stage-description",
          "Import counts and sample annotations, inspect sample quality, filter low-expression genes, and create the normalized expression matrices used downstream."
        )
      ),
      shiny::uiOutput(ns("stage_badge"))
    ),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        class = "preparation-card",
        bslib::card_header("Input data"),
        shiny::div(
          class = "control-section",
          shiny::h3("Counts"),
          shiny::selectInput(ns("counts_file"), "Session file", choices = NULL),
          shiny::fileInput(
            ns("counts_upload"),
            "Upload a counts table",
            accept = c(".csv", ".tsv", ".txt")
          ),
          shiny::div(
            class = "form-grid form-grid-2",
            shiny::selectInput(
              ns("counts_sep"),
              "Column separator",
              choices = separator_choices(),
              selected = "comma"
            ),
            shiny::selectizeInput(
              ns("gene_length_column"),
              "Gene-length column (optional)",
              choices = NULL,
              options = list(placeholder = "No gene-length column")
            )
          )
        ),
        shiny::hr(),
        shiny::div(
          class = "control-section",
          shiny::h3("Sample metadata"),
          shiny::selectInput(ns("metadata_file"), "Session file", choices = NULL),
          shiny::fileInput(
            ns("metadata_upload"),
            "Upload a metadata table",
            accept = c(".csv", ".tsv", ".txt")
          ),
          shiny::div(
            class = "form-grid form-grid-2",
            shiny::selectInput(
              ns("metadata_sep"),
              "Column separator",
              choices = separator_choices(),
              selected = "comma"
            ),
            shiny::textInput(ns("metadata_sample_id"), "Sample-ID column", value = "sample")
          ),
          shiny::selectizeInput(
            ns("qc_color_columns"),
            "Color/group QC plots by",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Optional: select one or more metadata columns",
              plugins = list("remove_button")
            )
          ),
          shiny::p(
            class = "field-help",
            "Selecting several columns creates a combined QC group without changing the source metadata file."
          )
        ),
        shiny::hr(),
        shiny::div(
          class = "control-section",
          shiny::h3("Traits (optional)"),
          shiny::selectInput(ns("traits_file"), "Session file", choices = NULL),
          shiny::fileInput(
            ns("traits_upload"),
            "Upload a traits table",
            accept = c(".csv", ".tsv", ".txt")
          ),
          shiny::div(
            class = "form-grid form-grid-2",
            shiny::selectInput(
              ns("traits_sep"),
              "Column separator",
              choices = separator_choices(),
              selected = "comma"
            ),
            shiny::textInput(ns("traits_sample_id"), "Sample-ID column", value = "sample")
          )
        )
      ),
      shiny::div(
        bslib::card(
          class = "preparation-card",
          bslib::card_header("Filtering and exclusions"),
          shiny::selectizeInput(
            ns("excluded_samples"),
            "Exclude samples",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Keep all samples",
              plugins = list("remove_button")
            )
          ),
          shiny::selectInput(
            ns("filtering_method"),
            "Low-expression rule",
            choices = c(
              "Minimum total reads across all samples" = "1",
              "Minimum reads in a number of samples" = "2",
              "Minimum reads in a proportion of samples" = "3"
            ),
            selected = "2"
          ),
          shiny::numericInput(
            ns("read_threshold"),
            "Minimum reads",
            value = 10,
            min = 0,
            step = 1
          ),
          shiny::uiOutput(ns("method_parameter")),
          shiny::uiOutput(ns("input_overview"))
        ),
        bslib::card(
          class = "run-card",
          bslib::card_header("Run preparation"),
          shiny::uiOutput(ns("dependency_status")),
          shiny::actionButton(
            ns("run"),
            "Run data preparation",
            icon = shiny::icon("play"),
            class = "btn-primary run-stage-button"
          ),
          shiny::p(
            class = "run-note",
            "Running again replaces this session's step-1 artifacts. Downstream outputs are left untouched and should be rerun when the prepared data changes."
          )
        )
      )
    ),
    shiny::uiOutput(ns("run_feedback")),
    shiny::uiOutput(ns("results_panel"))
  )
}

separator_choices <- function() {
  c("Comma (,)" = "comma", "Semicolon (;)" = "semicolon", "Tab" = "tab")
}

separator_value <- function(value) {
  switch(value, comma = ",", semicolon = ";", tab = "\t", ",")
}

separator_key <- function(value) {
  switch(value, `;` = "semicolon", `\t` = "tab", comma = "comma", "," = "comma", "comma")
}

detect_separator_key <- function(path, fallback = "comma") {
  if (is.null(path) || !file.exists(path)) {
    return(fallback)
  }
  header <- readLines(path, n = 1L, warn = FALSE)
  if (length(header) == 0L) {
    return(fallback)
  }
  separators <- c(comma = ",", semicolon = ";", tab = "\t")
  counts <- vapply(separators, function(separator) {
    matches <- gregexpr(separator, header, fixed = TRUE)[[1L]]
    if (length(matches) == 1L && matches[[1L]] == -1L) 0L else length(matches)
  }, integer(1))
  if (max(counts) == 0L) fallback else names(which.max(counts))
}

session_input_files <- function(session_path) {
  input_dir <- file.path(session_path, "in")
  if (!dir.exists(input_dir)) {
    return(character())
  }
  files <- list.files(
    input_dir,
    pattern = "[.](csv|tsv|txt)$",
    ignore.case = TRUE,
    full.names = FALSE
  )
  sort(files)
}

session_input_path <- function(session_path, filename, optional = FALSE) {
  if (is.null(filename) || !nzchar(filename)) {
    if (optional) {
      return(NULL)
    }
    stop("Select an input file.", call. = FALSE)
  }
  if (!identical(filename, basename(filename))) {
    stop("Invalid input filename.", call. = FALSE)
  }
  path <- file.path(session_path, "in", filename)
  if (!file.exists(path)) {
    stop("Input file not found: ", filename, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

read_table_columns <- function(path, sep) {
  if (is.null(path) || !file.exists(path)) {
    return(character())
  }
  names(utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    nrows = 2L,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
}

read_sample_ids <- function(path, sep, sample_id) {
  if (is.null(path) || !file.exists(path) || !nzchar(sample_id)) {
    return(character())
  }
  preview <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!sample_id %in% names(preview)) {
    return(character())
  }
  unique(as.character(preview[[sample_id]]))
}

uploaded_filename <- function(name) {
  name <- basename(name)
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", name)) {
    stop(
      "Uploaded filenames may contain only letters, numbers, dots, underscores, and hyphens.",
      call. = FALSE
    )
  }
  if (!grepl("[.](csv|tsv|txt)$", name, ignore.case = TRUE)) {
    stop("Upload a CSV, TSV, or TXT table.", call. = FALSE)
  }
  name
}

store_session_upload <- function(upload, session_path) {
  name <- uploaded_filename(upload$name)
  destination <- file.path(session_path, "in", name)
  if (file.exists(destination)) {
    stem <- tools::file_path_sans_ext(name)
    extension <- tools::file_ext(name)
    suffix <- format(Sys.time(), "%Y%m%d-%H%M%S")
    name <- paste0(stem, "-", suffix, ".", extension)
    destination <- file.path(session_path, "in", name)
  }
  if (!file.copy(upload$datapath, destination, overwrite = FALSE)) {
    stop("Could not store the uploaded file in this session.", call. = FALSE)
  }
  name
}

data_preparation_packages <- function() {
  c("WGCNA", "dendextend", "ggplot2", "ggrepel", "DESeq2", "SummarizedExperiment")
}

missing_data_preparation_packages <- function() {
  packages <- data_preparation_packages()
  packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
}

load_preparation_report <- function(session_path) {
  path <- file.path(session_path, "out", "1_preparation_report.rds")
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(error) NULL)
}

dataPreparationServer <- function(id, selected_session, on_complete = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    report <- shiny::reactiveVal(NULL)
    run_error <- shiny::reactiveVal(NULL)
    run_completed <- shiny::reactiveVal(FALSE)

    selected_path <- shiny::reactive({
      selected <- selected_session()
      if (is.null(selected)) NULL else selected$path
    })

    refresh_file_choices <- function(selected = list()) {
      path <- selected_path()
      if (is.null(path)) {
        return(invisible(NULL))
      }
      files <- session_input_files(path)
      required_choices <- stats::setNames(files, files)
      optional_choices <- c("None" = "", required_choices)

      if (length(files) == 0L) {
        shiny::updateSelectInput(
          session, "counts_file", choices = character(), selected = character()
        )
        shiny::updateSelectInput(
          session, "metadata_file", choices = character(), selected = character()
        )
        shiny::updateSelectInput(
          session, "traits_file", choices = c("None" = ""), selected = ""
        )
        return(invisible(NULL))
      }

      current_counts <- selected$counts %||% shiny::isolate(input$counts_file)
      current_metadata <- selected$metadata %||% shiny::isolate(input$metadata_file)
      current_traits <- selected$traits %||% shiny::isolate(input$traits_file)
      if (is.null(current_counts) || !current_counts %in% files) {
        current_counts <- files[grepl("count", files, ignore.case = TRUE)][1L] %||% files[1L]
      }
      if (is.null(current_metadata) || !current_metadata %in% files) {
        current_metadata <- files[grepl("metadata|sample", files, ignore.case = TRUE)][1L] %||% files[1L]
      }
      if (is.null(current_traits) || !current_traits %in% files) {
        trait_match <- files[grepl("trait", files, ignore.case = TRUE)]
        current_traits <- if (length(trait_match) > 0L) trait_match[[1L]] else ""
      }

      shiny::updateSelectInput(
        session, "counts_file", choices = required_choices, selected = current_counts
      )
      shiny::updateSelectInput(
        session, "metadata_file", choices = required_choices, selected = current_metadata
      )
      shiny::updateSelectInput(
        session, "traits_file", choices = optional_choices, selected = current_traits
      )
      invisible(NULL)
    }

    shiny::observeEvent(selected_path(), {
      report(load_preparation_report(selected_path()))
      run_error(NULL)
      run_completed(FALSE)
      saved <- saved_stage_parameters(selected_path(), "step1")
      saved_files <- if (is.null(saved$input_files)) list() else {
        lapply(saved$input_files, function(path) {
          if (is.null(path) || !nzchar(path)) "" else basename(path)
        })
      }
      refresh_file_choices(saved_files)
      if (!is.null(saved)) {
        shiny::updateSelectInput(
          session, "counts_sep", selected = separator_key(saved$counts_sep)
        )
        shiny::updateSelectInput(
          session, "metadata_sep", selected = separator_key(saved$metadata_sep)
        )
        shiny::updateSelectInput(
          session, "traits_sep", selected = separator_key(saved$traits_sep)
        )
        shiny::updateTextInput(
          session, "metadata_sample_id", value = saved$metadata_sample_id
        )
        shiny::updateTextInput(
          session, "traits_sample_id", value = saved$traits_sample_id
        )
        shiny::updateSelectInput(
          session, "filtering_method", selected = saved$filtering_method
        )
        shiny::updateNumericInput(
          session, "read_threshold", value = saved$n
        )
        shiny::updateSelectizeInput(
          session,
          "gene_length_column",
          selected = saved$gene_length_column %||% "",
          server = TRUE
        )
        shiny::updateSelectizeInput(
          session,
          "excluded_samples",
          selected = saved$excluded_samples,
          server = TRUE
        )
        qc_columns <- saved$qc_color_columns %||% saved$qc_color_column %||% character()
        shiny::updateSelectizeInput(
          session, "qc_color_columns", selected = qc_columns, server = TRUE
        )
      }
    }, ignoreNULL = TRUE)

    observe_upload <- function(input_id, target) {
      shiny::observeEvent(input[[input_id]], {
        tryCatch({
          name <- store_session_upload(input[[input_id]], selected_path())
          selection <- stats::setNames(list(name), target)
          refresh_file_choices(selection)
          shiny::showNotification(
            paste0("Stored ", name, " in this session."),
            type = "message"
          )
        }, error = function(error) {
          shiny::showNotification(conditionMessage(error), type = "error", duration = NULL)
        })
      }, ignoreInit = TRUE)
    }
    observe_upload("counts_upload", "counts")
    observe_upload("metadata_upload", "metadata")
    observe_upload("traits_upload", "traits")

    observe_separator <- function(file_input, separator_input, optional = FALSE) {
      shiny::observeEvent(
        list(selected_path(), input[[file_input]]),
        {
          path <- tryCatch(
            session_input_path(
              selected_path(), input[[file_input]], optional = optional
            ),
            error = function(error) NULL
          )
          if (!is.null(path)) {
            shiny::updateSelectInput(
              session,
              separator_input,
              selected = detect_separator_key(path)
            )
          }
        },
        ignoreInit = TRUE
      )
    }
    observe_separator("counts_file", "counts_sep")
    observe_separator("metadata_file", "metadata_sep")
    observe_separator("traits_file", "traits_sep", optional = TRUE)

    metadata_columns <- shiny::reactive({
      path <- tryCatch(
        session_input_path(selected_path(), input$metadata_file),
        error = function(error) NULL
      )
      tryCatch(
        read_table_columns(path, separator_value(input$metadata_sep)),
        error = function(error) character()
      )
    })

    sample_ids <- shiny::reactive({
      path <- tryCatch(
        session_input_path(selected_path(), input$metadata_file),
        error = function(error) NULL
      )
      tryCatch(
        read_sample_ids(
          path,
          separator_value(input$metadata_sep),
          input$metadata_sample_id
        ),
        error = function(error) character()
      )
    })

    shiny::observeEvent(metadata_columns(), {
      choices <- setdiff(metadata_columns(), input$metadata_sample_id)
      selected <- intersect(shiny::isolate(input$qc_color_columns), choices)
      shiny::updateSelectizeInput(
        session,
        "qc_color_columns",
        choices = choices,
        selected = selected,
        server = TRUE
      )
    })

    shiny::observeEvent(sample_ids(), {
      choices <- sample_ids()
      selected <- intersect(shiny::isolate(input$excluded_samples), choices)
      shiny::updateSelectizeInput(
        session,
        "excluded_samples",
        choices = choices,
        selected = selected,
        server = TRUE
      )
    })

    shiny::observeEvent(
      list(input$counts_file, input$counts_sep),
      {
        path <- tryCatch(
          session_input_path(selected_path(), input$counts_file),
          error = function(error) NULL
        )
        columns <- tryCatch(
          read_table_columns(path, separator_value(input$counts_sep)),
          error = function(error) character()
        )
        columns <- columns[nzchar(columns)]
        selected <- shiny::isolate(input$gene_length_column)
        if (is.null(selected) || !selected %in% columns) {
          length_match <- columns[grepl("length", columns, ignore.case = TRUE)]
          selected <- if (length(length_match) > 0L) length_match[[1L]] else ""
        }
        shiny::updateSelectizeInput(
          session,
          "gene_length_column",
          choices = c("None" = "", stats::setNames(columns, columns)),
          selected = selected,
          server = TRUE
        )
      },
      ignoreInit = TRUE
    )

    output$method_parameter <- shiny::renderUI({
      saved <- if (is.null(selected_path())) NULL else {
        saved_stage_parameters(selected_path(), "step1")
      }
      if (identical(input$filtering_method, "2")) {
        shiny::numericInput(
          session$ns("minimum_samples"),
          "Minimum number of samples",
          value = saved$m %||% 3,
          min = 1,
          step = 1
        )
      } else if (identical(input$filtering_method, "3")) {
        shiny::sliderInput(
          session$ns("minimum_proportion"),
          "Minimum proportion of samples",
          min = 0.05,
          max = 1,
          value = saved$p %||% 0.75,
          step = 0.05
        )
      }
    })

    output$input_overview <- shiny::renderUI({
      if (is.null(selected_path())) {
        return(NULL)
      }
      ids <- sample_ids()
      columns <- metadata_columns()
      shiny::div(
        class = "input-overview",
        shiny::div(
          shiny::span(class = "overview-value", length(ids)),
          shiny::span("samples detected")
        ),
        shiny::div(
          shiny::span(class = "overview-value", max(0L, length(columns) - 1L)),
          shiny::span("metadata fields")
        )
      )
    })

    output$dependency_status <- shiny::renderUI({
      missing <- missing_data_preparation_packages()
      if (length(missing) == 0L) {
        shiny::div(
          class = "dependency-ready",
          shiny::icon("circle-check"),
          "Analysis dependencies are available."
        )
      } else {
        shiny::div(
          class = "dependency-missing",
          shiny::icon("triangle-exclamation"),
          paste("Missing R packages:", paste(missing, collapse = ", "))
        )
      }
    })

    output$stage_badge <- shiny::renderUI({
      report()
      path <- selected_path()
      if (is.null(path)) return(status_badge("blocked"))
      session_stage_badge(path, "step1")
    })

    shiny::observeEvent(input$run, {
      run_error(NULL)
      run_completed(FALSE)
      missing <- missing_data_preparation_packages()
      if (length(missing) > 0L) {
        run_error(paste("Install the required R packages:", paste(missing, collapse = ", ")))
        return()
      }

      tryCatch({
        shiny::req(selected_path())
        counts_path <- session_input_path(selected_path(), input$counts_file)
        metadata_path <- session_input_path(selected_path(), input$metadata_file)
        traits_path <- session_input_path(
          selected_path(), input$traits_file, optional = TRUE
        )

        qc_columns <- input$qc_color_columns %||% character()
        qc_column <- if (length(qc_columns) == 1L) qc_columns else NULL
        metadata_transform <- NULL
        if (length(qc_columns) > 1L) {
          metadata_transform <- function(metadata) {
            missing_columns <- setdiff(qc_columns, colnames(metadata))
            if (length(missing_columns) > 0L) {
              stop(
                "QC grouping columns are missing from metadata: ",
                paste(missing_columns, collapse = ", "),
                call. = FALSE
              )
            }
            metadata$.qc_group <- do.call(
              interaction,
              c(
                unname(metadata[, qc_columns, drop = FALSE]),
                list(drop = TRUE, sep = " · ")
              )
            )
            metadata
          }
          qc_column <- ".qc_group"
        }

        m_value <- if (identical(input$filtering_method, "2")) {
          input$minimum_samples %||% 3
        } else {
          1
        }
        p_value <- if (identical(input$filtering_method, "3")) {
          input$minimum_proportion %||% 0.75
        } else {
          0.5
        }

        result <- shiny::withProgress(
          message = "Preparing expression data",
          detail = "Importing, checking, filtering, and normalizing…",
          value = 0.15,
          {
            value <- run_data_preparation(
              session_dir = selected_path(),
              counts_file = counts_path,
              metadata_file = metadata_path,
              gene_length_column = if (nzchar(input$gene_length_column %||% "")) {
                input$gene_length_column
              } else {
                NULL
              },
              traits_file = traits_path,
              excluded_samples = input$excluded_samples %||% character(),
              filtering_method = input$filtering_method,
              n = input$read_threshold,
              m = m_value,
              p = p_value,
              counts_sep = separator_value(input$counts_sep),
              metadata_sep = separator_value(input$metadata_sep),
              metadata_sample_id = input$metadata_sample_id,
              traits_sep = separator_value(input$traits_sep),
              traits_sample_id = input$traits_sample_id,
              qc_color_column = qc_column,
              metadata_transform = metadata_transform,
              create_plots = TRUE,
              save_artifact = TRUE,
              verbose = TRUE
            )
            shiny::incProgress(0.85, detail = "Saving preparation report…")
            value
          }
        )
        result$parameters$qc_color_columns <- qc_columns
        report(load_preparation_report(selected_path()))
        run_completed(TRUE)
        if (is.function(on_complete)) {
          on_complete(result, selected_path())
        }
        shiny::showNotification(
          "Data preparation completed successfully.",
          type = "message",
          duration = 6
        )
      }, error = function(error) {
        run_error(conditionMessage(error))
        shiny::showNotification(
          paste("Data preparation failed:", conditionMessage(error)),
          type = "error",
          duration = NULL
        )
      })
    })

    output$run_feedback <- shiny::renderUI({
      error <- run_error()
      if (!is.null(error)) {
        return(shiny::div(
          class = "run-feedback run-feedback-error",
          shiny::icon("circle-xmark"),
          shiny::div(shiny::strong("Preparation did not complete"), shiny::p(error))
        ))
      }
      if (isTRUE(run_completed())) {
        return(shiny::div(
          class = "run-feedback run-feedback-success",
          shiny::icon("circle-check"),
          shiny::div(
            shiny::strong("Prepared data saved"),
            shiny::p("Step 1 is complete. Review the metrics and QC plots below before continuing.")
          )
        ))
      }
      NULL
    })

    output$results_panel <- shiny::renderUI({
      if (is.null(report())) {
        return(shiny::div(
          class = "empty-results",
          shiny::icon("chart-area"),
          shiny::div(
            shiny::strong("No preparation report yet"),
            shiny::p("Configure the inputs and run step 1 to see its QC results here.")
          )
        ))
      }
      shiny::tagList(
        shiny::div(
          class = "results-heading",
          shiny::div(
            shiny::span(class = "eyebrow", "Latest result"),
            shiny::h2("Preparation summary")
          ),
          if (!is.null(report()$completed_at)) {
            shiny::span(
              class = "result-time",
              format(report()$completed_at, "%Y-%m-%d %H:%M")
            )
          }
        ),
        bslib::layout_columns(
          col_widths = c(4, 8),
          shiny::uiOutput(session$ns("summary_card")),
          shiny::div(
            class = "qc-card",
            bslib::navset_card_tab(
              title = "Quality-control plots",
              footer = shiny::uiOutput(session$ns("plot_downloads")),
              bslib::nav_panel(
                "PCA",
                shiny::plotOutput(session$ns("pca_plot"), height = "430px")
              ),
              bslib::nav_panel(
                "Sample clustering",
                shiny::plotOutput(
                  session$ns("clustering_plot"), height = "430px"
                )
              )
            )
          )
        )
      )
    })

    output$summary_card <- shiny::renderUI({
      current <- report()
      if (is.null(current)) {
        return(NULL)
      }
      summary <- current$summary
      metric <- function(value, label) {
        shiny::div(
          class = "summary-metric",
          shiny::strong(
            if (is.null(value) || length(value) == 0L || is.na(value)) {
              "—"
            } else {
              format(value, big.mark = ",", scientific = FALSE)
            }
          ),
          shiny::span(label)
        )
      }
      difference <- function(larger, smaller) {
        if (is.null(larger) || is.null(smaller)) return(NULL)
        max(0, as.numeric(larger) - as.numeric(smaller))
      }
      genes_after_missing_qc <- summary$genes_after_missing_qc %||%
        summary$retained_genes
      removed_missing_qc <- difference(
        summary$input_genes,
        genes_after_missing_qc
      )
      removed_expression_filter <- difference(
        genes_after_missing_qc,
        summary$retained_genes
      )
      excluded_sample_count <- difference(
        summary$input_samples,
        summary$retained_samples
      )
      excluded <- unique(c(
        summary$automatic_sample_exclusions,
        summary$configured_sample_exclusions
      ))
      bslib::card(
        class = "summary-card",
        bslib::card_header("Preparation metrics"),
        shiny::h3("Genes"),
        shiny::div(
          class = "summary-metrics",
          metric(summary$input_genes, "genes imported"),
          metric(genes_after_missing_qc, "after missing-value QC"),
          metric(summary$retained_genes, "genes retained"),
          metric(summary$removed_genes, "genes removed in total"),
          metric(removed_missing_qc, "removed by missing-value QC"),
          metric(removed_expression_filter, "removed by low-expression filter")
        ),
        shiny::hr(),
        shiny::h3("Samples"),
        shiny::div(
          class = "summary-metrics",
          metric(summary$input_samples, "samples imported"),
          metric(summary$retained_samples, "samples retained"),
          metric(excluded_sample_count, "samples excluded")
        ),
        shiny::hr(),
        shiny::h3("Sample exclusions"),
        if (length(excluded) == 0L) {
          shiny::p(class = "text-secondary", "No samples were excluded.")
        } else {
          shiny::div(
            class = "sample-tags",
            lapply(excluded, function(value) shiny::span(value))
          )
        },
        shiny::hr(),
        shiny::p(
          class = "artifact-path",
          shiny::strong("Prepared artifact"),
          shiny::code("out/1_prepared_data.RData")
        )
      )
    })

    output$pca_plot <- shiny::renderPlot({
      shiny::req(report(), report()$qc$pca)
      print(report()$qc$pca)
    }, res = 110)

    output$clustering_plot <- shiny::renderPlot({
      shiny::req(report(), report()$qc$clustering)
      graphics::plot(
        report()$qc$clustering,
        main = "Hierarchical clustering of samples",
        xlab = "Samples",
        sub = "",
        cex = 0.8
      )
    }, res = 110)

    output$plot_downloads <- shiny::renderUI({
      current <- report()
      if (is.null(current)) {
        return(NULL)
      }
      shiny::div(
        class = "plot-downloads",
        shiny::downloadButton(session$ns("download_pca"), "Download PCA PDF"),
        shiny::downloadButton(
          session$ns("download_clustering"), "Download clustering PDF"
        )
      )
    })

    output$download_pca <- shiny::downloadHandler(
      filename = function() "PCA_QC.pdf",
      content = function(file) {
        source <- report()$files$pca_qc
        shiny::req(source, file.exists(source))
        file.copy(source, file, overwrite = TRUE)
      }
    )

    output$download_clustering <- shiny::downloadHandler(
      filename = function() "clustering_QC.pdf",
      content = function(file) {
        source <- report()$files$clustering_qc
        shiny::req(source, file.exists(source))
        file.copy(source, file, overwrite = TRUE)
      }
    )
  })
}
