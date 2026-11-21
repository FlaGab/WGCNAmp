networkConstructionUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::div(
      class = "stage-heading",
      shiny::div(
        shiny::span(class = "stage-number", "Step 3"),
        shiny::h2("Network construction"),
        shiny::p(
          class = "stage-description",
          "Build the weighted co-expression network in a separate R process, preserving a responsive interface while WGCNA computes adjacency, TOM, and modules."
        )
      ),
      shiny::uiOutput(ns("stage_badge"))
    ),
    shiny::uiOutput(ns("prerequisite_notice")),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        class = "network-controls-card",
        bslib::card_header("Construction settings"),
        shiny::radioButtons(
          ns("method"),
          "Construction method",
          choices = c(
            "Blockwise modules (recommended)" = "blockwise",
            "Step-by-step adjacency and TOM" = "stepwise"
          ),
          selected = "blockwise"
        ),
        shiny::uiOutput(ns("method_controls")),
        shiny::div(
          class = "form-grid form-grid-2",
          shiny::selectInput(
            ns("correlation_type"),
            "Correlation",
            choices = c("Pearson" = "pearson", "Biweight midcorrelation" = "bicor"),
            selected = "pearson"
          ),
          shiny::selectInput(
            ns("tom_type"),
            "TOM type",
            choices = c(
              "Match network type" = "auto",
              "Signed" = "signed",
              "Unsigned" = "unsigned"
            ),
            selected = "auto"
          ),
          shiny::numericInput(
            ns("min_module_size"),
            "Minimum module size",
            value = 30,
            min = 2,
            step = 1
          ),
          shiny::selectInput(
            ns("deep_split"),
            "Module splitting sensitivity",
            choices = stats::setNames(0:4, c("0 — conservative", "1", "2 — balanced", "3", "4 — aggressive")),
            selected = 2
          ),
          shiny::sliderInput(
            ns("merge_cut_height"),
            "Module merge height",
            min = 0.05,
            max = 0.50,
            value = 0.30,
            step = 0.01
          ),
          shiny::numericInput(
            ns("random_seed"),
            "Random seed",
            value = 42,
            step = 1
          )
        ),
        shiny::checkboxInput(
          ns("kmeans_refinement"),
          "Apply the existing CoExpNets k-means module refinement",
          value = TRUE
        ),
        shiny::p(
          class = "field-help",
          "Disable refinement only if downstream analyses will use hierarchical-clustering modules exclusively."
        )
      ),
      shiny::div(
        bslib::card(
          class = "network-input-card",
          bslib::card_header("Inputs from previous steps"),
          shiny::uiOutput(ns("input_summary"))
        ),
        bslib::card(
          class = "run-card network-run-card",
          bslib::card_header("Background task"),
          shiny::uiOutput(ns("dependency_status")),
          shiny::uiOutput(ns("run_control")),
          shiny::p(
            class = "run-note",
            "The task continues in a separate R process while you use the interface. Running again replaces stage-3 network artifacts; downstream results should then be recomputed."
          )
        ),
        shiny::uiOutput(ns("task_status"))
      )
    ),
    shiny::uiOutput(ns("results_panel"))
  )
}

network_background_packages <- function(kmeans_refinement = TRUE) {
  packages <- c(
    "WGCNA", "dynamicTreeCut", "ggplot2", "future", "promises"
  )
  if (isTRUE(kmeans_refinement)) {
    packages <- c(packages, "CoExpNets")
  }
  packages
}

missing_network_background_packages <- function(kmeans_refinement = TRUE) {
  packages <- network_background_packages(kmeans_refinement)
  packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
}

network_input_summary <- function(session_path) {
  prepared <- file.path(session_path, "out", "1_prepared_data.RData")
  threshold <- file.path(session_path, "out", "2_pt.Rdata")
  if (!file.exists(prepared) || !file.exists(threshold)) {
    return(NULL)
  }
  environment <- new.env(parent = emptyenv())
  load(threshold, envir = environment)
  if (!exists("power_threshold", envir = environment, inherits = FALSE) ||
      !exists("network_type", envir = environment, inherits = FALSE)) {
    return(NULL)
  }
  preparation <- load_preparation_report(session_path)
  list(
    power = environment$power_threshold,
    network_type = environment$network_type,
    samples = if (is.null(preparation)) NULL else preparation$summary$retained_samples,
    genes = if (is.null(preparation)) NULL else preparation$summary$retained_genes
  )
}

networkConstructionServer <- function(id,
                                      selected_session,
                                      project_root,
                                      on_complete = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    report <- shiny::reactiveVal(NULL)
    run_error <- shiny::reactiveVal(NULL)
    running_job <- shiny::reactiveVal(NULL)
    handled_job <- shiny::reactiveVal(NULL)

    selected_path <- shiny::reactive({
      selected <- selected_session()
      if (is.null(selected)) NULL else selected$path
    })

    prerequisites_ready <- shiny::reactive({
      path <- selected_path()
      !is.null(path) &&
        file.exists(file.path(path, "out", "1_prepared_data.RData")) &&
        file.exists(file.path(path, "out", "2_pt.Rdata"))
    })

    task_state <- shiny::reactive({
      shiny::invalidateLater(1000, session)
      path <- selected_path()
      if (is.null(path)) NULL else read_network_task_state(path)
    })

    network_task <- shiny::ExtendedTask$new(function(config) {
      promises::future_promise({
        backend <- new.env(parent = globalenv())
        sys.source(config$backend_file, envir = backend)
        do.call(backend$run_network_construction_job, config$arguments)
      }, seed = TRUE)
    })
    bslib::bind_task_button(network_task, "run", session = session)

    shiny::observeEvent(selected_path(), {
      current <- load_network_report(selected_path())
      report(current)
      run_error(NULL)
      if (!is.null(current)) {
        parameters <- current$parameters
        shiny::updateRadioButtons(session, "method", selected = parameters$method)
        shiny::updateSelectInput(
          session, "correlation_type", selected = parameters$correlation_type
        )
        shiny::updateSelectInput(session, "tom_type", selected = parameters$tom_type)
        shiny::updateNumericInput(
          session, "min_module_size", value = parameters$min_module_size
        )
        shiny::updateSelectInput(session, "deep_split", selected = parameters$deep_split)
        shiny::updateSliderInput(
          session, "merge_cut_height", value = parameters$merge_cut_height
        )
        shiny::updateNumericInput(
          session, "random_seed", value = parameters$random_seed
        )
        shiny::updateCheckboxInput(
          session, "kmeans_refinement", value = parameters$kmeans_refinement
        )
      }
    }, ignoreNULL = TRUE)

    output$stage_badge <- shiny::renderUI({
      report()
      state <- task_state()
      if (!prerequisites_ready()) {
        return(status_badge("blocked"))
      }
      if (!is.null(state) && identical(state$status, "running")) {
        return(shiny::span(class = "status-badge status-running", "Running"))
      }
      complete <- file.exists(file.path(selected_path(), "out", "3_networks.rds"))
      status_badge(if (complete) "complete" else "ready")
    })

    output$prerequisite_notice <- shiny::renderUI({
      if (prerequisites_ready()) {
        return(NULL)
      }
      shiny::div(
        class = "run-feedback prerequisite-warning",
        shiny::icon("lock"),
        shiny::div(
          shiny::strong("Steps 1 and 2 are required"),
          shiny::p(
            "Prepare the expression matrix and save a soft-threshold power before constructing the network."
          )
        )
      )
    })

    output$method_controls <- shiny::renderUI({
      if (identical(input$method, "blockwise")) {
        shiny::tagList(
          shiny::numericInput(
            session$ns("max_block_size"),
            "Maximum genes per block",
            value = 20000,
            min = 2,
            step = 100
          ),
          shiny::p(
            class = "field-help",
            "WGCNA splits larger datasets into blocks and stores one TOM file per block."
          )
        )
      } else {
        shiny::tagList(
          shiny::checkboxInput(
            session$ns("save_adjacency"),
            "Save the adjacency matrix",
            value = TRUE
          ),
          shiny::div(
            class = "selection-guidance network-memory-warning",
            shiny::icon("triangle-exclamation"),
            shiny::span(
              "Step-by-step construction keeps full gene-by-gene matrices in memory. Use blockwise construction for large datasets."
            )
          )
        )
      }
    })

    output$input_summary <- shiny::renderUI({
      summary <- if (is.null(selected_path())) NULL else {
        network_input_summary(selected_path())
      }
      if (is.null(summary)) {
        return(shiny::p(class = "text-secondary", "Previous-stage inputs are not ready."))
      }
      metric <- function(value, label) {
        shiny::div(
          class = "network-input-metric",
          shiny::strong(if (is.null(value)) "—" else format(value, big.mark = ",")),
          shiny::span(label)
        )
      }
      shiny::div(
        class = "network-input-grid",
        metric(summary$power, "soft-threshold power"),
        metric(summary$network_type, "network type"),
        metric(summary$samples, "samples"),
        metric(summary$genes, "genes")
      )
    })

    output$dependency_status <- shiny::renderUI({
      missing <- missing_network_background_packages(input$kmeans_refinement)
      if (length(missing) == 0L) {
        shiny::div(
          class = "dependency-ready",
          shiny::icon("circle-check"),
          "Background execution dependencies are available."
        )
      } else {
        shiny::div(
          class = "dependency-missing",
          shiny::icon("triangle-exclamation"),
          paste("Missing R packages:", paste(missing, collapse = ", "))
        )
      }
    })

    output$run_control <- shiny::renderUI({
      button <- bslib::input_task_button(
        session$ns("run"),
        "Build network",
        icon = shiny::icon("diagram-project"),
        label_busy = "Building network…",
        type = "primary",
        class = "run-stage-button"
      )
      if (!prerequisites_ready() ||
          length(missing_network_background_packages(input$kmeans_refinement)) > 0L) {
        button$attribs$disabled <- "disabled"
      }
      button
    })

    shiny::observeEvent(input$run, {
      tryCatch({
        if (!prerequisites_ready()) {
          stop("Complete steps 1 and 2 before constructing the network.", call. = FALSE)
        }
        missing <- missing_network_background_packages(input$kmeans_refinement)
        if (length(missing) > 0L) {
          stop(
            "Install the required R packages: ", paste(missing, collapse = ", "),
            call. = FALSE
          )
        }
        path <- selected_path()
        job_id <- paste0(
          format(Sys.time(), "%Y%m%d-%H%M%S"), "-", sprintf("%04d", sample.int(9999, 1))
        )
        arguments <- list(
          session_dir = path,
          job_id = job_id,
          method = input$method,
          max_block_size = if (identical(input$method, "blockwise")) {
            input$max_block_size %||% 20000
          } else {
            20000
          },
          tom_type = input$tom_type,
          correlation_type = input$correlation_type,
          merge_cut_height = input$merge_cut_height,
          deep_split = as.integer(input$deep_split),
          min_module_size = input$min_module_size,
          random_seed = input$random_seed,
          kmeans_refinement = isTRUE(input$kmeans_refinement),
          save_adjacency = if (identical(input$method, "stepwise")) {
            isTRUE(input$save_adjacency)
          } else {
            FALSE
          },
          create_plots = TRUE,
          save_artifacts = TRUE,
          verbose = TRUE
        )
        run_error(NULL)
        handled_job(NULL)
        running_job(list(id = job_id, path = path))
        network_task$invoke(list(
          backend_file = file.path(project_root, "src", "network_construction_pipeline.R"),
          arguments = arguments
        ))
      }, error = function(error) {
        run_error(conditionMessage(error))
        shiny::showNotification(conditionMessage(error), type = "error", duration = NULL)
      })
    })

    shiny::observe({
      status <- network_task$status()
      job <- running_job()
      if (is.null(job) || identical(handled_job(), job$id)) {
        return()
      }
      if (identical(status, "success")) {
        value <- network_task$result()
        handled_job(job$id)
        if (identical(selected_path(), job$path)) {
          report(value)
          shiny::showNotification(
            paste0("Network construction completed with ", value$modules_hc, " HC modules."),
            type = "message",
            duration = 8
          )
        }
        if (is.function(on_complete)) {
          on_complete()
        }
      } else if (identical(status, "error")) {
        message <- tryCatch(
          {
            network_task$result()
            "Network construction failed."
          },
          error = function(error) conditionMessage(error)
        )
        handled_job(job$id)
        if (identical(selected_path(), job$path)) {
          run_error(message)
          shiny::showNotification(
            paste("Network construction failed:", message),
            type = "error",
            duration = NULL
          )
        }
      }
    })

    output$task_status <- shiny::renderUI({
      error <- run_error()
      if (!is.null(error)) {
        return(shiny::div(
          class = "network-task-state network-task-error",
          shiny::icon("circle-xmark"),
          shiny::div(shiny::strong("Task failed"), shiny::p(error))
        ))
      }
      state <- task_state()
      if (is.null(state)) {
        return(shiny::div(
          class = "network-task-state network-task-idle",
          shiny::icon("clock"),
          shiny::span("No network task has been started for this session.")
        ))
      }
      if (identical(state$status, "running")) {
        return(shiny::div(
          class = "network-task-state network-task-running",
          shiny::div(class = "network-spinner", role = "status"),
          shiny::div(
            shiny::strong("Network construction is running"),
            shiny::p(
              paste("Started", format(state$started_at, "%Y-%m-%d %H:%M:%S")),
              shiny::br(),
              "You may navigate to other stages while this process continues."
            )
          )
        ))
      }
      if (identical(state$status, "error")) {
        return(shiny::div(
          class = "network-task-state network-task-error",
          shiny::icon("circle-xmark"),
          shiny::div(shiny::strong("Last task failed"), shiny::p(state$error))
        ))
      }
      shiny::div(
        class = "network-task-state network-task-complete",
        shiny::icon("circle-check"),
        shiny::div(
          shiny::strong("Last task completed"),
          shiny::p(format(state$completed_at, "%Y-%m-%d %H:%M:%S"))
        )
      )
    })

    output$results_panel <- shiny::renderUI({
      current <- report()
      if (is.null(current)) {
        return(shiny::div(
          class = "empty-results",
          shiny::icon("diagram-project"),
          shiny::div(
            shiny::strong("No network report yet"),
            shiny::p("Start the background task to construct modules and generate diagnostics.")
          )
        ))
      }
      metric <- function(value, label) {
        shiny::div(
          class = "network-result-metric",
          shiny::strong(if (is.null(value)) "—" else format(value, big.mark = ",")),
          shiny::span(label)
        )
      }
      shiny::tagList(
        shiny::div(
          class = "results-heading",
          shiny::div(
            shiny::span(class = "eyebrow", "Latest result"),
            shiny::h2("Constructed network")
          ),
          shiny::span(
            class = "result-time",
            format(current$completed_at, "%Y-%m-%d %H:%M")
          )
        ),
        shiny::div(
          class = "network-result-grid",
          metric(current$genes, "genes"),
          metric(current$blocks, "blocks"),
          metric(current$modules_hc, "HC modules"),
          metric(current$modules_kmeans, "k-means modules")
        ),
        bslib::layout_columns(
          col_widths = c(8, 4),
          bslib::navset_card_tab(
            title = "Module sizes",
            bslib::nav_panel(
              "Hierarchical clustering",
              shiny::plotOutput(session$ns("hc_module_plot"), height = "480px")
            ),
            if (!is.null(current$modules_kmeans)) {
              bslib::nav_panel(
                "K-means refinement",
                shiny::plotOutput(session$ns("km_module_plot"), height = "480px")
              )
            },
            bslib::nav_panel(
              "Table",
              shiny::div(
                class = "threshold-table-wrap",
                shiny::tableOutput(session$ns("module_table"))
              )
            )
          ),
          bslib::card(
            class = "network-artifacts-card",
            bslib::card_header("Artifacts"),
            shiny::p(shiny::strong("Method"), shiny::br(), current$method),
            shiny::p(
              shiny::strong("TOM"),
              shiny::br(),
              paste(current$tom$method, current$tom$type, sep = " · ")
            ),
            shiny::p(
              shiny::strong("Network object"),
              shiny::br(),
              shiny::code("out/3_networks.rds")
            ),
            shiny::downloadButton(session$ns("download_log"), "Download task log")
          )
        )
      )
    })

    output$hc_module_plot <- shiny::renderPlot({
      shiny::req(report())
      counts <- subset(report()$module_counts, network == "Hierarchical clustering")
      print(network_module_count_plot(counts, "Genes per HC module"))
    }, res = 110)

    output$km_module_plot <- shiny::renderPlot({
      shiny::req(report(), report()$modules_kmeans)
      counts <- subset(report()$module_counts, network == "K-means refinement")
      print(network_module_count_plot(counts, "Genes per k-means module"))
    }, res = 110)

    output$module_table <- shiny::renderTable({
      shiny::req(report())
      report()$module_counts
    }, striped = TRUE, hover = TRUE, spacing = "s")

    output$download_log <- shiny::downloadHandler(
      filename = function() "3_network_construction.log",
      content = function(file) {
        source <- report()$files$log
        shiny::req(source, file.exists(source))
        file.copy(source, file, overwrite = TRUE)
      }
    )
  })
}
