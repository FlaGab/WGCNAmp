softThresholdUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::div(
      class = "stage-heading",
      shiny::div(
        shiny::span(class = "stage-number", "Step 2"),
        shiny::h2("Soft-threshold selection"),
        shiny::p(
          class = "stage-description",
          "Evaluate scale-free topology and mean connectivity across candidate powers, then save one explicit power for network construction."
        )
      ),
      shiny::uiOutput(ns("stage_badge"))
    ),
    shiny::uiOutput(ns("prerequisite_notice")),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        class = "threshold-controls-card",
        bslib::card_header("Candidate powers and network"),
        shiny::selectInput(
          ns("network_type"),
          "Network type",
          choices = c(
            "Signed" = "signed",
            "Unsigned" = "unsigned",
            "Signed hybrid" = "signed hybrid"
          ),
          selected = "signed"
        ),
        shiny::textInput(
          ns("candidate_powers"),
          "Candidate powers",
          value = "1:20, 22:50:2"
        ),
        shiny::p(
          class = "field-help",
          "Use comma-separated values or ranges. A third range value is the step, for example 22:50:2."
        ),
        shiny::sliderInput(
          ns("fit_cutoff"),
          "Scale-free fit target",
          min = 0.50,
          max = 0.95,
          value = 0.80,
          step = 0.01
        ),
        shiny::uiOutput(ns("candidate_summary"))
      ),
      shiny::div(
        bslib::card(
          class = "threshold-controls-card",
          bslib::card_header("Selection policy"),
          shiny::radioButtons(
            ns("selection_mode"),
            label = NULL,
            choices = c(
              "First power reaching the fit target" = "first_above_cutoff",
              "WGCNA automatic estimate" = "automatic",
              "Manual power" = "manual"
            ),
            selected = "first_above_cutoff"
          ),
          shiny::uiOutput(ns("manual_power_control")),
          shiny::div(
            class = "selection-guidance",
            shiny::icon("circle-info"),
            shiny::span(
              "The fit-target policy is reproducible and fails clearly when no evaluated power reaches the target."
            )
          )
        ),
        bslib::card(
          class = "run-card",
          bslib::card_header("Run selection"),
          shiny::uiOutput(ns("dependency_status")),
          shiny::uiOutput(ns("run_control")),
          shiny::p(
            class = "run-note",
            "Running again replaces the step-2 power, report, fit table, and PDF. Existing network results are not deleted, but should be rebuilt when the selected power changes."
          )
        )
      )
    ),
    shiny::uiOutput(ns("run_feedback")),
    shiny::uiOutput(ns("results_panel"))
  )
}

soft_threshold_packages <- function() {
  c("WGCNA", "ggplot2", "gridExtra")
}

missing_soft_threshold_packages <- function() {
  packages <- soft_threshold_packages()
  packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
}

load_soft_threshold_report <- function(session_path) {
  path <- file.path(session_path, "out", "2_soft_threshold_report.rds")
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(error) NULL)
}

softThresholdServer <- function(id, selected_session, on_complete = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    report <- shiny::reactiveVal(NULL)
    run_error <- shiny::reactiveVal(NULL)
    run_completed <- shiny::reactiveVal(FALSE)

    selected_path <- shiny::reactive({
      selected <- selected_session()
      if (is.null(selected)) NULL else selected$path
    })

    prepared_data_exists <- shiny::reactive({
      path <- selected_path()
      !is.null(path) && stage_dependencies_ready(path, "step2")
    })

    candidate_powers <- shiny::reactive({
      parse_soft_threshold_powers(input$candidate_powers)
    })

    shiny::observeEvent(selected_path(), {
      current <- load_soft_threshold_report(selected_path())
      report(current)
      run_error(NULL)
      run_completed(FALSE)
      if (!is.null(current)) {
        shiny::updateSelectInput(
          session, "network_type", selected = current$network_type
        )
        shiny::updateTextInput(
          session,
          "candidate_powers",
          value = paste(current$powers, collapse = ", ")
        )
        shiny::updateSliderInput(
          session, "fit_cutoff", value = current$fit_cutoff
        )
        shiny::updateRadioButtons(
          session, "selection_mode", selected = current$selection_mode
        )
      }
    }, ignoreNULL = TRUE)

    output$stage_badge <- shiny::renderUI({
      report()
      path <- selected_path()
      if (is.null(path)) return(status_badge("blocked"))
      session_stage_badge(path, "step2")
    })

    output$prerequisite_notice <- shiny::renderUI({
      if (prepared_data_exists()) {
        return(NULL)
      }
      shiny::div(
        class = "run-feedback prerequisite-warning",
        shiny::icon("lock"),
        shiny::div(
          shiny::strong("Step 1 is required"),
          shiny::p(
            "Complete data preparation first. Soft-threshold estimation uses the VST expression matrix saved in out/1_prepared_data.RData."
          )
        )
      )
    })

    output$candidate_summary <- shiny::renderUI({
      tryCatch({
        powers <- candidate_powers()
        shiny::div(
          class = "candidate-summary candidate-summary-valid",
          shiny::strong(length(powers)),
          shiny::span(
            paste0("powers from ", min(powers), " to ", max(powers))
          )
        )
      }, error = function(error) {
        shiny::div(
          class = "candidate-summary candidate-summary-invalid",
          shiny::icon("triangle-exclamation"),
          shiny::span(conditionMessage(error))
        )
      })
    })

    output$manual_power_control <- shiny::renderUI({
      if (!identical(input$selection_mode, "manual")) {
        return(NULL)
      }
      choices <- tryCatch(candidate_powers(), error = function(error) integer())
      current <- report()
      selected <- if (!is.null(current) &&
                      identical(current$selection_mode, "manual") &&
                      current$power_threshold %in% choices) {
        current$power_threshold
      } else if (length(choices) > 0L) {
        choices[[1L]]
      } else {
        NULL
      }
      shiny::selectInput(
        session$ns("manual_power"),
        "Manual power",
        choices = choices,
        selected = selected
      )
    })

    output$dependency_status <- shiny::renderUI({
      missing <- missing_soft_threshold_packages()
      if (length(missing) == 0L) {
        shiny::div(
          class = "dependency-ready",
          shiny::icon("circle-check"),
          "Soft-threshold dependencies are available."
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
      button <- shiny::actionButton(
        session$ns("run"),
        "Evaluate and save power",
        icon = shiny::icon("chart-line"),
        class = "btn-primary run-stage-button"
      )
      if (!prepared_data_exists()) {
        button$attribs$disabled <- "disabled"
      }
      button
    })

    shiny::observeEvent(input$run, {
      run_error(NULL)
      run_completed(FALSE)
      missing <- missing_soft_threshold_packages()
      if (length(missing) > 0L) {
        run_error(paste("Install the required R packages:", paste(missing, collapse = ", ")))
        return()
      }

      tryCatch({
        if (!prepared_data_exists()) {
          stop("Complete step 1 before estimating a soft threshold.", call. = FALSE)
        }
        powers <- candidate_powers()
        manual_power <- if (identical(input$selection_mode, "manual")) {
          suppressWarnings(as.numeric(input$manual_power))
        } else {
          NULL
        }

        result <- shiny::withProgress(
          message = "Estimating the soft threshold",
          detail = paste("Evaluating", length(powers), "candidate powers…"),
          value = 0.1,
          {
            value <- run_soft_threshold_selection(
              session_dir = selected_path(),
              network_type = input$network_type,
              powers = powers,
              fit_cutoff = input$fit_cutoff,
              selection_mode = input$selection_mode,
              manual_power = manual_power,
              save_artifact = TRUE,
              create_plot = TRUE,
              verbose = TRUE
            )
            shiny::incProgress(0.9, detail = "Saving the selected power and report…")
            value
          }
        )
        report(load_soft_threshold_report(selected_path()))
        run_completed(TRUE)
        if (is.function(on_complete)) {
          on_complete(result, selected_path())
        }
        shiny::showNotification(
          paste("Selected soft-threshold power", result$power_threshold),
          type = "message",
          duration = 6
        )
      }, error = function(error) {
        run_error(conditionMessage(error))
        shiny::showNotification(
          paste("Soft-threshold selection failed:", conditionMessage(error)),
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
          shiny::div(shiny::strong("Power selection did not complete"), shiny::p(error))
        ))
      }
      if (isTRUE(run_completed())) {
        return(shiny::div(
          class = "run-feedback run-feedback-success",
          shiny::icon("circle-check"),
          shiny::div(
            shiny::strong("Soft-threshold power saved"),
            shiny::p("Review the topology and connectivity diagnostics below before building the network.")
          )
        ))
      }
      NULL
    })

    output$results_panel <- shiny::renderUI({
      current <- report()
      if (is.null(current)) {
        return(shiny::div(
          class = "empty-results",
          shiny::icon("chart-line"),
          shiny::div(
            shiny::strong("No soft-threshold report yet"),
            shiny::p("Evaluate the candidate powers to see scale-free fit and mean connectivity.")
          )
        ))
      }

      automatic <- if (length(current$automatic_power) == 1L &&
                       is.finite(current$automatic_power)) {
        current$automatic_power
      } else {
        "Not available"
      }
      shiny::tagList(
        shiny::div(
          class = "results-heading",
          shiny::div(
            shiny::span(class = "eyebrow", "Latest result"),
            shiny::h2("Soft-threshold diagnostics")
          ),
          if (!is.null(current$completed_at)) {
            shiny::span(
              class = "result-time",
              format(current$completed_at, "%Y-%m-%d %H:%M")
            )
          }
        ),
        bslib::layout_columns(
          col_widths = c(4, 8),
          bslib::card(
            class = "threshold-summary-card",
            bslib::card_header("Selected power"),
            shiny::div(class = "selected-power", current$power_threshold),
            shiny::tags$dl(
              shiny::tags$dt("Network type"),
              shiny::tags$dd(current$network_type),
              shiny::tags$dt("Selection policy"),
              shiny::tags$dd(gsub("_", " ", current$selection_mode)),
              shiny::tags$dt("Fit target"),
              shiny::tags$dd(format(current$fit_cutoff, digits = 2)),
              shiny::tags$dt("WGCNA estimate"),
              shiny::tags$dd(automatic)
            ),
            shiny::p(
              class = "artifact-path",
              shiny::strong("Network input"),
              shiny::code("out/2_pt.Rdata")
            )
          ),
          bslib::navset_card_tab(
            title = "Diagnostics",
            footer = shiny::downloadButton(
              session$ns("download_plot"), "Download diagnostic PDF"
            ),
            bslib::nav_panel(
              "Scale independence",
              shiny::plotOutput(session$ns("fit_plot"), height = "440px")
            ),
            bslib::nav_panel(
              "Mean connectivity",
              shiny::plotOutput(session$ns("connectivity_plot"), height = "440px")
            ),
            bslib::nav_panel(
              "Fit table",
              shiny::div(
                class = "threshold-table-wrap",
                shiny::tableOutput(session$ns("fit_table"))
              )
            )
          )
        )
      )
    })

    output$fit_plot <- shiny::renderPlot({
      shiny::req(report(), report()$plots$scale_independence)
      print(report()$plots$scale_independence)
    }, res = 110)

    output$connectivity_plot <- shiny::renderPlot({
      shiny::req(report(), report()$plots$mean_connectivity)
      print(report()$plots$mean_connectivity)
    }, res = 110)

    output$fit_table <- shiny::renderTable({
      shiny::req(report())
      fit <- report()$fit_indices
      columns <- intersect(
        c("Power", "model_fit", "mean.k.", "median.k.", "max.k."),
        colnames(fit)
      )
      fit[, columns, drop = FALSE]
    }, digits = 3, striped = TRUE, hover = TRUE, spacing = "s")

    output$download_plot <- shiny::downloadHandler(
      filename = function() "sft.pdf",
      content = function(file) {
        source <- report()$files$plot
        shiny::req(source, file.exists(source))
        file.copy(source, file, overwrite = TRUE)
      }
    )
  })
}
