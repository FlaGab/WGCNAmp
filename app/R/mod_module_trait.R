moduleTraitUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::div(
      class = "stage-heading",
      shiny::div(
        shiny::span(class = "stage-number", "Step 4"),
        shiny::h2("Module–trait correlation"),
        shiny::p(
          class = "stage-description",
          "Choose traits and categorical contrasts from the current session. Controls update automatically when prepared data or networks change."
        )
      ),
      shiny::uiOutput(ns("stage_badge"))
    ),
    shiny::uiOutput(ns("prerequisite_notice")),
    bslib::layout_columns(
      col_widths = c(7, 5),
      shiny::div(
        bslib::card(
          class = "module-trait-controls-card",
          fill = FALSE,
          bslib::card_header("Analysis settings"),
          shiny::h3("Networks"),
          shiny::uiOutput(ns("network_selector")),
          shiny::hr(),
          shiny::h3("Continuous traits"),
          shiny::uiOutput(ns("continuous_empty")),
          shiny::uiOutput(ns("continuous_selector")),
          shiny::selectInput(
            ns("continuous_method"),
            "Correlation method",
            choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
            selected = "pearson"
          ),
          shiny::hr(),
          shiny::h3("Categorical metadata"),
          shiny::uiOutput(ns("categorical_selector")),
          shiny::uiOutput(ns("skipped_variables")),
          shiny::uiOutput(ns("reference_controls")),
          shiny::div(
            class = "form-grid form-grid-2 module-trait-options",
            shiny::checkboxInput(
              ns("include_pairwise"), "Pairwise contrasts", value = TRUE
            ),
            shiny::checkboxInput(
              ns("include_vs_all"), "Level versus all", value = TRUE
            ),
            shiny::numericInput(
              ns("min_count"),
              "Minimum observations per group",
              value = 1,
              min = 1,
              step = 1
            ),
            shiny::checkboxInput(
              ns("drop_redundant"), "Drop redundant contrasts", value = TRUE
            )
          ),
          shiny::hr(),
          shiny::h3("Categorical linear models"),
          shiny::checkboxInput(
            ns("fit_categorical_models"),
            "Fit contrast models with optional covariate adjustment",
            value = TRUE
          ),
          shiny::uiOutput(ns("model_covariate_selector")),
          shiny::uiOutput(ns("model_formula_preview"))
        )
      ),
      shiny::div(
        bslib::card(
          class = "module-trait-input-card",
          bslib::card_header("Available session data"),
          shiny::uiOutput(ns("input_summary"))
        ),
        bslib::card(
          class = "module-trait-contrast-card",
          bslib::card_header("Generated contrasts"),
          shiny::uiOutput(ns("contrast_status")),
          shiny::uiOutput(ns("contrast_selector"))
        ),
        bslib::card(
          class = "run-card module-trait-run-card",
          bslib::card_header("Run analysis"),
          shiny::uiOutput(ns("run_control")),
          shiny::p(
            class = "run-note",
            "Running again replaces the stage-4 report, correlation and model tables, and selected heatmaps."
          )
        ),
        shiny::uiOutput(ns("run_feedback"))
      )
    ),
    shiny::uiOutput(ns("results_panel"))
  )
}

module_trait_reference_id <- function(variable) {
  bytes <- as.integer(charToRaw(enc2utf8(variable)))
  paste0("reference_", paste(sprintf("%02x", bytes), collapse = ""))
}

module_trait_choices <- function(values, counts, suffix) {
  if (length(values) == 0L) return(character())
  labels <- paste0(values, " — ", unname(counts[values]), " ", suffix)
  stats::setNames(values, labels)
}

moduleTraitServer <- function(id, selected_session, on_complete = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    report <- shiny::reactiveVal(NULL)
    spec <- shiny::reactiveVal(NULL)
    spec_error <- shiny::reactiveVal(NULL)
    run_error <- shiny::reactiveVal(NULL)

    selected_path <- shiny::reactive({
      selected <- selected_session()
      if (is.null(selected)) NULL else selected$path
    })

    prerequisites_ready <- shiny::reactive({
      path <- selected_path()
      !is.null(path) && stage_dependencies_ready(path, "step4")
    })

    shiny::observeEvent(selected_path(), {
      path <- selected_path()
      current <- load_module_trait_report(path)
      report(current)
      run_error(NULL)
      current_spec <- tryCatch(
        module_trait_control_spec(path),
        error = function(error) {
          spec_error(conditionMessage(error))
          NULL
        }
      )
      spec(current_spec)
      if (is.null(current_spec)) return()
      spec_error(NULL)

      parameters <- if (is.null(current)) list() else current$parameters
      if (!is.null(parameters$continuous_method)) {
        shiny::updateSelectInput(
          session, "continuous_method", selected = parameters$continuous_method
        )
      }
      if (!is.null(parameters$include_pairwise)) {
        shiny::updateCheckboxInput(
          session, "include_pairwise", value = parameters$include_pairwise
        )
        shiny::updateCheckboxInput(
          session, "include_vs_all", value = parameters$include_vs_all
        )
        shiny::updateNumericInput(
          session, "min_count", value = parameters$min_count
        )
        shiny::updateCheckboxInput(
          session, "drop_redundant", value = parameters$drop_redundant
        )
      }
      if (!is.null(parameters$fit_categorical_models)) {
        shiny::updateCheckboxInput(
          session,
          "fit_categorical_models",
          value = parameters$fit_categorical_models
        )
      }
    }, ignoreNULL = TRUE)

    selected_from_report <- function(parameter, available) {
      current <- report()
      stored <- if (is.null(current)) NULL else current$parameters[[parameter]]
      intersect(stored %||% available, available)
    }

    output$network_selector <- shiny::renderUI({
      current <- spec()
      if (is.null(current)) return(NULL)
      available <- names(current$networks)
      shiny::checkboxGroupInput(
        session$ns("network_types"),
        label = NULL,
        choices = stats::setNames(available, unname(current$networks)),
        selected = selected_from_report("network_types", available)
      )
    })

    output$continuous_selector <- shiny::renderUI({
      current <- spec()
      if (is.null(current) || length(current$continuous) == 0L) return(NULL)
      shiny::checkboxGroupInput(
        session$ns("continuous_traits"),
        label = NULL,
        choices = module_trait_choices(
          current$continuous,
          current$continuous_counts,
          "distinct values"
        ),
        selected = selected_from_report("continuous_traits", current$continuous)
      )
    })

    output$categorical_selector <- shiny::renderUI({
      current <- spec()
      if (is.null(current) || length(current$categorical) == 0L) return(NULL)
      shiny::checkboxGroupInput(
        session$ns("categorical_variables"),
        label = NULL,
        choices = module_trait_choices(
          current$categorical,
          current$categorical_counts,
          "levels"
        ),
        selected = selected_from_report(
          "categorical_variables", current$categorical
        )
      )
    })

    available_model_covariates <- shiny::reactive({
      current <- spec()
      if (is.null(current) || nrow(current$covariates) == 0L) {
        return(data.frame(token = character(), label = character()))
      }
      selected_primary <- input$categorical_variables %||% character()
      current$covariates[
        !current$covariates$column %in% selected_primary,
        ,
        drop = FALSE
      ]
    })

    output$model_covariate_selector <- shiny::renderUI({
      if (!isTRUE(input$fit_categorical_models)) return(NULL)
      current <- spec()
      selected_categorical <- input$categorical_variables %||% character()
      if (is.null(current) || length(selected_categorical) == 0L) {
        return(shiny::p(
          class = "field-help",
          "Select categorical metadata to configure its regression model."
        ))
      }
      available <- available_model_covariates()
      stored <- if (is.null(report())) character() else {
        report()$parameters$model_covariates %||% character()
      }
      shiny::selectizeInput(
        session$ns("model_covariates"),
        "Adjustment covariates",
        choices = stats::setNames(available$token, available$label),
        selected = intersect(stored, available$token),
        multiple = TRUE,
        options = list(
          placeholder = "None — fit an unadjusted contrast model",
          plugins = list("remove_button")
        )
      )
    })

    output$model_formula_preview <- shiny::renderUI({
      if (!isTRUE(input$fit_categorical_models) ||
          length(input$categorical_variables %||% character()) == 0L) {
        return(NULL)
      }
      available <- available_model_covariates()
      selected <- input$model_covariates %||% character()
      labels <- available$label[match(selected, available$token)]
      labels <- labels[!is.na(labels)]
      formula <- paste(
        "module eigengene ~ binary contrast",
        if (length(labels) > 0L) {
          paste("+", paste(labels, collapse = " + "))
        } else {
          ""
        }
      )
      shiny::div(
        class = "selection-guidance module-model-formula",
        shiny::icon("calculator"),
        shiny::span(
          shiny::strong("Model: "),
          formula,
          shiny::br(),
          "Each contrast is fitted separately. Reported q-values use BH adjustment across all module–contrast tests within each network."
        )
      )
    })

    shiny::observeEvent(
      list(
        input$network_types,
        input$continuous_traits,
        input$categorical_variables,
        input$model_covariates,
        input$fit_categorical_models
      ),
      run_error(NULL),
      ignoreInit = TRUE
    )

    output$stage_badge <- shiny::renderUI({
      report()
      path <- selected_path()
      if (is.null(path)) return(status_badge("blocked"))
      session_stage_badge(path, "step4")
    })

    output$prerequisite_notice <- shiny::renderUI({
      if (prerequisites_ready()) return(NULL)
      shiny::div(
        class = "run-feedback prerequisite-warning",
        shiny::icon("lock"),
        shiny::div(
          shiny::strong("Prepared data and a network are required"),
          shiny::p("Complete steps 1 and 3 before configuring module–trait correlations.")
        )
      )
    })

    output$continuous_empty <- shiny::renderUI({
      current <- spec()
      if (is.null(current) || length(current$continuous) > 0L) return(NULL)
      shiny::p(
        class = "field-help module-trait-empty-note",
        "No informative numeric columns are available in the traits table."
      )
    })

    output$skipped_variables <- shiny::renderUI({
      current <- spec()
      if (is.null(current) || length(current$skipped_categorical) == 0L) {
        return(NULL)
      }
      shiny::div(
        class = "selection-guidance module-trait-skipped",
        shiny::icon("circle-info"),
        shiny::span(
          "Ignored because only one level is observed: ",
          paste(current$skipped_categorical, collapse = ", "),
          "."
        )
      )
    })

    output$reference_controls <- shiny::renderUI({
      current <- spec()
      selected <- input$categorical_variables %||% character()
      if (is.null(current) || length(selected) == 0L) return(NULL)
      stored <- if (is.null(report())) NULL else report()$parameters$reference_levels
      shiny::tagList(
        shiny::h4("Reference levels"),
        shiny::div(
          class = "form-grid form-grid-2 reference-level-grid",
          lapply(selected, function(variable) {
            levels <- current$categorical_levels[[variable]]
            selected_level <- levels[[1L]]
            if (!is.null(stored) && variable %in% names(stored) &&
                stored[[variable]] %in% levels) {
              selected_level <- stored[[variable]]
            }
            shiny::selectInput(
              session$ns(module_trait_reference_id(variable)),
              variable,
              choices = levels,
              selected = selected_level
            )
          })
        )
      )
    })

    selected_reference_levels <- shiny::reactive({
      current <- spec()
      selected <- input$categorical_variables %||% character()
      if (is.null(current) || length(selected) == 0L) return(NULL)
      references <- vapply(selected, function(variable) {
        input[[module_trait_reference_id(variable)]] %||%
          current$categorical_levels[[variable]][[1L]]
      }, character(1))
      stats::setNames(references, selected)
    })

    contrast_preview <- shiny::reactive({
      current <- spec()
      selected <- input$categorical_variables %||% character()
      if (is.null(current) || length(selected) == 0L) {
        return(list(matrix = NULL, error = NULL))
      }
      if (!isTRUE(input$include_pairwise) && !isTRUE(input$include_vs_all)) {
        return(list(
          matrix = NULL,
          error = "Enable pairwise or level-versus-all contrasts."
        ))
      }
      tryCatch(
        list(
          matrix = binary_contrast_matrix(
            data = current$metadata,
            variables = selected,
            reference_levels = selected_reference_levels(),
            vsAll = isTRUE(input$include_vs_all),
            includePairwise = isTRUE(input$include_pairwise),
            minCount = as.integer(input$min_count %||% 1L),
            drop_redundant = isTRUE(input$drop_redundant),
            verbose = FALSE
          ),
          error = NULL
        ),
        error = function(error) list(matrix = NULL, error = conditionMessage(error))
      )
    })

    output$contrast_status <- shiny::renderUI({
      preview <- contrast_preview()
      if (!is.null(preview$error)) {
        return(shiny::div(
          class = "candidate-summary candidate-summary-invalid",
          shiny::icon("triangle-exclamation"),
          shiny::span(preview$error)
        ))
      }
      if (is.null(preview$matrix)) {
        return(shiny::p(
          class = "text-secondary",
          "Select categorical variables to generate contrasts."
        ))
      }
      shiny::div(
        class = "candidate-summary candidate-summary-valid",
        shiny::strong(ncol(preview$matrix)),
        shiny::span("informative contrasts generated")
      )
    })

    output$contrast_selector <- shiny::renderUI({
      preview <- contrast_preview()
      if (is.null(preview$matrix)) return(NULL)
      available <- colnames(preview$matrix)
      stored <- if (is.null(report())) character() else {
        report()$parameters$selected_contrasts
      }
      selected <- intersect(stored, available)
      if (length(selected) == 0L) selected <- available
      shiny::checkboxGroupInput(
        session$ns("selected_contrasts"),
        "Contrasts to analyse",
        choices = available,
        selected = selected
      )
    })

    output$input_summary <- shiny::renderUI({
      current <- spec()
      if (is.null(current)) {
        return(shiny::p(
          class = "text-secondary",
          spec_error() %||% "Session inputs are not ready."
        ))
      }
      metric <- function(value, label) {
        shiny::div(
          class = "network-input-metric",
          shiny::strong(value),
          shiny::span(label)
        )
      }
      shiny::div(
        class = "network-input-grid",
        metric(current$samples, "samples"),
        metric(length(current$networks), "networks"),
        metric(length(current$continuous), "continuous traits"),
        metric(length(current$categorical), "categorical variables")
      )
    })

    output$run_control <- shiny::renderUI({
      button <- shiny::actionButton(
        session$ns("run"),
        "Run module–trait analysis",
        icon = shiny::icon("table-cells"),
        class = "btn-primary run-stage-button"
      )
      if (!prerequisites_ready() || is.null(spec())) {
        button$attribs$disabled <- "disabled"
      }
      button
    })

    shiny::observeEvent(input$run, {
      tryCatch({
        if (!prerequisites_ready()) {
          stop("Complete steps 1 and 3 before running this analysis.", call. = FALSE)
        }
        continuous <- input$continuous_traits %||% character()
        categorical <- input$categorical_variables %||% character()
        networks <- input$network_types %||% character()
        if (length(continuous) == 0L && length(categorical) == 0L) {
          stop("Select at least one continuous or categorical trait.", call. = FALSE)
        }
        if (length(networks) == 0L) {
          stop("Select at least one network.", call. = FALSE)
        }
        contrasts <- NULL
        if (length(categorical) > 0L) {
          preview <- contrast_preview()
          if (!is.null(preview$error)) stop(preview$error, call. = FALSE)
          contrasts <- intersect(
            input$selected_contrasts %||% character(),
            colnames(preview$matrix)
          )
          if (length(contrasts) == 0L) {
            stop("Select at least one generated contrast.", call. = FALSE)
          }
        }

        run_error(NULL)
        shiny::withProgress(message = "Calculating module–trait correlations", value = 0, {
          shiny::incProgress(0.2, detail = "Validating session data")
          result <- run_module_trait_analysis(
            session_dir = selected_path(),
            continuous_traits = continuous,
            categorical_variables = categorical,
            reference_levels = selected_reference_levels(),
            include_pairwise = isTRUE(input$include_pairwise),
            include_vs_all = isTRUE(input$include_vs_all),
            min_count = input$min_count,
            drop_redundant = isTRUE(input$drop_redundant),
            selected_contrasts = contrasts,
            network_types = networks,
            continuous_method = input$continuous_method,
            fit_categorical_models = isTRUE(input$fit_categorical_models),
            model_covariates = intersect(
              input$model_covariates %||% character(),
              available_model_covariates()$token
            ),
            create_plots = TRUE,
            save_artifacts = TRUE,
            verbose = FALSE
          )
          shiny::incProgress(0.8, detail = "Saving tables and heatmaps")
        })
        report(result)
        if (is.function(on_complete)) on_complete(result, selected_path())
        shiny::showNotification(
          "Module–trait analysis completed.", type = "message"
        )
      }, error = function(error) {
        run_error(conditionMessage(error))
        shiny::showNotification(
          paste("Module–trait analysis failed:", conditionMessage(error)),
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
          shiny::div(shiny::strong("Analysis failed"), shiny::p(error))
        ))
      }
      if (is.null(report())) return(NULL)
      shiny::div(
        class = "run-feedback run-feedback-success",
        shiny::icon("circle-check"),
        shiny::div(
          shiny::strong("Analysis complete"),
          shiny::p("Correlation tables and heatmaps are available below.")
        )
      )
    })

    module_trait_plot <- function(network, analysis) {
      shiny::renderPlot({
        current <- report()
        shiny::req(current, current$results[[network]][[analysis]]$plot)
        print(current$results[[network]][[analysis]]$plot)
      }, res = 110)
    }
    output$hc_continuous_plot <- module_trait_plot("hc", "continuous")
    output$hc_categorical_plot <- module_trait_plot("hc", "categorical")
    output$hc_categorical_model_plot <- module_trait_plot("hc", "categorical_model")
    output$kmeans_continuous_plot <- module_trait_plot("kmeans", "continuous")
    output$kmeans_categorical_plot <- module_trait_plot("kmeans", "categorical")
    output$kmeans_categorical_model_plot <- module_trait_plot(
      "kmeans", "categorical_model"
    )

    adjusted_model_table <- shiny::reactive({
      current <- report()
      if (is.null(current)) return(NULL)
      tables <- lapply(names(current$results), function(network) {
        model <- current$results[[network]]$categorical_model
        if (is.null(model)) return(NULL)
        data.frame(
          network = network,
          model$table,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      })
      tables <- Filter(Negate(is.null), tables)
      if (length(tables) == 0L) NULL else do.call(rbind, tables)
    })

    output$adjusted_model_table <- shiny::renderTable({
      table <- adjusted_model_table()
      shiny::req(table)
      columns <- intersect(
        c(
          "network", "module", "contrast", "estimate", "std_error",
          "conf_low", "conf_high", "p_value", "fdr", "n",
          "adjusted_r_squared", "omitted_covariates", "aliased_terms", "status"
        ),
        names(table)
      )
      utils::head(table[, columns, drop = FALSE], 250L)
    }, digits = 4, striped = TRUE, hover = TRUE, spacing = "s")

    output$results_panel <- shiny::renderUI({
      current <- report()
      if (is.null(current)) {
        return(shiny::div(
          class = "empty-results",
          shiny::icon("table-cells"),
          shiny::div(
            shiny::strong("No module–trait report yet"),
            shiny::p("Select traits and contrasts, then run the analysis.")
          )
        ))
      }
      panels <- list()
      for (network in names(current$results)) {
        for (analysis in c("continuous", "categorical", "categorical_model")) {
          if (is.null(current$results[[network]][[analysis]])) next
          output_id <- paste(network, analysis, "plot", sep = "_")
          network_label <- if (identical(network, "hc")) "HC" else "K-means"
          analysis_label <- switch(
            analysis,
            continuous = "continuous correlation",
            categorical = "categorical correlation",
            categorical_model = "categorical linear model"
          )
          panels[[length(panels) + 1L]] <- bslib::nav_panel(
            paste(network_label, analysis_label, sep = " · "),
            shiny::plotOutput(session$ns(output_id), height = "620px")
          )
        }
      }
      plots <- do.call(
        bslib::navset_card_tab,
        c(list(title = "Correlation heatmaps"), panels)
      )
      metric <- function(value, label) {
        shiny::div(
          class = "network-result-metric",
          shiny::strong(value),
          shiny::span(label)
        )
      }
      shiny::tagList(
        shiny::div(
          class = "results-heading",
          shiny::div(
            shiny::span(class = "eyebrow", "Latest result"),
            shiny::h2("Module–trait associations")
          ),
          shiny::span(
            class = "result-time",
            format(current$completed_at, "%Y-%m-%d %H:%M")
          )
        ),
        shiny::div(
          class = "network-result-grid",
          metric(current$samples, "samples"),
          metric(length(current$modules), "networks"),
          metric(length(current$parameters$continuous_traits), "continuous traits"),
          metric(current$contrasts, "categorical contrasts")
        ),
        plots,
        if (!is.null(adjusted_model_table())) {
          bslib::card(
            class = "module-model-table-card",
            bslib::card_header("Categorical linear-model estimates"),
            shiny::p(
              class = "field-help",
              "The coefficient is the model-estimated eigengene difference for contrast 1 minus contrast 0, adjusted when covariates are selected. The table preview is limited to 250 rows."
            ),
            shiny::div(
              class = "threshold-table-wrap module-model-table-wrap",
              shiny::tableOutput(session$ns("adjusted_model_table"))
            ),
            bslib::card_footer(
              shiny::downloadButton(
                session$ns("download_adjusted_models"),
                "Download complete model table"
              )
            )
          )
        },
        if (!is.null(current$files$binary_contrasts)) {
          shiny::div(
            class = "plot-downloads",
            shiny::downloadButton(
              session$ns("download_contrasts"), "Download binary contrast matrix"
            )
          )
        }
      )
    })

    output$download_contrasts <- shiny::downloadHandler(
      filename = function() "module_trait_binary_contrasts.tsv",
      content = function(file) {
        source <- report()$files$binary_contrasts
        shiny::req(source, file.exists(source))
        file.copy(source, file, overwrite = TRUE)
      }
    )

    output$download_adjusted_models <- shiny::downloadHandler(
      filename = function() "4_adjusted_categorical_models.tsv",
      content = function(file) {
        source <- report()$files$adjusted_models
        if (!is.null(source) && file.exists(source)) {
          file.copy(source, file, overwrite = TRUE)
        } else {
          utils::write.table(
            adjusted_model_table(),
            file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            na = "NA"
          )
        }
      }
    )

    list(
      spec = spec,
      contrast_preview = contrast_preview,
      report = report
    )
  })
}
