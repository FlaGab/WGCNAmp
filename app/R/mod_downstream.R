downstream_module_choices <- function(module_table, include_grey = FALSE) {
  if (is.null(module_table) || nrow(module_table) == 0L) return(character())
  if (!include_grey) {
    module_table <- module_table[module_table$module != "grey", , drop = FALSE]
  }
  labels <- paste0(module_table$module, " — ", module_table$genes, " genes")
  stats::setNames(module_table$module, labels)
}

downstream_selected <- function(report, parameter, available, default_count = 1L) {
  stored <- if (is.null(report)) NULL else report$parameters[[parameter]]
  selected <- intersect(stored %||% character(), available)
  if (length(selected) > 0L) return(selected)
  utils::head(available, default_count)
}

downstream_stage_heading <- function(step, title, description, ns) {
  shiny::div(
    class = "stage-heading",
    shiny::div(
      shiny::span(class = "stage-number", paste("Step", step)),
      shiny::h2(title),
      shiny::p(class = "stage-description", description)
    ),
    shiny::uiOutput(ns("stage_badge"))
  )
}

downstream_prerequisite_notice <- function(message) {
  shiny::div(
    class = "run-feedback prerequisite-warning",
    shiny::icon("lock"),
    shiny::div(shiny::strong("Previous-stage artifacts are required"), shiny::p(message))
  )
}

hubAnalysisUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    downstream_stage_heading(
      5,
      "Hub genes",
      "Rank genes within selected modules using WGCNA module membership and intramodular soft connectivity.",
      ns
    ),
    shiny::uiOutput(ns("prerequisite_notice")),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        class = "downstream-controls-card",
        fill = FALSE,
        bslib::card_header("Hub selection"),
        shiny::uiOutput(ns("network_selector")),
        shiny::uiOutput(ns("module_selector")),
        shiny::div(
          class = "form-grid form-grid-2",
          shiny::sliderInput(
            ns("min_mm"), "Minimum absolute module membership",
            min = 0, max = 1, value = 0.8, step = 0.05
          ),
          shiny::sliderInput(
            ns("connectivity_quantile"), "Within-module connectivity quantile",
            min = 0.5, max = 0.99, value = 0.9, step = 0.01
          ),
          shiny::uiOutput(ns("trait_selector")),
          shiny::selectInput(
            ns("trait_method"), "Gene-significance correlation",
            choices = c("Pearson" = "pearson", "Spearman" = "spearman")
          )
        ),
        shiny::div(
          class = "selection-guidance",
          shiny::icon("circle-info"),
          shiny::span(
            "Connectivity is calculated only inside selected modules with WGCNA softConnectivity(), avoiding a multi-gigabyte full TOM load."
          )
        )
      ),
      shiny::div(
        bslib::card(
          class = "downstream-summary-card",
          bslib::card_header("Available network"),
          shiny::uiOutput(ns("input_summary"))
        ),
        bslib::card(
          class = "run-card",
          bslib::card_header("Run analysis"),
          shiny::uiOutput(ns("run_control")),
          shiny::p(
            class = "run-note",
            "Larger modules take longer because connectivity is recalculated from expression data."
          )
        ),
        shiny::uiOutput(ns("run_feedback"))
      )
    ),
    shiny::uiOutput(ns("results_panel"))
  )
}

hubAnalysisServer <- function(id, selected_session, on_complete = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    spec <- shiny::reactiveVal(NULL)
    report <- shiny::reactiveVal(NULL)
    load_error <- shiny::reactiveVal(NULL)
    run_error <- shiny::reactiveVal(NULL)
    selected_path <- shiny::reactive({
      selected <- selected_session()
      if (is.null(selected)) NULL else selected$path
    })
    ready <- shiny::reactive({
      path <- selected_path()
      !is.null(path) && file.exists(file.path(path, "out", "1_prepared_data.RData")) &&
        !is.null(module_trait_network_artifact(path))
    })
    shiny::observeEvent(selected_path(), {
      report(load_hub_report(selected_path()))
      run_error(NULL)
      current <- tryCatch(
        downstream_control_spec(selected_path()),
        error = function(error) {
          load_error(conditionMessage(error))
          NULL
        }
      )
      spec(current)
      if (!is.null(current)) load_error(NULL)
    }, ignoreNULL = TRUE)

    output$stage_badge <- shiny::renderUI({
      report()
      if (!ready()) return(status_badge("blocked"))
      status_badge(if (stage_is_complete("step5", selected_path())) "complete" else "ready")
    })
    output$prerequisite_notice <- shiny::renderUI({
      if (ready()) return(NULL)
      downstream_prerequisite_notice("Complete data preparation and network construction first.")
    })
    output$network_selector <- shiny::renderUI({
      current <- spec()
      if (is.null(current)) return(NULL)
      available <- names(current$networks)
      shiny::selectInput(
        session$ns("network_type"), "Network",
        choices = stats::setNames(available, unname(current$networks)),
        selected = downstream_selected(report(), "network_type", available)
      )
    })
    selected_network <- shiny::reactive({
      current <- spec()
      shiny::req(current)
      input$network_type %||% names(current$networks)[[1L]]
    })
    output$module_selector <- shiny::renderUI({
      current <- spec()
      shiny::req(current)
      network <- selected_network()
      choices <- downstream_module_choices(current$module_tables[[network]])
      available <- unname(choices)
      shiny::checkboxGroupInput(
        session$ns("modules"), "Modules",
        choices = choices,
        selected = downstream_selected(report(), "modules", available, 1L)
      )
    })
    output$trait_selector <- shiny::renderUI({
      current <- spec()
      if (is.null(current)) return(NULL)
      choices <- c("None" = "", current$continuous_traits)
      stored <- if (is.null(report())) "" else report()$parameters$trait %||% ""
      shiny::selectInput(
        session$ns("trait"), "Optional continuous trait", choices, selected = stored
      )
    })
    output$input_summary <- shiny::renderUI({
      current <- spec()
      if (is.null(current)) return(shiny::p(class = "text-secondary", load_error()))
      table <- current$module_tables[[selected_network()]]
      shiny::div(
        class = "network-input-grid",
        downstream_metric(current$samples, "samples"),
        downstream_metric(current$genes, "genes"),
        downstream_metric(sum(table$module != "grey"), "non-grey modules"),
        downstream_metric(length(current$continuous_traits), "continuous traits")
      )
    })
    output$run_control <- shiny::renderUI({
      button <- shiny::actionButton(
        session$ns("run"), "Find hub genes",
        icon = shiny::icon("circle-nodes"), class = "btn-primary run-stage-button"
      )
      if (!ready() || is.null(spec())) button$attribs$disabled <- "disabled"
      button
    })
    shiny::observeEvent(input$run, {
      tryCatch({
        modules <- input$modules %||% character()
        if (length(modules) == 0L) stop("Select at least one module.", call. = FALSE)
        run_error(NULL)
        result <- shiny::withProgress(
          message = "Calculating hub-gene scores", value = 0.15,
          {
            result <- run_hub_analysis(
              session_dir = selected_path(),
              network_type = selected_network(),
              modules = modules,
              min_module_membership = input$min_mm,
              connectivity_quantile = input$connectivity_quantile,
              trait = input$trait,
              trait_method = input$trait_method,
              verbose = FALSE
            )
            shiny::incProgress(0.85, detail = "Saving hub tables and plot")
            result
          }
        )
        report(result)
        if (is.function(on_complete)) on_complete()
      }, error = function(error) {
        run_error(conditionMessage(error))
        shiny::showNotification(conditionMessage(error), type = "error", duration = NULL)
      })
    })
    output$run_feedback <- downstream_feedback(report, run_error, "Hub analysis complete")
    output$hub_plot <- shiny::renderPlot({
      shiny::req(report())
      print(report()$plot)
    }, res = 110)
    output$hub_table <- shiny::renderTable({
      shiny::req(report())
      table <- report()$hub_table
      table <- table[table$is_hub, , drop = FALSE]
      utils::head(table, 100L)
    }, digits = 4, striped = TRUE, hover = TRUE, spacing = "s")
    output$results_panel <- shiny::renderUI({
      current <- report()
      if (is.null(current)) return(downstream_empty("No hub-gene report yet."))
      shiny::tagList(
        downstream_results_heading("Hub-gene candidates", current$completed_at),
        shiny::div(
          class = "network-result-grid",
          downstream_metric(current$genes, "module genes scored"),
          downstream_metric(current$hubs, "hub candidates"),
          downstream_metric(length(current$parameters$modules), "modules"),
          downstream_metric(current$parameters$network_type, "network")
        ),
        bslib::navset_card_tab(
          title = "Hub results",
          bslib::nav_panel(
            "Connectivity plot", shiny::plotOutput(session$ns("hub_plot"), height = "600px")
          ),
          bslib::nav_panel(
            "Hub table",
            shiny::div(class = "threshold-table-wrap", shiny::tableOutput(session$ns("hub_table")))
          )
        )
      )
    })
  })
}

enrichmentUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    downstream_stage_heading(
      6,
      "GO enrichment",
      "Test selected WGCNA modules for over-represented Arabidopsis Gene Ontology terms.",
      ns
    ),
    shiny::uiOutput(ns("prerequisite_notice")),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        class = "downstream-controls-card",
        fill = FALSE,
        bslib::card_header("Enrichment settings"),
        shiny::uiOutput(ns("network_selector")),
        shiny::uiOutput(ns("module_selector")),
        shiny::div(
          class = "form-grid form-grid-2",
          shiny::selectInput(
            ns("ontology"), "Ontology",
            choices = c(
              "Biological process" = "BP", "Molecular function" = "MF",
              "Cellular component" = "CC", "All ontologies" = "ALL"
            )
          ),
          shiny::selectInput(
            ns("p_adjust"), "Multiple-testing correction",
            choices = c("Benjamini-Hochberg" = "BH", "Benjamini-Yekutieli" = "BY", "Bonferroni" = "bonferroni")
          ),
          shiny::numericInput(ns("p_cutoff"), "P-value cutoff", 0.05, min = 0, max = 1, step = 0.01),
          shiny::numericInput(ns("q_cutoff"), "Q-value cutoff", 0.2, min = 0, max = 1, step = 0.05),
          shiny::numericInput(ns("min_size"), "Minimum GO gene-set size", 10, min = 1, step = 1),
          shiny::numericInput(ns("show_categories"), "Terms shown per module", 15, min = 1, step = 1)
        ),
        shiny::checkboxInput(ns("simplify"), "Remove semantically redundant GO terms", TRUE),
        shiny::conditionalPanel(
          condition = sprintf("input['%s']", ns("simplify")),
          shiny::sliderInput(
            ns("similarity_cutoff"), "Semantic-similarity cutoff",
            min = 0.4, max = 0.95, value = 0.7, step = 0.05
          )
        )
      ),
      shiny::div(
        bslib::card(
          class = "downstream-summary-card",
          bslib::card_header("Dependencies"),
          shiny::uiOutput(ns("dependency_status"))
        ),
        bslib::card(
          class = "run-card",
          bslib::card_header("Run enrichment"),
          shiny::uiOutput(ns("run_control")),
          shiny::p(class = "run-note", "The full network gene set is used as the enrichment universe.")
        ),
        shiny::uiOutput(ns("run_feedback"))
      )
    ),
    shiny::uiOutput(ns("results_panel"))
  )
}

enrichmentServer <- function(id, selected_session, on_complete = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    spec <- shiny::reactiveVal(NULL)
    report <- shiny::reactiveVal(NULL)
    run_error <- shiny::reactiveVal(NULL)
    selected_path <- shiny::reactive({
      selected <- selected_session()
      if (is.null(selected)) NULL else selected$path
    })
    ready <- shiny::reactive({
      path <- selected_path()
      !is.null(path) && file.exists(file.path(path, "out", "1_prepared_data.RData")) &&
        !is.null(module_trait_network_artifact(path))
    })
    missing_packages <- shiny::reactive({
      packages <- c("clusterProfiler", "org.At.tair.db", "GO.db")
      available <- vapply(packages, function(package) {
        suppressPackageStartupMessages(requireNamespace(package, quietly = TRUE))
      }, logical(1))
      packages[!available]
    })
    shiny::observeEvent(selected_path(), {
      report(load_enrichment_report(selected_path()))
      run_error(NULL)
      spec(tryCatch(downstream_control_spec(selected_path()), error = function(error) NULL))
    }, ignoreNULL = TRUE)
    output$stage_badge <- shiny::renderUI({
      report()
      if (!ready()) return(status_badge("blocked"))
      status_badge(if (stage_is_complete("step6", selected_path())) "complete" else "ready")
    })
    output$prerequisite_notice <- shiny::renderUI({
      if (ready()) return(NULL)
      downstream_prerequisite_notice("Complete network construction before running enrichment.")
    })
    output$network_selector <- shiny::renderUI({
      current <- spec()
      if (is.null(current)) return(NULL)
      available <- names(current$networks)
      shiny::selectInput(
        session$ns("network_type"), "Network",
        choices = stats::setNames(available, unname(current$networks)),
        selected = downstream_selected(report(), "network_type", available)
      )
    })
    selected_network <- shiny::reactive({
      shiny::req(spec())
      input$network_type %||% names(spec()$networks)[[1L]]
    })
    output$module_selector <- shiny::renderUI({
      current <- spec()
      shiny::req(current)
      choices <- downstream_module_choices(current$module_tables[[selected_network()]])
      available <- unname(choices)
      shiny::checkboxGroupInput(
        session$ns("modules"), "Modules",
        choices = choices,
        selected = downstream_selected(report(), "modules", available, min(3L, length(available)))
      )
    })
    output$dependency_status <- shiny::renderUI({
      missing <- missing_packages()
      if (length(missing) == 0L) {
        return(shiny::div(class = "dependency-ready", shiny::icon("circle-check"), "GO annotation packages are available."))
      }
      shiny::div(
        class = "dependency-missing", shiny::icon("triangle-exclamation"),
        paste("Missing R packages:", paste(missing, collapse = ", "))
      )
    })
    output$run_control <- shiny::renderUI({
      button <- shiny::actionButton(
        session$ns("run"), "Run GO enrichment",
        icon = shiny::icon("flask"), class = "btn-primary run-stage-button"
      )
      if (!ready() || is.null(spec()) || length(missing_packages()) > 0L) {
        button$attribs$disabled <- "disabled"
      }
      button
    })
    shiny::observeEvent(input$run, {
      tryCatch({
        modules <- input$modules %||% character()
        if (length(modules) == 0L) stop("Select at least one module.", call. = FALSE)
        run_error(NULL)
        result <- shiny::withProgress(message = "Running GO enrichment", value = 0.1, {
          value <- run_module_enrichment(
            session_dir = selected_path(), network_type = selected_network(),
            modules = modules, ontology = input$ontology,
            p_adjust_method = input$p_adjust, p_value_cutoff = input$p_cutoff,
            q_value_cutoff = input$q_cutoff, min_gene_set_size = input$min_size,
            simplify_terms = isTRUE(input$simplify),
            similarity_cutoff = input$similarity_cutoff %||% 0.7,
            show_categories = input$show_categories, verbose = FALSE
          )
          shiny::incProgress(0.9, detail = "Saving enrichment results")
          value
        })
        report(result)
        if (is.function(on_complete)) on_complete()
      }, error = function(error) {
        run_error(conditionMessage(error))
        shiny::showNotification(conditionMessage(error), type = "error", duration = NULL)
      })
    })
    output$run_feedback <- downstream_feedback(report, run_error, "GO enrichment complete")
    output$enrichment_plot <- shiny::renderPlot({
      shiny::req(report())
      print(report()$plot)
    }, res = 110)
    output$enrichment_table <- shiny::renderTable({
      shiny::req(report())
      table <- report()$results
      columns <- intersect(c("module", "ID", "Description", "GeneRatio", "Count", "p.adjust"), names(table))
      utils::head(table[order(table$p.adjust), columns, drop = FALSE], 100L)
    }, digits = 4, striped = TRUE, hover = TRUE, spacing = "s")
    output$results_panel <- shiny::renderUI({
      current <- report()
      if (is.null(current)) return(downstream_empty("No enrichment report yet."))
      shiny::tagList(
        downstream_results_heading("GO enrichment", current$completed_at),
        shiny::div(
          class = "network-result-grid",
          downstream_metric(current$terms, "enriched terms"),
          downstream_metric(current$modules_with_terms, "modules with terms"),
          downstream_metric(length(current$parameters$modules), "modules tested"),
          downstream_metric(current$parameters$ontology, "ontology")
        ),
        bslib::navset_card_tab(
          title = "Enrichment results",
          bslib::nav_panel("Dot plot", shiny::plotOutput(session$ns("enrichment_plot"), height = "680px")),
          bslib::nav_panel("Terms", shiny::div(class = "threshold-table-wrap", shiny::tableOutput(session$ns("enrichment_table"))))
        )
      )
    })
  })
}

geneVisualizationUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    downstream_stage_heading(
      7,
      "Gene visualization",
      "Explore normalized expression, module summaries, or eigengene patterns without editing gene lists in a notebook.",
      ns
    ),
    shiny::uiOutput(ns("prerequisite_notice")),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        class = "downstream-controls-card",
        fill = FALSE,
        bslib::card_header("Visualization settings"),
        shiny::uiOutput(ns("network_selector")),
        shiny::uiOutput(ns("module_selector")),
        shiny::selectInput(
          ns("plot_type"), "View",
          choices = c(
            "Selected gene profiles" = "gene_profiles",
            "Module expression summaries" = "module_summary",
            "Module eigengene heatmap" = "eigengene_heatmap"
          )
        ),
        shiny::uiOutput(ns("gene_selector")),
        shiny::uiOutput(ns("group_selector")),
        shiny::checkboxInput(ns("scale_rows"), "Scale each gene to a z-score", TRUE)
      ),
      shiny::div(
        bslib::card(
          class = "downstream-summary-card",
          bslib::card_header("Current selection"),
          shiny::uiOutput(ns("selection_summary"))
        ),
        bslib::card(
          class = "run-card",
          bslib::card_header("Create visualization"),
          shiny::uiOutput(ns("run_control")),
          shiny::p(class = "run-note", "Plots and their underlying long-format data are saved with the session.")
        ),
        shiny::uiOutput(ns("run_feedback"))
      )
    ),
    shiny::uiOutput(ns("results_panel"))
  )
}

geneVisualizationServer <- function(id, selected_session, on_complete = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    spec <- shiny::reactiveVal(NULL)
    report <- shiny::reactiveVal(NULL)
    run_error <- shiny::reactiveVal(NULL)
    selected_path <- shiny::reactive({
      selected <- selected_session()
      if (is.null(selected)) NULL else selected$path
    })
    ready <- shiny::reactive({
      path <- selected_path()
      !is.null(path) && file.exists(file.path(path, "out", "1_prepared_data.RData")) &&
        !is.null(module_trait_network_artifact(path))
    })
    shiny::observeEvent(selected_path(), {
      report(load_visualization_report(selected_path()))
      run_error(NULL)
      spec(tryCatch(downstream_control_spec(selected_path()), error = function(error) NULL))
    }, ignoreNULL = TRUE)
    output$stage_badge <- shiny::renderUI({
      report()
      if (!ready()) return(status_badge("blocked"))
      status_badge(if (stage_is_complete("step7", selected_path())) "complete" else "ready")
    })
    output$prerequisite_notice <- shiny::renderUI({
      if (ready()) return(NULL)
      downstream_prerequisite_notice("Complete data preparation and network construction first.")
    })
    output$network_selector <- shiny::renderUI({
      current <- spec()
      if (is.null(current)) return(NULL)
      available <- names(current$networks)
      shiny::selectInput(
        session$ns("network_type"), "Network",
        choices = stats::setNames(available, unname(current$networks)),
        selected = downstream_selected(report(), "network_type", available)
      )
    })
    selected_network <- shiny::reactive({
      shiny::req(spec())
      input$network_type %||% names(spec()$networks)[[1L]]
    })
    output$module_selector <- shiny::renderUI({
      current <- spec()
      shiny::req(current)
      choices <- downstream_module_choices(current$module_tables[[selected_network()]])
      available <- unname(choices)
      shiny::checkboxGroupInput(
        session$ns("modules"), "Modules",
        choices = choices,
        selected = downstream_selected(report(), "modules", available, 1L)
      )
    })
    selected_modules <- shiny::reactive(input$modules %||% character())
    available_genes <- shiny::reactive({
      current <- spec()
      if (is.null(current) || length(selected_modules()) == 0L) return(character())
      colors <- current$colors[[selected_network()]]
      names(colors)[colors %in% selected_modules()]
    })
    output$gene_selector <- shiny::renderUI({
      if (!identical(input$plot_type %||% "gene_profiles", "gene_profiles")) return(NULL)
      genes <- available_genes()
      stored <- if (is.null(report())) character() else report()$parameters$genes
      selected <- intersect(stored, genes)
      if (length(selected) == 0L) selected <- utils::head(genes, 12L)
      shiny::selectizeInput(
        session$ns("genes"), "Genes (maximum 100)",
        choices = genes, selected = selected, multiple = TRUE,
        options = list(maxItems = 100, plugins = list("remove_button"))
      )
    })
    output$group_selector <- shiny::renderUI({
      current <- spec()
      if (is.null(current) || !identical(input$plot_type %||% "gene_profiles", "gene_profiles")) return(NULL)
      stored <- if (is.null(report())) "" else report()$parameters$group_column %||% ""
      shiny::selectInput(
        session$ns("group_column"), "Color samples by",
        choices = c("No grouping" = "", current$metadata_columns), selected = stored
      )
    })
    output$selection_summary <- shiny::renderUI({
      current <- spec()
      if (is.null(current)) return(shiny::p(class = "text-secondary", "Session inputs are not ready."))
      shiny::div(
        class = "network-input-grid",
        downstream_metric(length(selected_modules()), "modules selected"),
        downstream_metric(length(available_genes()), "genes available"),
        downstream_metric(current$samples, "samples"),
        downstream_metric(if (current$normalized_counts_available) "Yes" else "No", "normalized counts")
      )
    })
    output$run_control <- shiny::renderUI({
      button <- shiny::actionButton(
        session$ns("run"), "Create plot",
        icon = shiny::icon("chart-line"), class = "btn-primary run-stage-button"
      )
      if (!ready() || is.null(spec()) || !isTRUE(spec()$normalized_counts_available)) {
        button$attribs$disabled <- "disabled"
      }
      button
    })
    shiny::observeEvent(input$run, {
      tryCatch({
        modules <- selected_modules()
        if (length(modules) == 0L) stop("Select at least one module.", call. = FALSE)
        genes <- input$genes %||% character()
        if (identical(input$plot_type, "gene_profiles") && length(genes) == 0L) {
          stop("Select at least one gene.", call. = FALSE)
        }
        run_error(NULL)
        result <- run_gene_visualization(
          session_dir = selected_path(), network_type = selected_network(),
          modules = modules, genes = genes, plot_type = input$plot_type,
          group_column = input$group_column, scale_rows = isTRUE(input$scale_rows),
          verbose = FALSE
        )
        report(result)
        if (is.function(on_complete)) on_complete()
      }, error = function(error) {
        run_error(conditionMessage(error))
        shiny::showNotification(conditionMessage(error), type = "error", duration = NULL)
      })
    })
    output$run_feedback <- downstream_feedback(report, run_error, "Visualization created")
    output$visualization_plot <- shiny::renderPlot({
      shiny::req(report())
      print(report()$plot)
    }, res = 110)
    output$visualization_table <- shiny::renderTable({
      shiny::req(report())
      utils::head(report()$data, 100L)
    }, digits = 4, striped = TRUE, hover = TRUE, spacing = "s")
    output$results_panel <- shiny::renderUI({
      current <- report()
      if (is.null(current)) return(downstream_empty("No saved visualization yet."))
      shiny::tagList(
        downstream_results_heading("Expression visualization", current$completed_at),
        bslib::navset_card_tab(
          title = "Visualization result",
          bslib::nav_panel("Plot", shiny::plotOutput(session$ns("visualization_plot"), height = "680px")),
          bslib::nav_panel("Data", shiny::div(class = "threshold-table-wrap", shiny::tableOutput(session$ns("visualization_table"))))
        )
      )
    })
  })
}

downstream_metric <- function(value, label) {
  shiny::div(
    class = "network-input-metric",
    shiny::strong(if (is.null(value)) "—" else format(value, big.mark = ",")),
    shiny::span(label)
  )
}

downstream_feedback <- function(report, error, success_title) {
  shiny::renderUI({
    message <- error()
    if (!is.null(message)) {
      return(shiny::div(
        class = "run-feedback run-feedback-error",
        shiny::icon("circle-xmark"),
        shiny::div(shiny::strong("Analysis failed"), shiny::p(message))
      ))
    }
    if (is.null(report())) return(NULL)
    shiny::div(
      class = "run-feedback run-feedback-success",
      shiny::icon("circle-check"),
      shiny::div(shiny::strong(success_title), shiny::p("Results are saved in this session."))
    )
  })
}

downstream_empty <- function(message) {
  shiny::div(
    class = "empty-results",
    shiny::icon("chart-line"),
    shiny::div(shiny::strong(message), shiny::p("Configure this stage and run it to create results."))
  )
}

downstream_results_heading <- function(title, completed_at) {
  shiny::div(
    class = "results-heading",
    shiny::div(shiny::span(class = "eyebrow", "Latest result"), shiny::h2(title)),
    shiny::span(class = "result-time", format(completed_at, "%Y-%m-%d %H:%M"))
  )
}
