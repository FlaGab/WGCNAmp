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

format_storage_size <- function(bytes) {
  if (!is.numeric(bytes) || length(bytes) != 1L || !is.finite(bytes) || bytes <= 0) {
    return("0 B")
  }
  units <- c("B", "KB", "MB", "GB", "TB")
  unit <- min(floor(log(bytes, base = 1024)) + 1L, length(units))
  value <- bytes / 1024^(unit - 1L)
  paste0(format(round(value, if (value < 10) 1L else 0L), trim = TRUE), " ", units[[unit]])
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
      "Rank hub genes, inspect TOM connectivity, and export weighted networks for Cytoscape.",
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
        shiny::h3("Connectivity definition"),
        shiny::radioButtons(
          ns("connectivity_mode"),
          label = NULL,
          choices = c(
            "Scalable adjacency — recompute selected modules" = "adjacency",
            "TOM-based — reproduce the notebook definition" = "tom"
          ),
          selected = "adjacency"
        ),
        shiny::uiOutput(ns("connectivity_guidance")),
        shiny::h3("TOM visualization and Cytoscape export"),
        shiny::uiOutput(ns("tom_output_controls")),
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
            "TOM outputs load stored matrices even when scalable adjacency is used for hub ranking."
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
      !is.null(path) && stage_dependencies_ready(path, "step5")
    })
    shiny::observeEvent(selected_path(), {
      current_report <- load_hub_report(selected_path())
      report(current_report)
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
      if (!is.null(current_report)) {
        parameters <- current_report$parameters
        shiny::updateRadioButtons(
          session,
          "connectivity_mode",
          selected = parameters$connectivity_mode %||% "adjacency"
        )
        shiny::updateSliderInput(
          session, "min_mm", value = parameters$min_module_membership
        )
        shiny::updateSliderInput(
          session,
          "connectivity_quantile",
          value = parameters$connectivity_quantile
        )
        shiny::updateSelectInput(
          session, "trait_method", selected = parameters$trait_method
        )
      }
    }, ignoreNULL = TRUE)

    output$stage_badge <- shiny::renderUI({
      report()
      path <- selected_path()
      if (is.null(path)) return(status_badge("blocked"))
      session_stage_badge(path, "step5")
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
    output$connectivity_guidance <- shiny::renderUI({
      current <- spec()
      if (identical(input$connectivity_mode %||% "adjacency", "adjacency")) {
        return(shiny::div(
          class = "selection-guidance",
          shiny::icon("bolt"),
          shiny::span(
            shiny::strong("Scalable adjacency mode. "),
            "WGCNA softConnectivity() recomputes connectivity only for the selected modules. It is faster and avoids loading stored TOM matrices, but its kWithin values are not identical to the notebook's TOM-based values."
          )
        ))
      }
      if (is.null(current)) return(NULL)
      tom <- current$tom
      if (!isTRUE(tom$available)) {
        return(shiny::div(
          class = "selection-guidance network-memory-warning",
          shiny::icon("triangle-exclamation"),
          shiny::span(
            shiny::strong("TOM mode is unavailable. "),
            tom$message %||% "No compatible TOM artifact was found."
          )
        ))
      }
      shiny::div(
        class = "selection-guidance network-memory-warning",
        shiny::icon("hard-drive"),
        shiny::span(
          shiny::strong("Notebook-compatible TOM mode. "),
          if (identical(tom$method, "blockwise")) {
            paste0(
              "Stored TOM blocks are loaded sequentially, and only blocks containing selected modules are analysed. The manifest contains ",
              tom$blocks, " block(s), using ", format_storage_size(tom$bytes), " on disk. "
            )
          } else {
            paste0(
              "The full stored TOM matrix is loaded, using ",
              format_storage_size(tom$bytes), " on disk. "
            )
          },
          "Memory use can exceed the file size and the analysis can take substantially longer."
        )
      )
    })
    output$tom_output_controls <- shiny::renderUI({
      current <- spec()
      if (is.null(current)) return(NULL)
      if (!isTRUE(current$tom$available)) {
        return(shiny::div(
          class = "selection-guidance network-memory-warning",
          shiny::icon("triangle-exclamation"),
          shiny::span(
            shiny::strong("TOM outputs are unavailable. "),
            current$tom$message %||% "No compatible TOM artifact was found."
          )
        ))
      }
      parameters <- if (is.null(report())) list() else report()$parameters
      shiny::tagList(
        shiny::checkboxInput(
          session$ns("create_tom_outputs"),
          "Create a TOM heatmap and Cytoscape node/edge tables",
          parameters$create_tom_outputs %||% FALSE
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s']", session$ns("create_tom_outputs")),
          shiny::div(
            class = "form-grid form-grid-2",
            shiny::numericInput(
              session$ns("tom_max_genes"),
              "Maximum genes in TOM outputs",
              value = parameters$tom_max_genes %||% 250,
              min = 2, max = 1000, step = 10
            ),
            shiny::sliderInput(
              session$ns("tom_edge_quantile"),
              "Keep TOM edges at or above quantile",
              min = 0.5, max = 0.999,
              value = parameters$tom_edge_quantile %||% 0.95,
              step = 0.01
            )
          ),
          shiny::div(
            class = "selection-guidance",
            shiny::icon("circle-info"),
            shiny::span(
              "Hub candidates are selected first, followed by genes with the strongest module membership, while cycling across modules."
            )
          )
        )
      )
    })
    tom_outputs_requested <- shiny::reactive({
      isTRUE(input$create_tom_outputs) &&
        !is.null(spec()) && isTRUE(spec()$tom$available)
    })
    output$run_control <- shiny::renderUI({
      button <- shiny::actionButton(
        session$ns("run"), "Find hub genes",
        icon = shiny::icon("circle-nodes"), class = "btn-primary run-stage-button"
      )
      tom_unavailable <- identical(
        input$connectivity_mode %||% "adjacency", "tom"
      ) && (is.null(spec()) || !isTRUE(spec()$tom$available))
      if (!ready() || is.null(spec()) || tom_unavailable) {
        button$attribs$disabled <- "disabled"
      }
      button
    })
    shiny::observeEvent(input$run, {
      tryCatch({
        modules <- input$modules %||% character()
        if (length(modules) == 0L) stop("Select at least one module.", call. = FALSE)
        connectivity_mode <- input$connectivity_mode %||% "adjacency"
        if (identical(connectivity_mode, "tom") && !isTRUE(spec()$tom$available)) {
          stop(spec()$tom$message %||% "No compatible TOM artifact was found.", call. = FALSE)
        }
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
              connectivity_mode = connectivity_mode,
              trait = input$trait,
              trait_method = input$trait_method,
              create_tom_outputs = tom_outputs_requested(),
              tom_max_genes = input$tom_max_genes %||% 250,
              tom_edge_quantile = input$tom_edge_quantile %||% 0.95,
              verbose = FALSE
            )
            shiny::incProgress(0.85, detail = "Saving hub, TOM, and Cytoscape outputs")
            result
          }
        )
        report(result)
        if (is.function(on_complete)) on_complete(result, selected_path())
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
    output$tom_plot <- shiny::renderPlot({
      shiny::req(report(), report()$tom_network)
      print(report()$tom_network$plot)
    }, res = 110)
    output$cytoscape_edges <- shiny::renderTable({
      shiny::req(report(), report()$tom_network)
      utils::head(report()$tom_network$edge_table, 100L)
    }, digits = 4, striped = TRUE, hover = TRUE, spacing = "s")
    output$download_cytoscape_edges <- shiny::downloadHandler(
      filename = function() {
        shiny::req(report(), report()$files$cytoscape_edges)
        basename(report()$files$cytoscape_edges)
      },
      content = function(file) {
        shiny::req(report(), report()$files$cytoscape_edges)
        file.copy(report()$files$cytoscape_edges, file, overwrite = TRUE)
      }
    )
    output$download_cytoscape_nodes <- shiny::downloadHandler(
      filename = function() {
        shiny::req(report(), report()$files$cytoscape_nodes)
        basename(report()$files$cytoscape_nodes)
      },
      content = function(file) {
        shiny::req(report(), report()$files$cytoscape_nodes)
        file.copy(report()$files$cytoscape_nodes, file, overwrite = TRUE)
      }
    )
    output$results_panel <- shiny::renderUI({
      current <- report()
      if (is.null(current)) return(downstream_empty("No hub-gene report yet."))
      result_tabs <- list(
        bslib::nav_panel(
          "Hub score plot", shiny::plotOutput(session$ns("hub_plot"), height = "600px")
        ),
        bslib::nav_panel(
          "Hub table",
          shiny::div(class = "threshold-table-wrap", shiny::tableOutput(session$ns("hub_table")))
        )
      )
      if (!is.null(current$tom_network)) {
        result_tabs <- c(result_tabs, list(
          bslib::nav_panel(
            "TOM connectivity",
            shiny::plotOutput(session$ns("tom_plot"), height = "720px")
          ),
          bslib::nav_panel(
            "Cytoscape export",
            shiny::div(
              class = "selection-guidance",
              shiny::icon("circle-nodes"),
              shiny::div(
                shiny::strong(paste(
                  format(current$tom_network$genes, big.mark = ","), "nodes and",
                  format(current$tom_network$edges, big.mark = ","), "edges"
                )),
                shiny::p(paste0(
                  "Edges retain TOM weights at or above the ",
                  format(100 * current$tom_network$edge_quantile, trim = TRUE),
                  "th percentile (threshold ",
                  format(current$tom_network$threshold, digits = 4), ").",
                  if (identical(current$tom_network$block_method, "blockwise")) {
                    " Cross-block TOM values are unavailable, so cross-block edges are omitted."
                  } else {
                    ""
                  }
                ))
              )
            ),
            shiny::div(
              class = "form-grid form-grid-2",
              shiny::downloadButton(session$ns("download_cytoscape_nodes"), "Download node table"),
              shiny::downloadButton(session$ns("download_cytoscape_edges"), "Download edge table")
            ),
            shiny::h3("Edge preview"),
            shiny::div(
              class = "threshold-table-wrap",
              shiny::tableOutput(session$ns("cytoscape_edges"))
            )
          )
        ))
      }
      shiny::tagList(
        downstream_results_heading("Hub-gene candidates", current$completed_at),
        shiny::div(
          class = "network-result-grid",
          downstream_metric(current$genes, "module genes scored"),
          downstream_metric(current$hubs, "hub candidates"),
          downstream_metric(length(current$parameters$modules), "modules"),
          downstream_metric(
            if (identical(current$parameters$connectivity_mode, "tom")) {
              "TOM-based"
            } else {
              "Scalable adjacency"
            },
            "connectivity"
          )
        ),
        do.call(
          bslib::navset_card_tab,
          c(list(title = "Hub results"), result_tabs)
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
      "Functional enrichment",
      "Test selected WGCNA modules for over-represented Arabidopsis GO terms or KEGG pathways.",
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
        shiny::radioButtons(
          ns("database"), "Enrichment database",
          choices = c("Gene Ontology (GO)" = "GO", "KEGG pathways" = "KEGG"),
          selected = "GO", inline = TRUE
        ),
        shiny::uiOutput(ns("database_controls")),
        shiny::div(
          class = "form-grid form-grid-2",
          shiny::selectInput(
            ns("p_adjust"), "Multiple-testing correction",
            choices = c("Benjamini-Hochberg" = "BH", "Benjamini-Yekutieli" = "BY", "Bonferroni" = "bonferroni")
          ),
          shiny::numericInput(ns("p_cutoff"), "P-value cutoff", 0.05, min = 0, max = 1, step = 0.01),
          shiny::numericInput(ns("q_cutoff"), "Q-value cutoff", 0.2, min = 0, max = 1, step = 0.05),
          shiny::numericInput(ns("min_size"), "Minimum gene-set size", 10, min = 1, step = 1),
          shiny::numericInput(ns("max_size"), "Maximum gene-set size", 500, min = 1, step = 10),
          shiny::numericInput(ns("show_categories"), "Terms shown per module", 15, min = 1, step = 1)
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
          shiny::p(class = "run-note", "The full network gene set is used as the enrichment background.")
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
      !is.null(path) && stage_dependencies_ready(path, "step6")
    })
    selected_database <- shiny::reactive({
      database <- toupper(input$database %||% "GO")
      if (!database %in% c("GO", "KEGG")) "GO" else database
    })
    missing_packages <- shiny::reactive({
      packages <- if (identical(selected_database(), "GO")) {
        c("clusterProfiler", "org.At.tair.db", "GO.db")
      } else {
        c("clusterProfiler", "org.At.tair.db", "AnnotationDbi")
      }
      available <- vapply(packages, function(package) {
        suppressPackageStartupMessages(requireNamespace(package, quietly = TRUE))
      }, logical(1))
      packages[!available]
    })
    shiny::observeEvent(selected_path(), {
      current_report <- load_enrichment_report(selected_path())
      report(current_report)
      run_error(NULL)
      spec(tryCatch(downstream_control_spec(selected_path()), error = function(error) NULL))
      parameters <- if (is.null(current_report)) list() else current_report$parameters
      shiny::updateRadioButtons(
        session, "database", selected = parameters$database %||% "GO"
      )
      shiny::updateSelectInput(session, "p_adjust", selected = parameters$p_adjust_method %||% "BH")
      shiny::updateNumericInput(session, "p_cutoff", value = parameters$p_value_cutoff %||% 0.05)
      shiny::updateNumericInput(session, "q_cutoff", value = parameters$q_value_cutoff %||% 0.2)
      shiny::updateNumericInput(session, "min_size", value = parameters$min_gene_set_size %||% 10)
      shiny::updateNumericInput(session, "max_size", value = parameters$max_gene_set_size %||% 500)
      shiny::updateNumericInput(session, "show_categories", value = parameters$show_categories %||% 15)
    }, ignoreNULL = TRUE)
    output$stage_badge <- shiny::renderUI({
      report()
      path <- selected_path()
      if (is.null(path)) return(status_badge("blocked"))
      session_stage_badge(path, "step6")
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
    output$database_controls <- shiny::renderUI({
      current <- report()
      parameters <- if (is.null(current)) list() else current$parameters
      if (identical(selected_database(), "KEGG")) {
        return(shiny::div(
          class = "selection-guidance",
          shiny::icon("circle-info"),
          shiny::div(
            shiny::strong("Arabidopsis thaliana pathways (organism code: ath)"),
            shiny::p(
              "Network TAIR identifiers are converted to Entrez Gene IDs before analysis. ",
              "Unmapped genes are excluded and counted in the result. KEGG queries require internet access."
            )
          )
        ))
      }
      shiny::tagList(
        shiny::selectInput(
          session$ns("ontology"), "Ontology",
          choices = c(
            "Biological process" = "BP", "Molecular function" = "MF",
            "Cellular component" = "CC", "All ontologies" = "ALL"
          ),
          selected = parameters$ontology %||% "BP"
        ),
        shiny::checkboxInput(
          session$ns("simplify"),
          "Remove semantically redundant GO terms",
          parameters$simplify_terms %||% TRUE
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s']", session$ns("simplify")),
          shiny::sliderInput(
            session$ns("similarity_cutoff"), "Semantic-similarity cutoff",
            min = 0.4, max = 0.95,
            value = parameters$similarity_cutoff %||% 0.7, step = 0.05
          )
        )
      )
    })
    output$dependency_status <- shiny::renderUI({
      missing <- missing_packages()
      if (length(missing) == 0L) {
        label <- if (identical(selected_database(), "GO")) {
          "GO annotation packages are available."
        } else {
          "KEGG and identifier-mapping packages are available."
        }
        return(shiny::div(class = "dependency-ready", shiny::icon("circle-check"), label))
      }
      shiny::div(
        class = "dependency-missing", shiny::icon("triangle-exclamation"),
        paste("Missing R packages:", paste(missing, collapse = ", "))
      )
    })
    output$run_control <- shiny::renderUI({
      button <- shiny::actionButton(
        session$ns("run"), paste("Run", selected_database(), "enrichment"),
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
        if (input$max_size < input$min_size) {
          stop("Maximum gene-set size must be at least the minimum gene-set size.", call. = FALSE)
        }
        database <- selected_database()
        run_error(NULL)
        result <- shiny::withProgress(message = paste("Running", database, "enrichment"), value = 0.1, {
          value <- run_module_enrichment(
            session_dir = selected_path(), network_type = selected_network(),
            modules = modules, database = database,
            ontology = input$ontology %||% "BP",
            p_adjust_method = input$p_adjust, p_value_cutoff = input$p_cutoff,
            q_value_cutoff = input$q_cutoff, min_gene_set_size = input$min_size,
            max_gene_set_size = input$max_size,
            simplify_terms = isTRUE(input$simplify),
            similarity_cutoff = input$similarity_cutoff %||% 0.7,
            show_categories = input$show_categories, verbose = FALSE
          )
          shiny::incProgress(0.9, detail = "Saving enrichment results")
          value
        })
        report(result)
        if (is.function(on_complete)) on_complete(result, selected_path())
      }, error = function(error) {
        run_error(conditionMessage(error))
        shiny::showNotification(conditionMessage(error), type = "error", duration = NULL)
      })
    })
    output$run_feedback <- downstream_feedback(report, run_error, "Enrichment complete")
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
      database <- current$parameters$database %||% "GO"
      fourth_metric <- if (identical(database, "GO")) {
        downstream_metric(current$parameters$ontology, "ontology")
      } else {
        downstream_metric(database, "database")
      }
      mapping_note <- NULL
      if (identical(database, "KEGG") && !is.null(current$mapping)) {
        mapping_note <- shiny::div(
          class = "selection-guidance",
          shiny::icon("circle-info"),
          shiny::div(
            shiny::strong("TAIR to Entrez mapping"),
            shiny::p(
              paste0(
                format(current$mapping$mapped_genes, big.mark = ","), " of ",
                format(current$mapping$input_genes, big.mark = ","),
                " network genes mapped; ",
                format(current$mapping$unmapped_genes, big.mark = ","), " were excluded."
              )
            )
          )
        )
      }
      shiny::tagList(
        downstream_results_heading(paste(database, "enrichment"), current$completed_at),
        shiny::div(
          class = "network-result-grid",
          downstream_metric(current$terms, "enriched terms"),
          downstream_metric(current$modules_with_terms, "modules with terms"),
          downstream_metric(length(current$parameters$modules), "modules tested"),
          fourth_metric
        ),
        mapping_note,
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
      !is.null(path) && stage_dependencies_ready(path, "step7")
    })
    shiny::observeEvent(selected_path(), {
      report(load_visualization_report(selected_path()))
      run_error(NULL)
      spec(tryCatch(downstream_control_spec(selected_path()), error = function(error) NULL))
    }, ignoreNULL = TRUE)
    output$stage_badge <- shiny::renderUI({
      report()
      path <- selected_path()
      if (is.null(path)) return(status_badge("blocked"))
      session_stage_badge(path, "step7")
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
        if (is.function(on_complete)) on_complete(result, selected_path())
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
