required_packages <- c("shiny", "bslib")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_packages) > 0L) {
  stop(
    "Install the GUI dependencies before launching: ",
    paste(missing_packages, collapse = ", ")
  )
}

project_root <- wgcnamp_project_root()
max_upload_mb <- suppressWarnings(as.numeric(
  Sys.getenv("WGCNAMP_MAX_UPLOAD_MB", unset = "512")
))
if (is.na(max_upload_mb) || max_upload_mb <= 0) {
  max_upload_mb <- 512
}
options(shiny.maxRequestSize = max_upload_mb * 1024^2)

if (requireNamespace("future", quietly = TRUE) &&
    requireNamespace("promises", quietly = TRUE)) {
  background_workers <- suppressWarnings(as.integer(
    Sys.getenv("WGCNAMP_BACKGROUND_WORKERS", unset = "1")
  ))
  if (is.na(background_workers) || background_workers < 1L) {
    background_workers <- 1L
  }
  # I() prevents future from replacing a one-worker multisession plan with a
  # sequential plan; network construction must remain outside Shiny's process.
  # Stage 3 intentionally opens and closes PDF devices inside the worker.
  options(future.devices.onMisuse = "ignore")
  future::plan(future::multisession, workers = I(background_workers))
}

backend_files <- c(
  "src/importing_functions_v2.R",
  "src/quality_check.R",
  "src/filtering_functions.R",
  "src/normalization.R",
  "src/data_preparation_pipeline.R",
  "src/soft_threshold_pipeline.R",
  "src/network_construction_pipeline.R",
  "src/mod-trait_corr.R",
  "src/module_trait_pipeline.R",
  "src/downstream_pipeline.R"
)
backend_environment <- environment(dataPreparationServer)
for (backend_file in backend_files) {
  sys.source(file.path(project_root, backend_file), envir = backend_environment)
}

ui <- bslib::page_sidebar(
  window_title = "WGCNAmp",
  title = shiny::div(
    class = "brand",
    shiny::span(class = "brand-mark", "W"),
    shiny::div(
      shiny::div(class = "brand-title", "WGCNAmp"),
      shiny::div(class = "brand-subtitle", "Modular network analysis")
    )
  ),
  theme = bslib::bs_theme(
    version = 5,
    bg = "#f5f7fb",
    fg = "#172033",
    primary = "#176b87",
    secondary = "#667085"
  ),
  fillable = TRUE,
  sidebar = bslib::sidebar(
    width = 330,
    open = "always",
    sessionSelectorUI("session_selector"),
    shiny::hr(),
    shiny::div(class = "eyebrow", "Pipeline"),
    shiny::uiOutput("stage_navigation"),
    shiny::div(
      class = "sidebar-footer",
      shiny::icon("flask"),
      shiny::span("R / WGCNA workflow")
    )
  ),
  shiny::tags$head(
    shiny::tags$link(rel = "stylesheet", type = "text/css", href = "app.css")
  ),
  shiny::div(
    class = "app-content",
    shiny::uiOutput("session_header"),
    shiny::uiOutput("stage_content")
  )
)

server <- function(input, output, session) {
  selector <- sessionSelectorServer("session_selector", project_root)
  pipeline_revision <- shiny::reactiveVal(0L)

  stage_completed <- function(stage_id) {
    function(result = NULL, session_path = NULL) {
      path <- session_path
      if (is.null(path)) {
        selected <- shiny::isolate(selector$selected())
        path <- if (is.null(selected)) NULL else selected$path
      }
      if (!is.null(path)) record_stage_completion(path, stage_id, result)
      pipeline_revision(shiny::isolate(pipeline_revision()) + 1L)
    }
  }

  dataPreparationServer(
    "data_preparation",
    selected_session = selector$selected,
    on_complete = stage_completed("step1")
  )
  softThresholdServer(
    "soft_threshold",
    selected_session = selector$selected,
    on_complete = stage_completed("step2")
  )
  networkConstructionServer(
    "network_construction",
    selected_session = selector$selected,
    project_root = project_root,
    on_complete = stage_completed("step3")
  )
  moduleTraitServer(
    "module_trait",
    selected_session = selector$selected,
    on_complete = stage_completed("step4")
  )
  hubAnalysisServer(
    "hub_analysis",
    selected_session = selector$selected,
    on_complete = stage_completed("step5")
  )
  enrichmentServer(
    "enrichment",
    selected_session = selector$selected,
    on_complete = stage_completed("step6")
  )
  geneVisualizationServer(
    "gene_visualization",
    selected_session = selector$selected,
    on_complete = stage_completed("step7")
  )

  stage_table <- shiny::reactive({
    selected <- selector$selected()
    selector$revision()
    pipeline_revision()
    if (is.null(selected)) {
      return(NULL)
    }
    session_stage_status(selected$path)
  })

  shiny::observeEvent(selector$selected(), {
    selected_session <- selector$selected()
    stages <- stage_table()
    if (is.null(stages) || is.null(selected_session)) {
      return()
    }
    saved_stage <- read_session_state(
      selected_session$path, migrate = TRUE
    )$ui$active_stage
    incomplete <- which(stages$status != "complete")
    next_stage <- if (!is.null(saved_stage) && saved_stage %in% stages$id) {
      saved_stage
    } else if (length(incomplete) > 0L) {
      stages$id[[incomplete[[1L]]]]
    } else {
      stages$id[[nrow(stages)]]
    }
    shiny::updateRadioButtons(session, "active_stage", selected = next_stage)
  }, ignoreNULL = TRUE)

  shiny::observeEvent(input$active_stage, {
    selected <- shiny::isolate(selector$selected())
    if (!is.null(selected) && input$active_stage %in% stage_definitions()$id) {
      update_session_ui_state(selected$path, active_stage = input$active_stage)
    }
  }, ignoreInit = TRUE)

  output$stage_navigation <- shiny::renderUI({
    stages <- stage_table()
    if (is.null(stages)) {
      return(shiny::p(class = "sidebar-empty", "Select or create a session to view its pipeline."))
    }

    choice_labels <- paste(
      vapply(stages$status, status_icon, character(1)),
      stages$number,
      stages$label
    )
    choices <- stats::setNames(stages$id, choice_labels)
    selected_stage <- shiny::isolate(input$active_stage)
    if (is.null(selected_stage) || !selected_stage %in% stages$id) {
      incomplete <- which(stages$status != "complete")
      selected_stage <- if (length(incomplete) > 0L) {
        stages$id[[incomplete[[1L]]]]
      } else {
        stages$id[[nrow(stages)]]
      }
    }

    shiny::div(
      class = "stage-navigation",
      shiny::radioButtons(
        "active_stage",
        label = NULL,
        choices = choices,
        selected = selected_stage
      )
    )
  })

  output$session_header <- shiny::renderUI({
    selected <- selector$selected()
    stages <- stage_table()
    if (is.null(selected) || is.null(stages)) {
      return(shiny::div(
        class = "welcome-panel",
        shiny::div(class = "eyebrow", "WGCNA analysis workspace"),
        shiny::h1("Choose an analysis session"),
        shiny::p(
          "Select an existing session in the sidebar or create a new one. The GUI will expose each pipeline stage without requiring notebook edits."
        )
      ))
    }

    completed <- sum(stages$status == "complete")
    stale <- sum(stages$status == "stale")
    stale_stages <- stages[stages$status == "stale", , drop = FALSE]
    shiny::tagList(
      shiny::div(
        class = "session-header",
        shiny::div(
          shiny::div(class = "eyebrow", "Current analysis"),
          shiny::h1(selected$name),
          shiny::p(class = "session-path", selected$path)
        ),
        shiny::div(
          class = "progress-summary",
          shiny::strong(paste0(completed, " / ", nrow(stages))),
          shiny::span(
            if (stale > 0L) {
              paste("stages complete ·", stale, "stale")
            } else {
              "stages complete"
            }
          )
        )
      ),
      if (stale > 0L) {
        shiny::div(
          class = "stale-summary",
          shiny::icon("triangle-exclamation"),
          shiny::div(
            shiny::strong("Some results need to be rerun"),
            shiny::p(
              paste(stale_stages$label, collapse = ", "),
              shiny::br(),
              stale_stages$stale_reason[[1L]]
            )
          )
        )
      }
    )
  })

  output$stage_content <- shiny::renderUI({
    selected <- selector$selected()
    if (is.null(selected)) {
      return(NULL)
    }

    definitions <- stage_definitions()
    active <- input$active_stage %||% definitions$id[[1L]]
    if (!active %in% definitions$id) {
      active <- definitions$id[[1L]]
    }
    stage <- definitions[match(active, definitions$id), , drop = FALSE]

    if (identical(active, "step1")) {
      return(dataPreparationUI("data_preparation"))
    }
    if (identical(active, "step2")) {
      return(softThresholdUI("soft_threshold"))
    }
    if (identical(active, "step3")) {
      return(networkConstructionUI("network_construction"))
    }
    if (identical(active, "step4")) {
      return(moduleTraitUI("module_trait"))
    }
    if (identical(active, "step5")) {
      return(hubAnalysisUI("hub_analysis"))
    }
    if (identical(active, "step6")) {
      return(enrichmentUI("enrichment"))
    }
    if (identical(active, "step7")) {
      return(geneVisualizationUI("gene_visualization"))
    }

    bslib::layout_columns(
      col_widths = c(8, 4),
      bslib::card(
        class = "stage-card",
        bslib::card_header(
          shiny::div(
            class = "stage-card-header",
            shiny::div(
              shiny::span(class = "stage-number", paste0("Step ", stage$number)),
              shiny::h2(stage$label)
            ),
            status_badge(stage$status)
          )
        ),
        shiny::p(class = "stage-description", stage$description),
        shiny::div(
          class = "shell-notice",
          shiny::icon("screwdriver-wrench"),
          shiny::div(
            shiny::strong("Stage controls are coming next"),
            shiny::p(
              "This shell already tracks the session and pipeline state. The next implementation pass will move this stage's notebook parameters and Run action into this panel."
            )
          )
        ),
        bslib::card_footer(
          shiny::span(class = "text-secondary", "Current reference: "),
          shiny::code(stage$notebook)
        )
      ),
      bslib::card(
        class = "outputs-card",
        bslib::card_header("Recent artifacts"),
        if (length(outputs) == 0L) {
          shiny::p(class = "empty-artifacts", "No outputs have been generated for this session yet.")
        } else {
          shiny::tags$ul(
            class = "artifact-list",
            lapply(outputs, function(path) {
              shiny::tags$li(
                shiny::icon(if (grepl("[.]pdf$|[.]png$", path, ignore.case = TRUE)) "chart-line" else "file"),
                shiny::span(path)
              )
            })
          )
        }
      )
    )
  })
}

shiny::shinyApp(ui = ui, server = server)
