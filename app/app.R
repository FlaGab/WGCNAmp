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

status_icon <- function(status) {
  switch(status, complete = "✓", ready = "●", blocked = "○", "○")
}

status_label <- function(status) {
  switch(status, complete = "Complete", ready = "Ready", blocked = "Waiting", status)
}

status_badge <- function(status) {
  shiny::span(
    class = paste("status-badge", paste0("status-", status)),
    status_label(status)
  )
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

  stage_table <- shiny::reactive({
    selected <- selector$selected()
    selector$revision()
    if (is.null(selected)) {
      return(NULL)
    }
    session_stage_status(selected$path)
  })

  shiny::observeEvent(selector$selected(), {
    stages <- stage_table()
    if (is.null(stages)) {
      return()
    }
    incomplete <- which(stages$status != "complete")
    next_stage <- if (length(incomplete) > 0L) {
      stages$id[[incomplete[[1L]]]]
    } else {
      stages$id[[nrow(stages)]]
    }
    shiny::updateRadioButtons(session, "active_stage", selected = next_stage)
  }, ignoreNULL = TRUE)

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
        shiny::span("stages complete")
      )
    )
  })

  output$stage_content <- shiny::renderUI({
    selected <- selector$selected()
    stages <- stage_table()
    if (is.null(selected) || is.null(stages)) {
      return(NULL)
    }

    active <- input$active_stage %||% stages$id[[1L]]
    if (!active %in% stages$id) {
      active <- stages$id[[1L]]
    }
    stage <- stages[match(active, stages$id), , drop = FALSE]
    outputs <- session_output_files(selected$path)

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
