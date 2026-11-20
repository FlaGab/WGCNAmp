sessionSelectorUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::div(
    class = "session-selector",
    shiny::div(
      class = "eyebrow",
      "Analysis session"
    ),
    shiny::selectizeInput(
      ns("session_name"),
      label = NULL,
      choices = NULL,
      options = list(placeholder = "Select a session")
    ),
    shiny::div(
      class = "session-actions",
      shiny::actionButton(
        ns("refresh"),
        "Refresh",
        icon = shiny::icon("rotate"),
        class = "btn-outline-secondary"
      ),
      shiny::actionButton(
        ns("new_session"),
        "New session",
        icon = shiny::icon("plus"),
        class = "btn-primary"
      )
    )
  )
}

sessionSelectorServer <- function(id, project_root) {
  shiny::moduleServer(id, function(input, output, session) {
    sessions <- shiny::reactiveVal(discover_analysis_sessions(project_root))
    revision <- shiny::reactiveVal(0L)

    refresh_sessions <- function(select = NULL) {
      discovered <- discover_analysis_sessions(project_root)
      sessions(discovered)

      choices <- stats::setNames(discovered$name, discovered$name)
      current <- select %||% shiny::isolate(input$session_name)
      if (is.null(current) || length(current) == 0L || !current %in% discovered$name) {
        current <- if (nrow(discovered) > 0L) discovered$name[[1L]] else character()
      }

      shiny::updateSelectizeInput(
        session,
        "session_name",
        choices = choices,
        selected = current,
        server = TRUE
      )
      revision(revision() + 1L)
    }

    shiny::observeEvent(input$refresh, {
      refresh_sessions()
      shiny::showNotification("Session list refreshed.", type = "message")
    })

    shiny::observeEvent(input$new_session, {
      shiny::showModal(shiny::modalDialog(
        title = "Create analysis session",
        shiny::textInput(
          session$ns("new_session_name"),
          "Session name",
          placeholder = "e.g. drought_rnaseq"
        ),
        shiny::p(
          class = "text-secondary",
          "The app will create the input, output, plot, and log directories and copy the parameter template."
        ),
        footer = shiny::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton(
            session$ns("confirm_create"),
            "Create session",
            class = "btn-primary"
          )
        ),
        easyClose = TRUE
      ))
    })

    shiny::observeEvent(input$confirm_create, {
      tryCatch({
        name <- validate_session_name(input$new_session_name)
        create_analysis_session(name, project_root)
        shiny::removeModal()
        refresh_sessions(select = name)
        shiny::showNotification(
          paste0("Session '", name, "' created."),
          type = "message",
          duration = 5
        )
      }, error = function(error) {
        shiny::showNotification(
          conditionMessage(error),
          type = "error",
          duration = NULL
        )
      })
    })

    shiny::observeEvent(input$refresh, {
      refresh_sessions()
    }, ignoreInit = TRUE)

    session$onFlushed(function() {
      refresh_sessions()
    }, once = TRUE)

    selected_session <- shiny::reactive({
      revision()
      selected <- input$session_name
      table <- sessions()
      if (is.null(selected) || !nzchar(selected) || !selected %in% table$name) {
        return(NULL)
      }
      as.list(table[match(selected, table$name), , drop = FALSE])
    })

    list(
      selected = selected_session,
      revision = shiny::reactive(revision())
    )
  })
}

`%||%` <- function(left, right) {
  if (is.null(left) ||
      length(left) == 0L ||
      (length(left) == 1L && !nzchar(left))) {
    right
  } else {
    left
  }
}
