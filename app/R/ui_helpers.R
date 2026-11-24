status_icon <- function(status) {
  switch(status, complete = "✓", stale = "⚠", ready = "●", blocked = "○", "○")
}

status_label <- function(status) {
  switch(
    status,
    complete = "Complete",
    stale = "Stale — rerun",
    ready = "Ready",
    blocked = "Waiting",
    status
  )
}

status_badge <- function(status, reason = NULL) {
  shiny::span(
    class = paste("status-badge", paste0("status-", status)),
    title = reason,
    status_label(status)
  )
}

session_stage_badge <- function(session_path, stage_id) {
  stages <- session_stage_status(session_path)
  stage <- stages[match(stage_id, stages$id), , drop = FALSE]
  status <- stage$status[[1L]]
  reason <- if (identical(status, "stale")) {
    stage$stale_reason[[1L]]
  } else {
    NULL
  }
  status_badge(status, reason)
}
