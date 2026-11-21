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
