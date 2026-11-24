arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
project_root <- if (length(file_argument) == 1L) {
  normalizePath(dirname(sub("^--file=", "", file_argument)), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

user_arguments <- commandArgs(trailingOnly = TRUE)
launch_browser <- !("--no-browser" %in% user_arguments)
port_argument <- grep("^--port=", user_arguments, value = TRUE)
port <- if (length(port_argument) == 1L) {
  as.integer(sub("^--port=", "", port_argument))
} else {
  getOption("shiny.port", NULL)
}

Sys.setenv(WGCNAMP_PROJECT_ROOT = project_root)

shiny::runApp(
  appDir = file.path(project_root, "app"),
  launch.browser = launch_browser,
  port = port
)
