.microcyte_attach <- function(names, envir = parent.frame()) {
  if (!requireNamespace("MicroCyte", quietly = TRUE)) {
    stop("Install the MicroCyte package before sourcing this compatibility script.", call. = FALSE)
  }
  for (name in names) {
    assign(name, getExportedValue("MicroCyte", name), envir = envir)
  }
  invisible(names)
}

.microcyte_attach("ergodic_slicing")
rm(.microcyte_attach)
