.microcyte_attach <- function(names, envir = parent.frame()) {
  if (!requireNamespace("MicroCyte", quietly = TRUE)) {
    stop("Install the MicroCyte package before sourcing this compatibility script.", call. = FALSE)
  }
  for (name in names) {
    assign(name, getExportedValue("MicroCyte", name), envir = envir)
  }
  invisible(names)
}

.microcyte_attach(c(
  "sirmixaplot", "joiner", "grapho", "modthequad", "cull",
  "get_density", "color_density", "normalizer", "compensate", "re_map",
  "run_ca", "nearest_n", "cluster_count"
))
rm(.microcyte_attach)
