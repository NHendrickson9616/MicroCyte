ns <- loadNamespace("MicroCyte")

expected <- c(
  "foci_caller",
  "explore",
  "icellate",
  "sirmixaplot",
  "ternaplot"
)

stopifnot(all(vapply(expected, exists, logical(1L), envir = ns, inherits = FALSE)))

sample_data <- data.frame(
  x = c(1, 2, 3),
  y = c(2, 4, 8),
  Number = 1:3
)
culled <- MicroCyte::cull(sample_data, parameter = "x", direction = "b", value = 2)
stopifnot(nrow(culled) == 1L, identical(culled$x, 1))

analysis_formals <- lapply(expected, function(name) formals(get(name, envir = ns)))
ambient_defaults <- vapply(
  analysis_formals,
  function(args) {
    any(vapply(args, identical, logical(1L), quote(cells))) ||
      any(vapply(args, identical, logical(1L), quote(datum)))
  },
  logical(1L)
)
stopifnot(!any(ambient_defaults))
